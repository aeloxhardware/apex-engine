// Bindings pybind11 - ECHAFAUDAGE DE MIGRATION, pas un livrable.
//
// La cible est un moteur 100 % C++ ; ces bindings existent pour que les 305
// tests pytest servent d'ORACLE D'EQUIVALENCE pendant tout le portage (le
// Python et le C++ coexistent dans le meme processus, donc un test peut
// appeler les deux et les comparer). Ils seront retires a la bascule finale.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

#include "apex/bvh.hpp"
#include "apex/field_solver.hpp"
#include "apex/kdtree.hpp"
#include "apex/kinematics.hpp"
#include "apex/interpolate.hpp"
#include "apex/mesh_distance.hpp"
#include "apex/mesh_intersect.hpp"
#include "apex/mesh_repair.hpp"
#include "apex/printed_index.hpp"
#include "apex/support.hpp"
#include "apex/tet_topology.hpp"
#include "apex/types.hpp"

namespace py = pybind11;

namespace {

using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
using IntArray = py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;
using Int32Array = py::array_t<std::int32_t, py::array::c_style | py::array::forcecast>;

// Lit un point (3,) depuis un tableau numpy quelconque (converti si besoin).
apex::Vec3 read_vec3(const DoubleArray& array, const char* name) {
    if (array.size() != 3) {
        throw std::invalid_argument(std::string(name) + " doit contenir exactement 3 valeurs");
    }
    const double* data = array.data();
    return apex::Vec3{data[0], data[1], data[2]};
}

// Enveloppe Python de PrintedGrid : elle ajoute la conversion numpy et les
// deux formes de sortie attendues par les cinq sites d'appel.
class PyPrintedGrid {
public:
    PyPrintedGrid(double cell_size, double radius, int mode)
        : grid_(cell_size, radius,
                mode == 1 ? apex::QueryMode::Stamped : apex::QueryMode::Faithful) {}

    void append(const DoubleArray& start, const DoubleArray& end) {
        grid_.append(read_vec3(start, "start"), read_vec3(end, "end"));
    }

    // Insertion en LOT. Le chemin historique appelait `append` une fois par
    // capsule depuis Python - 744 324 franchissements de frontiere sur
    // death.stl, plus un `_segment_aabb` numpy par capsule : ~4 s mesurees
    // (1,79 + 1,19 + 1,01) pour une operation triviale. La granularite de
    // l'appel coutait plus cher que le calcul lui-meme.
    void append_batch(const DoubleArray& starts, const DoubleArray& ends) {
        if (starts.ndim() != 2 || starts.shape(1) != 3) {
            throw std::invalid_argument("starts doit etre de forme (N,3)");
        }
        if (ends.shape(0) != starts.shape(0) || ends.ndim() != 2 || ends.shape(1) != 3) {
            throw std::invalid_argument("ends doit avoir la meme forme que starts");
        }
        const auto count = static_cast<std::int64_t>(starts.shape(0));
        const double* s = starts.data();
        const double* e = ends.data();

        py::gil_scoped_release release;
        for (std::int64_t i = 0; i < count; ++i) {
            grid_.append(apex::Vec3{s[3 * i], s[3 * i + 1], s[3 * i + 2]},
                         apex::Vec3{e[3 * i], e[3 * i + 1], e[3 * i + 2]});
        }
    }

    // Motif A : step9 `_check_printed_part`, step8b `_path_collides`. Ils
    // bouclent capsule par capsule cote Python et ont besoin des objets, donc
    // des indices dans leur propre liste `_segments`.
    py::array_t<std::int32_t> candidates_indices(const DoubleArray& lo,
                                                 const DoubleArray& hi) {
        query(lo, hi);
        py::array_t<std::int32_t> out(static_cast<py::ssize_t>(buffer_.size()));
        if (!buffer_.empty()) {
            std::copy(buffer_.begin(), buffer_.end(), out.mutable_data());
        }
        return out;
    }

    // Motif B : step9b `_has_support_below` / `_has_support_in_any_direction`
    // et step5 `_layer_support_mask`. Ils construisaient deux tableaux (C,3)
    // par comprehension sur les objets candidats - 41 s de `numpy.array`
    // mesurees. On les leur rend directement.
    py::tuple candidates_endpoints(const DoubleArray& lo, const DoubleArray& hi) {
        query(lo, hi);
        const auto count = static_cast<py::ssize_t>(buffer_.size());

        py::array_t<double> starts({count, static_cast<py::ssize_t>(3)});
        py::array_t<double> ends({count, static_cast<py::ssize_t>(3)});
        double* s = starts.mutable_data();
        double* e = ends.mutable_data();

        for (py::ssize_t i = 0; i < count; ++i) {
            const apex::Vec3& a = grid_.start_of(buffer_[static_cast<std::size_t>(i)]);
            const apex::Vec3& b = grid_.end_of(buffer_[static_cast<std::size_t>(i)]);
            s[3 * i + 0] = a.x;
            s[3 * i + 1] = a.y;
            s[3 * i + 2] = a.z;
            e[3 * i + 0] = b.x;
            e[3 * i + 1] = b.y;
            e[3 * i + 2] = b.z;
        }
        return py::make_tuple(std::move(starts), std::move(ends));
    }

    // FUSION requete spatiale + test de support. C'est le point d'entree qui
    // rapporte : il ne materialise JAMAIS les ~690 candidats cote Python.
    // L'ancien chemin allouait deux tableaux (C,3) par appel (3,95 s mesurees
    // rien que dans `candidates_endpoints`), puis Python executait ~15
    // operations numpy dessus pour ne rendre qu'un booleen (11,97 s + ~10 s de
    // dispatch numpy). Ici, les extremites sont lues directement dans la
    // grille et la boucle s'arrete au premier appui trouve.
    bool has_support(const DoubleArray& lo, const DoubleArray& hi, const DoubleArray& position,
                     const DoubleArray& directions, double horizontal_radius,
                     double vertical_height) {
        query(lo, hi);
        if (buffer_.empty()) {
            return false;
        }
        if (directions.ndim() != 2 || directions.shape(1) != 3) {
            throw std::invalid_argument("directions doit etre de forme (D,3)");
        }

        // Extremites des candidats, compactees une fois par requete dans un
        // tampon reutilise (pas d'allocation par appel).
        starts_buffer_.resize(buffer_.size() * 3);
        ends_buffer_.resize(buffer_.size() * 3);
        for (std::size_t i = 0; i < buffer_.size(); ++i) {
            const apex::Vec3& a = grid_.start_of(buffer_[i]);
            const apex::Vec3& b = grid_.end_of(buffer_[i]);
            starts_buffer_[3 * i + 0] = a.x;
            starts_buffer_[3 * i + 1] = a.y;
            starts_buffer_[3 * i + 2] = a.z;
            ends_buffer_[3 * i + 0] = b.x;
            ends_buffer_[3 * i + 1] = b.y;
            ends_buffer_[3 * i + 2] = b.z;
        }

        const apex::Vec3 point = read_vec3(position, "position");
        const auto direction_count = static_cast<std::int64_t>(directions.shape(0));
        const auto candidate_count = static_cast<std::int64_t>(buffer_.size());
        const double* directions_data = directions.data();

        py::gil_scoped_release release;
        return apex::point_has_support_in_any_direction(
            point, directions_data, direction_count, starts_buffer_.data(), ends_buffer_.data(),
            candidate_count, horizontal_radius, vertical_height);
    }

    std::int64_t size() const { return grid_.size(); }
    void set_mode(int mode) {
        grid_.set_mode(mode == 1 ? apex::QueryMode::Stamped : apex::QueryMode::Faithful);
    }
    int mode() const { return grid_.mode() == apex::QueryMode::Stamped ? 1 : 0; }

    py::dict stats() const {
        py::dict d;
        d["visits"] = grid_.visit_count();
        d["uniques"] = grid_.unique_count();
        d["queries"] = grid_.query_count();
        const auto uniques = grid_.unique_count();
        d["duplication_factor"] =
            uniques == 0 ? 0.0
                         : static_cast<double>(grid_.visit_count()) / static_cast<double>(uniques);
        return d;
    }
    void reset_stats() const { grid_.reset_stats(); }

private:
    void query(const DoubleArray& lo, const DoubleArray& hi) {
        apex::AABB box;
        box.lo = read_vec3(lo, "min_corner");
        box.hi = read_vec3(hi, "max_corner");
        grid_.candidates(box, buffer_);
    }

    apex::PrintedGrid grid_;
    std::vector<std::int32_t> buffer_;
    // Tampons reutilises entre requetes - une allocation par appel serait le
    // cout meme que ce chemin existe pour supprimer.
    std::vector<double> starts_buffer_;
    std::vector<double> ends_buffer_;
};

// Enveloppe Python de KDTree : conversion numpy + les deux formes de sortie
// attendues (`(distances, indices)` comme `cKDTree.query`, liste d'indices
// comme `cKDTree.query_ball_point`).
class PyKDTree {
public:
    explicit PyKDTree(const DoubleArray& points) {
        if (points.ndim() != 2 || points.shape(1) != 3) {
            throw std::invalid_argument("points doit etre de forme (N,3)");
        }
        const auto count = static_cast<std::int64_t>(points.shape(0));
        py::gil_scoped_release release;
        tree_ = std::make_unique<apex::KDTree>(points.data(), count);
    }

    // Equivalent de `cKDTree.query(x)` avec k=1 : rend (distances (M,),
    // indices (M,)). L'appelant Python gere le cas x de forme (3,).
    py::tuple query_nn(const DoubleArray& queries) const {
        if (queries.ndim() != 2 || queries.shape(1) != 3) {
            throw std::invalid_argument("queries doit etre de forme (M,3)");
        }
        const auto count = static_cast<py::ssize_t>(queries.shape(0));
        py::array_t<double> distances(count);
        py::array_t<std::int64_t> indices(count);
        {
            py::gil_scoped_release release;
            tree_->query_nn(queries.data(), static_cast<std::int64_t>(count),
                            indices.mutable_data(), distances.mutable_data());
        }
        return py::make_tuple(std::move(distances), std::move(indices));
    }

    // Equivalent de `cKDTree.query_ball_point(point, r)` : indices dans un
    // rayon, ordre non specifie.
    std::vector<std::int32_t> query_radius(const DoubleArray& point, double radius) const {
        std::vector<std::int32_t> out;
        const apex::Vec3 p = read_vec3(point, "point");
        py::gil_scoped_release release;
        tree_->query_radius(p, radius, out);
        return out;
    }

    std::int64_t size() const { return tree_->size(); }

private:
    std::unique_ptr<apex::KDTree> tree_;
};

// Enveloppe Python du BVH de triangles. Rend le meme TRIPLET que
// `RayMeshIntersector.intersects_location` : (locations, index_ray,
// index_tri) - dans cet ordre, qui est celui de trimesh et non l'ordre
// naturel (index d'abord).
class PyTriangleBVH {
public:
    PyTriangleBVH(const DoubleArray& triangles, const DoubleArray& normals) {
        if (triangles.ndim() != 3 || triangles.shape(1) != 3 || triangles.shape(2) != 3) {
            throw std::invalid_argument("triangles doit etre de forme (N,3,3)");
        }
        if (normals.ndim() != 2 || normals.shape(1) != 3 ||
            normals.shape(0) != triangles.shape(0)) {
            throw std::invalid_argument("normals doit etre de forme (N,3), une par triangle");
        }
        const auto count = static_cast<std::int64_t>(triangles.shape(0));
        py::gil_scoped_release release;
        bvh_ = std::make_unique<apex::TriangleBVH>(triangles.data(), normals.data(), count);
    }

    py::tuple intersects_location(const DoubleArray& origins, const DoubleArray& directions) const {
        if (origins.ndim() != 2 || origins.shape(1) != 3) {
            throw std::invalid_argument("ray_origins doit etre de forme (M,3)");
        }
        if (directions.ndim() != 2 || directions.shape(1) != 3 ||
            directions.shape(0) != origins.shape(0)) {
            throw std::invalid_argument("ray_directions doit avoir la meme forme que ray_origins");
        }

        std::vector<apex::RayHit> hits;
        {
            py::gil_scoped_release release;
            hits = bvh_->intersects_location(origins.data(), directions.data(),
                                             static_cast<std::int64_t>(origins.shape(0)));
        }

        const auto count = static_cast<py::ssize_t>(hits.size());
        py::array_t<double> locations({count, static_cast<py::ssize_t>(3)});
        py::array_t<std::int64_t> index_ray(count);
        py::array_t<std::int64_t> index_tri(count);
        double* loc = locations.mutable_data();
        std::int64_t* iray = index_ray.mutable_data();
        std::int64_t* itri = index_tri.mutable_data();
        for (py::ssize_t i = 0; i < count; ++i) {
            const apex::RayHit& h = hits[static_cast<std::size_t>(i)];
            loc[3 * i + 0] = h.location.x;
            loc[3 * i + 1] = h.location.y;
            loc[3 * i + 2] = h.location.z;
            iray[i] = h.ray_index;
            itri[i] = h.triangle_index;
        }
        return py::make_tuple(std::move(locations), std::move(index_ray), std::move(index_tri));
    }

    std::int64_t size() const { return bvh_->size(); }

private:
    std::unique_ptr<apex::TriangleBVH> bvh_;
};

// Motif de sparsite FEM, expose comme OBJET plutot que recalcule a chaque
// matrice : les deux matrices de step4 (rigidite et terme de direction)
// partagent exactement le meme, qui ne depend que de la topologie du
// maillage. Le Python le reconstruisait deux fois via `coo_matrix.tocsr()`.
class PyFemSparsity {
public:
    PyFemSparsity(const IntArray& tetrahedra, std::int64_t node_count) {
        if (tetrahedra.ndim() != 2 || tetrahedra.shape(1) != 4) {
            throw std::invalid_argument("tetrahedra doit etre de forme (T,4)");
        }
        tet_count_ = static_cast<std::int64_t>(tetrahedra.shape(0));
        node_count_ = node_count;
        py::gil_scoped_release release;
        pattern_ = apex::build_fem_sparsity(tetrahedra.data(), tet_count_, node_count);
    }

    py::array_t<std::int32_t> indptr() const {
        py::array_t<std::int32_t> out(static_cast<py::ssize_t>(pattern_.indptr.size()));
        std::copy(pattern_.indptr.begin(), pattern_.indptr.end(), out.mutable_data());
        return out;
    }

    py::array_t<std::int32_t> indices() const {
        py::array_t<std::int32_t> out(static_cast<py::ssize_t>(pattern_.indices.size()));
        std::copy(pattern_.indices.begin(), pattern_.indices.end(), out.mutable_data());
        return out;
    }

    // `shape_gradients` et `scale` sont indexes COMME `active` quand il est
    // fourni (le Python les a deja sous-ensembles), sinon par tetraedre.
    py::array_t<double> assemble_gram(const IntArray& tetrahedra,
                                      const DoubleArray& shape_gradients, const DoubleArray& scale,
                                      const py::object& active) const {
        Int32Array active_array;
        const std::int32_t* active_ptr = nullptr;
        std::int64_t active_count = 0;
        if (!active.is_none()) {
            active_array = active.cast<Int32Array>();
            active_ptr = active_array.data();
            active_count = static_cast<std::int64_t>(active_array.size());
        }

        py::array_t<double> data(static_cast<py::ssize_t>(pattern_.indices.size()));
        {
            py::gil_scoped_release release;
            apex::assemble_gram_matrix(tetrahedra.data(), shape_gradients.data(), scale.data(),
                                       active_ptr, active_count, tet_count_, pattern_,
                                       data.mutable_data());
        }
        return data;
    }

    std::int64_t node_count() const { return node_count_; }

private:
    apex::SparsityPattern pattern_;
    std::int64_t tet_count_ = 0;
    std::int64_t node_count_ = 0;
};

}  // namespace

// Enveloppe Python de MeshDistanceIndex (T9, 2026-09-30) : ADAPTATEUR numpy ->
// vues brutes. Le coeur (mesh_distance.cpp) ignore pybind11 ; ce binding sera
// retire a la bascule, les etapes C++ amont remplissant MeshView directement.
class PyMeshDistance {
public:
    PyMeshDistance(const DoubleArray& vertices, const py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>& faces)
        : vertices_(vertices), faces_(faces) {
        if (vertices.ndim() != 2 || vertices.shape(1) != 3) throw std::invalid_argument("vertices : (V,3)");
        if (faces.ndim() != 2 || faces.shape(1) != 3) throw std::invalid_argument("faces : (F,3)");
        apex::MeshView view;
        view.vertices = vertices_.data();
        view.n_vertices = static_cast<std::size_t>(vertices_.shape(0));
        view.faces = faces_.data();
        view.n_faces = static_cast<std::size_t>(faces_.shape(0));
        index_ = std::make_unique<apex::MeshDistanceIndex>(view);
    }
    py::tuple query(const DoubleArray& points, int backend, int threads) const {
        if (points.ndim() != 2 || points.shape(1) != 3) throw std::invalid_argument("points : (N,3)");
        const auto n = static_cast<py::ssize_t>(points.shape(0));
        py::array_t<double> dist(n);
        py::array_t<double> closest({n, static_cast<py::ssize_t>(3)});
        {
            py::gil_scoped_release release;
            index_->distances(points.data(), static_cast<std::size_t>(n), dist.mutable_data(),
                              closest.mutable_data(), static_cast<apex::DistanceBackend>(backend), threads);
        }
        return py::make_tuple(dist, closest);
    }
private:
    DoubleArray vertices_;
    py::array_t<std::int64_t, py::array::c_style | py::array::forcecast> faces_;
    std::unique_ptr<apex::MeshDistanceIndex> index_;
};

PYBIND11_MODULE(apex_core, m) {
    m.doc() = "Noyau C++ du moteur Apex Slicer (portage en cours)";
    m.attr("__version__") = "0.1.0";

    m.def("selftest_roundtrip",
          [](const DoubleArray& a) {
              // Verification minimale du pas 0 : la chaine numpy -> C++ ->
              // numpy rend exactement le tableau d'entree.
              py::array_t<double> out(a.size());
              std::copy(a.data(), a.data() + a.size(), out.mutable_data());
              return out;
          },
          "Aller-retour numpy trivial, pour valider la chaine de build.");

    m.def(
        "resolve_axis_angles_batch",
        [](const DoubleArray& desired_directions, const DoubleArray& positions,
           const py::list& axes_spec, const DoubleArray& reference_direction,
           const py::object& platform_origin, double singularity_radius_epsilon,
           const py::object& initial_previous_angles) {
            if (desired_directions.ndim() != 2 || desired_directions.shape(1) != 3) {
                throw std::invalid_argument("desired_directions doit etre de forme (N,3)");
            }
            if (positions.ndim() != 2 || positions.shape(1) != 3) {
                throw std::invalid_argument("positions doit etre de forme (N,3)");
            }
            if (positions.shape(0) != desired_directions.shape(0)) {
                throw std::invalid_argument("positions et desired_directions de tailles differentes");
            }

            std::vector<apex::AxisSpec> axes;
            axes.reserve(axes_spec.size());
            for (const auto& item : axes_spec) {
                const auto tuple = item.cast<py::tuple>();
                const auto axis = tuple[0].cast<DoubleArray>();
                axes.push_back(apex::AxisSpec{read_vec3(axis, "axis"), tuple[1].cast<double>(),
                                              tuple[2].cast<double>()});
            }

            apex::Vec3 origin_value{};
            const apex::Vec3* origin_ptr = nullptr;
            if (!platform_origin.is_none()) {
                origin_value = read_vec3(platform_origin.cast<DoubleArray>(), "platform_origin");
                origin_ptr = &origin_value;
            }

            std::vector<double> previous_storage;
            const double* previous_ptr = nullptr;
            if (!initial_previous_angles.is_none()) {
                const auto arr = initial_previous_angles.cast<DoubleArray>();
                if (static_cast<std::size_t>(arr.size()) != axes.size()) {
                    throw std::invalid_argument(
                        "initial_previous_angles doit avoir un angle par axe");
                }
                previous_storage.assign(arr.data(), arr.data() + arr.size());
                previous_ptr = previous_storage.data();
            }

            const auto count = static_cast<std::int64_t>(desired_directions.shape(0));
            const apex::Vec3 reference_unit =
                read_vec3(reference_direction, "reference_direction");
            apex::IKBatchResult res;
            {
                // La resolution ne touche a aucun objet Python : liberer le GIL
                // permet a un autre thread de progresser pendant ce temps.
                py::gil_scoped_release release;
                res = apex::resolve_axis_angles_batch(
                    desired_directions.data(), positions.data(), count, axes, reference_unit,
                    origin_ptr, singularity_radius_epsilon, previous_ptr);
            }

            const auto axis_count = static_cast<py::ssize_t>(axes.size());
            py::array_t<double> angles({static_cast<py::ssize_t>(count), axis_count});
            py::array_t<double> residuals(static_cast<py::ssize_t>(count));
            if (!res.angles.empty()) {
                std::copy(res.angles.begin(), res.angles.end(), angles.mutable_data());
            }
            if (!res.residuals.empty()) {
                std::copy(res.residuals.begin(), res.residuals.end(), residuals.mutable_data());
            }
            return py::make_tuple(std::move(angles), std::move(residuals), res.fallback_index);
        },
        py::arg("desired_directions"), py::arg("positions"), py::arg("axes"),
        py::arg("reference_direction"), py::arg("platform_origin"),
        py::arg("singularity_radius_epsilon"), py::arg("initial_previous_angles") = py::none(),
        "Resolution IK N axes, sequentielle et REPRENABLE. Rend (angles (N,A), "
        "residuals (N,) en radians, fallback_index). fallback_index >= 0 signale que "
        "la paire externe n'a pas pu etre resolue analytiquement a CE point : les "
        "lignes precedentes sont valides, l'appelant resout ce point avec la grille "
        "Python, puis rappelle a partir du suivant en passant initial_previous_angles.");

    m.def(
        "tet_face_topology",
        [](const IntArray& tetrahedra) {
            if (tetrahedra.ndim() != 2 || tetrahedra.shape(1) != 4) {
                throw std::invalid_argument("tetrahedra doit etre de forme (N,4)");
            }
            const auto count = static_cast<std::int64_t>(tetrahedra.shape(0));
            apex::TetFaceTopology topo;
            {
                py::gil_scoped_release release;
                topo = apex::compute_tet_face_topology(tetrahedra.data(), count);
            }

            const auto boundary_count =
                static_cast<py::ssize_t>(topo.boundary_face_tets.size());
            py::array_t<std::int32_t> nodes({boundary_count, static_cast<py::ssize_t>(3)});
            py::array_t<std::int32_t> owners(boundary_count);
            py::array_t<std::int32_t> pairs(
                {static_cast<py::ssize_t>(topo.adjacency_pairs.size() / 2),
                 static_cast<py::ssize_t>(2)});
            if (!topo.boundary_face_nodes.empty()) {
                std::copy(topo.boundary_face_nodes.begin(), topo.boundary_face_nodes.end(),
                          nodes.mutable_data());
                std::copy(topo.boundary_face_tets.begin(), topo.boundary_face_tets.end(),
                          owners.mutable_data());
            }
            if (!topo.adjacency_pairs.empty()) {
                std::copy(topo.adjacency_pairs.begin(), topo.adjacency_pairs.end(),
                          pairs.mutable_data());
            }
            return py::make_tuple(std::move(nodes), std::move(owners), std::move(pairs));
        },
        py::arg("tetrahedra"),
        "Rend (boundary_face_nodes (K,3), boundary_face_tets (K,), adjacency_pairs (M,2)). "
        "Les faces frontiere sont dans l'ordre d'INSERTION, comme le dict Python.");

    m.def(
        "tet_unique_edges",
        [](const IntArray& tetrahedra) {
            if (tetrahedra.ndim() != 2 || tetrahedra.shape(1) != 4) {
                throw std::invalid_argument("tetrahedra doit etre de forme (N,4)");
            }
            const auto count = static_cast<std::int64_t>(tetrahedra.shape(0));
            std::vector<std::int32_t> edges;
            {
                py::gil_scoped_release release;
                edges = apex::compute_unique_edges(tetrahedra.data(), count);
            }
            py::array_t<std::int32_t> out(
                {static_cast<py::ssize_t>(edges.size() / 2), static_cast<py::ssize_t>(2)});
            if (!edges.empty()) {
                std::copy(edges.begin(), edges.end(), out.mutable_data());
            }
            return out;
        },
        py::arg("tetrahedra"), "Aretes uniques (a<b) du maillage tetraedrique, forme (M,2).");

    m.def("global_stats", []() {
        const apex::GlobalStats& s = apex::global_stats();
        py::dict d;
        d["visits"] = s.visits;
        d["uniques"] = s.uniques;
        d["queries"] = s.queries;
        d["appends"] = s.appends;
        // Rapport visites/uniques = facteur de duplication de la grille : une
        // capsule est indexee dans toutes les cellules que son AABB touche,
        // donc revisitee une fois par cellule a chaque requete. C'est la
        // mesure qui distingue le defaut ALGORITHMIQUE du seul cout du
        // langage (limite #27/#30 de CLAUDE.md, jamais chiffree jusqu'ici).
        d["duplication_factor"] =
            s.uniques == 0 ? 0.0 : static_cast<double>(s.visits) / static_cast<double>(s.uniques);
        return d;
    }, "Compteurs cumules sur toutes les instances de PrintedGrid du processus.");

    m.def("reset_global_stats", &apex::reset_global_stats);

    py::class_<PyPrintedGrid>(m, "PrintedGrid")
        .def(py::init<double, double, int>(), py::arg("cell_size"), py::arg("radius"),
             py::arg("mode") = 1)
        .def("append", &PyPrintedGrid::append, py::arg("start"), py::arg("end"))
        .def("append_batch", &PyPrintedGrid::append_batch, py::arg("starts"), py::arg("ends"),
             "Insere N capsules en un seul appel (voir append_batch, C++).")
        .def("candidates_indices", &PyPrintedGrid::candidates_indices, py::arg("min_corner"),
             py::arg("max_corner"))
        .def("candidates_endpoints", &PyPrintedGrid::candidates_endpoints,
             py::arg("min_corner"), py::arg("max_corner"))
        .def("has_support", &PyPrintedGrid::has_support, py::arg("min_corner"),
             py::arg("max_corner"), py::arg("position"), py::arg("directions"),
             py::arg("horizontal_radius"), py::arg("vertical_height"),
             "Requete spatiale ET test de support fusionnes, avec arret anticipe. "
             "Ne materialise jamais les candidats cote Python.")
        .def("stats", &PyPrintedGrid::stats)
        .def("reset_stats", &PyPrintedGrid::reset_stats)
        .def_property("mode", &PyPrintedGrid::mode, &PyPrintedGrid::set_mode)
        .def("__len__", &PyPrintedGrid::size);

    py::class_<PyMeshDistance>(m, "MeshDistance")
        .def(py::init<const DoubleArray&, const py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>&>(),
             py::arg("vertices"), py::arg("faces"))
        .def("query", &PyMeshDistance::query, py::arg("points"), py::arg("backend") = 1, py::arg("threads") = 0,
             "Distance exacte au maillage et point le plus proche (backend 0 serie, 1 threads, 2 Vulkan).");
    m.def("mesh_distance_backend", [](int requested) {
        std::string reason;
        const auto b = apex::MeshDistanceIndex::resolve_backend(static_cast<apex::DistanceBackend>(requested), &reason);
        return py::make_tuple(static_cast<int>(b), reason);
    });

    py::class_<PyKDTree>(m, "KDTree")
        .def(py::init<const DoubleArray&>(), py::arg("points"))
        .def("query_nn", &PyKDTree::query_nn, py::arg("queries"),
             "1-NN pour un lot de points (M,3). Rend (distances (M,), indices (M,)).")
        .def("query_radius", &PyKDTree::query_radius, py::arg("point"), py::arg("radius"),
             "Indices des points a distance euclidienne <= radius de point (3,).")
        .def("__len__", &PyKDTree::size);

    m.def(
        "blend_tet_directions",
        [](const DoubleArray& tet_centers, const DoubleArray& source_points,
           const Int32Array& source_offsets, const DoubleArray& source_directions,
           const py::object& platform_level, double radius, double power, double min_norm,
           double epsilon) {
            if (tet_centers.ndim() != 2 || tet_centers.shape(1) != 3) {
                throw std::invalid_argument("tet_centers doit etre de forme (T,3)");
            }
            if (source_directions.ndim() != 2 || source_directions.shape(1) != 3) {
                throw std::invalid_argument("source_directions doit etre de forme (S,3)");
            }
            const auto tet_count = static_cast<std::int64_t>(tet_centers.shape(0));
            const auto region_count = static_cast<std::int64_t>(source_offsets.size()) - 1;
            const bool has_platform = !platform_level.is_none();
            if (source_directions.shape(0) != region_count + (has_platform ? 1 : 0)) {
                throw std::invalid_argument(
                    "source_directions doit compter une direction par region, plateau inclus");
            }

            apex::BlendSources sources;
            sources.points = source_points.data();
            sources.offsets = source_offsets.data();
            sources.region_count = region_count;
            sources.directions = source_directions.data();
            sources.has_platform = has_platform;
            sources.platform_level = has_platform ? platform_level.cast<double>() : 0.0;

            apex::BlendParams blend_params;
            blend_params.radius = radius;
            blend_params.power = power;
            blend_params.min_norm = min_norm;
            blend_params.epsilon = epsilon;

            py::array_t<double> out(
                {static_cast<py::ssize_t>(tet_count), static_cast<py::ssize_t>(3)});
            {
                py::gil_scoped_release release;
                apex::blend_tet_directions(tet_centers.data(), tet_count, sources, blend_params,
                                           out.mutable_data());
            }
            return out;
        },
        py::arg("tet_centers"), py::arg("source_points"), py::arg("source_offsets"),
        py::arg("source_directions"), py::arg("platform_level"), py::arg("radius"),
        py::arg("power"), py::arg("min_norm"), py::arg("epsilon"),
        "Melange spatial des directions de gradient, avec elagage AABB par source.");

    m.def(
        "compute_shape_gradients",
        [](const DoubleArray& tet_vertices, const IntArray& tetrahedra) {
            if (tet_vertices.ndim() != 2 || tet_vertices.shape(1) != 3) {
                throw std::invalid_argument("tet_vertices doit etre de forme (P,3)");
            }
            if (tetrahedra.ndim() != 2 || tetrahedra.shape(1) != 4) {
                throw std::invalid_argument("tetrahedra doit etre de forme (T,4)");
            }
            const auto tet_count = static_cast<std::int64_t>(tetrahedra.shape(0));
            py::array_t<double> gradients({static_cast<py::ssize_t>(tet_count),
                                           static_cast<py::ssize_t>(3),
                                           static_cast<py::ssize_t>(4)});
            py::array_t<double> volumes(static_cast<py::ssize_t>(tet_count));
            {
                py::gil_scoped_release release;
                apex::compute_shape_gradients(tet_vertices.data(), tetrahedra.data(), tet_count,
                                              gradients.mutable_data(), volumes.mutable_data());
            }
            return py::make_tuple(std::move(gradients), std::move(volumes));
        },
        py::arg("tet_vertices"), py::arg("tetrahedra"),
        "Gradients des fonctions de forme P1 (T,3,4) et volumes (T,).");

    m.def(
        "solve_linear_system",
        [](const Int32Array& indptr, const Int32Array& indices, const DoubleArray& data,
           const DoubleArray& b, bool direct, double tolerance, int max_iterations,
           const py::object& initial_guess) {
            const auto node_count = static_cast<std::int64_t>(indptr.size()) - 1;
            if (b.size() != node_count) {
                throw std::invalid_argument("b doit avoir une valeur par noeud");
            }
            DoubleArray guess_array;
            const double* guess = nullptr;
            if (!initial_guess.is_none()) {
                guess_array = initial_guess.cast<DoubleArray>();
                if (guess_array.size() != node_count) {
                    throw std::invalid_argument("initial_guess doit avoir une valeur par noeud");
                }
                guess = guess_array.data();
            }

            apex::SolveResult res;
            {
                py::gil_scoped_release release;
                res = apex::solve_linear_system(
                    indptr.data(), indices.data(), data.data(), node_count, b.data(),
                    direct ? apex::SolverKind::Direct : apex::SolverKind::Iterative, tolerance,
                    max_iterations, guess);
            }

            py::array_t<double> values(static_cast<py::ssize_t>(node_count));
            std::copy(res.values.begin(), res.values.end(), values.mutable_data());
            py::dict info;
            info["converged"] = res.converged;
            info["iterations"] = res.iterations;
            info["error"] = res.error;
            return py::make_tuple(std::move(values), info);
        },
        py::arg("indptr"), py::arg("indices"), py::arg("data"), py::arg("b"), py::arg("direct"),
        py::arg("tolerance"), py::arg("max_iterations"), py::arg("initial_guess") = py::none(),
        "Resout A.phi = b (CSR, SPD). Rend (values, info). `converged=False` => "
        "l'appelant doit retomber sur le solveur Python.");

    m.def(
        "assemble_gradient_rhs",
        [](const IntArray& tetrahedra, const DoubleArray& shape_gradients,
           const DoubleArray& scale, const DoubleArray& directions, const Int32Array& active,
           std::int64_t node_count) {
            py::array_t<double> b(static_cast<py::ssize_t>(node_count));
            {
                py::gil_scoped_release release;
                apex::assemble_gradient_rhs(tetrahedra.data(), shape_gradients.data(),
                                            scale.data(), directions.data(), active.data(),
                                            static_cast<std::int64_t>(active.size()), node_count,
                                            b.mutable_data());
            }
            return b;
        },
        py::arg("tetrahedra"), py::arg("shape_gradients"), py::arg("scale"), py::arg("directions"),
        py::arg("active"), py::arg("node_count"),
        "Second membre du terme de direction : sum_e scale_e * G_e^T d_e.");

    py::class_<PyFemSparsity>(m, "FemSparsity")
        .def(py::init<const IntArray&, std::int64_t>(), py::arg("tetrahedra"),
             py::arg("node_count"))
        .def("indptr", &PyFemSparsity::indptr)
        .def("indices", &PyFemSparsity::indices)
        .def("assemble_gram", &PyFemSparsity::assemble_gram, py::arg("tetrahedra"),
             py::arg("shape_gradients"), py::arg("scale"), py::arg("active") = py::none(),
             "Valeurs CSR de sum_e scale_e * G_e^T G_e, dans ce motif.");

    py::class_<PyTriangleBVH>(m, "TriangleBVH")
        .def(py::init<const DoubleArray&, const DoubleArray&>(), py::arg("triangles"),
             py::arg("normals"))
        .def("intersects_location", &PyTriangleBVH::intersects_location, py::arg("ray_origins"),
             py::arg("ray_directions"),
             "Impacts dedoublonnes sur (location, index_ray). Rend le meme triplet que "
             "trimesh.ray : (locations (H,3), index_ray (H,), index_tri (H,)).")
        .def("__len__", &PyTriangleBVH::size);

    m.def(
        "triangles_cross",
        [](const DoubleArray& t0, const DoubleArray& t1) {
            if (t0.size() != 9 || t1.size() != 9) {
                throw std::invalid_argument("t0/t1 doivent contenir 3 sommets (9 doubles)");
            }
            const double* a = t0.data();
            const double* b = t1.data();
            apex::Triangle3 tri0{apex::Vec3{a[0], a[1], a[2]}, apex::Vec3{a[3], a[4], a[5]},
                                 apex::Vec3{a[6], a[7], a[8]}};
            apex::Triangle3 tri1{apex::Vec3{b[0], b[1], b[2]}, apex::Vec3{b[3], b[4], b[5]},
                                 apex::Vec3{b[6], b[7], b[8]}};
            return apex::triangles_cross(tri0, tri1);
        },
        py::arg("t0"), py::arg("t1"),
        "Test d'axe separateur (SAT) entre deux triangles (3,3) - transcription de "
        "_triangles_cross.");

    m.def(
        "interpolate_out_of_plane",
        [](const DoubleArray& query, const DoubleArray& sources, const DoubleArray& offsets,
           double power) {
            if (query.ndim() != 2 || query.shape(1) != 2) {
                throw std::invalid_argument("query doit etre de forme (K,2)");
            }
            if (sources.ndim() != 2 || sources.shape(1) != 2) {
                throw std::invalid_argument("sources doit etre de forme (M,2)");
            }
            if (offsets.ndim() != 1 || offsets.shape(0) != sources.shape(0)) {
                throw std::invalid_argument("offsets doit etre de forme (M,), aligne sur sources");
            }
            const auto query_count = static_cast<std::int64_t>(query.shape(0));
            const auto source_count = static_cast<std::int64_t>(sources.shape(0));
            py::array_t<double> out(static_cast<py::ssize_t>(query_count));
            {
                py::gil_scoped_release release;
                apex::interpolate_out_of_plane(query.data(), query_count, sources.data(),
                                               offsets.data(), source_count, power,
                                               out.mutable_data());
            }
            return out;
        },
        py::arg("query"), py::arg("sources"), py::arg("offsets"), py::arg("power"),
        "Ponderation inverse-distance 1/d^power de l'offset hors-plan (step6).");

    m.def(
        "triangle_triangle_segment",
        [](const DoubleArray& t0, const DoubleArray& t1) -> py::object {
            if (t0.size() != 9 || t1.size() != 9) {
                throw std::invalid_argument("t0/t1 doivent contenir 3 sommets (9 doubles)");
            }
            const double* pa = t0.data();
            const double* pb = t1.data();
            apex::Triangle3 tri0{apex::Vec3{pa[0], pa[1], pa[2]}, apex::Vec3{pa[3], pa[4], pa[5]},
                                 apex::Vec3{pa[6], pa[7], pa[8]}};
            apex::Triangle3 tri1{apex::Vec3{pb[0], pb[1], pb[2]}, apex::Vec3{pb[3], pb[4], pb[5]},
                                 apex::Vec3{pb[6], pb[7], pb[8]}};
            const auto result = apex::triangle_triangle_segment(tri0, tri1);
            if (!result.has_value()) {
                return py::none();
            }
            py::array_t<double> out({static_cast<py::ssize_t>(2), static_cast<py::ssize_t>(3)});
            double* o = out.mutable_data();
            o[0] = result->start.x;
            o[1] = result->start.y;
            o[2] = result->start.z;
            o[3] = result->end.x;
            o[4] = result->end.y;
            o[5] = result->end.z;
            return out;
        },
        py::arg("t0"), py::arg("t1"),
        "Segment d'intersection (2,3) de deux triangles non coplanaires, ou None.");

    m.def(
        "generalized_winding_number",
        [](const DoubleArray& points, const DoubleArray& triangles) {
            if (points.ndim() != 2 || points.shape(1) != 3) {
                throw std::invalid_argument("points doit etre de forme (P,3)");
            }
            if (triangles.ndim() != 3 || triangles.shape(1) != 3 || triangles.shape(2) != 3) {
                throw std::invalid_argument("triangles doit etre de forme (T,3,3)");
            }
            const auto point_count = static_cast<std::int64_t>(points.shape(0));
            const auto tri_count = static_cast<std::int64_t>(triangles.shape(0));
            py::array_t<double> out(static_cast<py::ssize_t>(point_count));
            {
                py::gil_scoped_release release;
                apex::generalized_winding_number(points.data(), point_count, triangles.data(),
                                                 tri_count, out.mutable_data());
            }
            return out;
        },
        py::arg("points"), py::arg("triangles"),
        "Nombre d'enroulement generalise (Van Oosterom & Strackee), bit-exact.");

    m.def(
        "find_intersection_segments",
        [](const DoubleArray& triangles, const IntArray& faces, double min_length) {
            if (triangles.ndim() != 3 || triangles.shape(1) != 3 || triangles.shape(2) != 3) {
                throw std::invalid_argument("triangles doit etre de forme (N,3,3)");
            }
            if (faces.ndim() != 2 || faces.shape(1) != 3 || faces.shape(0) != triangles.shape(0)) {
                throw std::invalid_argument("faces doit etre de forme (N,3), une par triangle");
            }
            const auto count = static_cast<std::int64_t>(triangles.shape(0));
            std::vector<apex::IntersectionSegment> segments;
            {
                py::gil_scoped_release release;
                segments = apex::find_intersection_segments(triangles.data(), faces.data(), count,
                                                             min_length);
            }
            const auto n = static_cast<py::ssize_t>(segments.size());
            py::array_t<std::int64_t> face_a(n);
            py::array_t<std::int64_t> face_b(n);
            py::array_t<double> points({n, static_cast<py::ssize_t>(2), static_cast<py::ssize_t>(3)});
            std::int64_t* fa = face_a.mutable_data();
            std::int64_t* fb = face_b.mutable_data();
            double* pts = points.mutable_data();
            for (py::ssize_t i = 0; i < n; ++i) {
                const apex::IntersectionSegment& s = segments[static_cast<std::size_t>(i)];
                fa[i] = s.face_a;
                fb[i] = s.face_b;
                pts[6 * i + 0] = s.segment.start.x;
                pts[6 * i + 1] = s.segment.start.y;
                pts[6 * i + 2] = s.segment.start.z;
                pts[6 * i + 3] = s.segment.end.x;
                pts[6 * i + 4] = s.segment.end.y;
                pts[6 * i + 5] = s.segment.end.z;
            }
            return py::make_tuple(std::move(face_a), std::move(face_b), std::move(points));
        },
        py::arg("triangles"), py::arg("faces"), py::arg("min_length"),
        "Segments d'intersection (face_a (S,), face_b (S,), points (S,2,3)) - "
        "l'appelant Python reconstruit le dict per_face.");

    m.def(
        "find_self_intersecting_faces",
        [](const DoubleArray& triangles, const IntArray& faces) {
            if (triangles.ndim() != 3 || triangles.shape(1) != 3 || triangles.shape(2) != 3) {
                throw std::invalid_argument("triangles doit etre de forme (N,3,3)");
            }
            if (faces.ndim() != 2 || faces.shape(1) != 3 || faces.shape(0) != triangles.shape(0)) {
                throw std::invalid_argument("faces doit etre de forme (N,3), une par triangle");
            }
            const auto count = static_cast<std::int64_t>(triangles.shape(0));
            std::vector<std::int32_t> bad;
            {
                py::gil_scoped_release release;
                bad = apex::find_self_intersecting_faces(triangles.data(), faces.data(), count);
            }
            py::array_t<std::int32_t> out(static_cast<py::ssize_t>(bad.size()));
            std::copy(bad.begin(), bad.end(), out.mutable_data());
            return out;
        },
        py::arg("triangles"), py::arg("faces"),
        "Indices (tries) des triangles impliques dans au moins une auto-intersection reelle. "
        "L'appelant Python les convertit en set (insensible a l'ordre).");
}
