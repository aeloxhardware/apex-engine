// Bindings pybind11 - ECHAFAUDAGE DE MIGRATION, pas un livrable.
//
// La cible est un moteur 100 % C++ ; ces bindings existent pour que les 305
// tests pytest servent d'ORACLE D'EQUIVALENCE pendant tout le portage (le
// Python et le C++ coexistent dans le meme processus, donc un test peut
// appeler les deux et les comparer). Ils seront retires a la bascule finale.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

#include "apex/bvh.hpp"
#include "apex/field_solver.hpp"
#include "apex/gcode_export.hpp"
#include "apex/gpu/vk_context.hpp"
#include "apex/geometry2d.hpp"
#include "apex/heat_method.hpp"
#include "apex/kdtree.hpp"
#include "apex/kinematics.hpp"
#include "apex/interpolate.hpp"
#include "apex/isosurface.hpp"
#include "apex/mesh_distance.hpp"
#include "apex/mesh_intersect.hpp"
#include "apex/mesh_repair.hpp"
#include "apex/overhang.hpp"
#include "apex/params.hpp"
#include "apex/past_material.hpp"
#include "apex/segmentation.hpp"
#include "apex/printed_index.hpp"
#include "apex/reorientation.hpp"
#include "apex/scheduling.hpp"
#include "apex/support.hpp"
#include "apex/tet_topology.hpp"
#include "apex/toolpath2d.hpp"
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
        "reorient_coupled_tiers",
        [](const DoubleArray& positions, const DoubleArray& extrusion, const DoubleArray& angles,
           const DoubleArray& axes, std::int64_t free_axis, const DoubleArray& grid,
           const DoubleArray& tiers, const DoubleArray& samples,
           const apex::PastMaterialIndex& material, double reach, const py::object& plate_origin,
           double plate_radius, double plate_half_height, double plate_center_offset,
           double safety_margin, double lo, double hi, double full_mm, double taper_mm,
           const IntArray& contacts) {
            apex::ReorientationInput in;
            in.n = static_cast<std::int64_t>(positions.shape(0));
            in.axis_count = static_cast<std::int64_t>(axes.shape(0));
            in.positions = positions.data();
            in.extrusion = extrusion.data();
            in.angles = angles.data();
            in.axes = axes.data();
            in.free_axis = free_axis;
            in.grid = grid.data();
            in.grid_count = static_cast<std::int64_t>(grid.size());
            in.tiers = tiers.data();
            in.tier_count = static_cast<std::int64_t>(tiers.size());
            in.samples = samples.data();
            in.sample_count = static_cast<std::int64_t>(samples.shape(0));
            in.reach = reach;
            if (!plate_origin.is_none()) {
                DoubleArray o = plate_origin.cast<DoubleArray>();
                in.has_plate = true;
                for (int i = 0; i < 3; ++i) in.plate_origin[i] = o.data()[i];
            }
            in.plate_radius = plate_radius;
            in.plate_half_height = plate_half_height;
            in.plate_center_offset = plate_center_offset;
            in.safety_margin = safety_margin;
            in.lo = lo; in.hi = hi; in.full_mm = full_mm; in.taper_mm = taper_mm;
            py::array_t<double> out({in.n, in.axis_count});
            std::vector<apex::ReorientationEntry> log;
            {
                py::gil_scoped_release release;
                log = apex::reorient_coupled_tiers(in, material, contacts.data(),
                                                   static_cast<std::int64_t>(contacts.size()),
                                                   out.mutable_data());
            }
            py::list entries;
            for (const auto& e : log) {
                entries.append(py::make_tuple(e.index, e.found, e.angle, e.tier));
            }
            return py::make_tuple(std::move(out), std::move(entries));
        },
        "Reorientation par paliers (etape 10) : (angles (n,A), [(indice, trouve, angle, palier)]).");

    m.def(
        "try_local_reorientation_coupled",
        [](const DoubleArray& positions, const DoubleArray& feedrates, const DoubleArray& angles,
           const DoubleArray& axes, std::int64_t point_index, const DoubleArray& inner_lo,
           const DoubleArray& inner_hi, std::int64_t samples_per_axis, double max_deviation,
           double max_velocity_deg_s, std::int64_t neighbor_margin, const DoubleArray& samples,
           const apex::PastMaterialIndex& material, double reach, const py::object& plate_origin,
           double plate_radius, double plate_half_height, double plate_center_offset,
           double safety_margin) {
            apex::LocalReorientationInput in;
            in.n = static_cast<std::int64_t>(positions.shape(0));
            in.axis_count = static_cast<std::int64_t>(axes.shape(0));
            in.point_index = point_index;
            in.positions = positions.data(); in.feedrates = feedrates.data();
            in.angles = angles.data(); in.axes = axes.data();
            in.inner_lo = inner_lo.data(); in.inner_hi = inner_hi.data();
            in.samples_per_axis = samples_per_axis;
            in.max_direction_deviation_deg = max_deviation;
            in.max_velocity_deg_s = max_velocity_deg_s;
            in.neighbor_margin = neighbor_margin;
            in.samples = samples.data();
            in.sample_count = static_cast<std::int64_t>(samples.shape(0));
            in.reach = reach;
            if (!plate_origin.is_none()) {
                DoubleArray o = plate_origin.cast<DoubleArray>();
                in.has_plate = true;
                for (int i = 0; i < 3; ++i) in.plate_origin[i] = o.data()[i];
            }
            in.plate_radius = plate_radius; in.plate_half_height = plate_half_height;
            in.plate_center_offset = plate_center_offset; in.safety_margin = safety_margin;
            py::array_t<double> out(in.axis_count);
            bool ok;
            {
                py::gil_scoped_release release;
                ok = apex::try_local_reorientation_coupled(in, material, out.mutable_data());
            }
            if (!ok) return py::object(py::none());
            return py::object(std::move(out));
        },
        "Reorientation locale (machine a premier axe couple) : angles retenus ou None.");

    m.def(
        "clamp_head_height_apply",
        [](const DoubleArray& z, const DoubleArray& current, const DoubleArray& grid,
           const DoubleArray& zlow, double plate_z, double safety_margin) {
            const auto n = static_cast<py::ssize_t>(z.size());
            py::array_t<double> out(n);
            py::array_t<std::uint8_t> moved(n);
            apex::clamp_head_height_apply(z.data(), current.data(), n, grid.data(), zlow.data(),
                                          static_cast<std::int64_t>(grid.size()), plate_z,
                                          safety_margin, out.mutable_data(), moved.mutable_data());
            return py::make_tuple(std::move(out), std::move(moved));
        },
        "Borne de B selon la hauteur de la pointe : (angles, masque des points ramenes).");

    m.def(
        "heat_method_distance",
        [](const DoubleArray& tet_vertices, const IntArray& tetrahedra,
           const DoubleArray& shape_gradients, const DoubleArray& volumes,
           const IntArray& source_nodes, double time_factor, double anchor_weight,
           double tolerance, int max_iterations, std::int64_t direct_threshold,
           const py::object& origin, const DoubleArray& up) {
            const auto tet_count = static_cast<std::int64_t>(tetrahedra.shape(0));
            const auto node_count = static_cast<std::int64_t>(tet_vertices.shape(0));
            apex::HeatParams p;
            p.time_factor = time_factor;
            p.anchor_weight = anchor_weight;
            p.tolerance = tolerance;
            p.max_iterations = max_iterations;
            p.direct_threshold = direct_threshold;
            if (!origin.is_none()) {
                DoubleArray o = origin.cast<DoubleArray>();
                p.project_reachable = true;
                for (int i = 0; i < 3; ++i) p.origin[i] = o.data()[i];
                for (int i = 0; i < 3; ++i) p.up[i] = up.data()[i];
            }
            apex::HeatResult res;
            {
                py::gil_scoped_release release;
                res = apex::heat_method_distance(
                    tet_vertices.data(), tetrahedra.data(), tet_count, node_count,
                    shape_gradients.data(), volumes.data(), source_nodes.data(),
                    static_cast<std::int64_t>(source_nodes.size()), p);
            }
            if (!res.converged) return py::object(py::none());
            py::array_t<double> phi(static_cast<py::ssize_t>(res.phi.size()));
            std::copy(res.phi.begin(), res.phi.end(), phi.mutable_data());
            return py::object(std::move(phi));
        },
        py::arg("tet_vertices"), py::arg("tetrahedra"), py::arg("shape_gradients"),
        py::arg("volumes"), py::arg("source_nodes"), py::arg("time_factor"),
        py::arg("anchor_weight"), py::arg("tolerance"), py::arg("max_iterations"),
        py::arg("direct_threshold"), py::arg("origin"), py::arg("up"),
        "Methode de la chaleur : distance au plateau par noeud, ou None si un solveur "
        "n'a pas converge (l'appelant retombe sur le Python).");

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

    // Pas 11 : sonde de l'en-tete genere apex/params.hpp (scripts/gen_params.py).
    // Prouve qu'il compile et que ses valeurs par defaut sont celles des
    // dataclasses Python (tests/test_config_schema.py compare).
    m.def("params_defaults_probe", []() {
        apex::params::PipelineParams p;
        apex::params::GcodeParams g;
        py::dict d;
        d["resolution"] = p.resolution;
        d["near_platform_resolution"] = p.near_platform_resolution.value_or(-1.0);
        d["layer_height_base"] = p.layers.layer_height_base;
        d["out_of_plane_interpolation_power"] = p.toolpath2d.out_of_plane_interpolation_power.value_or(-1.0);
        d["split_rings_into_runs"] = p.toolpath2d.split_rings_into_runs;
        d["infill_pattern"] = std::string(apex::params::to_string(p.toolpath2d.infill_pattern));
        d["direction_blend_radius"] = p.field_solver.direction_blend_radius;
        d["enforce_upward_growth"] = p.field_solver.enforce_upward_growth;
        d["rotation_axes_count"] = static_cast<int>(p.machine_profile.rotation_axes.size());
        d["b_max_deg"] = p.machine_profile.rotation_axes.size() > 1 ? p.machine_profile.rotation_axes[1].angle_range_deg[1] : -1.0;
        d["position_coupled_axes"] = p.machine_profile.position_coupled_axes;
        d["pivot_to_tip"] = g.pivot_to_tip;
        d["tip_at_bead_top"] = g.tip_at_bead_top;
        d["z_reference"] = g.z_reference;
        return d;
    });

    // Pas 13 : marching tetrahedra + soudure (step5), equivalence au bit.
    py::class_<apex::IsosurfaceExtractor>(m, "IsosurfaceExtractor")
        .def(py::init([](const DoubleArray& tet_vertices, const IntArray& tetrahedra, const DoubleArray& values) {
                 if (tet_vertices.ndim() != 2 || tet_vertices.shape(1) != 3) {
                     throw std::invalid_argument("tet_vertices doit etre de forme (N,3)");
                 }
                 if (tetrahedra.ndim() != 2 || tetrahedra.shape(1) != 4) {
                     throw std::invalid_argument("tetrahedra doit etre de forme (T,4)");
                 }
                 if (values.ndim() != 1 || values.shape(0) != tet_vertices.shape(0)) {
                     throw std::invalid_argument("values doit etre de forme (N,), aligne sur tet_vertices");
                 }
                 const auto n = static_cast<std::int64_t>(tet_vertices.shape(0));
                 const std::int64_t* q = tetrahedra.data();
                 const auto t = static_cast<std::int64_t>(tetrahedra.shape(0));
                 for (std::int64_t i = 0; i < t * 4; ++i) {
                     if (q[i] < 0 || q[i] >= n) {
                         throw std::invalid_argument("indice de tetraedre hors bornes");
                     }
                 }
                 return new apex::IsosurfaceExtractor(tet_vertices.data(), n, q, t, values.data());
             }),
             py::arg("tet_vertices"), py::arg("tetrahedra"), py::arg("values"))
        .def("extract",
             [](const apex::IsosurfaceExtractor& self, double isovalue, double weld_tolerance) {
                 apex::IsosurfaceResult r;
                 {
                     py::gil_scoped_release release;
                     r = self.extract(isovalue, weld_tolerance);
                 }
                 const auto nv = static_cast<py::ssize_t>(r.vertices.size() / 3);
                 const auto nf = static_cast<py::ssize_t>(r.source_tets.size());
                 py::array_t<double> verts({nv, static_cast<py::ssize_t>(3)});
                 std::copy(r.vertices.begin(), r.vertices.end(), verts.mutable_data());
                 py::array_t<std::int64_t> faces({nf, static_cast<py::ssize_t>(3)});
                 std::copy(r.faces.begin(), r.faces.end(), faces.mutable_data());
                 py::array_t<std::int64_t> src(nf);
                 std::copy(r.source_tets.begin(), r.source_tets.end(), src.mutable_data());
                 return py::make_tuple(verts, faces, src);
             },
             py::arg("isovalue"), py::arg("weld_tolerance"),
             "(vertices soudes (V,3), faces (F,3), tetraedre source (F,)) - "
             "_extract_isosurface_triangles puis _weld_vertices de step5.")
        .def_property_readonly("tet_count", &apex::IsosurfaceExtractor::tet_count);

    // Pas 14 : trace de bord (outline) natif, None = repli trimesh.
    m.def(
        "boundary_loops",
        [](const IntArray& faces, std::int64_t vertex_count) -> py::object {
            if (faces.ndim() != 2 || faces.shape(1) != 3) {
                throw std::invalid_argument("faces doit etre de forme (F,3)");
            }
            const std::int64_t* f = faces.data();
            const auto n = static_cast<std::int64_t>(faces.shape(0));
            for (std::int64_t i = 0; i < n * 3; ++i) {
                if (f[i] < 0 || f[i] >= vertex_count) {
                    throw std::invalid_argument("indice de face hors bornes");
                }
            }
            apex::BoundaryLoops r;
            {
                py::gil_scoped_release release;
                r = apex::boundary_loops(f, n, vertex_count);
            }
            if (!r.ok) {
                return py::none();
            }
            py::list loops;
            for (const auto& loop : r.loops) {
                py::array_t<std::int64_t> a(static_cast<py::ssize_t>(loop.size()));
                std::copy(loop.begin(), loop.end(), a.mutable_data());
                loops.append(a);
            }
            return loops;
        },
        py::arg("faces"), py::arg("vertex_count"),
        "Boucles de bord (indices, fermees) dans l'ordre de trimesh.outline(), "
        "ou None si le bord n'est pas fait de boucles simples (repli trimesh).");

    // Pas 12 : lancer de rayon 2D (core/geometry_utils.point_in_polygon_2d), au bit.
    m.def(
        "point_in_polygon_2d",
        [](double x, double y, const DoubleArray& polygon) {
            if (polygon.ndim() != 2 || polygon.shape(1) != 2) {
                throw std::invalid_argument("polygon doit etre de forme (N,2)");
            }
            return apex::point_in_polygon_2d(x, y, polygon.data(), static_cast<std::int64_t>(polygon.shape(0)));
        },
        py::arg("x"), py::arg("y"), py::arg("polygon"));

    // Pas 19 : matiere deja deposee (step9 rapide), minimum exact au bit.
    py::class_<apex::PastMaterialIndex>(m, "PastMaterialIndex")
        .def(py::init([](const DoubleArray& points, const IntArray& owners) {
                 if (points.ndim() != 2 || points.shape(1) != 3) {
                     throw std::invalid_argument("points doit etre de forme (N,3)");
                 }
                 if (owners.ndim() != 1 || owners.shape(0) != points.shape(0)) {
                     throw std::invalid_argument("owners doit etre de forme (N,)");
                 }
                 return new apex::PastMaterialIndex(points.data(), owners.data(),
                                                    static_cast<std::int64_t>(points.shape(0)));
             }),
             py::arg("points"), py::arg("owners"))
        .def("min_past_distance",
             [](const apex::PastMaterialIndex& self, const DoubleArray& queries, std::int64_t limit,
                double reach) {
                 if (queries.ndim() != 2 || queries.shape(1) != 3) {
                     throw std::invalid_argument("queries doit etre de forme (Q,3)");
                 }
                 py::gil_scoped_release release;
                 return self.min_past_distance(queries.data(), static_cast<std::int64_t>(queries.shape(0)),
                                               limit, reach);
             },
             py::arg("queries"), py::arg("limit"), py::arg("reach"))
        .def("min_past_distance_batch",
             [](const apex::PastMaterialIndex& self, const DoubleArray& queries, const IntArray& limits,
                double reach, int threads) {
                 if (queries.ndim() != 3 || queries.shape(2) != 3 || limits.shape(0) != queries.shape(0)) {
                     throw std::invalid_argument("queries (P,S,3) et limits (P,) attendus");
                 }
                 const auto P = static_cast<std::int64_t>(queries.shape(0));
                 py::array_t<double> out(P);
                 {
                     py::gil_scoped_release release;
                     self.min_past_distance_batch(queries.data(), P, static_cast<std::int64_t>(queries.shape(1)),
                                                  limits.data(), reach, out.mutable_data(), threads);
                 }
                 return out;
             },
             py::arg("queries"), py::arg("limits"), py::arg("reach"), py::arg("threads") = 0)
        .def_property_readonly("size", &apex::PastMaterialIndex::size);

    // Pas 24 : corps de l'export G-code, texte identique ligne a ligne.
    m.def(
        "gcode_header",
        [](double pivot_to_tip, double retract_length, double retract_speed, double nozzle_temp) {
            std::string out;
            if (!apex::gcode_header(pivot_to_tip, retract_length, retract_speed, nozzle_temp, out)) {
                return py::object(py::none());
            }
            return py::object(py::str(out));
        },
        py::arg("pivot_to_tip"), py::arg("retract_length"), py::arg("retract_speed"),
        py::arg("nozzle_temp"), "En-tete du G-code (None si le format Python diverge).");
    m.def("gcode_footer", []() { return apex::gcode_footer(); }, "Pied du G-code.");

    m.def(
        "export_gcode_body",
        [](const DoubleArray& positions, const DoubleArray& angles, const DoubleArray& extrusion,
           const py::array_t<std::uint8_t, py::array::c_style | py::array::forcecast>& skin,
           const DoubleArray& slow, const DoubleArray& zlow_grid, const DoubleArray& zlow_values,
           const py::dict& p, const py::object& travel) {
            const auto n = static_cast<std::int64_t>(positions.shape(0));
            if (positions.ndim() != 2 || positions.shape(1) != 3 || angles.ndim() != 2 || angles.shape(1) != 2 ||
                angles.shape(0) != n || extrusion.shape(0) != n || skin.shape(0) != n ||
                (n > 0 && slow.shape(0) != n - 1) || zlow_grid.shape(0) != zlow_values.shape(0)) {
                throw std::invalid_argument("formes incoherentes pour export_gcode_body");
            }
            apex::GcodeExportParams P{};
            P.layer_height = p["layer_height"].cast<double>();
            P.skin_thickness = p["skin_thickness"].cast<double>();
            P.line_width = p["line_width"].cast<double>();
            P.flow = p["flow"].cast<double>();
            P.print_speed = p["print_speed"].cast<double>();
            P.travel_speed = p["travel_speed"].cast<double>();
            P.rotation_feed = p["rotation_feed"].cast<double>();
            P.max_feed = p["max_feed"].cast<double>();
            P.retract_length = p["retract_length"].cast<double>();
            P.retract_speed = p["retract_speed"].cast<double>();
            P.retract_min_travel = p["retract_min_travel"].cast<double>();
            P.travel_clearance = p["travel_clearance"].cast<double>();
            P.pivot_to_tip = p["pivot_to_tip"].cast<double>();
            P.platform_margin = p["platform_margin"].cast<double>();
            P.singularity_radius = p["singularity_radius"].cast<double>();
            P.tilt_sign = p["tilt_sign"].cast<double>();
            P.b_sign = p["b_sign"].cast<double>();
            P.z_reference_tip = p["z_reference"].cast<std::string>() == "tip";
            P.decimals = p["decimals"].cast<int>();
            P.fil_area = p["fil_area"].cast<double>();
            if (p.contains("max_tip_deviation") && !p["max_tip_deviation"].is_none())
                P.max_tip_deviation = p["max_tip_deviation"].cast<double>();
            std::vector<double> zg(zlow_grid.data(), zlow_grid.data() + zlow_grid.shape(0));
            std::vector<double> zv(zlow_values.data(), zlow_values.data() + zlow_values.shape(0));
            // Controle des degagements (material_travel_check) : dict de DONNEES
            // {points (M,3), owners (M,), axis_t (K,), radius, axis_c (3,), axis_b (3,), max_factor}.
            std::unique_ptr<apex::PastMaterialIndex> material;
            apex::TravelCheck tc;
            const apex::TravelCheck* tcp = nullptr;
            if (!travel.is_none()) {
                const py::dict t = travel.cast<py::dict>();
                const auto pts = t["points"].cast<DoubleArray>();
                const auto own = t["owners"].cast<py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>>();
                const auto at = t["axis_t"].cast<DoubleArray>();
                const auto ac = t["axis_c"].cast<DoubleArray>();
                const auto ab = t["axis_b"].cast<DoubleArray>();
                if (pts.ndim() != 2 || pts.shape(1) != 3 || own.shape(0) != pts.shape(0) || ac.shape(0) != 3 ||
                    ab.shape(0) != 3) {
                    throw std::invalid_argument("formes incoherentes pour le controle des degagements");
                }
                material = std::make_unique<apex::PastMaterialIndex>(pts.data(), own.data(),
                                                                      static_cast<std::int64_t>(pts.shape(0)));
                tc.material = material.get();
                tc.axis_t.assign(at.data(), at.data() + at.shape(0));
                tc.radius = t["radius"].cast<double>();
                tc.max_factor = t["max_factor"].cast<double>();
                for (int k = 0; k < 3; ++k) { tc.axis_c[k] = ac.data()[k]; tc.axis_b[k] = ab.data()[k]; }
                tcp = &tc;
            }
            apex::GcodeExportStats st;
            std::string body;
            {
                py::gil_scoped_release release;
                body = apex::export_gcode_body(positions.data(), angles.data(), extrusion.data(), skin.data(),
                                               slow.data(), n, zg, zv, P, st, tcp);
            }
            py::dict s;
            s["moves"] = st.moves; s["travels"] = st.travels; s["retractions"] = st.retractions;
            s["subdivided_moves"] = st.subdivided_moves; s["max_tip_deviation_before"] = st.max_tip_deviation_before;
            s["clearance_Z"] = st.clearance_z; s["clearance_R"] = st.clearance_r; s["clearance_none"] = st.clearance_none;
            s["filament_mm"] = st.filament_mm; s["c_mismatch_deg"] = st.c_mismatch_deg;
            s["tip_z_min"] = st.tip_z_min; s["head_z_min"] = st.head_z_min;
            s["tip_below_plate"] = st.tip_below_plate; s["head_below_margin"] = st.head_below_margin;
            s["any"] = st.any;
            s["x"] = py::make_tuple(st.x_min, st.x_max); s["y"] = py::make_tuple(st.y_min, st.y_max);
            s["z"] = py::make_tuple(st.z_min, st.z_max);
            s["travel_collisions"] = st.travel_collisions;
            return py::make_tuple(py::str(body), s);
        },
        py::arg("positions"), py::arg("angles"), py::arg("extrusion"), py::arg("skin"), py::arg("slow"),
        py::arg("zlow_grid"), py::arg("zlow_values"), py::arg("params"), py::arg("travel") = py::none());
    // Capacite : export_gcode_body accepte le controle des degagements (2026-10-13).
    m.attr("EXPORT_TRAVEL_CHECK") = true;

    m.def(
        "feed_plan_caps",
        [](const DoubleArray& d, double vmax, double amax) {
            py::array_t<double> cap(d.shape(0));
            apex::feed_plan_caps(d.data(), static_cast<std::int64_t>(d.shape(0)), vmax, amax, cap.mutable_data());
            return cap;
        },
        py::arg("d"), py::arg("max_velocity_deg"), py::arg("max_accel_deg"));

    // Pas 15 : remplissage 2D d'une zone (step6), parois par Clipper2.
    m.def(
        "fill_region",
        [](const DoubleArray& outer, const std::vector<DoubleArray>& holes, double nozzle_width, int wall_count,
           double degenerate_area_ratio, double infill_density, int pattern, double sweep_angle_rad) {
            auto to_ring = [](const DoubleArray& a) {
                if (a.ndim() != 2 || a.shape(1) < 2) throw std::invalid_argument("anneau de forme (K,2) attendu");
                apex::Ring2 r(static_cast<std::size_t>(a.shape(0)));
                const auto cols = a.shape(1);
                for (py::ssize_t i = 0; i < a.shape(0); ++i) r[i] = {a.data()[i * cols], a.data()[i * cols + 1]};
                return r;
            };
            const apex::Ring2 o = to_ring(outer);
            std::vector<apex::Ring2> hs;
            for (const auto& h : holes) hs.push_back(to_ring(h));
            apex::FillRegionParams P;
            P.nozzle_width = nozzle_width;
            P.wall_count = wall_count;
            P.degenerate_area_ratio = degenerate_area_ratio;
            P.infill_density = infill_density;
            P.pattern = static_cast<apex::InfillKind>(pattern);
            P.sweep_angle_rad = sweep_angle_rad;
            apex::FillRegionResult res;
            {
                py::gil_scoped_release release;
                res = apex::fill_region(o, hs, P);
            }
            auto to_list = [](const std::vector<apex::Ring2>& rings) {
                py::list out;
                for (const auto& r : rings) {
                    py::array_t<double> a({static_cast<py::ssize_t>(r.size()), static_cast<py::ssize_t>(2)});
                    double* d = a.mutable_data();
                    for (std::size_t i = 0; i < r.size(); ++i) { d[2 * i] = r[i][0]; d[2 * i + 1] = r[i][1]; }
                    out.append(a);
                }
                return out;
            };
            return py::make_tuple(res.base_valid, to_list(res.infill), to_list(res.walls));
        },
        py::arg("outer"), py::arg("holes"), py::arg("nozzle_width"), py::arg("wall_count"),
        py::arg("degenerate_area_ratio"), py::arg("infill_density"), py::arg("pattern"),
        py::arg("sweep_angle_rad"));

    // Pas 21a : etape 2 entiere (surplombs, risque, ponts, portee).
    m.def(
        "overhang_map",
        [](const DoubleArray& vertices, const IntArray& faces, const DoubleArray& normals,
           double threshold_deg, int ray_density, double support_ratio, double max_gap,
           double horizontal_tol_deg, const DoubleArray& up) {
            apex::OverhangParams P;
            P.overhang_angle_threshold_deg = threshold_deg;
            P.bridge_detection_ray_density = ray_density;
            P.bridge_support_ratio_threshold = support_ratio;
            P.bridge_max_gap_distance = max_gap;
            P.bridge_horizontal_angle_tolerance_deg = horizontal_tol_deg;
            for (int k = 0; k < 3; ++k) P.up[k] = up.data()[k];
            const auto F = static_cast<std::int64_t>(faces.shape(0));
            apex::OverhangResult r;
            {
                py::gil_scoped_release release;
                r = apex::compute_overhang_map(vertices.data(), static_cast<std::int64_t>(vertices.shape(0)),
                                               faces.data(), normals.data(), F, P);
            }
            py::array_t<double> risk(F), span({F, static_cast<std::int64_t>(3)});
            py::array_t<bool> over(F), bridge(F);
            std::copy(r.risk_score.begin(), r.risk_score.end(), risk.mutable_data());
            std::copy(r.span_direction.begin(), r.span_direction.end(), span.mutable_data());
            for (std::int64_t f = 0; f < F; ++f) {
                over.mutable_data()[f] = r.is_overhang[f] != 0;
                bridge.mutable_data()[f] = r.is_bridge[f] != 0;
            }
            return py::make_tuple(risk, over, bridge, span);
        },
        py::arg("vertices"), py::arg("faces"), py::arg("normals"), py::arg("threshold_deg"),
        py::arg("ray_density"), py::arg("support_ratio"), py::arg("max_gap"), py::arg("horizontal_tol_deg"),
        py::arg("up"));

    // Pas 21b : etape 3, chemin par defaut.
    m.def(
        "segment_default",
        [](const DoubleArray& vertices, const IntArray& faces, const DoubleArray& normals, const DoubleArray& risk,
           const py::array_t<std::uint8_t, py::array::c_style | py::array::forcecast>& platform,
           const DoubleArray& directions, double max_angle_deg, double merge_angle_deg, const DoubleArray& up) {
            apex::SegmentationInput in{};
            in.vertices = vertices.data();
            in.faces = faces.data();
            in.face_normals = normals.data();
            in.risk = risk.data();
            in.platform = platform.data();
            in.face_count = static_cast<std::int64_t>(faces.shape(0));
            in.directions = directions.data();
            in.direction_count = static_cast<std::int64_t>(directions.shape(0));
            in.max_angle_deg = max_angle_deg;
            in.merge_angle_deg = merge_angle_deg;
            for (int k = 0; k < 3; ++k) in.up[k] = up.data()[k];
            apex::SegmentationResult r;
            {
                py::gil_scoped_release release;
                r = apex::segment_default(in);
            }
            auto vec = [](const std::vector<std::int64_t>& v) {
                py::array_t<std::int64_t> a(static_cast<py::ssize_t>(v.size()));
                std::copy(v.begin(), v.end(), a.mutable_data());
                return a;
            };
            return py::make_tuple(vec(r.labels), vec(r.direction_index), vec(r.depth), vec(r.dependencies));
        },
        py::arg("vertices"), py::arg("faces"), py::arg("normals"), py::arg("risk"), py::arg("platform"),
        py::arg("directions"), py::arg("max_angle_deg"), py::arg("merge_angle_deg"), py::arg("up"));

    // Pas 21c : etape 7 (couleur, trous temporaires, ordonnancement).
    // `layers` : liste de tuples (id, frame(9), [polylignes (K,C)], C, [regions], isovaleur ou None).
    // `face_color` : int64 par face, `painted` : bool par face (False = non peinte).
    m.def(
        "schedule_layers",
        [](const py::list& layers, const DoubleArray& vertices, const IntArray& faces, const DoubleArray& normals,
           const IntArray& face_color, const py::array_t<std::uint8_t, py::array::c_style | py::array::forcecast>& painted,
           double merge_tolerance, const DoubleArray& tet_vertices, const DoubleArray& field_values,
           const DoubleArray& escape_directions, double surface_color_distance, double reheat, double adhesion,
           double volume, double time_per, std::int64_t max_open_holes, std::int64_t ray_samples) {
            std::vector<apex::ScheduleLayer> L;
            for (const auto& item : layers) {
                const auto t = item.cast<py::tuple>();
                apex::ScheduleLayer layer;
                layer.id = t[0].cast<std::int64_t>();
                const auto frame = t[1].cast<DoubleArray>();
                std::copy(frame.data(), frame.data() + 9, layer.frame);
                layer.columns = t[3].cast<std::int64_t>();
                for (const auto& poly : t[2].cast<py::list>()) {
                    const auto a = poly.cast<DoubleArray>();
                    layer.toolpath.emplace_back(a.data(), a.data() + a.size());
                }
                layer.region_ids = t[4].cast<std::vector<std::int64_t>>();
                if (!t[5].is_none()) layer.isovalue = t[5].cast<double>();
                L.push_back(std::move(layer));
            }
            apex::ScheduleInput in;
            in.vertices = vertices.data(); in.vertex_count = static_cast<std::int64_t>(vertices.shape(0));
            in.faces = faces.data(); in.face_normals = normals.data(); in.face_count = static_cast<std::int64_t>(faces.shape(0));
            if (face_color.size() > 0) {
                in.face_color.resize(face_color.size());
                for (py::ssize_t f = 0; f < face_color.size(); ++f)
                    if (painted.data()[f]) in.face_color[f] = face_color.data()[f];
            }
            in.merge_tolerance = merge_tolerance;
            in.tet_vertices = tet_vertices.data(); in.tet_vertex_count = static_cast<std::int64_t>(tet_vertices.shape(0));
            in.field_values = field_values.data();
            in.escape_directions.assign(escape_directions.data(), escape_directions.data() + escape_directions.size());
            in.params.surface_color_distance = surface_color_distance;
            in.params.reheat_penalty = reheat; in.params.adhesion_risk_weight = adhesion;
            in.params.volume_per_purge = volume; in.params.time_per_purge = time_per;
            in.params.max_open_holes = max_open_holes; in.params.deadline_ray_samples = ray_samples;
            apex::ScheduleResult r;
            {
                py::gil_scoped_release release;
                r = apex::schedule_layers(L, in);
            }
            auto color = [](const apex::Color& c) -> py::object { return c ? py::object(py::int_(*c)) : py::none(); };
            py::list out_layers;
            for (std::size_t li = 0; li < L.size(); ++li) {
                py::list polys, colors;
                const auto C = static_cast<py::ssize_t>(L[li].columns);
                for (const auto& poly : r.toolpath[li]) {
                    py::array_t<double> a({static_cast<py::ssize_t>(poly.size()) / C, C});
                    std::copy(poly.begin(), poly.end(), a.mutable_data());
                    polys.append(a);
                }
                for (const auto& c : r.segment_colors[li]) colors.append(color(c));
                out_layers.append(py::make_tuple(polys, colors, color(r.dominant_color[li])));
            }
            py::list passes;
            for (const auto& p : r.passes) {
                py::list segs;
                for (const auto& s : p.segments) segs.append(py::make_tuple(s.first, s.second));
                passes.append(py::make_tuple(segs, color(p.color), p.is_hole_fill,
                                             py::make_tuple(p.layer_min, p.layer_max), p.deferred_regions));
            }
            return py::make_tuple(out_layers, passes);
        });

    // Socle Vulkan (phase 1 du plan Bird/Vulkan, 2026-10-13). Jamais d'exception
    // pour l'absence de GPU : gpu_status() rend la raison.
    m.def(
        "gpu_status",
        []() {
            const apex::gpu::GpuStatus s = apex::gpu::status();
            py::dict d;
            d["ok"] = s.ok; d["reason"] = s.reason; d["device"] = s.device; d["fp64"] = s.fp64;
            d["timestamp_period_ns"] = s.timestamp_period_ns;
            d["compiled"] = apex::gpu::compiled_with_vulkan();
            return d;
        },
        "Etat du backend Vulkan : {ok, reason, device, fp64, timestamp_period_ns, compiled}.");
    m.def(
        "gpu_saxpy",
        [](double a, const py::array& x, const py::array& y) -> py::tuple {
            // y <- a*x + y sur le GPU (fumee) ; float32 ou float64 selon le type de x. Rend (y, temps GPU s).
            if (py::isinstance<py::array_t<float>>(x)) {
                auto xf = py::array_t<float, py::array::c_style | py::array::forcecast>::ensure(x);
                auto yf = py::array_t<float, py::array::c_style | py::array::forcecast>::ensure(y);
                py::array_t<float> out(yf.size());
                std::copy(yf.data(), yf.data() + yf.size(), out.mutable_data());
                if (xf.size() != out.size()) throw std::invalid_argument("x et y de tailles differentes");
                double t;
                {
                    py::gil_scoped_release release;
                    t = apex::gpu::saxpy_f32(static_cast<float>(a), xf.data(), out.mutable_data(), xf.size());
                }
                return py::make_tuple(out, t);
            }
            auto xd = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(x);
            auto yd = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(y);
            py::array_t<double> out(yd.size());
            std::copy(yd.data(), yd.data() + yd.size(), out.mutable_data());
            if (xd.size() != out.size()) throw std::invalid_argument("x et y de tailles differentes");
            double t;
            {
                py::gil_scoped_release release;
                t = apex::gpu::saxpy_f64(a, xd.data(), out.mutable_data(), xd.size());
            }
            return py::make_tuple(out, t);
        },
        py::arg("a"), py::arg("x"), py::arg("y"));

}
