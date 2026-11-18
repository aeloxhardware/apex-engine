// Liaison Python de TetraNet (dossier tetranet/, futur dépôt) : `apex_core.tetranet_mesh`.
// Rend la surface nettoyée (soudée, orientée), le maillage et le diagnostic chiffré.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cstdint>
#include <stdexcept>
#include <string>

#include "tetranet/tetranet.hpp"

namespace py = pybind11;

namespace {

using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
using IntArray = py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;

// Réglages acceptés dans `options` ; toute autre clé est une erreur (jamais ignorée en silence).
tetranet::Params params_from(double cell_size, const py::dict& options) {
    tetranet::Params p;
    p.cell_size = cell_size;
    p.snap_alpha = 0.5;   // défaut historique de la liaison (le réseau n'en dépend qu'avec use_cuts = faux)
    for (const auto& item : options) {
        const std::string key = py::str(item.first);
        const py::handle v = item.second;
        if (key == "snap_alpha") p.snap_alpha = v.cast<double>();
        else if (key == "skin_tolerance") p.skin_tolerance = v.cast<double>();
        else if (key == "refuse_ratio") p.refuse_conflict_ratio = v.cast<double>();
        else if (key == "threads") p.threads = v.cast<std::uint32_t>();
        else if (key == "feature_angle_deg") p.feature_angle_deg = v.cast<double>();
        else if (key == "feature_snap_alpha") p.feature_snap_alpha = v.cast<double>();
        else if (key == "check_cover") p.check_cover = v.cast<bool>();
        else if (key == "cover_shallow_alpha") p.cover_shallow_alpha = v.cast<double>();
        else throw std::invalid_argument("tetranet_mesh : option inconnue « " + key + " »");
    }
    return p;
}

py::array_t<double> points(const std::vector<double>& flat) {
    py::array_t<double> a({static_cast<py::ssize_t>(flat.size() / 3), static_cast<py::ssize_t>(3)});
    std::copy(flat.begin(), flat.end(), a.mutable_data());
    return a;
}

py::array_t<std::int64_t> indices(const std::vector<std::int32_t>& flat, py::ssize_t width) {
    py::array_t<std::int64_t> a({static_cast<py::ssize_t>(flat.size()) / width, width});
    for (size_t i = 0; i < flat.size(); ++i) a.mutable_data()[i] = flat[i];
    return a;
}

void report_into(const tetranet::Report& r, py::dict& d) {
    d["status"] = static_cast<int>(r.status); d["reason"] = r.reason; d["seconds"] = r.seconds;
    d["tets"] = r.tets; d["sign_conflicts"] = r.sign_conflicts; d["ambiguous_points"] = r.ambiguous_points;
    d["band_vertices"] = r.band_vertices; d["conflict_ratio"] = r.conflict_ratio; d["ambiguous_ratio"] = r.ambiguous_ratio;
    d["volume_ratio"] = r.volume_ratio; d["min_dihedral_deg"] = r.min_dihedral_deg; d["skin_max_error"] = r.skin_max_error;
    d["surface_faces_count"] = r.surface_faces; d["surface_open_edges"] = r.surface_open_edges; d["surface_holes"] = r.surface_holes;
    d["surface_nonmanifold_edges"] = r.surface_nonmanifold_edges; d["surface_flipped_edges"] = r.surface_flipped_edges;
    d["surface_degenerate_faces"] = r.surface_degenerate_faces; d["surface_duplicate_faces"] = r.surface_duplicate_faces;
    d["surface_reoriented_faces"] = r.surface_reoriented_faces; d["warning"] = r.warning;
    d["skin_p99_error"] = r.skin_p99_error;
    d["skin_cover_max_error"] = r.skin_cover_max_error; d["skin_cover_p99_error"] = r.skin_cover_p99_error;
    d["cover_samples"] = r.cover_samples; d["cover_internal_samples"] = r.cover_internal_samples;
    d["t_checks"] = r.t_checks;
}

py::dict tetranet_mesh(const DoubleArray& vertices, const IntArray& faces, double cell_size, const py::dict& options,
                       const py::object& progress) {
    tetranet::SurfaceMesh in;
    in.vertices.assign(vertices.data(), vertices.data() + vertices.size());
    in.faces.resize(static_cast<size_t>(faces.size()));
    for (py::ssize_t i = 0; i < faces.size(); ++i) in.faces[static_cast<size_t>(i)] = static_cast<std::int32_t>(faces.data()[i]);
    tetranet::Params p = params_from(cell_size, options);
    if (!progress.is_none()) {
        // Appelée depuis les fils de calcul (le GIL est libéré pendant `mesh`) : on le reprend le temps du rappel.
        // `progress` reste vivant tout l'appel (référence du paramètre) ; une erreur Python est ignorée : une barre
        // ne doit jamais faire échouer le maillage.
        p.progress = [&progress](int phase, double fraction) {
            py::gil_scoped_acquire acquire;
            try { progress(phase, fraction); } catch (const py::error_already_set&) { PyErr_Clear(); }
        };
    }
    tetranet::TetMesh out; tetranet::SurfaceMesh cleaned; tetranet::Report r;
    {
        py::gil_scoped_release release;
        r = tetranet::mesh(in, p, out, &cleaned);
    }
    py::dict d;
    d["tet_vertices"] = points(out.vertices); d["tetrahedra"] = indices(out.tetrahedra, 4);
    d["surface_vertices"] = points(cleaned.vertices); d["surface_faces"] = indices(cleaned.faces, 3);
    report_into(r, d);
    return d;
}

}  // namespace

void bind_tetranet(py::module_& m) {
    m.def("tetranet_mesh", &tetranet_mesh,
          py::arg("vertices"), py::arg("faces"), py::arg("cell_size"), py::arg("options") = py::dict(),
          py::arg("progress") = py::none(),
          "Maillage tetraedrique TetraNet : dict (maillage, surface nettoyee, diagnostic). `options` : snap_alpha, "
          "skin_tolerance, refuse_ratio, threads, feature_angle_deg, feature_snap_alpha, check_cover, cover_shallow_alpha "
          "(cle inconnue : ValueError). "
          "`progress(phase, fraction)` : avancement reel (phase 0 preparation, 1 tetraedralisation, 2 verification ; "
          "fraction 0..1).");
}
