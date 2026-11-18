// TetraNet — réseau BCC, tri dedans/dehors, ramené à la peau (« isosurface stuffing »).
// Enchaînement des phases ; le détail est dans surface.cpp (entrée), stuffing.cpp (réseau, coupes, tétraèdres)
// et checks.cpp (contrôles).
#include "checks.hpp"
#include "stuffing.hpp"

namespace tetranet {

using namespace detail;

const char* version() { return "0.0.1"; }

namespace {

// Surface d'entrée soudée et orientée, copiée dans `cleaned` ; diagnostic de la surface dans le rapport.
bool prepare_surface(const SurfaceMesh& in, const Params& prm, Surface& surf, SurfaceMesh* cleaned, Report& rep) {
    if (!(prm.cell_size > 0) || !build_surface(in, surf, prm.fix_orientation)) {
        rep.status = Status::InvalidInput;
        rep.reason = "surface vide ou taille de maille invalide";
        return false;
    }
    if (cleaned) {
        cleaned->vertices.clear(); cleaned->faces.clear();
        cleaned->vertices.reserve(surf.p.size() * 3); cleaned->faces.reserve(surf.tri.size() * 3);
        for (const V3& q : surf.p) { cleaned->vertices.push_back(q.x); cleaned->vertices.push_back(q.y); cleaned->vertices.push_back(q.z); }
        for (const auto& t : surf.tri) { cleaned->faces.push_back(t[0]); cleaned->faces.push_back(t[1]); cleaned->faces.push_back(t[2]); }
    }
    rep.surface_faces = int64_t(surf.tri.size()); rep.surface_edges = surf.edges; rep.surface_open_edges = surf.open_edges; rep.surface_holes = surf.holes;
    rep.surface_nonmanifold_edges = surf.nonmanifold_edges; rep.surface_flipped_edges = surf.flipped_edges;
    rep.surface_degenerate_faces = surf.degenerate; rep.surface_reoriented_faces = surf.reoriented; rep.surface_duplicate_faces = surf.duplicates;
    return true;
}

}  // namespace

Report mesh(const SurfaceMesh& in, const Params& prm, TetMesh& out, SurfaceMesh* cleaned) {
    auto t_all = Clock::now();
    Report rep;
    out.vertices.clear();
    out.tetrahedra.clear();
    unsigned threads = prm.threads ? prm.threads : std::max(1u, std::thread::hardware_concurrency());
    Progress prog(prm.progress);
    prog.report(Prepare, 0.0);

    Surface surf;
    if (!prepare_surface(in, prm, surf, cleaned, rep)) return rep;
    prog.report(Prepare, 0.05);               // surface soudée, orientée
    auto t0 = Clock::now();
    Bvh bvh;
    bvh.build(surf);
    rep.t_bvh = since(t0);
    prog.report(Prepare, 0.10);               // arbre de recherche construit

    Stuffing st{prm, threads, prog, rep, surf, bvh};
    st.feature_radius = prm.feature_snap_alpha * prm.cell_size;
    st.features = make_features(surf, prm, st.feature_radius);
    if (!make_lattice(st)) return rep;
    t0 = Clock::now();
    classify_vertices(st);
    snap_to_skin(st);
    rep.t_classify = since(t0);
    prog.report(Prepare, 1.0);                // préparation terminée

    prog.report(Tetrahedralize, 0.0);
    t0 = Clock::now();
    emit_tetrahedra(st, out);
    rep.t_tets = since(t0);
    prog.report(Tetrahedralize, 1.0);         // tétraèdres assemblés, sommets compactés

    check_mesh(CheckContext{prm, threads, prog, rep, surf, bvh, prm.cell_size, prm.use_winding ? &st.wn : nullptr}, out);
    rep.seconds = since(t_all);
    final_status(prm, rep);
    return rep;
}

}  // namespace tetranet
