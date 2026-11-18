// TetraNet — contrôles du maillage produit.
#include "checks.hpp"
#include "cover.hpp"

namespace tetranet::detail {

namespace {

V3 vertex(const TetMesh& m, int32_t i) { return V3(m.vertices[3 * size_t(i)], m.vertices[3 * size_t(i) + 1], m.vertices[3 * size_t(i) + 2]); }

struct FK { int32_t a, b, c; bool operator==(const FK& o) const { return a == o.a && b == o.b && c == o.c; } };
struct FH { size_t operator()(const FK& k) const { return size_t((uint64_t(uint32_t(k.a)) * 0x9E3779B97F4A7C15ull) ^ (uint64_t(uint32_t(k.b)) << 21) ^ (uint64_t(uint32_t(k.c)) * 0xC2B2AE3D27D4EB4Full)); } };
using FaceCount = std::unordered_map<FK, int32_t, FH>;   // face triée -> nombre de tétraèdres

// Volume, volume minimal, angle dièdre minimal ; compte les faces de chaque tétraèdre.
void measure_tets(const CheckContext& c, const TetMesh& out, FaceCount& faces) {
    Report& rep = c.rep;
    Span verify_tets(c.prog, Verify, 0.0, 0.70, std::max<int64_t>(1, rep.tets));   // tétraèdres contrôlés
    double vol = 0, minvol = 1e300, mind = 180.0;
    {
    auto verify_tick = verify_tets.tick();
    for (int64_t t = 0; t < rep.tets; ++t) {
        verify_tick();
        const int32_t* q = &out.tetrahedra[4 * size_t(t)];
        V3 p[4] = {vertex(out, q[0]), vertex(out, q[1]), vertex(out, q[2]), vertex(out, q[3])};
        double v = tet_volume(p[0], p[1], p[2], p[3]);
        vol += v; minvol = std::min(minvol, v);
        mind = std::min(mind, min_dihedral(p));
        static const int fc[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
        for (auto& f : fc) {
            int32_t s[3] = {q[f[0]], q[f[1]], q[f[2]]};
            std::sort(s, s + 3);
            ++faces[FK{s[0], s[1], s[2]}];
        }
    }
    }
    rep.volume_ratio = c.surf.volume != 0 ? vol / c.surf.volume : 0;
    rep.min_volume = rep.tets ? minvol : 0;
    rep.min_dihedral_deg = mind;
}

// Faces de bord (un seul tétraèdre) ; arêtes du bord qui ne sont pas partagées par exactement deux de ces faces.
std::vector<std::array<int32_t, 3>> boundary(const FaceCount& faces, Report& rep) {
    std::unordered_map<uint64_t, int> bedge;
    std::vector<std::array<int32_t, 3>> bfaces;
    for (auto& kv : faces) if (kv.second == 1) {
        const FK& fk = kv.first;
        bfaces.push_back({fk.a, fk.b, fk.c});
        int32_t vv[3] = {fk.a, fk.b, fk.c};
        for (int k = 0; k < 3; ++k) ++bedge[ekey(vv[k], vv[(k + 1) % 3])];
    }
    rep.boundary_faces = int64_t(bfaces.size());
    for (auto& kv : bedge) if (kv.second != 2) ++rep.nonmanifold_edges;
    return bfaces;
}

// Écart de peau : centroïde et milieux d'arêtes de chaque face de bord, distance à la surface.
void skin_error(const CheckContext& c, const TetMesh& out, const std::vector<std::array<int32_t, 3>>& bfaces) {
    std::vector<float> errs(bfaces.size(), 0.f);
    Span skin_span(c.prog, Verify, 0.70, 0.80, int64_t(bfaces.size()));   // faces de peau comparées à la surface
    const double h = c.h;
    parallel_for(int64_t(bfaces.size()), c.threads, [&](int64_t b, int64_t e, unsigned) {
        auto tick = skin_span.tick();
        for (int64_t f = b; f < e; ++f) {
            tick();
            V3 a = vertex(out, bfaces[f][0]), bb = vertex(out, bfaces[f][1]), cc = vertex(out, bfaces[f][2]);
            V3 s[4] = {(a + bb + cc) * (1.0 / 3), (a + bb) * 0.5, (bb + cc) * 0.5, (cc + a) * 0.5};
            double m = 0;
            for (auto& x : s) { Hit hit = nearest(c.bvh, x, 10 * h); m = std::max(m, hit.found ? std::sqrt(hit.d2) : 10 * h); }
            errs[f] = float(m);
        }
    });
    if (!errs.empty()) {
        c.rep.skin_max_error = *std::max_element(errs.begin(), errs.end());
        size_t k = size_t(0.99 * (errs.size() - 1));
        std::nth_element(errs.begin(), errs.begin() + k, errs.end());
        c.rep.skin_p99_error = errs[k];
    }
}

// Écart surface -> bord du maillage (matière manquante, coins rentrants comblés), sur la seule peau EXTÉRIEURE du
// solide : un point d'une face interne (coques qui se chevauchent) n'a pas à être au bord (classement : cover.hpp).
// Le test d'enroulement n'est fait qu'au-delà de `cover_shallow_alpha` h : plus près du bord, l'écart vient d'un coin
// coupé ou comblé, pas d'une face interne (elles sont plus profondes) ; c'est l'essentiel du coût évité.
void cover_error(const CheckContext& c, const TetMesh& out, const std::vector<std::array<int32_t, 3>>& bfaces) {
    if (bfaces.empty()) return;
    std::vector<V3> points(out.vertices.size() / 3);
    for (size_t v = 0; v < points.size(); ++v) points[v] = vertex(out, int32_t(v));
    const SkinIndex skin(points, bfaces);
    const std::vector<Sample> samples = surface_samples(c.surf, 0.5 * c.h);
    const double tol = c.prm.skin_tolerance, far = 10 * c.h, shallow = std::max(tol, c.prm.cover_shallow_alpha * c.h);
    const OuterSkin outer(c.surf, c.wn, 0.5 * tol);
    std::vector<float> errs(samples.size(), -1.f);   // -1 : point d'une face interne, écarté
    Span span(c.prog, Verify, 0.80, 0.99, int64_t(samples.size()));   // points de surface comparés au bord
    parallel_for(int64_t(samples.size()), c.threads, [&](int64_t b, int64_t e, unsigned) {
        auto tick = span.tick();
        for (int64_t i = b; i < e; ++i) {
            tick();
            const Hit hit = nearest(skin.bvh, samples[i].p, far);
            const double d = hit.found ? std::sqrt(hit.d2) : far;
            if (d <= shallow || outer.outer(samples[i].tri)) errs[i] = float(d);   // près du bord : peau, sans test
        }
    });
    std::vector<float> kept;
    kept.reserve(errs.size());
    for (float e : errs) if (e >= 0) kept.push_back(e);
    c.rep.cover_samples = int64_t(samples.size());
    c.rep.cover_internal_samples = int64_t(samples.size() - kept.size());
    if (kept.empty()) return;
    c.rep.skin_cover_max_error = *std::max_element(kept.begin(), kept.end());
    const size_t k = size_t(0.99 * (kept.size() - 1));
    std::nth_element(kept.begin(), kept.begin() + k, kept.end());
    c.rep.skin_cover_p99_error = kept[k];
}

}  // namespace

void check_mesh(const CheckContext& c, const TetMesh& out) {
    c.prog.report(Verify, 0.0);
    const auto t0 = Clock::now();
    FaceCount faces;
    measure_tets(c, out, faces);
    const std::vector<std::array<int32_t, 3>> bfaces = boundary(faces, c.rep);
    skin_error(c, out, bfaces);
    if (c.prm.check_cover) cover_error(c, out, bfaces);
    c.rep.t_checks = since(t0);
    c.prog.report(Verify, 1.0);                 // contrôles terminés
}

void final_status(const Params& prm, Report& rep) {
    rep.conflict_ratio = rep.band_vertices ? double(rep.sign_conflicts) / double(rep.band_vertices) : 0.0;
    rep.ambiguous_ratio = rep.band_vertices ? double(rep.ambiguous_points) / double(rep.band_vertices) : 0.0;
    const double judged = prm.use_winding ? rep.ambiguous_ratio : rep.conflict_ratio;
    rep.warning = rep.conflict_ratio > prm.warn_conflict_ratio || rep.ambiguous_ratio > prm.warn_ambiguous_ratio || rep.surface_reoriented_faces > 0;
    if (judged > prm.refuse_conflict_ratio) { rep.status = Status::InconsistentSurface; rep.reason = "surface incohérente : l'intérieur ne se distingue pas de l'extérieur"; return; }
    if (rep.tets == 0) { rep.status = Status::InvalidInput; rep.reason = "aucun tétraèdre produit"; }
    else if (rep.nonmanifold_edges > 0) { rep.status = Status::OutOfTolerance; rep.reason = "bord du maillage non fermé ou non manifold"; }
    else if (rep.skin_max_error > prm.skin_tolerance) { rep.status = Status::OutOfTolerance; rep.reason = "écart de peau supérieur à la tolérance"; }
    else if (rep.skin_cover_max_error > prm.skin_tolerance) { rep.status = Status::OutOfTolerance; rep.reason = "surface loin du bord du maillage (matière manquante ou coin comblé)"; }
    else rep.status = Status::Ok;
}

}  // namespace tetranet::detail
