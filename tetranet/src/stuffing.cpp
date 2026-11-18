// TetraNet — réseau BCC, tri dedans/dehors, ramené à la peau, découpes de Labelle-Shewchuk.
#include "stuffing.hpp"

#include <numeric>

namespace tetranet::detail {

bool make_lattice(Stuffing& c) {
    Lattice& L = c.lat;
    L.h = c.prm.cell_size;
    const double h = L.h;
    V3 lo(1e300, 1e300, 1e300), hi(-1e300, -1e300, -1e300);
    for (const V3& q : c.surf.p) { lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)}; hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)}; }
    const double pad = 2.0 * h;
    lo = lo - V3(pad, pad, pad);
    L.lo = lo;
    L.nx = int(std::ceil((hi.x + pad - lo.x) / h)); L.ny = int(std::ceil((hi.y + pad - lo.y) / h)); L.nz = int(std::ceil((hi.z + pad - lo.z) / h));
    const int nx = L.nx, ny = L.ny, nz = L.nz;
    L.NA = int64_t(nx + 1) * (ny + 1) * (nz + 1); L.NB = int64_t(nx) * ny * nz; L.N = L.NA + L.NB;
    if (L.N > 400000000LL) { c.rep.status = Status::TooLarge; c.rep.reason = "réseau trop grand : augmenter cell_size"; return false; }
    c.rep.lattice_vertices = L.N;
    std::vector<V3>& pos = c.st.pos;
    pos.assign(L.N, V3());
    const int64_t NA = L.NA;
    parallel_for(L.NB, c.threads, [&](int64_t b, int64_t e, unsigned) {
        for (int64_t q = b; q < e; ++q) { int i = q % nx, j = (q / nx) % ny, k = q / (int64_t(nx) * ny); pos[NA + q] = lo + V3((i + .5) * h, (j + .5) * h, (k + .5) * h); }
    });
    parallel_for(NA, c.threads, [&](int64_t b, int64_t e, unsigned) {
        for (int64_t q = b; q < e; ++q) { int i = q % (nx + 1), j = (q / (nx + 1)) % (ny + 1), k = q / (int64_t(nx + 1) * (ny + 1)); pos[q] = lo + V3(i * h, j * h, k * h); }
    });
    c.prog.report(Prepare, 0.12);               // réseau posé
    return true;
}

namespace {

// Sans nombre d'enroulement : signe des sommets hors bande par composantes (union-find sur les arêtes du réseau).
void signs_without_winding(Stuffing& c) {
    const Lattice& L = c.lat;
    const int64_t N = L.N;
    std::vector<int8_t>& sgn = c.st.sgn;
    const std::vector<int8_t>& isband = c.st.isband;
    std::vector<int64_t> parent(N);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int64_t x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    auto unite = [&](int64_t a, int64_t b) { a = find(a); b = find(b); if (a != b) parent[std::max(a, b)] = std::min(a, b); };
    std::vector<std::pair<int64_t, int64_t>> mixed;
    auto edge = [&](int64_t u, int64_t v) {
        if (!isband[u] && !isband[v]) unite(u, v);
        else if (!isband[u]) mixed.emplace_back(u, v);
        else if (!isband[v]) mixed.emplace_back(v, u);
    };
    for (int k = 0; k <= L.nz; ++k) for (int j = 0; j <= L.ny; ++j) for (int i = 0; i <= L.nx; ++i) {
        int64_t a = L.idA(i, j, k);
        if (i < L.nx) edge(a, L.idA(i + 1, j, k));
        if (j < L.ny) edge(a, L.idA(i, j + 1, k));
        if (k < L.nz) edge(a, L.idA(i, j, k + 1));
    }
    for (int k = 0; k < L.nz; ++k) for (int j = 0; j < L.ny; ++j) for (int i = 0; i < L.nx; ++i) {
        int64_t bq = L.idB(i, j, k);
        for (int dk = 0; dk < 2; ++dk) for (int dj = 0; dj < 2; ++dj) for (int di = 0; di < 2; ++di) edge(bq, L.idA(i + di, j + dj, k + dk));
        if (i + 1 < L.nx) edge(bq, L.idB(i + 1, j, k));
        if (j + 1 < L.ny) edge(bq, L.idB(i, j + 1, k));
        if (k + 1 < L.nz) edge(bq, L.idB(i, j, k + 1));
    }
    std::vector<int8_t> csign(N, 0);
    for (auto& [u, v] : mixed) {
        int64_t r = find(u); int8_t sv = sgn[v];
        if (sv == 0) continue;
        if (csign[r] == 0) csign[r] = sv; else if (csign[r] != sv) ++c.rep.sign_conflicts;
    }
    int64_t root0 = find(0);
    if (csign[root0] == 0) csign[root0] = 1;
    for (int64_t v = 0; v < N; ++v) { if (isband[v]) continue; int8_t s = csign[find(v)]; sgn[v] = s ? s : 1; }
}

}  // namespace

void classify_vertices(Stuffing& c) {
    const int64_t N = c.lat.N;
    const double h = c.lat.h, band = 1.05 * h;
    VertexState& st = c.st;
    st.sgn.assign(N, 0); st.isband.assign(N, 0); st.qpt.assign(N, V3()); st.dist.assign(N, 0.f); st.wval.assign(N, 0.f);
    const bool use_winding = c.prm.use_winding;
    if (use_winding) { c.wn.build(c.bvh); c.wn.beta = c.prm.winding_beta; }
    std::vector<int8_t> ps_sign(N, 0);            // signe par pseudonormale, pour mesurer les désaccords
    std::vector<uint8_t> ambiguous(N, 0);
    Span classify_span(c.prog, Prepare, 0.12, 0.70, N);   // sommets classés dedans/dehors
    parallel_for(N, c.threads, [&](int64_t b, int64_t e, unsigned) {
        auto tick = classify_span.tick();
        for (int64_t v = b; v < e; ++v) {
            tick();
            Hit hit = nearest(c.bvh, st.pos[v], band);
            if (hit.found) {
                st.isband[v] = 1; st.qpt[v] = hit.q; st.dist[v] = float(std::sqrt(hit.d2));
                ps_sign[v] = hit.d2 == 0 ? 0 : sign_at(c.surf, hit, st.pos[v]);
            }
            if (use_winding) {
                double w = c.wn.query(st.pos[v], hit.found ? 0.0 : 1.2);   // loin de la peau : dipôle grossier, le signe a de la marge
                st.wval[v] = float(w);
                st.sgn[v] = w > 0.5 ? -1 : 1;
                if (hit.found && std::fabs(w - 0.5) < 0.25 && hit.d2 > 4e-4 * h * h) ambiguous[v] = 1;   // un sommet SUR la peau vaut 0,5 sans que rien ne soit cassé
            } else if (hit.found) st.sgn[v] = ps_sign[v];
        }
    });
    if (use_winding) {
        for (int64_t v = 0; v < N; ++v) if (st.isband[v]) { if (ps_sign[v] != 0 && ps_sign[v] != st.sgn[v]) ++c.rep.sign_conflicts; if (ambiguous[v]) ++c.rep.ambiguous_points; }
    } else {
        signs_without_winding(c);
    }
    for (int64_t v = 0; v < N; ++v) if (st.isband[v]) ++c.rep.band_vertices;
    c.prog.report(Prepare, 0.75);               // dedans/dehors tranché, bande comptée
}

namespace {

// Sans découpes : un sommet à moins de snap_alpha * h de la peau y est ramené.
void snap_simple(Stuffing& c) {
    VertexState& st = c.st;
    const double snap = c.prm.snap_alpha * c.lat.h;
    Span snap_span(c.prog, Prepare, 0.75, 0.98, c.lat.N);
    parallel_for(c.lat.N, c.threads, [&](int64_t b, int64_t e, unsigned) {
        auto tick = snap_span.tick();
        for (int64_t v = b; v < e; ++v) if (tick(), st.isband[v] && st.dist[v] <= snap) { st.pos[v] = snap_to_feature(c.features.get(), st.qpt[v], c.feature_radius); st.sgn[v] = 0; }
    });
}

struct VEdge { int64_t u, v; bool is_short; };

// Arêtes du réseau dont les deux bouts sont de signes opposés (« violées »).
std::vector<VEdge> violated_edges(const Stuffing& c) {
    const Lattice& L = c.lat;
    const VertexState& st = c.st;
    std::vector<VEdge> viol;
    auto consider = [&](int64_t u, int64_t v, bool is_short) {
        if (st.isband[u] && st.isband[v] && st.sgn[u] * st.sgn[v] < 0) viol.push_back({u, v, is_short});
    };
    for (int k = 0; k <= L.nz; ++k) for (int j = 0; j <= L.ny; ++j) for (int i = 0; i <= L.nx; ++i) {
        int64_t a = L.idA(i, j, k);
        if (i < L.nx) consider(a, L.idA(i + 1, j, k), false);
        if (j < L.ny) consider(a, L.idA(i, j + 1, k), false);
        if (k < L.nz) consider(a, L.idA(i, j, k + 1), false);
    }
    for (int k = 0; k < L.nz; ++k) for (int j = 0; j < L.ny; ++j) for (int i = 0; i < L.nx; ++i) {
        int64_t bq = L.idB(i, j, k);
        for (int dk = 0; dk < 2; ++dk) for (int dj = 0; dj < 2; ++dj) for (int di = 0; di < 2; ++di) consider(bq, L.idA(i + di, j + dj, k + dk), true);
        if (i + 1 < L.nx) consider(bq, L.idB(i + 1, j, k), false);
        if (j + 1 < L.ny) consider(bq, L.idB(i, j + 1, k), false);
        if (k + 1 < L.nz) consider(bq, L.idB(i, j, k + 1), false);
    }
    return viol;
}

// Découpes : point de coupe de chaque arête violée ; ramenage des sommets proches ; sommets de coupe restants.
void cut_edges(Stuffing& c) {
    VertexState& st = c.st;
    const double h = c.lat.h;
    const int64_t N = c.lat.N;
    const std::vector<VEdge> viol = violated_edges(c);
    c.prog.report(Prepare, 0.80);           // arêtes à couper recensées
    // point de coupe de chaque arête : interpolation des distances, ramenée sur la peau
    std::vector<V3> cuts(viol.size()), refs(viol.size());
    Span cut_span(c.prog, Prepare, 0.80, 0.97, int64_t(viol.size()));   // points de coupe calculés
    parallel_for(int64_t(viol.size()), c.threads, [&](int64_t b, int64_t e, unsigned) {
        auto tick = cut_span.tick();
        for (int64_t k = b; k < e; ++k) {
            tick();
            const auto& ed = viol[k];
            double du = st.dist[ed.u], dv = st.dist[ed.v];
            double t = (du + dv) > 0 ? du / (du + dv) : 0.5;
            V3 lin = st.pos[ed.u] + (st.pos[ed.v] - st.pos[ed.u]) * t;
            refs[k] = lin;
            Hit hit = nearest(c.bvh, lin, 2.0 * h);
            cuts[k] = snap_to_feature(c.features.get(), hit.found ? hit.q : lin, c.feature_radius);
        }
    });
    // ramenage : un sommet à moins de alpha*L d'une coupe de ses arêtes violées va sur la coupe la plus proche
    std::vector<float> ratio(N, 1.0f);
    std::vector<V3> target(N);
    for (size_t k = 0; k < viol.size(); ++k) {
        const auto& ed = viol[k];
        const double L = norm(st.pos[ed.u] - st.pos[ed.v]), alpha = ed.is_short ? c.prm.alpha_short : c.prm.alpha_long;
        for (int64_t w : {ed.u, ed.v}) {
            double r = norm(st.pos[w] - cuts[k]) / (alpha * L);
            if (r < ratio[w]) { ratio[w] = float(r); target[w] = cuts[k]; }
        }
    }
    for (int64_t v = 0; v < N; ++v) if (ratio[v] < 1.0f) { st.pos[v] = target[v]; st.sgn[v] = 0; }
    // arêtes encore violées : un sommet de coupe chacune
    for (size_t k = 0; k < viol.size(); ++k) {
        const auto& ed = viol[k];
        if (st.sgn[ed.u] == 0 || st.sgn[ed.v] == 0) continue;
        c.cuts.id[ekey(int(ed.u), int(ed.v))] = N + int64_t(c.cuts.pos.size());
        c.cuts.pos.push_back(cuts[k]); c.cuts.ref.push_back(refs[k]);
    }
}

}  // namespace

void snap_to_skin(Stuffing& c) {
    c.st.pos0 = c.st.pos;                          // positions de référence (réseau régulier)
    if (!c.prm.use_cuts) snap_simple(c);
    else cut_edges(c);
    for (int64_t v = 0; v < c.lat.N; ++v) if (c.st.isband[v] && c.st.sgn[v] == 0) ++c.rep.snapped_vertices;
}

namespace {

struct Local { std::vector<int32_t> tets; int64_t dropped = 0; };

// Tétraèdres d'une paire de centres voisins : intérieurs gardés, traversés découpés (gabarits de Labelle-Shewchuk).
class Emitter {
public:
    explicit Emitter(const Stuffing& c)
        : c_(c), st_(c.st), N_(c.lat.N), vmin_(c.prm.min_volume_ratio * c.lat.h * c.lat.h * c.lat.h / 12.0) {}
    void any(Local& L, int64_t a, int64_t b, int64_t cc, int64_t d) const { if (c_.prm.use_cuts) cut(L, a, b, cc, d); else plain(L, a, b, cc, d); }
    const V3& P(int64_t v) const { return v < N_ ? st_.pos[v] : c_.cuts.pos[size_t(v - N_)]; }

private:
    const V3& R(int64_t v) const { return v < N_ ? st_.pos0[v] : c_.cuts.ref[size_t(v - N_)]; }
    int64_t cutv(int64_t u, int64_t v) const { auto it = c_.cuts.id.find(ekey(int(u), int(v))); return it == c_.cuts.id.end() ? -1 : it->second; }

    bool centroid_inside(const std::array<int64_t, 4>& t) const {
        V3 c = (st_.pos[t[0]] + st_.pos[t[1]] + st_.pos[t[2]] + st_.pos[t[3]]) * 0.25;
        if (c_.prm.use_winding) { double m = 0; for (int k = 0; k < 4; ++k) m += st_.wval[t[k]]; m *= 0.25; if (m > 0.8) return true; if (m < 0.2) return false; return c_.wn.query(c) > 0.5; }   // enroulement interpolé, requête exacte seulement si indécis
        Hit hit = nearest(c_.bvh, c, 4 * c_.lat.h);
        if (!hit.found) return false;
        return sign_at(c_.surf, hit, c) < 0;
    }

    void plain(Local& L, int64_t a, int64_t b, int64_t c, int64_t d) const {
        int sg[4] = {st_.sgn[a], st_.sgn[b], st_.sgn[c], st_.sgn[d]};
        int nin = 0, nout = 0;
        for (int s : sg) { nin += s < 0; nout += s > 0; }
        if (nin == 0 && nout == 4) return;
        std::array<int64_t, 4> t{a, b, c, d};
        bool all_in = (nin == 4);
        if (!all_in) {
            bool keep;
            if (nin > 0 && nout == 0) keep = true;
            else if (nin == 0 && nout > 0) keep = false;
            else keep = centroid_inside(t);
            if (!keep) return;
        }
        double v0 = tet_volume(st_.pos0[a], st_.pos0[b], st_.pos0[c], st_.pos0[d]);
        if (v0 < 0) std::swap(t[2], t[3]);
        double v1 = tet_volume(st_.pos[t[0]], st_.pos[t[1]], st_.pos[t[2]], st_.pos[t[3]]);
        if (!(v1 > vmin_)) { ++L.dropped; return; }
        for (int k = 0; k < 4; ++k) L.tets.push_back(int32_t(t[k]));
    }

    void push(Local& L, int64_t p0, int64_t p1, int64_t p2, int64_t p3) const {
        if (p0 < 0 || p1 < 0 || p2 < 0 || p3 < 0) { ++L.dropped; return; }
        if (tet_volume(R(p0), R(p1), R(p2), R(p3)) < 0) std::swap(p2, p3);     // orientation de référence positive
        if (!(tet_volume(P(p0), P(p1), P(p2), P(p3)) > vmin_)) { ++L.dropped; return; }
        L.tets.push_back(int32_t(p0)); L.tets.push_back(int32_t(p1)); L.tets.push_back(int32_t(p2)); L.tets.push_back(int32_t(p3));
    }

    // prisme (p0 p1 p2 | q0 q1 q2) en trois tétraèdres : les diagonales des faces carrées passent par le plus petit indice
    void prism(Local& L, int64_t t1[3], int64_t t2[3]) const {
        int64_t* a = t1; int64_t* b = t2;
        int64_t mn = a[0]; int where = 0; bool in_second = false;
        for (int k = 0; k < 3; ++k) { if (a[k] < mn) { mn = a[k]; where = k; in_second = false; } if (b[k] < mn) { mn = b[k]; where = k; in_second = true; } }
        if (in_second) std::swap(a, b);
        const int64_t p0 = a[where], p1 = a[(where + 1) % 3], p2 = a[(where + 2) % 3];
        const int64_t q0 = b[where], q1 = b[(where + 1) % 3], q2 = b[(where + 2) % 3];
        const int64_t m12 = std::min({p1, p2, q1, q2});
        if (m12 == p1 || m12 == q2) { push(L, p0, p1, p2, q2); push(L, p0, p1, q1, q2); push(L, p0, q0, q1, q2); }
        else                        { push(L, p0, p2, p1, q1); push(L, p0, p2, q2, q1); push(L, p0, q0, q2, q1); }
    }

    void cut(Local& L, int64_t a, int64_t b, int64_t c, int64_t d) const {
        const int64_t v[4] = {a, b, c, d};
        int64_t ng[4], ps[4], zr[4]; int nn = 0, np = 0, nz = 0;
        for (int k = 0; k < 4; ++k) { if (st_.sgn[v[k]] < 0) ng[nn++] = v[k]; else if (st_.sgn[v[k]] > 0) ps[np++] = v[k]; else zr[nz++] = v[k]; }
        if (np == 0) { if (nn > 0) push(L, a, b, c, d); return; }          // tout dedans (ou sur la peau)
        if (nn == 0) return;                                               // tout dehors
        if (nn == 1) {
            int64_t q[3]; int m = 0;
            for (int k = 0; k < np; ++k) q[m++] = cutv(ng[0], ps[k]);
            for (int k = 0; k < nz; ++k) q[m++] = zr[k];
            push(L, ng[0], q[0], q[1], q[2]);
        } else if (nn == 2 && np == 2) {
            int64_t t1[3] = {ng[0], cutv(ng[0], ps[0]), cutv(ng[0], ps[1])};
            int64_t t2[3] = {ng[1], cutv(ng[1], ps[0]), cutv(ng[1], ps[1])};
            if (t1[1] < 0 || t1[2] < 0 || t2[1] < 0 || t2[2] < 0) { ++L.dropped; return; }
            prism(L, t1, t2);
        } else if (nn == 3 && np == 1) {
            int64_t t1[3] = {ng[0], ng[1], ng[2]};
            int64_t t2[3] = {cutv(ng[0], ps[0]), cutv(ng[1], ps[0]), cutv(ng[2], ps[0])};
            if (t2[0] < 0 || t2[1] < 0 || t2[2] < 0) { ++L.dropped; return; }
            prism(L, t1, t2);
        } else if (nn == 2 && np == 1 && nz == 1) {
            const int64_t c1 = cutv(ng[0], ps[0]), c2 = cutv(ng[1], ps[0]);
            if (c1 < 0 || c2 < 0) { ++L.dropped; return; }
            const int64_t q0 = ng[0], q1 = ng[1], q2 = c2, q3 = c1;       // base : a, b, bc, ac
            const int64_t m = std::min({q0, q1, q2, q3});
            if (m == q0 || m == q2) { push(L, zr[0], q0, q1, q2); push(L, zr[0], q0, q2, q3); }
            else                    { push(L, zr[0], q1, q2, q3); push(L, zr[0], q1, q3, q0); }
        } else { ++L.dropped; }
    }

    const Stuffing& c_;
    const VertexState& st_;
    const int64_t N_;
    const double vmin_;
};

// Les quatre tétraèdres de chaque face partagée par le centre (i,j,k) et son voisin +x, +y, +z.
void emit_cell(const Emitter& em, const Lattice& L, Local& out, int i, int j, int k) {
    int64_t B1 = L.idB(i, j, k);
    if (i + 1 < L.nx) {  // face x = i+1, coins en (j,k)
        int64_t B2 = L.idB(i + 1, j, k);
        int64_t c00 = L.idA(i + 1, j, k), c10 = L.idA(i + 1, j + 1, k), c11 = L.idA(i + 1, j + 1, k + 1), c01 = L.idA(i + 1, j, k + 1);
        em.any(out, B1, B2, c00, c10); em.any(out, B1, B2, c10, c11); em.any(out, B1, B2, c11, c01); em.any(out, B1, B2, c01, c00);
    }
    if (j + 1 < L.ny) {
        int64_t B2 = L.idB(i, j + 1, k);
        int64_t c00 = L.idA(i, j + 1, k), c10 = L.idA(i + 1, j + 1, k), c11 = L.idA(i + 1, j + 1, k + 1), c01 = L.idA(i, j + 1, k + 1);
        em.any(out, B1, B2, c00, c10); em.any(out, B1, B2, c10, c11); em.any(out, B1, B2, c11, c01); em.any(out, B1, B2, c01, c00);
    }
    if (k + 1 < L.nz) {
        int64_t B2 = L.idB(i, j, k + 1);
        int64_t c00 = L.idA(i, j, k + 1), c10 = L.idA(i + 1, j, k + 1), c11 = L.idA(i + 1, j + 1, k + 1), c01 = L.idA(i, j + 1, k + 1);
        em.any(out, B1, B2, c00, c10); em.any(out, B1, B2, c10, c11); em.any(out, B1, B2, c11, c01); em.any(out, B1, B2, c01, c00);
    }
}

// Concatène les tétraèdres des fils, ne garde que les sommets utilisés (dans l'ordre de première apparition).
void compact(const Emitter& em, Stuffing& c, std::vector<Local>& locals, TetMesh& out) {
    std::vector<int32_t> tets;
    for (auto& L : locals) { tets.insert(tets.end(), L.tets.begin(), L.tets.end()); c.rep.dropped_degenerate += L.dropped; }
    const int64_t NT = c.lat.N + int64_t(c.cuts.pos.size());
    std::vector<int32_t> remap(NT, -1);
    int32_t used = 0;
    for (int32_t& t : tets) { if (remap[t] < 0) remap[t] = used++; }
    out.vertices.resize(size_t(used) * 3);
    for (int64_t v = 0; v < NT; ++v) if (remap[v] >= 0) { const V3& q = em.P(v); size_t o = size_t(remap[v]) * 3; out.vertices[o] = q.x; out.vertices[o + 1] = q.y; out.vertices[o + 2] = q.z; }
    for (int32_t& t : tets) t = remap[t];
    out.tetrahedra = std::move(tets);
    c.rep.tets = int64_t(out.tetrahedra.size() / 4);
}

}  // namespace

void emit_tetrahedra(Stuffing& c, TetMesh& out) {
    const Lattice& L = c.lat;
    std::vector<Local> locals(c.threads);
    const Emitter em(c);
    const int64_t cells = int64_t(L.nz) * L.ny * L.nx;
    Span tet_span(c.prog, Tetrahedralize, 0.0, 0.90, cells);   // cellules du réseau découpées
    parallel_for(cells, c.threads, [&](int64_t b, int64_t e, unsigned th) {
        Local& local = locals[th];
        auto tick = tet_span.tick();
        for (int64_t q = b; q < e; ++q) {
            tick();
            int i = q % L.nx, j = (q / L.nx) % L.ny, k = q / (int64_t(L.nx) * L.ny);
            emit_cell(em, L, local, i, j, k);
        }
    });
    compact(em, c, locals, out);
}

}  // namespace tetranet::detail
