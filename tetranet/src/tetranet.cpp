// TetraNet — prototype : réseau BCC, tri dedans/dehors, ramené à la peau (snap).
#include "tetranet/tetranet.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <numeric>
#include <thread>
#include <unordered_map>

namespace tetranet {
namespace {

using Clock = std::chrono::steady_clock;
double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

struct V3 {
    double x = 0, y = 0, z = 0;
    V3() = default;
    V3(double a, double b, double c) : x(a), y(b), z(c) {}
    V3 operator+(const V3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    V3 operator-(const V3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    V3 operator*(double s) const { return {x * s, y * s, z * s}; }
};
double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
V3 cross(const V3& a, const V3& b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
double norm(const V3& a) { return std::sqrt(dot(a, a)); }
V3 unit(const V3& a) { double n = norm(a); return n > 0 ? a * (1.0 / n) : V3(); }

template <class F>
void parallel_for(int64_t n, unsigned threads, F fn) {  // fn(begin, end, thread)
    if (threads <= 1 || n < 2048) { fn(0, n, 0u); return; }
    std::vector<std::thread> pool;
    int64_t chunk = (n + threads - 1) / threads;
    for (unsigned t = 0; t < threads; ++t) {
        int64_t b = t * chunk, e = std::min<int64_t>(n, b + chunk);
        if (b >= e) break;
        pool.emplace_back([=, &fn] { fn(b, e, t); });
    }
    for (auto& th : pool) th.join();
}

// ---------------------------------------------------------------- surface + BVH

struct Surface {
    std::vector<V3> p;                       // sommets soudés
    std::vector<std::array<int, 3>> tri;     // triangles non dégénérés
    std::vector<V3> fn;                      // normales de face (unitaires)
    std::vector<V3> vn;                      // pseudonormales de sommet
    std::unordered_map<uint64_t, V3> en;     // pseudonormales d'arête
    double volume = 0.0;
    int64_t holes = 0, open_edges = 0, nonmanifold_edges = 0, flipped_edges = 0, edges = 0, degenerate = 0, duplicates = 0;
};

uint64_t ekey(int a, int b) { if (a > b) std::swap(a, b); return (uint64_t(a) << 32) | uint32_t(b); }

bool build_surface(const SurfaceMesh& in, Surface& s) {
    size_t n = in.vertices.size() / 3, m = in.faces.size() / 3;
    if (n == 0 || m == 0) return false;
    struct Key { uint64_t a, b, c; bool operator==(const Key& o) const { return a == o.a && b == o.b && c == o.c; } };
    struct H { size_t operator()(const Key& k) const { uint64_t h = k.a * 0x9E3779B97F4A7C15ull; h ^= (k.b + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2)); h ^= (k.c + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2)); return size_t(h); } };
    std::unordered_map<Key, int, H> weld;
    std::vector<int> map(n);
    for (size_t i = 0; i < n; ++i) {
        double c[3] = {in.vertices[3 * i] + 0.0, in.vertices[3 * i + 1] + 0.0, in.vertices[3 * i + 2] + 0.0};  // +0.0 : -0 == +0
        Key key;
        std::memcpy(&key.a, &c[0], 8); std::memcpy(&key.b, &c[1], 8); std::memcpy(&key.c, &c[2], 8);
        auto it = weld.find(key);
        if (it == weld.end()) { it = weld.emplace(key, int(s.p.size())).first; s.p.push_back({c[0], c[1], c[2]}); }
        map[i] = it->second;
    }
    for (size_t f = 0; f < m; ++f) {
        std::array<int, 3> t{map[in.faces[3 * f]], map[in.faces[3 * f + 1]], map[in.faces[3 * f + 2]]};
        if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) { ++s.degenerate; continue; }
        V3 nn = cross(s.p[t[1]] - s.p[t[0]], s.p[t[2]] - s.p[t[0]]);
        double a2 = norm(nn);
        if (!(a2 > 1e-14)) { ++s.degenerate; continue; }
        s.tri.push_back(t);
        s.fn.push_back(nn * (1.0 / a2));
        s.volume += dot(s.p[t[0]], cross(s.p[t[1]], s.p[t[2]])) / 6.0;
    }
    if (s.tri.empty()) return false;
    {   // topologie des arêtes : (nombre de faces, sens parcourus)
        struct EI { int count = 0; int fwd = 0; };
        std::unordered_map<uint64_t, EI> ed;
        ed.reserve(s.tri.size() * 2);
        std::unordered_map<uint64_t, int> faces;
        for (const auto& t : s.tri) {
            for (int k = 0; k < 3; ++k) { int a = t[k], b = t[(k + 1) % 3]; auto& e = ed[ekey(a, b)]; ++e.count; e.fwd += (a < b) ? 1 : 0; }
            int a = t[0], b = t[1], c = t[2];
            if (a > b) std::swap(a, b); if (b > c) std::swap(b, c); if (a > b) std::swap(a, b);
            if (++faces[(uint64_t(a) * 1000003ull + uint64_t(b)) * 1000033ull + uint64_t(c)] > 1) ++s.duplicates;
        }
        s.edges = int64_t(ed.size());
        std::unordered_map<int, int> uf;   // union-find sur les sommets des arêtes ouvertes
        auto ufind = [&](int x) { while (true) { auto it = uf.find(x); if (it == uf.end()) { uf[x] = x; return x; } if (it->second == x) return x; int r = it->second; auto j = uf.find(r); if (j != uf.end() && j->second != r) it->second = j->second; x = r; } };
        for (auto& kv : ed) if (kv.second.count == 1) { int a = int(kv.first >> 32), b = int(kv.first & 0xffffffffu); int ra = ufind(a), rb = ufind(b); if (ra != rb) uf[ra] = rb; }
        { std::unordered_map<int, int> roots; for (auto& kv : uf) roots[ufind(kv.first)] = 1; s.holes = int64_t(roots.size()); }
        for (auto& kv : ed) {
            if (kv.second.count == 1) ++s.open_edges;
            else if (kv.second.count > 2) ++s.nonmanifold_edges;
            else if (kv.second.fwd != 1) ++s.flipped_edges;   // 2 faces : un sens aller, un sens retour attendus
        }
    }
    s.vn.assign(s.p.size(), V3());
    for (size_t f = 0; f < s.tri.size(); ++f) {
        for (int k = 0; k < 3; ++k) {
            int a = s.tri[f][k], b = s.tri[f][(k + 1) % 3], c = s.tri[f][(k + 2) % 3];
            s.en[ekey(a, b)] = s.en[ekey(a, b)] + s.fn[f];
            V3 e1 = unit(s.p[b] - s.p[a]), e2 = unit(s.p[c] - s.p[a]);
            double ang = std::acos(std::max(-1.0, std::min(1.0, dot(e1, e2))));
            s.vn[a] = s.vn[a] + s.fn[f] * ang;
        }
    }
    return true;
}

struct Node { double lo[3], hi[3]; int left = -1, right = -1, first = 0, count = 0; };

struct Bvh {
    const Surface* s = nullptr;
    std::vector<Node> nodes;
    std::vector<int> order;

    void build(const Surface& surf) {
        s = &surf;
        order.resize(surf.tri.size());
        std::iota(order.begin(), order.end(), 0);
        nodes.reserve(surf.tri.size());
        rec(0, int(order.size()));
    }
    int rec(int b, int e) {
        int id = int(nodes.size());
        nodes.emplace_back();
        Node nd;
        for (int k = 0; k < 3; ++k) { nd.lo[k] = 1e300; nd.hi[k] = -1e300; }
        double clo[3] = {1e300, 1e300, 1e300}, chi[3] = {-1e300, -1e300, -1e300};
        for (int i = b; i < e; ++i) {
            const auto& t = s->tri[order[i]];
            double c[3] = {0, 0, 0};
            for (int v = 0; v < 3; ++v) {
                const V3& q = s->p[t[v]];
                double a[3] = {q.x, q.y, q.z};
                for (int k = 0; k < 3; ++k) { nd.lo[k] = std::min(nd.lo[k], a[k]); nd.hi[k] = std::max(nd.hi[k], a[k]); c[k] += a[k] / 3; }
            }
            for (int k = 0; k < 3; ++k) { clo[k] = std::min(clo[k], c[k]); chi[k] = std::max(chi[k], c[k]); }
        }
        if (e - b <= 4) { nd.first = b; nd.count = e - b; nodes[id] = nd; return id; }
        int ax = 0;
        for (int k = 1; k < 3; ++k) if (chi[k] - clo[k] > chi[ax] - clo[ax]) ax = k;
        int mid = (b + e) / 2;
        auto cen = [&](int ti) { const auto& t = s->tri[ti]; const V3 &a = s->p[t[0]], &bb = s->p[t[1]], &c = s->p[t[2]];
            return ax == 0 ? a.x + bb.x + c.x : ax == 1 ? a.y + bb.y + c.y : a.z + bb.z + c.z; };
        std::nth_element(order.begin() + b, order.begin() + mid, order.begin() + e, [&](int u, int v) { return cen(u) < cen(v); });
        nd.left = rec(b, mid);
        nd.right = rec(mid, e);
        nodes[id] = nd;
        return id;
    }
    static double box2(const Node& n, const V3& q) {
        double d = 0, a[3] = {q.x, q.y, q.z};
        for (int k = 0; k < 3; ++k) { double v = std::max({n.lo[k] - a[k], 0.0, a[k] - n.hi[k]}); d += v * v; }
        return d;
    }
};

// point le plus proche d'un triangle (Ericson) ; region : 0 face, 1..3 arêtes ab/bc/ca, 4..6 sommets
V3 closest_on_tri(const V3& p, const V3& a, const V3& b, const V3& c, int& region) {
    V3 ab = b - a, ac = c - a, ap = p - a;
    double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) { region = 4; return a; }
    V3 bp = p - b;
    double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) { region = 5; return b; }
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) { region = 1; return a + ab * (d1 / (d1 - d3)); }
    V3 cp = p - c;
    double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) { region = 6; return c; }
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) { region = 3; return a + ac * (d2 / (d2 - d6)); }
    double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) { region = 2; return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6))); }
    double den = 1.0 / (va + vb + vc);
    region = 0;
    return a + ab * (vb * den) + ac * (vc * den);
}

struct Hit { bool found = false; double d2 = 0; V3 q; int tri = -1, region = 0; };

Hit nearest(const Bvh& bvh, const V3& p, double maxd) {
    Hit h;
    h.d2 = maxd * maxd;
    int stack[128], sp = 0;
    stack[sp++] = 0;
    while (sp) {
        const Node& n = bvh.nodes[stack[--sp]];
        if (Bvh::box2(n, p) >= h.d2) continue;
        if (n.left < 0) {
            for (int i = n.first; i < n.first + n.count; ++i) {
                int ti = bvh.order[i];
                const auto& t = bvh.s->tri[ti];
                int reg;
                V3 q = closest_on_tri(p, bvh.s->p[t[0]], bvh.s->p[t[1]], bvh.s->p[t[2]], reg);
                V3 d = p - q;
                double d2 = dot(d, d);
                if (d2 < h.d2) { h.d2 = d2; h.q = q; h.tri = ti; h.region = reg; h.found = true; }
            }
        } else {
            double dl = Bvh::box2(bvh.nodes[n.left], p), dr = Bvh::box2(bvh.nodes[n.right], p);
            if (dl < dr) { stack[sp++] = n.right; stack[sp++] = n.left; }
            else { stack[sp++] = n.left; stack[sp++] = n.right; }
        }
    }
    return h;
}

// +1 dehors, -1 dedans (pseudonormale), 0 si p est sur la surface
int sign_at(const Surface& s, const Hit& h, const V3& p) {
    const auto& t = s.tri[h.tri];
    V3 n;
    switch (h.region) {
        case 0: n = s.fn[h.tri]; break;
        case 1: n = s.en.at(ekey(t[0], t[1])); break;
        case 2: n = s.en.at(ekey(t[1], t[2])); break;
        case 3: n = s.en.at(ekey(t[2], t[0])); break;
        case 4: n = s.vn[t[0]]; break;
        case 5: n = s.vn[t[1]]; break;
        default: n = s.vn[t[2]]; break;
    }
    double d = dot(p - h.q, n);
    return d > 0 ? 1 : -1;
}

double tet_volume(const V3& a, const V3& b, const V3& c, const V3& d) { return dot(cross(b - a, c - a), d - a) / 6.0; }

}  // namespace

const char* version() { return "0.0.1"; }

Report mesh(const SurfaceMesh& in, const Params& prm, TetMesh& out) {
    auto t_all = Clock::now();
    Report rep;
    out.vertices.clear();
    out.tetrahedra.clear();
    unsigned threads = prm.threads ? prm.threads : std::max(1u, std::thread::hardware_concurrency());

    Surface surf;
    if (!(prm.cell_size > 0) || !build_surface(in, surf)) {
        rep.status = Status::InvalidInput;
        rep.reason = "surface vide ou taille de maille invalide";
        return rep;
    }
    auto t0 = Clock::now();
    Bvh bvh;
    bvh.build(surf);
    rep.t_bvh = since(t0);
    rep.surface_faces = int64_t(surf.tri.size()); rep.surface_edges = surf.edges; rep.surface_open_edges = surf.open_edges; rep.surface_holes = surf.holes;
    rep.surface_nonmanifold_edges = surf.nonmanifold_edges; rep.surface_flipped_edges = surf.flipped_edges;
    rep.surface_degenerate_faces = surf.degenerate; rep.surface_duplicate_faces = surf.duplicates;

    // --- réseau BCC
    const double h = prm.cell_size;
    V3 lo(1e300, 1e300, 1e300), hi(-1e300, -1e300, -1e300);
    for (const V3& q : surf.p) { lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)}; hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)}; }
    const double pad = 2.0 * h;
    lo = lo - V3(pad, pad, pad);
    int nx = int(std::ceil((hi.x + pad - lo.x) / h)), ny = int(std::ceil((hi.y + pad - lo.y) / h)), nz = int(std::ceil((hi.z + pad - lo.z) / h));
    int64_t NA = int64_t(nx + 1) * (ny + 1) * (nz + 1), NB = int64_t(nx) * ny * nz, N = NA + NB;
    if (N > 400000000LL) { rep.status = Status::TooLarge; rep.reason = "réseau trop grand : augmenter cell_size"; return rep; }
    rep.lattice_vertices = N;
    auto idA = [&](int i, int j, int k) -> int64_t { return (int64_t(k) * (ny + 1) + j) * (nx + 1) + i; };
    auto idB = [&](int i, int j, int k) -> int64_t { return NA + (int64_t(k) * ny + j) * nx + i; };
    std::vector<V3> pos(N);
    parallel_for(NB, threads, [&](int64_t b, int64_t e, unsigned) {
        for (int64_t q = b; q < e; ++q) { int i = q % nx, j = (q / nx) % ny, k = q / (int64_t(nx) * ny); pos[NA + q] = lo + V3((i + .5) * h, (j + .5) * h, (k + .5) * h); }
    });
    parallel_for(NA, threads, [&](int64_t b, int64_t e, unsigned) {
        for (int64_t q = b; q < e; ++q) { int i = q % (nx + 1), j = (q / (nx + 1)) % (ny + 1), k = q / (int64_t(nx + 1) * (ny + 1)); pos[q] = lo + V3(i * h, j * h, k * h); }
    });

    // --- sommets de la bande : signe et point le plus proche
    t0 = Clock::now();
    const double band = 1.05 * h;
    std::vector<int8_t> sgn(N, 0), isband(N, 0);  // sgn : +1 dehors, -1 dedans, 0 sur la peau
    std::vector<V3> qpt(N);
    std::vector<float> dist(N, 0.f);
    parallel_for(N, threads, [&](int64_t b, int64_t e, unsigned) {
        for (int64_t v = b; v < e; ++v) {
            Hit hit = nearest(bvh, pos[v], band);
            if (!hit.found) continue;
            isband[v] = 1;
            sgn[v] = hit.d2 == 0 ? 0 : sign_at(surf, hit, pos[v]);
            qpt[v] = hit.q;
            dist[v] = float(std::sqrt(hit.d2));
        }
    });
    // --- composantes des sommets hors bande (union-find sur les arêtes du réseau)
    std::vector<int64_t> parent(N);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int64_t x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    auto unite = [&](int64_t a, int64_t b) { a = find(a); b = find(b); if (a != b) parent[std::max(a, b)] = std::min(a, b); };
    std::vector<std::pair<int64_t, int64_t>> mixed;  // (hors bande, bande)
    auto edge = [&](int64_t u, int64_t v) {
        if (!isband[u] && !isband[v]) unite(u, v);
        else if (!isband[u]) mixed.emplace_back(u, v);
        else if (!isband[v]) mixed.emplace_back(v, u);
    };
    for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i) {
        int64_t a = idA(i, j, k);
        if (i < nx) edge(a, idA(i + 1, j, k));
        if (j < ny) edge(a, idA(i, j + 1, k));
        if (k < nz) edge(a, idA(i, j, k + 1));
    }
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        int64_t bq = idB(i, j, k);
        for (int dk = 0; dk < 2; ++dk) for (int dj = 0; dj < 2; ++dj) for (int di = 0; di < 2; ++di) edge(bq, idA(i + di, j + dj, k + dk));
        if (i + 1 < nx) edge(bq, idB(i + 1, j, k));
        if (j + 1 < ny) edge(bq, idB(i, j + 1, k));
        if (k + 1 < nz) edge(bq, idB(i, j, k + 1));
    }
    std::vector<int8_t> csign(N, 0);
    for (auto& [u, v] : mixed) {
        int64_t r = find(u);
        int8_t sv = sgn[v];
        if (sv == 0) continue;
        if (csign[r] == 0) csign[r] = sv;
        else if (csign[r] != sv) ++rep.sign_conflicts;
    }
    int64_t root0 = find(0);
    if (csign[root0] == 0) csign[root0] = 1;  // le coin du réseau est dehors
    for (int64_t v = 0; v < N; ++v) {
        if (isband[v]) { ++rep.band_vertices; continue; }
        int8_t s = csign[find(v)];
        sgn[v] = s ? s : 1;
    }
    // --- ramener à la peau les sommets proches
    const double snap = prm.snap_alpha * h;
    std::vector<V3> pos0;  // positions d'origine des seuls sommets déplacés
    pos0 = pos;
    parallel_for(N, threads, [&](int64_t b, int64_t e, unsigned) {
        for (int64_t v = b; v < e; ++v) if (isband[v] && dist[v] <= snap) { pos[v] = qpt[v]; sgn[v] = 0; }
    });
    for (int64_t v = 0; v < N; ++v) if (isband[v] && sgn[v] == 0) ++rep.snapped_vertices;
    rep.t_classify = since(t0);

    // --- tétraèdres : paires de centres voisins + arête de la face partagée
    t0 = Clock::now();
    const double vmin = prm.min_volume_ratio * h * h * h / 12.0;
    struct Local { std::vector<int32_t> tets; int64_t dropped = 0; };
    std::vector<Local> locals(threads);
    auto centroid_inside = [&](const std::array<int64_t, 4>& t) {
        V3 c = (pos[t[0]] + pos[t[1]] + pos[t[2]] + pos[t[3]]) * 0.25;
        Hit hit = nearest(bvh, c, 4 * h);
        if (!hit.found) return false;
        return sign_at(surf, hit, c) < 0;
    };
    auto emit = [&](Local& L, int64_t a, int64_t b, int64_t c, int64_t d) {
        int sg[4] = {sgn[a], sgn[b], sgn[c], sgn[d]};
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
        double v0 = tet_volume(pos0[a], pos0[b], pos0[c], pos0[d]);
        if (v0 < 0) std::swap(t[2], t[3]);
        double v1 = tet_volume(pos[t[0]], pos[t[1]], pos[t[2]], pos[t[3]]);
        if (!(v1 > vmin)) { ++L.dropped; return; }
        for (int k = 0; k < 4; ++k) L.tets.push_back(int32_t(t[k]));
    };
    parallel_for(int64_t(nz) * ny * nx, threads, [&](int64_t b, int64_t e, unsigned th) {
        Local& L = locals[th];
        for (int64_t q = b; q < e; ++q) {
            int i = q % nx, j = (q / nx) % ny, k = q / (int64_t(nx) * ny);
            int64_t B1 = idB(i, j, k);
            if (i + 1 < nx) {  // face x = i+1, coins en (j,k)
                int64_t B2 = idB(i + 1, j, k);
                int64_t c00 = idA(i + 1, j, k), c10 = idA(i + 1, j + 1, k), c11 = idA(i + 1, j + 1, k + 1), c01 = idA(i + 1, j, k + 1);
                emit(L, B1, B2, c00, c10); emit(L, B1, B2, c10, c11); emit(L, B1, B2, c11, c01); emit(L, B1, B2, c01, c00);
            }
            if (j + 1 < ny) {
                int64_t B2 = idB(i, j + 1, k);
                int64_t c00 = idA(i, j + 1, k), c10 = idA(i + 1, j + 1, k), c11 = idA(i + 1, j + 1, k + 1), c01 = idA(i, j + 1, k + 1);
                emit(L, B1, B2, c00, c10); emit(L, B1, B2, c10, c11); emit(L, B1, B2, c11, c01); emit(L, B1, B2, c01, c00);
            }
            if (k + 1 < nz) {
                int64_t B2 = idB(i, j, k + 1);
                int64_t c00 = idA(i, j, k + 1), c10 = idA(i + 1, j, k + 1), c11 = idA(i + 1, j + 1, k + 1), c01 = idA(i, j + 1, k + 1);
                emit(L, B1, B2, c00, c10); emit(L, B1, B2, c10, c11); emit(L, B1, B2, c11, c01); emit(L, B1, B2, c01, c00);
            }
        }
    });
    std::vector<int32_t> tets;
    for (auto& L : locals) { tets.insert(tets.end(), L.tets.begin(), L.tets.end()); rep.dropped_degenerate += L.dropped; }
    pos0.clear(); pos0.shrink_to_fit();
    // compactage des sommets
    std::vector<int32_t> remap(N, -1);
    int32_t used = 0;
    for (int32_t& t : tets) { if (remap[t] < 0) remap[t] = used++; }
    out.vertices.resize(size_t(used) * 3);
    for (int64_t v = 0; v < N; ++v) if (remap[v] >= 0) { size_t o = size_t(remap[v]) * 3; out.vertices[o] = pos[v].x; out.vertices[o + 1] = pos[v].y; out.vertices[o + 2] = pos[v].z; }
    for (int32_t& t : tets) t = remap[t];
    out.tetrahedra = std::move(tets);
    rep.tets = int64_t(out.tetrahedra.size() / 4);
    rep.t_tets = since(t0);

    // --- contrôles
    t0 = Clock::now();
    auto P = [&](int32_t i) { return V3(out.vertices[3 * size_t(i)], out.vertices[3 * size_t(i) + 1], out.vertices[3 * size_t(i) + 2]); };
    double vol = 0, minvol = 1e300, mind = 180.0;
    struct FK { int32_t a, b, c; bool operator==(const FK& o) const { return a == o.a && b == o.b && c == o.c; } };
    struct FH { size_t operator()(const FK& k) const { return size_t((uint64_t(uint32_t(k.a)) * 0x9E3779B97F4A7C15ull) ^ (uint64_t(uint32_t(k.b)) << 21) ^ (uint64_t(uint32_t(k.c)) * 0xC2B2AE3D27D4EB4Full)); } };
    std::unordered_map<FK, int32_t, FH> faces;  // face triée -> nombre de tétraèdres
    for (int64_t t = 0; t < rep.tets; ++t) {
        const int32_t* q = &out.tetrahedra[4 * size_t(t)];
        V3 p[4] = {P(q[0]), P(q[1]), P(q[2]), P(q[3])};
        double v = tet_volume(p[0], p[1], p[2], p[3]);
        vol += v; minvol = std::min(minvol, v);
        static const int pr[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
        for (auto& e : pr) {
            int c = -1, d = -1;
            for (int k = 0; k < 4; ++k) if (k != e[0] && k != e[1]) (c < 0 ? c : d) = k;
            V3 ev = unit(p[e[1]] - p[e[0]]);
            V3 u = p[c] - p[e[0]]; u = u - ev * dot(u, ev);
            V3 w = p[d] - p[e[0]]; w = w - ev * dot(w, ev);
            double den = norm(u) * norm(w);
            if (den > 0) mind = std::min(mind, std::acos(std::max(-1.0, std::min(1.0, dot(u, w) / den))) * 180.0 / M_PI);
        }
        static const int fc[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
        for (auto& f : fc) {
            int32_t s[3] = {q[f[0]], q[f[1]], q[f[2]]};
            std::sort(s, s + 3);
            ++faces[FK{s[0], s[1], s[2]}];
        }
    }
    rep.volume_ratio = surf.volume != 0 ? vol / surf.volume : 0;
    rep.min_volume = rep.tets ? minvol : 0;
    rep.min_dihedral_deg = mind;
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
    // écart de peau : centroïde et milieux d'arêtes de chaque face de bord, distance à la surface
    std::vector<float> errs(bfaces.size(), 0.f);
    parallel_for(int64_t(bfaces.size()), threads, [&](int64_t b, int64_t e, unsigned) {
        for (int64_t f = b; f < e; ++f) {
            V3 a = P(bfaces[f][0]), bb = P(bfaces[f][1]), c = P(bfaces[f][2]);
            V3 s[4] = {(a + bb + c) * (1.0 / 3), (a + bb) * 0.5, (bb + c) * 0.5, (c + a) * 0.5};
            double m = 0;
            for (auto& x : s) { Hit hit = nearest(bvh, x, 10 * h); m = std::max(m, hit.found ? std::sqrt(hit.d2) : 10 * h); }
            errs[f] = float(m);
        }
    });
    if (!errs.empty()) {
        rep.skin_max_error = *std::max_element(errs.begin(), errs.end());
        size_t k = size_t(0.99 * (errs.size() - 1));
        std::nth_element(errs.begin(), errs.begin() + k, errs.end());
        rep.skin_p99_error = errs[k];
    }
    rep.t_checks = since(t0);
    rep.seconds = since(t_all);
    rep.conflict_ratio = rep.band_vertices ? double(rep.sign_conflicts) / double(rep.band_vertices) : 0.0;
    rep.warning = rep.conflict_ratio > prm.warn_conflict_ratio;
    if (rep.conflict_ratio > prm.refuse_conflict_ratio) { rep.status = Status::InconsistentSurface; rep.reason = "surface incohérente : l'intérieur ne se distingue pas de l'extérieur"; return rep; }
    if (rep.tets == 0) { rep.status = Status::InvalidInput; rep.reason = "aucun tétraèdre produit"; }
    else if (rep.nonmanifold_edges > 0) { rep.status = Status::OutOfTolerance; rep.reason = "bord du maillage non fermé ou non manifold"; }
    else if (rep.skin_max_error > prm.skin_tolerance) { rep.status = Status::OutOfTolerance; rep.reason = "écart de peau supérieur à la tolérance"; }
    else rep.status = Status::Ok;
    return rep;
}

}  // namespace tetranet
