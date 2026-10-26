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
    int64_t reoriented = 0, holes = 0, open_edges = 0, nonmanifold_edges = 0, flipped_edges = 0, edges = 0, degenerate = 0, duplicates = 0;
};

uint64_t ekey(int a, int b) { if (a > b) std::swap(a, b); return (uint64_t(a) << 32) | uint32_t(b); }

bool build_surface(const SurfaceMesh& in, Surface& s, bool fix_orientation) {
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
    if (fix_orientation) {
        // propagation de l'orientation sur les arêtes à deux faces (BFS), puis volume positif par composante
        const size_t T = s.tri.size();
        std::unordered_map<uint64_t, std::array<int, 2>> adj;   // arête -> faces (-1 si absente)
        std::unordered_map<uint64_t, int> cnt;
        adj.reserve(T * 2);
        for (size_t f = 0; f < T; ++f) for (int k = 0; k < 3; ++k) {
            uint64_t key = ekey(s.tri[f][k], s.tri[f][(k + 1) % 3]);
            int& c = cnt[key];
            auto& a = adj.try_emplace(key, std::array<int, 2>{-1, -1}).first->second;
            if (c < 2) a[c] = int(f);
            ++c;
        }
        auto dir_fwd = [&](size_t f, uint64_t key) {   // la face f parcourt-elle l'arête de petit vers grand indice ?
            for (int k = 0; k < 3; ++k) { int a = s.tri[f][k], b = s.tri[f][(k + 1) % 3]; if (ekey(a, b) == key) return a < b; }
            return true;
        };
        std::vector<int> comp(T, -1);
        std::vector<int> stack;
        int ncomp = 0;
        for (size_t seed = 0; seed < T; ++seed) {
            if (comp[seed] >= 0) continue;
            comp[seed] = ncomp; stack.assign(1, int(seed));
            std::vector<int> members;
            while (!stack.empty()) {
                int f = stack.back(); stack.pop_back(); members.push_back(f);
                for (int k = 0; k < 3; ++k) {
                    uint64_t key = ekey(s.tri[f][k], s.tri[f][(k + 1) % 3]);
                    if (cnt[key] != 2) continue;                 // seules les arêtes manifold propagent
                    const auto& a = adj[key];
                    int g = a[0] == f ? a[1] : a[0];
                    if (g < 0 || comp[g] >= 0) continue;
                    if (dir_fwd(size_t(g), key) == dir_fwd(size_t(f), key)) { std::swap(s.tri[g][1], s.tri[g][2]); ++s.reoriented; }
                    comp[g] = ncomp; stack.push_back(g);
                }
            }
            double vol = 0;
            for (int f : members) { const auto& t = s.tri[f]; vol += dot(s.p[t[0]], cross(s.p[t[1]], s.p[t[2]])) / 6.0; }
            if (vol < 0) for (int f : members) { std::swap(s.tri[f][1], s.tri[f][2]); }
            ++ncomp;
        }
        s.volume = 0;
        for (size_t f = 0; f < T; ++f) {
            const auto& t = s.tri[f];
            V3 nn = cross(s.p[t[1]] - s.p[t[0]], s.p[t[2]] - s.p[t[0]]);
            s.fn[f] = nn * (1.0 / norm(nn));
            s.volume += dot(s.p[t[0]], cross(s.p[t[1]], s.p[t[2]])) / 6.0;
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

// ---------------------------------------------------------------- nombre d'enroulement généralisé
// w(p) = (1/4π) Σ angle solide des triangles vus de p : 1 dedans, 0 dehors, entre les deux près des défauts.
// Les nœuds lointains sont remplacés par un dipôle (Barill et al., 2018) : normale pondérée par l'aire.
struct WData { V3 c, n; double area = 0, r = 0; };

struct Winding {
    const Bvh* bvh = nullptr;
    std::vector<WData> d;
    double beta = 3.0;

    void build(const Bvh& b) {
        bvh = &b;
        const Surface& s = *b.s;
        d.assign(b.nodes.size(), WData());
        for (int i = int(b.nodes.size()) - 1; i >= 0; --i) {   // les enfants ont un indice plus grand que leur parent
            const Node& nd = b.nodes[i];
            WData& w = d[i];
            if (nd.left < 0) {
                V3 csum; double asum = 0; V3 nsum;
                for (int k = nd.first; k < nd.first + nd.count; ++k) {
                    int ti = b.order[k]; const auto& t = s.tri[ti];
                    const V3 &a = s.p[t[0]], &bb = s.p[t[1]], &c = s.p[t[2]];
                    double ar = 0.5 * norm(cross(bb - a, c - a));
                    csum = csum + (a + bb + c) * (ar / 3.0); asum += ar; nsum = nsum + s.fn[ti] * ar;
                }
                w.area = asum; w.n = nsum; w.c = asum > 0 ? csum * (1.0 / asum) : V3();
                double r = 0;
                for (int k = nd.first; k < nd.first + nd.count; ++k) { const auto& t = s.tri[b.order[k]]; for (int v = 0; v < 3; ++v) r = std::max(r, norm(s.p[t[v]] - w.c)); }
                w.r = r;
            } else {
                const WData &l = d[nd.left], &rr = d[nd.right];
                w.area = l.area + rr.area; w.n = l.n + rr.n;
                w.c = w.area > 0 ? (l.c * l.area + rr.c * rr.area) * (1.0 / w.area) : V3();
                w.r = std::max(norm(l.c - w.c) + l.r, norm(rr.c - w.c) + rr.r);
            }
        }
    }
    static double solid(const V3& a, const V3& b, const V3& c) {   // Van Oosterom - Strakhov
        double la = norm(a), lb = norm(b), lc = norm(c);
        double den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
        return 2.0 * std::atan2(dot(a, cross(b, c)), den);
    }
    double query(const V3& p, double beta_override = 0.0) const {
        const double beta = beta_override > 0 ? beta_override : this->beta;
        const Surface& s = *bvh->s;
        double total = 0;
        int stack[256], sp = 0; stack[sp++] = 0;
        while (sp) {
            int i = stack[--sp]; const Node& nd = bvh->nodes[i]; const WData& w = d[i];
            V3 cp = w.c - p; double dist = norm(cp);
            if (dist > beta * w.r && w.area > 0) { total += dot(w.n, cp) / (dist * dist * dist); continue; }
            if (nd.left < 0) {
                for (int k = nd.first; k < nd.first + nd.count; ++k) {
                    const auto& t = s.tri[bvh->order[k]];
                    total += solid(s.p[t[0]] - p, s.p[t[1]] - p, s.p[t[2]] - p);
                }
            } else { stack[sp++] = nd.left; stack[sp++] = nd.right; }
        }
        return total / (4.0 * M_PI);
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

Report mesh(const SurfaceMesh& in, const Params& prm, TetMesh& out, SurfaceMesh* cleaned) {
    auto t_all = Clock::now();
    Report rep;
    out.vertices.clear();
    out.tetrahedra.clear();
    unsigned threads = prm.threads ? prm.threads : std::max(1u, std::thread::hardware_concurrency());

    Surface surf;
    if (!(prm.cell_size > 0) || !build_surface(in, surf, prm.fix_orientation)) {
        rep.status = Status::InvalidInput;
        rep.reason = "surface vide ou taille de maille invalide";
        return rep;
    }
    if (cleaned) {
        cleaned->vertices.clear(); cleaned->faces.clear();
        cleaned->vertices.reserve(surf.p.size() * 3); cleaned->faces.reserve(surf.tri.size() * 3);
        for (const V3& q : surf.p) { cleaned->vertices.push_back(q.x); cleaned->vertices.push_back(q.y); cleaned->vertices.push_back(q.z); }
        for (const auto& t : surf.tri) { cleaned->faces.push_back(t[0]); cleaned->faces.push_back(t[1]); cleaned->faces.push_back(t[2]); }
    }
    auto t0 = Clock::now();
    Bvh bvh;
    bvh.build(surf);
    rep.t_bvh = since(t0);
    rep.surface_faces = int64_t(surf.tri.size()); rep.surface_edges = surf.edges; rep.surface_open_edges = surf.open_edges; rep.surface_holes = surf.holes;
    rep.surface_nonmanifold_edges = surf.nonmanifold_edges; rep.surface_flipped_edges = surf.flipped_edges;
    rep.surface_degenerate_faces = surf.degenerate; rep.surface_reoriented_faces = surf.reoriented; rep.surface_duplicate_faces = surf.duplicates;

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

    // --- sommets de la bande : point le plus proche ; signe par nombre d'enroulement (tous les sommets)
    t0 = Clock::now();
    const double band = 1.05 * h;
    std::vector<int8_t> sgn(N, 0), isband(N, 0);  // sgn : +1 dehors, -1 dedans, 0 sur la peau
    std::vector<V3> qpt(N);
    std::vector<float> dist(N, 0.f);
    Winding wn;
    if (prm.use_winding) { wn.build(bvh); wn.beta = prm.winding_beta; }
    std::vector<int8_t> ps_sign(N, 0);            // signe par pseudonormale, pour mesurer les désaccords
    std::vector<uint8_t> ambiguous(N, 0);
    std::vector<float> wval(N, 0.f);               // nombre d'enroulement à la position d'origine du sommet
    parallel_for(N, threads, [&](int64_t b, int64_t e, unsigned) {
        for (int64_t v = b; v < e; ++v) {
            Hit hit = nearest(bvh, pos[v], band);
            if (hit.found) {
                isband[v] = 1; qpt[v] = hit.q; dist[v] = float(std::sqrt(hit.d2));
                ps_sign[v] = hit.d2 == 0 ? 0 : sign_at(surf, hit, pos[v]);
            }
            if (prm.use_winding) {
                double w = wn.query(pos[v], hit.found ? 0.0 : 1.2);   // loin de la peau : dipôle grossier, le signe a de la marge
                wval[v] = float(w);
                sgn[v] = w > 0.5 ? -1 : 1;
                if (hit.found && std::fabs(w - 0.5) < 0.25 && hit.d2 > 4e-4 * h * h) ambiguous[v] = 1;   // un sommet SUR la peau vaut 0,5 sans que rien ne soit cassé
            } else if (hit.found) sgn[v] = ps_sign[v];
        }
    });
    if (prm.use_winding) {
        for (int64_t v = 0; v < N; ++v) if (isband[v]) { if (ps_sign[v] != 0 && ps_sign[v] != sgn[v]) ++rep.sign_conflicts; if (ambiguous[v]) ++rep.ambiguous_points; }
    } else {
        // composantes des sommets hors bande (union-find sur les arêtes du réseau)
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
            int64_t r = find(u); int8_t sv = sgn[v];
            if (sv == 0) continue;
            if (csign[r] == 0) csign[r] = sv; else if (csign[r] != sv) ++rep.sign_conflicts;
        }
        int64_t root0 = find(0);
        if (csign[root0] == 0) csign[root0] = 1;
        for (int64_t v = 0; v < N; ++v) { if (isband[v]) continue; int8_t s = csign[find(v)]; sgn[v] = s ? s : 1; }
    }
    for (int64_t v = 0; v < N; ++v) if (isband[v]) ++rep.band_vertices;
    // --- ramener à la peau les sommets proches
    std::vector<V3> pos0 = pos;                         // positions de référence (réseau régulier)
    std::vector<V3> cut_pos, cut_ref;                   // sommets de coupe : position réelle (sur la peau) et position linéaire
    std::unordered_map<uint64_t, int64_t> cut_id;       // arête (u,v) -> indice du sommet de coupe
    if (!prm.use_cuts) {
        const double snap = prm.snap_alpha * h;
        parallel_for(N, threads, [&](int64_t b, int64_t e, unsigned) {
            for (int64_t v = b; v < e; ++v) if (isband[v] && dist[v] <= snap) { pos[v] = qpt[v]; sgn[v] = 0; }
        });
    } else {
        // 1. arêtes du réseau dont les deux bouts sont de signes opposés ("violées")
        struct VEdge { int64_t u, v; bool is_short; };
        std::vector<VEdge> viol;
        auto consider = [&](int64_t u, int64_t v, bool is_short) {
            if (isband[u] && isband[v] && sgn[u] * sgn[v] < 0) viol.push_back({u, v, is_short});
        };
        for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i) {
            int64_t a = idA(i, j, k);
            if (i < nx) consider(a, idA(i + 1, j, k), false);
            if (j < ny) consider(a, idA(i, j + 1, k), false);
            if (k < nz) consider(a, idA(i, j, k + 1), false);
        }
        for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            int64_t bq = idB(i, j, k);
            for (int dk = 0; dk < 2; ++dk) for (int dj = 0; dj < 2; ++dj) for (int di = 0; di < 2; ++di) consider(bq, idA(i + di, j + dj, k + dk), true);
            if (i + 1 < nx) consider(bq, idB(i + 1, j, k), false);
            if (j + 1 < ny) consider(bq, idB(i, j + 1, k), false);
            if (k + 1 < nz) consider(bq, idB(i, j, k + 1), false);
        }
        // 2. point de coupe de chaque arête : interpolation des distances, ramenée sur la peau
        std::vector<V3> cuts(viol.size()), refs(viol.size());
        parallel_for(int64_t(viol.size()), threads, [&](int64_t b, int64_t e, unsigned) {
            for (int64_t k = b; k < e; ++k) {
                const auto& ed = viol[k];
                double du = dist[ed.u], dv = dist[ed.v];
                double t = (du + dv) > 0 ? du / (du + dv) : 0.5;
                V3 lin = pos[ed.u] + (pos[ed.v] - pos[ed.u]) * t;
                refs[k] = lin;
                Hit hit = nearest(bvh, lin, 2.0 * h);
                cuts[k] = hit.found ? hit.q : lin;
            }
        });
        // 3. ramenage : un sommet à moins de alpha*L d'une coupe de ses arêtes violées va sur la coupe la plus proche
        std::vector<float> ratio(N, 1.0f);
        std::vector<V3> target(N);
        for (size_t k = 0; k < viol.size(); ++k) {
            const auto& ed = viol[k];
            const double L = norm(pos[ed.u] - pos[ed.v]), alpha = ed.is_short ? prm.alpha_short : prm.alpha_long;
            for (int64_t w : {ed.u, ed.v}) {
                double r = norm(pos[w] - cuts[k]) / (alpha * L);
                if (r < ratio[w]) { ratio[w] = float(r); target[w] = cuts[k]; }
            }
        }
        for (int64_t v = 0; v < N; ++v) if (ratio[v] < 1.0f) { pos[v] = target[v]; sgn[v] = 0; }
        // 4. arêtes encore violées : un sommet de coupe chacune
        for (size_t k = 0; k < viol.size(); ++k) {
            const auto& ed = viol[k];
            if (sgn[ed.u] == 0 || sgn[ed.v] == 0) continue;
            cut_id[ekey(int(ed.u), int(ed.v))] = N + int64_t(cut_pos.size());
            cut_pos.push_back(cuts[k]); cut_ref.push_back(refs[k]);
        }
    }
    for (int64_t v = 0; v < N; ++v) if (isband[v] && sgn[v] == 0) ++rep.snapped_vertices;
    rep.t_classify = since(t0);

    // --- tétraèdres : paires de centres voisins + arête de la face partagée
    t0 = Clock::now();
    const double vmin = prm.min_volume_ratio * h * h * h / 12.0;
    struct Local { std::vector<int32_t> tets; int64_t dropped = 0; };
    std::vector<Local> locals(threads);
    auto centroid_inside = [&](const std::array<int64_t, 4>& t) {
        V3 c = (pos[t[0]] + pos[t[1]] + pos[t[2]] + pos[t[3]]) * 0.25;
        if (prm.use_winding) { double m = 0; for (int k = 0; k < 4; ++k) m += wval[t[k]]; m *= 0.25; if (m > 0.8) return true; if (m < 0.2) return false; return wn.query(c) > 0.5; }   // enroulement interpolé, requête exacte seulement si indécis
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

    // ---- version avec découpes : sous-tétraèdres de la partie intérieure de chaque tétraèdre du réseau
    auto P = [&](int64_t v) -> const V3& { return v < N ? pos[v] : cut_pos[size_t(v - N)]; };
    auto R = [&](int64_t v) -> const V3& { return v < N ? pos0[v] : cut_ref[size_t(v - N)]; };
    auto cutv = [&](int64_t u, int64_t v) -> int64_t { auto it = cut_id.find(ekey(int(u), int(v))); return it == cut_id.end() ? -1 : it->second; };
    auto push = [&](Local& L, int64_t p0, int64_t p1, int64_t p2, int64_t p3) {
        if (p0 < 0 || p1 < 0 || p2 < 0 || p3 < 0) { ++L.dropped; return; }
        if (tet_volume(R(p0), R(p1), R(p2), R(p3)) < 0) std::swap(p2, p3);     // orientation de référence positive
        if (!(tet_volume(P(p0), P(p1), P(p2), P(p3)) > vmin)) { ++L.dropped; return; }
        L.tets.push_back(int32_t(p0)); L.tets.push_back(int32_t(p1)); L.tets.push_back(int32_t(p2)); L.tets.push_back(int32_t(p3));
    };
    // prisme (p0 p1 p2 | q0 q1 q2) en trois tétraèdres : les diagonales des faces carrées passent par le plus petit indice
    auto prism = [&](Local& L, int64_t t1[3], int64_t t2[3]) {
        int64_t* a = t1; int64_t* b = t2;
        int64_t mn = a[0]; int where = 0; bool in_second = false;
        for (int k = 0; k < 3; ++k) { if (a[k] < mn) { mn = a[k]; where = k; in_second = false; } if (b[k] < mn) { mn = b[k]; where = k; in_second = true; } }
        if (in_second) std::swap(a, b);
        const int64_t p0 = a[where], p1 = a[(where + 1) % 3], p2 = a[(where + 2) % 3];
        const int64_t q0 = b[where], q1 = b[(where + 1) % 3], q2 = b[(where + 2) % 3];
        const int64_t m12 = std::min({p1, p2, q1, q2});
        if (m12 == p1 || m12 == q2) { push(L, p0, p1, p2, q2); push(L, p0, p1, q1, q2); push(L, p0, q0, q1, q2); }
        else                        { push(L, p0, p2, p1, q1); push(L, p0, p2, q2, q1); push(L, p0, q0, q2, q1); }
    };
    auto emit_cut = [&](Local& L, int64_t a, int64_t b, int64_t c, int64_t d) {
        const int64_t v[4] = {a, b, c, d};
        int64_t ng[4], ps[4], zr[4]; int nn = 0, np = 0, nz = 0;
        for (int k = 0; k < 4; ++k) { if (sgn[v[k]] < 0) ng[nn++] = v[k]; else if (sgn[v[k]] > 0) ps[np++] = v[k]; else zr[nz++] = v[k]; }
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
    };
    auto emit_any = [&](Local& L, int64_t a, int64_t b, int64_t c, int64_t d) { if (prm.use_cuts) emit_cut(L, a, b, c, d); else emit(L, a, b, c, d); };
    parallel_for(int64_t(nz) * ny * nx, threads, [&](int64_t b, int64_t e, unsigned th) {
        Local& L = locals[th];
        for (int64_t q = b; q < e; ++q) {
            int i = q % nx, j = (q / nx) % ny, k = q / (int64_t(nx) * ny);
            int64_t B1 = idB(i, j, k);
            if (i + 1 < nx) {  // face x = i+1, coins en (j,k)
                int64_t B2 = idB(i + 1, j, k);
                int64_t c00 = idA(i + 1, j, k), c10 = idA(i + 1, j + 1, k), c11 = idA(i + 1, j + 1, k + 1), c01 = idA(i + 1, j, k + 1);
                emit_any(L, B1, B2, c00, c10); emit_any(L, B1, B2, c10, c11); emit_any(L, B1, B2, c11, c01); emit_any(L, B1, B2, c01, c00);
            }
            if (j + 1 < ny) {
                int64_t B2 = idB(i, j + 1, k);
                int64_t c00 = idA(i, j + 1, k), c10 = idA(i + 1, j + 1, k), c11 = idA(i + 1, j + 1, k + 1), c01 = idA(i, j + 1, k + 1);
                emit_any(L, B1, B2, c00, c10); emit_any(L, B1, B2, c10, c11); emit_any(L, B1, B2, c11, c01); emit_any(L, B1, B2, c01, c00);
            }
            if (k + 1 < nz) {
                int64_t B2 = idB(i, j, k + 1);
                int64_t c00 = idA(i, j, k + 1), c10 = idA(i + 1, j, k + 1), c11 = idA(i + 1, j + 1, k + 1), c01 = idA(i, j + 1, k + 1);
                emit_any(L, B1, B2, c00, c10); emit_any(L, B1, B2, c10, c11); emit_any(L, B1, B2, c11, c01); emit_any(L, B1, B2, c01, c00);
            }
        }
    });
    std::vector<int32_t> tets;
    for (auto& L : locals) { tets.insert(tets.end(), L.tets.begin(), L.tets.end()); rep.dropped_degenerate += L.dropped; }
    // compactage des sommets
    const int64_t NT = N + int64_t(cut_pos.size());
    std::vector<int32_t> remap(NT, -1);
    int32_t used = 0;
    for (int32_t& t : tets) { if (remap[t] < 0) remap[t] = used++; }
    out.vertices.resize(size_t(used) * 3);
    for (int64_t v = 0; v < NT; ++v) if (remap[v] >= 0) { const V3& q = P(v); size_t o = size_t(remap[v]) * 3; out.vertices[o] = q.x; out.vertices[o + 1] = q.y; out.vertices[o + 2] = q.z; }
    for (int32_t& t : tets) t = remap[t];
    out.tetrahedra = std::move(tets);
    rep.tets = int64_t(out.tetrahedra.size() / 4);
    rep.t_tets = since(t0);

    // --- contrôles
    t0 = Clock::now();
    auto OP = [&](int32_t i) { return V3(out.vertices[3 * size_t(i)], out.vertices[3 * size_t(i) + 1], out.vertices[3 * size_t(i) + 2]); };
    double vol = 0, minvol = 1e300, mind = 180.0;
    struct FK { int32_t a, b, c; bool operator==(const FK& o) const { return a == o.a && b == o.b && c == o.c; } };
    struct FH { size_t operator()(const FK& k) const { return size_t((uint64_t(uint32_t(k.a)) * 0x9E3779B97F4A7C15ull) ^ (uint64_t(uint32_t(k.b)) << 21) ^ (uint64_t(uint32_t(k.c)) * 0xC2B2AE3D27D4EB4Full)); } };
    std::unordered_map<FK, int32_t, FH> faces;  // face triée -> nombre de tétraèdres
    for (int64_t t = 0; t < rep.tets; ++t) {
        const int32_t* q = &out.tetrahedra[4 * size_t(t)];
        V3 p[4] = {OP(q[0]), OP(q[1]), OP(q[2]), OP(q[3])};
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
            V3 a = OP(bfaces[f][0]), bb = OP(bfaces[f][1]), c = OP(bfaces[f][2]);
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
    rep.ambiguous_ratio = rep.band_vertices ? double(rep.ambiguous_points) / double(rep.band_vertices) : 0.0;
    const double judged = prm.use_winding ? rep.ambiguous_ratio : rep.conflict_ratio;
    rep.warning = rep.conflict_ratio > prm.warn_conflict_ratio || rep.ambiguous_ratio > prm.warn_ambiguous_ratio || rep.surface_reoriented_faces > 0;
    if (judged > prm.refuse_conflict_ratio) { rep.status = Status::InconsistentSurface; rep.reason = "surface incohérente : l'intérieur ne se distingue pas de l'extérieur"; return rep; }
    if (rep.tets == 0) { rep.status = Status::InvalidInput; rep.reason = "aucun tétraèdre produit"; }
    else if (rep.nonmanifold_edges > 0) { rep.status = Status::OutOfTolerance; rep.reason = "bord du maillage non fermé ou non manifold"; }
    else if (rep.skin_max_error > prm.skin_tolerance) { rep.status = Status::OutOfTolerance; rep.reason = "écart de peau supérieur à la tolérance"; }
    else rep.status = Status::Ok;
    return rep;
}

}  // namespace tetranet
