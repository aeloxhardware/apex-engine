// TetraNet — surface d'entrée : soudure, topologie, orientation ; arbre de recherche, enroulement, point le plus proche.
#include "internal.hpp"

#include <cstring>
#include <numeric>

namespace tetranet::detail {

namespace {

// Soude les sommets confondus (coordonnées identiques au bit près, -0 == +0) ; `map` : sommet d'entrée -> soudé.
void weld_vertices(const SurfaceMesh& in, Surface& s, std::vector<int>& map) {
    const size_t n = in.vertices.size() / 3;
    struct Key { uint64_t a, b, c; bool operator==(const Key& o) const { return a == o.a && b == o.b && c == o.c; } };
    // Les doubles issus d'un STL (float32) ont leurs 29 bits de poids faible nuls : sans mélange final, les masques de
    // buckets ne voient presque que des collisions (soudure de 678 000 sommets : 149 s). splitmix64 sur chaque mot.
    struct H { size_t operator()(const Key& k) const { return size_t(mix64(mix64(mix64(k.a) ^ k.b) ^ k.c)); } };
    std::unordered_map<Key, int, H> weld;
    weld.reserve(n);
    map.assign(n, 0);
    for (size_t i = 0; i < n; ++i) {
        double c[3] = {in.vertices[3 * i] + 0.0, in.vertices[3 * i + 1] + 0.0, in.vertices[3 * i + 2] + 0.0};  // +0.0 : -0 == +0
        Key key;
        std::memcpy(&key.a, &c[0], 8); std::memcpy(&key.b, &c[1], 8); std::memcpy(&key.c, &c[2], 8);
        auto it = weld.find(key);
        if (it == weld.end()) { it = weld.emplace(key, int(s.p.size())).first; s.p.push_back({c[0], c[1], c[2]}); }
        map[i] = it->second;
    }
}

// Triangles non dégénérés, normales de face, volume signé.
void keep_triangles(const SurfaceMesh& in, const std::vector<int>& map, Surface& s) {
    const size_t m = in.faces.size() / 3;
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
}

// Topologie des arêtes : (nombre de faces, sens parcourus) ; trous, arêtes ouvertes / non manifold / retournées.
void measure_topology(Surface& s) {
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

// Propagation de l'orientation sur les arêtes à deux faces (BFS), puis volume positif par composante.
void fix_orientations(Surface& s) {
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

// Pseudonormales d'arête (somme des normales) et de sommet (pondérées par l'angle), pour le signe local.
void pseudonormals(Surface& s) {
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
}

}  // namespace

bool build_surface(const SurfaceMesh& in, Surface& s, bool fix_orientation) {
    if (in.vertices.size() / 3 == 0 || in.faces.size() / 3 == 0) return false;
    std::vector<int> map;
    weld_vertices(in, s, map);
    keep_triangles(in, map, s);
    if (s.tri.empty()) return false;
    measure_topology(s);
    if (fix_orientation) fix_orientations(s);
    pseudonormals(s);
    return true;
}

void Bvh::build(const Surface& surf) {
    s = &surf;
    order.resize(surf.tri.size());
    std::iota(order.begin(), order.end(), 0);
    nodes.reserve(surf.tri.size());
    rec(0, int(order.size()));
}

int Bvh::rec(int b, int e) {
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

void Winding::build(const Bvh& b) {
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

double Winding::solid(const V3& a, const V3& b, const V3& c) {   // Van Oosterom - Strakhov
    double la = norm(a), lb = norm(b), lc = norm(c);
    double den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
    return 2.0 * std::atan2(dot(a, cross(b, c)), den);
}

double Winding::query(const V3& p, double beta_override) const {
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

std::unique_ptr<FeatureEdges> make_features(const Surface& s, const Params& prm, double radius) {
    if (!(prm.feature_angle_deg > 0.0) || !(prm.feature_snap_alpha > 0.0)) return nullptr;
    std::vector<std::array<double, 3>> pts(s.p.size()), nrm(s.fn.size());
    for (size_t i = 0; i < s.p.size(); ++i) pts[i] = {s.p[i].x, s.p[i].y, s.p[i].z};
    for (size_t i = 0; i < s.fn.size(); ++i) nrm[i] = {s.fn[i].x, s.fn[i].y, s.fn[i].z};
    auto features = std::make_unique<FeatureEdges>(pts, s.tri, nrm, prm.feature_angle_deg, radius);
    if (features->empty()) return nullptr;
    return features;
}

V3 snap_to_feature(const FeatureEdges* features, const V3& q, double radius) {
    if (!features) return q;
    std::array<double, 3> out;
    return features->snap({q.x, q.y, q.z}, radius, out) ? V3(out[0], out[1], out[2]) : q;
}

}  // namespace tetranet::detail
