// TetraNet — points de la surface d'entrée, peau extérieure, bord interrogeable, angle dièdre.
#include "cover.hpp"

namespace tetranet::detail {

std::vector<Sample> surface_samples(const Surface& s, double spacing) {
    std::vector<Sample> out;
    out.reserve(s.p.size() + 2 * s.tri.size());
    std::vector<int32_t> owner(s.p.size(), -1);
    for (size_t f = 0; f < s.tri.size(); ++f) for (int v : s.tri[f]) if (owner[v] < 0) owner[v] = int32_t(f);
    for (size_t v = 0; v < s.p.size(); ++v) if (owner[v] >= 0) out.push_back({s.p[v], owner[v]});
    const double g = 1.32471795724474602596, a1 = 1.0 / g, a2 = 1.0 / (g * g);   // suite R2 (nombre plastique)
    for (size_t f = 0; f < s.tri.size(); ++f) {
        const auto& t = s.tri[f];
        const V3 &a = s.p[t[0]], &b = s.p[t[1]], &cc = s.p[t[2]];
        const double area = 0.5 * norm(cross(b - a, cc - a));
        const int n = std::clamp(int(std::ceil(2.0 * area / (spacing * spacing))), 1, 4096);
        for (int i = 0; i < n; ++i) {
            double u = std::fmod(0.5 + a1 * (i + 1), 1.0), w = std::fmod(0.5 + a2 * (i + 1), 1.0);
            if (u + w > 1.0) { u = 1.0 - u; w = 1.0 - w; }
            out.push_back({a + (b - a) * u + (cc - a) * w, int32_t(f)});
        }
        for (int k = 0; k < 3; ++k) {
            const int i0 = t[k], i1 = t[(k + 1) % 3];
            if (i0 > i1) continue;
            const V3 e = s.p[i1] - s.p[i0];
            const int m = int(norm(e) / spacing);
            for (int j = 1; j <= m; ++j) out.push_back({s.p[i0] + e * (double(j) / (m + 1)), int32_t(f)});
        }
    }
    return out;
}

OuterSkin::OuterSkin(const Surface& s, const Winding* wn, double eps)
    : s_(s), wn_(wn), eps_(eps), cls_(new std::atomic<int8_t>[s.tri.size()]) {
    for (size_t f = 0; f < s.tri.size(); ++f) cls_[f].store(0, std::memory_order_relaxed);
}

bool OuterSkin::outer(int32_t f) const {
    if (!wn_) return true;
    int8_t known = cls_[f].load(std::memory_order_relaxed);
    if (known == 0) {
        const auto& t = s_.tri[f];
        const V3 centre = (s_.p[t[0]] + s_.p[t[1]] + s_.p[t[2]]) * (1.0 / 3);
        known = wn_->query(centre + s_.fn[f] * eps_) < 0.5 ? 1 : 2;
        cls_[f].store(known, std::memory_order_relaxed);
    }
    return known == 1;
}

SkinIndex::SkinIndex(const std::vector<V3>& points, const std::vector<std::array<int32_t, 3>>& faces) {
    s.p = points;
    s.tri.reserve(faces.size());
    for (const auto& f : faces) s.tri.push_back({f[0], f[1], f[2]});
    if (!s.tri.empty()) bvh.build(s);
}

double min_dihedral(const V3 p[4]) {
    double mind = 180.0;
    static const int pr[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (auto& e : pr) {
        int cc = -1, d = -1;
        for (int k = 0; k < 4; ++k) if (k != e[0] && k != e[1]) (cc < 0 ? cc : d) = k;
        V3 ev = unit(p[e[1]] - p[e[0]]);
        V3 u = p[cc] - p[e[0]]; u = u - ev * dot(u, ev);
        V3 w = p[d] - p[e[0]]; w = w - ev * dot(w, ev);
        double den = norm(u) * norm(w);
        if (den > 0) mind = std::min(mind, std::acos(std::max(-1.0, std::min(1.0, dot(u, w) / den))) * 180.0 / M_PI);
    }
    return mind;
}

}  // namespace tetranet::detail
