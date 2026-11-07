#include "tetranet/features.hpp"

#include <algorithm>
#include <cmath>

namespace tetranet {

namespace {

using P3 = std::array<double, 3>;

double dot3(const P3& a, const P3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
P3 sub3(const P3& a, const P3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
P3 lerp3(const P3& a, const P3& b, double t) { return {a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t}; }

uint64_t edge_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (uint64_t(a) << 32) | uint32_t(b);
}

struct EdgeFaces { int f[2] = {-1, -1}; int count = 0; };

}  // namespace

int64_t FeatureEdges::cell_key(int64_t i, int64_t j, int64_t k) const {
    constexpr int64_t kBias = int64_t(1) << 20;
    return (((i + kBias) << 42) | ((j + kBias) << 21)) | (k + kBias);
}

FeatureEdges::FeatureEdges(const std::vector<P3>& points, const std::vector<std::array<int, 3>>& triangles,
                           const std::vector<P3>& normals, double angle_deg, double piece) {
    cell_ = 2.0 * piece;
    const double cos_limit = std::cos(angle_deg * 3.14159265358979323846 / 180.0);

    std::unordered_map<uint64_t, EdgeFaces> edges;
    edges.reserve(triangles.size() * 2);
    for (size_t t = 0; t < triangles.size(); ++t) {
        for (int e = 0; e < 3; ++e) {
            EdgeFaces& entry = edges[edge_key(triangles[t][e], triangles[t][(e + 1) % 3])];
            if (entry.count < 2) entry.f[entry.count] = int(t);
            ++entry.count;
        }
    }

    std::vector<uint64_t> keys;                                           // ordre trié : indépendant de la table de hachage
    for (const auto& [key, entry] : edges) {
        if (entry.count != 2) continue;                                   // bord ouvert ou arête non-variété : ignorée
        if (dot3(normals[entry.f[0]], normals[entry.f[1]]) >= cos_limit) continue;   // angle sous le seuil : arête lisse
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());
    for (uint64_t key : keys) {
        const P3& a = points[int(key >> 32)];
        const P3& b = points[int(key & 0xffffffffu)];
        const P3 d = sub3(b, a);
        const double len = std::sqrt(dot3(d, d));
        const int pieces = std::max(1, int(std::ceil(len / piece)));
        for (int k = 0; k < pieces; ++k) {
            Segment seg{lerp3(a, b, double(k) / pieces), lerp3(a, b, double(k + 1) / pieces)};
            segments_.push_back(seg);
        }
        ++edges_;
    }

    for (size_t s = 0; s < segments_.size(); ++s) {
        const P3 mid = lerp3(segments_[s].a, segments_[s].b, 0.5);
        grid_[cell_key(int64_t(std::floor(mid[0] / cell_)), int64_t(std::floor(mid[1] / cell_)),
                       int64_t(std::floor(mid[2] / cell_)))].push_back(int32_t(s));
    }
}

bool FeatureEdges::snap(const P3& p, double radius, P3& out) const {
    if (segments_.empty()) return false;
    const int64_t ci = int64_t(std::floor(p[0] / cell_)), cj = int64_t(std::floor(p[1] / cell_)),
                  ck = int64_t(std::floor(p[2] / cell_));
    double best = radius * radius;
    bool found = false;
    for (int64_t di = -1; di <= 1; ++di) for (int64_t dj = -1; dj <= 1; ++dj) for (int64_t dk = -1; dk <= 1; ++dk) {
        const auto it = grid_.find(cell_key(ci + di, cj + dj, ck + dk));
        if (it == grid_.end()) continue;
        for (int32_t s : it->second) {
            const Segment& seg = segments_[s];
            const P3 ab = sub3(seg.b, seg.a);
            const double denom = dot3(ab, ab);
            double t = denom > 0 ? dot3(sub3(p, seg.a), ab) / denom : 0.0;
            t = std::min(1.0, std::max(0.0, t));
            const P3 q = lerp3(seg.a, seg.b, t);
            const P3 r = sub3(p, q);
            const double d2 = dot3(r, r);
            if (d2 < best) { best = d2; out = q; found = true; }
        }
    }
    return found;
}

}  // namespace tetranet
