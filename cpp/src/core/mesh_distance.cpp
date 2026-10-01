#include "apex/mesh_distance.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>
#include <utility>

namespace apex {

namespace {

constexpr std::int32_t kLeafSize = 4;
constexpr double kInf = std::numeric_limits<double>::infinity();

inline Vec3 sub(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 add_scaled(const Vec3& a, const Vec3& d, double s) { return Vec3{a.x + s * d.x, a.y + s * d.y, a.z + s * d.z}; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double dist2(const Vec3& a, const Vec3& b) {
    const Vec3 d = sub(a, b);
    return dot(d, d);
}

inline double box_dist2(const AABB& b, const Vec3& p) {
    double d = 0.0;
    for (int a = 0; a < 3; ++a) {
        const double v = (&p.x)[a];
        const double lo = (&b.lo.x)[a];
        const double hi = (&b.hi.x)[a];
        if (v < lo) {
            d += (lo - v) * (lo - v);
        } else if (v > hi) {
            d += (v - hi) * (v - hi);
        }
    }
    return d;
}

}  // namespace

Vec3 closest_point_on_triangle(const Vec3& p, const Triangle3& t) {
    const Vec3 ab = sub(t.b, t.a), ac = sub(t.c, t.a), ap = sub(p, t.a);
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return t.a;
    const Vec3 bp = sub(p, t.b);
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return t.b;
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) return add_scaled(t.a, ab, d1 / (d1 - d3));
    const Vec3 cp = sub(p, t.c);
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return t.c;
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) return add_scaled(t.a, ac, d2 / (d2 - d6));
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        return add_scaled(t.b, sub(t.c, t.b), (d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    const double denom = va + vb + vc;
    if (denom == 0.0) {  // triangle degenere : plus proche des trois aretes
        Vec3 best = t.a;
        const Vec3 e[3][2] = {{t.a, t.b}, {t.b, t.c}, {t.c, t.a}};
        for (const auto& s : e) {
            const Vec3 d = sub(s[1], s[0]);
            const double l2 = dot(d, d);
            const double u = l2 > 0.0 ? std::clamp(dot(sub(p, s[0]), d) / l2, 0.0, 1.0) : 0.0;
            const Vec3 q = add_scaled(s[0], d, u);
            if (dist2(p, q) < dist2(p, best)) best = q;
        }
        return best;
    }
    const double v = vb / denom, w = vc / denom;
    return Vec3{t.a.x + ab.x * v + ac.x * w, t.a.y + ab.y * v + ac.y * w, t.a.z + ab.z * v + ac.z * w};
}

MeshDistanceIndex::MeshDistanceIndex(const MeshView& mesh) {
    tris_.reserve(mesh.n_faces);
    boxes_.reserve(mesh.n_faces);
    centroids_.reserve(mesh.n_faces);
    auto vertex = [&](std::int64_t i) {
        if (i < 0 || static_cast<std::size_t>(i) >= mesh.n_vertices) {
            throw std::out_of_range("MeshDistanceIndex : indice de sommet hors limites");
        }
        const double* v = mesh.vertices + 3 * i;
        return Vec3{v[0], v[1], v[2]};
    };
    for (std::size_t f = 0; f < mesh.n_faces; ++f) {
        const Triangle3 t{vertex(mesh.faces[3 * f]), vertex(mesh.faces[3 * f + 1]), vertex(mesh.faces[3 * f + 2])};
        tris_.push_back(t);
        AABB b;
        b.lo = Vec3{std::min({t.a.x, t.b.x, t.c.x}), std::min({t.a.y, t.b.y, t.c.y}), std::min({t.a.z, t.b.z, t.c.z})};
        b.hi = Vec3{std::max({t.a.x, t.b.x, t.c.x}), std::max({t.a.y, t.b.y, t.c.y}), std::max({t.a.z, t.b.z, t.c.z})};
        boxes_.push_back(b);
        centroids_.push_back(Vec3{(t.a.x + t.b.x + t.c.x) / 3.0, (t.a.y + t.b.y + t.c.y) / 3.0,
                                  (t.a.z + t.b.z + t.c.z) / 3.0});
    }
    order_.resize(tris_.size());
    for (std::size_t i = 0; i < order_.size(); ++i) order_[i] = static_cast<std::int32_t>(i);
    nodes_.reserve(tris_.size() * 2 + 1);
    root_ = build(0, static_cast<std::int32_t>(tris_.size()));
}

std::int32_t MeshDistanceIndex::build(std::int32_t lo, std::int32_t hi) {
    if (lo >= hi) return -1;
    const auto index = static_cast<std::int32_t>(nodes_.size());
    nodes_.push_back(Node{});
    AABB box{Vec3{kInf, kInf, kInf}, Vec3{-kInf, -kInf, -kInf}};
    Vec3 clo{kInf, kInf, kInf}, chi{-kInf, -kInf, -kInf};
    for (std::int32_t i = lo; i < hi; ++i) {
        const AABB& b = boxes_[order_[i]];
        const Vec3& c = centroids_[order_[i]];
        for (int a = 0; a < 3; ++a) {
            (&box.lo.x)[a] = std::min((&box.lo.x)[a], (&b.lo.x)[a]);
            (&box.hi.x)[a] = std::max((&box.hi.x)[a], (&b.hi.x)[a]);
            (&clo.x)[a] = std::min((&clo.x)[a], (&c.x)[a]);
            (&chi.x)[a] = std::max((&chi.x)[a], (&c.x)[a]);
        }
    }
    nodes_[index].box = box;
    if (hi - lo <= kLeafSize) {
        nodes_[index].first = lo;
        nodes_[index].count = hi - lo;
        return index;
    }
    int axis = 0;
    if (chi.y - clo.y > chi.x - clo.x) axis = 1;
    if (chi.z - clo.z > (&chi.x)[axis] - (&clo.x)[axis]) axis = 2;
    const std::int32_t mid = lo + (hi - lo) / 2;
    std::nth_element(order_.begin() + lo, order_.begin() + mid, order_.begin() + hi,
                     [this, axis](std::int32_t a, std::int32_t b) {
                         const double ca = (&centroids_[a].x)[axis], cb = (&centroids_[b].x)[axis];
                         return ca != cb ? ca < cb : a < b;
                     });
    const std::int32_t left = build(lo, mid);
    const std::int32_t right = build(mid, hi);
    nodes_[index].left = left;
    nodes_[index].right = right;
    return index;
}

double MeshDistanceIndex::query(const Vec3& p, Vec3* best_point, std::int64_t* best_face) const {
    double best = kInf;
    Vec3 best_q{};
    std::int64_t best_f = -1;
    if (root_ < 0) return kInf;
    // descente en profondeur, enfant le plus proche d'abord, elagage par la
    // distance a la boite (borne inferieure exacte)
    std::vector<std::int32_t> stack;
    stack.reserve(64);
    stack.push_back(root_);
    while (!stack.empty()) {
        const Node& n = nodes_[stack.back()];
        stack.pop_back();
        if (box_dist2(n.box, p) >= best) continue;
        if (n.count > 0) {
            for (std::int32_t i = 0; i < n.count; ++i) {
                const Vec3 q = closest_point_on_triangle(p, tris_[order_[n.first + i]]);
                const double d = dist2(p, q);
                if (d < best) {
                    best = d;
                    best_q = q;
                    best_f = order_[n.first + i];
                }
            }
            continue;
        }
        const double dl = n.left >= 0 ? box_dist2(nodes_[n.left].box, p) : kInf;
        const double dr = n.right >= 0 ? box_dist2(nodes_[n.right].box, p) : kInf;
        const std::int32_t first = dl <= dr ? n.left : n.right;
        const std::int32_t second = dl <= dr ? n.right : n.left;
        if (second >= 0) stack.push_back(second);
        if (first >= 0) stack.push_back(first);
    }
    if (best_point) *best_point = best_q;
    if (best_face) *best_face = best_f;
    return std::sqrt(best);
}

void MeshDistanceIndex::run_range(const double* points, std::size_t begin, std::size_t end, double* out,
                                  double* closest, std::int64_t* face_ids) const {
    for (std::size_t i = begin; i < end; ++i) {
        const Vec3 p{points[3 * i], points[3 * i + 1], points[3 * i + 2]};
        Vec3 q;
        out[i] = query(p, &q, face_ids ? face_ids + i : nullptr);
        if (closest) {
            closest[3 * i] = q.x;
            closest[3 * i + 1] = q.y;
            closest[3 * i + 2] = q.z;
        }
    }
}

DistanceBackend MeshDistanceIndex::resolve_backend(DistanceBackend requested, std::string* reason) {
    if (requested == DistanceBackend::Vulkan) {
        if (reason) *reason = "backend Vulkan non construit (APEX_WITH_VULKAN absent) : repli sur CpuParallel";
        return DistanceBackend::CpuParallel;
    }
    if (reason) reason->clear();
    return requested;
}

void MeshDistanceIndex::distances(const double* points, std::size_t n, double* out, double* closest,
                                  DistanceBackend backend, int threads, std::int64_t* face_ids) const {
    backend = resolve_backend(backend);
    if (backend == DistanceBackend::CpuSerial || n < 2048) {
        run_range(points, 0, n, out, closest, face_ids);
        return;
    }
    unsigned count = threads > 0 ? static_cast<unsigned>(threads) : std::max(1u, std::thread::hardware_concurrency());
    count = static_cast<unsigned>(std::min<std::size_t>(count, (n + 1023) / 1024));
    std::vector<std::thread> pool;
    const std::size_t chunk = (n + count - 1) / count;
    for (unsigned t = 0; t < count; ++t) {
        const std::size_t b = t * chunk, e = std::min(n, b + chunk);
        if (b >= e) break;
        pool.emplace_back([=, this] { run_range(points, b, e, out, closest, face_ids); });
    }
    for (auto& th : pool) th.join();
}

}  // namespace apex
