// Voir apex/past_material.hpp.

#include "apex/past_material.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <thread>

namespace apex {

PastMaterialIndex::PastMaterialIndex(const double* points, const std::int64_t* owners, std::int64_t count)
    : count_(count > 0 ? count : 0) {
    if (count_ == 0) return;
    points_.assign(points, points + count_ * 3);
    owners_.assign(owners, owners + count_);
    idx_.resize(static_cast<std::size_t>(count_));
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(count_); ++i) idx_[i] = i;
    nodes_.reserve(static_cast<std::size_t>(count_));
    root_ = build(0, count_);
}

std::int32_t PastMaterialIndex::build(std::int64_t lo, std::int64_t hi) {
    if (lo >= hi) return -1;
    double mn[3], mx[3];
    std::int64_t min_owner = std::numeric_limits<std::int64_t>::max();
    for (int a = 0; a < 3; ++a) {
        mn[a] = std::numeric_limits<double>::infinity();
        mx[a] = -std::numeric_limits<double>::infinity();
    }
    for (std::int64_t i = lo; i < hi; ++i) {
        const double* p = &points_[static_cast<std::size_t>(idx_[i]) * 3];
        for (int a = 0; a < 3; ++a) {
            mn[a] = std::min(mn[a], p[a]);
            mx[a] = std::max(mx[a], p[a]);
        }
        min_owner = std::min(min_owner, owners_[idx_[i]]);
    }
    int axis = 0;
    double span = mx[0] - mn[0];
    if (mx[1] - mn[1] > span) { span = mx[1] - mn[1]; axis = 1; }
    if (mx[2] - mn[2] > span) axis = 2;
    const std::int64_t mid = lo + (hi - lo) / 2;
    std::nth_element(idx_.begin() + lo, idx_.begin() + mid, idx_.begin() + hi,
                     [this, axis](std::int32_t a, std::int32_t b) {
                         const double va = points_[static_cast<std::size_t>(a) * 3 + axis];
                         const double vb = points_[static_cast<std::size_t>(b) * 3 + axis];
                         return va != vb ? va < vb : a < b;
                     });
    const auto node = static_cast<std::int32_t>(nodes_.size());
    nodes_.push_back(Node{});
    const std::int32_t point = idx_[mid];
    const std::int32_t left = build(lo, mid);
    const std::int32_t right = build(mid + 1, hi);
    Node& n = nodes_[static_cast<std::size_t>(node)];
    n.left = left; n.right = right; n.point = point; n.axis = static_cast<std::int8_t>(axis);
    n.min_owner = min_owner;
    for (int a = 0; a < 3; ++a) { n.lo[a] = mn[a]; n.hi[a] = mx[a]; }
    return node;
}

void PastMaterialIndex::search(std::int32_t ni, const double* q, std::int64_t limit, double& best_sq) const {
    if (ni < 0) return;
    const Node& n = nodes_[static_cast<std::size_t>(ni)];
    if (n.min_owner >= limit) return;                 // sous-arbre entierement futur
    double box_sq = 0.0;                              // distance a la boite
    for (int a = 0; a < 3; ++a) {
        const double d = q[a] < n.lo[a] ? n.lo[a] - q[a] : (q[a] > n.hi[a] ? q[a] - n.hi[a] : 0.0);
        box_sq += d * d;
    }
    if (box_sq > best_sq) return;
    if (owners_[n.point] < limit) {
        const double* p = &points_[static_cast<std::size_t>(n.point) * 3];
        const double dx = q[0] - p[0], dy = q[1] - p[1], dz = q[2] - p[2];
        const double sq = dx * dx + dy * dy + dz * dz;
        if (sq < best_sq) best_sq = sq;
    }
    const double* p = &points_[static_cast<std::size_t>(n.point) * 3];
    const bool go_left = q[n.axis] < p[n.axis];
    search(go_left ? n.left : n.right, q, limit, best_sq);
    search(go_left ? n.right : n.left, q, limit, best_sq);
}

double PastMaterialIndex::min_past_distance(const double* queries, std::int64_t query_count,
                                            std::int64_t limit, double reach) const {
    if (root_ < 0) return std::numeric_limits<double>::infinity();
    // seuil initial : reach^2 (un point exactement a `reach` reste admis,
    // comme query_ball_point ; best < reach est teste par l'appelant)
    double best_sq = reach * reach;
    bool found = false;
    for (std::int64_t i = 0; i < query_count; ++i) {
        double b = best_sq;
        search(root_, queries + i * 3, limit, b);
        if (b < best_sq) {      // `search` ne retient que sq < seuil, strictement
            best_sq = b;
            found = true;
        }
    }
    return found ? std::sqrt(best_sq) : std::numeric_limits<double>::infinity();
}

bool PastMaterialIndex::search_any(std::int32_t ni, const double* q, std::int64_t limit, double radius_sq) const {
    if (ni < 0) return false;
    const Node& n = nodes_[static_cast<std::size_t>(ni)];
    if (n.min_owner >= limit) return false;
    double box_sq = 0.0;
    for (int a = 0; a < 3; ++a) {
        const double d = q[a] < n.lo[a] ? n.lo[a] - q[a] : (q[a] > n.hi[a] ? q[a] - n.hi[a] : 0.0);
        box_sq += d * d;
    }
    if (box_sq > radius_sq) return false;
    if (owners_[n.point] < limit) {
        const double* p = &points_[static_cast<std::size_t>(n.point) * 3];
        const double dx = q[0] - p[0], dy = q[1] - p[1], dz = q[2] - p[2];
        if (dx * dx + dy * dy + dz * dz <= radius_sq) return true;
    }
    return search_any(n.left, q, limit, radius_sq) || search_any(n.right, q, limit, radius_sq);
}

bool PastMaterialIndex::any_past_within(const double* queries, std::int64_t query_count, std::int64_t limit,
                                        double radius) const {
    if (root_ < 0) return false;
    const double radius_sq = radius * radius;
    for (std::int64_t i = 0; i < query_count; ++i)
        if (search_any(root_, queries + i * 3, limit, radius_sq)) return true;
    return false;
}

void PastMaterialIndex::min_past_distance_batch(const double* queries, std::int64_t pose_count,
                                                std::int64_t samples, const std::int64_t* limits, double reach,
                                                double* out, int threads) const {
    auto work = [&](std::int64_t begin, std::int64_t end) {
        for (std::int64_t p = begin; p < end; ++p)
            out[p] = min_past_distance(queries + p * samples * 3, samples, limits[p], reach);
    };
    unsigned count = threads > 0 ? static_cast<unsigned>(threads) : std::max(1u, std::thread::hardware_concurrency());
    count = static_cast<unsigned>(std::min<std::int64_t>(count, std::max<std::int64_t>(1, pose_count / 16)));
    if (count <= 1) {
        work(0, pose_count);
        return;
    }
    std::vector<std::thread> pool;
    const std::int64_t chunk = (pose_count + count - 1) / count;
    for (unsigned t = 0; t < count; ++t) {
        const std::int64_t b = t * chunk, e = std::min(pose_count, b + chunk);
        if (b >= e) break;
        pool.emplace_back(work, b, e);
    }
    for (auto& th : pool) th.join();
}

}  // namespace apex
