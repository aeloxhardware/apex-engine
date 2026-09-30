#include "apex/kdtree.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace apex {

KDTree::KDTree(const double* points, std::int64_t point_count)
    : point_count_(point_count > 0 ? point_count : 0) {
    if (point_count_ == 0) {
        return;
    }
    points_.assign(points, points + static_cast<std::size_t>(point_count_) * 3);

    indices_.resize(static_cast<std::size_t>(point_count_));
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(point_count_); ++i) {
        indices_[static_cast<std::size_t>(i)] = i;
    }
    nodes_.reserve(static_cast<std::size_t>(point_count_));
    root_ = build(0, point_count_);
}

std::int32_t KDTree::build(std::int64_t lo, std::int64_t hi) {
    if (lo >= hi) {
        return -1;
    }

    // Axe de coupe : celui de plus grande etendue sur l'intervalle courant -
    // choix classique, meilleur equilibre qu'un axe cyclique fixe.
    double min_v[3] = {points_[static_cast<std::size_t>(indices_[static_cast<std::size_t>(lo)]) * 3 + 0],
                       points_[static_cast<std::size_t>(indices_[static_cast<std::size_t>(lo)]) * 3 + 1],
                       points_[static_cast<std::size_t>(indices_[static_cast<std::size_t>(lo)]) * 3 + 2]};
    double max_v[3] = {min_v[0], min_v[1], min_v[2]};
    for (std::int64_t i = lo + 1; i < hi; ++i) {
        const double* p = point_at(indices_[static_cast<std::size_t>(i)]);
        for (int a = 0; a < 3; ++a) {
            min_v[a] = std::min(min_v[a], p[a]);
            max_v[a] = std::max(max_v[a], p[a]);
        }
    }
    int axis = 0;
    double span = max_v[0] - min_v[0];
    if (max_v[1] - min_v[1] > span) {
        span = max_v[1] - min_v[1];
        axis = 1;
    }
    if (max_v[2] - min_v[2] > span) {
        axis = 2;
    }

    const std::int64_t mid = lo + (hi - lo) / 2;
    std::nth_element(
        indices_.begin() + lo, indices_.begin() + mid, indices_.begin() + hi,
        [this, axis](std::int32_t a, std::int32_t b) {
            const double va = point_at(a)[axis];
            const double vb = point_at(b)[axis];
            if (va != vb) {
                return va < vb;
            }
            return a < b;  // ordre total stable : evite toute dependance a l'ordre d'entree
        });

    const auto node_index = static_cast<std::int32_t>(nodes_.size());
    nodes_.push_back(Node{});
    // Node& n = nodes_[node_index] serait invalide apres les push_back des
    // sous-arbres (realloc) : on renseigne les champs apres construction.
    const std::int32_t point_index = indices_[static_cast<std::size_t>(mid)];
    const std::int32_t left = build(lo, mid);
    const std::int32_t right = build(mid + 1, hi);

    Node& n = nodes_[static_cast<std::size_t>(node_index)];
    n.axis = static_cast<std::int8_t>(axis);
    n.point = point_index;
    n.left = left;
    n.right = right;
    return node_index;
}

void KDTree::nearest(std::int32_t node_index, double qx, double qy, double qz,
                     std::int64_t& best_index, double& best_sqdist) const {
    if (node_index < 0) {
        return;
    }
    const Node& n = nodes_[static_cast<std::size_t>(node_index)];
    const double* p = point_at(n.point);
    const double dx = qx - p[0];
    const double dy = qy - p[1];
    const double dz = qz - p[2];
    const double sqdist = dx * dx + dy * dy + dz * dz;
    if (sqdist < best_sqdist || (sqdist == best_sqdist && n.point < best_index)) {
        best_sqdist = sqdist;
        best_index = n.point;
    }

    const double q[3] = {qx, qy, qz};
    const double diff = q[n.axis] - p[n.axis];
    const std::int32_t near_child = diff < 0.0 ? n.left : n.right;
    const std::int32_t far_child = diff < 0.0 ? n.right : n.left;

    nearest(near_child, qx, qy, qz, best_index, best_sqdist);
    // Cote oppose visite seulement si le plan de coupe est plus proche que le
    // meilleur candidat courant (`<=` : une egalite peut cacher un point de
    // plus petit indice a la meme distance).
    if (diff * diff <= best_sqdist) {
        nearest(far_child, qx, qy, qz, best_index, best_sqdist);
    }
}

void KDTree::within(std::int32_t node_index, double qx, double qy, double qz, double radius,
                    double radius_sq, std::vector<std::int32_t>& out) const {
    if (node_index < 0) {
        return;
    }
    const Node& n = nodes_[static_cast<std::size_t>(node_index)];
    const double* p = point_at(n.point);
    const double dx = qx - p[0];
    const double dy = qy - p[1];
    const double dz = qz - p[2];
    if (dx * dx + dy * dy + dz * dz <= radius_sq) {
        out.push_back(n.point);
    }

    const double q[3] = {qx, qy, qz};
    const double diff = q[n.axis] - p[n.axis];
    if (diff < 0.0) {
        within(n.left, qx, qy, qz, radius, radius_sq, out);
        if (-diff <= radius) {
            within(n.right, qx, qy, qz, radius, radius_sq, out);
        }
    } else {
        within(n.right, qx, qy, qz, radius, radius_sq, out);
        if (diff <= radius) {
            within(n.left, qx, qy, qz, radius, radius_sq, out);
        }
    }
}

void KDTree::query_nn(const double* queries, std::int64_t query_count, std::int64_t* out_indices,
                      double* out_distances) const {
    for (std::int64_t i = 0; i < query_count; ++i) {
        if (point_count_ == 0) {
            // Cas jamais atteint par le pipeline (tous les arbres sont
            // construits sur des ensembles non vides, cf. le garde-fou de
            // `_blend_tet_directions`) - protege seulement contre un OOB.
            out_indices[i] = -1;
            out_distances[i] = std::numeric_limits<double>::infinity();
            continue;
        }
        const double qx = queries[3 * i + 0];
        const double qy = queries[3 * i + 1];
        const double qz = queries[3 * i + 2];

        std::int64_t best_index = point_count_;  // sentinelle > tout indice valide
        double best_sqdist = std::numeric_limits<double>::infinity();
        nearest(root_, qx, qy, qz, best_index, best_sqdist);

        out_indices[i] = best_index;
        // Distance euclidienne : meme forme close que scipy (p=2), somme
        // accumulee dans l'ordre x, y, z pour un resultat identique au bit.
        const double* p = point_at(static_cast<std::int32_t>(best_index));
        const double dx = qx - p[0];
        const double dy = qy - p[1];
        const double dz = qz - p[2];
        out_distances[i] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
}

void KDTree::query_radius(const Vec3& point, double radius, std::vector<std::int32_t>& out) const {
    if (root_ < 0 || radius < 0.0) {
        return;
    }
    within(root_, point.x, point.y, point.z, radius, radius * radius, out);
}

}  // namespace apex
