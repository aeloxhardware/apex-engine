// BVH interne sur AABB de triangles - broad-phase generique, PARTAGE entre
// mesh_repair.cpp (SAT triangle-triangle, pas 8) et mesh_intersect.cpp
// (segments d'intersection + winding number, pas 9).
//
// Structure PUREMENT INTERNE, sans contrat d'equivalence : contrairement au
// BVH de rayons (bvh.hpp, pas 6), aucun des deux usages n'exige que ce
// broad-phase retourne EXACTEMENT les memes candidats qu'un rtree cote
// Python - la narrow-phase de chaque appelant (SAT, ou calcul de segment
// d'intersection) est deterministe et donne la meme reponse pour tout
// sur-ensemble conservateur de candidats. Voir mesh_repair.hpp pour le
// raisonnement complet sur pourquoi ce cas differe du pas 6.

#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "apex/types.hpp"

namespace apex::internal {

inline bool triangle_aabb_boxes_overlap(const AABB& a, const AABB& b) {
    return a.lo.x <= b.hi.x && a.hi.x >= b.lo.x && a.lo.y <= b.hi.y && a.hi.y >= b.lo.y &&
           a.lo.z <= b.hi.z && a.hi.z >= b.lo.z;
}

class TriangleAabbBvh {
public:
    explicit TriangleAabbBvh(const std::vector<AABB>& boxes) : boxes_(boxes) {
        const auto n = static_cast<std::int64_t>(boxes_.size());
        order_.resize(static_cast<std::size_t>(n));
        centroids_.resize(static_cast<std::size_t>(n));
        for (std::int64_t i = 0; i < n; ++i) {
            order_[static_cast<std::size_t>(i)] = static_cast<std::int32_t>(i);
            const AABB& b = boxes_[static_cast<std::size_t>(i)];
            centroids_[static_cast<std::size_t>(i)] =
                Vec3{(b.lo.x + b.hi.x) * 0.5, (b.lo.y + b.hi.y) * 0.5, (b.lo.z + b.hi.z) * 0.5};
        }
        nodes_.reserve(static_cast<std::size_t>(n) * 2);
        root_ = build(0, n);
    }

    // Candidats dont l'AABB recoupe `query`, ajoutes a `out` (non efface).
    void gather(const AABB& query, std::vector<std::int32_t>& out) const { gather(root_, query, out); }

private:
    // Feuille des que le noeud tient dans ce nombre de triangles - le test
    // d'AABB coute alors plus cher que les quelques tests qu'il economise.
    static constexpr std::int64_t kLeafSize = 4;

    struct Node {
        AABB box;
        std::int32_t left = -1;
        std::int32_t right = -1;
        std::int32_t first = 0;
        std::int32_t count = 0;
    };

    std::int32_t build(std::int64_t lo, std::int64_t hi) {
        if (lo >= hi) {
            return -1;
        }
        const auto node_index = static_cast<std::int32_t>(nodes_.size());
        nodes_.push_back(Node{});

        AABB box;
        box.lo = Vec3{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::infinity()};
        box.hi = Vec3{-std::numeric_limits<double>::infinity(),
                      -std::numeric_limits<double>::infinity(),
                      -std::numeric_limits<double>::infinity()};
        Vec3 centroid_lo = box.lo;
        Vec3 centroid_hi = box.hi;
        for (std::int64_t i = lo; i < hi; ++i) {
            const AABB& tb = boxes_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])];
            const Vec3& c = centroids_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])];
            for (int a = 0; a < 3; ++a) {
                (&box.lo.x)[a] = std::min((&box.lo.x)[a], (&tb.lo.x)[a]);
                (&box.hi.x)[a] = std::max((&box.hi.x)[a], (&tb.hi.x)[a]);
                (&centroid_lo.x)[a] = std::min((&centroid_lo.x)[a], (&c.x)[a]);
                (&centroid_hi.x)[a] = std::max((&centroid_hi.x)[a], (&c.x)[a]);
            }
        }
        nodes_[static_cast<std::size_t>(node_index)].box = box;

        const std::int64_t span = hi - lo;
        if (span <= kLeafSize) {
            Node& leaf = nodes_[static_cast<std::size_t>(node_index)];
            leaf.first = static_cast<std::int32_t>(lo);
            leaf.count = static_cast<std::int32_t>(span);
            return node_index;
        }

        int axis = 0;
        double extent = centroid_hi.x - centroid_lo.x;
        if (centroid_hi.y - centroid_lo.y > extent) {
            extent = centroid_hi.y - centroid_lo.y;
            axis = 1;
        }
        if (centroid_hi.z - centroid_lo.z > extent) {
            axis = 2;
        }

        const std::int64_t mid = lo + span / 2;
        std::nth_element(order_.begin() + lo, order_.begin() + mid, order_.begin() + hi,
                         [this, axis](std::int32_t a, std::int32_t b) {
                             const double ca = (&centroids_[static_cast<std::size_t>(a)].x)[axis];
                             const double cb = (&centroids_[static_cast<std::size_t>(b)].x)[axis];
                             if (ca != cb) {
                                 return ca < cb;
                             }
                             return a < b;
                         });

        const std::int32_t left = build(lo, mid);
        const std::int32_t right = build(mid, hi);
        Node& n = nodes_[static_cast<std::size_t>(node_index)];
        n.left = left;
        n.right = right;
        n.count = 0;
        return node_index;
    }

    void gather(std::int32_t node_index, const AABB& query, std::vector<std::int32_t>& out) const {
        if (node_index < 0) {
            return;
        }
        const Node& n = nodes_[static_cast<std::size_t>(node_index)];
        if (!triangle_aabb_boxes_overlap(n.box, query)) {
            return;
        }
        if (n.count > 0) {
            for (std::int32_t i = 0; i < n.count; ++i) {
                const std::int32_t t = order_[static_cast<std::size_t>(n.first + i)];
                if (triangle_aabb_boxes_overlap(boxes_[static_cast<std::size_t>(t)], query)) {
                    out.push_back(t);
                }
            }
            return;
        }
        gather(n.left, query, out);
        gather(n.right, query, out);
    }

    const std::vector<AABB>& boxes_;
    std::vector<Vec3> centroids_;
    std::vector<std::int32_t> order_;
    std::vector<Node> nodes_;
    std::int32_t root_ = -1;
};

}  // namespace apex::internal
