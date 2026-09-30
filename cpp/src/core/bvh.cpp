#include "apex/bvh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace apex {
namespace {

// Feuille des que le noeud tient dans ce nombre de triangles. Un seuil > 1
// evite une descente profonde pour un gain nul : le test d'AABB coute alors
// plus cher que les quelques tests de triangle qu'il economise.
constexpr std::int64_t kLeafSize = 4;

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

// Recouvrement d'AABB, bornes INCLUSES - c'est la semantique de
// `rtree.intersection`, que cette descente remplace.
inline bool boxes_overlap(const AABB& a, const AABB& b) {
    return a.lo.x <= b.hi.x && a.hi.x >= b.lo.x && a.lo.y <= b.hi.y && a.hi.y >= b.lo.y &&
           a.lo.z <= b.hi.z && a.hi.z >= b.lo.z;
}

// Cle de dedoublonnage : (location quantifiee a 8 decimales, index du rayon),
// transcription de `grouping.unique_rows` sur la pile
// `column_stack((locations, index_ray))`.
struct HitKey {
    std::int64_t x, y, z, ray;
    bool operator==(const HitKey& o) const {
        return x == o.x && y == o.y && z == o.z && ray == o.ray;
    }
};

struct HitKeyHash {
    std::size_t operator()(const HitKey& k) const noexcept {
        std::size_t h = static_cast<std::size_t>(k.x) * 0x9E3779B97F4A7C15ULL;
        h ^= static_cast<std::size_t>(k.y) + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
        h ^= static_cast<std::size_t>(k.z) + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
        h ^= static_cast<std::size_t>(k.ray) + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

// `trimesh.grouping.float_to_int` : multiplication par 10^digits puis
// troncature vers zero (`astype(np.int64)`), jamais un arrondi.
inline std::int64_t quantize(double value) {
    return static_cast<std::int64_t>(value * kMergeDigits);
}

}  // namespace

TriangleBVH::TriangleBVH(const double* triangles, const double* normals,
                         std::int64_t triangle_count)
    : triangle_count_(triangle_count > 0 ? triangle_count : 0) {
    if (triangle_count_ == 0) {
        return;
    }
    triangles_.assign(triangles, triangles + static_cast<std::size_t>(triangle_count_) * 9);
    normals_.assign(normals, normals + static_cast<std::size_t>(triangle_count_) * 3);

    boxes_.resize(static_cast<std::size_t>(triangle_count_));
    centroids_.resize(static_cast<std::size_t>(triangle_count_));
    order_.resize(static_cast<std::size_t>(triangle_count_));

    // AABB EXACTES, sans marge : ce sont celles que `triangles.bounds_tree`
    // indexe cote trimesh, et le jeu de candidats doit etre le meme.
    mesh_bounds_.lo = Vec3{std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity()};
    mesh_bounds_.hi = Vec3{-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()};

    for (std::int64_t t = 0; t < triangle_count_; ++t) {
        const double* v = triangles_.data() + static_cast<std::size_t>(t) * 9;
        AABB& box = boxes_[static_cast<std::size_t>(t)];
        for (int a = 0; a < 3; ++a) {
            const double c0 = v[a];
            const double c1 = v[3 + a];
            const double c2 = v[6 + a];
            const double lo = std::min(c0, std::min(c1, c2));
            const double hi = std::max(c0, std::max(c1, c2));
            (&box.lo.x)[a] = lo;
            (&box.hi.x)[a] = hi;
            (&centroids_[static_cast<std::size_t>(t)].x)[a] = (c0 + c1 + c2) / 3.0;
            (&mesh_bounds_.lo.x)[a] = std::min((&mesh_bounds_.lo.x)[a], lo);
            (&mesh_bounds_.hi.x)[a] = std::max((&mesh_bounds_.hi.x)[a], hi);
        }
        order_[static_cast<std::size_t>(t)] = static_cast<std::int32_t>(t);
    }

    nodes_.reserve(static_cast<std::size_t>(triangle_count_) * 2);
    root_ = build(0, triangle_count_);
}

std::int32_t TriangleBVH::build(std::int64_t lo, std::int64_t hi) {
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

    // Axe le plus etendu au sens des CENTROIDES (et non des AABB) : c'est ce
    // qui separe reellement les triangles, une AABB large pouvant deborder
    // des deux cotes du plan de coupe.
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
                         return a < b;  // ordre total : construction deterministe
                     });

    const std::int32_t left = build(lo, mid);
    const std::int32_t right = build(mid, hi);
    Node& n = nodes_[static_cast<std::size_t>(node_index)];
    n.left = left;
    n.right = right;
    n.count = 0;
    return node_index;
}

void TriangleBVH::gather(std::int32_t node_index, const AABB& query,
                         std::vector<std::int32_t>& out) const {
    if (node_index < 0) {
        return;
    }
    const Node& n = nodes_[static_cast<std::size_t>(node_index)];
    if (!boxes_overlap(n.box, query)) {
        return;
    }
    if (n.count > 0) {
        // Feuille : le test par noeud est conservateur, chaque triangle est
        // re-teste individuellement pour ne pas elargir le jeu de candidats.
        for (std::int32_t i = 0; i < n.count; ++i) {
            const std::int32_t t = order_[static_cast<std::size_t>(n.first + i)];
            if (boxes_overlap(boxes_[static_cast<std::size_t>(t)], query)) {
                out.push_back(t);
            }
        }
        return;
    }
    gather(n.left, query, out);
    gather(n.right, query, out);
}

AABB TriangleBVH::ray_bounds(const Vec3& origin, const Vec3& direction) const {
    // Transcription de `ray_triangle.ray_bounds`, un rayon a la fois.
    const double d[3] = {direction.x, direction.y, direction.z};
    const double o[3] = {origin.x, origin.y, origin.z};

    // Axe dominant du rayon (`np.abs(ray_directions).argmax`) - `argmax` rend
    // le PREMIER maximum en cas d'egalite, d'où les comparaisons strictes.
    int axis = 0;
    double best = std::abs(d[0]);
    if (std::abs(d[1]) > best) {
        best = std::abs(d[1]);
        axis = 1;
    }
    if (std::abs(d[2]) > best) {
        axis = 2;
    }

    const double axis_dir = d[axis];
    const double axis_ori = o[axis];
    const double axis_bound[2] = {(&mesh_bounds_.lo.x)[axis], (&mesh_bounds_.hi.x)[axis]};

    // `t = zeros; t[nonzero] = (bound - ori)/dir` puis `t[t < buffer] = buffer`.
    // Une direction nulle sur l'axe dominant (donc un vecteur nul) laisse les
    // deux t a zero, ensuite releves a buffer_dist - comportement conserve.
    double t[2] = {0.0, 0.0};
    if (axis_dir != 0.0) {
        t[0] = (axis_bound[0] - axis_ori) / axis_dir;
        t[1] = (axis_bound[1] - axis_ori) / axis_dir;
    }
    for (int i = 0; i < 2; ++i) {
        if (t[i] < kRayBufferDist) {
            t[i] = kRayBufferDist;
        }
    }

    AABB box;
    for (int a = 0; a < 3; ++a) {
        const double on_a = d[a] * t[0] + o[a];
        const double on_b = d[a] * t[1] + o[a];
        (&box.lo.x)[a] = std::min(on_a, on_b) - kRayBufferDist;
        (&box.hi.x)[a] = std::max(on_a, on_b) + kRayBufferDist;
    }
    return box;
}

std::vector<RayHit> TriangleBVH::intersects_location(const double* origins,
                                                     const double* directions,
                                                     std::int64_t ray_count) const {
    std::vector<RayHit> hits;
    if (triangle_count_ == 0 || ray_count <= 0) {
        return hits;
    }

    std::vector<std::int32_t> candidates;
    std::unordered_set<HitKey, HitKeyHash> seen;

    for (std::int64_t r = 0; r < ray_count; ++r) {
        const Vec3 origin{origins[3 * r], origins[3 * r + 1], origins[3 * r + 2]};
        const Vec3 direction{directions[3 * r], directions[3 * r + 1], directions[3 * r + 2]};

        candidates.clear();
        gather(root_, ray_bounds(origin, direction), candidates);
        if (candidates.empty()) {
            continue;
        }
        seen.clear();

        for (const std::int32_t t : candidates) {
            const double* v = triangles_.data() + static_cast<std::size_t>(t) * 9;
            const double* nrm = normals_.data() + static_cast<std::size_t>(t) * 3;

            // --- 1. intersections.planes_lines -------------------------------
            // plane_origin = premier sommet du triangle, comme trimesh.
            const Vec3 plane_origin{v[0], v[1], v[2]};
            const Vec3 plane_normal{nrm[0], nrm[1], nrm[2]};
            const Vec3 origin_vector{plane_origin.x - origin.x, plane_origin.y - origin.y,
                                     plane_origin.z - origin.z};

            const double projection_ori = dot(origin_vector, plane_normal);
            const double projection_dir = dot(direction, plane_normal);
            if (!(std::abs(projection_dir) > kPlaneParallelEps)) {
                continue;
            }
            const double distance = projection_ori / projection_dir;
            // `on_plane = line_direction * distance + line_origin` - meme
            // ordre (multiplication puis addition) que le numpy de trimesh.
            const Vec3 location{direction.x * distance + origin.x,
                                direction.y * distance + origin.y,
                                direction.z * distance + origin.z};

            // --- 2. triangles.points_to_barycentric, methode "cramer" --------
            const Vec3 edge0{v[3] - v[0], v[4] - v[1], v[5] - v[2]};
            const Vec3 edge1{v[6] - v[0], v[7] - v[1], v[8] - v[2]};
            const Vec3 w{location.x - v[0], location.y - v[1], location.z - v[2]};

            const double dot00 = dot(edge0, edge0);
            const double dot01 = dot(edge0, edge1);
            const double dot02 = dot(edge0, w);
            const double dot11 = dot(edge1, edge1);
            const double dot12 = dot(edge1, w);

            const double inverse_denominator = 1.0 / (dot00 * dot11 - dot01 * dot01);
            const double b2 = (dot00 * dot12 - dot01 * dot02) * inverse_denominator;
            const double b1 = (dot11 * dot02 - dot01 * dot12) * inverse_denominator;
            const double b0 = 1.0 - b1 - b2;

            // --- 3. appartenance au triangle ---------------------------------
            const bool inside = b0 > -kBarycentricTol && b1 > -kBarycentricTol &&
                                b2 > -kBarycentricTol && b0 < 1.0 + kBarycentricTol &&
                                b1 < 1.0 + kBarycentricTol && b2 < 1.0 + kBarycentricTol;
            if (!inside) {
                continue;
            }

            // --- 4. filtre "vers l'avant" ------------------------------------
            // Produit scalaire avec la direction telle quelle (non
            // normalisee), comme `util.diagonal_dot` cote trimesh.
            const Vec3 to_hit{location.x - origin.x, location.y - origin.y,
                              location.z - origin.z};
            if (!(dot(to_hit, direction) > kForwardTol)) {
                continue;
            }

            // --- 5. dedoublonnage (location, index_ray) ----------------------
            const HitKey key{quantize(location.x), quantize(location.y), quantize(location.z), r};
            if (!seen.insert(key).second) {
                continue;
            }

            hits.push_back(RayHit{r, t, location});
        }
    }
    return hits;
}

}  // namespace apex
