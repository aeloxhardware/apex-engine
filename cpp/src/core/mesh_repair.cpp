#include "apex/mesh_repair.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "apex/internal/triangle_aabb_bvh.hpp"

namespace apex {
namespace {

// Tolerance de separation - transcription exacte de `_COPLANAR_EPS`.
constexpr double kCoplanarEps = 1e-9;
// Norme d'axe en dessous de laquelle il est ignore (aretes paralleles ou
// face degeneree) - transcription de `norm < 1e-12`.
constexpr double kAxisNormEps = 1e-12;

inline Vec3 sub(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

using apex::internal::TriangleAabbBvh;

}  // namespace

bool triangles_cross(const Triangle3& t0, const Triangle3& t1) {
    const Vec3 edges_a[3] = {sub(t0.b, t0.a), sub(t0.c, t0.b), sub(t0.a, t0.c)};
    const Vec3 edges_b[3] = {sub(t1.b, t1.a), sub(t1.c, t1.b), sub(t1.a, t1.c)};

    Vec3 axes[11];
    int axis_count = 0;
    axes[axis_count++] = cross(edges_a[0], edges_a[1]);
    axes[axis_count++] = cross(edges_b[0], edges_b[1]);
    for (const Vec3& ea : edges_a) {
        for (const Vec3& eb : edges_b) {
            axes[axis_count++] = cross(ea, eb);
        }
    }

    const Vec3 va[3] = {t0.a, t0.b, t0.c};
    const Vec3 vb[3] = {t1.a, t1.b, t1.c};

    for (int i = 0; i < axis_count; ++i) {
        const double n = norm(axes[i]);
        if (n < kAxisNormEps) {
            continue;
        }
        const Vec3 axis{axes[i].x / n, axes[i].y / n, axes[i].z / n};

        double pa_min = dot(va[0], axis), pa_max = pa_min;
        for (int k = 1; k < 3; ++k) {
            const double p = dot(va[k], axis);
            pa_min = std::min(pa_min, p);
            pa_max = std::max(pa_max, p);
        }
        double pb_min = dot(vb[0], axis), pb_max = pb_min;
        for (int k = 1; k < 3; ++k) {
            const double p = dot(vb[k], axis);
            pb_min = std::min(pb_min, p);
            pb_max = std::max(pb_max, p);
        }

        if (pa_min >= pb_max - kCoplanarEps || pb_min >= pa_max - kCoplanarEps) {
            return false;
        }
    }
    return true;
}

std::vector<std::int32_t> find_self_intersecting_faces(const double* triangles,
                                                        const std::int64_t* faces,
                                                        std::int64_t triangle_count) {
    std::vector<std::int32_t> bad;
    if (triangle_count <= 0) {
        return bad;
    }

    std::vector<Triangle3> tris(static_cast<std::size_t>(triangle_count));
    std::vector<AABB> boxes(static_cast<std::size_t>(triangle_count));
    for (std::int64_t i = 0; i < triangle_count; ++i) {
        const double* v = triangles + static_cast<std::size_t>(i) * 9;
        Triangle3& t = tris[static_cast<std::size_t>(i)];
        t.a = Vec3{v[0], v[1], v[2]};
        t.b = Vec3{v[3], v[4], v[5]};
        t.c = Vec3{v[6], v[7], v[8]};

        AABB& box = boxes[static_cast<std::size_t>(i)];
        box.lo.x = std::min({t.a.x, t.b.x, t.c.x});
        box.lo.y = std::min({t.a.y, t.b.y, t.c.y});
        box.lo.z = std::min({t.a.z, t.b.z, t.c.z});
        box.hi.x = std::max({t.a.x, t.b.x, t.c.x});
        box.hi.y = std::max({t.a.y, t.b.y, t.c.y});
        box.hi.z = std::max({t.a.z, t.b.z, t.c.z});
    }

    const TriangleAabbBvh bvh(boxes);
    std::vector<bool> is_bad(static_cast<std::size_t>(triangle_count), false);
    std::vector<std::int32_t> candidates;

    for (std::int64_t i = 0; i < triangle_count; ++i) {
        candidates.clear();
        bvh.gather(boxes[static_cast<std::size_t>(i)], candidates);

        const std::int64_t fi0 = faces[i * 3 + 0];
        const std::int64_t fi1 = faces[i * 3 + 1];
        const std::int64_t fi2 = faces[i * 3 + 2];

        for (const std::int32_t j32 : candidates) {
            const std::int64_t j = j32;
            if (j <= i) {
                continue;
            }
            const std::int64_t fj0 = faces[j * 3 + 0];
            const std::int64_t fj1 = faces[j * 3 + 1];
            const std::int64_t fj2 = faces[j * 3 + 2];
            // Sommet partage => contact par construction, pas une
            // intersection - transcription de `np.intersect1d(...).size > 0`.
            if (fi0 == fj0 || fi0 == fj1 || fi0 == fj2 || fi1 == fj0 || fi1 == fj1 ||
                fi1 == fj2 || fi2 == fj0 || fi2 == fj1 || fi2 == fj2) {
                continue;
            }
            if (triangles_cross(tris[static_cast<std::size_t>(i)], tris[static_cast<std::size_t>(j)])) {
                is_bad[static_cast<std::size_t>(i)] = true;
                is_bad[static_cast<std::size_t>(j)] = true;
            }
        }
    }

    for (std::int64_t i = 0; i < triangle_count; ++i) {
        if (is_bad[static_cast<std::size_t>(i)]) {
            bad.push_back(static_cast<std::int32_t>(i));
        }
    }
    return bad;
}

}  // namespace apex
