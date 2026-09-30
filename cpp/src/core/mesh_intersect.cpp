#include "apex/mesh_intersect.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#include "apex/internal/triangle_aabb_bvh.hpp"

namespace apex {
namespace {

constexpr double kEps = 1e-10;  // transcription de _EPS
constexpr double kPi = std::numbers::pi;

inline Vec3 sub(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 add(const Vec3& a, const Vec3& b) { return Vec3{a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 scale(const Vec3& a, double s) { return Vec3{a.x * s, a.y * s, a.z * s}; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

inline const Vec3& vertex(const Triangle3& t, int index) {
    return index == 0 ? t.a : (index == 1 ? t.b : t.c);
}

inline double max_abs_component(const Vec3& v) {
    return std::max({std::abs(v.x), std::abs(v.y), std::abs(v.z)});
}

// Transcription de _interval_on_line : intervalle couvert par `triangle`
// sur la droite d'intersection des deux plans (parametree par `direction`),
// via les points ou le BORD du triangle rencontre l'AUTRE plan.
// `distances[i]` = distance signee du sommet i a l'autre plan.
std::optional<std::pair<double, double>> interval_on_line(const Triangle3& triangle,
                                                           const double distances[3],
                                                           const Vec3& direction, double scale_tol) {
    std::vector<Vec3> points;
    points.reserve(3);
    for (int index = 0; index < 3; ++index) {
        const int following = (index + 1) % 3;
        const double distance_index = distances[index];
        const double distance_following = distances[following];
        if (std::abs(distance_index) <= scale_tol) {
            points.push_back(vertex(triangle, index));
        }
        const bool crosses =
            (distance_index > scale_tol && distance_following < -scale_tol) ||
            (distance_index < -scale_tol && distance_following > scale_tol);
        if (crosses) {
            const double ratio = distance_index / (distance_index - distance_following);
            points.push_back(add(vertex(triangle, index),
                                 scale(sub(vertex(triangle, following), vertex(triangle, index)), ratio)));
        }
    }
    if (points.size() < 2) {
        return std::nullopt;
    }
    double lo = dot(direction, points[0]);
    double hi = lo;
    for (std::size_t i = 1; i < points.size(); ++i) {
        const double projection = dot(direction, points[i]);
        lo = std::min(lo, projection);
        hi = std::max(hi, projection);
    }
    return std::make_pair(lo, hi);
}

// Resolution d'un systeme lineaire 3x3 par la regle de Cramer (cofacteurs) -
// TOLERANCE EXPLICITE avec np.linalg.solve (LAPACK), pas bit-exact : voir
// mesh_intersect.hpp. Les LIGNES de `basis` sont les 3 equations
// (normal_a . x = offset_a, normal_b . x = offset_b, direction . x = 0),
// comme np.linalg.solve(basis, target) ou basis = [normal_a, normal_b, direction].
std::optional<Vec3> solve3x3(const Vec3 rows[3], const double target[3]) {
    const double det = rows[0].x * (rows[1].y * rows[2].z - rows[1].z * rows[2].y) -
                       rows[0].y * (rows[1].x * rows[2].z - rows[1].z * rows[2].x) +
                       rows[0].z * (rows[1].x * rows[2].y - rows[1].y * rows[2].x);
    if (std::abs(det) < 1e-300) {
        return std::nullopt;  // matrice singuliere - transcription de LinAlgError
    }
    // Cramer : remplace la colonne k par `target`, calcule le determinant,
    // divise par det. Les COLONNES de `basis` sont (rows[0][k], rows[1][k], rows[2][k]).
    auto column = [&](int k) {
        return Vec3{(&rows[0].x)[k], (&rows[1].x)[k], (&rows[2].x)[k]};
    };
    const Vec3 col0 = column(0);
    const Vec3 col1 = column(1);
    const Vec3 col2 = column(2);
    const Vec3 t{target[0], target[1], target[2]};

    auto det3 = [](const Vec3& c0, const Vec3& c1, const Vec3& c2) {
        return c0.x * (c1.y * c2.z - c1.z * c2.y) - c1.x * (c0.y * c2.z - c0.z * c2.y) +
               c2.x * (c0.y * c1.z - c0.z * c1.y);
    };

    const double x0 = det3(t, col1, col2) / det;
    const double x1 = det3(col0, t, col2) / det;
    const double x2 = det3(col0, col1, t) / det;
    return Vec3{x0, x1, x2};
}

}  // namespace

std::optional<Segment3> triangle_triangle_segment(const Triangle3& triangle_a, const Triangle3& triangle_b) {
    const Vec3 normal_a = cross(sub(triangle_a.b, triangle_a.a), sub(triangle_a.c, triangle_a.a));
    const Vec3 normal_b = cross(sub(triangle_b.b, triangle_b.a), sub(triangle_b.c, triangle_b.a));
    const double norm_a = norm(normal_a);
    const double norm_b = norm(normal_b);
    if (norm_a < kEps || norm_b < kEps) {
        return std::nullopt;
    }

    const double offset_a = dot(normal_a, triangle_a.a);
    const double offset_b = dot(normal_b, triangle_b.a);

    double distances_b[3] = {dot(normal_a, triangle_b.a) - offset_a, dot(normal_a, triangle_b.b) - offset_a,
                             dot(normal_a, triangle_b.c) - offset_a};
    double distances_a[3] = {dot(normal_b, triangle_a.a) - offset_b, dot(normal_b, triangle_a.b) - offset_b,
                             dot(normal_b, triangle_a.c) - offset_b};

    const double scale_a = kEps * norm_a * (std::max({max_abs_component(triangle_a.a),
                                                       max_abs_component(triangle_a.b),
                                                       max_abs_component(triangle_a.c)}) +
                                            1.0);
    const double scale_b = kEps * norm_b * (std::max({max_abs_component(triangle_b.a),
                                                       max_abs_component(triangle_b.b),
                                                       max_abs_component(triangle_b.c)}) +
                                            1.0);

    bool all_b_above = true, all_b_below = true;
    for (double d : distances_b) {
        all_b_above &= (d > scale_a);
        all_b_below &= (d < -scale_a);
    }
    if (all_b_above || all_b_below) {
        return std::nullopt;
    }
    bool all_a_above = true, all_a_below = true;
    for (double d : distances_a) {
        all_a_above &= (d > scale_b);
        all_a_below &= (d < -scale_b);
    }
    if (all_a_above || all_a_below) {
        return std::nullopt;
    }

    Vec3 direction = cross(normal_a, normal_b);
    const double direction_norm = norm(direction);
    if (direction_norm < kEps) {
        return std::nullopt;
    }
    direction = scale(direction, 1.0 / direction_norm);

    const auto interval_a = interval_on_line(triangle_a, distances_a, direction, scale_b);
    const auto interval_b = interval_on_line(triangle_b, distances_b, direction, scale_a);
    if (!interval_a.has_value() || !interval_b.has_value()) {
        return std::nullopt;
    }

    const double low = std::max(interval_a->first, interval_b->first);
    const double high = std::min(interval_a->second, interval_b->second);
    if (high - low <= kEps) {
        return std::nullopt;
    }

    const Vec3 rows[3] = {normal_a, normal_b, direction};
    const double target[3] = {offset_a, offset_b, 0.0};
    const auto reference = solve3x3(rows, target);
    if (!reference.has_value()) {
        return std::nullopt;
    }

    return Segment3{add(*reference, scale(direction, low)), add(*reference, scale(direction, high))};
}

void generalized_winding_number(const double* points, std::int64_t point_count,
                                const double* triangles, std::int64_t triangle_count, double* out) {
    for (std::int64_t p = 0; p < point_count; ++p) {
        const Vec3 point{points[3 * p], points[3 * p + 1], points[3 * p + 2]};
        double sum = 0.0;
        for (std::int64_t t = 0; t < triangle_count; ++t) {
            const double* tv = triangles + static_cast<std::size_t>(t) * 9;
            const Vec3 a{tv[0] - point.x, tv[1] - point.y, tv[2] - point.z};
            const Vec3 b{tv[3] - point.x, tv[4] - point.y, tv[5] - point.z};
            const Vec3 c{tv[6] - point.x, tv[7] - point.y, tv[8] - point.z};

            const double norm_a = norm(a);
            const double norm_b = norm(b);
            const double norm_c = norm(c);

            const double numerator = dot(a, cross(b, c));
            const double denominator = norm_a * norm_b * norm_c + dot(a, b) * norm_c +
                                       dot(b, c) * norm_a + dot(c, a) * norm_b;
            sum += std::atan2(numerator, denominator);
        }
        out[p] = sum / (2.0 * kPi);
    }
}

std::vector<IntersectionSegment> find_intersection_segments(const double* triangles,
                                                             const std::int64_t* faces,
                                                             std::int64_t triangle_count,
                                                             double min_length) {
    std::vector<IntersectionSegment> result;
    if (triangle_count <= 0) {
        return result;
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

    const internal::TriangleAabbBvh bvh(boxes);
    std::vector<std::int32_t> candidates;

    for (std::int64_t i = 0; i < triangle_count; ++i) {
        candidates.clear();
        bvh.gather(boxes[static_cast<std::size_t>(i)], candidates);
        // Tri par indice croissant : NE reproduit PAS l'ordre de
        // `tree.intersection(query)` cote Python (rtree/libspatialindex
        // rend un ordre spatial interne, verifie non trie par indice sur
        // aucune requete testee), mais rapproche legerement les deux
        // resultats sans cout. L'ORDRE des segments par face influence la
        // triangulation de sortie de `_subdivide_face` (voir
        // mesh_intersect.hpp pour la limite documentee, mesuree sur
        // turkish.stl) - un ordre parfaitement fidele exigerait de
        // reimplementer le parcours interne de libspatialindex, hors de
        // portee raisonnable pour un ecart deja mineur (~0.001% de volume).
        std::sort(candidates.begin(), candidates.end());

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
            // >=2 sommets partages => arete commune, pas une auto-intersection -
            // transcription de `np.intersect1d(faces[i], faces[j]).size >= 2`.
            int shared = 0;
            shared += (fi0 == fj0 || fi0 == fj1 || fi0 == fj2) ? 1 : 0;
            shared += (fi1 == fj0 || fi1 == fj1 || fi1 == fj2) ? 1 : 0;
            shared += (fi2 == fj0 || fi2 == fj1 || fi2 == fj2) ? 1 : 0;
            if (shared >= 2) {
                continue;
            }

            const auto segment = triangle_triangle_segment(tris[static_cast<std::size_t>(i)],
                                                            tris[static_cast<std::size_t>(j)]);
            if (!segment.has_value()) {
                continue;
            }
            if (norm(sub(segment->end, segment->start)) < min_length) {
                continue;
            }
            result.push_back(IntersectionSegment{static_cast<std::int32_t>(i),
                                                  static_cast<std::int32_t>(j), *segment});
        }
    }

    return result;
}

}  // namespace apex
