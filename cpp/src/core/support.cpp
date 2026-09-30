#include "apex/support.hpp"

#include <cmath>

namespace apex {
namespace {

constexpr double kParallelEps = 1e-9;

inline double dot3(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm3(const Vec3& a) { return std::sqrt(dot3(a, a)); }
inline Vec3 scale3(const Vec3& a, double s) { return Vec3{a.x * s, a.y * s, a.z * s}; }
inline Vec3 sub3(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }

inline Vec3 cross3(const Vec3& a, const Vec3& b) {
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

// Base locale (u, v, d) - MEME convention que `_perpendicular_basis` cote
// Python (reference fixe [1,0,0], repli [0,1,0] au-dela de 0.9). En changer
// donnerait une base differente, donc d'autres coordonnees locales.
struct LocalFrame {
    Vec3 u, v, d;
};

LocalFrame build_frame(const Vec3& direction_up) {
    Vec3 d = direction_up;
    const double n = norm3(d);
    d = (n > 1e-12) ? scale3(d, 1.0 / n) : Vec3{0.0, 0.0, 1.0};

    Vec3 reference{1.0, 0.0, 0.0};
    if (std::fabs(dot3(reference, d)) > 0.9) {
        reference = Vec3{0.0, 1.0, 0.0};
    }
    Vec3 u = sub3(reference, scale3(d, dot3(reference, d)));
    u = scale3(u, 1.0 / norm3(u));
    return LocalFrame{u, cross3(d, u), d};
}

// Le point du segment [a,b] le plus proche de `p`, exprime dans le repere
// local, puis les deux ecarts qui decident du support. Rend `true` des que la
// condition est remplie - c'est tout l'interet de ce portage.
inline bool candidate_supports(const LocalFrame& frame, double pu, double pv, double pd,
                               const Vec3& a, const Vec3& b, double horizontal_radius,
                               double vertical_height) {
    // Coordonnees locales des deux extremites.
    const double au = dot3(a, frame.u), av = dot3(a, frame.v), ad = dot3(a, frame.d);
    const double bu = dot3(b, frame.u), bv = dot3(b, frame.v), bd = dot3(b, frame.d);

    const double su = bu - au, sv = bv - av, sd = bd - ad;
    const double length_sq = su * su + sv * sv + sd * sd;

    // Segment degenere : `t = 0` cote Python (np.where(degenerate, 0.0, t)).
    double t = 0.0;
    if (length_sq >= 1e-18) {
        t = ((pu - au) * su + (pv - av) * sv + (pd - ad) * sd) / length_sq;
        if (t < 0.0) {
            t = 0.0;
        } else if (t > 1.0) {
            t = 1.0;
        }
    }

    const double delta_up = pd - (ad + t * sd);
    if (delta_up <= kParallelEps || delta_up > vertical_height) {
        return false;
    }
    const double dh_u = pu - (au + t * su);
    const double dh_v = pv - (av + t * sv);
    return std::hypot(dh_u, dh_v) <= horizontal_radius;
}

}  // namespace

bool point_has_support_along_direction(const Vec3& position, const Vec3& direction_up,
                                       const double* starts, const double* ends,
                                       std::int64_t candidate_count, double horizontal_radius,
                                       double vertical_height) {
    if (candidate_count <= 0) {
        return false;
    }

    const LocalFrame frame = build_frame(direction_up);
    const double pu = dot3(position, frame.u);
    const double pv = dot3(position, frame.v);
    const double pd = dot3(position, frame.d);

    for (std::int64_t i = 0; i < candidate_count; ++i) {
        const Vec3 a{starts[3 * i + 0], starts[3 * i + 1], starts[3 * i + 2]};
        const Vec3 b{ends[3 * i + 0], ends[3 * i + 1], ends[3 * i + 2]};
        if (candidate_supports(frame, pu, pv, pd, a, b, horizontal_radius, vertical_height)) {
            return true;  // ARRET ANTICIPE - le Python evaluait les ~690 candidats
        }
    }
    return false;
}

bool point_has_support_in_any_direction(const Vec3& position, const double* directions,
                                        std::int64_t direction_count, const double* starts,
                                        const double* ends, std::int64_t candidate_count,
                                        double horizontal_radius, double vertical_height) {
    if (candidate_count <= 0) {
        return false;
    }
    for (std::int64_t k = 0; k < direction_count; ++k) {
        const Vec3 direction{directions[3 * k + 0], directions[3 * k + 1], directions[3 * k + 2]};
        if (point_has_support_along_direction(position, direction, starts, ends, candidate_count,
                                              horizontal_radius, vertical_height)) {
            return true;
        }
    }
    return false;
}

}  // namespace apex
