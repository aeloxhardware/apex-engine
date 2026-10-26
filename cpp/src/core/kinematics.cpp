#include "apex/kinematics.hpp"

#include <cmath>
#include <optional>

namespace apex {
namespace {

constexpr double kPi = 3.14159265358979323846;

inline double deg2rad(double d) { return d * kPi / 180.0; }
inline double rad2deg(double r) { return r * 180.0 / kPi; }

inline double dot3(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm3(const Vec3& a) { return std::sqrt(dot3(a, a)); }

inline Vec3 cross3(const Vec3& a, const Vec3& b) {
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline Vec3 scale3(const Vec3& a, double s) { return Vec3{a.x * s, a.y * s, a.z * s}; }
inline Vec3 sub3(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }

inline Vec3 normalized(const Vec3& a) {
    const double n = norm3(a);
    return n > 0.0 ? scale3(a, 1.0 / n) : a;
}

// Transcription de `_rodrigues_rotate_single` - meme ordre d'operations.
Vec3 rodrigues_rotate_single(const Vec3& v, const Vec3& axis, double angle_rad) {
    const double cos_a = std::cos(angle_rad);
    const double sin_a = std::sin(angle_rad);
    const double d = v.x * axis.x + v.y * axis.y + v.z * axis.z;
    const double cross_x = axis.y * v.z - axis.z * v.y;
    const double cross_y = axis.z * v.x - axis.x * v.z;
    const double cross_z = axis.x * v.y - axis.y * v.x;
    const double factor = d * (1.0 - cos_a);
    return Vec3{v.x * cos_a + cross_x * sin_a + axis.x * factor,
                v.y * cos_a + cross_y * sin_a + axis.y * factor,
                v.z * cos_a + cross_z * sin_a + axis.z * factor};
}

// Transcription de `_angle_between` (vecteurs unitaires, resultat en degres).
double angle_between(const Vec3& a, const Vec3& b) {
    double cos_angle = a.x * b.x + a.y * b.y + a.z * b.z;
    if (cos_angle > 1.0) {
        cos_angle = 1.0;
    } else if (cos_angle < -1.0) {
        cos_angle = -1.0;
    }
    return rad2deg(std::acos(cos_angle));
}

// Transcription de `_perpendicular_basis` - MEME convention de reference
// (fixe [1,0,0], repli [0,1,0] au-dela de 0.9) : en changer produirait une
// base differente, donc des angles signes differents.
void perpendicular_basis(const Vec3& axis, Vec3& u, Vec3& v) {
    Vec3 reference{1.0, 0.0, 0.0};
    if (std::fabs(dot3(reference, axis)) > 0.9) {
        reference = Vec3{0.0, 1.0, 0.0};
    }
    u = sub3(reference, scale3(axis, dot3(reference, axis)));
    u = scale3(u, 1.0 / norm3(u));
    v = cross3(axis, u);
}

// Transcription de `_signed_angle_in_plane` - nullopt = angle indetermine.
std::optional<double> signed_angle_in_plane(const Vec3& vector_from, const Vec3& vector_to,
                                            const Vec3& axis) {
    Vec3 u, v;
    perpendicular_basis(axis, u, v);
    const Vec3 from_perp = sub3(vector_from, scale3(axis, dot3(vector_from, axis)));
    const Vec3 to_perp = sub3(vector_to, scale3(axis, dot3(vector_to, axis)));
    if (norm3(from_perp) < kParallelEps || norm3(to_perp) < kParallelEps) {
        return std::nullopt;
    }
    const double fu = dot3(from_perp, u), fv = dot3(from_perp, v);
    const double tu = dot3(to_perp, u), tv = dot3(to_perp, v);
    return rad2deg(std::atan2(tv * fu - tu * fv, tu * fu + tv * fv));
}

// Transcription de `_clamp_into_range`, y compris le clamp CIRCULAIRE (une
// plage bornee clampe a la borne la plus proche par enroulement, pas au sens
// lineaire naif).
double clamp_into_range(double angle_deg, double lo, double hi) {
    if (hi - lo >= 360.0) {
        // std::fmod garde le signe du dividende, contrairement au `%` de
        // Python : il faut remonter explicitement dans [0,360).
        double m = std::fmod(angle_deg - lo, 360.0);
        if (m < 0.0) {
            m += 360.0;
        }
        return lo + m;
    }
    if (lo <= angle_deg && angle_deg <= hi) {
        return angle_deg;
    }
    const double dist_lo = std::fmin(std::fabs(angle_deg - lo), 360.0 - std::fabs(angle_deg - lo));
    const double dist_hi = std::fmin(std::fabs(angle_deg - hi), 360.0 - std::fabs(angle_deg - hi));
    return dist_lo <= dist_hi ? lo : hi;
}

// Transcription de `_forward_direction` (composition la plus interne d'abord).
Vec3 forward_direction(const std::vector<AxisSpec>& axes, const std::vector<double>& angles_deg,
                       const Vec3& reference_direction) {
    Vec3 direction = reference_direction;
    for (std::size_t k = axes.size(); k-- > 0;) {
        const Vec3 axis = normalized(axes[k].axis);
        direction = rodrigues_rotate_single(direction, axis, deg2rad(angles_deg[k]));
    }
    return direction;
}

// Transcription de `_resolve_single_axis_analytic`.
double resolve_single_axis_analytic(const AxisSpec& spec, const Vec3& reference, const Vec3& target,
                                    double fallback_angle_deg) {
    if (spec.hi <= spec.lo) {
        return spec.lo;
    }
    const Vec3 axis = normalized(spec.axis);
    const std::optional<double> delta = signed_angle_in_plane(reference, target, axis);
    if (!delta.has_value()) {
        return clamp_into_range(fallback_angle_deg, spec.lo, spec.hi);
    }
    return clamp_into_range(*delta, spec.lo, spec.hi);
}

// Transcription de `_resolve_axis_pair_analytic` - nullopt = aucune racine
// bracketee, l'appelant doit basculer sur le repli grille.
std::optional<std::pair<double, double>> resolve_axis_pair_analytic(
    const AxisSpec& outer_spec, const AxisSpec& inner_spec, const Vec3& reference_effective,
    const Vec3& target, double fallback_outer, double fallback_inner) {
    const Vec3 outer_axis = normalized(outer_spec.axis);
    const Vec3 inner_axis = normalized(inner_spec.axis);

    const double outer_lo = outer_spec.lo, outer_hi = outer_spec.hi;
    const double inner_lo = inner_spec.lo, inner_hi = inner_spec.hi;

    const bool outer_locked = outer_hi <= outer_lo;
    const bool inner_locked = inner_hi <= inner_lo;

    if (outer_locked && inner_locked) {
        return std::make_pair(outer_lo, inner_lo);
    }

    if (inner_locked) {
        const Vec3 rotated_reference =
            rodrigues_rotate_single(reference_effective, inner_axis, deg2rad(inner_lo));
        const double angle =
            resolve_single_axis_analytic(outer_spec, rotated_reference, target, fallback_outer);
        return std::make_pair(angle, inner_lo);
    }

    if (outer_locked) {
        const Vec3 unwound_target =
            rodrigues_rotate_single(target, outer_axis, deg2rad(-outer_lo));
        const double angle =
            resolve_single_axis_analytic(inner_spec, reference_effective, unwound_target,
                                         fallback_inner);
        return std::make_pair(outer_lo, angle);
    }

    if (std::fabs(dot3(outer_axis, inner_axis)) > 1.0 - kAxesParallelEps) {
        // Axes paralleles : un seul DOF effectif, plage combinee.
        AxisSpec combined_spec;
        combined_spec.axis = inner_axis;
        combined_spec.lo = inner_lo + outer_lo;
        combined_spec.hi = inner_hi + outer_hi;
        const double combined = resolve_single_axis_analytic(
            combined_spec, reference_effective, target, fallback_outer + fallback_inner);
        const double inner_angle = clamp_into_range(combined - outer_lo, inner_lo, inner_hi);
        return std::make_pair(outer_lo, inner_angle);
    }

    const double target_angle_to_outer = angle_between(target, outer_axis);

    auto achieved_angle_to_outer = [&](double inner_angle_deg) {
        const Vec3 rotated =
            rodrigues_rotate_single(reference_effective, inner_axis, deg2rad(inner_angle_deg));
        return angle_between(rotated, outer_axis);
    };

    // Bracketing sur kBracketingSamples echantillons. Le Python vectorise
    // (np.linspace + rotation batch + arccos) ; ici une boucle scalaire suffit
    // et donne la meme suite de valeurs. `np.linspace` inclut les deux bornes,
    // avec le dernier point pose EXACTEMENT a `inner_hi` (pas lo + n*step, qui
    // accumulerait une derive d'arrondi) : on reproduit ce comportement.
    constexpr int samples = kBracketingSamples > 2 ? kBracketingSamples : 2;
    double sample_values[samples];
    double residuals[samples];
    for (int i = 0; i < samples; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(samples - 1);
        sample_values[i] = (i == samples - 1) ? inner_hi : inner_lo + t * (inner_hi - inner_lo);
        residuals[i] = achieved_angle_to_outer(sample_values[i]) - target_angle_to_outer;
    }

    auto sign_of = [](double x) { return x > 0.0 ? 1 : (x < 0.0 ? -1 : 0); };

    bool has_bracket = false;
    double lo_t = 0.0, hi_t = 0.0, lo_r = 0.0;
    for (int i = 0; i + 1 < samples; ++i) {
        if (sign_of(residuals[i + 1]) - sign_of(residuals[i]) != 0) {
            lo_t = sample_values[i];
            hi_t = sample_values[i + 1];
            lo_r = residuals[i];
            has_bracket = true;
            break;
        }
    }
    if (!has_bracket) {
        for (int i = 0; i < samples; ++i) {
            if (residuals[i] == 0.0) {
                lo_t = hi_t = sample_values[i];
                lo_r = residuals[i];
                has_bracket = true;
                break;
            }
        }
    }
    if (!has_bracket) {
        return std::nullopt;
    }

    double inner_angle = lo_t;
    for (int it = 0; it < kBisectionMaxIterations; ++it) {
        const double mid_t = 0.5 * (lo_t + hi_t);
        const double mid_r = achieved_angle_to_outer(mid_t) - target_angle_to_outer;
        inner_angle = mid_t;
        if (std::fabs(mid_r) < kBisectionToleranceDeg) {
            break;
        }
        if ((lo_r > 0.0) == (mid_r > 0.0)) {
            lo_t = mid_t;
            lo_r = mid_r;
        } else {
            hi_t = mid_t;
        }
    }

    const Vec3 intermediate =
        rodrigues_rotate_single(reference_effective, inner_axis, deg2rad(inner_angle));
    const std::optional<double> delta = signed_angle_in_plane(intermediate, target, outer_axis);
    const double outer_angle = delta.has_value()
                                   ? clamp_into_range(*delta, outer_lo, outer_hi)
                                   : clamp_into_range(fallback_outer, outer_lo, outer_hi);

    return std::make_pair(outer_angle, inner_angle);
}

// Transcription de `_axis_chain_world_lines` + `_detect_singular_axes`.
void detect_singular_axes(const Vec3& position, const Vec3& platform_origin,
                          const std::vector<AxisSpec>& axes,
                          const std::vector<double>& angles_deg, double epsilon,
                          std::vector<int>& out_singular) {
    out_singular.clear();
    if (axes.empty()) {
        return;
    }

    const Vec3 relative = sub3(position, platform_origin);
    std::vector<Vec3> accumulated_axes;
    std::vector<double> accumulated_angles;
    accumulated_axes.reserve(axes.size());
    accumulated_angles.reserve(axes.size());

    for (std::size_t axis_index = 0; axis_index < axes.size(); ++axis_index) {
        const Vec3 axis = normalized(axes[axis_index].axis);
        Vec3 world_axis = axis;
        // Rotations des DOF PRECEDENTS, appliquees en ordre inverse.
        for (std::size_t k = accumulated_axes.size(); k-- > 0;) {
            world_axis = rodrigues_rotate_single(world_axis, accumulated_axes[k],
                                                 deg2rad(accumulated_angles[k]));
        }

        const Vec3 projection = scale3(world_axis, dot3(relative, world_axis));
        if (norm3(sub3(relative, projection)) < epsilon) {
            out_singular.push_back(static_cast<int>(axis_index));
        }

        accumulated_axes.push_back(axis);
        accumulated_angles.push_back(angles_deg[axis_index]);
    }
}


// ---------------------------------------------------------------------------
// Repli par grille : transcription de `_resolve_axis_angles_via_grid_search`.
// np.radians = x * (pi / 180) : l'ordre des operations differe de `deg2rad` ci-dessus d'un ULP,
// et c'est la grille de numpy qui sert d'oracle.
// ---------------------------------------------------------------------------
inline double radians_np(double d) { return d * (kPi / 180.0); }

struct GridSolver {
    const std::vector<AxisSpec>& axes;
    Vec3 reference;
    std::vector<Vec3> axes_unit;          // axes normalises, ordre de la liste
    std::vector<Vec3> directions;         // grille produit, meshgrid 'ij' aplati
    std::vector<double> angles;           // directions.size() * axes.size(), degres
    double step0 = 0.0;
    int iterations = 0;

    GridSolver(const std::vector<AxisSpec>& a, const Vec3& ref, const GridFallback& g)
        : axes(a), reference(ref), iterations(g.refine_iterations) {
        const int samples = g.samples_per_axis > 1 ? g.samples_per_axis : 1;
        for (const auto& spec : axes) axes_unit.push_back(normalized(spec.axis));
        const std::size_t n_axes = axes.size();
        std::vector<std::vector<double>> per_axis(n_axes);
        for (std::size_t k = 0; k < n_axes; ++k) {
            const double lo = axes[k].lo, hi = axes[k].hi;
            const bool full_turn = (hi - lo) >= 360.0;
            const int n = hi <= lo ? 1 : samples;
            std::vector<double>& v = per_axis[k];
            v.resize(static_cast<std::size_t>(n));
            const bool endpoint = !full_turn;
            const int div = endpoint ? n - 1 : n;
            if (n == 1) {
                v[0] = lo;
            } else {
                const double step = (hi - lo) / div;       // numpy.linspace : i * step + start, derniere valeur = stop
                for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(i)] = i * step + lo;
                if (endpoint) v[static_cast<std::size_t>(n - 1)] = hi;
            }
        }
        std::size_t total = 1;
        for (const auto& v : per_axis) total *= v.size();
        directions.resize(total);
        angles.resize(total * n_axes);
        std::vector<std::size_t> idx(n_axes, 0);
        for (std::size_t p = 0; p < total; ++p) {
            Vec3 d = reference;
            for (std::size_t k = n_axes; k-- > 0;) {
                const double deg = per_axis[k][idx[k]];
                d = rodrigues_rotate_single(d, axes_unit[k], radians_np(deg));
            }
            directions[p] = d;
            for (std::size_t k = 0; k < n_axes; ++k) angles[p * n_axes + k] = per_axis[k][idx[k]];
            for (std::size_t k = n_axes; k-- > 0;) {                // 'ij' : le dernier axe varie le plus vite
                if (++idx[k] < per_axis[k].size()) break;
                idx[k] = 0;
            }
        }
        step0 = g.refine_multiplier * (360.0 / samples);
    }

    Vec3 direction_for(const std::vector<double>& a) const {
        Vec3 d = reference;
        for (std::size_t k = axes.size(); k-- > 0;) d = rodrigues_rotate_single(d, axes_unit[k], radians_np(a[k]));
        return d;
    }

    static double angular_error(const Vec3& d, const Vec3& t) {
        double c = d.x * t.x + d.y * t.y + d.z * t.z;
        if (c > 1.0) c = 1.0; else if (c < -1.0) c = -1.0;
        return std::acos(c);
    }

    // Renvoie l'erreur residuelle (radians) ; `out` recoit les angles affines.
    double solve(const Vec3& target, std::vector<double>& out) const {
        std::size_t best = 0;
        double best_sim = -2.0;
        for (std::size_t p = 0; p < directions.size(); ++p) {
            const double sim = directions[p].x * target.x + directions[p].y * target.y + directions[p].z * target.z;
            if (sim > best_sim) { best_sim = sim; best = p; }   // premier maximum, comme np.argmax
        }
        const std::size_t n_axes = axes.size();
        out.assign(angles.begin() + best * n_axes, angles.begin() + (best + 1) * n_axes);
        double best_error = angular_error(direction_for(out), target);
        double step = step0;
        std::vector<double> candidate;
        for (int it = 0; it < std::max(iterations, 0); ++it) {
            for (std::size_t k = 0; k < n_axes; ++k) {
                const double lo = axes[k].lo, hi = axes[k].hi;
                if (hi <= lo) continue;
                for (const double delta : {step, -step}) {
                    candidate = out;
                    double nudged = candidate[k] + delta;
                    if (nudged > hi) nudged = hi; else if (nudged < lo) nudged = lo;
                    candidate[k] = nudged;
                    const double error = angular_error(direction_for(candidate), target);
                    if (error < best_error) { out = candidate; best_error = error; }
                }
            }
            step *= 0.5;
        }
        return best_error;
    }
};

}  // namespace

IKBatchResult resolve_axis_angles_batch(const double* desired_directions, const double* positions,
                                        std::int64_t point_count,
                                        const std::vector<AxisSpec>& axes,
                                        const Vec3& reference_direction,
                                        const Vec3* platform_origin,
                                        double singularity_radius_epsilon,
                                        const double* initial_previous_angles, const GridFallback* grid) {
    const std::size_t axis_count = axes.size();
    std::optional<GridSolver> grid_solver;   // construit a la premiere utilisation (rare)
    std::vector<double> grid_angles;
    IKBatchResult result;
    result.angles.assign(static_cast<std::size_t>(point_count) * axis_count, 0.0);
    result.residuals.assign(static_cast<std::size_t>(point_count), 0.0);

    const Vec3 reference_unit = normalized(reference_direction);

    bool has_previous = initial_previous_angles != nullptr;
    std::vector<double> previous_angles(axis_count, 0.0);
    if (has_previous) {
        for (std::size_t k = 0; k < axis_count; ++k) {
            previous_angles[k] = initial_previous_angles[k];
        }
    }
    std::vector<double> angles(axis_count, 0.0);
    std::vector<double> inner_angles;
    std::vector<int> singular_axes;

    // Transcription de `resolve_axis_angles_analytic::_previous_for`.
    auto previous_for = [&](std::size_t index) -> double {
        if (has_previous && index < previous_angles.size()) {
            return previous_angles[index];
        }
        const double lo = axes[index].lo, hi = axes[index].hi;
        if (hi - lo >= 360.0) {
            return 0.0;
        }
        return (lo <= 0.0 && 0.0 <= hi) ? lo : 0.5 * (lo + hi);
    };

    for (std::int64_t i = 0; i < point_count; ++i) {
        Vec3 target{desired_directions[3 * i + 0], desired_directions[3 * i + 1],
                    desired_directions[3 * i + 2]};
        const double target_norm = norm3(target);
        if (target_norm > 1e-12) {
            target = scale3(target, 1.0 / target_norm);
        }

        if (axis_count == 0) {
            // `_resolve_axis_angles` court-circuite sur une liste d'axes vide
            // AVANT toute normalisation : residuelle nulle, pas l'angle a la
            // reference. Reproduit tel quel.
            result.residuals[static_cast<std::size_t>(i)] = 0.0;
            continue;
        }

        if (axis_count == 1) {
            angles[0] = resolve_single_axis_analytic(axes[0], reference_unit, target,
                                                     previous_for(0));
        } else {
            inner_angles.clear();
            for (std::size_t k = 2; k < axis_count; ++k) {
                inner_angles.push_back(previous_for(k));
            }

            Vec3 reference_effective = reference_unit;
            for (std::size_t k = axis_count; k-- > 2;) {
                const Vec3 axis = normalized(axes[k].axis);
                reference_effective = rodrigues_rotate_single(
                    reference_effective, axis, deg2rad(inner_angles[k - 2]));
            }

            const auto solved = resolve_axis_pair_analytic(axes[0], axes[1], reference_effective,
                                                           target, previous_for(0), previous_for(1));
            if (!solved.has_value()) {
                if (grid == nullptr || grid->samples_per_axis <= 0) {
                    // Repli grille non fourni : on rend la main a Python, qui resout CE point
                    // puis rappelle a partir du suivant (`initial_previous_angles`).
                    result.fallback_index = i;
                    return result;
                }
                if (!grid_solver) grid_solver.emplace(axes, reference_unit, *grid);
                const double residual = grid_solver->solve(target, grid_angles);
                ++result.grid_fallbacks;
                for (std::size_t k = 0; k < axis_count; ++k) angles[k] = grid_angles[k];
                result.residuals[static_cast<std::size_t>(i)] = residual;
                goto after_residual;
            }

            angles[0] = solved->first;
            angles[1] = solved->second;
            for (std::size_t k = 2; k < axis_count; ++k) {
                angles[k] = inner_angles[k - 2];
            }
        }

        {
            const Vec3 achieved = forward_direction(axes, angles, reference_unit);
            result.residuals[static_cast<std::size_t>(i)] = deg2rad(angle_between(achieved, target));
        }
        after_residual:

        // Heritage d'angle a la singularite (step8), AVANT de devenir
        // `previous_angles` du point suivant - l'ordre compte.
        if (platform_origin != nullptr) {
            const Vec3 position{positions[3 * i + 0], positions[3 * i + 1], positions[3 * i + 2]};
            detect_singular_axes(position, *platform_origin, axes, angles,
                                 singularity_radius_epsilon, singular_axes);
            if (!singular_axes.empty() && has_previous) {
                for (const int axis_index : singular_axes) {
                    if (static_cast<std::size_t>(axis_index) < previous_angles.size()) {
                        angles[static_cast<std::size_t>(axis_index)] =
                            previous_angles[static_cast<std::size_t>(axis_index)];
                    }
                }
            }
        }

        for (std::size_t k = 0; k < axis_count; ++k) {
            result.angles[static_cast<std::size_t>(i) * axis_count + k] = angles[k];
        }
        previous_angles = angles;
        has_previous = true;
    }

    return result;
}

}  // namespace apex
