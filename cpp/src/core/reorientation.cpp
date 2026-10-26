#include "apex/reorientation.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace apex {

namespace {

constexpr double kDeg2Rad = std::numbers::pi / 180.0;

// Rodrigues d'un point autour de `axis` unitaire : v*cos + (axis x v)*sin +
// axis*(v.axis)*(1-cos), meme expression que core.geometry_utils._rodrigues_rotate.
inline void rodrigues(double v[3], const double* axis, double c, double s) {
    const double dot = v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2];
    const double cx = axis[1] * v[2] - axis[2] * v[1];
    const double cy = axis[2] * v[0] - axis[0] * v[2];
    const double cz = axis[0] * v[1] - axis[1] * v[0];
    const double k = 1.0 - c;
    const double r0 = v[0] * c + cx * s + axis[0] * dot * k;
    const double r1 = v[1] * c + cy * s + axis[1] * dot * k;
    const double r2 = v[2] * c + cz * s + axis[2] * dot * k;
    v[0] = r0; v[1] = r1; v[2] = r2;
}

// Rotations en ordre INVERSE de la liste d'axes (axe interne d'abord).
void head_direction(const ReorientationInput& in, const double* angles, double out[3]) {
    out[0] = 0; out[1] = 0; out[2] = 1;
    for (std::int64_t k = in.axis_count - 1; k >= 0; --k) {
        const double a = angles[k] * kDeg2Rad;
        rodrigues(out, in.axes + 3 * k, std::cos(a), std::sin(a));
    }
}

}  // namespace

std::vector<ReorientationEntry> reorient_coupled_tiers(const ReorientationInput& in,
                                                       const PastMaterialIndex& material,
                                                       const std::int64_t* contacts,
                                                       std::int64_t contact_count,
                                                       double* new_angles) {
    const std::int64_t A = in.axis_count;
    std::copy(in.angles, in.angles + in.n * A, new_angles);
    std::vector<ReorientationEntry> log;
    log.reserve(static_cast<std::size_t>(contact_count));
    std::vector<double> verts(static_cast<std::size_t>(in.sample_count) * 3);
    std::vector<double> devs(static_cast<std::size_t>(in.grid_count));
    std::vector<std::int64_t> order(static_cast<std::size_t>(in.grid_count));
    std::vector<double> a(static_cast<std::size_t>(A));

    for (std::int64_t c = 0; c < contact_count; ++c) {
        const std::int64_t idx = contacts[c];
        const double* a0 = in.angles + idx * A;
        double d0[3];
        head_direction(in, a0, d0);

        for (std::int64_t g = 0; g < in.grid_count; ++g) order[static_cast<std::size_t>(g)] = g;
        const double cur = a0[in.free_axis];
        std::stable_sort(order.begin(), order.end(), [&](std::int64_t x, std::int64_t y) {
            return std::abs(in.grid[x] - cur) < std::abs(in.grid[y] - cur);
        });
        for (std::int64_t g = 0; g < in.grid_count; ++g) {
            for (std::int64_t k = 0; k < A; ++k) a[static_cast<std::size_t>(k)] = a0[k];
            a[static_cast<std::size_t>(in.free_axis)] = in.grid[g];
            double d[3];
            head_direction(in, a.data(), d);
            const double dot = std::clamp(d[0] * d0[0] + d[1] * d0[1] + d[2] * d0[2], -1.0, 1.0);
            devs[static_cast<std::size_t>(g)] = std::acos(dot) / kDeg2Rad;
        }

        // Controle d'un jeu d'angles : plateau puis matiere passee.
        auto clear_at = [&](const std::vector<double>& ang) {
            // Positionnement de la tete : rotations d'axe interne a externe, puis translation.
            const double* p = in.positions + 3 * idx;
            std::vector<double> cs(static_cast<std::size_t>(A)), sn(static_cast<std::size_t>(A));
            for (std::int64_t k = 0; k < A; ++k) {
                cs[static_cast<std::size_t>(k)] = std::cos(ang[static_cast<std::size_t>(k)] * kDeg2Rad);
                sn[static_cast<std::size_t>(k)] = std::sin(ang[static_cast<std::size_t>(k)] * kDeg2Rad);
            }
            for (std::int64_t s = 0; s < in.sample_count; ++s) {
                double v[3] = {in.samples[3 * s], in.samples[3 * s + 1], in.samples[3 * s + 2]};
                for (std::int64_t k = A - 1; k >= 0; --k) {
                    rodrigues(v, in.axes + 3 * k, cs[static_cast<std::size_t>(k)], sn[static_cast<std::size_t>(k)]);
                }
                verts[3 * static_cast<std::size_t>(s)] = v[0] + p[0];
                verts[3 * static_cast<std::size_t>(s) + 1] = v[1] + p[1];
                verts[3 * static_cast<std::size_t>(s) + 2] = v[2] + p[2];
            }
            if (in.has_plate) {
                double minimum = INFINITY;
                for (std::int64_t s = 0; s < in.sample_count; ++s) {
                    const double lx = verts[3 * static_cast<std::size_t>(s)] - in.plate_origin[0];
                    const double ly = verts[3 * static_cast<std::size_t>(s) + 1] - in.plate_origin[1];
                    const double lz = verts[3 * static_cast<std::size_t>(s) + 2] - in.plate_origin[2];
                    const double dr = std::max(std::hypot(lx, ly) - in.plate_radius, 0.0);
                    const double dz = std::max(std::abs(lz - in.plate_center_offset) - in.plate_half_height, 0.0);
                    minimum = std::min(minimum, std::hypot(dr, dz));
                }
                if (minimum < in.safety_margin) return false;
            }
            const double best = material.min_past_distance(verts.data(), in.sample_count, idx, in.reach);
            return !(best < in.reach);
        };

        bool found = false;
        double found_angle = 0.0;
        std::int64_t found_tier = 0;
        for (std::int64_t tier = 0; tier < in.tier_count && !found; ++tier) {
            const double tol = in.tiers[tier];
            for (std::int64_t oi = 0; oi < in.grid_count; ++oi) {
                const std::int64_t g = order[static_cast<std::size_t>(oi)];
                const double dv = devs[static_cast<std::size_t>(g)];
                if (dv > tol || (tier > 0 && dv <= in.tiers[tier - 1])) continue;
                for (std::int64_t k = 0; k < A; ++k) a[static_cast<std::size_t>(k)] = a0[k];
                a[static_cast<std::size_t>(in.free_axis)] = in.grid[g];
                if (clear_at(a)) {
                    found = true;
                    found_angle = in.grid[g];
                    found_tier = tier;
                    break;
                }
            }
        }
        log.push_back({idx, found, found_angle, found_tier});
        if (!found) continue;

        const double delta = found_angle - a0[in.free_axis];
        std::int64_t start = idx;
        while (start > 0 && in.extrusion[start] > 0) --start;
        std::int64_t stop = idx;
        while (stop + 1 < in.n && in.extrusion[stop + 1] > 0) ++stop;
        // abscisse curviligne cumulee sur [start, stop]
        std::vector<double> seg(static_cast<std::size_t>(stop - start + 1), 0.0);
        for (std::int64_t i = start + 1; i <= stop; ++i) {
            const double* p = in.positions + 3 * i;
            const double* q = in.positions + 3 * (i - 1);
            const double dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2];
            seg[static_cast<std::size_t>(i - start)] =
                seg[static_cast<std::size_t>(i - 1 - start)] + std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        const double s0 = seg[static_cast<std::size_t>(idx - start)];
        const double denom = std::max(in.taper_mm - in.full_mm, 1e-9);
        for (std::int64_t i = start; i <= stop; ++i) {
            const double s = std::abs(seg[static_cast<std::size_t>(i - start)] - s0);
            const double w = std::clamp((in.taper_mm - s) / denom, 0.0, 1.0);
            if (w > 0.0) {
                double& v = new_angles[i * A + in.free_axis];
                v = std::clamp(v + w * delta, in.lo, in.hi);
            }
        }
        new_angles[idx * A + in.free_axis] = found_angle;
    }
    return log;
}

}  // namespace apex
