#include "apex/step8_filters.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace apex {

void apply_per_point_normals(const double* local, const double* current,
                             const std::uint8_t* reversed_flag, std::int64_t n, double* desired,
                             std::int64_t* changed) {
    std::int64_t count = 0;
    for (std::int64_t i = 0; i < n; ++i) {
        double d[3] = {local[3 * i], local[3 * i + 1], local[3 * i + 2]};
        if (!(std::isfinite(d[0]) && std::isfinite(d[1]) && std::isfinite(d[2]))) continue;
        if (reversed_flag[i]) {
            d[0] = -d[0]; d[1] = -d[1]; d[2] = -d[2];
        }
        const double dot = d[0] * current[3 * i] + d[1] * current[3 * i + 1] + d[2] * current[3 * i + 2];
        if (dot <= 0.0) continue;
        desired[3 * i] = d[0]; desired[3 * i + 1] = d[1]; desired[3 * i + 2] = d[2];
        ++count;
    }
    if (changed != nullptr) *changed = count;
}

void median_filter_free_axes(double* angles, std::int64_t n, std::int64_t axis_count,
                             const std::int64_t* free_axes, std::int64_t free_count,
                             const std::int64_t* layer_ids, const std::int64_t* segments,
                             const std::uint8_t* travels, std::int64_t window) {
    const std::int64_t half = window / 2;
    std::vector<double> src, win(static_cast<std::size_t>(window));
    std::int64_t start = 0;
    for (std::int64_t i = 1; i <= n; ++i) {
        const bool boundary = i == n || layer_ids[i] != layer_ids[start] || segments[i] != segments[start] ||
                              travels[i] != travels[start] || travels[i];
        if (!boundary) continue;
        const std::int64_t len = i - start;
        if (len > 2 && !travels[start]) {
            for (std::int64_t f = 0; f < free_count; ++f) {
                const std::int64_t ax = free_axes[f];
                src.assign(static_cast<std::size_t>(len), 0.0);
                for (std::int64_t k = 0; k < len; ++k) src[static_cast<std::size_t>(k)] = angles[(start + k) * axis_count + ax];
                for (std::int64_t k = 0; k < len; ++k) {
                    for (std::int64_t w = 0; w < window; ++w) {
                        const std::int64_t j = std::clamp<std::int64_t>(k - half + w, 0, len - 1);
                        win[static_cast<std::size_t>(w)] = src[static_cast<std::size_t>(j)];
                    }
                    std::nth_element(win.begin(), win.begin() + half, win.end());
                    angles[(start + k) * axis_count + ax] = win[static_cast<std::size_t>(half)];
                }
            }
        }
        start = i;
    }
}

}  // namespace apex
