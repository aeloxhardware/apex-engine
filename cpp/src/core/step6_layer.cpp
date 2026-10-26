#include "apex/step6_layer.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "apex/geometry2d.hpp"
#include "apex/interpolate.hpp"

namespace apex {

namespace {

using Flat3 = std::vector<double>;   // K*3 (u, v, offset)
using Flat2 = std::vector<double>;   // K*2

bool allclose2(const double* a, const double* b) {
    for (int i = 0; i < 2; ++i) {
        if (!(std::abs(a[i] - b[i]) <= 1e-8 + 1e-5 * std::abs(b[i]))) return false;
    }
    return true;
}

double signed_area(const Flat2& poly) {
    const std::size_t n = poly.size() / 2;
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        sum += poly[2 * i] * poly[2 * j + 1] - poly[2 * j] * poly[2 * i + 1];
    }
    return 0.5 * sum;
}

Ring2 to_ring(const Flat2& v) {
    Ring2 r(v.size() / 2);
    for (std::size_t i = 0; i < r.size(); ++i) r[i] = {v[2 * i], v[2 * i + 1]};
    return r;
}

// --- troncons monotones ----------------------------------------------------

struct Run { std::int64_t start, stop; };

std::vector<Run> monotonic_runs(const std::vector<double>& off, double min_amp) {
    const auto n = static_cast<std::int64_t>(off.size());
    if (n <= 2) return {{0, n}};
    std::vector<Run> runs;
    std::int64_t start = 0, ext_index = 0;
    int direction = 0;
    double low = off[0], high = off[0], ext_value = off[0];
    for (std::int64_t i = 1; i < n; ++i) {
        const double value = off[static_cast<std::size_t>(i)];
        if (direction == 0) {
            low = std::min(low, value);
            high = std::max(high, value);
            if (value - low > min_amp) { direction = 1; ext_index = i; ext_value = value; }
            else if (high - value > min_amp) { direction = -1; ext_index = i; ext_value = value; }
        } else if (direction == 1) {
            if (value >= ext_value) { ext_index = i; ext_value = value; }
            else if (ext_value - value > min_amp) {
                runs.push_back({start, ext_index + 1});
                start = ext_index + 1;
                direction = -1; ext_index = i; ext_value = value;
            }
        } else {
            if (value <= ext_value) { ext_index = i; ext_value = value; }
            else if (value - ext_value > min_amp) {
                runs.push_back({start, ext_index + 1});
                start = ext_index + 1;
                direction = 1; ext_index = i; ext_value = value;
            }
        }
    }
    runs.push_back({start, n});
    return runs;
}

std::vector<Run> merge_short_runs(const std::vector<Run>& runs, const std::vector<double>& off) {
    std::vector<Run> merged;
    bool carried = false;
    std::int64_t carried_start = 0;
    for (Run r : runs) {
        if (carried) { r.start = carried_start; carried = false; }
        if (r.stop - r.start >= 2) { merged.push_back(r); continue; }
        if (!merged.empty() && off[static_cast<std::size_t>(r.start)] >= off[static_cast<std::size_t>(merged.back().stop - 1)]) {
            merged.back().stop = r.stop;
        } else {
            carried = true;
            carried_start = r.start;
        }
    }
    if (carried) {
        if (!merged.empty()) merged.back().stop = runs.back().stop;
        else merged.push_back({carried_start, runs.back().stop});
    }
    return merged;
}

bool is_closed_ring(const Flat3& poly) {
    const std::size_t n = poly.size() / 3;
    if (n < 3) return false;
    const double* first = poly.data();
    const double* last = poly.data() + 3 * (n - 1);
    if (allclose2(first, last)) return true;
    const double gap = std::hypot(last[0] - first[0], last[1] - first[1]);
    std::vector<double> steps(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        steps[i] = std::hypot(poly[3 * (i + 1)] - poly[3 * i], poly[3 * (i + 1) + 1] - poly[3 * i + 1]);
    }
    std::sort(steps.begin(), steps.end());
    const std::size_t m = steps.size();
    const double typical = (m % 2 == 1) ? steps[m / 2] : 0.5 * (steps[m / 2 - 1] + steps[m / 2]);
    return typical > 0.0 && gap <= 2.0 * typical;
}

void split_ring(const Flat3& ring, bool share, double min_amp, std::vector<Flat3>& out) {
    const std::size_t n = ring.size() / 3;
    if (n <= 2) { out.push_back(ring); return; }
    std::vector<double> off(n);
    for (std::size_t i = 0; i < n; ++i) off[i] = ring[3 * i + 2];
    const auto runs = monotonic_runs(off, min_amp);
    if (runs.size() <= 1) { out.push_back(ring); return; }
    const auto merged = merge_short_runs(runs, off);
    for (std::size_t k = 0; k < merged.size(); ++k) {
        std::int64_t start = merged[k].start;
        const std::int64_t stop = merged[k].stop;
        if (share && k > 0 && start > 0) --start;
        Flat3 seg(ring.begin() + 3 * start, ring.begin() + 3 * stop);
        const std::size_t len = seg.size() / 3;
        if (len >= 2 && seg[3 * (len - 1) + 2] < seg[2]) {
            Flat3 rev(seg.size());
            for (std::size_t i = 0; i < len; ++i) {
                for (int c = 0; c < 3; ++c) rev[3 * i + c] = seg[3 * (len - 1 - i) + c];
            }
            seg = std::move(rev);
        }
        out.push_back(std::move(seg));
    }
}

}  // namespace

std::vector<std::vector<double>> fill_layer(const LayerContours& layer, const Step6Config& cfg,
                                            const BridgeInfo& bridge) {
    const double* origin = layer.frame;
    const double* u = layer.frame + 3;
    const double* v = layer.frame + 6;
    const double normal[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};

    const std::size_t loops = layer.loops.size();
    std::vector<Flat2> loop2d(loops);       // projetees, point de fermeture compris
    std::vector<std::vector<double>> loop_off(loops);
    for (std::size_t k = 0; k < loops; ++k) {
        const auto n = static_cast<std::size_t>(layer.loop_sizes[k]);
        loop2d[k].resize(2 * n);
        loop_off[k].resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double* p = layer.loops[k] + 3 * i;
            const double r[3] = {p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]};
            loop2d[k][2 * i] = r[0] * u[0] + r[1] * u[1] + r[2] * u[2];
            loop2d[k][2 * i + 1] = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
            loop_off[k][i] = r[0] * normal[0] + r[1] * normal[1] + r[2] * normal[2];
        }
    }
    if (loops == 0) return {};

    // --- classification (_classify_contour_indices) ---
    std::vector<Flat2> open_loops(loops);
    for (std::size_t k = 0; k < loops; ++k) {
        const std::size_t n = loop2d[k].size() / 2;
        const bool closed = n > 0 && allclose2(loop2d[k].data(), loop2d[k].data() + 2 * (n - 1));
        open_loops[k].assign(loop2d[k].begin(), loop2d[k].begin() + 2 * (closed ? n - 1 : n));
    }
    std::vector<double> area(loops);
    for (std::size_t k = 0; k < loops; ++k) area[k] = std::abs(signed_area(open_loops[k]));
    std::vector<std::size_t> order(loops);
    for (std::size_t k = 0; k < loops; ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return area[a] > area[b]; });
    std::vector<std::vector<std::size_t>> containers(loops);
    for (std::size_t rank = 0; rank < loops; ++rank) {
        const std::size_t i = order[rank];
        const double* probe = open_loops[i].data();
        for (std::size_t r2 = 0; r2 < rank; ++r2) {
            const std::size_t j = order[r2];
            if (point_in_polygon_2d(probe[0], probe[1], open_loops[j].data(),
                                    static_cast<std::int64_t>(open_loops[j].size() / 2))) {
                containers[i].push_back(j);
            }
        }
    }
    std::vector<char> is_hole(loops);
    for (std::size_t k = 0; k < loops; ++k) is_hole[k] = containers[k].size() % 2 == 1;
    std::vector<std::vector<std::size_t>> groups(loops);
    for (std::size_t k = 0; k < loops; ++k) {
        if (!is_hole[k]) continue;
        long parent = -1;
        for (const std::size_t p : containers[k]) {
            if (is_hole[p]) continue;
            if (parent < 0 || area[p] < area[static_cast<std::size_t>(parent)]) parent = static_cast<long>(p);
        }
        if (parent >= 0) groups[static_cast<std::size_t>(parent)].push_back(k);
    }

    std::vector<Flat3> toolpath;
    for (std::size_t outer_idx = 0; outer_idx < loops; ++outer_idx) {
        if (is_hole[outer_idx]) continue;
        const auto& hole_idxs = groups[outer_idx];
        std::vector<char> member(loops, 0);
        member[outer_idx] = 1;
        for (const auto h : hole_idxs) member[h] = 1;

        Flat2 source2d;
        std::vector<double> source_off;
        for (std::size_t k = 0; k < loops; ++k) {
            if (!member[k]) continue;
            const std::size_t open_n = open_loops[k].size() / 2;
            source2d.insert(source2d.end(), open_loops[k].begin(), open_loops[k].end());
            source_off.insert(source_off.end(), loop_off[k].begin(), loop_off[k].begin() + static_cast<long>(open_n));
        }
        const std::size_t outer_n = open_loops[outer_idx].size() / 2;

        // --- detection de pont sur le contour exterieur ---
        bool has_bridge = false;
        double span2d[2] = {0, 0};
        if (bridge.tree != nullptr) {
            const auto K = static_cast<std::int64_t>(layer.loop_sizes[outer_idx]);
            std::vector<std::int64_t> faces(static_cast<std::size_t>(K));
            std::vector<double> dist(static_cast<std::size_t>(K));
            bridge.tree->query_nn(layer.loops[outer_idx], K, faces.data(), dist.data());
            std::int64_t bridged = 0;
            for (const auto f : faces) bridged += bridge.is_bridge[f] ? 1 : 0;
            const double ratio = K ? static_cast<double>(bridged) / static_cast<double>(K) : 0.0;
            if (!(ratio < cfg.bridge_contour_ratio_threshold)) {
                double sum[3] = {0, 0, 0};
                std::int64_t valid = 0;
                for (const auto f : faces) {
                    if (!bridge.is_bridge[f]) continue;
                    const double* s = bridge.span + 3 * f;
                    if (std::isnan(s[0]) || std::isnan(s[1]) || std::isnan(s[2])) continue;
                    sum[0] += s[0]; sum[1] += s[1]; sum[2] += s[2];
                    ++valid;
                }
                if (valid > 0) {
                    double m[3] = {sum[0] / static_cast<double>(valid), sum[1] / static_cast<double>(valid), sum[2] / static_cast<double>(valid)};
                    const double mn = std::sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
                    if (!(mn < 1e-12)) {
                        for (double& c : m) c /= mn;
                        const double a = m[0] * u[0] + m[1] * u[1] + m[2] * u[2];
                        const double b = m[0] * v[0] + m[1] * v[1] + m[2] * v[2];
                        const double sn = std::sqrt(a * a + b * b);
                        if (!(sn < 1e-12)) {
                            span2d[0] = a / sn; span2d[1] = b / sn;
                            has_bridge = true;
                        }
                    }
                }
            }
        }

        // --- parois + remplissage ---
        const Ring2 outer_ring = to_ring(open_loops[outer_idx]);
        std::vector<Ring2> hole_rings;
        for (const auto h : hole_idxs) hole_rings.push_back(to_ring(open_loops[h]));
        FillRegionParams fp;
        fp.nozzle_width = cfg.nozzle_width;
        fp.degenerate_area_ratio = cfg.degenerate_area_ratio;
        fp.infill_density = cfg.infill_density;
        fp.pattern = cfg.pattern;
        fp.wall_count = cfg.wall_count;
        fp.sweep_angle_rad = (cfg.pattern == InfillKind::CrossedZigzag) ? layer.sweep_angle : 0.0;
        // categorie : 0 infill, 1 paroi n>=1, 2 contour exterieur, 3 trou
        struct Generated { Ring2 ring; int kind; std::size_t hole_loop; };
        std::vector<Generated> generated;
        if (has_bridge) {
            fp.wall_count = cfg.wall_count + cfg.bridge_anchor_extra_walls;
            fp.sweep_angle_rad = std::atan2(span2d[1], span2d[0]) + std::numbers::pi / 2.0;
            FillRegionResult r = fill_region(outer_ring, hole_rings, fp);
            if (!r.base_valid) {
                generated.push_back({outer_ring, 2, 0});
            } else {
                for (auto& s : r.infill) generated.push_back({std::move(s), 0, 0});
                for (auto& w : r.walls) generated.push_back({std::move(w), 1, 0});
                generated.push_back({outer_ring, 2, 0});
                for (std::size_t h = 0; h < hole_rings.size(); ++h) generated.push_back({hole_rings[h], 3, hole_idxs[h]});
            }
        } else {
            FillRegionResult r = fill_region(outer_ring, hole_rings, fp);
            if (!r.base_valid) continue;
            for (auto& s : r.infill) generated.push_back({std::move(s), 0, 0});
            for (auto& w : r.walls) generated.push_back({std::move(w), 1, 0});
            generated.push_back({outer_ring, 2, 0});
            for (std::size_t h = 0; h < hole_rings.size(); ++h) generated.push_back({hole_rings[h], 3, hole_idxs[h]});
        }

        // --- offset hors-plan, format (u, v, offset), fermeture des contours d'origine ---
        for (auto& g : generated) {
            const std::size_t n = g.ring.size();
            std::vector<double> off(n, 0.0);
            if (g.kind == 2) {
                for (std::size_t i = 0; i < n && i < outer_n; ++i) off[i] = loop_off[outer_idx][i];
            } else if (cfg.has_power && !source_off.empty() && n > 0) {
                std::vector<double> q(2 * n);
                for (std::size_t i = 0; i < n; ++i) { q[2 * i] = g.ring[i][0]; q[2 * i + 1] = g.ring[i][1]; }
                interpolate_out_of_plane(q.data(), static_cast<std::int64_t>(n), source2d.data(), source_off.data(),
                                         static_cast<std::int64_t>(source_off.size()), cfg.power, off.data());
            }
            Flat3 arr(3 * n);
            for (std::size_t i = 0; i < n; ++i) {
                arr[3 * i] = g.ring[i][0]; arr[3 * i + 1] = g.ring[i][1]; arr[3 * i + 2] = off[i];
            }
            if (cfg.continuous_rings && (g.kind == 2 || g.kind == 3) && n >= 3 &&
                !(arr[0] == arr[3 * (n - 1)] && arr[1] == arr[3 * (n - 1) + 1])) {
                arr.insert(arr.end(), arr.begin(), arr.begin() + 3);
            }
            toolpath.push_back(std::move(arr));
        }
    }

    if (cfg.split_rings_into_runs) {
        std::vector<Flat3> result;
        for (auto& poly : toolpath) {
            if (is_closed_ring(poly) || cfg.split_open_polylines) {
                split_ring(poly, cfg.continuous_rings, cfg.run_split_min_amplitude, result);
            } else {
                result.push_back(std::move(poly));
            }
        }
        toolpath = std::move(result);
    }
    return toolpath;
}

}  // namespace apex
