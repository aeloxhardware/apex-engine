// Voir apex/scheduling.hpp : transcription de pipeline/step7_scheduling.py.

#include "apex/scheduling.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>

#include "apex/bvh.hpp"
#include "apex/kdtree.hpp"
#include "apex/mesh_distance.hpp"

namespace apex {

namespace {

struct Frame {
    double o[3], u[3], v[3], n[3];
    explicit Frame(const double* f) {
        for (int k = 0; k < 3; ++k) { o[k] = f[k]; u[k] = f[3 + k]; v[k] = f[6 + k]; }
        // step6._frame_normal
        n[0] = u[1] * v[2] - u[2] * v[1];
        n[1] = u[2] * v[0] - u[0] * v[2];
        n[2] = u[0] * v[1] - u[1] * v[0];
    }
    // step6._lift_to_3d : origin + a*u + b*v (+ c*n), meme ordre d'operations
    void lift(const double* row, std::int64_t columns, double* out) const {
        for (int k = 0; k < 3; ++k) {
            double x = o[k] + row[0] * u[k];
            x = x + row[1] * v[k];
            if (columns >= 3) x = x + row[2] * n[k];
            out[k] = x;
        }
    }
};

struct Group {
    std::int64_t layer_index, layer_id;
    Color color;
    std::vector<std::int64_t> segments;
    std::vector<std::int64_t> regions;
};

struct Run {
    std::size_t start, stop;
    Color color;
};

std::vector<Run> color_runs(const std::vector<Color>& colors) {
    std::vector<Run> runs;
    std::size_t start = 0;
    for (std::size_t i = 1; i <= colors.size(); ++i)
        if (i == colors.size() || colors[i] != colors[start]) {
            runs.push_back({start, i, colors[start]});
            start = i;
        }
    return runs;
}

}  // namespace

ScheduleResult schedule_layers(const std::vector<ScheduleLayer>& layers, const ScheduleInput& in) {
    ScheduleResult out;
    const std::size_t L = layers.size();
    const bool painted = !in.face_color.empty() &&
                         std::any_of(in.face_color.begin(), in.face_color.end(), [](const Color& c) { return c.has_value(); });

    std::unique_ptr<MeshDistanceIndex> surface;
    if (painted) {
        MeshView view{in.vertices, static_cast<std::size_t>(in.vertex_count), in.faces,
                      static_cast<std::size_t>(in.face_count)};
        surface = std::make_unique<MeshDistanceIndex>(view);
    }

    // --- 1. couleur par arete, lissage, decoupage (_apply_colors_to_layers) ---
    out.toolpath.resize(L);
    out.segment_colors.resize(L);
    out.dominant_color.resize(L);
    for (std::size_t li = 0; li < L; ++li) {
        const auto& layer = layers[li];
        const std::int64_t C = layer.columns;
        if (!painted) {
            out.toolpath[li] = layer.toolpath;
            out.segment_colors[li].assign(layer.toolpath.size(), std::nullopt);
            continue;
        }
        const Frame frame(layer.frame);
        for (const auto& poly : layer.toolpath) {
            const std::size_t K = poly.size() / C;
            if (K < 2) continue;      // ni arete ni sous-polyligne (len < 2 rejetee)
            std::vector<double> mids((K - 1) * 3), dist(K - 1), lengths(K - 1);
            std::vector<double> row(C);
            for (std::size_t e = 0; e + 1 < K; ++e) {
                double sq = 0.0;
                for (std::int64_t c = 0; c < C; ++c) {
                    row[c] = 0.5 * (poly[e * C + c] + poly[(e + 1) * C + c]);
                    const double d = poly[(e + 1) * C + c] - poly[e * C + c];
                    sq += d * d;
                }
                lengths[e] = std::sqrt(sq);
                frame.lift(row.data(), C, mids.data() + e * 3);
            }
            std::vector<std::int64_t> face(K - 1);
            surface->distances(mids.data(), K - 1, dist.data(), nullptr, DistanceBackend::CpuSerial, 0, face.data());
            std::vector<Color> colors(K - 1);
            for (std::size_t e = 0; e + 1 < K; ++e) {
                if (dist[e] > in.params.surface_color_distance || face[e] < 0) colors[e] = std::nullopt;
                else colors[e] = in.face_color[face[e]];
            }
            // _smooth_edge_colors
            if (in.merge_tolerance > 0) {
                while (true) {
                    const auto runs = color_runs(colors);
                    if (runs.size() <= 1) break;
                    std::vector<double> rl(runs.size(), 0.0);
                    for (std::size_t r = 0; r < runs.size(); ++r)
                        for (std::size_t e = runs[r].start; e < runs[r].stop; ++e) rl[r] += lengths[e];
                    std::size_t shortest = 0;
                    for (std::size_t r = 1; r < rl.size(); ++r)
                        if (rl[r] < rl[shortest]) shortest = r;
                    if (rl[shortest] >= in.merge_tolerance) break;
                    const double prev = shortest > 0 ? rl[shortest - 1] : -1.0;
                    const double next = shortest + 1 < runs.size() ? rl[shortest + 1] : -1.0;
                    const Color absorbing = prev >= next ? runs[shortest - 1].color : runs[shortest + 1].color;
                    for (std::size_t e = runs[shortest].start; e < runs[shortest].stop; ++e) colors[e] = absorbing;
                }
            }
            for (const auto& r : color_runs(colors)) {
                std::vector<double> sub(poly.begin() + r.start * C, poly.begin() + (r.stop + 1) * C);
                out.toolpath[li].push_back(std::move(sub));
                out.segment_colors[li].push_back(r.color);
            }
        }
        // _dominant_color : premiere des plus frequentes, ordre de premiere apparition
        std::vector<std::pair<std::int64_t, std::int64_t>> counts;
        for (const auto& c : out.segment_colors[li]) {
            if (!c) continue;
            auto it = std::find_if(counts.begin(), counts.end(), [&](const auto& p) { return p.first == *c; });
            if (it == counts.end()) counts.push_back({*c, 1});
            else it->second++;
        }
        if (!counts.empty()) {
            auto best = counts.begin();
            for (auto it = counts.begin(); it != counts.end(); ++it)
                if (it->second > best->second) best = it;
            out.dominant_color[li] = best->first;
        }
    }

    // --- 2. groupes (couche, couleur), ordre de premiere apparition ---
    std::vector<std::vector<Group>> groups(L);
    for (std::size_t li = 0; li < L; ++li) {
        for (std::size_t s = 0; s < out.segment_colors[li].size(); ++s) {
            const Color c = out.segment_colors[li][s];
            auto it = std::find_if(groups[li].begin(), groups[li].end(), [&](const Group& g) { return g.color == c; });
            if (it == groups[li].end()) {
                groups[li].push_back({static_cast<std::int64_t>(li), layers[li].id, c, {}, layers[li].region_ids});
                it = groups[li].end() - 1;
            }
            it->segments.push_back(static_cast<std::int64_t>(s));
        }
    }

    // --- 3. deadlines d'accessibilite (compute_hole_deadline) ---
    bool has_iso = L > 0;
    std::vector<double> sorted_iso;
    for (const auto& l : layers) {
        if (!l.isovalue) has_iso = false;
        else sorted_iso.push_back(*l.isovalue);
    }
    std::sort(sorted_iso.begin(), sorted_iso.end());
    std::unique_ptr<TriangleBVH> bvh;
    std::unique_ptr<KDTree> field_tree;
    auto ensure_tools = [&]() {
        if (bvh) return;
        std::vector<double> tri(in.face_count * 9);
        for (std::int64_t f = 0; f < in.face_count; ++f)
            for (int k = 0; k < 3; ++k)
                for (int c = 0; c < 3; ++c) tri[f * 9 + k * 3 + c] = in.vertices[in.faces[f * 3 + k] * 3 + c];
        bvh = std::make_unique<TriangleBVH>(tri.data(), in.face_normals, in.face_count);
        field_tree = std::make_unique<KDTree>(in.tet_vertices, in.tet_vertex_count);
    };
    const std::int64_t layer_count = static_cast<std::int64_t>(L);
    std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t> cache;   // (couche, couleur ou -inf)
    auto deadline_of = [&](const Group& g) -> std::int64_t {
        if (!has_iso) return layer_count - 1;
        const auto key = std::make_pair(g.layer_index, g.color ? *g.color : std::numeric_limits<std::int64_t>::min());
        auto found = cache.find(key);
        if (found != cache.end()) return found->second;
        ensure_tools();
        // _group_position_3d : centroide des segments, releve en 3D
        const auto& layer = layers[g.layer_index];
        const auto& tp = out.toolpath[g.layer_index];
        const std::int64_t C = layer.columns;
        std::vector<double> mean(C, 0.0);
        std::size_t count = 0;
        for (std::int64_t s : g.segments) {
            if (s >= static_cast<std::int64_t>(tp.size())) continue;
            const auto& poly = tp[s];
            for (std::size_t p = 0; p < poly.size() / C; ++p, ++count)
                for (std::int64_t c = 0; c < C; ++c) mean[c] += poly[p * C + c];
        }
        double origin[3];
        const Frame frame(layer.frame);
        if (count == 0) {
            for (int k = 0; k < 3; ++k) origin[k] = frame.o[k];
        } else {
            for (auto& m : mean) m /= double(count);
            frame.lift(mean.data(), C, origin);
        }
        const std::int64_t D = static_cast<std::int64_t>(in.escape_directions.size() / 3);
        const std::int64_t samples = in.params.deadline_ray_samples;
        const double eps = 1e-6;
        std::int64_t best_blocking = -1;
        for (std::int64_t d = 0; d < D; ++d) {
            const double* dir = in.escape_directions.data() + d * 3;
            const auto hits = bvh->intersects_location(origin, dir, 1);
            std::int64_t blocking = layer_count;
            double exit = std::numeric_limits<double>::infinity();
            for (const auto& h : hits) {
                const double dd = std::sqrt((h.location.x - origin[0]) * (h.location.x - origin[0]) +
                                            (h.location.y - origin[1]) * (h.location.y - origin[1]) +
                                            (h.location.z - origin[2]) * (h.location.z - origin[2]));
                if (dd > eps) exit = std::min(exit, dd);
            }
            if (std::isfinite(exit) && samples > 0) {
                std::vector<double> pts(samples * 3);
                const double step = samples > 1 ? (exit - eps) / double(samples - 1) : 0.0;
                for (std::int64_t k = 0; k < samples; ++k) {
                    const double t = (k == samples - 1 && samples > 1) ? exit : eps + double(k) * step;
                    for (int c = 0; c < 3; ++c) pts[k * 3 + c] = origin[c] + dir[c] * t;
                }
                std::vector<std::int64_t> idx(samples);
                std::vector<double> nd(samples);
                field_tree->query_nn(pts.data(), samples, idx.data(), nd.data());
                std::int64_t future = std::numeric_limits<std::int64_t>::max();
                for (std::int64_t k = 0; k < samples; ++k) {
                    const double phi = in.field_values[idx[k]];
                    const std::int64_t li = std::lower_bound(sorted_iso.begin(), sorted_iso.end(), phi) - sorted_iso.begin();
                    if (li > g.layer_index) future = std::min(future, li);
                }
                if (future != std::numeric_limits<std::int64_t>::max()) blocking = future;
            }
            best_blocking = std::max(best_blocking, blocking);
        }
        const std::int64_t deadline = std::clamp(best_blocking - 1, g.layer_index, layer_count - 1);
        cache[key] = deadline;
        return deadline;
    };

    // --- 4. ordonnancement (_schedule_passes) ---
    auto make_pass = [&](const std::vector<const Group*>& gs, Color color, bool fill) -> std::optional<SchedulePass> {
        SchedulePass p;
        for (const Group* g : gs)
            for (std::int64_t s : g->segments) p.segments.push_back({g->layer_id, s});
        if (p.segments.empty()) return std::nullopt;
        p.color = color;
        p.is_hole_fill = fill;
        p.layer_min = p.layer_max = p.segments[0].first;
        for (const auto& s : p.segments) {
            p.layer_min = std::min(p.layer_min, s.first);
            p.layer_max = std::max(p.layer_max, s.first);
        }
        if (fill) {
            for (const Group* g : gs) p.deferred_regions.insert(p.deferred_regions.end(), g->regions.begin(), g->regions.end());
            std::sort(p.deferred_regions.begin(), p.deferred_regions.end());
            p.deferred_regions.erase(std::unique(p.deferred_regions.begin(), p.deferred_regions.end()), p.deferred_regions.end());
        }
        return p;
    };
    struct Pending { const Group* group; std::int64_t deadline; };
    std::vector<Pending> pending;
    Color active;
    const double purge_cost = in.params.time_per_purge + in.params.volume_per_purge;
    auto flush = [&](std::int64_t up_to, bool force_all) {
        std::stable_sort(pending.begin(), pending.end(),
                         [](const Pending& a, const Pending& b) { return a.deadline < b.deadline; });
        std::vector<Pending> remaining;
        for (const auto& h : pending) {
            if (!(force_all || h.deadline <= up_to)) { remaining.push_back(h); continue; }
            if (auto p = make_pass({h.group}, h.group->color, true)) {
                out.passes.push_back(std::move(*p));
                active = h.group->color;
            }
        }
        pending = std::move(remaining);
    };
    for (std::size_t li = 0; li < L; ++li) {
        std::vector<const Group*> free_groups, colored;
        for (const auto& g : groups[li]) (g.color ? colored : free_groups).push_back(&g);
        if (colored.empty()) {
            if (auto p = make_pass(free_groups, active, false)) out.passes.push_back(std::move(*p));
            flush(static_cast<std::int64_t>(li), false);
            continue;
        }
        const Group* chosen = nullptr;
        for (const Group* g : colored)
            if (active && g->color == active) chosen = g;   // dict : la derniere occurrence d'une couleur gagne
        if (!chosen) {
            chosen = colored[0];
            for (const Group* g : colored)
                if (g->segments.size() > chosen->segments.size()) chosen = g;
        }
        std::vector<const Group*> printed{chosen};
        printed.insert(printed.end(), free_groups.begin(), free_groups.end());
        if (auto p = make_pass(printed, chosen->color, false)) out.passes.push_back(std::move(*p));
        active = chosen->color;
        for (const Group* g : colored) {
            if (g == chosen) continue;
            const std::int64_t deadline = deadline_of(*g);
            const std::int64_t span = deadline - static_cast<std::int64_t>(li);
            const bool defer = deadline >= static_cast<std::int64_t>(li) + 1 &&
                               static_cast<std::int64_t>(pending.size()) < in.params.max_open_holes &&
                               in.params.reheat_penalty + in.params.adhesion_risk_weight * double(std::max<std::int64_t>(span, 0)) <= purge_cost;
            if (defer) {
                pending.push_back({g, deadline});
            } else {
                if (auto p = make_pass({g}, g->color, false)) out.passes.push_back(std::move(*p));
                active = g->color;
            }
        }
        flush(static_cast<std::int64_t>(li), false);
    }
    flush(layer_count, true);
    return out;
}

}  // namespace apex
