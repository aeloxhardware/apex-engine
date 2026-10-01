// Voir apex/segmentation.hpp : transcription du chemin par defaut de
// pipeline/step3_segmentation.py, ordres compris.

#include "apex/segmentation.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <numeric>
#include <unordered_map>

namespace apex {

namespace {

constexpr double kRadToDeg = 180.0 / 3.14159265358979323846;

double angle_deg(const double* a, const double* b) {
    double d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    d = std::clamp(d, -1.0, 1.0);
    return std::acos(d) * kRadToDeg;
}

}  // namespace

SegmentationResult segment_default(const SegmentationInput& in) {
    const std::int64_t F = in.face_count, D = in.direction_count;
    SegmentationResult out;

    // --- compatibilite face / direction (_face_direction_compatibility) ---
    std::vector<std::uint8_t> compat(F * D);
    for (std::int64_t f = 0; f < F; ++f) {
        const double* n = in.face_normals + f * 3;
        double norm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (norm == 0.0) norm = 1.0;
        const double u[3] = {n[0] / norm, n[1] / norm, n[2] / norm};
        for (std::int64_t d = 0; d < D; ++d) {
            const double* dir = in.directions + d * 3;
            double c = u[0] * -dir[0] + u[1] * -dir[1] + u[2] * -dir[2];
            c = std::clamp(c, -1.0, 1.0);
            const double a = std::acos(c) * kRadToDeg;
            compat[f * D + d] = a < in.max_angle_deg ? 0 : 1;
        }
    }

    // --- adjacence de faces (build_face_adjacency), ordre CSR (ligne, colonne) ---
    std::map<std::pair<std::int64_t, std::int64_t>, std::vector<std::int64_t>> edge_faces;
    for (std::int64_t f = 0; f < F; ++f)
        for (int i = 0; i < 3; ++i) {
            const std::int64_t a = in.faces[f * 3 + i], b = in.faces[f * 3 + (i + 1) % 3];
            edge_faces[{std::min(a, b), std::max(a, b)}].push_back(f);
        }
    std::vector<std::vector<std::int64_t>> nbr(F);
    for (const auto& [edge, fs] : edge_faces) {
        if (fs.size() != 2) continue;
        nbr[fs[0]].push_back(fs[1]);
        nbr[fs[1]].push_back(fs[0]);
    }
    std::vector<std::int64_t> row, col;
    for (std::int64_t f = 0; f < F; ++f) {
        auto& v = nbr[f];
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());   // doublons sommes par scipy
        for (std::int64_t g : v) { row.push_back(f); col.push_back(g); }
    }
    const std::size_t E = row.size();

    // --- croissance de regions (_region_growing), tri STABLE ---
    std::vector<double> edge_risk(E);
    for (std::size_t e = 0; e < E; ++e) edge_risk[e] = in.risk[row[e]] + in.risk[col[e]];
    std::vector<std::size_t> order(E);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return -edge_risk[a] < -edge_risk[b]; });

    std::vector<std::int64_t> region(F);
    std::iota(region.begin(), region.end(), 0);
    std::unordered_map<std::int64_t, std::vector<std::uint8_t>> rcompat;
    std::unordered_map<std::int64_t, std::vector<double>> rresidual;
    const bool gate = in.merge_angle_deg >= 0.0;
    for (std::int64_t f = 0; f < F; ++f) {
        rcompat[f].assign(compat.begin() + f * D, compat.begin() + (f + 1) * D);
        if (gate) {
            auto& r = rresidual[f];
            r.resize(D);
            for (std::int64_t d = 0; d < D; ++d) r[d] = (compat[f * D + d] ? 0.0 : 1.0) * in.risk[f];
        }
    }
    auto find = [&](std::int64_t f) {
        while (region[f] != f) {
            region[f] = region[region[f]];
            f = region[f];
        }
        return f;
    };
    auto argmin = [&](const std::vector<double>& v) {
        std::int64_t best = 0;
        for (std::int64_t d = 1; d < static_cast<std::int64_t>(v.size()); ++d)
            if (v[d] < v[best]) best = d;
        return best;
    };
    std::vector<std::uint8_t> shared(D);
    for (std::size_t idx : order) {
        const std::int64_t r1 = find(row[idx]), r2 = find(col[idx]);
        if (r1 == r2) continue;
        const auto& c1 = rcompat[r1];
        const auto& c2 = rcompat[r2];
        bool any = false;
        for (std::int64_t d = 0; d < D; ++d) {
            shared[d] = c1[d] & c2[d];
            any = any || shared[d];
        }
        if (!any) continue;
        if (gate) {
            const auto& a = rresidual[r1];
            const auto& b = rresidual[r2];
            std::int64_t ic = -1;
            double best = std::numeric_limits<double>::infinity();
            for (std::int64_t d = 0; d < D; ++d) {
                const double v = shared[d] ? a[d] + b[d] : std::numeric_limits<double>::infinity();
                if (ic < 0 || v < best) { best = v; ic = d; }
            }
            const double* dc = in.directions + ic * 3;
            const double a1 = angle_deg(in.directions + argmin(a) * 3, dc);
            const double a2 = angle_deg(in.directions + argmin(b) * 3, dc);
            if (std::max(a1, a2) > in.merge_angle_deg) continue;
        }
        region[r2] = r1;
        rcompat[r1] = shared;
        if (gate) {
            auto& a = rresidual[r1];
            const auto& b = rresidual[r2];
            for (std::int64_t d = 0; d < D; ++d) a[d] = a[d] + b[d];
            rresidual.erase(r2);
        }
        rcompat.erase(r2);
    }
    std::vector<std::int64_t> roots(F);
    for (std::int64_t f = 0; f < F; ++f) roots[f] = find(f);
    std::vector<std::int64_t> unique_roots(roots);
    std::sort(unique_roots.begin(), unique_roots.end());
    unique_roots.erase(std::unique(unique_roots.begin(), unique_roots.end()), unique_roots.end());
    const std::int64_t R = static_cast<std::int64_t>(unique_roots.size());
    out.labels.resize(F);
    for (std::int64_t f = 0; f < F; ++f)
        out.labels[f] = std::lower_bound(unique_roots.begin(), unique_roots.end(), roots[f]) - unique_roots.begin();

    // --- direction par region (_select_region_directions) ---
    std::vector<std::vector<double>> residual(R, std::vector<double>(D, 0.0));
    for (std::int64_t f = 0; f < F; ++f) {
        auto& r = residual[out.labels[f]];
        for (std::int64_t d = 0; d < D; ++d) r[d] = r[d] + (compat[f * D + d] ? 0.0 : 1.0) * in.risk[f];
    }
    out.direction_index.resize(R);
    for (std::int64_t r = 0; r < R; ++r) out.direction_index[r] = argmin(residual[r]);

    // --- adjacence entre regions (networkx : premiere apparition) ---
    std::vector<std::vector<std::int64_t>> radj(R);
    for (std::size_t e = 0; e < E; ++e) {
        const std::int64_t a = out.labels[row[e]], b = out.labels[col[e]];
        if (a == b) continue;
        if (std::find(radj[a].begin(), radj[a].end(), b) == radj[a].end()) {
            radj[a].push_back(b);
            radj[b].push_back(a);
        }
    }

    // --- profondeurs (_compute_print_depths) : BFS multi-sources par composante ---
    std::vector<std::uint8_t> is_root(R, 0);
    for (std::int64_t f = 0; f < F; ++f)
        if (in.platform[f]) is_root[out.labels[f]] = 1;
    const double un = std::sqrt(in.up[0] * in.up[0] + in.up[1] * in.up[1] + in.up[2] * in.up[2]);
    const double up[3] = {in.up[0] / un, in.up[1] / un, in.up[2] / un};
    std::vector<double> lowest(R, std::numeric_limits<double>::infinity());
    for (std::int64_t f = 0; f < F; ++f) {
        double c[3];
        for (int k = 0; k < 3; ++k)
            c[k] = (in.vertices[in.faces[f * 3] * 3 + k] + in.vertices[in.faces[f * 3 + 1] * 3 + k] +
                    in.vertices[in.faces[f * 3 + 2] * 3 + k]) / 3.0;
        lowest[out.labels[f]] = std::min(lowest[out.labels[f]], c[0] * up[0] + c[1] * up[1] + c[2] * up[2]);
    }
    out.depth.assign(R, -1);
    std::vector<std::int64_t> component(R, -1);
    for (std::int64_t s = 0; s < R; ++s) {
        if (component[s] >= 0) continue;
        std::vector<std::int64_t> members{s};
        component[s] = s;
        for (std::size_t k = 0; k < members.size(); ++k)
            for (std::int64_t v : radj[members[k]])
                if (component[v] < 0) { component[v] = s; members.push_back(v); }
        std::deque<std::int64_t> queue;
        for (std::int64_t m : members)
            if (is_root[m]) { out.depth[m] = 0; queue.push_back(m); }
        if (queue.empty()) {
            std::int64_t best = members[0];
            for (std::int64_t m : members)
                if (lowest[m] < lowest[best] || (lowest[m] == lowest[best] && m < best)) best = m;
            out.depth[best] = 0;
            queue.push_back(best);
        }
        while (!queue.empty()) {
            const std::int64_t u = queue.front();
            queue.pop_front();
            for (std::int64_t v : radj[u])
                if (out.depth[v] < 0) { out.depth[v] = out.depth[u] + 1; queue.push_back(v); }
        }
    }

    // --- dependances, ordre de networkx (_build_dependency_graph) ---
    std::vector<std::vector<std::int64_t>> succ(R);
    std::vector<std::uint8_t> seen(R, 0);
    for (std::int64_t u = 0; u < R; ++u) {
        for (std::int64_t v : radj[u]) {
            if (seen[v]) continue;
            if (out.depth[u] < out.depth[v]) succ[u].push_back(v);
            else if (out.depth[v] < out.depth[u]) succ[v].push_back(u);
        }
        seen[u] = 1;
    }
    for (std::int64_t u = 0; u < R; ++u)
        for (std::int64_t v : succ[u]) { out.dependencies.push_back(u); out.dependencies.push_back(v); }
    return out;
}

}  // namespace apex
