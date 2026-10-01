// Voir apex/overhang.hpp : transcription de pipeline/step2_overhang.py.

#include "apex/overhang.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>

#include <Eigen/Dense>

#include "apex/bvh.hpp"

namespace apex {

namespace {

constexpr double kRadToDeg = 180.0 / 3.14159265358979323846;

std::int64_t find_root(std::vector<std::int64_t>& parent, std::int64_t i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

}  // namespace

OverhangResult compute_overhang_map(const double* V, std::int64_t vertex_count, const std::int64_t* F,
                                    const double* N, std::int64_t face_count, const OverhangParams& P) {
    (void)vertex_count;
    OverhangResult out;
    out.risk_score.resize(face_count);
    out.is_overhang.resize(face_count);
    out.is_bridge.assign(face_count, 0);
    out.span_direction.assign(face_count * 3, std::numeric_limits<double>::quiet_NaN());

    // --- angles a la verticale descendante (_compute_overhang_angles) ---
    const double up_norm = std::sqrt(P.up[0] * P.up[0] + P.up[1] * P.up[1] + P.up[2] * P.up[2]);
    const double down[3] = {-P.up[0] / up_norm, -P.up[1] / up_norm, -P.up[2] / up_norm};
    std::vector<double> angle(face_count);
    for (std::int64_t f = 0; f < face_count; ++f) {
        const double* n = N + f * 3;
        double norm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (norm == 0.0) norm = 1.0;
        double dot = (n[0] / norm) * down[0] + (n[1] / norm) * down[1] + (n[2] / norm) * down[2];
        dot = std::clamp(dot, -1.0, 1.0);
        angle[f] = std::acos(dot) * kRadToDeg;
    }

    // --- surplomb et score de risque (sigmoide centree sur le seuil) ---
    const double threshold = P.overhang_angle_threshold_deg;
    const double width = std::max(threshold * 0.25, 1e-3);
    for (std::int64_t f = 0; f < face_count; ++f) {
        out.is_overhang[f] = angle[f] < threshold ? 1 : 0;
        const double x = (threshold - angle[f]) / width;
        out.risk_score[f] = 1.0 / (1.0 + std::exp(-x));
    }

    // --- ponts : rayons verticaux depuis une grille barycentrique ---
    std::vector<std::int64_t> candidates;
    for (std::int64_t f = 0; f < face_count; ++f) {
        if (angle[f] < P.bridge_horizontal_angle_tolerance_deg) candidates.push_back(f);
    }
    if (!candidates.empty()) {
        std::vector<double> triangles(face_count * 9);
        for (std::int64_t f = 0; f < face_count; ++f)
            for (int k = 0; k < 3; ++k)
                for (int c = 0; c < 3; ++c) triangles[f * 9 + k * 3 + c] = V[F[f * 3 + k] * 3 + c];
        TriangleBVH bvh(triangles.data(), N, face_count);

        const int n = std::max(P.bridge_detection_ray_density, 1);
        std::vector<std::array<double, 3>> bary;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n - i; ++j) bary.push_back({double(i) / n, double(j) / n, double(n - i - j) / n});
        const double eps = 1e-6;
        std::vector<double> origins(bary.size() * 3), directions(bary.size() * 3, 0.0);
        for (std::size_t p = 0; p < bary.size(); ++p) directions[p * 3 + 2] = -1.0;

        for (std::int64_t f : candidates) {
            const double* a = V + F[f * 3 + 0] * 3;
            const double* b = V + F[f * 3 + 1] * 3;
            const double* c = V + F[f * 3 + 2] * 3;
            for (std::size_t p = 0; p < bary.size(); ++p) {
                for (int k = 0; k < 3; ++k) {
                    const double point = bary[p][0] * a[k] + bary[p][1] * b[k] + bary[p][2] * c[k];
                    origins[p * 3 + k] = point + (k == 2 ? -eps : 0.0);
                }
            }
            const auto hits = bvh.intersects_location(origins.data(), directions.data(),
                                                      static_cast<std::int64_t>(bary.size()));
            std::vector<double> best(bary.size(), std::numeric_limits<double>::infinity());
            for (const auto& h : hits) {
                const double* o = origins.data() + h.ray_index * 3;
                const double d = std::sqrt((h.location.x - o[0]) * (h.location.x - o[0]) +
                                           (h.location.y - o[1]) * (h.location.y - o[1]) +
                                           (h.location.z - o[2]) * (h.location.z - o[2]));
                best[h.ray_index] = std::min(best[h.ray_index], d);
            }
            std::size_t unsupported = 0;
            for (double d : best) unsupported += (d > P.bridge_max_gap_distance) ? 1 : 0;   // inf compris
            const double ratio = double(unsupported) / double(bary.size());
            if (ratio >= P.bridge_support_ratio_threshold) out.is_bridge[f] = 1;
        }
    }

    // --- composantes connexes (graphe de TOUTES les faces) et axe principal ---
    bool any_bridge = std::any_of(out.is_bridge.begin(), out.is_bridge.end(), [](std::uint8_t v) { return v; });
    if (!any_bridge) return out;
    std::map<std::pair<std::int64_t, std::int64_t>, std::vector<std::int64_t>> edge_faces;
    for (std::int64_t f = 0; f < face_count; ++f)
        for (int i = 0; i < 3; ++i) {
            std::int64_t a = F[f * 3 + i], b = F[f * 3 + (i + 1) % 3];
            edge_faces[{std::min(a, b), std::max(a, b)}].push_back(f);
        }
    std::vector<std::int64_t> parent(face_count);
    std::iota(parent.begin(), parent.end(), 0);
    for (const auto& [edge, fs] : edge_faces) {
        if (fs.size() != 2) continue;
        const std::int64_t ra = find_root(parent, fs[0]), rb = find_root(parent, fs[1]);
        if (ra != rb) parent[ra] = rb;
    }
    std::map<std::int64_t, std::vector<std::int64_t>> components;
    for (std::int64_t f = 0; f < face_count; ++f)
        if (out.is_bridge[f]) components[find_root(parent, f)].push_back(f);

    for (const auto& [root, fs] : components) {
        std::vector<std::int64_t> verts;
        for (std::int64_t f : fs)
            for (int k = 0; k < 3; ++k) verts.push_back(F[f * 3 + k]);
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();
        for (std::int64_t v : verts) mean += Eigen::Vector3d(V[v * 3], V[v * 3 + 1], V[v * 3 + 2]);
        mean /= double(verts.size());
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        for (std::int64_t v : verts) {
            const Eigen::Vector3d d = Eigen::Vector3d(V[v * 3], V[v * 3 + 1], V[v * 3 + 2]) - mean;
            cov += d * d.transpose();
        }
        Eigen::Vector3d axis(1.0, 0.0, 0.0);
        if (verts.size() >= 2) {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
            if (solver.info() == Eigen::Success) axis = solver.eigenvectors().col(2);   // plus grande valeur propre
        }
        // signe canonique, comme le Python : plus grande composante positive
        int big = 0;
        for (int k = 1; k < 3; ++k)
            if (std::fabs(axis[k]) > std::fabs(axis[big])) big = k;
        if (axis[big] < 0) axis = -axis;
        for (std::int64_t f : fs)
            for (int k = 0; k < 3; ++k) out.span_direction[f * 3 + k] = axis[k];
    }
    return out;
}

}  // namespace apex
