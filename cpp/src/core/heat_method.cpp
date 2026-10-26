#include "apex/heat_method.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "apex/field_solver.hpp"

namespace apex {

namespace {

double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    const std::size_t n = v.size();
    std::nth_element(v.begin(), v.begin() + n / 2, v.end());
    const double hi = v[n / 2];
    if (n % 2 == 1) return hi;
    const double lo = *std::max_element(v.begin(), v.begin() + n / 2);
    return 0.5 * (lo + hi);  // np.median
}

struct Solved {
    std::vector<double> x;
    bool ok = false;
};

// Equivalent de `solve` du Python : direct sous le seuil, CG diagonal sinon,
// avec UN repli a max(20000, 10 * maxiter) iterations.
Solved solve(const SparsityPattern& pattern, const std::vector<double>& data,
             std::int64_t node_count, const std::vector<double>& b, const HeatParams& p) {
    Solved out;
    if (node_count < p.direct_threshold) {
        SolveResult r = solve_linear_system(pattern.indptr.data(), pattern.indices.data(),
                                            data.data(), node_count, b.data(),
                                            SolverKind::Direct, p.tolerance, p.max_iterations,
                                            nullptr);
        out.ok = r.converged;
        out.x = std::move(r.values);
        return out;
    }
    SolveResult r = solve_linear_system(pattern.indptr.data(), pattern.indices.data(),
                                        data.data(), node_count, b.data(), SolverKind::Iterative,
                                        p.tolerance, p.max_iterations, nullptr);
    if (!r.converged) {
        r = solve_linear_system(pattern.indptr.data(), pattern.indices.data(), data.data(),
                                node_count, b.data(), SolverKind::Iterative, p.tolerance,
                                std::max(20000, 10 * p.max_iterations), nullptr);
    }
    out.ok = r.converged;
    out.x = std::move(r.values);
    return out;
}

}  // namespace

HeatResult heat_method_distance(const double* tet_vertices, const std::int64_t* tetrahedra,
                                std::int64_t tet_count, std::int64_t node_count,
                                const double* shape_gradients, const double* volumes,
                                const std::int64_t* source_nodes, std::int64_t source_count,
                                const HeatParams& p) {
    HeatResult result;
    if (source_count == 0) {
        result.phi.assign(static_cast<std::size_t>(node_count),
                          std::numeric_limits<double>::infinity());
        result.converged = true;
        return result;
    }
    const auto N = static_cast<std::size_t>(node_count);

    // Rigidite L (volume comme poids) et masse condensee M (diagonale).
    const SparsityPattern pattern = build_fem_sparsity(tetrahedra, tet_count, node_count);
    std::vector<double> L(pattern.indices.size());
    assemble_gram_matrix(tetrahedra, shape_gradients, volumes, nullptr, 0, tet_count, pattern,
                         L.data());
    std::vector<double> mass(N, 0.0);
    for (std::int64_t t = 0; t < tet_count; ++t) {
        const double c = std::abs(volumes[t]) / 4.0;
        for (int k = 0; k < 4; ++k) mass[static_cast<std::size_t>(tetrahedra[4 * t + k])] += c;
    }

    // Temps de diffusion : t = facteur * h^2, h = arete (0,1) mediane.
    std::vector<double> lengths(static_cast<std::size_t>(tet_count));
    for (std::int64_t t = 0; t < tet_count; ++t) {
        const double* a = tet_vertices + 3 * tetrahedra[4 * t];
        const double* b = tet_vertices + 3 * tetrahedra[4 * t + 1];
        const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
        lengths[static_cast<std::size_t>(t)] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    double h = median(std::move(lengths));
    if (h == 0.0) h = 1.0;
    const double t_diff = p.time_factor * h * h;

    // (M + t L) u = source
    std::vector<double> A(L.size());
    for (std::size_t i = 0; i < L.size(); ++i) A[i] = t_diff * L[i];
    for (std::size_t row = 0; row < N; ++row) {
        for (std::int32_t s = pattern.indptr[row]; s < pattern.indptr[row + 1]; ++s) {
            if (static_cast<std::size_t>(pattern.indices[static_cast<std::size_t>(s)]) == row) {
                A[static_cast<std::size_t>(s)] += mass[row];
            }
        }
    }
    std::vector<double> source(N, 0.0);
    std::vector<char> is_source(N, 0);
    for (std::int64_t i = 0; i < source_count; ++i) {
        source[static_cast<std::size_t>(source_nodes[i])] = 1.0;
        is_source[static_cast<std::size_t>(source_nodes[i])] = 1;
    }
    Solved u = solve(pattern, A, node_count, source, p);
    if (!u.ok) return result;

    // Flux normalise par tetraedre.
    std::vector<double> flow(static_cast<std::size_t>(tet_count) * 3, 0.0);
    std::vector<char> usable(static_cast<std::size_t>(tet_count), 0);
    for (std::int64_t t = 0; t < tet_count; ++t) {
        const double* g = shape_gradients + 12 * t;
        double grad[3];
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int j = 0; j < 4; ++j) {
                s += g[i * 4 + j] * u.x[static_cast<std::size_t>(tetrahedra[4 * t + j])];
            }
            grad[i] = s;
        }
        const double n = std::sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
        if (n > 1e-12) {
            usable[static_cast<std::size_t>(t)] = 1;
            for (int i = 0; i < 3; ++i) flow[3 * static_cast<std::size_t>(t) + i] = -grad[i] / n;
        }
    }

    if (p.project_reachable) {
        const double un = std::sqrt(p.up[0] * p.up[0] + p.up[1] * p.up[1] + p.up[2] * p.up[2]);
        const double up[3] = {p.up[0] / un, p.up[1] / un, p.up[2] / un};
        for (std::int64_t t = 0; t < tet_count; ++t) {
            if (!usable[static_cast<std::size_t>(t)]) continue;
            double c[3] = {0, 0, 0};
            for (int k = 0; k < 4; ++k) {
                const double* v = tet_vertices + 3 * tetrahedra[4 * t + k];
                for (int i = 0; i < 3; ++i) c[i] += v[i];
            }
            for (int i = 0; i < 3; ++i) c[i] = c[i] / 4.0 - p.origin[i];
            const double along = c[0] * up[0] + c[1] * up[1] + c[2] * up[2];
            double radial[3];
            for (int i = 0; i < 3; ++i) radial[i] = c[i] - along * up[i];
            const double rn = std::sqrt(radial[0] * radial[0] + radial[1] * radial[1] +
                                        radial[2] * radial[2]);
            if (!(rn > 1e-6)) continue;  // sur l'axe : direction laissee intacte
            double* d = &flow[3 * static_cast<std::size_t>(t)];
            double d_r = (d[0] * radial[0] + d[1] * radial[1] + d[2] * radial[2]) / rn;
            const double d_z = d[0] * up[0] + d[1] * up[1] + d[2] * up[2];
            d_r = std::max(d_r, 0.0);
            double q[3];
            for (int i = 0; i < 3; ++i) q[i] = d_r * (radial[i] / rn) + d_z * up[i];
            double qn = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
            if (qn < 1e-9) {
                for (int i = 0; i < 3; ++i) q[i] = up[i];
                qn = 1.0;
            }
            for (int i = 0; i < 3; ++i) d[i] = q[i] / qn;
        }
    }

    // div(X) : rhs[noeud i] += V * sum_k G[k][i] * X[k]
    std::vector<double> rhs(N, 0.0);
    for (std::int64_t t = 0; t < tet_count; ++t) {
        const double* g = shape_gradients + 12 * t;
        const double* x = &flow[3 * static_cast<std::size_t>(t)];
        for (int i = 0; i < 4; ++i) {
            const double s = g[0 * 4 + i] * x[0] + g[1 * 4 + i] * x[1] + g[2 * 4 + i] * x[2];
            rhs[static_cast<std::size_t>(tetrahedra[4 * t + i])] += volumes[t] * s;
        }
    }

    // (L + ancrage sur la source) phi = rhs
    std::vector<double> B = L;
    for (std::size_t row = 0; row < N; ++row) {
        if (!is_source[row]) continue;
        for (std::int32_t s = pattern.indptr[row]; s < pattern.indptr[row + 1]; ++s) {
            if (static_cast<std::size_t>(pattern.indices[static_cast<std::size_t>(s)]) == row) {
                B[static_cast<std::size_t>(s)] += p.anchor_weight;
            }
        }
    }
    Solved phi = solve(pattern, B, node_count, rhs, p);
    if (!phi.ok) return result;

    double mean_all = 0.0, mean_src = 0.0, min_src = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < N; ++i) mean_all += phi.x[i];
    mean_all /= static_cast<double>(N);
    for (std::int64_t i = 0; i < source_count; ++i) {
        const double v = phi.x[static_cast<std::size_t>(source_nodes[i])];
        mean_src += v;
    }
    mean_src /= static_cast<double>(source_count);
    if (mean_src > mean_all) {
        for (double& v : phi.x) v = -v;
    }
    for (std::int64_t i = 0; i < source_count; ++i) {
        min_src = std::min(min_src, phi.x[static_cast<std::size_t>(source_nodes[i])]);
    }
    for (double& v : phi.x) v -= min_src;
    result.phi = std::move(phi.x);
    result.converged = true;
    return result;
}

}  // namespace apex
