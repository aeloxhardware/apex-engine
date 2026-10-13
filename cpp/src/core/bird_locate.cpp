// Voir apex/bird_locate.hpp. Derive de S4 Slicer (Joshua Bird, GPL-3.0).
#include "apex/bird_locate.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>

#include "apex/internal/parallel_pool.hpp"

namespace apex {
namespace {

internal::ParallelPool& locate_pool() {
    static internal::ParallelPool* pool =
        new internal::ParallelPool(std::max(1u, std::thread::hardware_concurrency()));
    return *pool;
}

// Resout A x = b (A 3x3, lignes) par LU a pivot partiel, comme dgesv. Rend false
// si un pivot est nul (systeme singulier).
bool solve3(double a[3][3], double b[3], double x[3]) {
    for (int j = 0; j < 3; ++j) {
        int p = j;
        for (int i = j + 1; i < 3; ++i)
            if (std::fabs(a[i][j]) > std::fabs(a[p][j])) p = i;
        if (a[p][j] == 0.0) return false;
        if (p != j) {
            for (int c = 0; c < 3; ++c) std::swap(a[j][c], a[p][c]);
            std::swap(b[j], b[p]);
        }
        const double r = 1.0 / a[j][j];
        for (int i = j + 1; i < 3; ++i) {
            const double l = a[i][j] * r;
            for (int c = j + 1; c < 3; ++c) a[i][c] -= l * a[j][c];
            b[i] -= l * b[j];
        }
    }
    for (int i = 2; i >= 0; --i) {
        double s = b[i];
        for (int c = i + 1; c < 3; ++c) s -= a[i][c] * x[c];
        x[i] = s / a[i][i];
    }
    return true;
}

void one(const double* p, const double* d, const double* o, double* score, double* mapped) {
    double a[3][3], b[3], lam[3];
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) a[r][c] = d[(c + 1) * 3 + r] - d[r];
        b[r] = p[r] - d[r];
    }
    double bary[4];
    if (solve3(a, b, lam)) {
        bary[1] = lam[0], bary[2] = lam[1], bary[3] = lam[2];
        bary[0] = 1.0 - ((lam[0] + lam[1]) + lam[2]);
    } else {
        bary[0] = bary[1] = bary[2] = bary[3] = std::numeric_limits<double>::quiet_NaN();
    }
    double m = bary[0];
    bool nan = std::isnan(m);
    for (int v = 1; v < 4; ++v) {
        nan = nan || std::isnan(bary[v]);
        m = std::min(m, bary[v]);
    }
    *score = nan ? -std::numeric_limits<double>::infinity() : m;
    for (int c = 0; c < 3; ++c)
        mapped[c] = bary[0] * o[c] + bary[1] * o[3 + c] + bary[2] * o[6 + c] + bary[3] * o[9 + c];
}

}  // namespace

void tet_barycentric_candidates(const double* points, std::int64_t n, const double* deformed,
                                const double* original, std::int64_t t, const std::int64_t* cand,
                                std::int64_t k, double* score, double* mapped, int threads) {
    for (std::int64_t i = 0; i < n * k; ++i)
        if (cand[i] < 0 || cand[i] >= t) throw std::out_of_range("indice de tetraedre candidat hors limites");
    auto block = [&](std::int64_t b, std::int64_t e) {
        for (std::int64_t i = b; i < e; ++i)
            for (std::int64_t j = 0; j < k; ++j) {
                const std::int64_t c = cand[i * k + j];
                one(points + i * 3, deformed + c * 12, original + c * 12, score + i * k + j,
                    mapped + (i * k + j) * 3);
            }
    };
    if (threads == 1) {
        block(0, n);
        return;
    }
    locate_pool().run(n, block, 2048);
}

}  // namespace apex
