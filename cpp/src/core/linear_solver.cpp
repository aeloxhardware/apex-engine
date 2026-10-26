// Resolution du systeme lineaire de step4, via Eigen.
//
// Fichier SEPARE de field_solver.cpp uniquement pour le temps de
// compilation : les en-tetes Eigen sont lourds, et seule cette unite doit
// les payer.
//
// EQUIVALENCE : c'est le SEUL portage du projet qui ne peut pas etre
// bit-a-bit, et c'est assume. `scipy.sparse.linalg.cg` et le gradient
// conjugue d'Eigen resolvent le meme systeme par le meme algorithme, avec le
// meme preconditionneur diagonal et le meme critere d'arret relatif - mais
// leurs produits scalaires n'accumulent pas dans le meme ordre (numpy passe
// par BLAS, Eigen par ses propres noyaux vectorises). A `rtol = 1e-8`, les
// deux champs concordent a ~1e-8 en relatif ; exiger le bit reviendrait a
// reimplementer BLAS.
//
// Correspondance des deux chemins, verifiee ligne a ligne cote Python :
//   - `node_count < direct_solve_node_threshold` -> solveur DIRECT. Le
//     Python appelle `spsolve` (SuperLU, matrice generale) ; le systeme
//     etant symetrique defini positif par construction (rigidite FEM +
//     penalites quadratiques), `SimplicialLDLT` est l'equivalent adapte et
//     n'exige aucune hypothese que le Python ne fasse deja.
//   - sinon -> gradient conjugue preconditionne par la DIAGONALE
//     (`DiagonalPreconditioner` = le `x / diagonal` du Python), `rtol` et
//     `maxiter` repris tels quels.
//
// Le Python protege sa diagonale (`diagonal[diagonal == 0.0] = 1.0`) ; le
// preconditionneur d'Eigen fait la meme chose en interne (il neutralise les
// entrees nulles), donc aucune divergence de comportement sur une diagonale
// degeneree.

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include <algorithm>
#include <cmath>
#include <thread>

#include "apex/field_solver.hpp"
#include "apex/internal/parallel_pool.hpp"

namespace apex {

namespace {

// Gradient conjugue a preconditionneur de Jacobi, parallele ET reproductible :
// produits matrice-vecteur par lignes (chaque ligne est sommee seule, dans
// l'ordre du CSR), produits scalaires par blocs de taille FIXE dont les sommes
// partielles sont additionnees dans l'ordre -- le resultat ne depend donc ni du
// nombre de fils ni de leur ordonnancement. Meme algorithme et meme critere
// d'arret que scipy.sparse.linalg.cg(M=diag) : ||r|| <= rtol * ||b||.
constexpr std::int64_t kBlock = 4096;

internal::ParallelPool& solver_pool() {
    static internal::ParallelPool* pool =
        new internal::ParallelPool(std::max(1u, std::thread::hardware_concurrency()));
    return *pool;
}

struct ParallelCg {
    const std::int32_t* indptr;
    const std::int32_t* indices;
    const double* data;
    std::int64_t n;
    std::vector<double> partial;

    template <class F>
    void for_rows(F&& f) {
        solver_pool().run(n, [&](std::int64_t b, std::int64_t e) { f(b, e); }, 8192);
    }

    void spmv(const double* x, double* y) {
        for_rows([&](std::int64_t b, std::int64_t e) {
            for (std::int64_t row = b; row < e; ++row) {
                double sum = 0.0;
                for (std::int32_t k = indptr[row]; k < indptr[row + 1]; ++k) sum += data[k] * x[indices[k]];
                y[row] = sum;
            }
        });
    }

    double dot(const double* a, const double* b) {
        const std::int64_t blocks = (n + kBlock - 1) / kBlock;
        partial.assign(static_cast<std::size_t>(blocks), 0.0);
        solver_pool().run(blocks, [&](std::int64_t b0, std::int64_t b1) {
            for (std::int64_t blk = b0; blk < b1; ++blk) {
                const std::int64_t lo = blk * kBlock, hi = std::min(n, lo + kBlock);
                double sum = 0.0;
                for (std::int64_t i = lo; i < hi; ++i) sum += a[i] * b[i];
                partial[static_cast<std::size_t>(blk)] = sum;
            }
        }, 4);
        double total = 0.0;
        for (const double v : partial) total += v;
        return total;
    }
};

SolveResult parallel_jacobi_cg(const std::int32_t* indptr, const std::int32_t* indices,
                               const double* data, std::int64_t n, const double* b,
                               double tolerance, int max_iterations, const double* guess) {
    SolveResult result;
    result.values.assign(static_cast<std::size_t>(n), 0.0);
    ParallelCg cg{indptr, indices, data, n, {}};
    std::vector<double> inv_diag(static_cast<std::size_t>(n), 1.0);
    for (std::int64_t row = 0; row < n; ++row) {
        for (std::int32_t k = indptr[row]; k < indptr[row + 1]; ++k) {
            if (indices[k] == row && data[k] != 0.0) inv_diag[static_cast<std::size_t>(row)] = 1.0 / data[k];
        }
    }
    std::vector<double>& x = result.values;
    std::vector<double> r(static_cast<std::size_t>(n)), z(static_cast<std::size_t>(n)),
        p(static_cast<std::size_t>(n)), q(static_cast<std::size_t>(n));
    if (guess != nullptr) {
        std::copy(guess, guess + n, x.begin());
        cg.spmv(x.data(), q.data());
        cg.for_rows([&](std::int64_t lo, std::int64_t hi) {
            for (std::int64_t i = lo; i < hi; ++i) r[static_cast<std::size_t>(i)] = b[i] - q[static_cast<std::size_t>(i)];
        });
    } else {
        std::copy(b, b + n, r.begin());
    }
    const double bnorm = std::sqrt(cg.dot(b, b));
    if (bnorm == 0.0) {
        result.converged = true;
        return result;
    }
    const double threshold = tolerance * bnorm;
    double rz_old = 0.0;
    for (int it = 0; it < max_iterations; ++it) {
        const double rnorm = std::sqrt(cg.dot(r.data(), r.data()));
        result.iterations = it;
        result.error = rnorm / bnorm;
        if (rnorm <= threshold) {
            result.converged = true;
            return result;
        }
        cg.for_rows([&](std::int64_t lo, std::int64_t hi) {
            for (std::int64_t i = lo; i < hi; ++i) z[static_cast<std::size_t>(i)] = r[static_cast<std::size_t>(i)] * inv_diag[static_cast<std::size_t>(i)];
        });
        const double rz = cg.dot(r.data(), z.data());
        if (it == 0) {
            p = z;
        } else {
            const double beta = rz / rz_old;
            cg.for_rows([&](std::int64_t lo, std::int64_t hi) {
                for (std::int64_t i = lo; i < hi; ++i) p[static_cast<std::size_t>(i)] = z[static_cast<std::size_t>(i)] + beta * p[static_cast<std::size_t>(i)];
            });
        }
        rz_old = rz;
        cg.spmv(p.data(), q.data());
        const double pq = cg.dot(p.data(), q.data());
        if (pq == 0.0) break;
        const double alpha = rz / pq;
        cg.for_rows([&](std::int64_t lo, std::int64_t hi) {
            for (std::int64_t i = lo; i < hi; ++i) {
                x[static_cast<std::size_t>(i)] += alpha * p[static_cast<std::size_t>(i)];
                r[static_cast<std::size_t>(i)] -= alpha * q[static_cast<std::size_t>(i)];
            }
        });
    }
    const double rnorm = std::sqrt(cg.dot(r.data(), r.data()));
    result.error = rnorm / bnorm;
    result.converged = rnorm <= threshold;
    return result;
}

}  // namespace

SolveResult solve_linear_system(const std::int32_t* indptr, const std::int32_t* indices,
                                const double* data, std::int64_t node_count, const double* b,
                                SolverKind kind, double tolerance, int max_iterations,
                                const double* initial_guess) {
    using SparseMatrix = Eigen::SparseMatrix<double>;

    // Vue SANS COPIE sur le CSR fourni par scipy. `Map` impose que les
    // indices soient tries par ligne, ce que garantit `.tocsr()` cote Python
    // (csr_matrix trie ses indices a la construction).
    const Eigen::Map<const SparseMatrix> A(node_count, node_count, indptr[node_count], indptr,
                                           indices, data);
    const Eigen::Map<const Eigen::VectorXd> rhs(b, node_count);

    SolveResult result;
    result.values.assign(static_cast<std::size_t>(node_count), 0.0);
    Eigen::Map<Eigen::VectorXd> x(result.values.data(), node_count);

    if (kind == SolverKind::Direct) {
        Eigen::SimplicialLDLT<SparseMatrix> solver;
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            return result;  // converged reste false : l'appelant retombe sur Python
        }
        x = solver.solve(rhs);
        result.converged = solver.info() == Eigen::Success;
        return result;
    }

    // Iteratif : CG de Jacobi parallele et reproductible (voir ParallelCg).
    // (L'ancien CG d'Eigen, mono-fil et mesure a x1,01 contre scipy, est retire.)
    return parallel_jacobi_cg(indptr, indices, data, node_count, b, tolerance, max_iterations, initial_guess);
}

}  // namespace apex
