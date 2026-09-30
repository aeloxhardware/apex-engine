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

#include "apex/field_solver.hpp"

namespace apex {

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

    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper,
                             Eigen::DiagonalPreconditioner<double>>
        solver;
    solver.setTolerance(tolerance);
    solver.setMaxIterations(max_iterations);
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        return result;
    }

    if (initial_guess != nullptr) {
        const Eigen::Map<const Eigen::VectorXd> guess(initial_guess, node_count);
        x = solver.solveWithGuess(rhs, guess);
    } else {
        x = solver.solve(rhs);
    }

    result.iterations = static_cast<int>(solver.iterations());
    result.error = solver.error();
    result.converged = solver.info() == Eigen::Success;
    return result;
}

}  // namespace apex
