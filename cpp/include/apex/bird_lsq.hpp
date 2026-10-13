// Moindres carres non lineaires de la voie Bird (pipeline/bird/deformation.py).
//
// Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer) :
// les deux problemes (rotation par tetraedre, deformation du maillage) sont ceux
// de son notebook. Le solveur est une TRANSCRIPTION de scipy.optimize.least_squares
// (method="trf", sans bornes, tr_solver="lsmr", regularize=True) et de
// scipy.sparse.linalg.lsmr (atol = btol = 1e-6, conlim = 1e8, maxiter = min(m, n)).
//
// EQUIVALENCE : TOLERANCE DECLAREE, PAS LE BIT. La jacobienne est arrondie en
// float32 (comme dans l'original, `.astype(np.float32)`) et l'echelle `jac` est
// sommee en float32 comme scipy ; la QR a deux colonnes et les racines du
// quartique (np.roots) passent par LAPACK cote Python, par une QR de Householder et
// les valeurs propres d'Eigen ici. Le resultat est compare sur la deformation finale
// (ecart de positions) et la decision de la boucle (zones de repli).
//
// Interface : la jacobienne a une STRUCTURE FIXE (CSR, colonnes triees) ; seules
// ses valeurs changent d'une evaluation a l'autre.
#pragma once

#include <cstdint>
#include <vector>

namespace apex {

// Compteurs du dernier appel a bird_lsq_deform / bird_lsq_rotation (mesure de la
// phase 0 du plan Bird C++/Vulkan). Aucun effet sur le calcul.
struct LsqProfile {
    double total_s = 0.0, matvec_s = 0.0, rmatvec_s = 0.0, fun_s = 0.0, jac_s = 0.0;
    std::int64_t matvec_calls = 0, rmatvec_calls = 0, fun_calls = 0, jac_calls = 0;
    std::int64_t lsmr_calls = 0, lsmr_iterations = 0, nfev = 0;
    std::int64_t rows = 0, cols = 0, nnz = 0;
};
LsqProfile bird_lsq_last_profile();

// Nombre de fils des produits matrice-vecteur, residus et jacobienne (0 = tous les
// coeurs, 1 = sequentiel). Resultat identique au bit quel que soit ce nombre.
// `min_parallel` : taille de boucle en dessous de laquelle on reste sur un bloc
// (4096 par defaut ; 1 dans les tests, pour exercer le chemin parallele).
void bird_lsq_set_threads(int threads, std::int64_t min_parallel = 4096);

struct LsqResult {
    std::vector<double> x;
    double cost = 0.0;
    int nfev = 0;
    int status = 0;           // codes de scipy : 0 plafond, 1 gtol, 2 ftol, 3 xtol, 4 ftol+xtol
};

// Deformation (cellule 7 de Bird) : residu par tetraedre ||N V_c - T_c||^2,
// T_c = R_c N V0_c (fourni : old_transformed (T,3,4)). Noeuds epingles :
// column[i] < 0 ; inconnues = base_vertices des noeuds libres, dans l'ordre de column.
LsqResult bird_lsq_deform(const double* base_vertices, std::int64_t n_points,
                          const std::int64_t* cells, std::int64_t n_cells,
                          const std::int64_t* column, std::int64_t n_free,
                          const double* old_transformed, int max_nfev,
                          double ftol, double xtol, double gtol);

// Rotations (optimize_rotations) : residus W (f[c1] - f[c2])^2 pour chaque paire
// (pairs (P,2)), puis (f[t] - initial[t])^2 pour chaque cible ; cellules figees a 0 :
// column[c] < 0. x0 = 0, x_scale = 1.
LsqResult bird_lsq_rotation(const std::int64_t* pairs, std::int64_t n_pairs,
                            const std::int64_t* targets, std::int64_t n_targets,
                            const double* initial, std::int64_t n_cells,
                            const std::int64_t* column, std::int64_t n_free,
                            double weight, int max_nfev, double ftol, double xtol, double gtol);

}  // namespace apex

namespace apex {

// Retour d'une polyligne de l'espace deforme vers la piece (deformed_slicing.
// DeformedToOriginal.map_polyline) : programmation dynamique (Viterbi) qui choisit,
// pour chaque point, l'un de ses k candidats (`mapped` (n,k,3), `valid` (n,k)) de
// sorte que la somme des distances entre points ramenes consecutifs soit minimale.
// Transcription terme a terme (norme = sqrt de la somme x->y->z, premier minimum
// en cas d'egalite) : IDENTIQUE AU BIT a la version numpy.
void viterbi_choice(const double* mapped, const std::uint8_t* valid, std::int64_t n, std::int64_t k,
                    std::int64_t* choice);

}  // namespace apex
