// Retour de l'espace deforme vers la piece d'origine (voie Bird,
// pipeline/bird/deformed_slicing.py, `_parallel_barycentric`).
//
// Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer).
//
// Pour chaque point i et chacun de ses k tetraedres candidats c = cand[i, j] :
//   - coordonnees barycentriques du point dans le tetraedre DEFORME c ;
//   - score = plus petite coordonnee (NaN -> -inf, comme np.nan_to_num) ;
//   - point ramene = meme combinaison des sommets du tetraedre d'ORIGINE c.
// Les candidats (k plus proches centres) restent calcules par scipy cote Python :
// leurs indices sont donc identiques au bit.
//
// EQUIVALENCE : TOLERANCE DECLAREE, PAS LE BIT. Le Python resout le systeme 3x3
// par LAPACK (OpenBLAS, dgesv) ; une LU a pivot partiel transcrite ne coincide
// au bit que sur ~53 % des cas (mesure du 2026-10-13, Benchy passe 4) : OpenBLAS
// a ses propres noyaux. Ecart attendu : quelques ULP. La consequence se mesure en
// bout de chaine (seuils de validite a 1e-6, Viterbi, collisions, G-code).
//
// Systeme singulier (tetraedre plat) : le Python leve LinAlgError ; ici les
// coordonnees sont NaN, donc score -inf et point ramene NaN (candidat jamais
// retenu). Aucun cas mesure sur les modeles de reference.
#pragma once

#include <cstdint>

namespace apex {

// points (n,3) ; deformed (t,4,3) et original (t,4,3) : sommets des tetraedres
// deformes et d'origine ; cand (n,k) indices dans [0, t). Sorties : score (n,k),
// mapped (n,k,3). threads <= 0 : tous les coeurs.
void tet_barycentric_candidates(const double* points, std::int64_t n, const double* deformed,
                                const double* original, std::int64_t t, const std::int64_t* cand,
                                std::int64_t k, double* score, double* mapped, int threads = 0);

}  // namespace apex
