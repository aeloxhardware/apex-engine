// Interpolation inverse-distance 2D - portage de
// `pipeline/step6_toolpath2d.py::_interpolate_out_of_plane`.
//
// CONTEXTE. Apres le pas 9, le reprofilage du pipeline (step1->step8,
// death.stl, resolution 0.8, C++ actif pas 1-9) donne 43,1 s profilees dont
// 22,4 s d'attente TetGen (52 %, sous-processus, plancher dur). Hors TetGen,
// PLUS AUCUN poste Python isole ne depasse 2,1 s - et celui-ci est le
// dernier : `_interpolate_out_of_plane`, 2,10 s de temps propre sur 2853
// appels. C'est donc le dernier portage "cible par le profil" ; ce qui
// restera ensuite est du portage SYSTEMATIQUE (etapes entieres a 0,1-0,7 s
// par poste), pas de l'optimisation guidee.
//
// CE QUE FAIT LA FONCTION : ponderation inverse-distance `1/d^power` d'un
// scalaire (l'offset hors-plan d'une couche) en des points 2D GENERES
// (parois offsettees, infill) a partir des offsets REELS connus sur les
// boucles de contour. Le cout vient de la matrice de distances complete
// (K requetes x M sources), materialisee cote numpy.
//
// EQUIVALENCE : TOLERANCE SERREE, PAS LE BIT - mesure, pas supposee.
// Il n'y a ici ni solve ni appel BLAS (contrairement aux pas 7 et 9),
// uniquement des operations elementaires ; la transcription est donc
// terme a terme. Mais les deux reductions (`weights.sum(axis=1)` et
// `(weights * offsets).sum(axis=1)`) sont accumulees SEQUENTIELLEMENT ici,
// la ou numpy peut sommer par paires/vectorise. MESURE sur 200 tirages
// aleatoires (tailles et coincidences variees) : ecart max **8,88e-16**,
// et 36/200 cas seulement bit-identiques - exactement le meme phenomene
// qu'au pas 9 pour `generalized_winding_number`. Un ULP de bruit, sans
// consequence : cet offset est une position hors-plan en millimetres, dont
// les seuils d'usage en aval sont ~1e-9 (voir `_COINCIDENT_EPS`) a 1e-3.
//
// GAIN mesure : x2,8 sur une taille realiste (400 requetes x 300 sources,
// 300 repetitions : 1,64 s -> 0,58 s). Il vient de ne jamais materialiser
// la matrice de distances (K,M) que numpy construit en entier.
//
// Les trois cas limites du Python sont transcrits tels quels, car ils ne
// sont PAS cosmetiques :
//   - `distances < 1e-9` (_COINCIDENT_EPS) : `0**-power` exploserait ; la
//     ligne concernee rend la valeur source EXACTE (celle du plus proche) ;
//   - `total <= 0` : sous-debordement numerique sur de tres grandes
//     distances - repli sur la source la plus proche, jamais une division
//     par zero ;
//   - `argmin` sur les distances : premier minimum en cas d'egalite (meme
//     convention que numpy), d'ou une comparaison STRICTE dans la boucle.

#pragma once

#include <cstdint>

namespace apex {

// `query` : K*2 doubles (u,v). `sources` : M*2. `offsets` : M.
// `out` : K, rempli par l'appelant.
void interpolate_out_of_plane(const double* query, std::int64_t query_count,
                              const double* sources, const double* offsets,
                              std::int64_t source_count, double power, double* out);

}  // namespace apex
