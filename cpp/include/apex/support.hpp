// Test de support d'un point contre des capsules - portage de
// `core/geometry_utils.py::point_has_support_along_direction`.
//
// POURQUOI CE MODULE EST LE TROISIEME PORTE : apres les pas 1 (index spatial)
// et 2 (cinematique inverse), le re-profilage de death.stl @ 0.8 la place en
// tete des postes PORTABLES avec 11,97 s de temps propre pour 273 146 appels,
// plus ~10 s de dispatch numpy qu'elle declenche (`cross`,
// `normalize_axis_tuple`, `moveaxis`) - soit ~22 s des 98,2 s profilees. Le
// seul poste plus gros est TetGen lui-meme (16,6 s), qui est deja natif et
// incompressible.
//
// LE VRAI LEVIER EST L'ARRET ANTICIPE, PAS LE LANGAGE. La fonction ne rend
// qu'un BOOLEEN ("existe-t-il un appui ?") et execute pourtant ~15 operations
// numpy sur les ~690 candidats a chaque appel, sans jamais court-circuiter.
// Mesure a l'appui : la grille rend 690,6 candidats uniques par requete en
// moyenne, et ce nombre ne depend PAS de la structure spatiale (un BVH
// parfait rendrait les memes). C'est donc ici, et non dans l'index, que se
// trouve le plancher identifie au pas 1.
//
// EQUIVALENCE : meme construction de base orthonormee que
// `_perpendicular_basis` (reference fixe [1,0,0], repli [0,1,0] au-dela de
// 0.9), meme projection scalaire clampee, memes comparaisons. L'arret
// anticipe ne change pas le resultat : `np.any(...)` est vrai des qu'un
// candidat satisfait les trois conditions, donc s'arreter au premier rend
// exactement la meme reponse.

#pragma once

#include <cstdint>

#include "apex/types.hpp"

namespace apex {

// `position` a-t-il un appui parmi les segments `starts[i]`-`ends[i]` le long
// de `direction_up` ? Un candidat supporte le point si son point le plus
// proche est STRICTEMENT en dessous (`delta_up > 1e-9`), a moins de
// `vertical_height` en dessous, et a moins de `horizontal_radius`
// transversalement.
bool point_has_support_along_direction(const Vec3& position, const Vec3& direction_up,
                                       const double* starts, const double* ends,
                                       std::int64_t candidate_count, double horizontal_radius,
                                       double vertical_height);

// Variante multi-directions : les MEMES candidats sont testes pour chaque
// direction, avec court-circuit des la premiere qui donne un appui.
// Transcription de `step9b._has_support_in_any_direction` et de la boucle de
// `step5._layer_support_mask`, qui font exactement cela cote Python - mais en
// repayant tout le cout vectoriel a chaque direction.
bool point_has_support_in_any_direction(const Vec3& position, const double* directions,
                                        std::int64_t direction_count, const double* starts,
                                        const double* ends, std::int64_t candidate_count,
                                        double horizontal_radius, double vertical_height);

}  // namespace apex
