// Cinematique inverse N axes - portage de la chaine analytique de
// `core/geometry_utils.py` et de la boucle point-par-point de
// `pipeline/step8_machine_toolpath.py::run`.
//
// POURQUOI CE MODULE EST LE DEUXIEME PORTE : apres le portage de l'index
// spatial (pas 1), le re-profilage de death.stl @ 0.8 place step8 en tete
// avec 86,1 s sur 168,0 s profilees (51 %). Les fonctions dominantes sont
// toutes ici : `_rodrigues_rotate_single` (13,7 s / 5,74 M appels),
// `_resolve_axis_pair_analytic` (7,5 s / 298 k), `_angle_between`
// (4,1 s / 5,4 M) - plus ~25 s de PUR DISPATCH NUMPY (`cross`,
// `normalize_axis_tuple`, `moveaxis`, `norm`) sur des vecteurs de dimension 3.
// C'est le cas d'ecole ou l'interpreteur coute plus cher que le calcul.
//
// LA BOUCLE N'EST PAS PARALLELISABLE TELLE QUELLE, contrairement a ce que
// supposait le plan de portage : `previous_angles` porte une dependance
// sequentielle reelle d'un point au suivant (repliement des axes internes
// axes[2:], repli en cas de singularite locale, et heritage d'angle a la
// singularite). Le portage est donc SEQUENTIEL, et le gain vient du seul
// changement de langage. Une parallelisation exigerait de casser cette
// chaine, ce qui est un changement de COMPORTEMENT, pas de langage - hors
// de portee d'un portage a equivalence stricte.
//
// REPLI GRILLE, et pourquoi ce module est REPRENABLE.
// `_resolve_axis_pair_analytic` peut echouer a bracketer une racine. Mesure
// sur la cinematique reelle (C tour complet + B [0,185]) : 0 echec sur 2000
// directions aleatoires, MAIS echec systematique quand la cible est colineaire
// a l'axe externe (typiquement +Z, un depot vertical). Le residu du bracketing
// y est une TANGENCE a zero, jamais un changement de signe - le premier point
// d'une trajectoire a fond plat tombe pile dessus.
//
// La grille historique n'est deliberement PAS portee : elle sert 0,001 % des
// points, la porter serait du code a maintenir et une occasion de diverger
// pour rien. Ce module s'arrete au point fautif, remplit tout ce qui precede,
// et rend son index. L'appelant Python resout CE point avec son implementation
// d'origine, puis rappelle ce module a partir du suivant en passant les angles
// obtenus via `initial_previous_angles`.
//
// La reprise est ce qui rend le mecanisme viable : une premiere version
// abandonnait le lot ENTIER au premier repli, si bien qu'UN point sur 297 994
// annulait tout le gain (mesure : aucune acceleration). La chaine
// `previous_angles` est preservee exactement, donc l'equivalence l'est aussi.

#pragma once

#include <cstdint>
#include <vector>

#include "apex/types.hpp"

namespace apex {

// Equivalent de core.mesh_types.RotationAxisSpec.
struct AxisSpec {
    Vec3 axis;      // pas necessairement unitaire : normalise a l'usage
    double lo = 0.0;   // angle_range_deg[0]
    double hi = 0.0;   // angle_range_deg[1]
};

// Repli par grille + affinement local (core.geometry_utils.composed_rotation_angles_and_directions
// + refine_axis_angles_local_search), applique par le natif quand la resolution analytique echoue.
// `samples_per_axis == 0` : repli non porte, le natif rend la main a Python (fallback_index).
struct GridFallback {
    int samples_per_axis = 0;
    double refine_multiplier = 0.0;   // Step8Params.ik_refine_step_multiplier
    int refine_iterations = 0;        // Step8Params.ik_refine_iterations
};

struct IKBatchResult {
    std::vector<double> angles;     // point_count * axes.size(), degres
    std::vector<double> residuals;  // point_count, RADIANS (meme contrat que le Python)
    // Index du premier point exigeant le repli grille, -1 si aucun.
    std::int64_t fallback_index = -1;
    std::int64_t grid_fallbacks = 0;   // points resolus par le repli grille natif
};

// Constantes reprises telles quelles de core/geometry_utils.py - toute
// divergence casserait l'equivalence numerique.
inline constexpr double kParallelEps = 1e-9;
inline constexpr double kAxesParallelEps = 1e-9;
inline constexpr double kBisectionToleranceDeg = 1e-5;
inline constexpr int kBisectionMaxIterations = 60;
inline constexpr int kBracketingSamples = 64;

// Resout toute la trajectoire, en conservant la dependance sequentielle.
// `platform_origin` a nullptr desactive la detection de singularite (cas
// `platform_origin is None` cote Python).
// `initial_previous_angles` a nullptr = premier point de la trajectoire
// (aucun angle precedent) ; sinon, angles du point precedent, pour reprendre
// apres un repli traite cote Python.
IKBatchResult resolve_axis_angles_batch(const double* desired_directions,
                                        const double* positions,
                                        std::int64_t point_count,
                                        const std::vector<AxisSpec>& axes,
                                        const Vec3& reference_direction,
                                        const Vec3* platform_origin,
                                        double singularity_radius_epsilon,
                                        const double* initial_previous_angles,
                                        const GridFallback* grid = nullptr);

}  // namespace apex
