// Etape 10 : reorientation par paliers d'une machine a axe couple a la
// position (`step10_adjustment.reorient_coupled_tiers`, limite #59).
//
// Port du coeur de la boucle (candidats tries par ecart, paliers de tolerance
// sur l'ecart de direction de tete, controle plateau + matiere passee,
// propagation ponderee aux voisins de la polyligne). Meme controle que
// `FastCollisionChecker.penetrations` : le minimum de distance a la matiere
// DEJA deposee vient du `PastMaterialIndex` (pas 19), la tete est positionnee
// par rotations de Rodrigues en ordre inverse de la liste d'axes.
//
// EQUIVALENCE : memes operations flottantes que le Python sauf le produit
// scalaire `vecteurs @ axe` (BLAS) ; les angles trouves sont des valeurs de
// la grille, donc identiques tant qu'aucun candidat ne tombe a ~1e-15 d'un
// seuil de palier.

#pragma once

#include <cstdint>
#include <vector>

#include "apex/past_material.hpp"

namespace apex {

struct ReorientationInput {
    std::int64_t n = 0;              // points de trajectoire
    std::int64_t axis_count = 0;     // axes de rotation
    const double* positions = nullptr;   // n*3
    const double* extrusion = nullptr;   // n
    const double* angles = nullptr;      // n*axis_count (degres)
    const double* axes = nullptr;        // axis_count*3 (unitaires)
    std::int64_t free_axis = 0;
    const double* grid = nullptr;        // valeurs candidates de l'axe libre
    std::int64_t grid_count = 0;
    const double* tiers = nullptr;       // paliers de tolerance (degres)
    std::int64_t tier_count = 0;
    const double* samples = nullptr;     // S*3, echantillons de tete (repere local)
    std::int64_t sample_count = 0;
    double reach = 0.0;                  // safety_margin + rayon de buse
    bool has_plate = false;
    double plate_origin[3] = {0, 0, 0};
    double plate_radius = 150.0;
    double plate_half_height = 0.5;
    double plate_center_offset = -0.5;
    double safety_margin = 1.0;
    double lo = 0.0, hi = 185.0;         // plage de l'axe libre
    double full_mm = 1.0, taper_mm = 3.0;
};

struct ReorientationEntry {
    std::int64_t index;
    bool found;
    double angle;
    std::int64_t tier;
};

// `contacts` : indices tries, uniques. `new_angles` (n*axis_count) est ecrit
// (copie de `angles` puis corrections).
std::vector<ReorientationEntry> reorient_coupled_tiers(const ReorientationInput& in,
                                                       const PastMaterialIndex& material,
                                                       const std::int64_t* contacts,
                                                       std::int64_t contact_count,
                                                       double* new_angles);

// `clamp_head_height` (meme fichier Python) : pour chaque point, ramene l'angle
// de l'axe libre a la valeur autorisee la plus proche de la grille si la tete
// passerait sous la marge du plateau. `zlow[g]` : point le plus bas de la tete
// pour grid[g]. Ecrit `out_angle` (n) et `moved` (n, 0/1).
void clamp_head_height_apply(const double* z, const double* current, std::int64_t n,
                             const double* grid, const double* zlow, std::int64_t grid_count,
                             double plate_z, double safety_margin, double* out_angle,
                             std::uint8_t* moved);

// `_try_local_reorientation` (etape 10, niveau 1) pour une machine a premier
// axe couple a la position : C fige, grille produit sur les axes internes,
// candidats par ecart de direction croissant (tri stable), seuil de surplomb,
// limite de vitesse angulaire vis-a-vis des voisins, controle rapide (plateau
// + matiere passee). Rend les angles retenus dans `out` (axis_count) et vrai,
// sinon faux.
struct LocalReorientationInput {
    std::int64_t n = 0, axis_count = 0, point_index = 0;
    const double* positions = nullptr;    // n*3
    const double* feedrates = nullptr;    // n
    const double* angles = nullptr;       // n*axis_count, etat COURANT
    const double* axes = nullptr;         // axis_count*3, unitaires
    const double* inner_lo = nullptr;     // axis_count-1
    const double* inner_hi = nullptr;
    std::int64_t samples_per_axis = 48;
    double max_direction_deviation_deg = 45.0;
    double max_velocity_deg_s = 0.0;      // <= 0 : pas de limite
    std::int64_t neighbor_margin = 1;
    const double* samples = nullptr;      // S*3
    std::int64_t sample_count = 0;
    double reach = 0.0;
    bool has_plate = false;
    double plate_origin[3] = {0, 0, 0};
    double plate_radius = 150.0, plate_half_height = 0.5, plate_center_offset = -0.5;
    double safety_margin = 1.0;
};

bool try_local_reorientation_coupled(const LocalReorientationInput& in,
                                     const PastMaterialIndex& material, double* out);

// Controle rapide de step9 pour un LOT de poses (`FastCollisionChecker.penetrations` x m),
// en parallele : tete positionnee (rotations d'axe interne a externe, translation),
// distance au plateau (disque d'axe Z) puis a la matiere DEJA deposee
// (`limits[k]` = indice du point de trajectoire). `plate_pen`/`part_pen` (m) :
// penetration, ou NaN si aucune.
struct FastBatchInput {
    std::int64_t axis_count = 0;
    const double* axes = nullptr;         // axis_count*3 unitaires
    const double* samples = nullptr;      // S*3
    std::int64_t sample_count = 0;
    double reach = 0.0;
    bool has_plate = false;
    double plate_origin[3] = {0, 0, 0};
    double plate_radius = 150.0, plate_half_height = 0.5, plate_center_offset = -0.5;
    double safety_margin = 1.0;
};

void fast_penetrations_batch(const FastBatchInput& in, const PastMaterialIndex& material,
                             const double* positions, const double* angles, const std::int64_t* limits,
                             std::int64_t count, double* plate_pen, double* part_pen);

}  // namespace apex
