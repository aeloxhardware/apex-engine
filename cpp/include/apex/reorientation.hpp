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

}  // namespace apex
