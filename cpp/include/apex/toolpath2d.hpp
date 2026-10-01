// Pas 15 - remplissage 2D d'une zone de couche (step6 : `_offset_walls` +
// `_generate_infill`), parois par Clipper2 au lieu de shapely/GEOS.
//
// CONTRAT : PAS au bit pres, et c'est assume. Clipper2 et GEOS calculent le
// meme offset en onglet (meme limite d'onglet, 5 fois la distance) mais avec
// une autre arithmetique (Clipper2 : entiers, 1e-6 mm ici) et un autre point
// de depart d'anneau. L'equivalence est verifiee sur des grandeurs physiques
// (aire de chaque niveau de paroi, longueur de remplissage) a tolerance
// declaree, puis sur la metrique step9b du pipeline (CLAUDE.md, limite #56).
// Opt-in cote Python (`step6_toolpath2d._NATIVE_FILL`).
//
// Reproduit, a ces tolerances pres :
//   - orientation des anneaux de GEOS : exterieur HORAIRE, trous ANTI-HORAIRES,
//     anneaux fermes (dernier point = premier) ;
//   - filtrage des morceaux d'offset par `degenerate_area_ratio` ;
//   - motifs LINES, GRID (deux balayages a densite / 2) et CROSSED_ZIGZAG
//     (zone retrecie d'une demi-buse, raccords qui suivent le bord, a defaut
//     raccord droit s'il reste dans la zone) ;
//   - espacement nozzle_width / min(density, 1), borne a 5000 lignes.
// Les contours d'origine (paroi n = 0) restent ajoutes cote Python.

#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace apex {

using Point2 = std::array<double, 2>;
using Ring2 = std::vector<Point2>;

enum class InfillKind : int { Lines = 0, Grid = 1, CrossedZigzag = 2 };

struct FillRegionParams {
    double nozzle_width = 0.4;
    int wall_count = 2;
    double degenerate_area_ratio = 0.0;
    double infill_density = 0.2;
    InfillKind pattern = InfillKind::CrossedZigzag;
    double sweep_angle_rad = 0.0;
    double miter_limit = 5.0;     // mitre_limit par defaut de shapely/GEOS
    double scale = 1e6;           // coordonnees entieres de Clipper2 : 1e-6 mm
};

struct FillRegionResult {
    bool base_valid = false;      // faux : aire nulle, la zone est ignoree (aucune paroi, meme n = 0)
    std::vector<Ring2> infill;    // polylignes ouvertes
    std::vector<Ring2> walls;     // anneaux fermes, du centre vers l'exterieur (n = wall_count .. 1)
};

FillRegionResult fill_region(const Ring2& outer, const std::vector<Ring2>& holes, const FillRegionParams& params);

}  // namespace apex
