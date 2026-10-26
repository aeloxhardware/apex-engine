// Etape 6 - remplissage 2D d'une couche ENTIERE en natif (`_fill_layer`,
// pipeline/step6_toolpath2d.py) : projection dans le repere local, classification
// des contours par imbrication, detection de pont (kd-tree des centres de faces),
// parois + remplissage (`fill_region`, Clipper2), offset hors-plan (exact pour le
// contour n=0, interpole sinon), fermeture des anneaux, decoupage en troncons
// monotones. Les motifs OSCILLATING_SPIRAL restent en Python.
//
// EQUIVALENCE : celle de `fill_region` (tolerance, pas le bit : Clipper2 contre
// GEOS), plus quelques ULP de BLAS sur les projections `points @ axe`.
// Les ponts passent ici aussi par Clipper2 (le Python utilise shapely pour eux).
#pragma once

#include <cstdint>
#include <vector>

#include "apex/kdtree.hpp"
#include "apex/toolpath2d.hpp"

namespace apex {

struct Step6Config {
    double nozzle_width = 0.4;
    int wall_count = 2;
    double degenerate_area_ratio = 1e-6;
    double infill_density = 0.2;
    InfillKind pattern = InfillKind::CrossedZigzag;
    double bridge_contour_ratio_threshold = 0.05;
    int bridge_anchor_extra_walls = 2;
    bool has_power = false;
    double power = 2.0;
    bool continuous_rings = false;
    bool split_rings_into_runs = false;
    bool split_open_polylines = false;
    double run_split_min_amplitude = 1e-12;
};

struct BridgeInfo {
    const KDTree* tree = nullptr;               // centres de faces ; nullptr = pas de detection de pont
    const std::uint8_t* is_bridge = nullptr;    // par face
    const double* span = nullptr;               // par face, 3 doubles (NaN si non applicable)
};

struct LayerContours {
    const double* frame = nullptr;                   // origine, axe_u, axe_v (9)
    std::vector<const double*> loops;                // chaque boucle : K*3
    std::vector<std::int64_t> loop_sizes;
    double sweep_angle = 0.0;
};

// Polylignes (K*3 : u, v, offset hors-plan) de la couche, dans l'ordre de `_fill_layer`.
std::vector<std::vector<double>> fill_layer(const LayerContours& layer, const Step6Config& config,
                                            const BridgeInfo& bridge);

}  // namespace apex
