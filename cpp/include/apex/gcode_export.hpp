// Pas 24 - corps de l'export G-code (file_io/gcode_export.py::export_gcode).
//
// CONTRAT : TEXTE IDENTIQUE LIGNE A LIGNE a la version Python (test
// tests/test_cpp_equivalence.py). Transcriptions qui ne sont PAS des details :
//   - le modulo Python sur les flottants (`a % b` prend le signe du DIVISEUR,
//     pas celui du dividende comme fmod) : `py_mod` ;
//   - `f"{v:.4f}".rstrip("0").rstrip(".")`, "-0" et "" -> "0" : `fmt`, via
//     std::to_chars (arrondi EXACT, comme le dtoa de Python) ;
//   - norme de (X, Y, Z, B) firmware a QUATRE composantes pour l'avance ;
//   - accumulation sequentielle de E, dans l'ordre des mouvements.
// En-tete, pied de fichier et appels de rappel Python (`travel_check`) restent
// cote Python ; `head_zlow` est passe sous forme de TABLE (grille de B, valeur).
// Precalculs numpy passes tels quels : positions decalees au sommet du cordon,
// facteur de ralentissement `feed_plan`, section du filament.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace apex {

struct GcodeExportParams {
    double layer_height, skin_thickness, line_width, flow;
    double print_speed, travel_speed, rotation_feed, max_feed;
    double retract_length, retract_speed, retract_min_travel, travel_clearance;
    double pivot_to_tip, platform_margin, singularity_radius;
    double tilt_sign;
    bool z_reference_tip;   // true = "tip", false = "pivot"
    int decimals;
    double fil_area;        // pi * (d/2)^2, calcule cote Python
};

struct GcodeExportStats {
    std::int64_t moves = 0, travels = 0, retractions = 0;
    std::int64_t clearance_z = 0, clearance_r = 0, clearance_none = 0;
    double filament_mm = 0.0, c_mismatch_deg = 0.0;
    double tip_z_min = 0.0, head_z_min = 0.0;   // +inf si aucun mouvement / pas de table
    std::int64_t tip_below_plate = 0, head_below_margin = 0;
    double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;
    bool any = false;
};

// `positions` n*3 (deja decalees), `angles` n*2 (C, B), `extrusion` n,
// `skin` n (0/1), `slow` n-1 (facteur feed_plan). `zlow_grid`/`zlow_values`
// : table de head_zlow (vides = pas de table). Rend le CORPS (lignes separees
// par '\n', sans '\n' final) et remplit `stats`.
std::string export_gcode_body(const double* positions, const double* angles, const double* extrusion,
                              const std::uint8_t* skin, const double* slow, std::int64_t n,
                              const std::vector<double>& zlow_grid, const std::vector<double>& zlow_values,
                              const GcodeExportParams& params, GcodeExportStats& stats);

// `feed_plan` (meme fichier Python) : plafonds de vitesse angulaire par passes
// avant puis arriere, remise a zero a chaque changement de signe (np.sign).
// `d` : variations d'angle (degres) de chaque segment ; ecrit `cap`.
void feed_plan_caps(const double* d, std::int64_t count, double max_velocity_deg, double max_accel_deg,
                    double* cap);

}  // namespace apex
