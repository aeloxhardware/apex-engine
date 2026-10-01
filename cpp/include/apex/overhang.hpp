// Pas 21a - etape 2 entiere (pipeline/step2_overhang.py::run) : surplombs,
// score de risque, ponts internes par lancer de rayons verticaux, direction
// de portee des ponts (axe principal par composante connexe).
//
// Gain de temps nul (0,02 s en Python) : c'est une piece de l'orchestrateur
// natif, pas une optimisation.
//
// EQUIVALENCE, par sortie :
//   - is_overhang, is_bridge : identiques (les rayons passent par le BVH du
//     pas 6, au bit pres ; les angles ne servent qu'a des comparaisons a un
//     seuil, ou un ULP d'ecart sur acos ne peut basculer qu'une face pile au
//     seuil) ;
//   - risk_score : tolerance 1e-12 (acos et exp du CRT MSVC contre numpy) ;
//   - bridge_span_direction : au SIGNE pres (vecteur propre d'une matrice de
//     covariance 3x3 contre SVD de LAPACK) ; le seul consommateur (step6)
//     balaie perpendiculairement a cette direction, insensible au signe.
//
// Composantes connexes : comme le Python, sur le graphe d'adjacence de TOUTES
// les faces (aretes partagees par exactement deux faces), puis restreintes aux
// faces-ponts - deux ponts relies par des faces ordinaires forment donc UNE
// composante. Transcrit tel quel, ce n'est pas corrige ici.

#pragma once

#include <cstdint>
#include <vector>

namespace apex {

struct OverhangParams {
    double overhang_angle_threshold_deg = 45.0;
    int bridge_detection_ray_density = 4;
    double bridge_support_ratio_threshold = 0.7;
    double bridge_max_gap_distance = 2.0;
    double bridge_horizontal_angle_tolerance_deg = 10.0;
    double up[3] = {0.0, 0.0, 1.0};
};

struct OverhangResult {
    std::vector<double> risk_score;          // F
    std::vector<std::uint8_t> is_overhang;   // F
    std::vector<std::uint8_t> is_bridge;     // F
    std::vector<double> span_direction;      // F*3, NaN hors ponts
};

OverhangResult compute_overhang_map(const double* vertices, std::int64_t vertex_count,
                                    const std::int64_t* faces, const double* face_normals,
                                    std::int64_t face_count, const OverhangParams& params);

}  // namespace apex
