// Pas 21c - etape 7 (pipeline/step7_scheduling.py) : couleur par arete de
// trajectoire, lissage le long du chemin, decoupage "une polyligne = une
// couleur", deadline d'accessibilite des trous temporaires, arbitrage trou /
// purge et ordonnancement earliest-deadline-first.
//
// Non porte : `build_color_dependency_graph`, garde-fou dont le Python ne lit
// pas le resultat (aucune sortie n'en depend).
//
// EQUIVALENCE : identique hors deux sources de tolerance, toutes deux aux
// frontieres :
//   - distance a la surface et triangle le plus proche (MeshDistance contre
//     trimesh.proximity, quelques ULP ; un point a egale distance de deux
//     faces de couleurs differentes peut changer de couleur) ;
//   - sommes de longueurs d'arc (numpy somme par paires au-dela de 8 termes),
//     qui ne comptent que comparees a `merge_tolerance`.
// Ordres transcrits : groupes par couleur dans l'ordre de premiere apparition
// (dict Python), couleur dominante = premiere des plus frequentes, `max` sur
// les groupes = premier du plus grand nombre de segments, tri EDF stable.

#pragma once

#include <cstdint>
#include <optional>
#include <vector>

namespace apex {

using Color = std::optional<std::int64_t>;

struct ScheduleLayer {
    std::int64_t id = 0;
    double frame[9] = {0};                         // origine, axe u, axe v
    std::vector<std::vector<double>> toolpath;     // polylignes, `columns` colonnes par point
    std::int64_t columns = 2;                      // 2 (u, v) ou 3 (u, v, offset)
    std::vector<std::int64_t> region_ids;
    std::optional<double> isovalue;
};

struct ScheduleParams {
    double surface_color_distance = 0.5;
    double reheat_penalty = 1.0, adhesion_risk_weight = 1.0;
    double volume_per_purge = 1.0, time_per_purge = 1.0;
    std::int64_t max_open_holes = 4;
    std::int64_t deadline_ray_samples = 24;
};

struct ScheduleInput {
    const double* vertices = nullptr; std::int64_t vertex_count = 0;
    const std::int64_t* faces = nullptr; const double* face_normals = nullptr; std::int64_t face_count = 0;
    // couleur peinte par face ; vide = aucune peinture (cas par defaut)
    std::vector<Color> face_color;
    double merge_tolerance = 0.0;
    const double* tet_vertices = nullptr; std::int64_t tet_vertex_count = 0;
    const double* field_values = nullptr;
    std::vector<double> escape_directions;         // D*3
    ScheduleParams params;
};

struct SchedulePass {
    std::vector<std::pair<std::int64_t, std::int64_t>> segments;   // (layer_id, indice)
    Color color;
    bool is_hole_fill = false;
    std::int64_t layer_min = 0, layer_max = 0;
    std::vector<std::int64_t> deferred_regions;
};

struct ScheduleResult {
    std::vector<std::vector<std::vector<double>>> toolpath;   // par couche, recoupees par couleur
    std::vector<std::vector<Color>> segment_colors;           // par couche
    std::vector<Color> dominant_color;                        // par couche
    std::vector<SchedulePass> passes;
};

ScheduleResult schedule_layers(const std::vector<ScheduleLayer>& layers, const ScheduleInput& input);

}  // namespace apex
