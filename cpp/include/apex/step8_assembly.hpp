// Etape 8 - assemblage des points bruts de la trajectoire
// (`_iter_schedule_points`, pipeline/step8_machine_toolpath.py), en natif.
//
// Parcourt les segments du planning DANS L'ORDRE, releve chaque polyligne 2D du
// repere local vers le 3D (`_lift_to_3d`, memes operations flottantes dans le
// meme ordre : au bit), et insere un point de trajet a vide devant chaque
// polyligne dont la cle (couche, segment) differe de la precedente.
#pragma once

#include <cstdint>
#include <vector>

namespace apex {

struct RawPointsOutput {
    std::vector<double> position;            // N*3
    std::vector<double> desired_direction;   // N*3, normale de couche (signe de deposition_reversed)
    std::vector<std::int64_t> layer_id;      // N (identifiant de couche, pas l'indice de table)
    std::vector<std::int64_t> segment_index; // N
    std::vector<std::int64_t> pass_index;    // N : planning -> couleur de la passe
    std::vector<std::uint8_t> is_travel;     // N
};

// Table de couches (L entrees) : `frames` L*9 (origine, u, v), `reversed` L,
// `layer_ids` L, `poly_begin` L+1 (CSR sur les polylignes), `point_begin` P+1 (CSR
// sur `points`), `points` M*columns (columns = 2 ou 3).
// Segments du planning (S) : indice dans la table des couches, indice de polyligne
// DANS sa couche, indice de passe.
RawPointsOutput build_raw_points(std::int64_t layer_count, const double* frames,
                                 const std::uint8_t* reversed, const std::int64_t* layer_ids,
                                 const std::int64_t* poly_begin, const std::int64_t* point_begin,
                                 const double* points, std::int64_t columns, std::int64_t segment_count,
                                 const std::int64_t* seg_layer, const std::int64_t* seg_index,
                                 const std::int64_t* seg_pass);

}  // namespace apex
