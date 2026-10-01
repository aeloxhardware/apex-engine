// Pas 21b - etape 3, CHEMIN PAR DEFAUT (pipeline/step3_segmentation.py) :
// compatibilite face/direction, croissance de regions, direction par region,
// profondeurs de croissance depuis le plateau, graphe de dependances.
// Les mecanismes opt-in (`branch_cut_threshold`, `anchor_growth_weight`)
// restent en Python.
//
// EQUIVALENCE : identique, ordre des dependances compris, a condition que le
// Python trie les aretes de facon STABLE (adopte le 2026-10-01, voir
// `_region_growing` : le tri non stable de numpy rendait la segmentation
// dependante de son implementation - 11 ou 24 regions sur turkish).
// Ordres transcrits, parce que l'aval en depend :
//   - aretes de faces dans l'ordre (ligne, colonne) de la matrice CSR de
//     `build_face_adjacency`, triees de facon stable par risque decroissant ;
//   - regions renumerotees par identifiant de racine croissant
//     (`np.unique(..., return_inverse=True)`) ;
//   - dependances dans l'ordre ou networkx les rend : aretes non orientees du
//     `nx.Graph` (noeuds croissants, voisins par premiere apparition), puis
//     `DiGraph.edges()` (noeuds croissants, successeurs par ordre d'ajout).
// Seul ecart possible : `acos` (CRT contre numpy) dans la compatibilite, pour
// une face EXACTEMENT au seuil de 45 deg d'une direction candidate.

#pragma once

#include <cstdint>
#include <vector>

namespace apex {

struct SegmentationInput {
    const double* vertices;           // V*3
    const std::int64_t* faces;        // F*3
    const double* face_normals;       // F*3
    const double* risk;               // F (risque EFFECTIF : plateau a 0)
    const std::uint8_t* platform;     // F
    std::int64_t face_count;
    const double* directions;         // D*3 (la verticale en premier)
    std::int64_t direction_count;
    double max_angle_deg;
    double merge_angle_deg;           // < 0 : desactive (None)
    double up[3];
};

struct SegmentationResult {
    std::vector<std::int64_t> labels;           // F, regions 0..R-1
    std::vector<std::int64_t> direction_index;  // R, indice dans `directions`
    std::vector<std::int64_t> depth;            // R
    std::vector<std::int64_t> dependencies;     // E*2, (avant, apres)
};

SegmentationResult segment_default(const SegmentationInput& in);

}  // namespace apex
