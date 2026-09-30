// Tables topologiques d'un maillage tetraedrique - portage de
// `step4_field_solver._tet_face_table` et `step1_ingest._build_tet_adjacency`.
//
// POURQUOI CE MODULE EST LE QUATRIEME PORTE : apres les pas 1-3, le
// re-profilage de death.stl @ 0.8 (66,5 s profilees) designe ces deux
// fonctions comme la plus grosse grappe portable restante, TetGen mis a part
// (15,7 s, deja natif et incompressible) :
//   - `_tet_face_table`      : 1,79 s + 1,78 s de generateur sur 12,5 M
//                              appels, plus une part des 2,37 s de `sorted`
//                              et des 1,43 s de `setdefault` -> ~6-7 s ;
//   - `_build_tet_adjacency` : 2,91 s + 1,35 s de `set.add` -> ~4,3 s.
// Les deux construisent la meme chose - une table de sous-simplexes avec un
// dictionnaire a cles tuple - d'ou un seul module pour les deux.
//
// CONTRAT D'EQUIVALENCE : l'ordre des faces frontiere suit l'ordre
// d'INSERTION (premiere occurrence en balayant les tetraedres dans l'ordre,
// puis les 4 faces locales dans l'ordre (0,1,2), (0,1,3), (0,2,3), (1,2,3)),
// et non l'ordre trie des cles. C'est l'ordre d'iteration du dict Python, dont
// depend l'ordre des lignes de `boundary_face_nodes` - donc celui de la
// requete kdtree en aval dans `_map_boundary_faces_to_surface`. Un ordre
// different produirait des tableaux permutes et des comparaisons faussees.
//
// Pour les ARETES (adjacence noeud-noeud) l'ordre est en revanche sans
// importance : `scipy.sparse.csr_matrix` trie ses indices a la construction,
// donc la matrice finale est identique quel que soit l'ordre fourni. Verifie,
// pas suppose.

#pragma once

#include <cstdint>
#include <vector>

namespace apex {

struct TetFaceTopology {
    // Faces n'appartenant qu'a UN seul tetraedre, dans l'ordre d'insertion.
    std::vector<int32_t> boundary_face_nodes;  // K*3, triplets tries
    std::vector<int32_t> boundary_face_tets;   // K
    // Paires de tetraedres partageant une face interne.
    std::vector<int32_t> adjacency_pairs;      // M*2
};

// `tetrahedra` : N*4 indices de noeuds, ligne-majeur.
TetFaceTopology compute_tet_face_topology(const int64_t* tetrahedra, int64_t tet_count);

// Aretes uniques (a < b) du maillage, 6 par tetraedre avant dedoublonnage.
// Rend un tableau plat M*2.
std::vector<int32_t> compute_unique_edges(const int64_t* tetrahedra, int64_t tet_count);

}  // namespace apex
