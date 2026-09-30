// Detection d'auto-intersections triangle-triangle - portage de
// `pipeline/step1_ingest.py::_triangles_cross` / `_find_self_intersecting_faces`.
//
// PERIMETRE. step1_ingest.py est domine (~90%) par de la logique de
// reparation de maillage qui appelle trimesh/pymeshfix directement
// (slice_plane, booleans, triangulation de polygones, MeshFix) - pas du code
// maison portable, et pas souhaitable a porter (bibliotheques tierces).
// `_build_tet_adjacency` est deja natif (`tet_unique_edges`). Le seul autre
// point chaud reel, identifie par lecture, est le test SAT triangle-triangle
// et la boucle qui l'appelle sur toutes les paires de triangles candidates :
//   - `_find_self_intersecting_faces` (broad-phase `bounds_tree` + narrow-
//     phase `_triangles_cross`), utilisee par le rapport d'ingestion ET pour
//     amorcer la reparation locale (`_repair_local_self_intersections`) ;
//   - `_fill_loop_min_weight` (triangulation Klincsek du rebouchage), qui
//     appelle `_triangles_cross` en boucle interne contre les obstacles et
//     les triangles deja acceptes du patch - fonction elle-meme NON portee
//     (etat Python riche : cache memoire, DP O(n^3)), mais qui beneficie de
//     la primitive SAT unitaire exposee ici.
//
// Chemin nominal (death.stl) : CE MODULE N'EST JAMAIS EXERCE - aucune
// auto-intersection SAT sur ce modele nettoye. Le gain mesurable porte sur
// les modeles qui declenchent la cascade de reparation (turkish.stl : 50
// faces auto-intersectantes detectees par ce test, cf. CLAUDE.md).
//
// EQUIVALENCE. Transcription terme a terme de `_triangles_cross` (theoreme
// de l'axe separateur, 11 axes : 2 normales de face + 9 produits croises
// arete x arete), meme tolerance `_COPLANAR_EPS = 1e-9`, memes axes de norme
// < 1e-12 ignores (paralleles ou face degeneree). Bit-exact par construction
// (memes operations flottantes, meme ordre).
//
// Le broad-phase (`trimesh.triangles.bounds_tree`, un RTree) n'a PAS besoin
// d'etre reproduit a l'identique ici, contrairement au cas du BVH de rayons
// (voir bvh.hpp) : le broad-phase ne fait ici qu'ELARGIR ou RESTREINDRE
// l'ensemble de paires soumises au test EXACT, qui est lui deterministe et
// symetrique - un AABB exact (comme celui indexe par bounds_tree) est un
// sur-ensemble conservateur suffisant, aucune paire vraiment intersectante
// n'est perdue tant que l'AABB reste exact (pas de marge ajoutee).

#pragma once

#include <cstdint>
#include <vector>

#include "apex/types.hpp"

namespace apex {

// Test d'axe separateur (SAT) entre deux triangles - transcription exacte de
// `_triangles_cross`. Ne teste PAS le partage de sommet : c'est a l'appelant
// de filtrer les paires adjacentes AVANT d'appeler cette fonction (meme
// convention que le Python, ou le filtrage est fait par `_find_self_intersecting_faces`
// avant l'appel).
bool triangles_cross(const Triangle3& t0, const Triangle3& t1);

// Indices (dans `faces`) des triangles impliques dans au moins une
// auto-intersection reelle - transcription de `_find_self_intersecting_faces`.
//
// `triangles` : `triangle_count * 9` doubles (v0,v1,v2 par triangle,
// ligne-majeur, meme layout que `mesh.triangles`).
// `faces` : `triangle_count * 3` entiers (indices de sommets), utilises
// UNIQUEMENT pour detecter le partage de sommet entre deux triangles (memes
// discipline que le Python : `np.intersect1d(faces[i], faces[j]).size > 0`
// exclut la paire).
//
// Broad-phase : BVH maison (meme structure que TriangleBVH, AABB exactes),
// construit une fois puis descendu pour chaque triangle - equivalent
// fonctionnel de `bounds_tree(triangles)` + `tree.intersection(query)`, sans
// dependance a `rtree`. Chaque paire (i,j) avec j>i n'est visitee qu'une
// fois, comme le Python.
//
// Retourne les indices en ordre croissant (un `std::vector` trie, pas un
// `set` non ordonne - l'appelant Python le convertit en `set`, insensible a
// l'ordre).
std::vector<std::int32_t> find_self_intersecting_faces(const double* triangles,
                                                        const std::int64_t* faces,
                                                        std::int64_t triangle_count);

}  // namespace apex
