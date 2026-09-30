// Segments d'intersection triangle-triangle + nombre d'enroulement
// generalise - portage de `pipeline/_self_intersection_resolve.py`
// (triangle_triangle_segment, generalized_winding_number,
// _intersection_segments).
//
// CONTEXTE. Mesure au profil (pipeline complet step1->step8, death.stl,
// resolution 0.8) : `_finish_geometric_resolution` / `resolve_self_intersections`
// pese ~13.5s sur 40.4s totaux (le plus gros poste Python restant apres le
// pas 8). Le fichier source (953 lignes) est majoritairement du calcul
// numpy maison (102 refs np. contre 9 trimesh.) - profil inverse de
// step1_ingest.py, donc un bon candidat.
//
// PERIMETRE. Seules les fonctions CHAUDES et SANS dependance shapely sont
// portees : `triangle_triangle_segment` (15856 appels mesures, 2.59s
// cumulees), `generalized_winding_number` (formule fermee vectorisee,
// 0.92s), et le broad+narrow-phase de `_intersection_segments` (3.22s).
// `_subdivide_face`, `_split_keyhole_polygon`, `_weld_line_endpoints`,
// `_close_small_boundary_loops` NE SONT PAS portees : elles dependent
// fortement de shapely (polygonize, unary_union, Polygon) et networkx,
// avec beaucoup d'etat/cas particuliers - meme diagnostic que la partie
// non portable de step1_ingest.py (pas 8).
//
// EQUIVALENCE : TOLERANCE EXPLICITE, PAS BIT-EXACTE - decision explicite,
// comme step4 (field_solver.hpp), et pour une raison analogue.
// `triangle_triangle_segment` resout un systeme lineaire 3x3
// (`np.linalg.solve`, LAPACK dgesv a pivot partiel) pour localiser un point
// de reference sur la droite d'intersection des deux plans porteurs. Le pas
// 7 a deja mesure qu'une inversion 3x3 par cofacteurs (la seule alternative
// sans nouvelle dependance) diverge de LAPACK a ~1e-14 relatif - jamais au
// bit. Ce module est par ailleurs deja construit sur des tolerances
// EXPLICITES largement superieures a ce bruit (_EPS=1e-10, snap a distance
// 1e-9*scale) : l'oracle compare donc avec une tolerance documentee, pas
// `array_equal`.
//
// `generalized_winding_number`, en revanche, est une formule fermee
// (arctan2 sur des produits scalaires/vectoriels) sans solve, donc
// QUASI-bit-exacte : chaque terme de la somme sur les triangles suit le
// meme ordre d'operations que le numpy source (Van Oosterom & Strackee).
// MESURE (pas supposee) : ecart maximal 1.1e-15 sur 300 points x 1280
// triangles - un ULP de bruit, du a la reduction `np.sum` de numpy
// (potentiellement pairwise/vectorisee) qui n'accumule pas terme a terme
// comme la boucle C++ ci-dessous, meme phenomene que le solveur CG de
// step4. Sans consequence pratique : ce nombre sert a un seuil
// `winding < 0.5` (voir resolve_self_intersections), largement au-dessus
// de ce bruit.
//
// Le broad-phase de `_intersection_segments` reutilise TriangleAabbBvh
// (apex/internal/triangle_aabb_bvh.hpp, deja partage avec mesh_repair.cpp,
// pas 8).
//
// LIMITE ASSUMEE ET MESUREE, DIFFERENTE DU PAS 8 : ici, CONTRAIREMENT au
// SAT (qui ne rend qu'un booleen), l'ORDRE des segments retournes PAR FACE
// est significatif en aval. `_subdivide_face` (non porte, cf. ci-dessus)
// itere `per_face[i]` dans l'ordre recu et ce dernier influence quels
// points de la chaine de decoupe sont vus en premier, donc la
// triangulation exacte produite. Le Python recoit cet ordre de
// `tree.intersection(query)` (rtree/libspatialindex), qui n'est PAS trie
// par indice de triangle (verifie : aucune requete testee sur turkish.stl
// ne rend un ordre trie) - un ordre spatial interne a la structure R-tree,
// non reproductible sans reimplementer libspatialindex. Le tri par indice
// croissant retenu ici NE reproduit PAS cet ordre, il le rapproche
// seulement sans cout.
//
// MESURE bout-en-bout sur `res/test_models/turkish.stl`
// (`resolve_self_intersections`, le seul appelant reel) :
//   - vitesse : 29.98s (Python/rtree) -> 2.60s (natif), soit **x11.5** ;
//   - resultat : 10270 faces (Python) contre 10264 (natif), volume
//     identique a 1.13e-5 RELATIF pres, mais `euler_number` -34 contre -35
//     - une vraie divergence TOPOLOGIQUE mineure (pas seulement du bruit
//     flottant), signature d'une boucle de bord fermee differemment par
//     `_close_small_boundary_loops` en aval.
// Decision explicite de l'utilisateur : le gain (x11.5, sur le plus gros
// poste Python restant apres le pas 8) l'emporte sur cette divergence
// mineure, GARDEE et DOCUMENTEE plutot que masquee - `resolve_self_intersections`
// est deja un mecanisme best-effort qui rend `None` dans de nombreux cas
// (cf. sa docstring), pas un calcul cense produire une reponse canonique
// unique.

#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "apex/types.hpp"

namespace apex {

struct Segment3 {
    Vec3 start;
    Vec3 end;
};

// Segment d'intersection de deux triangles non coplanaires - transcription
// de `triangle_triangle_segment`. `std::nullopt` couvre tous les cas None
// du Python (triangle degenere, pas d'intersection, contact ponctuel,
// plans paralleles, systeme 3x3 singulier).
std::optional<Segment3> triangle_triangle_segment(const Triangle3& a, const Triangle3& b);

// Nombre d'enroulement generalise de `points` (P,3) par rapport a la
// surface fermee `triangles` (T,3,3) - transcription bit-exacte de
// `generalized_winding_number` (Van Oosterom & Strackee).
void generalized_winding_number(const double* points, std::int64_t point_count,
                                const double* triangles, std::int64_t triangle_count,
                                double* out);

// Un segment d'intersection pour une PAIRE de faces (i,j), i<j, du maillage.
struct IntersectionSegment {
    std::int32_t face_a = 0;
    std::int32_t face_b = 0;
    Segment3 segment;
};

// Transcription du broad+narrow-phase de `_intersection_segments` :
// exclut les paires de triangles partageant >=2 sommets (arete commune,
// contact par construction), filtre les segments plus courts que
// `min_length`. `faces` : `triangle_count*3` entiers (indices de sommets,
// meme role que dans mesh_repair.hpp - detection du partage de sommet).
std::vector<IntersectionSegment> find_intersection_segments(const double* triangles,
                                                             const std::int64_t* faces,
                                                             std::int64_t triangle_count,
                                                             double min_length);

}  // namespace apex
