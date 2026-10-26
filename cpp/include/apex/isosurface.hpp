// Pas 13 - marching tetrahedra + soudure des sommets (step5).
//
// Portage de `pipeline/step5_layers.py` : `_extract_isosurface_triangles`
// (qui appelle `_march_tetrahedra_batch`) puis `_weld_vertices`. Le trace de
// bord (`_extract_closed_contours`, trimesh.outline) reste en Python a ce pas.
//
// EQUIVALENCE VISEE : AU BIT PRES, et c'est atteignable ici, contrairement aux
// pas 7/9/10 : aucune reduction, aucun solve. Chaque sommet est
// `p_from + t * (p_to - p_from)` avec `t = (iso - v_from) / (v_to - v_from)`,
// composante par composante, exactement la formule numpy ; l'orientation ne
// depend que du SIGNE d'un produit mixte.
//
// L'ORDRE DE SORTIE est celui de la version numpy par lots, et il compte :
// l'indice des triangles se propage (tetraedre source par triangle) et la
// soudure numerote les sommets. Donc :
//   1. tous les tetraedres "cap" (1 ou 3 sommets au-dessus), dans l'ordre des
//      candidats, 1 triangle chacun ;
//   2. puis tous les "quad" (2 au-dessus), dans l'ordre des candidats,
//      2 triangles consecutifs chacun (tri1 puis tri2).
// Candidats = tetraedres tels que min <= iso <= max (bornes incluses).
// Conventions numpy reprises : `above = v >= iso` ; argmax/argmin rendent le
// PREMIER extremum en cas d'egalite ; `lone` = premier sommet dont `above`
// vaut la valeur minoritaire.
//
// SOUDURE (`_weld_vertices`) : cle = `np.round(v / tol)` en int64 - arrondi
// au pair le plus proche (std::nearbyint, mode FE_TONEAREST par defaut) ;
// `np.unique(axis=0, return_index, return_inverse)` = tri lexicographique
// STABLE des cles (mergesort quand return_index est demande) : les sommets
// soudes sortent dans l'ordre des cles croissantes, chacun represente par sa
// PREMIERE occurrence. Transcrit tel quel.

#pragma once

#include <cstdint>
#include <vector>

namespace apex {

struct IsosurfaceResult {
    std::vector<double> vertices;        // V*3 (soudes si tol > 0)
    std::vector<std::int64_t> faces;     // F*3, indices dans `vertices`
    std::vector<std::int64_t> source_tets;  // F, tetraedre source (indice global)
};

class IsosurfaceExtractor {
public:
    // `tet_vertices` : N*3 ; `tetrahedra` : T*4 ; `values` : N (champ aux noeuds).
    // Les tableaux sont COPIES (l'extracteur vit plus longtemps que l'appel).
    IsosurfaceExtractor(const double* tet_vertices, std::int64_t vertex_count,
                        const std::int64_t* tetrahedra, std::int64_t tet_count,
                        const double* values);

    // Triangles de l'isosurface `isovalue`, soudes a `weld_tolerance` pres
    // (<= 0 : pas de soudure, sommets bruts 3 par triangle).
    IsosurfaceResult extract(double isovalue, double weld_tolerance) const;

    std::int64_t tet_count() const { return static_cast<std::int64_t>(tet_min_.size()); }

private:
    std::vector<double> vertices_;
    std::vector<std::int64_t> tets_;
    std::vector<double> values_;
    std::vector<double> tet_min_;
    std::vector<double> tet_max_;
};

// Pas 14 - trace de bord (`_extract_closed_contours`, trimesh.outline()).
//
// Reproduit EXACTEMENT outline() dans le cas courant : triangles degeneres
// (deux indices egaux) retires, aretes de bord = aretes (triees) presentes
// UNE seule fois, et chaque sommet du bord de degre 2 (boucles simples). Dans
// ce cas, le parcours DFS de scipy (graph.traversals + fill_traversals) part
// du PLUS PETIT indice de chaque composante, suit son PLUS PETIT voisin, fait
// le tour et se referme ; les boucles sortent par indice de depart croissant.
// Verifie empiriquement contre trimesh (tests/test_cpp_equivalence.py).
//
// Sinon (sommet de degre != 2 : bord non-variete, boucle ouverte), rend
// `ok = false` et l'appelant retombe sur trimesh pour CETTE couche - reprise
// a grain fin, jamais tout ou rien (lecon du pas 2).
struct BoundaryLoops {
    bool ok = false;
    std::vector<std::vector<std::int64_t>> loops;  // indices de sommets, fermes (dernier == premier)
};

BoundaryLoops boundary_loops(const std::int64_t* faces, std::int64_t face_count, std::int64_t vertex_count);

// Couche par defaut de step5 (chemin sans clustering) : isosurface soudee, contours
// fermes, region_ids (regions >= 0 des tetraedres sources, tries) et repere local
// (`_compute_local_frame`), en un appel. status : 0 = aucune couche (pas de face ou
// pas de contour), 1 = ok, 2 = bord non simple (l'appelant retombe sur trimesh).
struct LayerExtraction {
    int status = 0;
    std::vector<std::vector<double>> contours;   // chaque contour : K*3
    std::vector<std::int64_t> region_ids;
    double origin[3] = {0, 0, 0}, axis_u[3] = {0, 0, 0}, axis_v[3] = {0, 0, 0};
    std::vector<std::int64_t> unique_tets;
};

LayerExtraction extract_layer(const IsosurfaceExtractor& extractor, double isovalue,
                              double weld_tolerance, const double* gradients,
                              const double* volumes, const std::int64_t* tet_regions);

// Volume de chaque tetraedre (`_tet_volumes`) : |det(aretes)| / 6.
void tet_volumes(const double* vertices, const std::int64_t* tetrahedra, std::int64_t tet_count,
                 double* out);

}  // namespace apex
