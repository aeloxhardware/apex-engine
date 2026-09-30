// BVH de triangles + intersection rayon/maillage - portage de l'usage de
// `trimesh.ray.ray_triangle.RayMeshIntersector` dans le pipeline.
//
// POURQUOI CE MODULE, ET CE QU'IL N'EST PAS. Comme le pas 5 (k-d tree), c'est
// un portage ARCHITECTURAL : il retire `trimesh.ray` + `rtree` du coeur et
// debloque le portage de step1/step2/step7 (et, plus tard, de
// `ProximityQuery`, qui partagera cette meme structure). Contrairement au
// k-d tree, la reference Python n'est PAS native : `pyembree` n'est pas
// installe dans ce venv (verifie, pas suppose), donc `RayMeshIntersector`
// tombe sur l'implementation pure numpy de trimesh. Un gain de vitesse est
// donc possible ici - mesure : 27 400 rayons/s sur death.stl (2418 faces),
// mais 2 000 rayons/s sur death_hi-res.stl (226 156 faces). Le cout suit le
// nombre de faces, ce que le broad phase de trimesh ne borne pas.
//
// PERIMETRE VOLONTAIREMENT MINIMAL, verrouille par lecture des sites d'appel.
// UN SEUL point d'entree est utilise par tout le pipeline :
// `intersects_location`, appele par `step2_overhang` (detection de ponts,
// rayons verticaux) et `step7_scheduling` (deadline d'accessibilite, un rayon
// a la fois). `intersects_first`, `intersects_any` et `contains_points` ne
// sont appeles NULLE PART - non portes, comme la grille de repli du pas 2.
//
// LES DEUX CONSOMMATEURS IGNORENT `index_tri` (`_` et `_index_tri` dans leur
// code) et sont INSENSIBLES A L'ORDRE : tous deux reduisent par `.min()` sur
// les distances. Le contrat reellement porte est donc l'ENSEMBLE des couples
// (index_ray, location) dedoublonnes - meme discipline de comparaison que le
// pas 1 pour `_PrintedGeometryIndex`, dont l'ordre venait deja d'un `set`.
//
// EQUIVALENCE : LA NARROW PHASE EST UNE TRANSCRIPTION, PAS UN ALGORITHME
// MEILLEUR. Il aurait ete tentant d'ecrire un Moller-Trumbore, plus court et
// plus rapide ; il donne d'autres flottants, donc d'autres resultats aux
// bords de tolerance, donc un oracle inutilisable. Ce module reproduit la
// chaine exacte de trimesh, operation par operation :
//   1. `intersections.planes_lines` - intersection LIGNE INFINIE / plan du
//      triangle, valide si |dot(direction, normale)| > 1e-5 (seuil ABSOLU,
//      pas relatif : le reproduire tel quel) ;
//   2. `triangles.points_to_barycentric` methode "cramer" (le defaut), avec
//      le meme ordre d'operations et la meme forme `1 - b1 - b2` ;
//   3. test d'appartenance `-tol.zero < b < 1 + tol.zero` avec
//      tol.zero = 1e-13 ;
//   4. filtre "vers l'avant" `dot(location - origine, direction) > -1e-6`
//      (produit scalaire avec la direction NON normalisee, comme
//      `util.diagonal_dot`) ;
//   5. dedoublonnage `grouping.unique_rows` sur (location, index_ray),
//      quantifie a 8 decimales (`tol.merge = 1e-8`). C'est lui qui fait
//      qu'un impact PILE SUR UNE ARETE, touchant deux triangles au meme
//      point, ne compte qu'une fois.
//
// LE BROAD PHASE FAIT PARTIE DU COMPORTEMENT OBSERVABLE - LECON PAYEE ICI.
// Une premiere version descendait le BVH par un test rayon/AABB classique
// (tranches), en supposant que le broad phase n'etait qu'une optimisation :
// deux sur-ensembles differents, la narrow phase tranchant ensuite de la
// meme facon. C'ETAIT FAUX, et la mesure sur le pipeline l'a montre la ou
// 3 000 rayons aleatoires ne voyaient rien - `step2` passait de 49 ponts
// detectes a 0.
//
// Cause : `ray_bounds` clippe le parametre `t` a `buffer_dist = 1e-5` avant
// de construire la boite du rayon, ce qui EXCLUT DES CANDIDATS la face d'ou
// part le rayon. Le filtre "vers l'avant" de la narrow phase, lui, la
// laisserait passer : `step2` decale son origine de 1e-6 sous la face, donc
// l'impact avec sa propre face est a `dot = -1e-6`, et le calcul flottant
// reel donne -9.999999994736442e-07, strictement superieur au seuil -1e-6.
// Autrement dit : **c'est le broad phase de trimesh, plus strict que sa
// propre narrow phase, qui evite l'auto-intersection.** Le reproduire n'est
// donc pas facultatif.
//
// D'ou le choix retenu : le BVH remplace **rtree comme INDEX**, pas la
// REQUETE. La boite du rayon est celle de `ray_bounds`, transcrite a
// l'identique, et le BVH sert a trouver vite les triangles dont l'AABB la
// recoupe - exactement ce que faisait `tree.intersection(bounds)`. Les AABB
// de triangles sont donc EXACTES et sans marge, comme celles indexees par
// `triangles.bounds_tree`. Le jeu de candidats est identique par
// construction, et non par chance.
//
// Verification prealable conservee, car elle reste vraie et utile : compare
// a un Moller-Trumbore exhaustif sur les 2418 triangles de death.stl
// (300 rayons), le broad phase de trimesh ne RATE aucun impact reel ; les
// ecarts de comptage observes s'expliquaient entierement par le
// dedoublonnage d'arete du point 5 (7 triangles touches pour 6 points
// distincts).

#pragma once

#include <cstdint>
#include <vector>

#include "apex/types.hpp"

namespace apex {

// Constantes reprises telles quelles de trimesh - toute divergence casserait
// l'equivalence numerique. Voir l'en-tete pour l'origine de chacune.
inline constexpr double kPlaneParallelEps = 1e-5;   // intersections.planes_lines
inline constexpr double kBarycentricTol = 1e-13;    // trimesh.constants.tol.zero
inline constexpr double kForwardTol = -1e-6;        // ray_triangle_id, filtre "forward"
inline constexpr double kMergeDigits = 1e8;         // tol.merge = 1e-8 -> 8 decimales
inline constexpr double kRayBufferDist = 1e-5;      // ray_bounds(buffer_dist=...)

// Un impact, tel que le pipeline le consomme.
struct RayHit {
    std::int64_t ray_index = 0;
    std::int32_t triangle_index = 0;
    Vec3 location;
};

// BVH statique sur un ensemble de triangles fixe (median-split sur l'axe le
// plus etendu, meme strategie que `kdtree.cpp` - coherent avec le reste du
// noyau, et suffisant : les maillages du projet sont statiques et construits
// une seule fois par execution).
class TriangleBVH {
public:
    // `triangles` : `triangle_count * 9` doubles (v0, v1, v2 par triangle,
    // ligne-majeur). `normals` : `triangle_count * 3`, les normales de face
    // du maillage - PASSEES et non recalculees, car trimesh utilise
    // `mesh.face_normals` quand elles sont disponibles (c'est le cas ici) et
    // une normale recalculee donnerait d'autres flottants.
    TriangleBVH(const double* triangles, const double* normals, std::int64_t triangle_count);

    // Tous les impacts de `ray_count` rayons, DEDOUBLONNES sur
    // (location quantifiee, index_ray) comme `intersects_location`.
    // Les rayons sont parcourus dans l'ordre ; l'ordre a l'interieur d'un
    // rayon suit la descente du BVH et n'est pas significatif (voir en-tete).
    std::vector<RayHit> intersects_location(const double* origins, const double* directions,
                                            std::int64_t ray_count) const;

    std::int64_t size() const { return triangle_count_; }

private:
    struct Node {
        AABB box;
        std::int32_t left = -1;    // -1 sur une feuille
        std::int32_t right = -1;
        std::int32_t first = 0;    // indice dans `order_`
        std::int32_t count = 0;    // 0 sur un noeud interne
    };

    std::int32_t build(std::int64_t lo, std::int64_t hi);
    // Candidats dont l'AABB recoupe `query` - equivalent exact de
    // `tree.intersection(bounds)`. `out` n'est pas efface : l'appelant
    // fournit un tampon reutilise.
    void gather(std::int32_t node_index, const AABB& query,
                std::vector<std::int32_t>& out) const;
    // Transcription de `ray_triangle.ray_bounds` pour UN rayon.
    AABB ray_bounds(const Vec3& origin, const Vec3& direction) const;

    std::int64_t triangle_count_ = 0;
    std::vector<double> triangles_;  // triangle_count_ * 9
    std::vector<double> normals_;    // triangle_count_ * 3
    std::vector<AABB> boxes_;        // AABB EXACTE par triangle (comme bounds_tree)
    std::vector<Vec3> centroids_;
    std::vector<std::int32_t> order_;
    std::vector<Node> nodes_;
    AABB mesh_bounds_;               // union des AABB, = `tree.bounds` cote rtree
    std::int32_t root_ = -1;
};

}  // namespace apex
