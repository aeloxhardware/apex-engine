// Arbre k-d 3D statique - portage de l'usage de `scipy.spatial.cKDTree` dans
// step1 (`mesh_data.kdtree`), step3, step4, step5 et step7.
//
// POURQUOI CE MODULE EST PORTE : contrairement aux pas 1-4, ce n'est PAS un
// point chaud du profil. `scipy.spatial.cKDTree` est deja natif (C), donc le
// gain de vitesse attendu est nul a legerement negatif. Le but est
// ARCHITECTURAL : `cKDTree` est la derniere primitive spatiale de scipy
// utilisee dans le coeur, et la retirer rapproche le moteur de la cible
// 100 % C++ (limite #56). Les gains de vitesse ne sont pas cherches ici.
//
// PERIMETRE VOLONTAIREMENT MINIMAL. Les 9 sites d'appel `.query(...)` du
// pipeline ont tous ete relus : TOUS passent un tableau `(M,3)` et demandent
// le PLUS PROCHE VOISIN unique (`k=1` par defaut de scipy). Aucun ne demande
// k>1, ni une requete par lot de rayons, ni un arbre-contre-arbre. Ce module
// n'expose donc que :
//   - `query_nn`      : 1-NN pour un lot de points (step3/4/5/6/7) ;
//   - `query_radius`  : indices dans un rayon autour d'UN point
//                       (`core/constraints.py::query_ball_point`, seul site).
// Porter davantage serait du code a maintenir et une occasion de diverger de
// scipy pour un cas que personne n'appelle.
//
// EQUIVALENCE.
//   - `query_nn` rend la distance EUCLIDIENNE (p=2, comme le defaut de
//     `cKDTree.query`), calculee `sqrt(dx*dx + dy*dy + dz*dz)` avec
//     accumulation x -> y -> z, sous `/fp:precise` : meme somme, meme
//     `sqrt` correctement arrondi que scipy, donc meme bit.
//   - EGALITE DE DISTANCE (deux points a distance identique du point
//     requete) : ce module tranche par PLUS PETIT INDICE (comparaison
//     stricte `<` en balayant). Le comportement de scipy dans ce cas n'est
//     pas specifie, mais les donnees reelles du pipeline (centres de
//     tetraedres, centres de faces, sommets en virgule flottante) ne
//     produisent jamais d'egalite exacte. L'oracle d'equivalence utilise
//     des donnees non degenerees pour cette raison.
//   - `query_radius` : rayon INCLUSIF (`distance <= radius`), comme
//     `cKDTree.query_ball_point`. L'ordre des indices rendus n'est pas
//     specifie - le seul appelant re-filtre et re-pondere derriere, il est
//     insensible a l'ordre (verifie).

#pragma once

#include <cstdint>
#include <vector>

#include "apex/types.hpp"

namespace apex {

// Arbre k-d 3D sur un nuage de points fixe (median-split sur l'axe le plus
// etendu). Construit une fois, interroge N fois.
class KDTree {
public:
    // `points` : `point_count * 3` doubles, ligne-majeur. Copies a la
    // construction (l'arbre possede ses donnees, comme `cKDTree`).
    KDTree(const double* points, std::int64_t point_count);

    // 1-NN pour `query_count` points (`query_count * 3` doubles, ligne-majeur).
    // Ecrit l'indice du plus proche point dans `out_indices[i]` et la distance
    // euclidienne correspondante dans `out_distances[i]`.
    void query_nn(const double* queries, std::int64_t query_count,
                  std::int64_t* out_indices, double* out_distances) const;

    // Indices de tous les points a distance euclidienne <= `radius` de
    // `point`. Ajoutes a `out` (jamais efface au prealable - l'appelant
    // fournit un tampon vide ou reutilise).
    void query_radius(const Vec3& point, double radius, std::vector<std::int32_t>& out) const;

    std::int64_t size() const { return point_count_; }
    // Acces aux points copies, pour l'equivalent de `cKDTree.data`.
    const std::vector<double>& data() const { return points_; }

private:
    struct Node {
        std::int32_t left = -1;   // indice de noeud, -1 = absent
        std::int32_t right = -1;
        std::int32_t point = -1;  // indice du point porte par ce noeud
        std::int8_t axis = 0;     // 0/1/2 : axe de coupe
    };

    // Construit recursivement sur `indices_[lo, hi)`, rend l'indice du noeud
    // cree (ou -1 si l'intervalle est vide).
    std::int32_t build(std::int64_t lo, std::int64_t hi);

    void nearest(std::int32_t node_index, double qx, double qy, double qz,
                 std::int64_t& best_index, double& best_sqdist) const;
    void within(std::int32_t node_index, double qx, double qy, double qz,
                double radius, double radius_sq, std::vector<std::int32_t>& out) const;

    inline const double* point_at(std::int32_t i) const {
        return points_.data() + static_cast<std::size_t>(i) * 3;
    }

    std::int64_t point_count_ = 0;
    std::vector<double> points_;          // point_count_ * 3
    std::vector<Node> nodes_;
    std::vector<std::int32_t> indices_;   // permutation de points, scratch de build
    std::int32_t root_ = -1;
};

}  // namespace apex
