// Pas 19 - controle de collision rapide de step9 (`FastCollisionChecker`).
//
// La matiere deposee est densifiee le long des cordons ; chaque point porte
// l'indice du point de trajectoire qui le depose (`owner`). Pour un point de
// trajectoire d'indice `limit`, la question est : distance MINIMALE entre les
// sommets de tete positionnes et la matiere DEJA deposee (owner < limit), dans
// le rayon `reach`.
//
// La version Python (K plus proches voisins scipy, puis requete par boule
// quand la liste est saturee de matiere future) calcule exactement ce
// minimum. Ici, un arbre k-d dont chaque noeud connait la BOITE et le PLUS
// PETIT `owner` de son sous-arbre : un sous-arbre entierement futur ou plus
// loin que le meilleur candidat est elague d'un coup. Meme minimum, distance
// en meme forme close (sqrt(dx*dx + dy*dy + dz*dz), ordre x, y, z) : resultat
// identique au bit (tests/test_cpp_equivalence.py).

#pragma once

#include <cstdint>
#include <vector>

namespace apex {

class PastMaterialIndex {
public:
    PastMaterialIndex(const double* points, const std::int64_t* owners, std::int64_t count);

    // Minimum, sur les `query_count` requetes, de la distance a un point
    // d'owner < `limit` situe a distance <= `reach`. +inf si aucun.
    double min_past_distance(const double* queries, std::int64_t query_count, std::int64_t limit,
                             double reach) const;

    // `pose_count` poses de `samples` points chacune (pose_count*samples*3
    // doubles) : `out[p]` = min_past_distance de la pose p avec la limite
    // `limits[p]`. Meme calcul, pose par pose, reparti sur `threads` fils
    // (<= 0 : nombre de coeurs) -- resultat identique au bit a l'appel unitaire.
    void min_past_distance_batch(const double* queries, std::int64_t pose_count, std::int64_t samples,
                                 const std::int64_t* limits, double reach, double* out, int threads = 0) const;

    // Vrai si un point d'owner < `limit` est a distance <= `radius` d'une des
    // `query_count` requetes (meme test que scipy query_ball_point puis
    // `(owner < limit).any()`). Arret au premier point trouve.
    bool any_past_within(const double* queries, std::int64_t query_count, std::int64_t limit,
                         double radius) const;

    std::int64_t size() const { return count_; }

private:
    bool search_any(std::int32_t node, const double* q, std::int64_t limit, double radius_sq) const;
    struct Node {
        std::int32_t left = -1, right = -1, point = -1;
        std::int8_t axis = 0;
        double lo[3], hi[3];          // boite du sous-arbre
        std::int64_t min_owner = 0;   // plus petit owner du sous-arbre
    };
    std::int32_t build(std::int64_t lo, std::int64_t hi);
    void search(std::int32_t node, const double* q, std::int64_t limit, double& best_sq) const;

    std::int64_t count_ = 0;
    std::vector<double> points_;
    std::vector<std::int64_t> owners_;
    std::vector<std::int32_t> idx_;
    std::vector<Node> nodes_;
    std::int32_t root_ = -1;
};

}  // namespace apex
