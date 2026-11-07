// TetraNet — arêtes vives de la surface d'entrée (piste C de R-030 point 8).
//
// Le réseau régulier ramené à la peau n'a aucune raison de placer un nœud SUR une arête vive : les triangles du
// bord de maillage coupent alors le coin (mesure du 2026-11-07 : jusqu'à 0,55 mm sur death et turkish à cellule
// 1 mm). Ici : repérer les arêtes dont les deux faces font un angle supérieur au seuil, puis ramener un point de
// coupe voisin d'une telle arête SUR l'arête.
//
// Interface volontairement indépendante des types internes de tetranet.cpp (tableaux plats).
#pragma once
#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace tetranet {

class FeatureEdges {
public:
    // points : 3*N ; triangles : 3*M (indices) ; normals : 3*M (normales de face unitaires).
    // Une arête partagée par exactement deux faces dont l'angle entre normales dépasse `angle_deg` est vive.
    // `piece` : longueur maximale d'un morceau d'arête dans l'index (les arêtes longues sont découpées).
    FeatureEdges(const std::vector<std::array<double, 3>>& points, const std::vector<std::array<int, 3>>& triangles,
                 const std::vector<std::array<double, 3>>& normals, double angle_deg, double piece);

    int64_t edge_count() const { return edges_; }
    bool empty() const { return segments_.empty(); }

    // Point de l'arête vive le plus proche de `p`, s'il est à moins de `radius` ; sinon faux.
    bool snap(const std::array<double, 3>& p, double radius, std::array<double, 3>& out) const;

private:
    struct Segment { std::array<double, 3> a, b; };
    int64_t cell_key(int64_t i, int64_t j, int64_t k) const;
    std::vector<Segment> segments_;
    std::unordered_map<int64_t, std::vector<int32_t>> grid_;
    double cell_ = 1.0;
    int64_t edges_ = 0;
};

}  // namespace tetranet
