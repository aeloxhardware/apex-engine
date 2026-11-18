// TetraNet — points de la surface d'entrée et classement peau extérieure / face interne, partagés par les contrôles
// (écart surface -> bord) et le raffinement de la peau.
#pragma once
#include "internal.hpp"

namespace tetranet::detail {

// Point de la surface d'entrée, rattaché à un triangle qui le porte (son classement peau extérieure / face interne).
struct Sample { V3 p; int32_t tri; };

// Sommets (les coins vifs comptent) ; dans chaque triangle, ceil(2 aire / pas²) points d'une suite à faible discrépance
// (R2) ; le long des arêtes plus longues que le pas, un point par pas (une seule fois par arête manifold : dans le
// triangle qui la parcourt du petit vers le grand indice). Un triangle long et fin n'est donc pas sur-échantillonné.
std::vector<Sample> surface_samples(const Surface& s, double spacing);

// Un triangle est de la peau EXTÉRIEURE du solide si, juste au-dessus de son centre (le long de sa normale), on est
// dehors (enroulement < 1/2) ; sinon c'est une face interne (coques qui se chevauchent), qui n'a pas à être au bord.
// Calcul paresseux, une fois par triangle, sûr entre fils (course bénigne : même valeur). Sans enroulement : tout est
// peau extérieure.
class OuterSkin {
public:
    OuterSkin(const Surface& s, const Winding* wn, double eps);
    bool outer(int32_t tri) const;
private:
    const Surface& s_;
    const Winding* wn_;
    double eps_;
    std::unique_ptr<std::atomic<int8_t>[]> cls_;   // 0 inconnu, 1 peau extérieure, 2 face interne
};

// Bord d'un maillage tétraédrique comme surface interrogeable (point le plus proche) : faces de bord, sommets partagés.
struct SkinIndex {
    Surface s;
    Bvh bvh;
    SkinIndex(const std::vector<V3>& points, const std::vector<std::array<int32_t, 3>>& faces);
    SkinIndex(const SkinIndex&) = delete;
    SkinIndex& operator=(const SkinIndex&) = delete;
};

// Angle dièdre minimal (degrés) d'un tétraèdre.
double min_dihedral(const V3 p[4]);

}  // namespace tetranet::detail
