// Types geometriques de base du moteur C++.
//
// Volontairement reduit au strict necessaire du portage en cours (pas 1 du
// plan : _PrintedGeometryIndex). Les 40 dataclasses de core/mesh_types.py ne
// sont PAS transcrites d'avance : on ne construit que ce qui sert, sinon on
// fige un contrat de donnees avant d'avoir porte le code qui le consomme.
//
// Aucune dependance tierce a ce stade - pas meme Eigen. Le pas 1 ne manipule
// que des AABB et des distances point-segment, pour lesquelles trois doubles
// suffisent.

#pragma once

#include <algorithm>

namespace apex {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct AABB {
    Vec3 lo;
    Vec3 hi;
};

// Un triangle, 3 sommets consecutifs - partage entre mesh_repair.hpp
// (pas 8, SAT) et mesh_intersect.hpp (pas 9, segments d'intersection).
struct Triangle3 {
    Vec3 a, b, c;
};

// AABB d'une capsule (segment dilate de `radius`) - transcription exacte de
// `pipeline/step9_collision_sim.py::_segment_aabb`.
inline AABB segment_aabb(const Vec3& start, const Vec3& end, double radius) {
    AABB box;
    box.lo.x = std::min(start.x, end.x) - radius;
    box.lo.y = std::min(start.y, end.y) - radius;
    box.lo.z = std::min(start.z, end.z) - radius;
    box.hi.x = std::max(start.x, end.x) + radius;
    box.hi.y = std::max(start.y, end.y) + radius;
    box.hi.z = std::max(start.z, end.z) + radius;
    return box;
}

}  // namespace apex
