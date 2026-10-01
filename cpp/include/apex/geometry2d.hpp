// Pas 12 - primitives 2D de core/geometry_utils.py.
//
// `point_in_polygon_2d` : lancer de rayon standard, transcrit terme a terme
// (meme test `(y1 > y) != (y2 > y)`, meme formule d'intersection
// `(x2 - x1) * (y - y1) / (y2 - y1) + x1`, meme comparaison stricte), donc
// AU BIT PRES : le Python opere deja en doubles IEEE sur les memes valeurs.
// Appele par step5 (connexite de boucles) et step6 (classification
// matiere/trou, limite #11) : 0,52 s sur death @ 0,8 en Python.

#pragma once

#include <cstdint>

namespace apex {

// `polygon` : n*2 doubles (x, y), non ferme ou ferme (le Python ne distingue pas).
bool point_in_polygon_2d(double x, double y, const double* polygon, std::int64_t n);

}  // namespace apex
