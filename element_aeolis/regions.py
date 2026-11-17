"""Région d'une couche : masque booléen sur la grille -> polygone shapely (contourpy)."""
from __future__ import annotations

import contourpy
import numpy as np
import shapely
from shapely.geometry import Polygon


def mask_to_geometry(xs: np.ndarray, ys: np.ndarray, mask: np.ndarray):
    """Union des zones vraies du masque (contour à mi-chemin entre nœuds). Vide si aucune."""
    if not mask.any():
        return Polygon()
    return level_geometry(xs, ys, mask.astype(float), 0.5, outside=0.0)


def level_geometry(xs: np.ndarray, ys: np.ndarray, values: np.ndarray, level: float, outside: float | None = None):
    """Union des zones où `values` (grille, interpolation bilinéaire) est ≥ `level`. Vide si aucune. Hors de la
    grille et aux nœuds non finis, la valeur est `outside` (défaut : level − 1)."""
    finite = np.isfinite(values)
    if not (finite & (values >= level)).any():
        return Polygon()
    lowest = level - 1.0 if outside is None else outside
    v = np.where(finite, values, lowest)
    top = max(float(v.max()), level) + 1.0
    padded = np.pad(v, 1, constant_values=lowest)
    dx, dy = xs[1] - xs[0], ys[1] - ys[0]
    gx = np.r_[xs[0] - dx, xs, xs[-1] + dx]
    gy = np.r_[ys[0] - dy, ys, ys[-1] + dy]
    gen = contourpy.contour_generator(gx, gy, padded, fill_type=contourpy.FillType.OuterOffset)
    pts_list, off_list = gen.filled(level, top)
    polys = []
    for pts, off in zip(pts_list, off_list):
        rings = [pts[a:b] for a, b in zip(off[:-1], off[1:])]
        if len(rings[0]) >= 4:
            polys.append(Polygon(rings[0], [r for r in rings[1:] if len(r) >= 4]))
    if not polys:
        return Polygon()
    return shapely.union_all([p.buffer(0) for p in polys])
