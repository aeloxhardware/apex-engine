"""Surfaces de tranchage : interpolation entre le dessous B et le dessus T.

Couche 1 : plane, à B + t1 (adhérence). Couches k = 2..N : interpolation linéaire de la couche 1
au dessus T, S_k = S_1 + ((k−1)/(N−1))·(T − S_1). L'épaisseur locale d'une couche k ≥ 2 vaut
(T − S_1)/(N − 1) : fine là où la pièce est basse, plus épaisse sous les crêtes.
Constat sur la sortie d'IceSL (plaque 5) : première couche plane de 0,2 à 0,3 mm, puis
interpolation linéaire jusqu'au dessus ; les deux premiers points sont repris ici.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .heightfield import Heightfield


@dataclass
class LayerStack:
    field: Heightfield
    count: int                # N
    first_layer: float = 0.2  # t1

    def z(self, k: int, x, y) -> np.ndarray:
        """Hauteur de dépôt de la couche k (1..N) aux points (x, y)."""
        b = self.field.bottom_at(x, y)
        t = self.field.top_at(x, y)
        s1 = np.minimum(b + self.first_layer, t)
        if self.count == 1 or k <= 1:
            return s1 if self.count > 1 else t
        return s1 + ((k - 1) / (self.count - 1)) * (t - s1)

    def thickness(self, k: int, x, y) -> np.ndarray:
        """Épaisseur locale de la couche k aux points (x, y)."""
        return self.z(k, x, y) - (self.field.bottom_at(x, y) if k <= 1 else self.z(k - 1, x, y))


def build(field: Heightfield, layer_height: float, first_layer: float | None = None) -> LayerStack:
    """N est choisi pour que l'épaisseur maximale des couches 2..N ne dépasse pas `layer_height`."""
    t1 = layer_height if first_layer is None else first_layer
    hmax = float(np.nanmax(field.thickness))
    if hmax <= 0:
        raise ValueError("pièce d'épaisseur nulle")
    rest = max(hmax - t1, 0.0)
    return LayerStack(field, 1 + max(1, math.ceil(rest / layer_height - 1e-9)) if rest > 0 else 1, t1)


def max_slope_deg(field: Heightfield, k: int, count: int, first_layer: float = 0.2) -> float:
    """Pente maximale (degrés) de la surface S_k, mesurée sur la grille."""
    s1 = np.minimum(field.bottom + first_layer, field.top)
    s = s1 + ((k - 1) / max(count - 1, 1)) * (field.top - s1)
    gy, gx = np.gradient(s, field.ys, field.xs)
    g = np.hypot(gx, gy)
    return float(np.degrees(np.arctan(np.nanmax(g))))
