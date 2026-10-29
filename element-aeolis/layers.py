"""Surfaces de tranchage : interpolation entre le dessous B et le dessus T, épaisseur bornée.

Plancher `min_layer_height` et plafond `max_layer_height` viennent du profil d'imprimante
(défauts 0,1 et 0,25 mm).

Couche 1 : plane, à B + t1 (adhérence). Couches k ≥ 2 : S_k = S_1 + (k−1)·d, où d est
l'épaisseur locale d = max((T − S_1)/(N − 1), min_layer_height), avec N tel que d ≤ plafond
partout. Là où la pièce est basse, d est relevée au plancher : la colonne atteint T avant la
couche N et ne reçoit plus de matière ensuite (couche « inactive » en ce point). Le reste
r < plancher d'une colonne est absorbé par sa dernière couche (étirée jusqu'à T), de sorte
que toutes les couches actives respectent plancher et plafond et que la dernière est sur T.

Inter-couches (`intermediate_layers`, activé par défaut) : N est dimensionné sur la colonne
typique (hauteur médiane) et l'épaisseur nominale. Là où l'espacement d entre deux couches
principales dépasse le plafond, chaque intervalle est subdivisé en m = ceil(d / plafond) :
les m − 1 couches intermédiaires ne couvrent QUE les colonnes concernées. Désactivé, N est
dimensionné sur la colonne la plus haute et le plafond borne d partout.

Constat sur la sortie d'IceSL (plaque 5) : première couche plane de 0,2 à 0,3 mm, épaisseurs
fines en bas (≈ 0,07 mm, sous le plancher ici), mais la fonction exacte d'IceSL n'est pas
linéaire en k et n'est pas reproduite.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .heightfield import Heightfield

DEFAULT_MIN_LAYER = 0.1
DEFAULT_MAX_LAYER = 0.25


@dataclass
class LayerStack:
    field: Heightfield
    count: int                       # nombre de couches (maximum sur la pièce)
    first_layer: float = 0.2         # t1
    min_layer: float = DEFAULT_MIN_LAYER
    max_layer: float = DEFAULT_MAX_LAYER
    intermediate: bool = False       # subdivise localement les intervalles trop grands

    # -- géométrie d'une colonne -------------------------------------------------------
    def _column(self, b, t):
        """(s1, d, last) : première surface, épaisseur locale, dernière couche active."""
        b = np.asarray(b, float)
        t = np.asarray(t, float)
        s1 = np.minimum(b + self.first_layer, t)
        rem = t - s1
        n = max(self.count - 1, 1)
        d = np.maximum(rem / n, self.min_layer)
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.where(rem > 1e-9, rem / d, 0.0)
        nfull = np.floor(m + 1e-9)
        r = rem - nfull * d
        stretch = (r > 1e-9) & (r < self.min_layer) & (d + r <= self.max_layer + 1e-9)
        last = 1 + nfull + np.where((r > 1e-9) & ~stretch, 1, 0)
        last = np.where(rem > 1e-9, last, 1)
        return s1, d, last

    def _split(self, d) -> np.ndarray:
        """Nombre de couches par intervalle principal (1 sans inter-couches)."""
        if not self.intermediate:
            return np.ones_like(d)
        return np.maximum(1.0, np.ceil(d / self.max_layer - 1e-9))

    def split_at(self, x, y) -> np.ndarray:
        s1, d, last = self._column(self.field.bottom_at(x, y), self.field.top_at(x, y))
        return self._split(d)

    def last_layer(self, x, y) -> np.ndarray:
        """Indice de la dernière couche active aux points (x, y)."""
        s1, d, last = self._column(self.field.bottom_at(x, y), self.field.top_at(x, y))
        return last

    def z(self, k: int, x, y) -> np.ndarray:
        """Hauteur de dépôt de la couche k aux points ; T au-delà de la dernière couche active."""
        b = self.field.bottom_at(x, y)
        t = self.field.top_at(x, y)
        s1, d, last = self._column(b, t)
        z = np.where(k >= last, t, s1 + (k - 1) * d)
        return np.where(last <= 1, t, z) if k == 1 else z

    def z_sub(self, k: int, j: int, x, y) -> np.ndarray:
        """Couche intermédiaire j (1..m−1) entre les couches principales k et k+1."""
        m = self.split_at(x, y)
        zk, zn = self.z(k, x, y), self.z(k + 1, x, y)
        return zk + (j / m) * (zn - zk)

    def thickness(self, k: int, x, y) -> np.ndarray:
        """Épaisseur locale de la couche principale k (0 là où elle n'est plus active)."""
        if k <= 1:
            return self.z(1, x, y) - self.field.bottom_at(x, y)
        zk, zp = self.z(k, x, y), self.z(k - 1, x, y)
        m = self.split_at(x, y)
        return (zk - zp) / m

    def thickness_sub(self, k: int, j: int, x, y) -> np.ndarray:
        m = self.split_at(x, y)
        return (self.z(k + 1, x, y) - self.z(k, x, y)) / m

    def sequence(self) -> list[tuple[int, int]]:
        """Ordre de dépôt : (k, 0) principale, puis (k, j) intermédiaires de l'intervalle k→k+1."""
        out = []
        mmax = int(np.nanmax(self.split_grid())) if self.intermediate else 1
        for k in range(1, self.count + 1):
            out.append((k, 0))
            if k < self.count:
                out.extend((k, j) for j in range(1, mmax))
        return out

    def split_grid(self) -> np.ndarray:
        f = self.field
        return self._split(self._column(f.bottom, f.top)[1])

    # -- région active d'une couche ------------------------------------------------------
    def last_grid(self) -> np.ndarray:
        """Dernière couche active aux nœuds de la grille (NaN hors empreinte)."""
        f = self.field
        return self._column(f.bottom, f.top)[2]

    def active_mask(self, k: int, j: int = 0, last: np.ndarray | None = None,
                    split: np.ndarray | None = None) -> np.ndarray:
        """Nœuds où la couche (k, j) existe ; vrai hors empreinte (le contour de la pièce coupe)."""
        last = self.last_grid() if last is None else last
        if j == 0:
            return np.where(np.isnan(last), True, last >= k)
        split = self.split_grid() if split is None else split
        return np.where(np.isnan(last), True, (last >= k + 1) & (split > j))


def build(field: Heightfield, layer_height: float | None = None, first_layer: float | None = None,
          *, min_layer: float = DEFAULT_MIN_LAYER, max_layer: float = DEFAULT_MAX_LAYER,
          intermediate: bool = True) -> LayerStack:
    """Sans inter-couches, `layer_height` (cible, par défaut le plafond) dimensionne N sur la
    colonne la plus haute. Avec, elle le dimensionne sur la colonne médiane (défaut 0,2 mm) et le
    plafond est tenu par subdivision locale."""
    if min_layer <= 0 or max_layer < min_layer:
        raise ValueError("hauteurs de couche : 0 < plancher ≤ plafond requis")
    target = min(layer_height or (0.2 if intermediate else max_layer), max_layer)
    target = max(target, min_layer)
    t1 = first_layer if first_layer is not None else min(max(0.2, min_layer), max_layer)
    hmax = float(np.nanmax(field.thickness))
    if hmax <= 0:
        raise ValueError("pièce d'épaisseur nulle")
    href = float(np.nanmedian(field.thickness)) if intermediate else hmax
    rest = max(href - t1, 0.0)
    n = 1 + (max(1, math.ceil(rest / target - 1e-9)) if rest > 0 else 0)
    stack = LayerStack(field, n, t1, min_layer, max_layer, intermediate)
    stack.count = int(np.nanmax(stack.last_grid()))
    return stack


def max_slope_deg(field: Heightfield, stack: LayerStack, k: int) -> float:
    """Pente maximale (degrés) de la surface S_k, mesurée sur la grille."""
    xx, yy = np.meshgrid(field.xs, field.ys)
    s = stack.z(k, xx.ravel(), yy.ravel()).reshape(xx.shape)
    s = np.where(np.isnan(field.top), np.nan, s)
    gy, gx = np.gradient(s, field.ys, field.xs)
    return float(np.degrees(np.arctan(np.nanmax(np.hypot(gx, gy)))))
