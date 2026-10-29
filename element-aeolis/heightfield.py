"""Modèle « surface en hauteur » : dessous B(x,y), dessus T(x,y) et empreinte en XY.

Un maillage entre dans cette classe si chaque rayon vertical qui le traverse le coupe en
exactement deux endroits (un dessous, un dessus) : plaque, relief, coque sans surplomb.
Les maillages qui ne le sont pas sont refusés (`NotAHeightfieldError`) plutôt que
tranchés faux.

Le dessous et le dessus sont rastérisés directement depuis les triangles (barycentriques),
sans lancer de rayons : exact pour des faces planes, indépendant de rtree/embree.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import shapely
import trimesh
from shapely.geometry import MultiPolygon, Polygon

_EDGE_EPS = 1e-9


class NotAHeightfieldError(ValueError):
    """Le maillage n'est pas un dessous + un dessus au-dessus d'une empreinte XY."""


@dataclass
class Heightfield:
    xs: np.ndarray            # (nx,) abscisses des nœuds de la grille
    ys: np.ndarray            # (ny,) ordonnées
    top: np.ndarray           # (ny, nx) T, NaN hors empreinte
    bottom: np.ndarray        # (ny, nx) B, NaN hors empreinte
    footprint: Polygon | MultiPolygon

    @property
    def thickness(self) -> np.ndarray:
        return self.top - self.bottom

    def _interp(self, grid: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Interpolation bilinéaire (les points hors grille sont ramenés au bord)."""
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        dx = self.xs[1] - self.xs[0]
        dy = self.ys[1] - self.ys[0]
        fx = np.clip((x - self.xs[0]) / dx, 0, len(self.xs) - 1 - 1e-9)
        fy = np.clip((y - self.ys[0]) / dy, 0, len(self.ys) - 1 - 1e-9)
        ix = fx.astype(int)
        iy = fy.astype(int)
        tx = fx - ix
        ty = fy - iy
        g00 = grid[iy, ix]
        g10 = grid[iy, ix + 1]
        g01 = grid[iy + 1, ix]
        g11 = grid[iy + 1, ix + 1]
        # Un nœud NaN (bord de l'empreinte) retombe sur ses voisins valides.
        stack = np.stack([g00, g10, g01, g11])
        w = np.stack([(1 - tx) * (1 - ty), tx * (1 - ty), (1 - tx) * ty, tx * ty])
        valid = ~np.isnan(stack)
        wsum = (w * valid).sum(axis=0)
        val = (np.where(valid, stack, 0.0) * w).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(wsum > 0, val / wsum, np.nan)

    def top_at(self, x, y):
        return self._interp(self.top, x, y)

    def bottom_at(self, x, y):
        return self._interp(self.bottom, x, y)


def _rasterize(tris: np.ndarray, xs: np.ndarray, ys: np.ndarray, keep_max: bool) -> np.ndarray:
    """Hauteur z des triangles aux nœuds (ys × xs) ; max (dessus) ou min (dessous)."""
    out = np.full((len(ys), len(xs)), np.nan)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    for t in tris:
        (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) = t
        det = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
        if abs(det) < 1e-14:
            continue
        i0 = max(int(np.ceil((min(x0, x1, x2) - xs[0]) / dx - _EDGE_EPS)), 0)
        i1 = min(int(np.floor((max(x0, x1, x2) - xs[0]) / dx + _EDGE_EPS)), len(xs) - 1)
        j0 = max(int(np.ceil((min(y0, y1, y2) - ys[0]) / dy - _EDGE_EPS)), 0)
        j1 = min(int(np.floor((max(y0, y1, y2) - ys[0]) / dy + _EDGE_EPS)), len(ys) - 1)
        if i0 > i1 or j0 > j1:
            continue
        gx, gy = np.meshgrid(xs[i0:i1 + 1], ys[j0:j1 + 1])
        a = ((y1 - y2) * (gx - x2) + (x2 - x1) * (gy - y2)) / det
        b = ((y2 - y0) * (gx - x2) + (x0 - x2) * (gy - y2)) / det
        c = 1.0 - a - b
        eps = 1e-9
        inside = (a >= -eps) & (b >= -eps) & (c >= -eps)
        if not inside.any():
            continue
        z = a * z0 + b * z1 + c * z2
        sub = out[j0:j1 + 1, i0:i1 + 1]
        cand = np.where(inside, z, np.nan)
        if keep_max:
            sub[:] = np.fmax(sub, cand)
        else:
            sub[:] = np.fmin(sub, cand)
    return out


def _footprint(mesh: trimesh.Trimesh) -> Polygon | MultiPolygon:
    """Silhouette en XY : union des triangles tournés vers le haut."""
    up = mesh.face_normals[:, 2] > 1e-9
    tris = mesh.triangles[up][:, :, :2]
    polys = shapely.polygons(tris)
    polys = polys[shapely.area(polys) > 1e-12]
    merged = shapely.union_all(shapely.buffer(polys, 1e-7))
    merged = shapely.buffer(merged, -1e-7)
    return shapely.simplify(merged, 1e-6)


def from_mesh(mesh: trimesh.Trimesh, cell: float = 0.25) -> Heightfield:
    """Construit le modèle depuis un maillage ; `cell` est le pas de la grille en mm."""
    n = mesh.face_normals
    lo, hi = mesh.bounds[0, :2], mesh.bounds[1, :2]
    nx = int(round((hi[0] - lo[0]) / cell)) + 1
    ny = int(round((hi[1] - lo[1]) / cell)) + 1
    xs = np.linspace(lo[0], hi[0], nx)
    ys = np.linspace(lo[1], hi[1], ny)
    top = _rasterize(mesh.triangles[n[:, 2] > 1e-9], xs, ys, keep_max=True)
    bottom = _rasterize(mesh.triangles[n[:, 2] < -1e-9], xs, ys, keep_max=False)
    covered = ~np.isnan(top)
    if not covered.any() or (covered != ~np.isnan(bottom)).any():
        raise NotAHeightfieldError("dessus et dessous ne couvrent pas la même empreinte")
    if (top[covered] - bottom[covered] < -1e-6).any():
        raise NotAHeightfieldError("le dessus passe sous le dessous (surplomb ou cavité)")
    field = Heightfield(xs, ys, top, bottom, _footprint(mesh))
    if mesh.is_watertight and mesh.volume > 0:
        # Un surplomb ou une cavité laisse du vide entre B et T : le volume intégré dépasse
        # alors celui du maillage.
        dx, dy = xs[1] - xs[0], ys[1] - ys[0]
        w = np.ones_like(top)
        w[[0, -1], :] *= 0.5
        w[:, [0, -1]] *= 0.5
        integral = float(np.nansum((top - bottom) * w) * dx * dy)
        if abs(integral - mesh.volume) > 0.02 * mesh.volume + 1e-6:
            raise NotAHeightfieldError(
                f"volume intégré {integral:.1f} mm³ ≠ volume du maillage {mesh.volume:.1f} mm³ "
                "(surplomb ou cavité)")
    return field
