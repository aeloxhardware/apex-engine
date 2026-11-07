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

from core import progress as _progress

_EDGE_EPS = 1e-9
_FOOTPRINT_CHUNKS = 16
_FOOTPRINT_MIN_CHUNK = 500


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


_PAIR_BUDGET = 2_000_000    # (triangle, nœud) évalués par paquet : borne la mémoire
_MIN_DET = 1e-14


def _rasterize(tris: np.ndarray, xs: np.ndarray, ys: np.ndarray, keep_max: bool,
               span: tuple[float, float] | None = None) -> np.ndarray:
    """Hauteur z des triangles aux nœuds (ys × xs) ; max (dessus) ou min (dessous).

    Vectorisé sur les couples (triangle, nœud de sa boîte englobante) : mêmes formules barycentriques
    que le balayage triangle par triangle, donc mêmes valeurs ; le max/min ne dépend pas de l'ordre.
    `span` : part [lo, hi] de la barre de l'étape que ce balayage occupe (progression réelle)."""
    out = np.full((len(ys), len(xs)), np.nan)
    if len(tris) == 0:
        return out
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    t = np.asarray(tris, float)
    x0, y0, z0 = t[:, 0, 0], t[:, 0, 1], t[:, 0, 2]
    x1, y1, z1 = t[:, 1, 0], t[:, 1, 1], t[:, 1, 2]
    x2, y2, z2 = t[:, 2, 0], t[:, 2, 1], t[:, 2, 2]
    det = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
    xmin, xmax = np.minimum(np.minimum(x0, x1), x2), np.maximum(np.maximum(x0, x1), x2)
    ymin, ymax = np.minimum(np.minimum(y0, y1), y2), np.maximum(np.maximum(y0, y1), y2)
    i0 = np.maximum(np.ceil((xmin - xs[0]) / dx - _EDGE_EPS), 0).astype(int)
    i1 = np.minimum(np.floor((xmax - xs[0]) / dx + _EDGE_EPS), len(xs) - 1).astype(int)
    j0 = np.maximum(np.ceil((ymin - ys[0]) / dy - _EDGE_EPS), 0).astype(int)
    j1 = np.minimum(np.floor((ymax - ys[0]) / dy + _EDGE_EPS), len(ys) - 1).astype(int)
    ok = (np.abs(det) >= _MIN_DET) & (i0 <= i1) & (j0 <= j1)
    idx = np.flatnonzero(ok)
    w = (i1 - i0 + 1)[idx]
    h = (j1 - j0 + 1)[idx]
    cost = np.cumsum(w * h)
    reduce_at = np.fmax.at if keep_max else np.fmin.at
    eps = 1e-9
    start = 0
    while start < len(idx):
        base = cost[start - 1] if start else 0
        stop = max(int(np.searchsorted(cost, base + _PAIR_BUDGET, side="right")), start + 1)
        sel = idx[start:stop]
        ww, hh = w[start:stop], h[start:stop]
        n_pairs = ww * hh
        tri = np.repeat(sel, n_pairs)
        first = np.cumsum(n_pairs) - n_pairs
        local = np.arange(int(n_pairs.sum())) - np.repeat(first, n_pairs)
        wr = np.repeat(ww, n_pairs)
        ix = i0[tri] + local % wr
        iy = j0[tri] + local // wr
        gx, gy = xs[ix], ys[iy]
        d = det[tri]
        a = ((y1[tri] - y2[tri]) * (gx - x2[tri]) + (x2[tri] - x1[tri]) * (gy - y2[tri])) / d
        b = ((y2[tri] - y0[tri]) * (gx - x2[tri]) + (x0[tri] - x2[tri]) * (gy - y2[tri])) / d
        c = 1.0 - a - b
        inside = (a >= -eps) & (b >= -eps) & (c >= -eps)
        z = a * z0[tri] + b * z1[tri] + c * z2[tri]
        reduce_at(out, (iy[inside], ix[inside]), z[inside])
        start = stop
        if span is not None:
            _progress.report(start, len(idx), span)
    return out


def _footprint(mesh: trimesh.Trimesh, span: tuple[float, float] | None = None) -> Polygon | MultiPolygon:
    """Silhouette en XY : union des triangles tournés vers le haut.

    Les triangles du dessus d'une surface en hauteur se partagent leurs arêtes sans se recouvrir : l'union
    « de couverture » de GEOS les réunit en ~0,1 s, contre ~6 s pour l'union générale sur 54 000
    triangles. Elle ne vaut que pour une couverture valide : on le vérifie (l'aire de l'union égale la
    somme des aires, et le résultat est valide) et, sinon, on retombe sur l'union robuste (triangles élargis de 1e-7 mm), faite par
    paquets pour garder un compte rendu de progression réel."""
    up = mesh.face_normals[:, 2] > 1e-9
    tris = mesh.triangles[up][:, :, :2]
    polys = shapely.polygons(tris)
    polys = polys[shapely.area(polys) > 1e-12]
    merged = None
    try:
        candidate = shapely.coverage_union_all(polys)
        if candidate.is_valid and \
                abs(float(shapely.area(candidate)) - float(shapely.area(polys).sum())) <= 1e-9 * max(1.0, candidate.area):
            merged = candidate
    except shapely.errors.GEOSException:
        pass
    if merged is None:
        polys = shapely.buffer(polys, 1e-7)
        n_chunks = max(1, min(_FOOTPRINT_CHUNKS, len(polys) // _FOOTPRINT_MIN_CHUNK))
        parts = []
        for i, chunk in enumerate(np.array_split(polys, n_chunks)):
            parts.append(shapely.union_all(chunk))
            if span is not None:
                _progress.report(i + 1, n_chunks + 1, span)
        merged = shapely.buffer(shapely.union_all(parts), -1e-7)
    return shapely.simplify(merged, 1e-6)


def from_mesh(mesh: trimesh.Trimesh, cell: float = 0.25) -> Heightfield:
    """Construit le modèle depuis un maillage ; `cell` est le pas de la grille en mm."""
    n = mesh.face_normals
    lo, hi = mesh.bounds[0, :2], mesh.bounds[1, :2]
    nx = int(round((hi[0] - lo[0]) / cell)) + 1
    ny = int(round((hi[1] - lo[1]) / cell)) + 1
    xs = np.linspace(lo[0], hi[0], nx)
    ys = np.linspace(lo[1], hi[1], ny)
    top = _rasterize(mesh.triangles[n[:, 2] > 1e-9], xs, ys, keep_max=True, span=(0.0, 0.1))
    bottom = _rasterize(mesh.triangles[n[:, 2] < -1e-9], xs, ys, keep_max=False, span=(0.1, 0.2))
    covered = ~np.isnan(top)
    if not covered.any() or (covered != ~np.isnan(bottom)).any():
        raise NotAHeightfieldError("dessus et dessous ne couvrent pas la même empreinte")
    if (top[covered] - bottom[covered] < -1e-6).any():
        raise NotAHeightfieldError("le dessus passe sous le dessous (surplomb ou cavité)")
    field = Heightfield(xs, ys, top, bottom, _footprint(mesh, span=(0.2, 0.95)))
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
