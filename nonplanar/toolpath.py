"""Chemins de dépôt d'une couche : parois et remplissage en 2D, puis relevés sur la surface.

Une surface en hauteur se projette sans ambiguïté sur XY : parois et remplissage se calculent
dans le plan, puis chaque point reçoit z = S_k(x,y). Les segments sont redécoupés pour que Z
suive la surface.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon

from .layers import LayerStack


@dataclass
class ToolpathParams:
    line_width: float = 0.4
    wall_count: int = 2
    infill_density: float = 0.2
    solid_bottom_mm: float = 0.6      # profondeur remplie à 100 % depuis le dessous
    solid_top_mm: float = 0.8         # idem sous le dessus (mesurée en z, pas en nombre de couches)
    max_segment: float = 0.5          # longueur maximale d'un segment 3D relevé


@dataclass
class Polyline:
    xy: np.ndarray                    # (K, 2)
    kind: str                         # "wall_outer", "wall_inner", "infill"


@dataclass
class LayerPath:
    index: int                        # k, de 1 à N
    polylines: list[Polyline] = field(default_factory=list)


def _rings(geom) -> list[np.ndarray]:
    polys = [geom] if isinstance(geom, Polygon) else list(getattr(geom, "geoms", []))
    out = []
    for p in polys:
        if not isinstance(p, Polygon) or p.is_empty:
            continue
        out.append(np.asarray(p.exterior.coords))
        out.extend(np.asarray(r.coords) for r in p.interiors)
    return out


def _zigzag(region, spacing: float, angle_deg: float) -> list[np.ndarray]:
    """Lignes parallèles balayant `region`, raccordées en zigzag quand le raccord reste dedans."""
    if region.is_empty:
        return []
    c = np.array(region.centroid.coords[0])
    a = np.radians(angle_deg)
    d = np.array([np.cos(a), np.sin(a)])
    n = np.array([-d[1], d[0]])
    minx, miny, maxx, maxy = region.bounds
    reach = float(np.hypot(maxx - minx, maxy - miny)) / 2 + spacing
    lines = []
    for i, off in enumerate(np.arange(-reach, reach + spacing, spacing)):
        p0 = c + n * off - d * reach
        p1 = c + n * off + d * reach
        inter = region.intersection(LineString([p0, p1]))
        segs = [inter] if isinstance(inter, LineString) else list(getattr(inter, "geoms", []))
        for s in segs:
            if isinstance(s, LineString) and s.length > 1e-6:
                pts = np.asarray(s.coords)
                lines.append(pts if (len(lines) % 2 == 0) else pts[::-1])
    # raccords : on chaîne les lignes consécutives dont le raccord est couvert par la région
    chains: list[list[np.ndarray]] = []
    safe = region.buffer(1e-6)
    for pts in lines:
        if chains:
            end = chains[-1][-1][-1]
            if np.linalg.norm(end - pts[0]) < 3 * spacing and safe.contains(LineString([end, pts[0]])):
                chains[-1].append(pts)
                continue
        chains.append([pts])
    return [np.vstack(c) for c in chains]


def layer_polylines(region, solid_region, index: int, params: ToolpathParams) -> LayerPath:
    """Parois (de l'intérieur vers l'extérieur) puis remplissage, en 2D.

    `region` : zone où la couche existe ; `solid_region` : partie de `region` à remplir à 100 %
    (couches du bas et du haut de chaque colonne).
    """
    w = params.line_width
    path = LayerPath(index)
    n_walls = params.wall_count
    inner = region.buffer(-(n_walls * w))
    angle = 45.0 if index % 2 else 135.0
    if solid_region is None or solid_region.is_empty:
        solid, sparse = inner.buffer(0).difference(inner.buffer(0)), inner
    else:
        solid = inner.intersection(solid_region.buffer(0))
        sparse = inner.difference(solid)
    sparse_spacing = w / max(min(params.infill_density, 1.0), 1e-3)
    for geom, spacing in ((sparse, sparse_spacing), (solid, w)):
        for line in _zigzag(geom, spacing, angle):
            path.polylines.append(Polyline(line, "infill"))
    for i in range(n_walls - 1, -1, -1):
        for ring in _rings(region.buffer(-(i + 0.5) * w)):
            path.polylines.append(Polyline(ring, "wall_outer" if i == 0 else "wall_inner"))
    return path


def densify(xy: np.ndarray, stack: LayerStack, k: int, max_segment: float, j: int = 0) -> np.ndarray:
    """Relève (K,2) sur la surface de la couche (k, j) et redécoupe : retourne (M,3)."""
    out = [xy[0:1]]
    for a, b in zip(xy[:-1], xy[1:]):
        length = float(np.linalg.norm(b - a))
        m = max(1, int(np.ceil(length / max_segment)))
        t = np.linspace(0, 1, m + 1)[1:, None]
        out.append(a + (b - a) * t)
    pts = np.vstack(out)
    z = stack.z(k, pts[:, 0], pts[:, 1]) if j == 0 else stack.z_sub(k, j, pts[:, 0], pts[:, 1])
    return np.column_stack([pts, z])
