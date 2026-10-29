# Porte de pipeline/bird/skin/surface_patterns.py (branche bird_old, 2026-10-28) pour la voie chaleur : imports seulement.
# Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer).
# Porte de research/S4_Slicer/surface_patterns.py (2026-10-07, plan bird-into-apex-cpp) : code de calcul
# copie TEL QUEL, seuls les imports changent (egalite au bit avec la recherche).
# Fait partie de la voie Bird d'Apex : derive de S4 Slicer (Joshua Bird, GPL-3.0).
"""
Les 10 motifs de trace des couches de lissage non-planaire (et de la
premiere couche), calcules en 2D sur une zone plane : la zone de surface
DEPLIEE (surface_skin.py) ou le plan de la premiere couche.

Interface commune : generate(pattern, zone, spacing, angle=0.0) -> liste de
polylignes (K,2), dans l'ordre de trace. `zone` est un Polygon/MultiPolygon
shapely, `spacing` l'ecart entre passes voisines (largeur de buse pour une
couche pleine), `angle` l'orientation (radians) des motifs orientables.

Definitions (validees par l'utilisateur le 2026-09-27) :
  concentric           cercles imbriques autour du centre de la zone
  rectilinear          lignes paralleles reliees en zigzag
  monotonic_line       lignes paralleles toutes tracees dans le MEME sens,
                       dans l'ordre du balayage (pas de raccord)
  aligned_rectilinear  comme rectilinear, mais selon une direction GLOBALE
                       fixe (`angle`), identique pour toutes les zones, au
                       lieu de l'axe principal propre a chaque zone
  hilbert              courbe de Hilbert, pas = spacing
  archimedean          spirale d'Archimede, ecart entre spires = spacing
  concentric_square    carres imbriques autour du centre
  concentric_conform   contours de parois imbriques (decalages successifs du
                       bord) jusqu'au centre
  concentric_triangle  triangles equilateraux imbriques autour du centre
  octagram_spiral      spirale suivant une etoile a 8 branches

Centre des motifs imbriques/spirales : pole d'inaccessibilite de la zone
(shapely.ops.polylabel), toujours A L'INTERIEUR, contrairement au centroide
d'une zone concave. Toutes les formes sont rognees par la zone ; les bouts
plus courts que `spacing` sont ecartes.

Venv d'Apex (reutilise les balayages de pipeline/step6_toolpath2d).
"""
import math
import os
import sys

import numpy as np
from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from shapely.ops import polylabel

from pipeline.step6_toolpath2d import _clipped_scanlines, _connect_zigzag, _polygon_pieces  # noqa: E402

PATTERNS = (
    "concentric", "rectilinear", "monotonic_line", "aligned_rectilinear", "hilbert",
    "archimedean", "concentric_square", "concentric_conform", "concentric_triangle", "octagram_spiral",
)


def generate(pattern, zone, spacing, angle=0.0):
    if pattern not in PATTERNS:
        raise ValueError(f"motif inconnu : {pattern} (disponibles : {', '.join(PATTERNS)})")
    if zone is None or zone.is_empty or spacing <= 0:
        return []
    pieces = _polygon_pieces(zone)
    paths = []
    for piece in pieces:                   # chaque morceau a son propre centre
        paths.extend(_GENERATORS[pattern](piece, spacing, angle))
    return [p for p in paths if len(p) >= 2 and _length(p) >= spacing]


# ---------------------------------------------------------------------------
# outils communs
# ---------------------------------------------------------------------------

def _length(path):
    return float(np.linalg.norm(np.diff(path, axis=0), axis=1).sum())


def _clip(path, zone):
    """Rogne une polyligne par la zone, en gardant l'ordre de parcours."""
    clipped = LineString(path).intersection(zone)
    if clipped.is_empty:
        return []
    parts = clipped.geoms if hasattr(clipped, "geoms") else [clipped]
    out = [np.asarray(g.coords) for g in parts if g.geom_type == "LineString" and len(g.coords) >= 2]
    start = path[0]
    out.sort(key=lambda c: _arc_position(path, c[0], start))
    return out


def _arc_position(path, point, _start):
    """Abscisse curviligne approchee du point le plus proche sur la polyligne."""
    return LineString(path).project(Point(point))


def _center(zone):
    return np.asarray(polylabel(zone, tolerance=max(zone.length * 1e-4, 1e-6)).coords[0])


def _max_radius(zone, center):
    x0, y0, x1, y1 = zone.bounds
    corners = np.array([[x0, y0], [x0, y1], [x1, y0], [x1, y1]])
    return float(np.linalg.norm(corners - center, axis=1).max())


def _principal_angle(zone):
    """Axe le plus long du rectangle d'emprise minimal : le balayage suit la
    zone, ce qui limite le nombre de demi-tours."""
    rect = np.asarray(zone.minimum_rotated_rectangle.exterior.coords)[:4]
    edges = np.diff(np.vstack([rect, rect[:1]]), axis=0)
    longest = edges[np.argmax(np.linalg.norm(edges, axis=1))]
    return math.atan2(longest[1], longest[0])


def _nested(zone, spacing, outline):
    """Formes imbriquees autour du centre : outline(center, size) -> boucle
    fermee (M,2) de 'taille' size ; tailles espacees de `spacing`, rognees par
    la zone, de l'exterieur vers l'interieur."""
    center = _center(zone)
    reach = _max_radius(zone, center)
    paths = []
    for size in np.arange(reach, spacing * 0.5, -spacing):
        loop = outline(center, size)
        paths.extend(_clip(np.vstack([loop, loop[:1]]), zone))
    return paths


# ---------------------------------------------------------------------------
# les 10 motifs
# ---------------------------------------------------------------------------

def _concentric(zone, spacing, _angle):
    def circle(c, r):
        n = max(24, int(2 * math.pi * r / (spacing * 0.5)))
        t = np.linspace(0, 2 * math.pi, n, endpoint=False)
        return c + r * np.column_stack([np.cos(t), np.sin(t)])
    return _nested(zone, spacing, circle)


def _concentric_square(zone, spacing, angle):
    def square(c, half):
        u = np.array([math.cos(angle), math.sin(angle)])
        v = np.array([-u[1], u[0]])
        return np.array([c + half * (sx * u + sy * v) for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))])
    return _nested(zone, spacing, square)


def _concentric_triangle(zone, spacing, angle):
    # "taille" = rayon inscrit : deux triangles voisins ont leurs cotes
    # paralleles distants de exactement `spacing`.
    def triangle(c, inradius):
        circum = 2.0 * inradius
        t = angle + math.pi / 2 + np.arange(3) * 2 * math.pi / 3
        return c + circum * np.column_stack([np.cos(t), np.sin(t)])
    return _nested(zone, spacing, triangle)


def _concentric_conform(zone, spacing, _angle):
    paths, level = [], zone
    while not level.is_empty:
        for piece in _polygon_pieces(level):
            for ring in (piece.exterior, *piece.interiors):
                paths.append(np.asarray(ring.coords))
        level = level.buffer(-spacing, join_style="mitre")
    return paths


def _rectilinear(zone, spacing, angle, aligned=False):
    sweep = angle if aligned else _principal_angle(zone)
    scans, used_spacing = _clipped_scanlines(zone, 1.0, spacing, sweep)
    return _connect_zigzag(scans, zone, used_spacing)


def _aligned_rectilinear(zone, spacing, angle):
    return _rectilinear(zone, spacing, angle, aligned=True)


def _monotonic_line(zone, spacing, _angle):
    scans, _ = _clipped_scanlines(zone, 1.0, spacing, _principal_angle(zone))
    # _clipped_scanlines oriente deja chaque segment dans le sens du balayage :
    # on les garde tous dans ce sens, dans l'ordre des lignes.
    return [segment for scan in scans for segment in scan]


def _hilbert(zone, spacing, _angle):
    x0, y0, x1, y1 = zone.bounds
    extent = max(x1 - x0, y1 - y0)
    order = max(1, math.ceil(math.log2(max(extent / spacing, 1.0))))
    order = min(order, 10)                 # 1 048 576 points au plus
    n = 2 ** order
    step = extent / n if extent / n >= spacing * 0.999 else spacing
    d = np.arange(n * n)
    xs, ys = _hilbert_d2xy(n, d)
    # Grille calee sur le coin de l'emprise. Limite connue : une courbe sur
    # grille laisse une frange d'environ spacing/2 le long des bords (la ou
    # les autres motifs suivent le bord) -- ~95 % de couverture sur un disque
    # de 10 mm, 94,9 % sur un anneau. Centrer la grille a ete mesure PIRE
    # (94,1 %). La combler exigerait un contour ajoute, hors definition du motif.
    path = np.column_stack([x0 + (xs + 0.5) * step, y0 + (ys + 0.5) * step])
    return _clip(path, zone)


def _hilbert_d2xy(n, d):
    """Indice le long de la courbe -> (x, y) sur une grille n x n (vectorise)."""
    x = np.zeros_like(d)
    y = np.zeros_like(d)
    t = d.copy()
    s = 1
    while s < n:
        rx = 1 & (t // 2)
        ry = 1 & (t ^ rx)
        swap = ry == 0
        flip = swap & (rx == 1)
        x = np.where(flip, s - 1 - x, x)
        y = np.where(flip, s - 1 - y, y)
        x, y = np.where(swap, y, x), np.where(swap, x, y)
        x = x + s * rx
        y = y + s * ry
        t //= 4
        s *= 2
    return x, y


def _spiral(zone, spacing, profile):
    """Spirale partant du centre : rayon = spacing * tours * profile(theta),
    echantillonnee a ~spacing/2 le long de la courbe."""
    center = _center(zone)
    turns = _max_radius(zone, center) / spacing + 1
    theta_max = 2 * math.pi * turns
    n = max(64, int(math.pi * turns * turns * 2 * 2) + 64)   # longueur ~ pi*R^2/spacing
    theta = np.linspace(0, theta_max, n)
    r = spacing * theta / (2 * math.pi) * profile(theta)
    path = center + np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    return _clip(path, zone)


def _archimedean(zone, spacing, _angle):
    return _spiral(zone, spacing, lambda t: np.ones_like(t))


def _octagram_spiral(zone, spacing, angle):
    # profil en etoile a 8 branches : alterne entre 1 (pointes) et cos(pi/8)
    # (creux), lisse ; l'ecart minimal entre spires reste >= spacing * cos(pi/8).
    dip = 1.0 - math.cos(math.pi / 8)
    return _spiral(zone, spacing, lambda t: 1.0 - dip * (0.5 - 0.5 * np.cos(8 * (t - angle))))


_GENERATORS = {
    "concentric": _concentric,
    "rectilinear": _rectilinear,
    "monotonic_line": _monotonic_line,
    "aligned_rectilinear": _aligned_rectilinear,
    "hilbert": _hilbert,
    "archimedean": _archimedean,
    "concentric_square": _concentric_square,
    "concentric_conform": _concentric_conform,
    "concentric_triangle": _concentric_triangle,
    "octagram_spiral": _octagram_spiral,
}
