"""
element_ardor/step6_solid.py

Couches PLEINES de dessus et de dessous (etape 6, apres le remplissage).

Constat (Benchy, 2026-11-18) : Element Ardor ne connaissait aucune couche pleine. Chaque couche recevait ses parois et
un remplissage clairseme (20 %), y compris celles qui FORMENT une surface visible : le quadrillage du remplissage
etait a nu sur le pont, le toit, le dessus de la cheminee et du porte-canne, et sous le bateau. Un trancheur plan
remplit a 100 % les N couches sous chaque dessus et au-dessus de chaque dessous ; c'est ce que fait ce module, sur des
couches COURBES.

Pour une couche L de rang r (rang de son isovaleur) :
  - ses voisines de rang r + j (j = 1..N dessus, -1..-N dessous) sont projetees dans le repere local de L, comme L
    elle-meme. Une couche courbe est proche de son plan moyen ; deux couches empilees se ressemblent, la projection
    orthogonale sur le plan de L approche donc la projection le long de la normale aux couches ;
  - un morceau de couche voisine (une couche de l'etape 5 est un morceau d'isosurface) n'est retenu que s'il est
    PROCHE de L en 3D : boites englobantes qui se touchent a l'ecart d'isovaleur pres, et ecart au plan de L dans la
    plage des ecarts de L elle-meme. Sans ce filtre, une autre branche de la piece a la meme isovaleur, ailleurs,
    pourrait se projeter sur L et masquer un dessus ;
  - zone couverte = intersection, sur tous les j, des voisines retenues ; zone pleine = L moins la zone couverte. Un
    rang sans voisine (premiere couche, sommet) rend toute la couche pleine.

Dans la zone pleine (prise a l'interieur de la paroi la plus interieure), le remplissage clairseme est retire et
remplace par des lignes jointives (densite 1, en zigzag, orientation alternee d'une couche a l'autre comme le
remplissage). Les ponts
(type BRIDGE) ne sont pas touches. Les deux chemins de l'etape 6 (natif et Python) passent par ici : le resultat ne
depend pas du commutateur natif au-dela de ce qu'il en dependait deja.

Desactive (0 couche) dans `Step6Params` ; le pipeline l'active (`main.PipelineParams.toolpath2d`).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import shapely
import shapely.affinity
from shapely.geometry import LineString, Polygon

from core import line_types
from core.mesh_types import Layer
from element_ardor.step6_geometry import (
    _classify_contours,
    _interpolate_out_of_plane,
    _out_of_plane_offset,
    _polygon_pieces,
    _project_to_local_2d,
    _to_toolpath_3col,
)
from element_ardor.step6_infill import _generate_infill, _infill_angle_for_layer
from element_ardor.step6_bridge_runs import _split_rings_with_kinds
from element_ardor.step6_params import InfillPattern, Step6Params

_PROXIMITY_MARGIN = 1.0       # mm : marge du filtre de proximite 3D, au-dela de l'ecart d'isovaleur
_SAME_AXES = float(np.cos(np.radians(2.0)))   # deux reperes dont les axes sont a moins de 2 deg : meme projection
_MIN_SOLID_AREA_WIDTHS = 1.0  # une zone pleine plus etroite que cette fraction de buse est ignoree (eclats)


@dataclass(frozen=True)
class _Piece:
    """Morceau de couche de l'etape 5, vu par ce module : boucles 3D, boite englobante, isovaleur, rang."""
    index: int
    loops: list
    lo: np.ndarray
    hi: np.ndarray
    isovalue: float
    rank: int


def apply(layers: list[Layer], params: Step6Params, isovalue_rank: dict) -> list[Layer]:
    """Couches avec leurs zones pleines remplies a 100 % (sans effet si `solid_top_layers` et `solid_bottom_layers`
    sont nuls, ou sans isovaleurs), puis sans cordon de remplissage plus court que `min_fill_length`."""
    return [_drop_short_fill(layer, params.min_fill_length) for layer in _apply_solid(layers, params, isovalue_rank)]


def _drop_short_fill(layer: Layer, min_length: float) -> Layer:
    """Retire les cordons de remplissage (clairseme, plein, pont) plus courts que `min_length` mm (bouts de
    quelques dixiemes de millimetre entre deux parois, anneaux concentriques minuscules) : un trancheur classique ne
    les depose pas, la buse y laisse un pate. Les parois ne sont pas touchees."""
    if min_length <= 0 or not layer.toolpath_2d or layer.toolpath_kinds is None:
        return layer
    fill = (line_types.INFILL, line_types.BRIDGE)
    keep = [k not in fill or len(p) < 2 or float(np.hypot(*np.diff(np.asarray(p)[:, :2], axis=0).T).sum()) >= min_length
            for p, k in zip(layer.toolpath_2d, layer.toolpath_kinds)]
    if all(keep):
        return layer
    return replace(layer, toolpath_2d=[p for p, ok in zip(layer.toolpath_2d, keep) if ok],
                   toolpath_kinds=[k for k, ok in zip(layer.toolpath_kinds, keep) if ok])


def _apply_solid(layers: list[Layer], params: Step6Params, isovalue_rank: dict) -> list[Layer]:
    n_top, n_bottom = int(params.solid_top_layers), int(params.solid_bottom_layers)
    if (n_top <= 0 and n_bottom <= 0) or not layers:
        return layers
    pieces = _pieces(layers, isovalue_rank)
    regions = _Regions(layers)
    by_rank: dict[int, list[_Piece]] = {}
    for p in pieces:
        if p is not None:
            by_rank.setdefault(p.rank, []).append(p)
    out = []
    for layer, piece in zip(layers, pieces):
        if piece is None or not layer.toolpath_2d:
            out.append(layer)
            continue
        offsets = [j for j in range(1, n_top + 1)] + [-j for j in range(1, n_bottom + 1)]
        solid = _solid_region(layer, piece, by_rank, offsets, params, regions)
        out.append(_fill_solid(layer, solid, params, isovalue_rank, piece.rank == 0) if solid is not None else layer)
    return out


def _pieces(layers: list[Layer], isovalue_rank: dict) -> list:
    out = []
    for k, layer in enumerate(layers):
        loops = [np.asarray(c, dtype=float).reshape(-1, 3) for c in layer.contour_3d if len(c) >= 3]
        if layer.isovalue is None or layer.isovalue not in isovalue_rank or not loops:
            out.append(None)
            continue
        pts = np.concatenate(loops, axis=0)
        out.append(_Piece(k, loops, pts.min(0), pts.max(0), float(layer.isovalue), int(isovalue_rank[layer.isovalue])))
    return out


def _region(loops: list, frame) -> object:
    """Matiere d'un ensemble de boucles 3D projetees dans `frame` (imbrication -> trous), en un objet shapely."""
    groups = _classify_contours([_project_to_local_2d(c, frame) for c in loops])
    polys = []
    for outer, holes in groups:
        if len(outer) < 3:
            continue
        poly = Polygon(outer, holes=[h for h in holes if len(h) >= 3])
        if not poly.is_valid:
            poly = poly.buffer(0)
        polys.append(poly)
    return _valid(shapely.unary_union(polys)) if polys else Polygon()


def _valid(geom):
    """Geometrie valide et sans morceaux degeneres (lignes, points) : les operations suivantes n'echouent plus sur
    un contour qui se recoupe (\"side location conflict\" de GEOS, turkish)."""
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    if geom.geom_type == "GeometryCollection":
        geom = shapely.unary_union([g for g in geom.geoms if g.geom_type in ("Polygon", "MultiPolygon")])
    return geom


def _union(geoms: list):
    try:
        return shapely.unary_union(geoms)
    except shapely.errors.GEOSException:
        return shapely.unary_union([_valid(g.buffer(0)) for g in geoms])


def _intersection(a, b):
    try:
        return a.intersection(b)
    except shapely.errors.GEOSException:
        return _valid(a.buffer(0)).intersection(_valid(b.buffer(0)))


class _Regions:
    """Matiere de chaque couche dans son PROPRE repere, calculee une fois. Vue depuis une autre couche dont le repere a
    les memes axes (a 2 deg pres : toutes les couches planes, et la plupart des couches voisines courbes ; erreur de
    position < 0,2 mm pour 5 mm d'ecart hors plan, au bord des zones pleines seulement), c'est la meme region translatee : la projection et le classement des boucles (en O(n2) sur une couche de 240 morceaux) ne
    sont plus refaits pour chacune des voisines."""

    def __init__(self, layers: list[Layer]):
        self.layers = layers
        self.own: dict[int, object] = {}

    def own_region(self, k: int):
        if k not in self.own:
            loops = [np.asarray(c, dtype=float).reshape(-1, 3) for c in self.layers[k].contour_3d if len(c) >= 3]
            self.own[k] = _region(loops, self.layers[k].local_frame)
        return self.own[k]

    def seen_from(self, k: int, frame, loops: list):
        o, u, v = (np.asarray(a, dtype=float) for a in self.layers[k].local_frame)
        o2, u2, v2 = (np.asarray(a, dtype=float) for a in frame)
        if float(u @ u2) > _SAME_AXES and float(v @ v2) > _SAME_AXES:
            d = o - o2
            return shapely.affinity.translate(self.own_region(k), float(d @ u2), float(d @ v2))
        return _region(loops, frame)


def _near(piece: _Piece, other: _Piece, frame, own_range: tuple) -> bool:
    """`other` est-il la continuation de `piece` dans l'empilement (et non une autre branche de la piece) ?"""
    reach = abs(other.isovalue - piece.isovalue) + _PROXIMITY_MARGIN
    if np.any(other.lo > piece.hi + reach) or np.any(other.hi < piece.lo - reach):
        return False
    o, u, v = (np.asarray(x, dtype=float) for x in frame)
    mid = float((0.5 * (other.lo + other.hi) - o) @ np.cross(u, v))     # ecart du centre au plan de la couche
    return own_range[0] - reach <= mid <= own_range[1] + reach


def _solid_region(layer: Layer, piece: _Piece, by_rank: dict, offsets: list, params: Step6Params,
                  regions: "_Regions"):
    """Zone pleine de la couche (shapely, repere local), ou None si elle est entierement couverte."""
    frame = layer.local_frame
    own = regions.own_region(piece.index)
    if own.is_empty:
        return None
    own_off = np.concatenate([_out_of_plane_offset(c, frame) for c in piece.loops])
    own_range = (float(own_off.min()), float(own_off.max()))
    covered = own
    for j in offsets:
        neighbours = [o for o in by_rank.get(piece.rank + j, []) if _near(piece, o, frame, own_range)]
        if not neighbours:
            covered = Polygon()
            break
        covered = _intersection(covered, _union([regions.seen_from(o.index, frame, o.loops) for o in neighbours]))
        if covered.is_empty:
            break
    w = params.nozzle_width
    if covered.is_empty:
        solid = own
    else:
        try:
            solid = own.difference(covered)
        except shapely.errors.GEOSException:
            solid = _valid(own.buffer(0)).difference(_valid(covered.buffer(0)))
    inner = own.buffer(-w * params.wall_count, join_style="mitre")
    solid = solid.intersection(inner)
    # ouverture : retire les eclats plus etroits qu'une buse (paroi en pente douce, bruit de projection)
    solid = solid.buffer(-0.5 * w * _MIN_SOLID_AREA_WIDTHS).buffer(0.5 * w * _MIN_SOLID_AREA_WIDTHS).intersection(inner)
    return None if solid.is_empty or solid.area < w * w else solid


def _sources(layer: Layer):
    """Points de contour 2D et leurs ecarts hors-plan REELS : source de l'interpolation des points generes."""
    pts = [np.asarray(c, dtype=float).reshape(-1, 3) for c in layer.contour_3d if len(c) >= 2]
    if not pts:
        return None, None
    src2 = np.concatenate([_project_to_local_2d(c, layer.local_frame) for c in pts])
    off = np.concatenate([_out_of_plane_offset(c, layer.local_frame) for c in pts])
    return src2, off


_FLAT_OFFSET_TOL = 0.05    # mm : ecart hors-plan des contours sous lequel une couche est plane (rien a redecouper)


def _resample(line: np.ndarray, step: float) -> np.ndarray:
    """Polyligne 2D dont aucun segment ne depasse `step` (sommets d'origine gardes) ; inchangee si `step` <= 0."""
    if step <= 0 or len(line) < 2:
        return line
    seg = line[1:] - line[:-1]
    n = np.maximum(1, np.ceil(np.hypot(seg[:, 0], seg[:, 1]) / step).astype(int))
    t = np.concatenate([np.arange(k) / k for k in n])
    start = np.repeat(np.arange(len(seg)), n)
    return np.vstack([line[start] + t[:, None] * seg[start], line[-1:]])


def _clip_outside(poly: np.ndarray, solid) -> list:
    """Morceaux d'une polyligne (K,3) hors de `solid` ; l'ecart hors-plan est interpole le long de la polyligne."""
    line = LineString(poly[:, :2])
    if line.length <= 0 or not line.intersects(solid):
        return [poly]
    try:
        rest = line.difference(solid)
    except shapely.errors.GEOSException:
        return [poly]                       # coupe impossible : la ligne clairsemee reste entiere
    seg = np.hypot(*np.diff(poly[:, :2], axis=0).T)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    out = []
    for part in getattr(rest, "geoms", [rest]):
        if part.is_empty or part.geom_type != "LineString" or part.length < 1e-6:
            continue
        xy = np.asarray(part.coords)
        s = np.array([line.project(shapely.Point(p)) for p in xy])
        out.append(np.column_stack([xy, np.interp(s, cum, poly[:, 2])]))
    return out


def _piece_angle(piece, layer_angle: float) -> float:
    """Orientation des lignes pleines d'une zone : le long de son grand axe si elle est allongee (rapport > 2 : bande
    d'un dessus en pente, rebord), comme un trancheur classique ; sinon celle de la couche (alternee d'une couche a
    l'autre). En travers d'une bande de 1 mm, le zigzag faisait un virage tous les 0,4 mm (toit du Benchy)."""
    rect = piece.minimum_rotated_rectangle
    if rect.geom_type != "Polygon":
        return layer_angle
    c = np.asarray(rect.exterior.coords)[:4]
    e1, e2 = c[1] - c[0], c[2] - c[1]
    l1, l2 = float(np.hypot(*e1)), float(np.hypot(*e2))
    if min(l1, l2) <= 1e-9 or max(l1, l2) < 2.0 * min(l1, l2):
        return layer_angle
    e = e1 if l1 >= l2 else e2
    return float(np.arctan2(e[1], e[0]))


def _narrow_ring(piece, nozzle_width: float) -> bool:
    """Zone pleine etroite et non allongee (anneau du dessus d'une cheminee, d'un tube) : largeur moyenne 2A/P sous
    4 cordons et pas de grand axe. Un zigzag y laisse des bouts de quelques dixiemes de millimetre (dessus en
    « labyrinthe ») ; des anneaux concentriques la remplissent d'un seul trait par tour."""
    rect = piece.minimum_rotated_rectangle
    if rect.geom_type == "Polygon":
        c = np.asarray(rect.exterior.coords)[:4]
        l1, l2 = float(np.hypot(*(c[1] - c[0]))), float(np.hypot(*(c[2] - c[1])))
        if min(l1, l2) > 1e-9 and max(l1, l2) >= 2.0 * min(l1, l2):
            return False
    return piece.length > 0 and 2.0 * piece.area / piece.length < 4.0 * nozzle_width


def _concentric_lines(piece, nozzle_width: float) -> list:
    """Anneaux concentriques fermes, du bord vers l'interieur, a un cordon d'ecart (premier a un demi-cordon du bord)."""
    out = []
    for k in range(200):
        inner = _valid(piece.buffer(-nozzle_width * (0.5 + k)))
        if inner.is_empty:
            break
        for part in _polygon_pieces(inner):
            for ring in [part.exterior, *part.interiors]:
                coords = np.asarray(ring.coords)[:, :2]
                if len(coords) >= 4:
                    out.append(coords)
    return out


def _dense_lines(piece, nozzle_width: float, angle: float) -> list:
    """Lignes jointives d'une zone pleine. Une zone que GEOS ne sait pas balayer (segments presque confondus : « non-noded
    intersection », l-triple) est d'abord simplifiee de 1 um ; si elle echoue encore, elle reste vide
    plutot qu'un plantage du tranchage."""
    if _narrow_ring(piece, nozzle_width):
        return _concentric_lines(piece, nozzle_width)
    for zone in (piece, piece.simplify(1e-3).buffer(0)):
        try:
            return _generate_infill(zone, 1.0, InfillPattern.CROSSED_ZIGZAG, nozzle_width, angle)
        except shapely.errors.GEOSException:
            continue
    return []


def _fill_solid(layer: Layer, solid, params: Step6Params, isovalue_rank: dict, first: bool = False) -> Layer:
    """Retire le remplissage clairseme de `solid` et y ajoute des lignes jointives."""
    kinds = list(layer.toolpath_kinds) if layer.toolpath_kinds is not None else [line_types.NONE] * len(
        layer.toolpath_2d)
    # Premiere couche : posee sur le plateau, rien n'y est un pont (l'etape 2 marque pont toute face horizontale
    # tournee vers le bas sans rien dessous, le fond compris) ; ses cordons de pont sont remplaces comme le remplissage.
    sparse = (line_types.INFILL, line_types.BRIDGE) if first else (line_types.INFILL,)
    bridges = [] if first else [np.asarray(p)[:, :2] for p, k in zip(layer.toolpath_2d, kinds)
                                if k == line_types.BRIDGE and len(p) >= 2]
    if bridges:
        solid = solid.difference(shapely.unary_union([LineString(b).buffer(params.nozzle_width) for b in bridges]))
        if solid.is_empty:
            return layer
    kept, kept_kinds = [], []
    for poly, kind in zip(layer.toolpath_2d, kinds):
        poly = np.asarray(poly, dtype=float)
        if kind in sparse and len(poly) >= 2:
            parts = _clip_outside(poly, solid)
            kept += parts
            kept_kinds += [kind] * len(parts)
        else:
            kept.append(poly)
            kept_kinds.append(kind)
    angle = _infill_angle_for_layer(layer, isovalue_rank)
    src2, src_off = _sources(layer)
    power = params.out_of_plane_interpolation_power
    dense = []
    step = params.solid_line_resample if src_off is not None and np.ptp(src_off) > _FLAT_OFFSET_TOL else 0.0
    for piece in _polygon_pieces(_valid(solid.buffer(0))):
        for line in _dense_lines(piece, params.nozzle_width, _piece_angle(piece, angle)):
            line = _resample(np.asarray(line, dtype=float), step)
            off = (_interpolate_out_of_plane(line, src2, src_off, power)
                   if power is not None and src2 is not None else None)
            dense.append(_to_toolpath_3col(line, off))
    dense_kinds = [line_types.INFILL] * len(dense)
    if params.split_rings_into_runs and dense:
        dense, dense_kinds = _split_rings_with_kinds(
            dense, dense_kinds, include_open_polylines=params.split_open_polylines,
            share_vertices=params.continuous_rings, min_amplitude=params.run_split_min_amplitude,
            min_slope_deg=params.run_split_min_slope_deg)
    return replace(layer, toolpath_2d=list(dense) + kept, toolpath_kinds=list(dense_kinds) + kept_kinds)
