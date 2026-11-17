"""Tranchage GÉNÉRAL d'Element Aeolis : n'importe quelle géométrie, lissage surfacique dans la limite de la buse.

Le moteur « surface en hauteur » (`heightfield.py`, `layers.py`) ne sait trancher que les maillages dont chaque
rayon vertical coupe exactement un dessous et un dessus. Celui-ci prend tout le reste (surplombs, cavités, pièces
creuses, plusieurs corps) : il TRANCHE PAR SECTIONS PLANAIRES DÉFORMÉES.

Principe. La couche k est la section du maillage par la surface

    z = F_k(x, y) = c_k + δ_k(x, y)

où c_k est un plan (c_k = z_bas + t1 + (k−1)·h) et δ_k le RELIEF LOCAL de la couche (`relief.py`), qui descend ou
monte la couche. Sans relief
(δ = 0), c'est un tranchage planaire classique ; avec lui, CHAQUE surface tournée vers le haut (pont, toit, dessus
d'une caisse), à n'importe quelle hauteur, est une couche à elle seule au lieu de marches d'escalier. Comme toute
section planaire, cela vaut pour n'importe quelle topologie (une colonne verticale peut traverser plusieurs
intervalles de matière).

Limites de la buse. Deux réglages de l'imprimante :

  * `nozzle_angle_deg` : pente maximale, depuis l'horizontale, que la surface déjà imprimée peut avoir autour de la
    pointe (0° = épaule plate, 90° = aiguille) ;
  * `nozzle_reach` : longueur de buse sous le conduit de ventilation, donc dénivelé maximal toléré entre la pointe
    et la matière voisine plus haute (mm).

Le relief vérifie, dans chaque couche, une pente ≤ tan α (et ≤ 70°), et max δ − min δ ≤ longueur de buse. Comme la
matière déjà posée en q est sous F_{k−1}(q) ≤ F_k(q) et que la pointe est en F_k(p), la matière voisine ne peut jamais
dépasser la pointe de plus que min(r·tan α, longueur) : la buse ne touche rien, PAR CONSTRUCTION (le contrôle de
collision ne fait que le vérifier).

Épaisseur. Elle vaut h + δ_k − δ_{k−1}, dans [0, 1,25·h] : sous une surface épousée, seules les couches de sa peau la
suivent, celles du dessous s'écrasent (épaisseur nulle) au lieu de se tordre. Une couche écrasée ou plus mince que
le plancher n'est pas imprimée ; sa tranche revient à la couche suivante (`GeneralStack.base`), sans dépasser le
plafond.

Historique. Jusqu'au 2026-11-14 le relief était UNIQUE : D(x, y), calé sur les seules colonnes les plus hautes et
appliqué à toute la hauteur par une rampe λ(c). Mesuré sur un Benchy (80°, 10 mm) : marches sur toutes les surfaces
plus basses que le toit, et cônes à 80° creusés sous le toit à travers toute la pièce (parois jusqu'à 1,5 mm de la
surface du modèle).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import shapely
import trimesh
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon

from core import progress as _progress
from . import headcheck as hc
from . import heightfield as hf
from .regions import level_geometry
from . import toolpath as tp
from .relief import Relief, amplitude, build_relief, cone_minorant  # noqa: F401  (`cone_minorant` : nom historique)
from .relief import _AMP_FRACTION, _CHAMFER_SAFETY

RELIEF_CELL = 0.25           # pas de la grille du relief et de la transformée de distance (mm)
_MAX_EDGE = 1.0              # côté xy maximal des triangles qu'on déforme (mm)
_MAX_FACES = 700_000


class NotSliceableError(ValueError):
    """Le maillage n'a donné aucune couche (vide, ouvert, dégénéré)."""


# --------------------------------------------------------------------------------------- piles de couches

class GeneralStack:
    """Même interface que `LayerStack` pour le G-code : `sample(x, y)` → colonne avec `layer_z`, `layer_thickness`."""

    def __init__(self, relief: Relief, z_bottom: float, planes: np.ndarray, min_thickness: float = 0.0,
                 max_thickness: float = math.inf):
        self.relief = relief
        self.z_bottom = z_bottom
        self.c = planes
        self.count = len(planes)
        self.base = self._bases(min_thickness, max_thickness)

    def _bases(self, min_thickness: float, max_thickness: float) -> np.ndarray:
        """Dessus de la matière posée sous chaque couche, aux nœuds : la surface de la dernière couche IMPRIMÉE en
        dessous. Une couche plus mince que `min_thickness` (écrasée sous une surface épousée) n'est pas imprimée ; sa
        tranche revient à la couche suivante, qui la dépose en plus de la sienne. Sans ce report, chaque tranche trop
        mince laissait un vide, visible en encoches sur les parois (flanc d'un Benchy sous le plat-bord). Le report ne
        doit pas porter la couche suivante au-delà du plafond `max_thickness` : la tranche est alors imprimée mince
        (comptée dans `fraction_below_min`), une buse trop haute sur la matière étant pire qu'une couche trop fine."""
        base = np.empty_like(self.relief.d)
        self.kept = np.zeros(self.relief.d.shape, dtype=bool)
        floor = np.full(self.relief.d.shape[1:], self.z_bottom, dtype=self.relief.d.dtype)
        for k in range(self.count):
            base[k] = floor
            zk = self.c[k] + self.relief.d[k]
            nxt = self.c[k + 1] + self.relief.d[k + 1] if k + 1 < self.count else zk
            keep = (zk - floor >= min_thickness) | ((zk - floor > 1e-6) & (nxt - floor > max_thickness))
            self.kept[k] = keep
            floor = np.where(keep, zk, floor)
        return base

    def _printed(self, g: np.ndarray, k: int, where: tuple) -> np.ndarray:
        """Interpolation bilinéaire de la grille `g` de la couche k, limitée aux nœuds où la couche est IMPRIMÉE.

        À côté d'un petit détail épousé, la même couche est écrasée un ou deux millimètres plus bas (pont autour de
        l'anneau incliné d'un porte-canne de Benchy, 0,8 mm de large, soit trois nœuds) : une interpolation ordinaire
        mêlait les deux hauteurs et les cordons du bord de l'anneau ondulaient en bosses. Sans nœud imprimé autour du
        point, interpolation ordinaire."""
        ix, iy, tx, ty = where
        corners = ((iy, ix, (1 - tx) * (1 - ty)), (iy, ix + 1, tx * (1 - ty)),
                   (iy + 1, ix, (1 - tx) * ty), (iy + 1, ix + 1, tx * ty))
        m = self.kept[min(max(k, 1), self.count) - 1]
        num = np.zeros(np.shape(tx)); den = np.zeros(np.shape(tx)); plain = np.zeros(np.shape(tx))
        for jy, jx, w in corners:
            v = g[jy, jx]
            plain += w * v
            wm = w * m[jy, jx]
            num += wm * v
            den += wm
        return np.where(den > 1e-9, num / np.maximum(den, 1e-12), plain)

    def base_xy(self, k: int, where: tuple) -> np.ndarray:
        """Dessus de la matière posée sous la couche k (à partir de 1) aux points déjà localisés."""
        return self._printed(self.base[k - 1], k, where)

    def z_xy(self, k: int, x=None, y=None, *, where: tuple | None = None) -> np.ndarray:
        """Hauteur de la couche k (à partir de 1) aux points (x, y), ou aux points déjà localisés (`where`)."""
        where = where if where is not None else self.relief.locate(x, y)
        return self.c[k - 1] + self._printed(self.relief.d[min(max(k, 1), self.count) - 1], k, where)

    def z_grid(self, k: int) -> np.ndarray:
        """Hauteur de la couche k aux nœuds de la grille du relief."""
        return self.c[k - 1] + self.relief.d[k - 1]

    def sample(self, x, y) -> "_Column":
        return _Column(self, self.relief.locate(x, y))


class _Column:
    def __init__(self, stack: GeneralStack, where: tuple):
        self.stack, self.where = stack, where

    def layer_z(self, k: int, j: int = 0) -> np.ndarray:
        return self.stack.z_xy(k, where=self.where)

    def layer_thickness(self, k: int, j: int = 0) -> np.ndarray:
        """Épaisseur DÉPOSÉE : jusqu'au dessus de la dernière couche imprimée dessous (`GeneralStack.base`)."""
        return self.layer_z(k) - self.stack.base_xy(k, self.where)


# --------------------------------------------------------------------------------------- sections

def _refine(vertices: np.ndarray, faces: np.ndarray, max_edge: float) -> tuple[np.ndarray, np.ndarray]:
    """Coupe en deux, au milieu, chaque arête dont la longueur EN PROJECTION XY dépasse `max_edge` (la déformation ne
    dépend que de x, y : un triangle haut et mince est déjà exact), et redécoupe chaque triangle selon ses arêtes
    coupées (1, 2 ou 3). Les deux triangles d'une arête partagent son milieu : le maillage reste étanche, sans
    jonction en T (`trimesh.remesh.subdivide(face_index=...)` en laisse, ce qui ouvre les contours des sections).
    Borné à `_MAX_FACES` : le côté maximal double sinon."""
    v = np.asarray(vertices, float)
    f = np.asarray(faces, dtype=np.int64)
    for _ in range(16):
        edges = np.stack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]], axis=1)             # (F, 3, 2)
        length = np.linalg.norm(v[edges[..., 0], :2] - v[edges[..., 1], :2], axis=2)     # (F, 3)
        mark = length > max_edge
        if not mark.any():
            break
        if len(f) + 3 * int(mark.any(axis=1).sum()) > _MAX_FACES:
            max_edge *= 2.0
            continue
        n = len(v)
        key = edges.min(axis=2) * n + edges.max(axis=2)
        uniq = np.unique(key[mark])
        mid = np.where(mark, n + np.searchsorted(uniq, key), -1)
        v = np.vstack([v, 0.5 * (v[uniq // n] + v[uniq % n])])
        cnt = mark.sum(axis=1)
        parts = [f[cnt == 0]]

        def rot(sel, r):
            return f[sel, r], f[sel, (r + 1) % 3], f[sel, (r + 2) % 3]

        for r in range(3):                       # une arête coupée : e_r = (a, b)
            sel = (cnt == 1) & mark[:, r]
            a, b, c = rot(sel, r)
            m = mid[sel, r]
            parts += [np.column_stack([a, m, c]), np.column_stack([m, b, c])]
        for r in range(3):                       # deux arêtes coupées : e_r = (a, b) et e_{r+1} = (b, c)
            sel = (cnt == 2) & mark[:, r] & mark[:, (r + 1) % 3]
            a, b, c = rot(sel, r)
            m0, m1 = mid[sel, r], mid[sel, (r + 1) % 3]
            parts += [np.column_stack([m0, b, m1]), np.column_stack([a, m0, m1]), np.column_stack([a, m1, c])]
        sel = cnt == 3
        a, b, c = f[sel, 0], f[sel, 1], f[sel, 2]
        m0, m1, m2 = mid[sel, 0], mid[sel, 1], mid[sel, 2]
        parts += [np.column_stack([a, m0, m2]), np.column_stack([m0, b, m1]),
                  np.column_stack([m2, m1, c]), np.column_stack([m0, m1, m2])]
        f = np.vstack(parts)
    return v, f


def _edge_point(v: np.ndarray, f: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Point d'annulation de f sur l'arête (i, j), calculé dans un ordre canonique : les deux triangles qui
    partagent l'arête obtiennent EXACTEMENT le même point, donc les boucles se referment sans tolérance."""
    lo, hi = np.minimum(i, j), np.maximum(i, j)
    fl, fh = f[lo], f[hi]
    t = fl / (fl - fh)
    return v[lo, :2] + (v[hi, :2] - v[lo, :2]) * t[:, None]


def section_region(v: np.ndarray, faces: np.ndarray, f: np.ndarray):
    """Région (shapely) où la fonction f (par sommet) est ≤ 0… plus exactement le contour {f = 0} rempli en
    parité paire-impaire : une cavité reste un trou, un îlot dans une cavité redevient de la matière."""
    s = f[faces] > 0
    cnt = s.sum(axis=1)
    sel = np.flatnonzero((cnt == 1) | (cnt == 2))
    if sel.size == 0:
        return Polygon()
    tri, sm = faces[sel], s[sel]
    odd = np.where(cnt[sel] == 1, sm.argmax(axis=1), (~sm).argmax(axis=1))
    r = np.arange(len(sel))
    a, b, c = tri[r, odd], tri[r, (odd + 1) % 3], tri[r, (odd + 2) % 3]
    p0, p1 = _edge_point(v, f, a, b), _edge_point(v, f, a, c)
    keep = np.any(p0 != p1, axis=1)
    if not keep.any():
        return Polygon()
    lines = shapely.linestrings(np.stack([p0[keep], p1[keep]], axis=1))
    noded = shapely.union_all(lines)
    faces_ = shapely.get_parts(shapely.polygonize(shapely.get_parts(noded)))
    if len(faces_) == 0:
        return Polygon()
    filled = np.array([Polygon(p.exterior) for p in faces_ if p.area > 1e-9], dtype=object)
    if len(filled) == 0:
        return Polygon()
    region = filled[0]
    for poly in filled[1:]:                 # parité paire-impaire : ou exclusif de toutes les boucles remplies
        region = region.symmetric_difference(poly)
    return _polygons_only(region)


def _polygons_only(geom):
    if geom.is_empty:
        return Polygon()
    if isinstance(geom, (Polygon, MultiPolygon)):
        return geom if geom.is_valid else geom.buffer(0)
    parts = [g for g in getattr(geom, "geoms", []) if isinstance(g, (Polygon, MultiPolygon)) and not g.is_empty]
    return shapely.union_all(parts) if parts else Polygon()


# --------------------------------------------------------------------------------------- plan

@dataclass
class GeneralPlan:
    relief: Relief
    stack: GeneralStack
    vertices: np.ndarray
    faces: np.ndarray
    thickness: float                          # h : espacement nominal des plans
    field: hf.Heightfield                     # dessus/dessous rastérisés (étendue pour le contrôle de la tête)
    z_bottom: float
    regions: list = field(default_factory=list)
    sections: list = field(default_factory=list)  # sections brutes, couches écrasées comprises (peau pleine)
    slope_deg: float = 0.0
    min_thickness: float = 0.0                # plancher : une couche plus mince (écrasée) n'est pas imprimée

    @property
    def count(self) -> int:
        return self.stack.count


def analyze(mesh: trimesh.Trimesh, *, nozzle_angle_deg: float, nozzle_reach: float, layer_height: float | None,
            first_layer: float | None, min_layer: float, max_layer: float) -> GeneralPlan:
    """Dessus de la pièce, relief doux et plans des couches. (Étape « heightfield » : analyse du modèle.)"""
    if len(mesh.faces) == 0:
        raise NotSliceableError("maillage vide")
    lo, hi = mesh.bounds[0, :2], mesh.bounds[1, :2]
    nx = int(round((hi[0] - lo[0]) / RELIEF_CELL)) + 1
    ny = int(round((hi[1] - lo[1]) / RELIEF_CELL)) + 1
    xs, ys = np.linspace(lo[0], hi[0], max(nx, 2)), np.linspace(lo[1], hi[1], max(ny, 2))
    n = mesh.face_normals
    top = hf._rasterize(mesh.triangles[n[:, 2] > 1e-9], xs, ys, keep_max=True, span=(0.0, 0.3))
    bottom = hf._rasterize(mesh.triangles[n[:, 2] < -1e-9], xs, ys, keep_max=False, span=(0.3, 0.6))
    if np.isnan(top).all():
        raise NotSliceableError("aucune face tournée vers le haut")
    z_bottom = float(mesh.bounds[0, 2])
    z_top = float(mesh.bounds[1, 2])
    height = z_top - z_bottom
    amp = amplitude(height, nozzle_angle_deg, nozzle_reach)
    # Épaisseur nominale : avec relief, elle varie de ±25 % ; elle doit rester dans [plancher, plafond].
    spread = 1.0 + (_AMP_FRACTION if amp > 0 else 0.0)
    h = layer_height if layer_height else 0.2
    h = min(max(h, min_layer * spread), max_layer / spread) if max_layer / spread >= min_layer * spread else max(min_layer, min(h, max_layer))
    t1 = first_layer if first_layer is not None else min(max(0.2, min_layer), max_layer)
    rest = max(height - t1, 0.0)
    n_main = 1 + (max(1, int(round(rest / h))) if rest > 1e-9 else 0)
    h_eff = rest / (n_main - 1) if n_main > 1 else h
    # Plans jusqu'à z_top + amp : une colonne dont le dessus est sous la surface du relief finit plus tard.
    n_total = n_main + int(math.ceil(amp / max(h_eff, 1e-9))) + 1
    planes = z_bottom + t1 + h_eff * np.arange(n_total)
    relief = build_relief(mesh, (xs, ys), planes, nozzle_angle_deg, nozzle_reach, min_layer)
    stack = GeneralStack(relief, z_bottom, planes, min_layer, max_layer)
    v, f = _refine(mesh.vertices, mesh.faces, _MAX_EDGE)
    field_ = hf.Heightfield(xs, ys, top, bottom, None)   # type: ignore[arg-type]
    slope = math.degrees(math.atan(relief.max_gradient)) if relief.amp > 0 else 0.0
    return GeneralPlan(relief, stack, v, f, h_eff, field_, z_bottom, slope_deg=slope, min_thickness=min_layer)


def compute_regions(plan: GeneralPlan) -> None:
    """Région de chaque couche : section du maillage par la surface F_k (à mi-épaisseur), privée des zones où la
    couche est écrasée (épaisseur sous `min_thickness`, sous une surface épousée) : rien n'y est déposé. Progression
    réelle."""
    st, rel = plan.stack, plan.relief
    where = rel.locate(plan.vertices[:, 0], plan.vertices[:, 1])
    z = plan.vertices[:, 2]
    below = np.zeros(len(z))                      # δ de la couche du dessous aux sommets (plateau : 0)
    plan.regions, plan.sections = [], []
    for i in range(st.count):
        _progress.report(i, st.count)
        c_mid = st.c[i] - (plan.stack.c[0] - plan.z_bottom if i == 0 else plan.thickness) / 2.0
        here = st._printed(rel.d[i], i + 1, where)      # nœuds imprimés seulement (voir `_printed`)
        f = z - c_mid - 0.5 * (here + below)      # surface à mi-épaisseur : moyenne des deux couches
        below = here
        region = section_region(plan.vertices, plan.faces, f)
        plan.sections.append(region if region.area >= 1e-3 else Polygon())
        thick = st.z_grid(i + 1) - st.base[i]
        kept = (st.base[i + 1] > st.base[i]) if i + 1 < st.count else thick > 1e-6
        if not region.is_empty and not kept.all():
            # Bord tracé sur l'épaisseur interpolée, pas sur le masque des nœuds : un masque coupe en escalier au pas
            # de la grille (encoches de 0,1 mm sur les parois du plat-bord d'un Benchy).
            level = np.where(kept, np.maximum(thick, plan.min_thickness), np.minimum(thick, plan.min_thickness - 1e-4))
            region = _polygons_only(region.intersection(level_geometry(rel.xs, rel.ys, level, plan.min_thickness - 5e-5)))
        plan.regions.append(region if region.area >= 1e-3 else Polygon())
    while plan.regions and plan.regions[-1].is_empty:
        plan.regions.pop()
        plan.sections.pop()
        st.count -= 1
        st.c = st.c[:-1]
    if not plan.regions:
        raise NotSliceableError("aucune section : le maillage est ouvert ou dégénéré")


def layer_paths(plan: GeneralPlan, params, speed_first: float, speed: float):
    """Parois et remplissage de chaque couche, relevés sur F_k. Rend ([(k, 0, vitesse, polylignes 3D)], [types de chaque polyligne par couche]).

    Peau pleine : la part de la région absente d'une des `nb` couches du dessous ou d'une des `nt` du dessus. Présence
    lue sur les sections brutes : une couche écrasée sous une surface épousée n'est pas imprimée, mais sa matière l'est
    (par la couche suivante) ; la compter absente rendait pleines toutes les couches sous chaque peau (+13 % de
    filament sur un prisme en pente)."""
    tpar: tp.ToolpathParams = params.toolpath
    st, rel = plan.stack, plan.relief
    h = plan.thickness
    nb = max(1, math.ceil(tpar.solid_bottom_mm / h))
    nt = max(1, math.ceil(tpar.solid_top_mm / h))
    regs = plan.regions
    total = len(regs)
    out, kinds = [], []
    for i, region in enumerate(regs):
        _progress.report(i, total)
        if region.is_empty:
            continue
        below = _intersection(plan.sections, range(i - nb, i))
        above = _intersection(plan.sections, range(i + 1, i + 1 + nt))
        core = below.intersection(above) if not (below.is_empty or above.is_empty) else Polygon()
        solid = region.difference(core) if not core.is_empty else region
        if not core.is_empty:
            solid = _polygons_only(solid.buffer(-0.1).buffer(0.1))
        path = tp.layer_polylines(region, solid, i + 1, tpar, region, None)
        xy_list = [tp.densify_xy(pl.xy, tpar.max_segment) for pl in path.polylines]
        polys3d, layer_kinds = [], []
        for pl, xy in zip(path.polylines, xy_list):
            if len(xy) < 2:
                continue
            z = st.z_xy(i + 1, xy[:, 0], xy[:, 1])
            polys3d.append(np.column_stack([xy, z]))
            layer_kinds.append(pl.kind)
        out.append((i + 1, 0, speed_first if i == 0 else speed, polys3d))
        kinds.append(layer_kinds)
    return out, kinds


def _intersection(regs, idx):
    acc = None
    for j in idx:
        if j < 0 or j >= len(regs) or regs[j].is_empty:
            return Polygon()
        acc = regs[j] if acc is None else acc.intersection(regs[j])
        if acc.is_empty:
            return Polygon()
    return acc if acc is not None else Polygon()


# --------------------------------------------------------------------------------------- contrôles

def surfaces_violation(stack, field_: hf.Heightfield, angle_deg: float, reach: float, *, cell: float = 1.0) -> float:
    """Pire dépassement (mm) de la condition « la matière voisine ne dépasse pas la pointe de plus que φ » sur les
    surfaces du moteur « surface en hauteur ». 0 ou moins : la buse passe partout. Sert à décider si ce moteur respecte
    les limites de la buse ou s'il faut passer au tranchage général.

    Pour une couche, p parcourt les nœuds où elle dépose et q les nœuds déjà imprimés ; on compare la hauteur déjà posée
    en q à la pointe en p, plus φ(|p − q|). Calcul par transformée en chanfrein (une par couche), en partant des couches
    du haut, les plus penchées, pour sortir au premier dépassement."""
    stride = max(1, int(round(cell / (field_.xs[1] - field_.xs[0]))))
    xs, ys = field_.xs[::stride], field_.ys[::stride]
    gx, gy = np.meshgrid(xs, ys)
    col = stack.sample(gx.ravel(), gy.ravel())
    valid = ~np.isnan(field_.top[::stride, ::stride])
    last = stack.last_grid()[::stride, ::stride]
    split = stack.split_grid()[::stride, ::stride]
    slope = math.tan(math.radians(min(max(angle_deg, 0.0), 89.0))) * _CHAMFER_SAFETY
    step = slope * (xs[1] - xs[0])
    worst = -np.inf
    for k, j in reversed(stack.sequence()):
        active = stack.active_mask(k, j, last, split) & valid
        if not active.any():
            continue
        tip = col.layer_z(k, j).reshape(gx.shape)
        printed = tip - col.layer_thickness(k, j).reshape(gx.shape)
        src = np.where(active, tip, np.inf)
        bound = np.minimum(cone_minorant(src, step), float(src.min()) + reach)
        worst = max(worst, float(np.nanmax(np.where(valid, printed - bound, -np.inf))))
        if worst > 0.0:
            break
    return worst


def head_check(plan: GeneralPlan, layers_data, env: hc.HeadEnvelope, margin: float, report: hc.HeadReport,
               stride: int = 6) -> None:
    """Contrôle de la tête contre la hauteur DÉJÀ imprimée, entretenue couche après couche sur la grille du relief."""
    rel = plan.relief
    gx, gy = np.meshgrid(rel.xs, rel.ys)
    printed_h = np.full(gx.shape, plan.z_bottom)
    grid = hf.Heightfield(rel.xs, rel.ys, printed_h, printed_h, None)   # type: ignore[arg-type]
    total = sum(len(pts[::stride]) for _, _, _, polys in layers_data for pts in polys) or 1
    done = 0
    for k, j, _speed, polys in layers_data:
        region = plan.regions[k - 1]
        if polys:
            samp = np.vstack([pts[::stride] for pts in polys])
            snapshot = printed_h.copy()

            def printed(x, y, snap=snapshot):
                return grid._interp(snap, x, y)

            hit, pen = hc.check_points(env, plan.field, None, None, samp[:, :2], samp[:, 2], margin, printed_fn=printed)
            report.points += len(samp)
            report.colliding_points += hit
            report.max_penetration = max(report.max_penetration, pen)
            done += len(samp)
            _progress.report(done, total)
        if not region.is_empty:
            shapely.prepare(region)
            inside = shapely.contains_xy(region, gx.ravel(), gy.ravel()).reshape(gx.shape)
            zk = plan.stack.z_grid(k)
            printed_h[:] = np.where(inside, np.maximum(printed_h, zk), printed_h)
