"""Relief LOCAL du tranchage général d'Element Aeolis : chaque surface tournée vers le haut est épousée par une PEAU.

La couche k est la section du maillage par z = F_k(x, y) = c_k + δ_k(x, y) (une grille δ par couche).

Construction (révision du 2026-11-15, « peau et écrasement »).
  1. CIBLES. Les faces tournées vers le haut que la buse peut suivre (pente ≤ min(angle de buse, `_FOLLOW_ANGLE_DEG`))
     forment des surfaces connexes ; une surface entière est une bande (coupée seulement au-delà de la longueur de
     buse). Sa couche de tête E est le plan le plus proche de son sommet ; la peau (`_skin`, `_SKIN_LAYERS` couches et
     leur fond) lui est parallèle (`_skin_targets`). Ignorées : les surfaces plus étroites que `_MIN_WIDTH`, et le bord
     d'une surface au pied d'une paroi qui la domine (`_clear_of_walls`).
  2. LIMITES. Dans chaque couche, pente ≤ L = tan(min(angle de buse, `_RAMP_ANGLE_DEG`)) ; en chaque point, épaisseur
     dans [0, (1 + a)·h] : une couche peut S'ÉCRASER (épaisseur nulle) mais pas s'épaissir de plus de 25 % ; première
     couche plate ; étendue par couche ≤ longueur de buse. δ est l'extension lipschitzienne 3D des cibles sous ces
     limites (`_bounds`, séparable) ; δ = 0 (le plan) partout où c'est permis (`_within_nozzle`). Sous une peau, les
     couches s'écrasent donc sur son fond ; ailleurs elles restent planes. La matière déjà posée ne dépasse jamais la
     pointe de plus que min(r·tan α, longueur) : la buse ne touche rien, par construction.
  3. CONFLITS. Cibles acceptées de la plus haute à la plus basse, surfaces en pente avant surfaces plates
     (`_accept_from_the_top`).

Historique. Relief unique calé sur le sommet (jusqu'au 2026-11-14) : marches sur toutes les surfaces plus basses.
Relief local à épaisseur h·(1 ± 25 %) (2026-11-14/15) : chaque décalage de couche ne s'éteignait qu'en 4 fois sa
hauteur, à travers les parois, d'où des lignes tordues sur les flancs, la face arrière et autour des hublots d'un
Benchy ; le couper en bandes laissait un raccord visible par bande. Voir `docs/limites/087.md`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import shapely
import trimesh
from scipy import ndimage

from . import heightfield as hf

_CHAMFER_SAFETY = 0.93       # cos(22,5°) : la transformée en chanfrein tolère jusqu'à 1/0,924 fois la pente
_AMP_FRACTION = 0.25         # taux a : épaisseur au plus (1 + a)·h
_SKIN_LAYERS = 3             # couches parallèles à une surface épousée (peau) ; dessous, elles s'écrasent
_MIN_PATCH_AREA = 1.0        # surface minimale (mm², en projection) d'une surface épousée
_FINE = 4                    # finesse des cibles et de leurs cônes : cellule du relief / 4 (voir `_bounds`)
_ACCEPT_SLICE = 0.5          # tranche de hauteur de surface (mm) pour accepter les cibles de haut en bas
_SMALL_AREA = 10.0           # aire projetée (mm²) d'un petit détail, à peau plus épaisse (`_skin_depth`)
_EXTEND = 0.5                # largeur (mm) du prolongement des cibles dans le vide voisin (`_extensions`)
_MIN_WIDTH = 0.6             # largeur minimale (mm, diamètre inscrit en projection) d'une surface épousée
_RAMP_ANGLE_DEG = 70.0       # pente maximale des raccords (plafonnée par l'angle de la buse) : voir `build_relief`
_FOLLOW_ANGLE_DEG = 45.0     # pente maximale d'une surface épousée : au-delà, une marche (h / tan pente) est plus
                             # étroite que le cordon, donc invisible, et la couche penchée n'apporte plus rien.
                             # Les faces plus raides SÉPARENT les surfaces : le pont et le plat-bord d'une coque ne
                             # forment pas une seule surface reliée par le pavois.


def cone_minorant(f: np.ndarray, step: float) -> np.ndarray:
    """min_q f(q) + step·|p − q|/cellule, par transformée en chanfrein (poids 1 et √2 × `step`).

    `f` vaut +inf aux nœuds sans source. Deux balayages ; dans chaque ligne, la propagation gauche-droite
    puis droite-gauche est faite d'un coup par cumul de minimum."""
    out = np.array(f, dtype=float)
    ny, nx = out.shape
    wd = step * math.sqrt(2.0)
    ramp = np.arange(nx) * step
    for rows in (range(ny), range(ny - 1, -1, -1)):
        prev = None
        for j in rows:
            row = out[j]
            if prev is not None:
                row = np.minimum(row, prev + step)
                row[1:] = np.minimum(row[1:], prev[:-1] + wd)
                row[:-1] = np.minimum(row[:-1], prev[1:] + wd)
            row = np.minimum.accumulate(row - ramp) + ramp
            row = np.minimum.accumulate((row + ramp)[::-1])[::-1] - ramp
            out[j] = row
            prev = row
    return out


def _cone(f: np.ndarray, step: float, span: float) -> np.ndarray:
    """`cone_minorant` restreint à la boîte des sources élargie de `span / step` cellules, +inf ailleurs. EXACT pour
    l'usage qui en est fait : δ reste dans un intervalle de largeur `span`, donc un cône plus loin de ses sources que
    span / pente dépasse toute valeur possible de δ et ne contraint rien."""
    rows = np.flatnonzero(np.isfinite(f).any(axis=1))
    out = np.full(f.shape, np.inf, dtype=np.float32)
    if not rows.size:
        return out
    cols = np.flatnonzero(np.isfinite(f).any(axis=0))
    m = int(math.ceil(span / max(step, 1e-12))) + 1
    win = (slice(max(rows[0] - m, 0), rows[-1] + m + 1), slice(max(cols[0] - m, 0), cols[-1] + m + 1))
    out[win] = cone_minorant(f[win], step)
    return out


def amplitude(height: float, angle_deg: float, reach: float) -> float:
    """Amplitude du relief : la longueur de buse. 0 : planaire. (Jusqu'au 2026-11-15, aussi bornée au quart de la
    hauteur de la pièce : les couches ne pouvant varier que de 25 %, on ne pouvait creuser plus. Elles s'écrasent
    désormais sous la peau ; la profondeur est bornée par la première couche, voir `_band_targets`.)"""
    amp = min(reach, height)
    return amp if amp >= 1e-3 and angle_deg > 0.0 else 0.0


@dataclass
class Relief:
    xs: np.ndarray
    ys: np.ndarray
    d: np.ndarray                 # (K, ny, nx) : δ_k aux nœuds, dans [−amp, 0]
    amp: float
    max_gradient: float           # max |∇δ_k| (pente), sans unité
    patches: int = 0              # nombre de bandes épousées

    def _cells(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        fx = np.clip((x - self.xs[0]) / (self.xs[1] - self.xs[0]), 0, len(self.xs) - 1 - 1e-9)
        fy = np.clip((y - self.ys[0]) / (self.ys[1] - self.ys[0]), 0, len(self.ys) - 1 - 1e-9)
        ix, iy = fx.astype(int), fy.astype(int)
        return ix, iy, fx - ix, fy - iy

    def locate(self, x, y) -> tuple:
        """Cellules et poids bilinéaires des points : à calculer une fois, à réutiliser pour chaque couche."""
        return self._cells(x, y)

    def at(self, k: int, x=None, y=None, *, where: tuple | None = None) -> np.ndarray:
        """δ de la couche k (à partir de 1) aux points (x, y), ou aux points déjà localisés (`where`)."""
        ix, iy, tx, ty = where if where is not None else self._cells(x, y)
        g = self.d[min(max(k, 1), len(self.d)) - 1]
        return ((1 - tx) * (1 - ty) * g[iy, ix] + tx * (1 - ty) * g[iy, ix + 1]
                + (1 - tx) * ty * g[iy + 1, ix] + tx * ty * g[iy + 1, ix + 1])


def _patches(mesh: trimesh.Trimesh, min_nz: float) -> list[np.ndarray]:
    """Surfaces connexes de faces tournées vers le haut (n_z ≥ `min_nz`), en indices de faces ; les surfaces trop
    petites (`_MIN_PATCH_AREA` en projection) sont ignorées."""
    nz = mesh.face_normals[:, 2]
    keep = nz >= min_nz
    if not keep.any():
        return []
    adj = mesh.face_adjacency
    adj = adj[keep[adj[:, 0]] & keep[adj[:, 1]]]
    comps = trimesh.graph.connected_components(adj, nodes=np.flatnonzero(keep), min_len=1)
    proj = mesh.area_faces * nz
    return [np.asarray(c) for c in comps if proj[np.asarray(c)].sum() >= _MIN_PATCH_AREA]


@dataclass
class _Band:
    k: int                        # couche de tête : le plan le plus proche du sommet de la surface
    idx: np.ndarray               # pixels de la grille FINE (indices à plat)
    val: np.ndarray               # δ cible en ces pixels
    depth: int = _SKIN_LAYERS     # couches de peau (`_skin_depth`)
    ext: bool = False             # prolongement dans le vide autour d'une surface (`_extensions`)


def _skin(k: int, depth: int = _SKIN_LAYERS) -> range:
    """Surfaces de couche qui épousent une surface dont la couche de tête est k (indices à partir de 0) : k et les
    `_SKIN_LAYERS` du dessous, parallèles. Les `_SKIN_LAYERS` couches du haut ont ainsi l'épaisseur nominale ; la plus
    basse est le FOND de la peau, sur lequel s'écrasent les couches du dessous. (Sans ce fond, elles s'écrasaient sur
    la dernière couche de peau, qui tombait à une épaisseur nulle.) Jamais la première couche, qui reste plate.
    `depth` : nombre de couches de peau (`_skin_depth`)."""
    return range(max(k - depth, 1), k + 1) if k >= 1 else range(0, 1)


def _skin_depth(area: float, rise: float, h: float) -> int:
    """Couches de peau d'une surface d'aire projetée `area` et de dénivelé `rise`. Sous la peau, les couches planes
    s'arrêtent le long d'une ligne inclinée, visible sur la paroi qui porte la surface. Sur un PETIT détail (aire ≤
    `_SMALL_AREA`), la paroi est courte et ces arrêts, serrés, faisaient des bosses (cylindre du porte-canne d'un
    Benchy, anneau incliné de 6 mm², 1 mm de dénivelé) : la peau couvre alors tout le dénivelé, le détail est incliné
    sur toute sa hauteur et les arrêts tombent à son pied. Sur une grande surface, une peau aussi épaisse tordrait
    les parois sur toute la hauteur du dénivelé."""
    if area > _SMALL_AREA:
        return _SKIN_LAYERS
    return max(_SKIN_LAYERS, int(math.ceil(rise / h)) + 2)


def _vertical(dc: float) -> tuple[float, float]:
    """Écarts (au-dessus, au-dessous) tolérés entre le δ d'une couche et celui d'une autre plus basse de dc (mm ; dc < 0 :
    plus haute), au même point. L'épaisseur reste dans [0, (1 + a)·h] : en montant d'une couche, δ gagne au plus a·h
    et perd au plus h (couche écrasée)."""
    return (_AMP_FRACTION * dc, dc) if dc >= 0 else (-dc, -_AMP_FRACTION * dc)


def _band_targets(mesh: trimesh.Trimesh, patches: list, fine: tuple, planes: np.ndarray, amp: float,
                  slope: float) -> list:
    """Surfaces à épouser, rastérisées sur la grille FINE. Une surface entière est UNE bande (coupée seulement si son
    dénivelé dépasse l'amplitude) ; sa couche de tête E est le plan le plus proche de son sommet, donc δ = t − c_E
    ≤ h/2 : sous ses parties basses, les couches du dessous s'écrasent au lieu de se tordre. δ est borné pour que la
    peau (`_skin`) reste au-dessus de la première couche.

    Une surface plus étroite que `_MIN_WIDTH` n'est pas suivie : ses marches, plus courtes qu'un cordon, ne se voient
    pas, alors que la suivre tord la paroi qui la porte (dessus de l'encadrement d'un hublot de Benchy, 0,38 mm de
    large et 1,7 mm de dénivelé : lignes ondulées tout autour du hublot). L'anneau du porte-canne (0,8 mm) reste suivi.

    Le bord d'une surface qui touche une paroi plus haute n'est pas suivi (`_clear_of_walls`)."""
    xf, yf = fine
    sections: dict = {}
    cell = float(xf[1] - xf[0]) * float(yf[1] - yf[0])
    h = float(planes[1] - planes[0])
    rise = max(amp - h, h)
    bands = []
    for faces in patches:
        t = hf._rasterize(mesh.triangles[faces], xf, yf, keep_max=True)
        rest = ~np.isnan(t)
        if 2.0 * float(ndimage.distance_transform_edt(rest).max()) * float(xf[1] - xf[0]) < _MIN_WIDTH:
            continue
        while rest.any():
            top = float(t[rest].max())
            band = rest & (t > top - rise - 1e-9)
            rest &= ~band
            if band.sum() * cell < _MIN_PATCH_AREA and bands:
                continue
            k = int(np.argmin(np.abs(planes - top)))
            if k >= 1:
                band = _clear_of_walls(mesh, band, t, (float(planes[k - 1]), slope), fine, sections)
                if not band.any():
                    continue
            depth = _skin_depth(float(band.sum()) * cell, float(t[band].max() - t[band].min()), h)
            low = float(planes[_skin(k, depth)[0]])
            lo_d, up_d = float(planes[0]) - low, _AMP_FRACTION * (low - float(planes[0]))
            idx = np.flatnonzero(band)
            bands.append(_Band(k, idx, np.clip(t.ravel()[idx] - planes[k], lo_d, up_d).astype(np.float32), depth))
    return bands


@dataclass
class _Accepted:
    tgt: dict                     # couche -> cibles acceptées (grille FINE à plat, +inf hors empreinte)
    up_c: dict                    # couche -> cône majorant des cibles
    lo_c: dict                    # couche -> cône minorant
    planes: np.ndarray
    fine_shape: tuple
    slope: float
    span: float
    min_thickness: float          # plancher d'épaisseur du profil (voir `_skin_targets`)


def _clear_of_walls(mesh: trimesh.Trimesh, band: np.ndarray, t: np.ndarray, level: tuple, fine: tuple,
                    sections: dict) -> np.ndarray:
    """Pixels de la bande assez loin des parois qui la dominent. `level` = (c_{E−1}, pente des raccords).

    Une paroi qui monte au-dessus d'une surface épousée (cheminée sur un toit, bastingage au bord d'un pont) est déjà
    imprimée jusqu'à c_{E−1} quand la peau de la surface est posée en E. Pour que la buse passe, il faut
    c_{E−1} ≤ t(p) + pente·distance. Si la surface est épousée jusqu'au pied de la paroi, la paroi doit s'enfoncer, et
    ses couches ne regagnent leur plan qu'à 25 % d'épaisseur en plus par couche : lignes ondulées sur 4 fois
    l'enfoncement (pied de la cheminée d'un Benchy, face intérieure du bastingage). On renonce donc à épouser la
    surface à moins de (c_{E−1} − t)/pente de la paroi : moins de 0,5 mm autour de la cheminée. Paroi : matière à la
    hauteur c_{E−1} (section plane du maillage), hors de l'empreinte de la bande elle-même."""
    plane, slope = level
    need = np.where(band, (plane - t) / slope, 0.0)
    if not (need > 0.0).any():
        return band
    xf, yf = fine
    cell = float(xf[1] - xf[0])
    r = int(math.ceil(float(need.max()) / cell)) + 1
    iy, ix = np.nonzero(band)
    win = (slice(max(int(iy.min()) - r, 0), int(iy.max()) + r + 1), slice(max(int(ix.min()) - r, 0), int(ix.max()) + r + 1))
    poly = _solid_at(mesh, plane, sections)
    if poly.is_empty:
        return band
    gx, gy = np.meshgrid(xf[win[1]], yf[win[0]])
    solid = shapely.contains_xy(poly, gx, gy) & ~band[win]
    if not solid.any():
        return band
    dist = ndimage.distance_transform_edt(~solid) * cell
    out = band.copy()
    out[win] &= dist >= need[win]
    return out


def _skin_targets(k: int, val: np.ndarray, planes: np.ndarray, min_thickness: float, depth: int) -> list:
    """δ cibles des couches de la peau (`_skin`) d'une surface de couche de tête k, aux pixels de valeurs `val`.

    La peau a m couches au-dessus de son fond. Fond nominal : t − m·h. Sous lui, les couches s'écrasent sur le fond ;
    la plus haute qui reste a l'épaisseur fond − c_p (c_p : plan juste dessous). Si elle est plus mince que le plancher
    `min_thickness`, on descend le fond sur c_p et les m couches de peau se partagent ce reste (au plus plancher / m de
    plus chacune, sous le plafond de 25 %). Sans cela, la tranche trop mince était reportée sur la couche suivante,
    jusqu'à h + plancher : au-dessus du plafond du profil."""
    layers = list(_skin(k, depth))
    if len(layers) == 1:
        return [(k, val)]
    h = float(planes[1] - planes[0])
    m = len(layers) - 1
    top = planes[k] + val
    floor = top - m * h
    below = planes[np.clip(np.searchsorted(planes, floor + 1e-6, side="right") - 1, 0, len(planes) - 1)]
    floor = np.where((floor - below < min_thickness) & (below >= planes[0]), below, floor)
    floor = np.maximum(floor, planes[0])
    step = (top - floor) / m
    return [(j, (top - (k - j) * step - planes[j]).astype(np.float32)) for j in layers]


def _solid_at(mesh: trimesh.Trimesh, plane: float, sections: dict):
    """Section plane du maillage à la hauteur `plane` (matière à cette hauteur), mise en cache dans `sections`."""
    if plane not in sections:
        from .general import section_region      # import tardif : general importe ce module
        sections[plane] = section_region(mesh.vertices, mesh.faces, mesh.vertices[:, 2] - plane - 1e-6)
    return sections[plane]


def _extensions(mesh: trimesh.Trimesh, bands: list, fine: tuple, planes: np.ndarray) -> list:
    """Prolongements des surfaces épousées dans le VIDE voisin, sur `_EXTEND`, à la valeur du pixel le plus proche.

    Rien n'y est déposé, mais la hauteur d'un point de trajet ou de coupe est interpolée sur la grille du relief (0,25
    mm) : au bord d'une surface étroite, elle mêlait la surface et la couche voisine restée plane dans le vide. Cas
    réel : l'anneau incliné du porte-canne d'un Benchy (0,8 mm, trois nœuds) bordé par son trou, où la couche reste à
    son plan ; sur tout le côté bas de l'anneau, la coupe passait au-dessus de lui et la couche de tête disparaissait
    (dessus bosselé). Vide : pas de matière à c_{E−1}, et aucune autre surface. Acceptés en dernier."""
    xf, yf = fine
    shape = (len(yf), len(xf))
    occupied = np.zeros(shape[0] * shape[1], dtype=bool)
    for b in bands:
        occupied[b.idx] = True
    occupied = occupied.reshape(shape)
    r = int(round(_EXTEND / float(xf[1] - xf[0])))
    sections: dict = {}
    out = []
    for b in bands:
        if b.k < 1 or not b.idx.size:
            continue
        iy, ix = np.divmod(b.idx, shape[1])
        y0, x0 = max(int(iy.min()) - r, 0), max(int(ix.min()) - r, 0)
        win = (slice(y0, int(iy.max()) + r + 1), slice(x0, int(ix.max()) + r + 1))
        vals = np.full(occupied[win].shape, np.nan, dtype=np.float32)
        vals[iy - y0, ix - x0] = b.val
        dist, (ny, nx) = ndimage.distance_transform_edt(np.isnan(vals), return_indices=True)
        ring = (dist > 0) & (dist <= r) & ~occupied[win]
        ry, rx = np.nonzero(ring)
        if not ry.size:
            continue
        poly = _solid_at(mesh, float(planes[b.k - 1]), sections)
        air = ~shapely.contains_xy(poly, xf[rx + x0], yf[ry + y0]) if not poly.is_empty else np.ones(len(ry), bool)
        ry, rx = ry[air], rx[air]
        if ry.size:
            out.append(_Band(b.k, (ry + y0) * shape[1] + rx + x0, vals[ny[ry, rx], nx[ry, rx]], b.depth, ext=True))
    return out


def _accept_pass(bands: list, acc: _Accepted, n: int) -> None:
    """Une passe d'acceptation de haut en bas (voir `_accept_from_the_top`) sur les bandes `bands`."""
    planes = acc.planes
    z = [planes[b.k] + b.val for b in bands]
    if not z:
        return
    level = max(float(zz.max()) for zz in z)
    low = min(float(zz.min()) for zz in z)
    while level >= low - _ACCEPT_SLICE:
        changed = set()
        for b, zz in zip(bands, z):
            sl = np.flatnonzero((zz <= level + 1e-9) & (zz > level - _ACCEPT_SLICE))
            if not sl.size:
                continue
            idx, val = b.idx[sl], b.val[sl]
            ok = np.ones(len(idx), bool)
            for j in acc.up_c:
                ru, rl = _vertical(float(planes[b.k] - planes[j]))
                ok &= (val <= acc.up_c[j][idx] + ru + 1e-4) & (val >= acc.lo_c[j][idx] - rl - 1e-4)
            if ok.any():
                for k, v in _skin_targets(b.k, val[ok], planes, acc.min_thickness, b.depth):
                    cur = acc.tgt.setdefault(k, np.full(n, np.inf, dtype=np.float32))
                    np.minimum.at(cur, idx[ok], v)
                    changed.add(k)
        for k in changed:
            t = acc.tgt[k].reshape(acc.fine_shape)
            acc.up_c[k] = _cone(t, acc.slope, acc.span).ravel()
            acc.lo_c[k] = -_cone(np.where(np.isinf(t), np.inf, -t), acc.slope, acc.span).ravel()
        level -= _ACCEPT_SLICE


def _accept_from_the_top(bands: list, planes: np.ndarray, fine_shape: tuple, slope: float, span: float,
                         min_thickness: float = 0.0) -> dict:
    """Cibles de chaque couche (grille FINE, NaN hors des empreintes), choisies de la surface la plus HAUTE à la plus
    basse, par tranches de hauteur de surface (`_ACCEPT_SLICE`), pixel par pixel. Une cible vaut pour toute la peau
    de sa surface (`_skin`).

    Deux cibles sont incompatibles si l'écart de leurs δ dépasse pente·|p − q| plus l'écart vertical toléré
    (`_vertical`) ; l'extension céderait alors en abaissant la plus haute, la surface visible. Un pixel n'est donc
    gardé que s'il est entre les bornes 3D des tranches déjà acceptées (plus hautes). Cas réel : le fond du trou d'un
    porte-canne, sous et à côté de l'anneau incliné qui le borde, écrasait l'anneau en marches. L'ordre se fait par
    pixel et non par bande : une bande dont le sommet est haut peut passer, ailleurs, sous une surface plus haute. Les
    bornes ne sont lues qu'aux pixels de la tranche ; les cônes ne sont recalculés que pour les couches modifiées.

    Deux passes : les surfaces en pente (dénivelé > h) d'abord, les plates ensuite. Une surface plate n'a presque rien
    à gagner à être épousée (au pire une demi-couche d'écart), une surface en pente évite une marche par couche : sur
    un Benchy, l'appui de la fenêtre avant, plat et plus haut, refusait 31 % du pont avant, alors en marches."""
    n = fine_shape[0] * fine_shape[1]
    acc = _Accepted({}, {}, {}, planes, fine_shape, slope, span, min_thickness)
    h = float(planes[1] - planes[0])
    real = [b for b in bands if b.val.size and not b.ext]
    sloped = [b for b in real if float(b.val.max() - b.val.min()) > h]
    flat = [b for b in real if float(b.val.max() - b.val.min()) <= h]
    for group in (sloped, flat, [b for b in bands if b.ext]):
        _accept_pass(group, acc, n)
    tgt = acc.tgt
    return {k: np.where(np.isinf(v), np.nan, v).reshape(fine_shape) for k, v in tgt.items()}


def _bounds(fine_tgt: dict, shape: tuple, slope: float, h: float, span: float) -> tuple[np.ndarray, np.ndarray]:
    """Bornes U (haute) et L (basse) de l'extension lipschitzienne 3D des cibles, sur la pile du relief `shape` :
    U(p, k) = min sur les cibles (q, j) de t(q, j) + pente·|p − q| + écart vertical toléré (`_vertical`), L de même
    avec max et −. La première couche est une cible δ = 0 partout. Séparable : cônes en chanfrein sur la grille FINE
    de chaque couche cible (`slope` par cellule fine), relus aux nœuds du relief, puis cumuls en hauteur, vers le haut
    puis vers le bas (taux a·h et h par couche). Seules les EMPREINTES contraignent.

    Pourquoi la grille fine : les cônes partent des pixels du bord des empreintes, et un bord oblique à la grille est en
    escalier (une rangée de pixels tous les cellule / pente du bord). Entre deux empreintes proches (pont et plat-bord
    d'un Benchy), les couches qui traversent le flanc n'ont aucune liberté : l'escalier y fixe leur hauteur, et la
    section du flanc évasé ondule (mesuré : p99 0,049 mm à 0,25 mm de cellule, 0,017 mm avec des cônes 4 fois plus
    fins). La grille fine n'existe qu'en 2D, pour les seules couches cibles : la pile reste à la cellule du relief."""
    up = np.full(shape, np.inf, dtype=np.float32)
    lo = np.full(shape, -np.inf, dtype=np.float32)
    for k, t in fine_tgt.items():
        on = ~np.isnan(t)
        up[k] = _cone(np.where(on, t, np.inf), slope, span)[::_FINE, ::_FINE]
        lo[k] = -_cone(np.where(on, -t, np.inf), slope, span)[::_FINE, ::_FINE]
    up[0] = np.minimum(up[0], 0.0)
    lo[0] = np.maximum(lo[0], 0.0)
    k = np.arange(shape[0], dtype=np.float32)[:, None, None]
    rise, drop = k * np.float32(_AMP_FRACTION * h), k * np.float32(h)
    up = np.minimum.accumulate(up - rise, axis=0) + rise
    up = np.minimum.accumulate((up + drop)[::-1], axis=0)[::-1] - drop
    lo = np.maximum.accumulate(lo + drop, axis=0) - drop
    lo = np.maximum.accumulate((lo - rise)[::-1], axis=0)[::-1] + rise
    return lo, up


def _within_nozzle(fine_tgt: dict, planes: np.ndarray, shape: tuple, slope: float, cell: float, span: float) -> np.ndarray:
    """Champ δ : 0 (le plan) ramené entre les bornes de la buse. Min et max de champs qui respectent la pente et les
    écarts verticaux : le résultat les respecte, vaut la cible sur les empreintes compatibles et 0 partout où c'est
    permis. Sous une surface, la peau la suit et les couches du dessous s'écrasent sur la peau (épaisseur nulle) au
    lieu de se tordre ; à côté, elles reprennent leur plan à la pente `slope`. `cell` : cellule FINE."""
    h = float(planes[1] - planes[0]) if len(planes) > 1 else 0.0
    lo, up = _bounds(fine_tgt, shape, slope * cell, h, span)
    return np.minimum(np.maximum(0.0, lo), up).astype(np.float32)


def build_relief(mesh: trimesh.Trimesh, grid: tuple, planes: np.ndarray, angle_deg: float, reach: float,
                 min_thickness: float = 0.0) -> Relief:
    """Champ δ des couches `planes` (pas uniforme) sur la grille `grid` = (xs, ys). Voir l'en-tête du module."""
    xs, ys = grid
    shape = (len(planes), len(ys), len(xs))
    height = float(mesh.bounds[1, 2] - mesh.bounds[0, 2])
    amp = amplitude(height, angle_deg, reach)
    flat = Relief(xs, ys, np.zeros(shape, dtype=np.float32), 0.0, 0.0)
    if amp == 0.0:
        return flat
    # Surfaces suivies : pente ≤ min(buse, 45°). Raccords : pente ≤ min(buse, 70°). Un raccord raide reste COURT : à
    # 45°, le bastingage d'un Benchy, à côté d'un pont épousé 5 mm plus bas, gardait ses couches enfoncées sur toute
    # sa hauteur, en lignes penchées terminées par un pli sur le flanc ; à 70°, le flanc est plat. À 79°, des stries
    # apparaissaient au pied de la cheminée. Les raccords raides tombent surtout dans le vide ou dans les couches
    # écrasées ; la buse les tolère (contrôle de la tête : 0 collision).
    limit = math.tan(math.radians(min(angle_deg, _FOLLOW_ANGLE_DEG)))
    patches = _patches(mesh, 1.0 / math.hypot(1.0, limit))
    if not patches:
        return flat
    fine = (np.linspace(xs[0], xs[-1], (len(xs) - 1) * _FINE + 1),
            np.linspace(ys[0], ys[-1], (len(ys) - 1) * _FINE + 1))
    fine_shape = (len(fine[1]), len(fine[0]))
    bands = _band_targets(mesh, patches, fine, planes, amp, math.tan(math.radians(min(angle_deg, _RAMP_ANGLE_DEG))))
    count = len(bands)
    bands += _extensions(mesh, bands, fine, planes)
    slope = math.tan(math.radians(min(angle_deg, _RAMP_ANGLE_DEG))) * _CHAMFER_SAFETY
    cell = float(xs[1] - xs[0]) / _FINE
    span = amp + float(planes[1] - planes[0])           # δ ∈ [−amp, h/2]
    fine_tgt = _accept_from_the_top(bands, planes, fine_shape, slope * cell, span, min_thickness)
    d = _within_nozzle(fine_tgt, planes, shape, slope, cell, span)
    return Relief(xs, ys, d, amp, _max_slope(d, xs, ys), count)


def _max_slope(d: np.ndarray, xs: np.ndarray, ys: np.ndarray) -> float:
    """Plus forte pente entre deux nœuds voisins (côtés et diagonales) d'une même couche. Une différence centrée
    mélangerait, à une arête du relief, les pentes de deux versants et surestimerait jusqu'à √2."""
    dx, dy = float(xs[1] - xs[0]), float(ys[1] - ys[0])
    diag = math.hypot(dx, dy)
    return float(max(np.abs(np.diff(d, axis=2)).max() / dx, np.abs(np.diff(d, axis=1)).max() / dy,
                     np.abs(d[:, 1:, 1:] - d[:, :-1, :-1]).max() / diag,
                     np.abs(d[:, 1:, :-1] - d[:, :-1, 1:]).max() / diag))
