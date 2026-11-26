"""
element_ardor/step4_planar.py

Couches PLANES partout ou l'impression plane suffit (etape 4, apres la methode de la chaleur).

Constat (Benchy droit, 2026-11-18) : la distance au plateau a travers la matiere courbe les couches PARTOUT ou la
matiere ne monte pas a la verticale, meme la ou un trancheur classique imprime sans peine. La premiere couche
allait de z = 0 a 0,55 mm (fond tranche en escalier, motif en labyrinthe sous la coque, premiere couche qui
n'adhere pas) ; les dessus plats (pont, toit, bastingage) etaient coupes en biais par des dizaines de couches,
chacune laissant ses parois en boucles serpentines sur la surface. Or le Benchy est concu pour s'imprimer a plat.

Principe : la direction de croissance X d'un tetraedre est un melange de la verticale et de la direction de la
methode de la chaleur, `X = normaliser((1 - w) * haut + w * d_chaleur)`, puis le champ est reintegre
(`L phi = div X`, ancre a 0 sur le plateau). Le poids w vaut :
  - 0 (couches planes, trancheur classique) loin de tout surplomb que l'impression plane ne sait pas faire ;
  - 1 (champ de la chaleur, inchange) a moins de `core` mm d'un tel surplomb, avec une rampe douce (`ramp` mm)
    entre les deux : les couches tournent progressivement, sans ligne de separation entre zone plane et zone
    inclinee (l'influence d'un surplomb se propage, elle n'est pas coupee) ;
  - 0 sous `flat_height` mm (premieres couches planes et pleines sur le plateau), rampe comprise.

Surplomb « a incliner » (critere d'un trancheur classique, `tilt_points`) : pente au-dela de `overhang_deg` depuis
la verticale ; plafond qui n'est ni une levre ni un pont (un plafond porte des deux cotes s'imprime en pont, aucune
inclinaison ne l'aide) ; en morceaux d'au moins `min_patch_area` mm2 (un trou de 2 mm ou une gravure s'impriment a
plat avec un leger affaissement, comme dans un trancheur classique ; incliner pour eux deformerait toutes les couches
voisines).

Desactive (None) dans `FieldSolverParams.planar`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from core.mesh_types import MeshData, TetMesh
from element_ardor.step4_conform import ConformParams


@dataclass
class PlanarParams:
    overhang_deg: float = 60.0       # surplomb en pente imprimable a plat (depuis la verticale) ; au-dela, on incline
    ceiling_deg: float = 80.0        # au-dela : plafond (porte dans toutes les directions : levre ou pont)
    lip: float = 1.5                 # mm : un surplomb a moins de cela d'un appui s'imprime a plat (levre, rebord)
    slope_reach: float = 2.0         # mm : une pente portee par un appui a moins de cela, dans le sens de la descente
    bridge_span: float = 25.0        # mm : portee d'un pont imprime a plat entre deux appuis (Benchy : cabine, 18 mm)
    min_patch_area: float = 4.0      # mm2 : un surplomb plus petit s'imprime a plat
    core: float = 2.0                # mm : distance au surplomb sous laquelle le champ de la chaleur est garde tel quel
    ramp: float = 8.0                # mm : largeur de la transition chaleur -> plan
    flat_height: float = 0.6         # mm au-dessus du plateau : couches toujours planes (plus une rampe de meme hauteur)
    smooth_iterations: int = 20      # lissages de la direction (moyenne aux noeuds puis aux tetraedres) dans la zone inclinee
    # couches conformes sous les grandes surfaces en pente (`element_ardor/step4_conform.py`, 2026-11-25, limite #92) ;
    # None = aucune (couches planes jusqu'au dessus : toit du Benchy en escalier)
    conform: Optional[ConformParams] = field(default_factory=ConformParams)


@dataclass(frozen=True)
class FieldInputs:
    """Ce dont le melange a besoin : maillages, gradients de forme, volumes, noeuds du plateau, faces en contact."""
    tet_mesh: TetMesh
    mesh_data: MeshData
    shape_gradients: Optional[np.ndarray]
    volumes: Optional[np.ndarray]
    platform_nodes: np.ndarray
    platform_faces: Optional[np.ndarray]


def _smoothstep(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def _area(v: np.ndarray, faces: np.ndarray) -> np.ndarray:
    return 0.5 * np.linalg.norm(np.cross(v[faces[:, 1]] - v[faces[:, 0]], v[faces[:, 2]] - v[faces[:, 0]]), axis=1)


def _overhang_faces(mesh_data: MeshData, platform_faces, params: PlanarParams, up: np.ndarray):
    """(faces en surplomb au-dela de `overhang_deg`, sinus de l'angle de surplomb par face). Ignore le contact
    plateau, les faces a moins de 2 x `flat_height` du plateau (premieres couches toujours planes) et les morceaux
    connexes de moins de `min_patch_area` mm2."""
    faces = np.asarray(mesh_data.faces)
    v = np.asarray(mesh_data.vertices, dtype=float)
    n = np.asarray(mesh_data.face_normals, dtype=float)
    n = n / np.maximum(np.linalg.norm(n, axis=1), 1e-12)[:, None]
    down = -(n @ up)
    over = down > math.sin(math.radians(params.overhang_deg))
    over &= v[faces].mean(axis=1) @ up > float((v @ up).min()) + 2.0 * params.flat_height
    if platform_faces is not None and len(platform_faces) == len(over):
        over &= ~np.asarray(platform_faces, dtype=bool)
    idx = np.flatnonzero(over)
    if idx.size:
        rows = np.repeat(np.arange(len(idx)), 3)
        inc = sp.csr_matrix((np.ones(len(rows)), (rows, faces[idx].ravel())), shape=(len(idx), len(v)))
        count, label = connected_components(inc @ inc.T, directed=False)
        patch_area = np.bincount(label, weights=_area(v, faces[idx]), minlength=count)
        over[idx[patch_area[label] < params.min_patch_area]] = False
    return over, down


def _sample_faces(v: np.ndarray, faces: np.ndarray, spacing: float) -> tuple:
    """Points repartis sur chaque face (grille barycentrique, pas <= `spacing`) : (points, indice de leur face)."""
    tri = v[faces]
    longest = np.max(np.linalg.norm(tri - np.roll(tri, 1, axis=1), axis=2), axis=1)
    div = np.clip(np.ceil(longest / spacing).astype(int), 1, 64)
    pts, owner = [], []
    for d in np.unique(div):
        k = np.flatnonzero(div == d)
        i, j = np.meshgrid(np.arange(d + 1), np.arange(d + 1), indexing="ij")
        ok = i + j <= d
        a, b = i[ok] / d, j[ok] / d
        w = np.stack([1.0 - a - b, a, b], axis=1)
        pts.append(np.einsum("sk,fkd->fsd", w, tri[k]).reshape(-1, 3))
        owner.append(np.repeat(k, len(w)))
    return np.concatenate(pts), np.concatenate(owner)


class _Solid:
    """Test « dans la matiere » sur le maillage tetraedrique (robuste aux defauts du STL : jonctions en T, trous)."""

    def __init__(self, tm: TetMesh):
        t = tm.tet_vertices[tm.tetrahedra]
        self.origin = t[:, 0]
        m = np.stack([t[:, 1] - t[:, 0], t[:, 2] - t[:, 0], t[:, 3] - t[:, 0]], axis=2)
        ok = np.abs(np.linalg.det(m)) > 1e-15
        self.inv = np.zeros_like(m)
        self.inv[ok] = np.linalg.inv(m[ok])
        self.tree = cKDTree(t.mean(axis=1))

    def contains(self, p: np.ndarray, k: int = 12) -> np.ndarray:
        if not len(p):
            return np.zeros(0, dtype=bool)
        _, idx = self.tree.query(p, k=k)
        bc = np.einsum("nkij,nkj->nki", self.inv[idx], p[:, None, :] - self.origin[idx])
        inside = (bc >= -1e-6).all(axis=2) & (bc.sum(axis=2) <= 1.0 + 1e-6)
        return inside.any(axis=1)


def _printable_surface(v: np.ndarray, faces: np.ndarray, over: np.ndarray, down: np.ndarray):
    """Distance exacte aux faces imprimables (parois et pentes sous `overhang_deg`, pas les dessus) : requete de l'etape
    5c (native si active). Exacte, pas par centres de faces : les grands triangles coplanaires d'une union (flanc en L
    d'un porte-a-faux) passaient pour une paroi sous tout le bras."""
    from element_ardor.step5c_skin_projection import _Surface
    rest = np.flatnonzero(~over & (down > -0.999))
    sub = faces[rest]
    used, remap = np.unique(sub, return_inverse=True)
    n = np.zeros((len(rest), 3))
    return _Surface(v[used], remap.reshape(-1, 3), n, 200_000)


def _supports(v: np.ndarray, faces: np.ndarray, up: np.ndarray, solid: _Solid, printable) -> np.ndarray:
    """Points d'appui de faces en surplomb : sur leurs aretes de bord, la ou une sonde placee juste au-dela du bord et
    un peu plus bas est dans la matiere (paroi qui descend, meme sans sommet partage : jonction en T) ou contre une face
    imprimable (`printable` : flanc d'une arche, sous lequel il n'y a que du vide). Le cote d'un porte-a-faux, ou la
    paroi monte, n'est pas un appui."""
    spacing = 0.5
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    third = np.concatenate([faces[:, 2], faces[:, 0], faces[:, 1]])
    _, inv, cnt = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
    border = cnt[inv.ravel()] == 1
    a, b, o = v[edges[border, 0]], v[edges[border, 1]], v[third[border]]
    if not len(a):
        return np.zeros((0, 3))
    m = np.maximum(1, np.ceil(np.linalg.norm(b - a, axis=1) / spacing).astype(int))
    e = np.repeat(np.arange(len(a)), m)
    t = (np.arange(m.sum()) - np.repeat(np.cumsum(m) - m, m) + 0.5) / np.repeat(m, m)
    q = a[e] + (b - a)[e] * t[:, None]
    out = np.cross((b - a)[e], up)
    out /= np.maximum(np.linalg.norm(out, axis=1), 1e-12)[:, None]
    out[np.einsum("nd,nd->n", out, o[e] - q) > 0.0] *= -1.0     # vers l'exterieur de la face (loin du 3e sommet)
    # deux sondes : presque sous le bord (paroi qui descend) et plus au-dela (flanc incline d'une arche)
    held = np.zeros(len(q), dtype=bool)
    for side, drop in ((0.1, 0.3), (0.3, 0.2)):
        _c, dist, _f = printable.closest(q + side * out - drop * up)
        held |= dist <= 0.15
    return q[held | solid.contains(q + 0.3 * out - 0.4 * up)]


def _anchored(p: np.ndarray, sup: np.ndarray, params: PlanarParams) -> np.ndarray:
    """Points portes : levre (appui a moins de `lip`) ou pont (appuis opposes, portee <= `bridge_span`). L'appui
    oppose est cherche dans toute la boule de rayon `bridge_span` : les appuis sont serres le long d'un bord, les k plus
    proches peuvent tous etre du meme cote."""
    if not len(sup):
        return np.zeros(len(p), dtype=bool)
    tree = cKDTree(sup)
    d1, i1 = tree.query(p)
    held = d1 <= params.lip
    todo = np.flatnonzero(~held & (d1 < params.bridge_span))
    if todo.size:
        balls = tree.query_ball_point(p[todo], r=params.bridge_span)
        for n, (j, ball) in enumerate(zip(todo, balls)):
            if not ball:
                continue
            rel = sup[ball] - p[j]
            dist = np.linalg.norm(rel, axis=1)
            ok = (rel @ (sup[i1[j]] - p[j]) < 0.0) & (dist + d1[j] <= params.bridge_span)
            held[j] = bool(ok.any())
    return held


def tilt_points(inputs: FieldInputs, params: PlanarParams, up: np.ndarray) -> np.ndarray:
    """Points de surface ou le surplomb exige d'incliner les couches (critere d'un trancheur classique). Une face en
    surplomb au-dela de `overhang_deg` (depuis la verticale) s'imprime quand meme a plat si le point est a moins de
    `lip` mm d'un appui (levre, chanfrein, rebord : ce qui compte est la distance en porte-a-faux, pas l'aire), ou
    s'il a des appuis de part et d'autre a moins de `bridge_span` mm (pont : plafond tendu, haut d'une arche ou d'un
    trou). Un porte-a-faux n'a d'appui que d'un cote : il reste a incliner. Appui : bord du surplomb sous lequel une paroi descend ; la ou elle monte, le bord est libre (bout d'un
    porte-a-faux). `OverhangMap.is_bridge` (rien a moins de 2 mm sous la face) ne distingue pas un pont d'un
    porte-a-faux : il n'est pas utilise."""
    md = inputs.mesh_data
    faces = np.asarray(md.faces)
    v = np.asarray(md.vertices, dtype=float)
    over, down = _overhang_faces(md, inputs.platform_faces, params, up)
    idx = np.flatnonzero(over)
    if not idx.size:
        return np.zeros((0, 3))
    sup = _supports(v, faces[idx], up, _Solid(inputs.tet_mesh), _printable_surface(v, faces, over, down))
    pts, owner = _sample_faces(v, faces[idx], 0.5)
    face = idx[owner]
    ceiling = down[face] > math.sin(math.radians(params.ceiling_deg))
    held = np.zeros(len(pts), dtype=bool)
    held[ceiling] = _anchored(pts[ceiling], sup, params)
    n = np.asarray(md.face_normals, dtype=float)[face[~ceiling]]
    held[~ceiling] = _slope_held(pts[~ceiling], n, sup, params, up)
    return pts[~held]


def _slope_held(p: np.ndarray, normals: np.ndarray, sup: np.ndarray, params: PlanarParams, up: np.ndarray):
    """Points d'une PENTE portes : un appui a moins de `lip` mm (levre etroite), ou a moins de `slope_reach` mm dans le
    cone de descente (+-60 deg autour de la ligne de plus grande pente). Une pente n'est portee que par ce qui est sous elle : les parois qui la bordent sur
    les cotes (paroi arriere d'une cabine inclinee entre ses deux flancs) ne la portent pas, alors qu'elles feraient
    passer pour un pont tout plafond tendu entre elles. Le haut d'une arche reste porte : ses flancs, plus raides,
    sont dans le sens de la descente."""
    if not len(p) or not len(sup):
        return np.zeros(len(p), dtype=bool)
    n = normals / np.maximum(np.linalg.norm(normals, axis=1), 1e-12)[:, None]
    desc = -(up[None, :] - (n @ up)[:, None] * n)
    desc /= np.maximum(np.linalg.norm(desc, axis=1), 1e-12)[:, None]
    tree = cKDTree(sup)
    d1, _ = tree.query(p)
    held = d1 <= params.lip                 # levre etroite (rebord arrondi), portee de tous cotes
    for j, ball in enumerate(tree.query_ball_point(p, r=params.slope_reach)):
        if ball and not held[j]:
            rel = sup[ball] - p[j]
            dist = np.maximum(np.linalg.norm(rel, axis=1), 1e-9)
            held[j] = bool(((rel @ desc[j]) >= 0.5 * dist).any())
    return held


def tilt_weight(inputs: FieldInputs, params: PlanarParams, up: np.ndarray) -> np.ndarray:
    """Poids w par tetraedre (0 = couches planes, 1 = champ de la chaleur)."""
    tm = inputs.tet_mesh
    centers = tm.tet_vertices[tm.tetrahedra].mean(axis=1)
    need = tilt_points(inputs, params, up)
    if len(need):
        # recherche bornee a core + ramp (au-dela, w = 0) et parallele : 4,6 s -> une fraction sur turkish
        dist, _ = cKDTree(need).query(centers, distance_upper_bound=params.core + params.ramp, workers=-1)
        w = np.where(np.isfinite(dist), 1.0 - _smoothstep((dist - params.core) / max(params.ramp, 1e-9)), 0.0)
    else:
        w = np.zeros(len(centers))
    plate = float(np.min(tm.tet_vertices[inputs.platform_nodes] @ up)) if len(inputs.platform_nodes) else 0.0
    w *= _smoothstep((centers @ up - plate - params.flat_height) / max(params.flat_height, 1e-9))
    return w


def blend(inputs: FieldInputs, values: np.ndarray, params: PlanarParams, solver_params):
    """(champ, champ d'impression, niveaux) : champ reintegre a partir de la direction melangee, `values` etant le
    champ de la chaleur (par noeud). Le champ d'impression est celui d'avant les couches conformes (il donne l'axe de
    la tete) ; None sans couches conformes. Niveaux : `step4_conform.offset` (vide sans couches conformes)."""
    from element_ardor import step4_field_solver as s4
    up = np.asarray(solver_params.up_vector if solver_params.up_vector is not None else [0.0, 0.0, 1.0], float)
    up = up / np.linalg.norm(up)
    tm = inputs.tet_mesh
    grad = s4._compute_per_tet_gradient(tm, inputs.shape_gradients, values)
    norm = np.linalg.norm(grad, axis=1)
    heat = np.where(norm[:, None] > 1e-12, grad / np.maximum(norm, 1e-12)[:, None], up[None, :])
    w = tilt_weight(inputs, params, up)
    flow = (1.0 - w)[:, None] * up[None, :] + w[:, None] * heat
    fn = np.linalg.norm(flow, axis=1)
    flow = np.where(fn[:, None] > 1e-9, flow / np.maximum(fn, 1e-9)[:, None], up[None, :])
    if params.smooth_iterations > 0:
        flow = _smooth_directions(tm, flow, w, inputs.volumes, int(params.smooth_iterations))
    if solver_params.heat_project_reachable and solver_params.radial_tilt_origin is not None:
        # meme projection que le flux de la chaleur : le gradient du champ obtenu n'est pas exactement dans le
        # demi-plan atteignable par la tete (C impose par la position), la direction melangee y est ramenee
        flow, _changed = s4._project_to_reachable_radial_tilt(tm, flow, solver_params.radial_tilt_origin, up)
    L = s4._assemble_stiffness_matrix(tm, inputs.shape_gradients, inputs.volumes)
    rhs = s4._assemble_divergence_rhs(tm, inputs.shape_gradients, inputs.volumes, flow)
    # Ancrage : chaque noeud du plateau a sa HAUTEUR au-dessus du plateau, pas 0. Les faces de contact comprennent les
    # chanfreins du fond (Benchy : noeuds jusqu'a 0,48 mm) ; tous a 0, ils tordaient les premieres couches et
    # decalaient tout le champ (phi - z de -0,12 +- 0,05 mm sur un Benchy entierement plan).
    pn = inputs.platform_nodes
    height = tm.tet_vertices @ up
    target = np.zeros(tm.node_count)
    target[pn] = height[pn] - height[pn].min()
    anchor = np.zeros(tm.node_count)
    anchor[pn] = solver_params.platform_anchor_weight
    A = (L + sp.diags(anchor, format="csr")).tocsr()
    phi = _solve(A, rhs + anchor * target, solver_params, values)
    if params.conform is None:
        return _place_unanchored(tm, pn, phi, up), None, np.zeros(0)
    from element_ardor import step4_conform
    shaped, levels = step4_conform.offset(inputs, phi, params.conform, solver_params, up)
    return _place_unanchored(tm, pn, shaped, up), _place_unanchored(tm, pn, phi, up), levels


def _smooth_directions(tm: TetMesh, flow: np.ndarray, w: np.ndarray, volumes: np.ndarray, iterations: int):
    """Directions lissees dans la zone inclinee (w > 0) : moyenne ponderee par le volume aux noeuds, puis moyenne des
    quatre noeuds par tetraedre, `iterations` fois. La ou deux fronts de la chaleur se rejoignent (au-dessus d'une
    ouverture, contournee des deux cotes), les couches se croisaient en angle, une ligne de separation ; lissees, elles
    passent en courbe. Les tetraedres plans (w = 0) gardent la verticale exacte."""
    T = tm.tetrahedra
    vol = np.abs(volumes) + 1e-12
    inc = sp.csr_matrix((np.ones(T.size), (np.repeat(np.arange(len(T)), 4), T.ravel())), shape=(len(T), tm.node_count))
    to_node = sp.diags(1.0 / np.maximum(inc.T @ vol, 1e-30)) @ inc.T @ sp.diags(vol)    # moyenne ponderee aux noeuds
    to_tet = inc * 0.25                                                                  # moyenne des 4 noeuds
    out = flow.copy()
    tilted = (w > 0.0)[:, None]
    for _ in range(iterations):
        tet = to_tet @ (to_node @ out)
        tet /= np.maximum(np.linalg.norm(tet, axis=1), 1e-12)[:, None]
        out = np.where(tilted, tet, flow)
    return out


def _solve(A: sp.csr_matrix, b: np.ndarray, solver_params, guess: np.ndarray) -> np.ndarray:
    """A.phi = b : noyau natif (`solve_linear_system`, direct sous `direct_solve_node_threshold` noeuds, sinon gradient
    conjugue) quand il est actif ; repli sur le solveur Python s'il est absent ou ne converge pas."""
    from core import _cpp
    from element_ardor import step4_field_solver as s4
    if _cpp.enabled() and hasattr(_cpp.module(), "solve_linear_system"):
        A = A.tocsr()
        A.sort_indices()
        values, info = _cpp.module().solve_linear_system(
            np.ascontiguousarray(A.indptr, dtype=np.int32), np.ascontiguousarray(A.indices, dtype=np.int32),
            np.ascontiguousarray(A.data, dtype=np.float64), np.ascontiguousarray(b, dtype=np.float64),
            A.shape[0] < solver_params.direct_solve_node_threshold, float(solver_params.convergence_tolerance),
            int(solver_params.max_solver_iterations), np.ascontiguousarray(guess, dtype=np.float64))
        if info.get("converged", False):
            return np.asarray(values, dtype=np.float64)
    return s4._solve_linear_system(A, b, solver_params, initial_guess=guess)


def _place_unanchored(tm: TetMesh, platform_nodes: np.ndarray, phi: np.ndarray, up: np.ndarray) -> np.ndarray:
    """Un morceau du maillage sans noeud du plateau (eclat detache : 15 noeuds a z = 11,5 en haut du porte-canne d'un
    Benchy TetraNet, element plus fin que la maille) n'a pas d'ancrage : sa valeur resolue est arbitraire. Il prend celle
    du noeud ancre le plus proche plus l'ecart de hauteur (couches planes a travers l'eclat). Avant le 2026-11-25, il
    gardait le champ de la chaleur, qui vaut 0 la ou elle n'arrive pas : l'eclat tombait dans la premiere couche et ses
    lignes de niveau sur le STL (etape 5e) faisaient du porte-canne un emmelement de boucles."""
    rows = np.repeat(np.arange(len(tm.tetrahedra)), 4)
    inc = sp.csr_matrix((np.ones(len(rows)), (rows, tm.tetrahedra.ravel())), shape=(len(tm.tetrahedra), tm.node_count))
    count, label = connected_components(inc.T @ inc, directed=False)
    if count <= 1:
        return phi
    anchored = np.zeros(count, dtype=bool)
    anchored[np.unique(label[platform_nodes])] = True
    fixed = anchored[label]
    loose = np.flatnonzero(~fixed)
    if not len(loose) or not fixed.any():
        return phi
    keep = np.flatnonzero(fixed)
    _d, near = cKDTree(tm.tet_vertices[keep]).query(tm.tet_vertices[loose], workers=-1)
    out = phi.copy()
    out[loose] = phi[keep[near]] + (tm.tet_vertices[loose] - tm.tet_vertices[keep[near]]) @ up
    return out
