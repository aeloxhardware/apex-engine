"""Lissage non planaire : motif de Hilbert (dépliage ARAP, localisation UV, grille géodésique).

Découpé de `surface_paths.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `surface_paths`.
"""

import numpy as np
from element_ardor.skin.surface_mesh import (  # noqa: F401
    NOZZLE_WIDTH,
    PatternResult,
    ZoneMesh,
    boundary_edges,
)
from element_ardor.skin.surface_patterns import (  # noqa: F401
    first_level,
    line_direction,
    signed_line_distance,
    zone_center,
)


# --- Hilbert : depliage ARAP (etape 4d) ----------------------------------------
def arap_unfold(zm, iterations=30):
    """Depliage au plus rigide (ARAP, libigl, MPL-2.0) initialise par LSCM
    (deux sommets de bord les plus eloignes epingles). ARAP vise des
    triangles plans ISOMETRIQUES aux triangles 3D : l'echelle est conservee,
    l'ecart restant est mesure par les valeurs singulieres (arap_stretch)."""
    import igl
    V, F = np.ascontiguousarray(zm.V, float), np.ascontiguousarray(zm.F, np.int64)
    bnd = np.unique(boundary_edges(zm.F))
    i0 = int(bnd[np.argmax(np.linalg.norm(zm.V[bnd] - zm.V[bnd].mean(0), axis=1))])
    i1 = int(bnd[np.argmax(np.linalg.norm(zm.V[bnd] - zm.V[i0], axis=1))])
    b = np.array([i0, i1], np.int64)
    bc = np.array([[0.0, 0.0], [np.linalg.norm(zm.V[i1] - zm.V[i0]), 0.0]])
    uv0, _q = igl.lscm(V, F, b, bc)
    data = igl.ARAPData()
    data.max_iter = iterations
    # un sommet epingle (sa position LSCM) : sans aucune contrainte, la
    # precalculation echoue sur le demi-cylindre (systeme singulier) ; deux
    # epingles imposeraient la corde 3D au lieu de la distance sur la surface
    igl.arap_precomputation(V, F, 2, np.array([i0], np.int32), data)
    return np.asarray(igl.arap_solve(np.ascontiguousarray(uv0[[i0]]), data, np.ascontiguousarray(uv0)))


def arap_stretch(zm, uv):
    """Valeurs singulieres (T,2) du passage triangle 3D -> plan (1 = isometrie)
    et aire 3D par triangle."""
    tri = zm.V[zm.F]
    e1, e2 = tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]
    n = np.cross(e1, e2)
    area = 0.5 * np.linalg.norm(n, axis=1)
    x = e1 / np.linalg.norm(e1, axis=1)[:, None]
    y = np.cross(n / (2 * area)[:, None], x)
    dq = np.stack([np.stack([(e1 * x).sum(1), (e1 * y).sum(1)], 1),
                   np.stack([(e2 * x).sum(1), (e2 * y).sum(1)], 1)], axis=2)       # (T,2,2)
    t = uv[zm.F]
    du = np.stack([t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]], axis=2)
    return np.linalg.svd(du @ np.linalg.inv(dq), compute_uv=False), area


def locate_uv(zm, uv, pts2d, ks=(16, 64, 256, 1024)):
    """Triangle et barycentriques de chaque point 2D dans le depliage ; -1 si
    hors. Candidats = triangles aux centres les plus proches, recherche
    elargie (ks) pour les points non trouves : sur un maillage d'aiguilles
    (demi-cylindre), 12 candidats laissaient 8 % de points non localises en
    PLEIN milieu de zone, donc des trous dans la courbe."""
    from scipy.spatial import cKDTree
    tri = uv[zm.F]
    cen = tri.mean(axis=1)
    tree = cKDTree(cen)
    face = np.full(len(pts2d), -1)
    bary = np.zeros((len(pts2d), 3))
    # hors du depliage : plus loin de tout centre que le plus grand rayon de triangle
    rmax = np.linalg.norm(tri - cen[:, None, :], axis=2).max()
    outside = tree.query(pts2d)[0] > rmax + 1e-9
    done = 0
    for k in ks:
        todo_idx = np.flatnonzero((face < 0) & ~outside)
        if not len(todo_idx):
            break
        k = min(k, len(zm.F))
        _d, cand = tree.query(pts2d[todo_idx], k=k)
        cand = cand.reshape(len(todo_idx), -1)
        for c in range(done, cand.shape[1]):
            sub = face[todo_idx] < 0
            if not sub.any():
                break
            ii = todo_idx[sub]
            f = cand[sub, c]
            a, b_, c_ = tri[f, 0], tri[f, 1], tri[f, 2]
            v0, v1, v2 = b_ - a, c_ - a, pts2d[ii] - a
            den = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
            den = np.where(np.abs(den) < 1e-30, 1e-30, den)
            l1 = (v2[:, 0] * v1[:, 1] - v1[:, 0] * v2[:, 1]) / den
            l2 = (v0[:, 0] * v2[:, 1] - v2[:, 0] * v0[:, 1]) / den
            inside = (l1 >= -1e-9) & (l2 >= -1e-9) & (l1 + l2 <= 1 + 1e-9)
            face[ii[inside]] = f[inside]
            bary[ii[inside]] = np.column_stack([1 - l1[inside] - l2[inside], l1[inside], l2[inside]])
        done = cand.shape[1]
        if k == len(zm.F):
            break
    return face, bary


def hilbert(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=0.0, d_edge=None, uv=None, **_):
    """Courbe de Hilbert au pas w tracee dans le depliage ARAP de la zone,
    ramenee sur la surface par barycentriques (points densifies a w/4 : la
    courbe suit la surface entre deux noeuds de grille), rognee au premier
    niveau. Le pas reel sur la surface = w x etirement local du depliage (mesure,
    pas corrige : decision de l'utilisateur apres mesure)."""
    import math
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    uv = arap_unfold(zm) if uv is None else uv
    a = np.radians(angle_deg)
    rot = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    q = uv @ rot                                            # grille tournee de angle_deg
    lo, hi = q.min(axis=0), q.max(axis=0)
    extent = float(max(hi - lo))
    order = max(1, math.ceil(math.log2(max(extent / width, 1.0))))
    n = 2 ** order
    xs, ys = _hilbert_d2xy(n, np.arange(n * n))
    grid = lo + (np.column_stack([xs, ys]) + 0.5) * width   # pas exact w (grille plus grande que la zone)
    # densification a w/4 le long de la courbe
    sub = 4
    tt = np.arange(sub) / sub
    dense = (grid[:-1, None, :] + tt[None, :, None] * (grid[1:] - grid[:-1])[:, None, :]).reshape(-1, 2)
    dense = np.vstack([dense, grid[-1:]])
    face, bary = locate_uv(zm, q, dense)
    ok = face >= 0
    pts = np.zeros((len(dense), 3))
    de = np.full(len(dense), -1.0)
    fv = zm.F[face[ok]]
    pts[ok] = (bary[ok, :, None] * zm.V[fv]).sum(axis=1)
    de[ok] = (bary[ok] * d_edge[fv]).sum(axis=1)
    keep = de >= first_level(width, overlap)
    paths, start = [], None
    for i in range(len(dense) + 1):
        if i < len(dense) and keep[i]:
            start = i if start is None else start
        elif start is not None:
            if i - start >= 2:
                paths.append(pts[start:i])
            start = None
    return PatternResult(paths, [], d_edge)


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


# --- Hilbert : avertissement et remedes (decision utilisateur 2026-09-28) ------
HILBERT_STRETCH_TOLERANCE = 0.05     # etirement local admis (5 % = tolerance d'espacement)


def hilbert_stretch(zm, uv=None):
    """Depliage ARAP et part de surface etiree au-dela de la tolerance."""
    uv = arap_unfold(zm) if uv is None else uv
    sv, area = arap_stretch(zm, uv)
    over = np.abs(sv - 1).max(axis=1) > HILBERT_STRETCH_TOLERANCE
    return uv, float(area[over].sum() / area.sum())


def geodesic_grid_coords(zm, angle_deg, center):
    """Remede C : coordonnees (U, V) = distances signees MESUREES SUR LA
    SURFACE a deux lignes perpendiculaires passant par le centre (meme graine
    que le rectiligne). Chaque famille de lignes de grille est alors exactement
    espacee de w ; les deux familles ne sont plus perpendiculaires sur une
    surface bombee : les carres du motif se deforment."""
    t, normal = line_direction(zm, angle_deg)
    m = np.cross(t, normal)
    m /= np.linalg.norm(m)
    p = zm.V[center]
    return np.column_stack([signed_line_distance(zm, t, p), signed_line_distance(zm, m, p)])


def hilbert_geodesic(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=0.0, d_edge=None, **_):
    """Hilbert, remede C. Les triangles retournes dans (U, V) (les deux champs
    de distance deviennent paralleles loin du centre) sont ecartes de la
    localisation : aucun point n'y est pose deux fois."""
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    uv = geodesic_grid_coords(zm, angle_deg, zone_center(zm, d_edge))
    t = uv[zm.F]
    a, b = t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]
    det = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    good = det * np.sign(np.median(det)) > 0
    sub = ZoneMesh.__new__(ZoneMesh)
    sub.V, sub.F = zm.V, zm.F[good]
    return hilbert(sub, width, overlap, 0.0, d_edge=d_edge, uv=uv)
