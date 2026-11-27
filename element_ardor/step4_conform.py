"""
Couches CONFORMES : les couches du haut d'une grande surface en pente (toit, pont) s'inclinent jusqu'a lui etre
paralleles (2026-11-25, limite #92 ; idee de l'utilisateur, 2026-11-19 : « influencer l'inclinaison progressive des
couches pour qu'elles s'inclinent jusqu'a etre paralleles aux surfaces »).

En couches planes (`step4_planar`), un dessus en pente sort en escalier : chaque couche y laisse une marche et ses
parois. Ici, le champ plan phi recoit un decalage u harmonique (L u = 0, meme matrice de rigidite que le champ) :
  - sur chaque surface en pente retenue (faces inclinees de `tilt_min_deg` a `tilt_max_deg`, en morceaux connexes de
    largeur 2A/P >= `min_width` mm : un toit, pas le dessus d'un bastingage de 2 mm), aux noeuds a moins de
    `surface_reach` mm d'elle : u = lambda (c - h), h la hauteur du point de la surface le plus proche, c la hauteur
    moyenne (ponderee par l'aire) du morceau. Avec lambda = 1, phi + u = c sur toute la surface : elle EST une couche ;
  - u = 0 a plus de `depth` mm de toute surface retenue, sur le plateau, et pres des surfaces tournees vers le bas
    (`ceiling_tilt_deg`) : le reste de la piece garde ses couches, un plafond reste plat (pont accroche aux parois) ;
  - entre les deux, u passe en douceur (harmonique) : les couches s'inclinent progressivement.
lambda = min(1, `ratio` * depth / demi-hauteur du morceau) : une surface haute de plus de `depth` / `ratio` n'est
pas rendue entierement plane (les couches s'y resserreraient trop d'un cote).

Desactive (None) dans `PlanarParams.conform`. Les contours des couches courbes sont repris sur le STL par l'etape 5e.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


@dataclass
class ConformParams:
    depth: float = 10.0            # mm : au-dela de cette distance d'une surface retenue, les couches ne bougent pas
    tilt_min_deg: float = 2.0      # pente minimale d'une face retenue (sous : dessus plat, deja une couche plane)
    tilt_max_deg: float = 45.0     # pente maximale (au-dela, la surface est une paroi, pas un dessus)
    min_width: float = 5.0         # mm : largeur 2A/P minimale d'un morceau (un dessus etroit garde ses couches planes)
    ratio: float = 0.5             # lambda = min(1, ratio * depth / demi-hauteur du morceau)
    surface_reach: float = 0.75    # mm : noeuds imposes autour de la surface
    sample_spacing: float = 0.4    # mm : echantillonnage des faces retenues
    # Faces tournees vers le bas a plus de cet angle de la verticale (plafonds, dessous) : u = 0 a moins de
    # `surface_reach`. Sans cela, le toit du Benchy (2 mm d'epaisseur) etait decale en entier : ses couches, paralleles
    # au dessus en pente, coupaient en biais le plafond horizontal de la cabine et le refermaient en plusieurs couches,
    # chacune posee dans le vide (longueur deposee sans appui 41 -> 150 mm). 180 = aucune.
    ceiling_tilt_deg: float = 135.0
    # Dessus plats (pente < `tilt_min_deg`) de largeur 2A/P >= ceci (mm) : leur valeur du champ est aussi un niveau
    # (voir `offset`).
    flat_min_width: float = 1.0
    # mm : ecart maximal entre la cible exacte pres d'une surface rendue entierement parallele et l'ancienne pour que
    # la premiere remplace la seconde (voir `_flush`, limite #96).
    flush_tolerance: float = 0.25


def _face_patches(vertices: np.ndarray, faces: np.ndarray, normals: np.ndarray, platform_faces, up: np.ndarray,
                  params: ConformParams) -> tuple[np.ndarray, np.ndarray]:
    """(faces retenues, etiquette de morceau par face retenue, de 0 a n-1)."""
    tilt = np.degrees(np.arccos(np.clip(normals @ up, -1.0, 1.0)))
    sel = (tilt >= params.tilt_min_deg) & (tilt <= params.tilt_max_deg)
    if platform_faces is not None and len(platform_faces):
        sel[np.asarray(platform_faces)] = False
    idx = np.flatnonzero(sel)
    if not len(idx):
        return idx, idx
    tri = faces[idx]
    edges = np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(idx)), 3)
    key = edges[:, 0].astype(np.int64) * (len(vertices) + 1) + edges[:, 1]
    order = np.argsort(key, kind="stable")
    ks, own, edges = key[order], owner[order], edges[order]
    pair = np.flatnonzero(ks[1:] == ks[:-1])
    graph = sp.coo_matrix((np.ones(len(pair)), (own[pair], own[pair + 1])), shape=(len(idx), len(idx)))
    count, label = connected_components(graph, directed=False)
    # largeur 2A/P de chaque morceau ; une arete de bord n'apparait qu'une fois dans le morceau
    first = np.searchsorted(ks, ks)
    border = np.bincount(first, minlength=len(ks))[first] == 1
    length = np.linalg.norm(vertices[edges[:, 0]] - vertices[edges[:, 1]], axis=1)
    perimeter = np.bincount(label[own[border]], length[border], count)
    area = np.bincount(label, _area(vertices, tri), count)
    wide = 2.0 * area / np.maximum(perimeter, 1e-9) >= params.min_width
    keep = wide[label]
    _u, relabel = np.unique(label[keep], return_inverse=True)
    return idx[keep], relabel


def _area(vertices: np.ndarray, tri: np.ndarray) -> np.ndarray:
    a, b, c = vertices[tri[:, 0]], vertices[tri[:, 1]], vertices[tri[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def _targets(vertices: np.ndarray, tri: np.ndarray, label: np.ndarray, up: np.ndarray, params: ConformParams):
    """(points echantillonnes, cible u par point : lambda * (c - h) de son morceau, morceau de chaque point, lambda
    par morceau)."""
    from element_ardor.step4_planar import _sample_faces
    pts, owner = _sample_faces(vertices, tri, params.sample_spacing)
    count = int(label.max()) + 1
    corner_h = vertices[tri] @ up                                   # (F, 3)
    area = _area(vertices, tri)
    centre = np.bincount(label, area * corner_h.mean(axis=1), count) / np.maximum(np.bincount(label, area, count), 1e-12)
    lo = np.full(count, np.inf)
    hi = np.full(count, -np.inf)
    np.minimum.at(lo, label, corner_h.min(axis=1))
    np.maximum.at(hi, label, corner_h.max(axis=1))
    half = np.maximum(centre - lo, hi - centre)
    lam = np.minimum(1.0, params.ratio * params.depth / np.maximum(half, 1e-9))
    patch = label[owner]
    return pts, lam[patch] * (centre[patch] - pts @ up), patch, lam


def _near_ceilings(md, nodes: np.ndarray, up: np.ndarray, params: ConformParams) -> np.ndarray:
    """Noeuds a moins de `surface_reach` d'une face tournee vers le bas (au-dela de `ceiling_tilt_deg`)."""
    from element_ardor.step4_planar import _sample_faces
    normals = np.asarray(md.face_normals, dtype=float)
    down = np.flatnonzero(normals @ up < np.cos(np.radians(params.ceiling_tilt_deg)))
    if not len(down) or params.ceiling_tilt_deg >= 180.0:
        return np.zeros(len(nodes), dtype=bool)
    pts, _owner = _sample_faces(np.asarray(md.vertices, dtype=float), np.asarray(md.faces)[down], params.sample_spacing)
    dist, _i = cKDTree(pts).query(nodes, distance_upper_bound=params.surface_reach, workers=-1)
    return np.isfinite(dist)


def wide_patch_points(mesh_data, up: np.ndarray, params: ConformParams) -> np.ndarray:
    """Points (echantillonnes a `sample_spacing`) des surfaces que les couches conformes rendent paralleles."""
    from element_ardor.step4_planar import _sample_faces
    vertices, faces = np.asarray(mesh_data.vertices, dtype=float), np.asarray(mesh_data.faces)
    idx, _label = _face_patches(vertices, faces, np.asarray(mesh_data.face_normals, dtype=float), None, up, params)
    if not len(idx):
        return np.zeros((0, 3))
    return _sample_faces(vertices, faces[idx], params.sample_spacing)[0]


def offset(inputs, phi: np.ndarray, params: ConformParams, solver_params, up: np.ndarray
           ) -> tuple[np.ndarray, np.ndarray]:
    """(phi + u (voir le module), niveaux). phi est inchange s'il n'y a aucune surface retenue. Niveaux : valeur du
    champ obtenu sur chaque surface parallele aux couches, des plus grandes aux plus petites : surfaces rendues
    entierement paralleles (lambda = 1), dessus plats et plafonds plats. L'etape 5 place chacun a mi-chemin entre
    deux isovaleurs : une isovaleur qui frolait une telle surface (bruit du champ de quelques centiemes) n'en couvrait
    qu'une partie, en dentelle (toit du Benchy a la finesse 0,8, sommet de la cheminee ; 2026-11-26, limite #95), ou,
    sous un plafond, refermait celui-ci en deux couches posees en partie dans le vide (limite #96)."""
    md = inputs.mesh_data
    vertices, faces = np.asarray(md.vertices, dtype=float), np.asarray(md.faces)
    normals = np.asarray(md.face_normals, dtype=float)
    idx, label = _face_patches(vertices, faces, normals, inputs.platform_faces, up, params)
    shaped, levels = phi, []
    if len(idx):
        shaped, levels = _shift(inputs, phi, (vertices, faces[idx], label), params, solver_params, up)
    for low, high in ((0.0, params.tilt_min_deg), (180.0 - params.tilt_min_deg, 180.0)):   # dessus, plafonds
        flat = replace(params, tilt_min_deg=low, tilt_max_deg=high, min_width=params.flat_min_width)
        fidx, flabel = _face_patches(vertices, faces, normals, inputs.platform_faces, up, flat)
        if len(fidx):
            pts, _t, patch, _lam = _targets(vertices, faces[fidx], flabel, up, flat)
            levels += _patch_levels(shaped, inputs.tet_mesh, pts, patch, np.ones(int(flabel.max()) + 1), flat)
    levels.sort(key=lambda item: -item[0])
    return shaped, np.asarray([level for _area_, *level in levels], dtype=float).reshape(-1, 2)


def _shift(inputs, phi: np.ndarray, patches: tuple, params: ConformParams, solver_params,
           up: np.ndarray) -> tuple[np.ndarray, list]:
    """(phi + u, niveaux [(aire, valeur)] des morceaux rendus entierement plans) ; `patches` : (sommets, faces
    retenues, morceau par face)."""
    from element_ardor import step4_field_solver as s4
    from element_ardor.step4_planar import _solve
    md, tm = inputs.mesh_data, inputs.tet_mesh
    pts, target_pts, patch, lam = _targets(*patches, up, params)
    dist, nearest = cKDTree(pts).query(tm.tet_vertices, workers=-1)
    free = np.zeros(tm.node_count, dtype=bool)
    free[inputs.platform_nodes] = True
    free |= dist >= params.depth
    free |= _near_ceilings(md, tm.tet_vertices, up, params)
    near = (dist <= params.surface_reach) & ~free
    target = np.zeros(tm.node_count)
    target[near] = target_pts[nearest[near]]
    weight = np.zeros(tm.node_count)
    weight[near | free] = solver_params.platform_anchor_weight
    L = s4._assemble_stiffness_matrix(tm, inputs.shape_gradients, inputs.volumes)
    A = (L + sp.diags(weight, format="csr")).tocsr()
    u = _solve(A, weight * target, solver_params, np.zeros(tm.node_count))
    clean = _clean_patches(phi + u, tm, pts, target_pts + pts @ up, patch, lam, params)
    mask = near & clean[patch[nearest]]
    if mask.any() and _flush(phi, target, (mask, nearest), (pts, target_pts), inputs, patches, params, up):
        u = _solve(A, weight * target, solver_params, u)
    shaped = phi + u
    return shaped, _patch_levels(shaped, tm, pts, patch, lam, params)


def _clean_patches(values: np.ndarray, tm, pts: np.ndarray, level: np.ndarray, patch: np.ndarray, lam: np.ndarray,
                   params: ConformParams) -> np.ndarray:
    """Par morceau : rendu entierement parallele (lambda = 1) et champ a moins de `flush_tolerance` de c en au moins
    95 % de ses points (au plus `_LEVEL_SAMPLES`). Sur un morceau que le decalage n'atteint pas partout (death : 5 %
    des points a plus d'1 mm), la cible exacte de `_flush` tordait les couches (1 431 -> 1 478 mm sans appui)."""
    clean = np.zeros(len(lam), dtype=bool)
    for j in np.flatnonzero(lam >= 1.0):
        own = np.flatnonzero(patch == j)
        own = own[::max(1, len(own) // _LEVEL_SAMPLES)]
        if len(own):
            error = _field_at(pts[own], tm, values) - level[own]
            clean[j] = np.mean(np.abs(error) <= params.flush_tolerance) >= 0.95
    return clean


def _field_at(points: np.ndarray, tm, values: np.ndarray) -> np.ndarray:
    from core.mesh_types import ScalarField
    from element_ardor.step5e_surface_levels import phi_at
    return phi_at(points, tm, ScalarField(values=values))


def _flush(phi: np.ndarray, target: np.ndarray, links: tuple, samples: tuple, inputs, patches: tuple,
           params: ConformParams, up: np.ndarray) -> bool:
    """Cible exacte pres des surfaces `links[0]` (rendues entierement paralleles et propres, voir `_clean_patches`),
    modifie `target` en place ; vrai si elle a change.

    La cible lambda (c - h), ajoutee a phi, en recopie le bruit : sur le toit du Benchy basse resolution, phi + u
    variait de 0,36 mm, presque deux couches, et le toit sortait en dentelle (limite #96). Le champ vise y est
    c - profondeur sous le plan de la face la plus proche, au gradient median de phi pres. Un noeud ne prend cette
    cible que si elle s'ecarte de moins de `flush_tolerance` de l'ancienne (bruit) : au-dela, le champ a sa propre
    structure (surplomb, plafond) et la remplacer tordait les couches (death : 1 431 -> 1 470 mm sans appui).
    `links` : (noeuds concernes, point de surface le plus proche de chaque noeud) ; `samples` : (points, cible)."""
    from element_ardor import step4_field_solver as s4
    mask, nearest = links
    pts, target_pts = samples
    nodes = np.flatnonzero(mask)
    tm = inputs.tet_mesh
    k = nearest[nodes]
    grad = np.median(np.linalg.norm(s4._compute_per_tet_gradient(tm, inputs.shape_gradients, phi), axis=1))
    exact = target_pts[k] + pts[k] @ up - grad * _depth(patches, tm.tet_vertices[nodes], up) - phi[nodes]
    keep = np.abs(exact - target[nodes]) <= params.flush_tolerance
    target[nodes[keep]] = exact[keep]
    return bool(keep.any())


def _depth(patches: tuple, nodes: np.ndarray, up: np.ndarray) -> np.ndarray:
    """Profondeur de chaque noeud sous le plan de la face retenue la plus proche (`patches` : sommets, faces, morceau
    par face)."""
    import trimesh
    mesh = trimesh.Trimesh(patches[0], patches[1], process=False)
    closest, _distance, face = trimesh.proximity.closest_point(mesh, nodes)
    normal = mesh.face_normals[face] * np.sign(mesh.face_normals[face] @ up)[:, None]
    return np.einsum("ij,ij->i", closest - nodes, normal)


_LEVEL_SAMPLES = 4000


def _patch_levels(values: np.ndarray, tm, pts: np.ndarray, patch: np.ndarray, lam: np.ndarray,
                  params: ConformParams) -> list:
    """[(aire approchee, valeur mediane du champ, etendue p5-p95)] du champ aux points (au plus `_LEVEL_SAMPLES`)
    de chaque morceau de lambda = 1. L'etape 5 juge de l'etendue acceptable (une epaisseur de couche)."""
    from core.mesh_types import ScalarField
    from element_ardor.step5e_surface_levels import phi_at
    out = []
    for j in np.flatnonzero(lam >= 1.0):
        own = np.flatnonzero(patch == j)
        if not len(own):
            continue
        v = phi_at(pts[own[::max(1, len(own) // _LEVEL_SAMPLES)]], tm, ScalarField(values=values))
        low, mid, high = np.percentile(v, [5, 50, 95])
        out.append((len(own) * params.sample_spacing ** 2, float(mid), float(high - low)))
    return out
