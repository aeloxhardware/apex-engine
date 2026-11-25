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

from dataclasses import dataclass

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
    """(points echantillonnes, cible u par point) : lambda * (c - h) du morceau de chaque point."""
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
    return pts, lam[patch] * (centre[patch] - pts @ up)


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


def offset(inputs, phi: np.ndarray, params: ConformParams, solver_params, up: np.ndarray) -> np.ndarray:
    """phi + u (voir le module) ; phi inchange s'il n'y a aucune surface retenue."""
    from element_ardor import step4_field_solver as s4
    from element_ardor.step4_planar import _solve
    md, tm = inputs.mesh_data, inputs.tet_mesh
    vertices, faces = np.asarray(md.vertices, dtype=float), np.asarray(md.faces)
    idx, label = _face_patches(vertices, faces, np.asarray(md.face_normals, dtype=float), inputs.platform_faces,
                               up, params)
    if not len(idx):
        return phi
    pts, target_pts = _targets(vertices, faces[idx], label, up, params)
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
    return phi + _solve(A, weight * target, solver_params, np.zeros(tm.node_count))
