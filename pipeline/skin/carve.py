"""
P3 sur la voie chaleur : creuser le corps sous les zones de lissage retenues.

Port de `pipeline/bird/skin/skin_carve.py` (branche `bird_old`), refait pour des couches qui sont des isosurfaces
COURBES de phi (la voie Bird tranchait des plans dans un espace deforme) :

  1. psi aux nœuds du maillage tetraedrique = distance a la plus proche face retenue, moins l'epaisseur `t` du
     lissage. psi < 0 : la matiere appartient a la couche de surface, pas au corps. psi est interpole LINEAIREMENT
     dans chaque tetraedre : exact quand la peau retenue est plane (la distance a un plan est lineaire), a
     O(h^2 / R) pres sur une paroi courbe de rayon R, par exces de creusement (la distance a plusieurs faces est
     concave : l'interpolant lui est inferieur) ;
  2. une couche n'est recalculee que si un tetraedre qu'elle coupe a un nœud a psi < 0 ; les autres gardent
     exactement leurs contours de l'etape 5 ;
  3. pour une couche touchee : patch triangule de l'isosurface (comme l'etape 5), psi aux sommets par coordonnees
     barycentriques dans le tetraedre source, DECOUPE des triangles le long de psi = 0 (exacte : psi est lineaire
     sur chaque triangle), puis trace de bord : la boucle fermee de la couche creusee est le bord du patch conserve.

Les points de coupe sont calcules depuis l'arete canonique (indice le plus petit d'abord) : deux triangles qui
partagent une arete produisent le MEME sommet, bit pour bit, donc le patch reste etanche.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np
import trimesh
from scipy.spatial import cKDTree

from core.mesh_types import Layer, ScalarField, TetMesh
from pipeline import step5_layers as s5
from pipeline.skin.zones import barycentric


@dataclass
class CarveParams:
    thickness: float = 0.2          # epaisseur de la couche de lissage (mm)
    band: float = 8.0               # au-dela de cette distance a la peau retenue, psi n'est plus calcule (mm)
    refine_edge: float = 0.25       # arete visee des triangles autour de la ligne de coupe (mm)
    max_levels: int = 4             # subdivisions successives au plus (1,5 mm -> 0,1 mm)
    weld_tolerance_factor: float = 1e-6


@dataclass
class CarveReport:
    touched_layers: int = 0
    dropped_layers: int = 0         # couches entierement dans la couche de surface
    nodes_in_skin: int = 0          # nœuds a psi < 0
    seconds: dict = field(default_factory=dict)


class ZoneDistance:
    """Distance EXACTE (3D) a la peau retenue, moins l'epaisseur : psi < 0 dans la couche de surface. Au-dela de
    `band` de la peau, psi vaut `band - t` sans calcul."""

    def __init__(self, zone_vertices: np.ndarray, zone_faces: np.ndarray, params: CarveParams):
        self.params = params
        self.empty = len(zone_faces) == 0
        self.calls = 0
        if self.empty:
            return
        self.mesh = trimesh.Trimesh(zone_vertices, zone_faces, process=False)
        centers = zone_vertices[zone_faces].mean(axis=1)
        self.tree = cKDTree(centers)
        self.reach = params.band + float(np.linalg.norm(zone_vertices[zone_faces] - centers[:, None, :], axis=2).max())

    def psi(self, points: np.ndarray) -> np.ndarray:
        t, band = self.params.thickness, self.params.band
        out = np.full(len(points), band - t)
        if self.empty or len(points) == 0:
            return out
        idx = np.flatnonzero(self.tree.query(points)[0] <= self.reach)
        for s in range(0, len(idx), 50_000):
            part = idx[s:s + 50_000]
            _c, dist, _t = trimesh.proximity.closest_point(self.mesh, points[part])
            out[part] = np.minimum(dist, band) - t
            self.calls += len(part)
        return out


def skin_psi(tet_mesh: TetMesh, zone_vertices: np.ndarray, zone_faces: np.ndarray, params: CarveParams) -> np.ndarray:
    """psi aux nœuds du maillage (distance exacte a la peau retenue, moins l'epaisseur)."""
    return ZoneDistance(zone_vertices, zone_faces, params).psi(tet_mesh.tet_vertices)


def refine_patch(vertices: np.ndarray, faces: np.ndarray, marked: np.ndarray):
    """Subdivision CONFORME : les triangles `marked` sont coupes en 4 (milieux de leurs aretes) ; un voisin non
    marque dont une arete est coupee est coupe en 2 ou 3 (fermeture), jamais de sommet pendant. Rend
    (vertices, faces, parent, red) : triangle d'origine de chaque nouveau triangle, et `red` vrai pour les 4 enfants
    d'un triangle marque (les seuls a pouvoir etre raffines encore)."""
    if not marked.any():
        return vertices, faces, np.arange(len(faces)), np.ones(len(faces), bool)
    n = len(vertices)
    edges = np.concatenate([faces[marked][:, [0, 1]], faces[marked][:, [1, 2]], faces[marked][:, [2, 0]]])
    edges = np.sort(edges, axis=1)
    keys, first = np.unique(edges[:, 0].astype(np.int64) * n + edges[:, 1], return_index=True)
    lo, hi = edges[first, 0], edges[first, 1]
    mid_vertices = 0.5 * (vertices[lo] + vertices[hi])
    mid_of = {int(k): n + i for i, k in enumerate(keys)}
    sorted_keys = keys
    def has(a, b):
        k = np.minimum(a, b).astype(np.int64) * n + np.maximum(a, b)
        i = np.searchsorted(sorted_keys, k)
        return (i < len(sorted_keys)) & (sorted_keys[np.minimum(i, len(sorted_keys) - 1)] == k)

    split = has(faces[:, 0], faces[:, 1]) | has(faces[:, 1], faces[:, 2]) | has(faces[:, 2], faces[:, 0])
    out_faces = [faces[~split]]
    out_parent = [np.flatnonzero(~split)]
    out_red = [np.ones((~split).sum(), bool)]

    def mid(a, b):
        return mid_of.get(int(min(a, b)) * n + int(max(a, b)))

    new_f, new_p, new_r = [], [], []
    for i in np.flatnonzero(split):
        a, b, c = (int(x) for x in faces[i])
        mab, mbc, mca = mid(a, b), mid(b, c), mid(c, a)
        m = (mab is not None) + (mbc is not None) + (mca is not None)
        if m == 3:
            for t in ((a, mab, mca), (mab, b, mbc), (mca, mbc, c), (mab, mbc, mca)):
                new_f.append(t); new_p.append(i); new_r.append(True)
        else:
            for rot in range(3):   # fermeture : tourner le triangle jusqu'a ce que (x, y) soit une arete coupee
                x, y, z = (a, b, c) if rot == 0 else (b, c, a) if rot == 1 else (c, a, b)
                mxy, myz = mid(x, y), mid(y, z)
                if m == 1 and mxy is not None:
                    new_f += [(x, mxy, z), (mxy, y, z)]; new_p += [i, i]; new_r += [False, False]
                    break
                if m == 2 and mxy is not None and myz is not None:
                    new_f += [(mxy, y, myz), (x, mxy, myz), (x, myz, z)]; new_p += [i] * 3; new_r += [False] * 3
                    break
    out_faces.append(np.array(new_f, dtype=np.int64).reshape(-1, 3))
    out_parent.append(np.array(new_p, dtype=np.int64))
    out_red.append(np.array(new_r, dtype=bool))
    return (np.concatenate([vertices, mid_vertices]), np.concatenate(out_faces), np.concatenate(out_parent),
            np.concatenate(out_red))


def clip_patch(vertices: np.ndarray, faces: np.ndarray, psi: np.ndarray):
    """Garde la partie psi >= 0 d'un patch triangule. Rend (vertices, faces, source) ou `source` donne, par
    triangle conserve, l'indice du triangle d'origine."""
    pos = psi[faces] >= 0.0
    count = pos.sum(axis=1)
    keep_all = count == 3
    one = np.flatnonzero(count == 1)      # 1 sommet conserve -> 1 triangle
    two = np.flatnonzero(count == 2)      # 2 sommets conserves -> quadrilatere -> 2 triangles

    new_vertices = [vertices]
    next_id = len(vertices)
    cut_ids: dict = {}

    def cut(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """Indices des sommets de coupe sur les aretes (a, b) (psi de signes opposes), canoniques et partages."""
        nonlocal next_id
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        key = lo.astype(np.int64) * len(vertices) + hi
        uniq, first, inverse = np.unique(key, return_index=True, return_inverse=True)
        ids = np.empty(len(uniq), dtype=np.int64)
        fresh = []
        for n, k in enumerate(uniq):
            existing = cut_ids.get(int(k))
            if existing is None:
                cut_ids[int(k)] = next_id + len(fresh)
                ids[n] = next_id + len(fresh)
                fresh.append(n)
            else:
                ids[n] = existing
        if fresh:
            l, h = lo[first[fresh]], hi[first[fresh]]
            t = psi[l] / (psi[l] - psi[h])
            new_vertices.append(vertices[l] + t[:, None] * (vertices[h] - vertices[l]))
            next_id += len(fresh)
        return ids[inverse]

    out_faces = [faces[keep_all]]
    out_source = [np.flatnonzero(keep_all)]
    if len(one):
        f = faces[one]
        p = pos[one]
        k = p.argmax(axis=1)                       # le sommet conserve
        r = np.arange(len(one))
        a = f[r, k]
        b = f[r, (k + 1) % 3]
        c = f[r, (k + 2) % 3]
        ab, ac = cut(a, b), cut(a, c)
        out_faces.append(np.stack([a, ab, ac], axis=1))      # meme sens que (a, b, c)
        out_source.append(one)
    if len(two):
        f = faces[two]
        p = pos[two]
        k = (~p).argmax(axis=1)                    # le sommet supprime
        r = np.arange(len(two))
        d = f[r, k]
        a = f[r, (k + 1) % 3]
        b = f[r, (k + 2) % 3]
        da, db = cut(d, a), cut(d, b)              # sur (d, a) et (d, b)
        out_faces.append(np.stack([a, b, db], axis=1))
        out_faces.append(np.stack([a, db, da], axis=1))
        out_source.append(two)
        out_source.append(two)
    return np.concatenate(new_vertices), np.concatenate(out_faces).astype(np.int64), np.concatenate(out_source)


def carve_layers(
    layers: list[Layer],
    tet_mesh: TetMesh,
    scalar_field: ScalarField,
    zone: ZoneDistance,
    params: Optional[CarveParams] = None,
) -> tuple[list[Layer], CarveReport]:
    """Couches du corps creuse : les couches non touchees sont rendues telles quelles (meme objet), les touchees
    recoivent des `contour_3d` et un `local_frame` recalcules. Une couche entierement creusee est retiree. Les `id`
    sont inchanges (l'appelant les renumerote s'il le veut).

    Une couche est touchee si un tetraedre qu'elle coupe est PROCHE de la peau retenue (un nœud a moins de t + son
    arete maximale : majorant sur, pas une estimation). Pour chaque couche touchee, le patch de l'isosurface est
    RAFFINE (subdivision conforme) autour de la ligne psi = 0 jusqu'a `refine_edge`, avec psi EXACT a chaque
    sommet, puis decoupe : l'erreur de la ligne de coupe ne depend plus de la taille des tetraedres."""
    params = params or CarveParams()
    report = CarveReport()
    t0 = time.time()
    values = scalar_field.values
    tets = tet_mesh.tetrahedra
    tet_min, tet_max = values[tets].min(axis=1), values[tets].max(axis=1)
    tet_vertices = tet_mesh.tet_vertices[tets]                                   # (T,4,3)
    longest = np.linalg.norm(tet_vertices[:, :, None, :] - tet_vertices[:, None, :, :], axis=3).max(axis=(1, 2))
    node_psi = zone.psi(tet_mesh.tet_vertices)
    near_tet = node_psi[tets].min(axis=1) <= longest                             # psi_min <= arete max : majorant
    volumes = s5._tet_volumes(tet_mesh)
    scale = float(np.ptp(tet_mesh.tet_vertices, axis=0).max()) or 1.0
    weld = params.weld_tolerance_factor * scale
    report.nodes_in_skin = int((node_psi < 0).sum())
    out: list[Layer] = []
    for layer in layers:
        iso = float(layer.isovalue)
        cut_tets = (tet_min <= iso) & (iso <= tet_max)
        if not (cut_tets & near_tet).any():
            out.append(layer)
            continue
        report.touched_layers += 1
        vertices, faces, source_tets = s5._extract_isosurface_triangles(tet_mesh, values, iso, tet_min, tet_max)
        if len(faces) == 0:
            out.append(layer)
            continue
        w_vertices, w_faces = s5._weld_vertices(vertices, faces, weld)
        ok = (w_faces[:, 0] != w_faces[:, 1]) & (w_faces[:, 1] != w_faces[:, 2]) & (w_faces[:, 0] != w_faces[:, 2])
        w_faces, w_source = w_faces[ok], source_tets[ok]
        near_face = near_tet[w_source]

        # psi exact aux sommets des triangles proches ; les autres : loin de la peau (jamais coupes)
        psi = np.full(len(w_vertices), params.band - params.thickness)
        used = np.unique(w_faces[near_face])
        psi[used] = zone.psi(w_vertices[used])
        parent = np.arange(len(w_faces))
        red = np.ones(len(w_faces), bool)
        for _level in range(params.max_levels):
            tri = w_vertices[w_faces]
            edge = np.linalg.norm(tri - np.roll(tri, 1, axis=1), axis=2).max(axis=1)
            pf = psi[w_faces]
            # 1-Lipschitz : psi ne peut s'annuler dans le triangle que si min(psi) <= arete et max(psi) >= -arete
            mark = red & (edge > params.refine_edge) & (pf.min(axis=1) <= edge) & (pf.max(axis=1) >= -edge)
            if not mark.any():
                break
            n_before = len(w_vertices)
            w_vertices, w_faces, new_parent, new_red = refine_patch(w_vertices, w_faces, mark)
            w_source = w_source[new_parent]
            psi = np.concatenate([psi, np.full(len(w_vertices) - n_before, params.band - params.thickness)])
            fresh = np.arange(n_before, len(w_vertices))
            psi[fresh] = zone.psi(w_vertices[fresh])
            red = new_red & np.isin(new_parent, np.flatnonzero(mark)) | (new_red & red[new_parent] & ~mark[new_parent])
        c_vertices, c_faces, c_source = clip_patch(w_vertices, w_faces, psi)
        c_vertices, c_faces = s5._weld_vertices(c_vertices, c_faces, weld)
        contours = s5._extract_closed_contours(c_vertices, c_faces) if len(c_faces) else []
        if not contours:
            report.dropped_layers += 1
            continue
        used_tets = np.unique(w_source[c_source])
        frame = s5._compute_local_frame(contours, scalar_field.gradient[used_tets], volumes[used_tets])
        out.append(replace(layer, contour_3d=contours, local_frame=frame))
    report.seconds["carve"] = time.time() - t0
    report.seconds["distance_calls"] = zone.calls
    return out, report
