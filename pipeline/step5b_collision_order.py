"""
Etape 5b (opt-in, 2026-10-01) - ORDRE DES COUCHES CONTRAINT PAR LA TETE.

Constat qui la motive (CLAUDE.md, limites #64 et #67) : sur death et turkish,
le pipeline par defaut fait toucher la vraie tete sur ~40 % des points
echantillonnes, et ce contre de la matiere deposee BIEN PLUS TOT et plus
haut que la couche en cours. L'ordre des couches (step5) ne connait que
l'isovaleur et la continuite spatiale : il ignore totalement la geometrie de
la tete. C'est une boucle ouverte, au sens de la section « COMMENT LIRE CE
FICHIER » : le critere qui compte (la tete passe-t-elle ?) n'est evalue
qu'a l'etape 9, apres tout, sans rien reinjecter.

Mecanisme : si la tete, posee sur la couche L, touche la matiere de la couche
M, alors L doit etre imprimee AVANT M (M n'existe pas encore quand la tete
passe sur L). Ces contraintes s'ajoutent a l'ordre de step5 par un tri
topologique STABLE (Kahn, priorite = rang dans l'ordre d'origine) : sans
contrainte, l'ordre est inchange. Un cycle (L avant M ET M avant L) est
rompu en prenant la couche la mieux classee ; le nombre de ruptures est
rendu, jamais masque.

Approximations assumees, toutes du cote de la detection et non du depot :
  - pose de tete sur des points du CONTOUR de L tous les `pose_spacing` mm,
    orientation par l'IK de step8 (meme cinematique, meme convention) avec la
    normale de couche pour direction ;
  - matiere de M = contour densifie + grille interieure (`material_spacing`)
    dans le plan local, sans l'offset hors-plan ;
  - tete echantillonnee tous les `head_sample_spacing` mm, zone de contact
    de la buse retiree comme a l'etape 9.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np

from core.mesh_types import Layer, MachineProfile


@dataclass
class CollisionOrderParams:
    pose_spacing: float = 3.0           # mm entre deux poses de tete sur un contour
    material_spacing: float = 1.5       # mm, echantillonnage de la matiere d'une couche
    head_sample_spacing: float = 6.0    # mm, echantillonnage de la surface de la tete
    margin: float = 1.0                 # mm, distance tete / matiere consideree comme contact
    nozzle_contact_clearance: float = 3.0   # meme defaut que CollisionSimParams
    neighbours: int = 4                 # matiere la plus proche consultee par point de tete
    # Contraintes d'APPUI : si la matiere de M est juste sous un point du
    # contour de L (a `support_depth` le long de -normale, a moins de
    # `support_radius`), alors M avant L - sauf si la tete impose l'inverse
    # (la contrainte de tete prime : une collision est un echec certain).
    # Mesure (limite #67) : sans elles, death 0,227 -> 0,311 % de matiere sans
    # appui ; avec, 0,222 % et collisions 1,2 -> 0,13 %. 0,3 = 1,5 hauteur de
    # couche par defaut. None = desactive.
    support_depth: Optional[float] = 0.3
    support_radius: float = 0.6


@dataclass
class CollisionOrderReport:
    constraints: int        # couples (L avant M) detectes
    moved_layers: int       # couches dont la position change
    cycle_breaks: int       # contraintes violees faute d'ordre compatible


def _densify(polyline: np.ndarray, spacing: float) -> np.ndarray:
    """Points tous les `spacing` au plus le long de la polyligne, extremites
    comprises (vectorise : un appel par polyligne, pas par segment)."""
    polyline = np.asarray(polyline, dtype=float)
    if len(polyline) < 2:
        return polyline.copy()
    seg = np.diff(polyline, axis=0)
    k = np.maximum(1, np.ceil(np.linalg.norm(seg, axis=1) / spacing)).astype(np.int64)
    which = np.repeat(np.arange(len(seg)), k)
    offset = np.arange(len(which)) - np.repeat(np.cumsum(k) - k, k)
    t = ((offset + 1) / np.repeat(k, k))[:, None]
    return np.vstack([polyline[:1], polyline[which] + seg[which] * t])


def _layer_material(layer: Layer, spacing: float) -> np.ndarray:
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    from shapely import contains_xy

    origin, u, v = layer.local_frame
    contours = [np.asarray(c, dtype=float) for c in layer.contour_3d if len(c) >= 3]
    parts = [_densify(c, spacing) for c in contours]
    loops = [np.column_stack([(c - origin) @ u, (c - origin) @ v]) for c in contours]
    try:
        area = unary_union([Polygon(l).buffer(0) for l in loops])
        if not area.is_empty:
            # pair-impair : une grille dans l'union des boucles moins les trous
            minx, miny, maxx, maxy = area.bounds
            gx, gy = np.meshgrid(np.arange(minx, maxx, spacing), np.arange(miny, maxy, spacing))
            gx, gy = gx.ravel(), gy.ravel()
            inside = np.zeros(len(gx), bool)
            for loop in loops:
                inside ^= contains_xy(Polygon(loop), gx, gy)
            if inside.any():
                parts.append(origin + np.outer(gx[inside], u) + np.outer(gy[inside], v))
    except Exception:
        pass
    return np.vstack(parts) if parts else np.zeros((0, 3))


def _layer_poses(layer: Layer, spacing: float) -> np.ndarray:
    poses = []
    for c in layer.contour_3d:
        c = np.asarray(c, dtype=float)
        if len(c) < 2:
            continue
        lengths = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(c, axis=0), axis=1))]
        if lengths[-1] <= 0:
            continue
        s = np.arange(0.0, lengths[-1], spacing)
        poses.append(np.column_stack([np.interp(s, lengths, c[:, k]) for k in range(3)]))
    return np.vstack(poses) if poses else np.zeros((0, 3))


def _rotation_matrices(axes, angles: np.ndarray) -> np.ndarray:
    """(P,3,3) : composition des rotations de chaque pose, axe le plus interne
    applique en premier - meme convention que
    step9._positioned_toolhead_vertices, calculee par lot."""
    count = len(angles)
    total = np.broadcast_to(np.eye(3), (count, 3, 3)).copy()
    for k, spec in enumerate(axes):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        theta = np.radians(angles[:, k])
        c, s_ = np.cos(theta)[:, None, None], np.sin(theta)[:, None, None]
        cross = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        rot = c * np.eye(3) + s_ * cross + (1 - c) * np.outer(axis, axis)
        total = total @ rot          # R0 R1 ... RN : RN agit en premier sur le vecteur
    return total


def _stable_topological_order(n: int, edges: set[tuple[int, int]]) -> tuple[list[int], int]:
    """Kahn avec priorite = indice d'origine ; cycle rompu par la couche au
    plus petit degre entrant restant. Rend (ordre, contraintes violees)."""
    import heapq

    succ: list[list[int]] = [[] for _ in range(n)]
    indeg = [0] * n
    for a, b in edges:
        succ[a].append(b)
        indeg[b] += 1
    ready = [i for i in range(n) if indeg[i] == 0]
    heapq.heapify(ready)
    placed = [False] * n
    order: list[int] = []
    remaining = set(range(n))
    while remaining:
        if not ready:
            # cycle : la couche dont le placement viole le MOINS de contraintes
            # (degre entrant restant minimal), a egalite la mieux classee -
            # heuristique classique d'ensemble d'arcs de retour
            forced = min(remaining, key=lambda r: (indeg[r], r))
            ready = [forced]
        i = heapq.heappop(ready)
        if placed[i]:
            continue
        placed[i] = True
        remaining.discard(i)
        order.append(i)
        for j in succ[i]:
            indeg[j] -= 1
            if indeg[j] == 0 and not placed[j]:
                heapq.heappush(ready, j)
    position = {layer: k for k, layer in enumerate(order)}
    broken = sum(1 for a, b in edges if position[a] > position[b])
    return order, broken


def run(layers: list[Layer], machine_profile: MachineProfile,
        params: CollisionOrderParams) -> tuple[list[Layer], CollisionOrderReport]:
    """Reordonne `layers` (ids renumerotes dans le nouvel ordre)."""
    import trimesh
    from scipy.spatial import cKDTree

    from pipeline import step8_machine_toolpath as s8
    from pipeline import step9_collision_sim as s9

    head = machine_profile.nozzle_model
    if not isinstance(head, trimesh.Trimesh) or len(layers) < 2:
        return list(layers), CollisionOrderReport(0, 0, 0)
    samples = s9._filter_nozzle_contact_zone(
        s9._surface_sample_points(head, params.head_sample_spacing), params.nozzle_contact_clearance)

    material, owner = [], []
    for k, layer in enumerate(layers):
        m = _layer_material(layer, params.material_spacing)
        material.append(m)
        owner.append(np.full(len(m), k))
    material_all = np.vstack(material)
    owner_all = np.concatenate(owner)
    tree = cKDTree(material_all)

    axes = machine_profile.rotation_axes or []
    edges: set[tuple[int, int]] = set()
    support_edges: set[tuple[int, int]] = set()
    for k, layer in enumerate(layers):
        poses = _layer_poses(layer, params.pose_spacing)
        if not len(poses):
            continue
        normal = s8._layer_normal(layer.local_frame)
        if getattr(layer, "deposition_reversed", False):
            normal = -normal
        raw = [s8._RawPoint(position=p, desired_direction=normal, layer_id=k, segment_index=0,
                            color_id=None, feedrate=1.0, extrusion=1.0, is_travel=False) for p in poses]
        if getattr(machine_profile, "position_coupled_axes", None) and axes:
            angles, _ = s8._resolve_all_axis_angles_coupled(
                raw, axes, np.array([0.0, 0.0, 1.0]), machine_profile, s8.Step8Params())
        else:
            angles, _ = s8._resolve_all_axis_angles(
                raw, axes, np.array([0.0, 0.0, 1.0]), machine_profile.platform_origin, s8.Step8Params())
        hit_layers: set[int] = set()
        rotations = _rotation_matrices(axes, np.asarray(angles, dtype=float).reshape(len(poses), len(axes)))
        chunk = max(1, 200_000 // max(len(samples), 1))     # <= ~200 000 points par requete
        for start in range(0, len(poses), chunk):
            rot = rotations[start:start + chunk]
            pts = np.matmul(samples[None, :, :], rot.transpose(0, 2, 1)) + poses[start:start + chunk, None, :]
            hit_layers |= _hits(tree, owner_all, pts.reshape(-1, 3), params)
        hit_layers.discard(k)
        edges.update((k, m) for m in hit_layers)
        if params.support_depth is not None:
            below = poses - normal * params.support_depth
            dist, idx = tree.query(below, k=1, distance_upper_bound=params.support_radius, workers=-1)
            ok = np.isfinite(dist)
            supporters = set(np.unique(owner_all[idx[ok]]).tolist()) - {k}
            support_edges.update((m, k) for m in supporters)

    edges |= {e for e in support_edges if (e[1], e[0]) not in edges}

    order, broken = _stable_topological_order(len(layers), edges)
    moved = sum(1 for new, old in enumerate(order) if new != old)
    reordered = [replace(layers[old], id=new) for new, old in enumerate(order)]
    return reordered, CollisionOrderReport(len(edges), moved, broken)


def _hits(tree, owner_all: np.ndarray, points: np.ndarray, params: CollisionOrderParams) -> set[int]:
    dist, idx = tree.query(points, k=params.neighbours, distance_upper_bound=params.margin, workers=-1)
    idx = np.atleast_2d(idx)
    ok = np.isfinite(np.atleast_2d(dist))
    return set(np.unique(owner_all[idx[ok]]).tolist())
