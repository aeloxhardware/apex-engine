"""
Etape 5d : contours exacts des couches PLANES, recoupes dans le STL d'origine (2026-11-19).

Une couche d'Ardor est une isosurface de phi sur le maillage tetraedrique. La ou rien n'oblige les couches a pencher
(limite #88), phi = z a ~0,02 mm pres et la couche est un plan ; son contour reste pourtant celui du maillage
tetraedrique, qui n'epouse la surface qu'a l'echelle de ses mailles (limite #86). L'etape 5c ramene les points du
contour sur la surface en se deplacant DANS le plan de la couche, ce qui est impossible sur une surface presque
parallele a la couche : sur le Benchy, la derniere couche de la cheminee (z = 47,9, dessus a 48,0) coupait des faces
de bord du maillage placees a 47,91-47,92 et sortait en polygone crenele (dessus « en labyrinthe »).

Ici, comme un trancheur classique, boucle par boucle : une boucle plane (etendue en z < `flat_tol`) dont la couche
ne traverse, a moins de `neighbourhood` d'elle, que des tetraedres a gradient de phi vertical (a `angle_tol_deg`
pres) est recoupee par le plan z = sa mediane dans la surface d'origine (`MeshData`), et remplacee par la boucle
FERMEE de la section qui lui correspond dans les deux sens (p99 des distances < `match_tol`), simplifiee a
`simplify_tol` pres. Une boucle sans correspondante (STL non etanche : section ouverte) reste celle de l'isosurface.
Le repere local et l'isovaleur ne changent pas.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np
import trimesh
from scipy.spatial import cKDTree

from core.mesh_types import Layer, ScalarField, TetMesh
from element_ardor import step5e_surface_levels as s5e
from element_ardor.step5e_surface_levels import SurfaceLevelParams

_UP = np.array([0.0, 0.0, 1.0])


@dataclass
class PlanarSectionParams:
    angle_tol_deg: float = 0.5      # gradient de phi vertical a cet angle pres, dans tous les tetraedres traverses
    flat_tol: float = 0.05          # etendue en z maximale du contour d'une couche plane (mm)
    match_tol: float = 0.5          # p99 des distances, dans les deux sens, entre une boucle et sa section (mm)
    # Pres d'un dessus, l'isosurface crenele (mailles coupees de biais) : Benchy, derniere couche de la cheminee, p99
    # 0,53 et 0,69 mm, boucle suivante a plus de 1,9 mm. Une correspondante sans rivale est acceptee jusqu'a ceci (mm).
    loose_match_tol: float = 1.0
    ambiguity_ratio: float = 2.5    # la deuxieme boucle la plus proche doit etre au moins ce multiple plus loin
    close_tol: float = 1e-3         # une section dont les extremites sont plus loin que ceci est ouverte (mm)
    sample_step: float = 0.05       # pas d'echantillonnage des boucles pour la comparaison (mm)
    # simplification (Douglas-Peucker) des boucles recoupees : le STL du Benchy (225 000 faces) donnait des contours
    # a un sommet tous les quelques microns (+28 % de points de trajectoire) ; ecart au STL borne par cette valeur (mm)
    simplify_tol: float = 0.01
    neighbourhood: float = 2.0      # rayon autour d'une boucle ou tous les tetraedres traverses doivent etre verticaux (mm)
    # Etape 5e (2026-11-25, limite #92) : les AUTRES boucles (couches courbes, ou planes pres d'une matiere inclinee)
    # prises sur le STL comme ligne de niveau de phi. None = contours de l'isosurface (et de l'etape 5c).
    surface_levels: Optional[SurfaceLevelParams] = field(default_factory=SurfaceLevelParams)


@dataclass
class PlanarSectionReport:
    planar_layers: int = 0          # couches ayant au moins une boucle plane
    replaced_loops: int = 0
    kept_loops: int = 0
    surface_loops: int = 0          # etape 5e : boucles examinees, puis remplacees par leur ligne de niveau sur le STL
    surface_replaced: int = 0


def _densify(loop: np.ndarray, step: float) -> np.ndarray:
    seg = np.diff(loop, axis=0)
    n = np.maximum(1, np.ceil(np.linalg.norm(seg, axis=1) / step).astype(int))
    i = np.repeat(np.arange(len(seg)), n)
    t = (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)) / np.repeat(n, n)
    return np.vstack([loop[:-1][i] + seg[i] * t[:, None], loop[-1:]])


def _planar_loops(layers: list[Layer], tet_mesh: TetMesh, scalar: ScalarField, params: PlanarSectionParams) -> dict:
    """{(couche, boucle): hauteur} des boucles PLANES : etendue en z < `flat_tol` et gradient de phi vertical dans les
    tetraedres que la couche traverse a moins de `neighbourhood` de la boucle. Une couche peut etre plane d'un cote et
    inclinee de l'autre (Benchy incline : coque plane, cabine inclinee)."""
    if scalar.gradient is None:
        return {}
    tets = tet_mesh.tet_vertices[tet_mesh.tetrahedra]
    centres = tets.mean(axis=1)
    radius = np.linalg.norm(tets - centres[:, None, :], axis=2).max(axis=1)     # un tetraedre touche a radius de son centre
    phi = scalar.values[tet_mesh.tetrahedra]
    lo, hi = phi.min(axis=1), phi.max(axis=1)
    g = scalar.gradient / np.maximum(np.linalg.norm(scalar.gradient, axis=1), 1e-12)[:, None]
    vertical = g @ _UP > np.cos(np.radians(params.angle_tol_deg))
    out = {}
    for k, layer in enumerate(layers):
        if not layer.contour_3d or layer.isovalue is None:
            continue
        cut = np.flatnonzero((lo <= layer.isovalue) & (layer.isovalue <= hi))
        if not len(cut):
            continue
        whole = bool(vertical[cut].all())
        reach = 0.0 if whole else params.neighbourhood + float(radius[cut[~vertical[cut]]].max())
        tilted = None if whole else cKDTree(centres[cut[~vertical[cut]]])
        for j, loop in enumerate(layer.contour_3d):
            loop = np.asarray(loop, dtype=float)
            if np.ptp(loop[:, 2]) >= params.flat_tol:
                continue
            if tilted is not None and np.isfinite(tilted.query(loop, distance_upper_bound=reach)[0]).any():
                continue
            out[(k, j)] = float(np.median(loop[:, 2]))
    return out


def _closed_loops(segments: np.ndarray, transform: np.ndarray, close_tol: float) -> list[np.ndarray]:
    if len(segments) == 0:
        return []
    seg3 = np.concatenate([segments, np.zeros((len(segments), 2, 1))], axis=2)
    seg3 = trimesh.transform_points(seg3.reshape(-1, 3), transform).reshape(-1, 2, 3)
    try:
        path = trimesh.load_path(seg3)
    except Exception:                       # section degeneree : aucune boucle sure
        return []
    out = []
    for loop in path.discrete:
        loop = np.asarray(loop, dtype=float)
        if len(loop) >= 4 and np.linalg.norm(loop[0] - loop[-1]) <= close_tol:
            out.append(loop)
    return out


def _simplify(loop: np.ndarray, tol: float) -> np.ndarray:
    """Boucle fermee simplifiee a `tol` pres dans son plan (z constant), au moins 4 sommets."""
    if tol <= 0:
        return loop
    import shapely
    ring = shapely.LinearRing(loop[:, :2]).simplify(tol, preserve_topology=True)
    xy = np.asarray(ring.coords)
    if len(xy) < 4:
        return loop
    return np.column_stack([xy, np.full(len(xy), float(np.median(loop[:, 2])))])


def _match(loop: np.ndarray, candidates: list[np.ndarray], dense: list[np.ndarray], trees: list, params):
    """Boucle de la section qui correspond a `loop` dans les deux sens, ou None : la plus proche (p99 des distances,
    pire des deux sens), acceptee sous `match_tol`, ou sous `loose_match_tol` si elle est sans rivale (la suivante au
    moins `ambiguity_ratio` fois plus loin)."""
    own = _densify(np.asarray(loop, dtype=float), params.sample_step)
    own_tree = cKDTree(own[:, :2])
    scores = []
    for pts, tree in zip(dense, trees):
        score = np.percentile(tree.query(own[:, :2])[0], 99)
        if score <= params.loose_match_tol:
            score = max(score, np.percentile(own_tree.query(pts[:, :2])[0], 99))
        scores.append(float(score))
    if not scores:
        return None
    order = np.argsort(scores)
    best = scores[order[0]]
    second = scores[order[1]] if len(order) > 1 else np.inf
    if best <= params.match_tol or (best <= params.loose_match_tol and second >= params.ambiguity_ratio * best):
        return candidates[order[0]]
    return None


def run(layers: list[Layer], mesh_data, tet_mesh: TetMesh, scalar_field: ScalarField,
        params: PlanarSectionParams = PlanarSectionParams()) -> tuple[list[Layer], PlanarSectionReport]:
    report = PlanarSectionReport()
    layers, replaced = _run_planar(layers, mesh_data, tet_mesh, scalar_field, params, report)
    if params.surface_levels is not None:
        layers, surface = s5e.run(layers, mesh_data, tet_mesh, scalar_field, replaced, params.surface_levels)
        report.surface_loops, report.surface_replaced = surface.loops, surface.replaced
    return layers, report


def _run_planar(layers, mesh_data, tet_mesh, scalar_field, params, report) -> tuple[list[Layer], set]:
    """Etape 5d proprement dite : (couches, boucles (k, j) remplacees par leur section plane)."""
    planar = _planar_loops(layers, tet_mesh, scalar_field, params)
    report.planar_layers = len({k for k, _j in planar})
    replaced = set()
    if not planar:
        return layers, replaced
    mesh = trimesh.Trimesh(np.asarray(mesh_data.vertices), np.asarray(mesh_data.faces), process=False)
    key = {kj: round(z, 4) for kj, z in planar.items()}
    heights = sorted(set(key.values()))
    lines, transforms, _faces = trimesh.intersections.mesh_multiplane(mesh, [0.0, 0.0, 0.0], _UP, np.asarray(heights))
    sections = {}
    for h, segs, transform in zip(heights, lines, transforms):
        loops = _closed_loops(segs, transform, params.close_tol)
        dense = [_densify(q, params.sample_step) for q in loops]
        sections[h] = (loops, dense, [cKDTree(q[:, :2]) for q in dense])
    out = list(layers)
    for k in sorted({k for k, _j in planar}):
        contours = list(layers[k].contour_3d)
        for j, loop in enumerate(contours):
            if (k, j) not in planar:
                continue
            loops, dense, trees = sections[key[(k, j)]]
            found = _match(loop, loops, dense, trees, params) if loops else None
            if found is None:
                report.kept_loops += 1
            else:
                contours[j] = _simplify(found, params.simplify_tol)
                report.replaced_loops += 1
                replaced.add((k, j))
        out[k] = replace(layers[k], contour_3d=contours)
    return out, replaced
