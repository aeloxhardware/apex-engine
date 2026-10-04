"""
Ecart entre la SURFACE du modele et la frontiere du volume tetraedrise.

Un TetMesh dont la frontiere s'ecarte de la peau du modele (tetraedres
absents sur une coque fine, par exemple) laisse des trous dans les couches :
les contours de step5 sont les sections des tetraedres, donc identiques des
deux cotes de toute comparaison avec elles (check_step5_sections) et un manque
du maillage passe inapercu. Seul un test contre la surface D'ORIGINE le voit.

Mesure : distance de points de la surface a la frontiere du maillage (faces de
tetraedre qui n'appartiennent qu'a un seul tetraedre), plus le rapport des
volumes. Un test « point decale vers l'interieur dans un tetraedre » a ete
essaye puis ABANDONNE : trop fragile pres des aretes (0,5 % de faux positifs
sur un maillage dont la frontiere epouse la peau a 0,04 mm).
Fonction pure, sans dependance au reste du pipeline.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import trimesh
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


@dataclass
class OutsideCluster:
    area: float            # mm2 de peau a plus de `tolerance` de la frontiere du maillage
    center: np.ndarray     # (3,)


@dataclass
class TetCoverageReport:
    sample_count: int
    surface_area: float
    tolerance: float
    max_distance: float            # mm, pire ecart surface -> frontiere du maillage
    far_area: float                # mm2 estimes de peau a plus de `tolerance`
    far_fraction: float
    boundary_area: float           # aire de la frontiere du maillage (a comparer a surface_area)
    volume_ratio: float            # volume des tetraedres / volume du modele (1.0 = identiques)
    clusters: list[OutsideCluster] = field(default_factory=list)  # plus grands d'abord
    far_points: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))


def tet_boundary_faces(tetrahedra: np.ndarray) -> np.ndarray:
    """Faces (indices de sommets) n'appartenant qu'a un seul tetraedre."""
    faces = np.concatenate([tetrahedra[:, [0, 1, 2]], tetrahedra[:, [0, 1, 3]],
                            tetrahedra[:, [0, 2, 3]], tetrahedra[:, [1, 2, 3]]])
    _, first, count = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_counts=True)
    return faces[first[count == 1]]


def surface_outside_volume(
    vertices: np.ndarray,
    faces: np.ndarray,
    tet_vertices: np.ndarray,
    tetrahedra: np.ndarray,
    tolerance: float | None = None,
    samples_per_mm2: float = 2.0,
    max_samples: int = 60_000,
    cluster_link: float = 2.0,
    seed: int = 0,
) -> TetCoverageReport:
    """
    Part de la surface du modele (`vertices`, `faces`) a plus de `tolerance`
    (mm) de la frontiere du maillage tetraedrique.

    `tolerance` : defaut 5e-3 de la diagonale de la boite englobante (0,35 mm
    pour un modele de 70 mm) -- bien au-dessus de l'ecart d'un maillage
    conforme (0,04 mm mesure sur le Benchy), bien en dessous d'une coque
    manquante.
    """
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    area = float(mesh.area)
    count = int(min(max_samples, max(5_000, area * samples_per_mm2)))
    pts, _ = trimesh.sample.sample_surface(mesh, count, seed=seed)
    if tolerance is None:
        tolerance = 5e-3 * float(np.linalg.norm(np.ptp(vertices, axis=0)))

    boundary = trimesh.Trimesh(vertices=tet_vertices, faces=tet_boundary_faces(tetrahedra), process=False)
    _, dist, _ = trimesh.proximity.closest_point(boundary, pts)

    corners = tet_vertices[tetrahedra]
    tet_volume = float(np.abs(np.einsum("ij,ij->i", np.cross(corners[:, 1] - corners[:, 0],
                                                            corners[:, 2] - corners[:, 0]),
                                        corners[:, 3] - corners[:, 0])).sum() / 6.0)
    far = np.flatnonzero(dist > tolerance)
    per_sample = area / count
    report = TetCoverageReport(
        sample_count=count, surface_area=area, tolerance=tolerance, max_distance=float(dist.max()),
        far_area=len(far) * per_sample, far_fraction=len(far) / count,
        boundary_area=float(boundary.area), volume_ratio=tet_volume / abs(float(mesh.volume)),
        far_points=pts[far],
    )
    if len(far):
        pairs = cKDTree(pts[far]).query_pairs(cluster_link, output_type="ndarray")
        graph = sp.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(far),) * 2)
        _, label = connected_components(graph, directed=False)
        sizes = np.bincount(label)
        for c in np.argsort(-sizes):
            report.clusters.append(OutsideCluster(float(sizes[c] * per_sample), pts[far][label == c].mean(axis=0)))
    return report
