"""
Juge des replis de la boucle Bird (porte de research/S4_Slicer/fold_judge.py).

Definition de l'utilisateur (dessin du 2026-09-26) : une barre qui se courbe en
U est saine ; une barre dont la SURFACE passe a travers elle-meme est un repli.
Mesure : triangles de la surface deformee qui se traversent, regroupes en zones
d'un seul tenant ; le critere est l'AIRE DE LA PLUS GRANDE ZONE (Benchy : zones
saines <= 13,2 mm2, plus petit repli de l'utilisateur 69,7 mm2 ; seuil 30 mm2).
Criteres ecartes, tous mesures contre les verdicts de l'utilisateur : tetraedres
retournes, recouvrements entre tetraedres non voisins, enfoncement median,
distance geodesique entre triangles croises.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import trimesh
from scipy.sparse.csgraph import connected_components

from pipeline.step1_ingest import _find_self_intersecting_faces


def crossing_zones(mesh: trimesh.Trimesh) -> list[tuple[float, float, int]]:
    """Zones de croisement, triees par aire decroissante : (aire mm2, etendue mm, faces)."""
    bad = np.array(sorted(_find_self_intersecting_faces(mesh)), dtype=np.int64)
    if len(bad) == 0:
        return []
    index = -np.ones(len(mesh.faces), dtype=np.int64)
    index[bad] = np.arange(len(bad))
    adj = mesh.face_adjacency
    both = (index[adj] >= 0).all(axis=1)
    graph = sp.coo_matrix((np.ones(both.sum()), (index[adj[both, 0]], index[adj[both, 1]])),
                          shape=(len(bad), len(bad)))
    n_zones, label = connected_components(graph, directed=False)
    centers = mesh.triangles_center[bad]
    areas = np.bincount(label, weights=mesh.area_faces[bad], minlength=n_zones)
    counts = np.bincount(label, minlength=n_zones)
    zones = [(float(areas[z]), float(np.ptp(centers[label == z], axis=0).max()), int(counts[z]))
             for z in range(n_zones)]
    return sorted(zones, reverse=True)


def judge(mesh: trimesh.Trimesh) -> dict:
    zones = crossing_zones(mesh)
    return {"largest_zone_mm2": zones[0][0] if zones else 0.0,
            "largest_zone_extent_mm": zones[0][1] if zones else 0.0,
            "zones": len(zones),
            "total_mm2": sum(a for a, _, _ in zones)}


def as_judged_mesh(points: np.ndarray, faces: np.ndarray) -> trimesh.Trimesh:
    """Meme maillage que celui que la recherche jugeait : surface exportee en STL
    (coordonnees en float32) puis relue par trimesh (fusion des sommets)."""
    return trimesh.Trimesh(np.asarray(points, dtype=np.float32).astype(float), faces, process=True)
