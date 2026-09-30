"""
scripts/diagnose_contour_vs_rendered_mesh.py

Diagnostic ad hoc, correction d'un defaut du diagnostic precedent
(diagnose_contour_overflow.py) : ce dernier comparait Layer.contour_3d au
STL SOURCE charge depuis le disque (res/test_models/death.stl), alors que
le GLB affiche result.mesh_data - le mesh INTERNE potentiellement REPARE
par pymeshfix (limite #43 du CLAUDE.md : "0.27% d'ecart de volume, 1014
faces retirees" sur death.stl specifiquement). Comparer contre le mauvais
mesh de reference peut masquer un vrai debordement visuel : un point de
contour peut deborder du mesh REPARE/AFFICHE tout en restant proche du STL
BRUT (si la reparation a legerement deplace la surface a cet endroit),
d'ou le resultat "quasi 0 debordement" (max 0.008mm) du diagnostic
precedent malgre un ecart visuel net rapporte par l'utilisateur sur trois
captures d'ecran distinctes.

Usage :
    .venv/Scripts/python.exe scripts/diagnose_contour_vs_rendered_mesh.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import trimesh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"
OUTSIDE_TOL_MM = 0.005


def main() -> None:
    params = PipelineParams()
    result = run_slicing_only(MODEL_PATH, params)
    layers = result.colored_layers

    # Le mesh REELLEMENT affiche en gris dans le GLB
    # (file_io/visualization_export.py ligne 177-181 :
    # trimesh.Trimesh(vertices=result.mesh_data.vertices, faces=result.mesh_data.faces, ...))
    rendered_mesh = trimesh.Trimesh(
        vertices=result.mesh_data.vertices, faces=result.mesh_data.faces, process=False
    )
    print(f"mesh_data (rendu GLB)  : {len(rendered_mesh.vertices)} sommets, "
          f"{len(rendered_mesh.faces)} faces, volume={rendered_mesh.volume:.3f}")

    source_mesh = trimesh.load(MODEL_PATH, force="mesh")
    print(f"STL source (disque)    : {len(source_mesh.vertices)} sommets, "
          f"{len(source_mesh.faces)} faces, volume={source_mesh.volume:.3f}")
    print(f"Ecart volume : {abs(rendered_mesh.volume - source_mesh.volume) / source_mesh.volume * 100:.3f}%")

    all_points: list[np.ndarray] = []
    layer_of_point: list[int] = []
    for layer in layers:
        for loop in layer.contour_3d:
            pts = np.asarray(loop)
            if len(pts) >= 2 and np.allclose(pts[0], pts[-1]):
                pts = pts[:-1]
            all_points.append(pts)
            layer_of_point.extend([layer.id] * len(pts))

    all_points = np.concatenate(all_points, axis=0)
    layer_of_point = np.array(layer_of_point)

    for label, mesh in [("mesh_data (rendu GLB)", rendered_mesh), ("STL source (disque)", source_mesh)]:
        prox = trimesh.proximity.ProximityQuery(mesh)
        closest, dist, triangle_id = prox.on_surface(all_points)
        face_normals = mesh.face_normals[triangle_id]
        to_point = all_points - closest
        signed = np.einsum("ij,ij->i", to_point, face_normals)
        outside_mask = (dist > OUTSIDE_TOL_MM) & (signed > 0)
        n_outside = int(outside_mask.sum())

        print(f"\n=== Contre {label} : {n_outside}/{len(all_points)} points hors solide "
              f"(tol={OUTSIDE_TOL_MM}mm) ===")

        outside_only_dist = dist[signed > 0]
        if len(outside_only_dist) > 0:
            print(f"  max distance (cote exterieur) = {outside_only_dist.max():.4f}mm  "
                  f"p99.9={np.percentile(outside_only_dist, 99.9):.4f}mm")

        if n_outside > 0:
            outside_dist = dist[outside_mask]
            outside_layers = layer_of_point[outside_mask]
            order = np.argsort(-outside_dist)
            print(f"  Top 10 points les plus eloignes :")
            for i in order[:10]:
                p = all_points[outside_mask][i]
                print(f"    ({p[0]:8.3f}, {p[1]:8.3f}, {p[2]:8.3f})  dist={outside_dist[i]:.4f}mm  "
                      f"layer={outside_layers[i]}")


if __name__ == "__main__":
    main()
