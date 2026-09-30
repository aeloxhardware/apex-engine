"""
scripts/diagnose_toolpath_color.py

Diagnostic ad hoc : inspecte les points de ToolpathMachine autour des plus
longs segments (deja identifies par scripts/diagnose_toolpath_jumps.py),
et reproduit EXACTEMENT la logique de couleur de
file_io/visualization_export.py::_toolpath_path3d (is_extruding = les DEUX
extremites du segment ont extrusion > 0), pour verifier si le rendu GLB
colorie ces longs sauts en vert (EXTRUDING_COLOR) ou cyan (TRAVEL_COLOR).

L'utilisateur rapporte visuellement un pic VERT qui depasse du solide -
si le code dit "extrusion=0 => cyan" mais que le pic observe est vert,
soit (a) le point pense etre un trajet a vide a en realite extrusion>0
(bug dans _iter_schedule_points, pas dans le rendu), soit (b) le point
observe visuellement n'est pas celui identifie par le diagnostic precedent.

Usage :
    .venv/Scripts/python.exe scripts/diagnose_toolpath_color.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"
CONTEXT = 4  # points avant/apres a afficher


def main() -> None:
    params = PipelineParams()
    result = run_slicing_only(MODEL_PATH, params)

    points = result.toolpath.points
    positions = np.array([p.position for p in points])
    extrusions = np.array([p.extrusion for p in points])

    deltas = positions[1:] - positions[:-1]
    seg_lengths = np.linalg.norm(deltas, axis=1)

    order = np.argsort(-seg_lengths)
    top_indices = order[:10]

    for seg_idx in top_indices:
        print(f"\n=== Segment idx={seg_idx}  longueur={seg_lengths[seg_idx]:.3f}mm ===")
        lo = max(0, seg_idx - CONTEXT)
        hi = min(len(points), seg_idx + CONTEXT + 2)
        for i in range(lo, hi):
            p = points[i]
            marker = ""
            if i == seg_idx:
                marker = "  <-- segment start (i-1)"
            elif i == seg_idx + 1:
                marker = "  <-- segment end (i)"
            print(f"  i={i:6d}  pos=({p.position[0]:8.2f},{p.position[1]:8.2f},{p.position[2]:8.2f})  "
                  f"extrusion={p.extrusion:.3f}  color_id={p.color_id}{marker}")

        # Reproduit EXACTEMENT la logique de _toolpath_path3d
        i = seg_idx + 1
        is_extruding = points[i].extrusion > 0.0 and points[i - 1].extrusion > 0.0
        color = "VERT (EXTRUDING_COLOR)" if is_extruding else "CYAN (TRAVEL_COLOR)"
        print(f"  => couleur rendue par _toolpath_path3d : {color}")


if __name__ == "__main__":
    main()
