"""
scripts/extract_nested_hole_case.py

Extrait un cas reel (contours 2D) ou le plus petit conteneur PAR AIRE d'un
trou est lui-meme un trou - la condition qui faisait planter l'ancienne
version de `_classify_contours` (KeyError sur res/test_models/death.stl).
Instrumente `_classify_contours` pour DETECTER cette condition sans
planter (le fix est deja applique), afin de fabriquer un test de
non-regression fidele au vrai cas plutot qu'une geometrie synthetique
devinee.

Usage :
    .venv/Scripts/python.exe scripts/extract_nested_hole_case.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402
from core.geometry_utils import point_in_polygon_2d  # noqa: E402
from pipeline.step6_toolpath2d import _signed_area_2d  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"


def _find_pathological_case(contours_2d: list[np.ndarray]):
    if not contours_2d:
        return None
    open_loops = [c[:-1] if np.allclose(c[0], c[-1]) else c for c in contours_2d]
    areas = [abs(_signed_area_2d(loop)) for loop in open_loops]
    order = sorted(range(len(open_loops)), key=lambda i: areas[i], reverse=True)
    containers = {i: [] for i in range(len(open_loops))}
    for rank, i in enumerate(order):
        probe_point = open_loops[i][0]
        for j in order[:rank]:
            if point_in_polygon_2d(probe_point, open_loops[j]):
                containers[i].append(j)
    is_hole = {i: (len(containers[i]) % 2 == 1) for i in range(len(open_loops))}

    for i in range(len(open_loops)):
        if not is_hole[i] or not containers[i]:
            continue
        smallest = min(containers[i], key=lambda p: areas[p])
        if is_hole[smallest]:
            return open_loops, is_hole, containers, i, smallest
    return None


def main() -> None:
    import pipeline.step6_toolpath2d as st6

    found = {}

    orig = st6._classify_contours

    def wrapped(contours_2d):
        if "case" not in found:
            result = _find_pathological_case(contours_2d)
            if result is not None:
                found["case"] = (contours_2d, result)
        return orig(contours_2d)

    st6._classify_contours = wrapped

    params = PipelineParams(resolution=1.5)
    run_slicing_only(MODEL_PATH, params)

    if "case" not in found:
        print("Aucun cas pathologique trouve sur ce modele avec ces parametres.")
        return

    contours_2d, (open_loops, is_hole, containers, hole_idx, smallest_container) = found["case"]
    print(f"Trouve : boucle {hole_idx} (trou) a pour plus petit conteneur la boucle "
          f"{smallest_container} (elle-meme un trou)")
    print(f"Nombre de boucles dans cette couche : {len(contours_2d)}")
    for i, loop in enumerate(open_loops):
        marker = ""
        if i == hole_idx:
            marker = "  <- TROU PROBLEMATIQUE"
        elif i == smallest_container:
            marker = "  <- SON PLUS PETIT CONTENEUR (un trou)"
        print(f"  boucle {i}: is_hole={is_hole[i]}, containers={containers[i]}, "
              f"n_points={len(loop)}{marker}")

    np.savez(
        "out/failing_contours.npz",
        **{f"loop_{i}": loop for i, loop in enumerate(open_loops)},
        hole_idx=hole_idx,
        smallest_container=smallest_container,
    )
    print("Sauvegarde dans out/failing_contours.npz")


if __name__ == "__main__":
    main()
