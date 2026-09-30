"""
scripts/measure_travel_routing_effect.py

Mesure directe de l'effet du routage lift-and-traverse (CLAUDE.md limite
#46) sur death.stl : compare, AVANT et APRES activation de
PipelineParams.travel_routing, la distance en ligne droite entre le point
de fin d'un trajet a vide et le point precedent, pour tout trajet dont la
longueur d'origine depassait length_threshold.

Le routage remplace chaque tel trajet par plusieurs points (lift/traverse/
descente) : ce script verifie que le SEGMENT le plus long entre deux points
CONSECUTIFS de la trajectoire routee (jamais la distance totale start->end
d'origine) chute nettement, confirmant que le routage decompose bien le
saut en ligne droite plutot que de le laisser intact.

Usage :
    .venv/Scripts/python.exe scripts/measure_travel_routing_effect.py
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from pipeline.step8b_travel_routing import TravelRoutingParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

import sys as _sys
MODEL_PATH = _sys.argv[1] if len(_sys.argv) > 1 else "res/test_models/death.stl"


def _max_consecutive_travel_gap(points) -> float:
    max_gap = 0.0
    for i in range(1, len(points)):
        if points[i].extrusion <= 0.0 or points[i - 1].extrusion <= 0.0:
            gap = float(np.linalg.norm(points[i].position - points[i - 1].position))
            max_gap = max(max_gap, gap)
    return max_gap


def main() -> None:
    base_params = PipelineParams(resolution=1.5)

    t0 = time.time()
    before = run_slicing_only(MODEL_PATH, base_params)
    elapsed_before = time.time() - t0
    max_gap_before = _max_consecutive_travel_gap(before.toolpath.points)

    routed_params = PipelineParams(
        resolution=1.5,
        travel_routing=TravelRoutingParams(length_threshold=5.0, lift_height=3.0),
    )
    t0 = time.time()
    after = run_slicing_only(MODEL_PATH, routed_params)
    elapsed_after = time.time() - t0
    max_gap_after = _max_consecutive_travel_gap(after.toolpath.points)

    print("\n=== AVANT routage ===")
    print(f"  points trajectoire : {len(before.toolpath.points)}")
    print(f"  plus grand saut a vide (ligne droite) : {max_gap_before:.2f}mm")
    print(f"  temps : {elapsed_before:.1f}s")

    print("\n=== APRES routage (lift_height=3.0, traverse_samples=1) ===")
    print(f"  points trajectoire : {len(after.toolpath.points)}")
    print(f"  plus grand saut a vide (ligne droite) : {max_gap_after:.2f}mm")
    print(f"  trajets encore en collision apres routage : {len(after.travel_routing_diagnostics)}")
    print(f"  temps : {elapsed_after:.1f}s")

    print(f"\n=== COMPARAISON ===")
    print(f"  plus grand saut : {max_gap_before:.2f}mm -> {max_gap_after:.2f}mm")


if __name__ == "__main__":
    main()
