"""
scripts/measure_support_check.py

Mesure directe, sur un modele reel, de la matiere extrudee sans rien en
dessous au moment du depot - le defaut que tout le mecanisme de support du
slicer existe pour reduire.

Execute step1->step8 + step9b (pipeline/step9b_support_check.py) et rapporte :
  - unsupported_fraction : fraction de la matiere extrudee sans support au
    moment T ;
  - ordering_defect_count : depots dont le support EXISTE mais arrive plus
    tard dans la trajectoire - un defaut d'ORDRE, corrigible ;
  - true_overhang_count : surplombs geometriques reels (aucune orientation
    de tete atteignable ne les resout) ;
  - tilt_resolved_count : points sans support vertical strict, mais resolus
    par une orientation de tete atteignable.

Rapporte AUSSI la distance totale de trajet a vide et le plus grand saut :
les mecanismes de step6 qui reduisent les defauts d'ordre le font en
decoupant et en reordonnant les polylignes d'une couche, ce qui peut se
payer en trajets a vide supplementaires. Un gain de support obtenu au prix
d'une explosion du trajet a vide n'est pas un gain - la seule facon de le
savoir est de mesurer les deux ensemble, jamais l'un sans l'autre.

Usage :
    .venv/Scripts/python.exe scripts/measure_support_check.py [chemin_stl]
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from pipeline.step9b_support_check import SupportCheckParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

args = [a for a in sys.argv[1:] if not a.startswith("--")]
MODEL_PATH = args[0] if args else "res/test_models/death.stl"


def _travel_stats(points: list) -> tuple[float, float, int]:
    """(distance totale a vide, plus grand saut a vide, nombre de trajets)."""
    total = 0.0
    longest = 0.0
    count = 0
    for i in range(1, len(points)):
        if points[i].extrusion > 0.0:
            continue
        length = float(np.linalg.norm(points[i].position - points[i - 1].position))
        total += length
        longest = max(longest, length)
        count += 1
    return total, longest, count


def measure(label: str) -> None:
    params = PipelineParams()
    params.support_check = SupportCheckParams()

    t0 = time.time()
    result = run_slicing_only(MODEL_PATH, params)
    elapsed = time.time() - t0

    report = result.support_check_report
    travel_total, travel_longest, travel_count = _travel_stats(result.toolpath.points)

    print(f"\n=== {label} ===")
    print(f"  temps total step1->8+9b : {elapsed:.1f}s")
    print(f"  couches: {len(result.layers)} | points de trajectoire: {len(result.toolpath.points)}")
    print(f"  polylignes totales      : {sum(len(l.toolpath_2d or []) for l in result.colored_layers)}")
    print(f"  points testes (en extrusion) : {report.tested_extruded_points}")
    print(f"  depots sans support     : {len(report.events)} ({report.unsupported_fraction:.2%})")
    print(f"  dont defaut d'ordre     : {report.ordering_defect_count}")
    print(f"  dont resolu par inclinaison de tete : {report.tilt_resolved_count}")
    print(f"  dont surplomb reel      : {report.true_overhang_count}")
    print(f"  trajets a vide          : {travel_count} | total {travel_total:.1f}mm | "
          f"plus grand saut {travel_longest:.2f}mm")


if __name__ == "__main__":
    measure(f"support_check sur {MODEL_PATH}")
