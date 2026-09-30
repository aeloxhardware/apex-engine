"""
scripts/diagnose_toolpath_jumps.py

Diagnostic ad hoc : cherche les segments de ToolpathMachine (sortie step8)
les plus longs sur death.stl, distingue extrusion vs trajet a vide, pour
verifier l'hypothese que l'artefact visuel signale par l'utilisateur (pics
verts/cyan qui depassent nettement du mesh dans out/death_preview.glb) est
un trajet a vide insere en ligne droite (step8._iter_schedule_points, pas
de planification - voir CLAUDE.md) entre deux passes eloignees du
PrintSchedule (DFS spatial de continuite, limite #40), et NON un
debordement geometrique de Layer.contour_3d (deja mesure et ECARTE :
0/64698 points hors solide, voir scripts/diagnose_contour_overflow.py).

Usage :
    .venv/Scripts/python.exe scripts/diagnose_toolpath_jumps.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"
TOP_N = 30


def main() -> None:
    params = PipelineParams()
    result = run_slicing_only(MODEL_PATH, params)

    points = result.toolpath.points
    positions = np.array([p.position for p in points])
    extrusions = np.array([p.extrusion for p in points])

    deltas = positions[1:] - positions[:-1]
    seg_lengths = np.linalg.norm(deltas, axis=1)
    is_travel = ~((extrusions[1:] > 0.0) & (extrusions[:-1] > 0.0))

    order = np.argsort(-seg_lengths)
    print(f"Total segments : {len(seg_lengths)}")
    print(f"Longueur segment : max={seg_lengths.max():.3f}mm "
          f"median={np.median(seg_lengths):.3f}mm mean={seg_lengths.mean():.3f}mm")
    print(f"Nombre de trajets a vide (is_travel) : {int(is_travel.sum())}/{len(is_travel)}")

    print(f"\nTop {TOP_N} segments les plus longs (idx, longueur_mm, is_travel, p_start, p_end):")
    for i in order[:TOP_N]:
        p0 = positions[i]
        p1 = positions[i + 1]
        print(f"  idx={i:6d}  len={seg_lengths[i]:8.3f}mm  travel={bool(is_travel[i])!s:5}  "
              f"start=({p0[0]:7.2f},{p0[1]:7.2f},{p0[2]:7.2f})  "
              f"end=({p1[0]:7.2f},{p1[1]:7.2f},{p1[2]:7.2f})")

    # Correle longueur de segment et is_travel : si les segments longs sont
    # presque tous des trajets a vide, ca confirme l'hypothese "insertion en
    # ligne droite entre passes eloignees", pas un artefact de couche.
    long_threshold = np.percentile(seg_lengths, 99.5)
    long_mask = seg_lengths > long_threshold
    print(f"\nSegments > percentile 99.5 ({long_threshold:.2f}mm) : {int(long_mask.sum())}, "
          f"dont travel={int((long_mask & is_travel).sum())} "
          f"extrusion={int((long_mask & ~is_travel).sum())}")


if __name__ == "__main__":
    main()
