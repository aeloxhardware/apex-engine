"""
scripts/diagnose_contour_fragmentation.py

Diagnostic ad hoc : les couches 187, 220-222 (zone identifiee par
scripts/diagnose_contour_overflow.py comme celle des quelques points les
plus proches de la surface exterieure, meme si le debordement mesure reste
infime, <0.01mm) sont examinees pour compter le nombre de boucles de
contour et leur taille - verifie l'hypothese que le motif visuel de zigzag
signale par l'utilisateur (screenshots) vient d'une FRAGMENTATION en
nombreuses petites boucles sur une zone a geometrie tres detaillee
(motif festonne), pas d'un vrai debordement geometrique (deja mesure et
quasi nul).

Usage :
    .venv/Scripts/python.exe scripts/diagnose_contour_fragmentation.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"


def main() -> None:
    params = PipelineParams()
    result = run_slicing_only(MODEL_PATH, params)
    layers = result.colored_layers

    loop_counts = np.array([len(layer.contour_3d) for layer in layers])
    print(f"Nombre de boucles par couche : min={loop_counts.min()} max={loop_counts.max()} "
          f"median={np.median(loop_counts):.1f} mean={loop_counts.mean():.2f}")

    order = np.argsort(-loop_counts)
    print("\nTop 15 couches avec le plus de boucles de contour :")
    for i in order[:15]:
        layer = layers[i]
        loop_sizes = [len(loop) for loop in layer.contour_3d]
        z_values = [float(np.mean(np.asarray(loop)[:, 2])) for loop in layer.contour_3d]
        print(f"  layer.id={layer.id:4d}  isovalue={layer.isovalue:.4f}  "
              f"n_loops={loop_counts[i]:3d}  loop_sizes={sorted(loop_sizes)[:10]}"
              f"{'...' if len(loop_sizes) > 10 else ''}  z_mean~={np.mean(z_values):.2f}")

    # Zoom sur les couches specifiquement signalees par le diagnostic
    # precedent (187, 220, 221, 222).
    print("\n=== Detail des couches 187, 220-222 (zone du debordement infime mesure) ===")
    by_id = {layer.id: layer for layer in layers}
    for lid in [187, 220, 221, 222]:
        layer = by_id.get(lid)
        if layer is None:
            print(f"  layer.id={lid} : absent")
            continue
        loop_sizes = [len(loop) for loop in layer.contour_3d]
        print(f"  layer.id={lid}  isovalue={layer.isovalue:.4f}  n_loops={len(layer.contour_3d)}  "
              f"loop_sizes={sorted(loop_sizes)}")


if __name__ == "__main__":
    main()
