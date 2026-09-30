"""
scripts/diagnose_global_resolution_effect.py

Test de l'hypothese alternative au raffinement local (limite #48) : sur
death.stl, la zone de jonction fusionnee couvre la quasi-totalite du volume
(rayon 47mm sur un modele 30x61x50mm) - le mecanisme point+rayon
(junction_refinement_resolution) ne peut structurellement pas s'y
appliquer (aucun seed 'rest' possible hors de la zone). Teste si un
raffinement GLOBAL (resolution plus fine, sans aucun mecanisme add_region)
resout la meme non-monotonie du champ (limite #47), mesuree via
LayerParams.detect_junction_variance.

Usage :
    .venv/Scripts/python.exe scripts/diagnose_global_resolution_effect.py
"""

from __future__ import annotations

import os
import sys
import time
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"


def _run_and_measure(params: PipelineParams, label: str) -> dict:
    t0 = time.time()
    result = run_slicing_only(MODEL_PATH, params)
    elapsed = time.time() - t0

    layers = result.colored_layers
    anomalous = [layer for layer in layers if layer.junction_variance_score is not None]

    print(f"\n=== {label} ===")
    print(f"  temps total : {elapsed:.1f}s")
    print(f"  tetraedres  : {result.tet_mesh.tet_count}")
    print(f"  couches     : {len(layers)}")
    print(f"  anomalies   : {len(anomalous)}/{len(layers)}")

    return {
        "elapsed": elapsed,
        "tet_count": result.tet_mesh.tet_count,
        "n_layers": len(layers),
        "n_anomalies": len(anomalous),
    }


def main() -> None:
    base_params = PipelineParams()
    base_params.layers = replace(base_params.layers, detect_junction_variance=True)

    baseline = _run_and_measure(base_params, "BASELINE resolution=1.5 (defaut)")

    results = [baseline]
    for resolution in [1.0, 0.7, 0.5, 0.3]:
        params = replace(base_params, resolution=resolution)
        results.append(_run_and_measure(params, f"resolution={resolution}"))

    print("\n=== COMPARAISON ===")
    for r, label in zip(results, ["1.5 (defaut)", "1.0", "0.7", "0.5", "0.3"]):
        print(f"  resolution={label:12s} tetraedres={r['tet_count']:7d}  "
              f"anomalies={r['n_anomalies']:3d}/{r['n_layers']}  temps={r['elapsed']:.1f}s")


if __name__ == "__main__":
    main()
