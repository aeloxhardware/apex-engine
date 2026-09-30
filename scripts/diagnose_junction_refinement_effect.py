"""
scripts/diagnose_junction_refinement_effect.py

Mesure finale du Livrable 2 (limite #48) : compare le diagnostic de
non-monotonie (LayerParams.detect_junction_variance, limite #47) AVANT et
APRES activation de junction_refinement_resolution sur death.stl, pour
verifier que le raffinement local aux jonctions multi-region reduit
effectivement le nombre de couches anormales (27/228 mesure avant tout
fix).

Usage :
    .venv/Scripts/python.exe scripts/diagnose_junction_refinement_effect.py
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

    before = _run_and_measure(base_params, "AVANT junction_refinement_resolution")

    after_params = replace(
        base_params,
        junction_refinement_resolution=1.0,
        junction_refinement_angle_deg=30.0,
        junction_refinement_margin=3.0,
        # near_platform_resolution desactive pour ce test : la zone de
        # jonction detectee sur death.stl (rayon ~47mm) englobe le seed du
        # raffinement plateau sur ce modele precis - conflit geometrique
        # reel (garde-fou explicite, pas un bug), hors perimetre de cette
        # mesure isolee de l'effet du raffinement de jonction seul.
        near_platform_resolution=None,
    )
    after = _run_and_measure(after_params, "APRES junction_refinement_resolution=1.0")

    print("\n=== COMPARAISON ===")
    print(f"  anomalies   : {before['n_anomalies']} -> {after['n_anomalies']}")
    print(f"  tetraedres  : {before['tet_count']} -> {after['tet_count']}")
    print(f"  temps       : {before['elapsed']:.1f}s -> {after['elapsed']:.1f}s")


if __name__ == "__main__":
    main()
