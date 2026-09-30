"""
scripts/diagnose_global_resolution_single.py

Version ciblee de diagnose_global_resolution_effect.py (le balayage 5
resolutions etait trop long) : teste UNE SEULE resolution globale plus
fine que le defaut (1.5), pour verifier rapidement si elle reduit les
anomalies de non-monotonie (27/228 a resolution=1.5) sur death.stl.

Usage :
    .venv/Scripts/python.exe scripts/diagnose_global_resolution_single.py [resolution]
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


def main() -> None:
    resolution = float(sys.argv[1]) if len(sys.argv) > 1 else 0.7

    params = PipelineParams(resolution=resolution)
    params.layers = replace(params.layers, detect_junction_variance=True)

    t0 = time.time()
    result = run_slicing_only(MODEL_PATH, params)
    elapsed = time.time() - t0

    layers = result.colored_layers
    anomalous = [layer for layer in layers if layer.junction_variance_score is not None]

    print(f"\n=== resolution={resolution} ===")
    print(f"  temps total : {elapsed:.1f}s")
    print(f"  tetraedres  : {result.tet_mesh.tet_count}")
    print(f"  couches     : {len(layers)}")
    print(f"  anomalies   : {len(anomalous)}/{len(layers)}")


if __name__ == "__main__":
    main()
