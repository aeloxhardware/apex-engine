"""Pas 21a : step2 natif contre Python sur un modele reel (memes surplombs et
ponts, risque a 1e-12, portee au signe pres)."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from core import _cpp  # noqa: E402
from pipeline import step1_ingest, step2_overhang as s2  # noqa: E402
from main import PipelineParams  # noqa: E402


def main():
    p = PipelineParams(resolution=0.8)
    mesh_data, _tet, _rep = step1_ingest.run(sys.argv[1], p.resolution)
    _cpp.set_enabled(False)
    ref = s2.run(mesh_data, p.overhang)
    _cpp.set_enabled(True)
    nat = s2.run(mesh_data, p.overhang)
    assert np.array_equal(ref.is_overhang, nat.is_overhang)
    assert np.array_equal(ref.is_bridge, nat.is_bridge)
    dr = np.abs(ref.risk_score - nat.risk_score).max()
    m = ref.is_bridge
    nan_ok = np.array_equal(np.isnan(ref.bridge_span_direction), np.isnan(nat.bridge_span_direction))
    ds = np.abs(ref.bridge_span_direction[m] - nat.bridge_span_direction[m]).max() if m.any() else 0.0
    print(f"{os.path.basename(sys.argv[1])} : {m.sum()} faces-ponts, {ref.is_overhang.sum()} surplombs IDENTIQUES ; "
          f"risque ecart max {dr:.1e} ; portee ecart max (signe compris) {ds:.1e}, NaN identiques {nan_ok}")


if __name__ == "__main__":
    main()
