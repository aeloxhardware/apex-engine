"""
Equivalence de step5 natif (pas 13 + 14) contre Python sur un modele REEL :
meme champ scalaire, step5 rejoue deux fois (interrupteur force), couches et
contours compares au bit, temps mesures d'affilee, taux de repli du trace de
bord natif rapporte.

Usage : APEX_CPP_LIB_DIR=... .venv/Scripts/python.exe -u scripts/check_step5_native.py <modele.stl> [--resolution=0.8]
"""
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import main as apex_main  # noqa: E402
from core import _cpp  # noqa: E402
import pipeline.step5_layers as s5  # noqa: E402


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    model = args[0]
    params = apex_main.PipelineParams(resolution=float(flags.get("--resolution", 0.8)))
    captured = {}
    real_run = s5.run

    def spy(*a, **k):
        captured["args"], captured["kwargs"] = a, k
        return real_run(*a, **k)

    s5.run = spy
    apex_main.layers_run = spy if hasattr(apex_main, "layers_run") else None
    from scripts.export_slicing_preview import run_slicing_only
    import scripts.export_slicing_preview as esp
    for mod in (apex_main, esp):
        for name, val in list(vars(mod).items()):
            if val is real_run:
                setattr(mod, name, spy)
    _cpp.set_enabled(True)
    run_slicing_only(model, params)
    s5.run = real_run
    a, k = captured["args"], captured["kwargs"]

    out = {}
    for native in (False, True, False, True):
        _cpp.set_enabled(native)
        before = s5._NATIVE_OUTLINE_FALLBACKS
        t0 = time.time()
        layers = real_run(*a, **k)
        dt = time.time() - t0
        out.setdefault(native, []).append(dt)
        out[("layers", native)] = layers
        if native:
            out["fallbacks"] = s5._NATIVE_OUTLINE_FALLBACKS - before
    lp, ln = out[("layers", False)], out[("layers", True)]
    assert len(lp) == len(ln), (len(lp), len(ln))
    same = all(
        p.isovalue == n.isovalue and len(p.contour_3d) == len(n.contour_3d)
        and all(np.array_equal(x, y) for x, y in zip(p.contour_3d, n.contour_3d))
        and np.array_equal(p.local_frame[0], n.local_frame[0])
        for p, n in zip(lp, ln)
    )
    print(f"{os.path.basename(model)} : {len(lp)} couches, contours identiques au bit : {same}")
    print(f"  step5 Python {min(out[False]):.2f} s, natif {min(out[True]):.2f} s "
          f"(x{min(out[False]) / max(min(out[True]), 1e-9):.2f}) ; replis trimesh : {out['fallbacks']}")
    if not same:
        sys.exit(1)


if __name__ == "__main__":
    main()
