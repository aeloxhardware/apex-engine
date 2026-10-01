"""
Equivalence et gain d'une etape rejouee en SERIE puis en PARALLELE sur un
modele reel (memes entrees, capturees pendant un vrai pipeline).
Usage : .venv/Scripts/python.exe -u scripts/check_step_parallel.py <modele.stl> [--resolution=0.8]
"""
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import main as apex_main  # noqa: E402
import pipeline.step6_toolpath2d as s6  # noqa: E402
import scripts.export_slicing_preview as esp  # noqa: E402


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    model = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    params = apex_main.PipelineParams(resolution=float(flags.get("--resolution", 0.8)))
    captured = {}
    real = s6.run

    def spy(*a, **k):
        captured["a"], captured["k"] = a, k
        return real(*a, **k)

    for mod in (apex_main, esp):
        for name, val in list(vars(mod).items()):
            if val is real:
                setattr(mod, name, spy)
    esp.run_slicing_only(model, params)
    a, k = captured["a"], captured["k"]
    out = {}
    for threads in (1, None, 1, None):
        s6._STEP6_THREADS = threads
        t0 = time.time()
        layers = real(*a, **k)
        out.setdefault(threads, []).append(time.time() - t0)
        out[("l", threads)] = layers
    ls, lp = out[("l", 1)], out[("l", None)]
    same = len(ls) == len(lp) and all(
        len(x.toolpath_2d) == len(y.toolpath_2d) and all(np.array_equal(p, q) for p, q in zip(x.toolpath_2d, y.toolpath_2d))
        for x, y in zip(ls, lp))
    print(f"{os.path.basename(model)} step6 : serie {min(out[1]):.2f} s, parallele {min(out[None]):.2f} s, "
          f"identique au bit : {same}")
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()
