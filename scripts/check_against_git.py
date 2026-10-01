"""
Compare une etape du pipeline a sa version COMMITEE (git HEAD) sur un modele
reel : memes entrees, capturees pendant un vrai pipeline, sorties comparees
au bit. Sert a prouver qu'une optimisation Python ne change rien.
Usage : .venv/Scripts/python.exe -u scripts/check_against_git.py <modele.stl> <out/oldmods/stepX_old.py>
"""
import importlib.util
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
    model, old_path = sys.argv[1], sys.argv[2]
    spec = importlib.util.spec_from_file_location("step6_old", old_path)
    old = importlib.util.module_from_spec(spec)
    sys.modules["step6_old"] = old
    spec.loader.exec_module(old)
    # Les CLASSES (enums, dataclasses de parametres) doivent etre celles du
    # module courant : sinon `pattern is InfillPattern.X` compare deux enums
    # distincts et l'ancienne version prend silencieusement une autre branche
    # (piege rencontre le 2026-10-01, faux ecart de 33 contre 35 polylignes).
    for name, val in list(vars(old).items()):
        if isinstance(val, type) and getattr(val, "__module__", "") == "step6_old" and hasattr(s6, name):
            setattr(old, name, getattr(s6, name))
    params = apex_main.PipelineParams(resolution=0.8)
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
    t0 = time.time()
    new_layers = real(*a, **k)
    t_new = time.time() - t0
    t0 = time.time()
    old_layers = old.run(*a, **k)
    t_old = time.time() - t0
    same = len(new_layers) == len(old_layers) and all(
        len(x.toolpath_2d) == len(y.toolpath_2d) and all(np.array_equal(p, q) for p, q in zip(x.toolpath_2d, y.toolpath_2d))
        for x, y in zip(new_layers, old_layers))
    shown = 0
    for li, (x, y) in enumerate(zip(new_layers, old_layers)):
        if len(x.toolpath_2d) != len(y.toolpath_2d):
            print(f"  couche {li}: {len(x.toolpath_2d)} polylignes (actuel) contre {len(y.toolpath_2d)} (git)")
            shown += 1
        else:
            for pi, (p, q) in enumerate(zip(x.toolpath_2d, y.toolpath_2d)):
                if not np.array_equal(p, q):
                    d = np.abs(p - q).max() if p.shape == q.shape else None
                    print(f"  couche {li} polyligne {pi}: formes {p.shape} {q.shape}, ecart max {d}")
                    shown += 1
                    break
        if shown >= 6:
            break
    print(f"{os.path.basename(model)} : version git {t_old:.2f} s, version actuelle {t_new:.2f} s, identique au bit : {same}")
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()
