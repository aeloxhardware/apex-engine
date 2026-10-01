"""
Pas 15 : step6 avec parois/remplissage natifs (Clipper2) contre shapely, sur
les MEMES couches step5 d'un modele reel. Pas au bit : compare par couche le
nombre de polylignes, l'aire enclose par les parois et la longueur de
remplissage, et le temps de step6 dans les deux modes.
Usage : APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/check_native_fill.py <stl>
"""
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

import scripts.export_slicing_preview as esp  # noqa: E402
from main import PipelineParams  # noqa: E402
from pipeline import step6_toolpath2d as s6  # noqa: E402


def length(polys):
    return float(sum(np.linalg.norm(np.diff(p[:, :2], axis=0), axis=1).sum() for p in polys if len(p) > 1))


def main():
    captured = {}
    original = esp.toolpath2d_run

    def capture(layers, mesh_data, overhang_map, params):
        captured.update(layers=layers, args=(mesh_data, overhang_map, params))
        raise KeyboardInterrupt  # arret apres step5 : seules les couches comptent

    esp.toolpath2d_run = capture
    try:
        esp.run_slicing_only([a for a in sys.argv[1:] if not a.startswith('--')][0], PipelineParams(resolution=0.8))
    except KeyboardInterrupt:
        pass
    esp.toolpath2d_run = original
    layers, args = captured["layers"], captured["args"]
    if "--continuous" in sys.argv:
        import dataclasses
        args = (args[0], args[1], dataclasses.replace(args[2], continuous_rings=True))

    out = {}
    for native in (False, True):
        s6._NATIVE_FILL = native
        t0 = time.time()
        res = s6.run(layers, *args)
        out[native] = (res, time.time() - t0)
    (ref, t_ref), (nat, t_nat) = out[False], out[True]
    dl, dp = [], 0
    for a, b in zip(ref, nat):
        la, lb = length(a.toolpath_2d or []), length(b.toolpath_2d or [])
        if la > 0:
            dl.append((lb - la) / la)
        dp += abs(len(a.toolpath_2d or []) - len(b.toolpath_2d or []))
    dl = np.abs(np.array(dl))
    tot_a = sum(length(l.toolpath_2d or []) for l in ref)
    tot_b = sum(length(l.toolpath_2d or []) for l in nat)
    print(f"{os.path.basename([a for a in sys.argv[1:] if not a.startswith('--')][0])} : {len(layers)} couches ; step6 shapely {t_ref:.2f} s, natif {t_nat:.2f} s")
    print(f"  longueur totale {tot_a:.1f} / {tot_b:.1f} mm ({100 * (tot_b - tot_a) / tot_a:+.3f} %)")
    print(f"  ecart relatif de longueur par couche : median {np.median(dl):.2e}, max {dl.max():.2e}, "
          f"couches > 1 % : {(dl > 0.01).sum()}")
    print(f"  ecart cumule du nombre de polylignes : {dp}")


if __name__ == "__main__":
    main()
