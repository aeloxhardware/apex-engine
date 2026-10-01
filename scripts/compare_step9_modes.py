"""
step9 EXACT (capsules) contre RAPIDE (arbre k-d, matiere densifiee) sur un
modele : memes points touches ? Ecart de penetration ? Le rapide surestime la
distance d'au plus un demi-pas de densification (limite #62).
Usage : APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/compare_step9_modes.py <stl> [--noradial]
"""
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from main import PipelineParams  # noqa: E402
from pipeline import step9_collision_sim as s9  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402


def main():
    model = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    params = PipelineParams(resolution=0.8)
    if "--noradial" in sys.argv:
        params.reachable_tilt_from_machine = False
    tp = run_slicing_only(model, params).toolpath
    out = {}
    for fast in (True, False):
        p = s9.CollisionSimParams(**{**params.collision.__dict__, "fast_kdtree": fast})
        t0 = time.time()
        rep = s9.run(tp, params.machine_profile, p)
        out[fast] = ({(e.point_index, e.obstacle_type): e.penetration_depth for e in rep.events}, time.time() - t0)
        print(f"{'rapide' if fast else 'exact '} : {len(out[fast][0])} evenements en {out[fast][1]:.0f} s", flush=True)
    f, e = out[True][0], out[False][0]
    both = set(f) & set(e)
    print(f"communs {len(both)}, rapide seul {len(set(f) - set(e))}, exact seul {len(set(e) - set(f))}")
    if both:
        diffs = [e[k] - f[k] for k in both]
        print(f"penetration exact - rapide : min {min(diffs):.3f} max {max(diffs):.3f} mm")
    only_e = sorted(set(e) - set(f), key=lambda k: -e[k])[:5]
    print("exact seul, plus profondes :", [(k, round(e[k], 3)) for k in only_e])


if __name__ == "__main__":
    main()
