"""
Temps et evenements de step9 (collision avec la vraie tete) sur un modele,
mode rapide (arbre k-d) et, en option, mode exact sur un sous-echantillon.
Usage : APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/measure_step9.py <stl> [--resolution=0.8] [--exact-limit=0]
"""
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from main import PipelineParams  # noqa: E402
from pipeline import step9_collision_sim as s9  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    model = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    params = PipelineParams(resolution=float(flags.get("--resolution", 0.8)))
    if "--noradial" in sys.argv:
        params.reachable_tilt_from_machine = False
    if "--ppn" in sys.argv:
        params.machine_toolpath.per_point_normal = True
        params.machine_toolpath.per_point_normal_nodal = True
        params.machine_toolpath.free_axis_median_window = 9
    if "--median" in flags:
        # filtre median des axes libres (limite #63), applique dans step8 donc AVANT step9
        params.machine_toolpath.free_axis_median_window = int(flags["--median"])
    if "--no-ppn" in sys.argv:
        params.machine_toolpath.per_point_normal = False
    result = run_slicing_only(model, params)
    tp = result.toolpath
    p = s9.CollisionSimParams(**{**params.collision.__dict__, "fast_kdtree": True})
    t0 = time.time()
    rep = s9.run(tp, params.machine_profile, p)
    dt = time.time() - t0
    kinds = {}
    for e in rep.events:
        kinds[e.obstacle_type] = kinds.get(e.obstacle_type, 0) + 1
    n_sampled = len(s9._sample_point_indices(tp.points, p.path_sampling_step))
    print(f"{os.path.basename(model)} : {len(tp.points)} points, {n_sampled} echantillonnes ; "
          f"step9 rapide {dt:.1f} s ; evenements {kinds}")
    ev = [e for e in rep.events if e.obstacle_type == "printed_part"]
    if ev:
        b = np.array([tp.points[e.point_index].axis_angles[1] for e in ev])
        z = np.array([tp.points[e.point_index].position[2] for e in ev])
        pen = np.array([e.penetration_depth for e in ev])
        zmin = min(p_.position[2] for p_ in tp.points)
        allb = np.array([tp.points[i].axis_angles[1] for i in s9._sample_point_indices(tp.points, p.path_sampling_step)])
        print(f"  contacts piece : B median {np.median(b):.1f} (tous points : {np.median(allb):.1f}), "
              f"B>30 : {100*(b>30).mean():.0f} % (tous : {100*(allb>30).mean():.0f} %)")
        print(f"  hauteur au-dessus du bas de piece : mediane {np.median(z - zmin):.1f} mm ; "
              f"penetration mediane {np.median(pen):.3f}, >0,5 mm : {(pen > 0.5).sum()}")
    print(f"  penetration max { max((e.penetration_depth for e in rep.events), default=0):.3f} mm")


if __name__ == "__main__":
    main()
