"""
step9 RAPIDE (arbre k-d, matiere densifiee) contre EXACT (capsules) sur un
ECHANTILLON de points : le controle exact de toute la trajectoire prend des
heures sur turkish/death, alors que l'exact en quelques points est bon marche
(`step10._collisions_at_points`, index exact reconstruit jusqu'au point).
Echantillon : N points signales par le rapide + N points echantillonnes par
step9 et non signales. Mesure : faux positifs, faux negatifs, ecart de
penetration.
Usage : APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/compare_step9_sampled.py <stl> [--n=1500]
"""
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from main import PipelineParams  # noqa: E402
from pipeline import step9_collision_sim as s9  # noqa: E402
from pipeline.step10_adjustment import _collisions_at_points  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    model = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    n = int(flags.get("--n", 1500))
    params = PipelineParams(resolution=0.8)
    tp = run_slicing_only(model, params).toolpath
    mp = params.machine_profile
    fast_p = s9.CollisionSimParams(**{**params.collision.__dict__, "fast_kdtree": True})
    exact_p = s9.CollisionSimParams(**{**params.collision.__dict__, "fast_kdtree": False})
    t0 = time.time()
    fast = s9.run(tp, mp, fast_p)
    t_fast = time.time() - t0
    fast_depth = {}
    for e in fast.events:
        key = (e.point_index, e.obstacle_type)
        fast_depth[key] = max(fast_depth.get(key, -1.0), e.penetration_depth)
    flagged = sorted({e.point_index for e in fast.events})
    sampled = s9._sample_point_indices(tp.points, params.collision.path_sampling_step)
    clean = sorted(set(sampled) - set(flagged))
    rng = np.random.default_rng(0)
    pick_f = sorted(rng.choice(flagged, min(n, len(flagged)), replace=False).tolist()) if flagged else []
    pick_c = sorted(rng.choice(clean, min(n, len(clean)), replace=False).tolist()) if clean else []
    t0 = time.time()
    exact = _collisions_at_points(tp, mp, exact_p, pick_f + pick_c)
    t_exact = time.time() - t0
    exact_depth = {}
    for e in exact.__iter__():
        key = (e.point_index, e.obstacle_type)
        exact_depth[key] = max(exact_depth.get(key, -1.0), e.penetration_depth)
    exact_pts = {k[0] for k in exact_depth}
    fp = [i for i in pick_f if i not in exact_pts]
    fn = [i for i in pick_c if i in exact_pts]
    common = set(fast_depth) & set(exact_depth)
    diffs = [exact_depth[k] - fast_depth[k] for k in common]
    print(f"{os.path.basename(model)} : rapide {len(flagged)} points signales sur {len(sampled)} echantillonnes ({t_fast:.0f} s)")
    print(f"  exact sur {len(pick_f)} signales + {len(pick_c)} non signales ({t_exact:.0f} s) : "
          f"faux positifs {len(fp)}, faux negatifs {len(fn)}")
    if diffs:
        print(f"  penetration exact - rapide : min {min(diffs):.3f} max {max(diffs):.3f} mm ; "
              f"faux positifs, profondeur rapide max {max([fast_depth.get((i, t), 0) for i in fp for t in ('printed_part', 'platform')] or [0]):.3f} mm ; "
              f"faux negatifs, profondeur exacte max {max([exact_depth.get((i, t), 0) for i in fn for t in ('printed_part', 'platform')] or [0]):.3f} mm")


if __name__ == "__main__":
    main()
