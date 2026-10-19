import os, sys
ROOT = r"R:\APEX_NEW_VERSION\slicer_engine"
sys.path.insert(0, ROOT); os.chdir(ROOT); os.environ.setdefault("APEX_USE_CPP","1")
import numpy as np
from main import PipelineParams
from pipeline.step9b_support_check import SupportCheckParams
from scripts.export_slicing_preview import run_slicing_only
from scripts.measure_matrix import apply_variant
model, variant = sys.argv[1], sys.argv[2]
p = PipelineParams(resolution=0.8); p.support_check = SupportCheckParams(locate_supporting_point=True)
apply_variant(p, variant)
r = run_slicing_only(os.path.join(r"D:\apex_multibranch\models", model + ".stl"), p)
ev = r.support_check_report.events
print("events", len(ev), "tested", r.support_check_report.tested_extruded_points, "fields:", [k for k in vars(ev[0])] if ev else None)
if ev:
    pos = np.array([e.position for e in ev])
    fut = np.array([e.has_future_support for e in ev]); tilt = np.array([e.supported_by_tilt for e in ev])
    print("future support (ordering defects):", int(fut.sum()), "tilt resolved:", int(tilt.sum()), "true:", int((~fut & ~tilt).sum()))
    for lo,hi in [(i, i+4) for i in range(0, 60, 4)]:
        m=(pos[:,2]>=lo)&(pos[:,2]<hi)
        if m.any(): print(f"z[{lo},{hi}) n={int(m.sum())} y-range {pos[m,1].min():.1f}..{pos[m,1].max():.1f}")
    idx=np.array([e.point_index for e in ev]); print("point_index range", idx.min(), idx.max(), "of", len(r.toolpath.points))
if ev:
    gap = np.array([e.supporting_point_index - e.point_index if e.supporting_point_index is not None and e.supporting_point_index >= 0 else -1 for e in ev])
    g = gap[(gap > 0)]
    print("ordering gap (points) percentiles 10/50/90:", np.percentile(g, [10, 50, 90]) if len(g) else None, "of", len(r.toolpath.points))
    print("z hist of events:", np.histogram(pos[:, 2], bins=10)[0].tolist(), "z range", round(float(pos[:,2].min()),1), round(float(pos[:,2].max()),1))
    print("x,y centroid of events", np.round(pos[:, :2].mean(0), 1), "extents of model z", round(float(max(p.position[2] for p in r.toolpath.points)),1))
    # layers : portee des evenements en indice de couche
    import collections
    print("point_index decile hist:", np.histogram(idx, bins=10)[0].tolist())
if ev:
    ext = np.array([p.extrusion > 0 for p in r.toolpath.points]); run = np.cumsum(~ext)
    sup = np.array([e.supporting_point_index if e.supporting_point_index is not None else -1 for e in ev])
    fut = np.array([e.has_future_support for e in ev]); ok = fut & (sup >= 0)
    same_run = (run[idx[ok]] == run[sup[ok]])
    print("ordering defects:", int(fut.sum()), "localised:", int(ok.sum()), "support in SAME extruded run:", int(same_run.sum()), f"({100*same_run.mean():.0f}%)")
    far = ok.copy(); far[ok] = ~same_run
    print("gap percentiles of other-run defects:", np.percentile((sup - idx)[far], [10, 50, 90]) if far.any() else None)
