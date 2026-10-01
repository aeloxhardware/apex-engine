"""
step8 en tableaux (`_RawPoints`, 2026-10-01) contre la version a objets
(out/oldmods/step8_before_arrays.py), sur les MEMES couches et le meme
ordonnancement : trajectoire et diagnostics identiques au bit, et temps.
Usage : APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/check_step8_arrays.py <stl> [--ppn] [--median9]
"""
import importlib.util
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

import scripts.export_slicing_preview as esp  # noqa: E402
from main import PipelineParams  # noqa: E402
from pipeline import step8_machine_toolpath as new8  # noqa: E402


def load_old():
    spec = importlib.util.spec_from_file_location("step8_old", os.path.join(ROOT, "out", "oldmods", "step8_before_arrays.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["step8_old"] = mod
    spec.loader.exec_module(mod)
    mod.Step8Params = new8.Step8Params      # memes classes que le module courant
    return mod


def main():
    model = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    params = PipelineParams(resolution=0.8)
    if "--median9" in sys.argv:
        params.machine_toolpath.free_axis_median_window = 9
    cap = {}

    def capture(*args, **kwargs):
        cap.update(args=args, kwargs=kwargs)
        return new8.run(*args, **kwargs)

    esp.machine_toolpath_run = capture
    if "--ppn" in sys.argv:
        params.machine_toolpath.per_point_normal = True
    esp.run_slicing_only(model, params)
    old8 = load_old()
    out = {}
    for name, mod in (("objets", old8), ("tableaux", new8)):
        t0 = time.time()
        out[name] = (mod.run(*cap["args"], **cap["kwargs"]), time.time() - t0)
    (ta, da), t_a = out["objets"]
    (tb, db), t_b = out["tableaux"]
    assert len(ta.points) == len(tb.points)
    for pa, pb in zip(ta.points, tb.points):
        assert np.array_equal(pa.position, pb.position)
        assert pa.axis_angles == pb.axis_angles
        assert pa.extrusion == pb.extrusion and pa.feedrate == pb.feedrate and pa.color_id == pb.color_id
    assert [(d.point_index, d.layer_id, d.segment_index, d.residual_angle_error_deg) for d in da] == \
           [(d.point_index, d.layer_id, d.segment_index, d.residual_angle_error_deg) for d in db]
    print(f"{os.path.basename(model)} : {len(ta.points)} points, {len(da)} diagnostics IDENTIQUES au bit ; "
          f"step8 objets {t_a:.2f} s, tableaux {t_b:.2f} s")


if __name__ == "__main__":
    main()
