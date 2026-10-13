"""Longueur pivot-pointe courante (file_io/pivot_length.py)."""
import math

import pytest

from file_io import pivot_length as pl
from file_io.gcode_export import GcodeParams, machine_to_firmware


def test_default_then_last_recorded_value(tmp_path):
    path = str(tmp_path / "l.json")
    assert pl.current_pivot_to_tip(path)[0] == pl.DEFAULT_PIVOT_TO_TIP
    pl.record_pivot_length(46.0, "4.0 mm", path=path)
    pl.record_pivot_length(45.6, "3.6 mm", path=path)
    assert [r.l_mm for r in pl.load_history(path)] == [46.0, 45.6]
    assert pl.current_pivot_to_tip(path)[0] == pytest.approx(45.6)


def test_wrong_pivot_length_drifts_the_tip():
    # G-code calcule avec L faux de 0,4 mm, execute par une machine a L juste :
    # la pointe derive de dL sin(theta) en XY et dL (1 - cos theta) en Z.
    true_l, wrong_l, theta = 45.6, 46.0, 45.0
    tip = [60.0, 0.0, 40.0]
    p = GcodeParams()
    x, _, z, _ = machine_to_firmware(tip, theta, wrong_l, None, p.singularity_radius, p.tilt_sign, p.z_reference, p.b_sign)
    x_ok, _, z_ok, _ = machine_to_firmware(tip, theta, true_l, None, p.singularity_radius, p.tilt_sign, p.z_reference, p.b_sign)
    dl, t = wrong_l - true_l, math.radians(theta)
    assert abs(x - x_ok) == pytest.approx(dl * math.sin(t), abs=1e-9)
    assert abs(z - z_ok) == pytest.approx(dl * (1 - math.cos(t)), abs=1e-9)
