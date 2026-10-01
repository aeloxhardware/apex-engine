import numpy as np
import pytest

import main
from file_io.gcode_export import (GcodeParams, clearance_move, export_gcode, feed_plan, firmware_to_machine,
                                  head_axis_vectors, machine_to_firmware, parse_gcode)
from pipeline.step8_machine_toolpath import coupled_axis_angle
from pipeline.step9_collision_sim import _positioned_toolhead_vertices


@pytest.fixture(scope="module")
def machine():
    m = main._default_machine_profile()
    m.nozzle_model = None
    return m


def _c(machine, p):
    return coupled_axis_angle(p, machine.platform_origin, machine.rotation_axes[0].axis,
                              machine.radial_direction, 1e-6)


def test_pivot_offset_matches_step9_head_orientation(machine):
    """Le pivot est sur le CORPS de la tete, dans la direction que step9 donne
    a l'axe +Z du modele de tete : meme convention que la collision."""
    p = np.array([30.0, 40.0, 10.0])
    for b in (0.0, 30.0, 90.0, 150.0, 185.0):
        body = _positioned_toolhead_vertices(np.array([[0.0, 0.0, 45.0]]), p, [_c(machine, p), b],
                                             machine.rotation_axes)[0]
        fx, fy, fz, fb = machine_to_firmware(p, b, 45.0)
        r_body = float(np.hypot(body[0], body[1]))
        assert fx == pytest.approx(r_body, abs=1e-9)
        assert fz == pytest.approx(body[2] - 45.0, abs=1e-9)     # reference : pointe a B = 0
        assert fb == b


def test_known_cases(machine):
    fx, fy, fz, fb = machine_to_firmware(np.array([10.0, 0.0, 5.0]), 0.0, 45.0)
    assert (fx, fy, fz) == pytest.approx((10.0, 0.0, 5.0))
    fx, fy, fz, _ = machine_to_firmware(np.array([10.0, 0.0, 5.0]), 90.0, 45.0)
    assert (fx, fz) == pytest.approx((55.0, 5.0 - 45.0))
    _, fy, _, _ = machine_to_firmware(np.array([0.0, 10.0, 0.0]), 0.0, 45.0)
    assert fy == pytest.approx(90.0)


def test_firmware_y_is_c_minus_90(machine):
    rng = np.random.default_rng(0)
    for p in rng.uniform(-100, 100, (50, 3)):
        y = machine_to_firmware(p, 20.0, 45.0)[1]
        assert (y - (_c(machine, p) - 90.0) + 180.0) % 360.0 - 180.0 == pytest.approx(0.0, abs=1e-9)


def test_round_trip(machine):
    rng = np.random.default_rng(1)
    for p, b in zip(rng.uniform(-100, 100, (50, 3)), rng.uniform(0, 185, 50)):
        fw = machine_to_firmware(p, b, 45.0)
        assert firmware_to_machine(*fw, 45.0) == pytest.approx(p, abs=1e-9)


def test_y_unwrapped():
    prev = 179.0
    fw = machine_to_firmware(np.array([-10.0, -0.1, 0.0]), 0.0, 45.0, previous_y=prev)
    assert fw[1] == pytest.approx(180.0 + np.degrees(np.arctan2(0.1, 10.0)), abs=1e-6)


def test_clearance_single_axis():
    params = GcodeParams(travel_clearance=1.0)
    p = np.array([20.0, 0.0, 5.0])
    axis, q = clearance_move(p, 10.0, params)                 # tete quasi verticale : +Z
    assert axis == "Z" and q == pytest.approx(p + [0, 0, 1])
    axis, q = clearance_move(p, 175.0, params)                # tete retournee : -Z
    assert axis == "Z" and q == pytest.approx(p - [0, 0, 1])
    axis, q = clearance_move(p, 80.0, params)                 # tete couchee : R, vers le corps
    assert axis == "R" and q == pytest.approx(p + [1, 0, 0])
    p_low = np.array([20.0, 0.0, 1.5])                        # -Z interdit pres du plateau
    axis, q = clearance_move(p_low, 175.0, params, head_zlow=lambda b: 0.0)
    assert axis == "R"
    assert np.count_nonzero(q - p_low) == 1


def _toolpath(machine):
    pos = np.array([[20, 0, 0.2], [30, 0, 0.2], [30, 10, 0.2], [0, 30, 0.2], [0, 40, 0.2]], dtype=float)
    ang = np.array([[_c(machine, p), 0.0] for p in pos])
    return {"position": pos, "axis_angles": ang, "extrusion": np.array([0, 1, 1, 0, 1], dtype=float),
            "is_skin": np.array([0, 0, 0, 0, 1], dtype=bool)}


def test_export_extrusion_and_retraction(machine):
    tp = _toolpath(machine)
    params = GcodeParams()
    lines, rep = export_gcode(tp, machine, params)
    assert lines[:7][2:] == ["G21", "G90", "M82", "G92 E0", "M109 S195"]
    area = np.pi * (params.filament_diameter / 2) ** 2
    expected = (10 + 10 + 10) * 0.2 * 0.4 / area
    assert rep.filament_mm == pytest.approx(expected)
    assert rep.retractions == 1 and rep.travels == 2
    assert rep.c_mismatch_deg < 1e-9
    es = [float(w[1:]) for ln in lines for w in ln.split(";")[0].split() if w.startswith("E")]
    assert es[-1] == pytest.approx(expected - params.retract_length, abs=1e-4)   # 4 decimales
    assert all(b >= a - params.retract_length - 1e-9 for a, b in zip(es, es[1:]))


def test_parse_round_trip_positions(machine):
    tp = _toolpath(machine)
    params = GcodeParams(decimals=6)
    lines, _ = export_gcode(tp, machine, params)
    moves = parse_gcode(lines)
    tips = np.array([firmware_to_machine(x, y, z, b, 45.0) for x, y, z, b, _ in moves])
    # la pointe est au SOMMET du cordon : point de trajectoire (milieu de couche) + h/2 selon l'axe de tete
    half = np.where(tp["is_skin"], params.skin_thickness, params.layer_height) / 2
    for p, h in zip(tp["position"], half):
        assert np.min(np.linalg.norm(tips - (p + [0, 0, h]), axis=1)) < 1e-5
    extruded = [m for m in moves if m[4]]
    assert len(extruded) == 3


def test_feed_plan_slows_fast_b_changes():
    pts = np.column_stack([np.arange(5.0), np.zeros(5), np.zeros(5)])
    angles = np.column_stack([np.zeros(5), [0, 60, 0, 60, 0]])
    t_nom, t_real, factor = feed_plan(pts, angles, 750.0, 300.0, [0], 25.0)
    assert t_real > t_nom and np.all(factor <= 1.0 + 1e-12)
    angles[:, 1] = 0.0
    t_nom, t_real, factor = feed_plan(pts, angles, 750.0, 300.0, [0], 25.0)
    assert t_real == pytest.approx(t_nom)


def test_tip_at_bead_top_along_head_axis(machine):
    """Couche plate a z = h/2 (milieu de couche, step5) -> pointe emise a z = h ;
    tete inclinee -> decalage de h/2 le long de l'axe de la tete, pas en Z seul."""
    pos = np.array([[20.0, 0.0, 0.1], [30.0, 0.0, 0.1]])
    ang = np.array([[_c(machine, p), 0.0] for p in pos])
    tp = {"position": pos, "axis_angles": ang, "extrusion": np.array([0.0, 1.0])}
    lines, rep = export_gcode(tp, machine, GcodeParams(decimals=6))
    tips = np.array([firmware_to_machine(x, y, z, b, 45.0) for x, y, z, b, _ in parse_gcode(lines)])
    assert tips[:, 2].min() == pytest.approx(0.2, abs=1e-9)
    assert rep.tip_z_min == pytest.approx(0.2, abs=1e-9)
    lines, rep = export_gcode(tp, machine, GcodeParams(decimals=6, tip_at_bead_top=False))
    assert rep.tip_z_min == pytest.approx(0.1, abs=1e-9)
    v = head_axis_vectors(np.array([[10.0, 0.0, 5.0]]), 90.0)[0]
    assert v == pytest.approx([1.0, 0.0, 0.0], abs=1e-12)


def test_head_geometry_variants_round_trip():
    rng = np.random.default_rng(2)
    for sign in (1.0, -1.0):
        for ref in ("tip", "pivot"):
            for p, b in zip(rng.uniform(-100, 100, (20, 3)), rng.uniform(0, 185, 20)):
                fw = machine_to_firmware(p, b, 45.0, tilt_sign=sign, z_reference=ref)
                assert firmware_to_machine(*fw, 45.0, tilt_sign=sign, z_reference=ref) == pytest.approx(p, abs=1e-9)
    assert machine_to_firmware(np.array([10.0, 0.0, 5.0]), 0.0, 45.0, z_reference="pivot")[2] == pytest.approx(50.0)
    with pytest.raises(ValueError):
        machine_to_firmware(np.array([10.0, 0.0, 5.0]), 0.0, 45.0, z_reference="autre")


def test_guards_and_strict(machine, tmp_path):
    pos = np.array([[20.0, 0.0, 1.0], [30.0, 0.0, -0.5], [40.0, 0.0, 1.0]])
    ang = np.array([[_c(machine, p), 0.0] for p in pos])
    tp = {"position": pos, "axis_angles": ang, "extrusion": np.array([0.0, 1.0, 1.0])}
    out = tmp_path / "x.gcode"
    _, rep = export_gcode(tp, machine, GcodeParams(), str(out), head_zlow=lambda b: -2.0)
    assert rep.tip_below_plate == 1 and rep.tip_z_min == pytest.approx(-0.4)
    assert rep.head_z_min == pytest.approx(-2.4) and rep.head_below_margin >= 1
    assert rep.fw_x_range[0] <= 20.0 and rep.fw_x_range[1] >= 40.0 and rep.written and out.exists()
    out2 = tmp_path / "y.gcode"
    _, rep = export_gcode(tp, machine, GcodeParams(strict=True), str(out2))
    assert not rep.written and not out2.exists()
