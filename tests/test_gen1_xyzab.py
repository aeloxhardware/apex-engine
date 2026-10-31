"""Aetheris-4 Gen 1 : XYZ + A (autour de X) + B (autour de Y) portes par la tete."""
import numpy as np
import pytest

from file_io.xyzab_kinematics import HeadChain, body_axis, firmware_to_tip, tip_to_firmware
from file_io.gcode_export_xyzab import GcodeParamsXYZAB, export_gcode_xyzab, chain_for
from file_io import kinematics_description as kin


def _machine(order="AB"):
    import main
    m = main.aetheris4_gen1_profile(order=order)
    m.nozzle_model = None
    return m


def test_zero_angles_firmware_is_the_tip():
    ch = HeadChain()
    tip = np.array([[10.0, 20.0, 5.0]])
    assert np.allclose(tip_to_firmware(ch, tip, [[0.0, 0.0]]), tip)


def test_body_axis_follows_iso_conventions():
    ch = HeadChain()
    # A autour de X, main droite : +Z -> (0, -sin a, cos a) ; B autour de Y : +Z -> (sin b, 0, cos b)
    assert np.allclose(body_axis(ch, [[90.0, 0.0]])[0], [0, -1, 0], atol=1e-12)
    assert np.allclose(body_axis(ch, [[0.0, 90.0]])[0], [1, 0, 0], atol=1e-12)


@pytest.mark.parametrize("angles", [[30.0, 0.0], [0.0, -40.0], [25.0, 60.0], [-70.0, 35.0]])
def test_pivot_compensation_round_trips_and_matches_single_pivot_formula(angles):
    ch = HeadChain(pivot_to_tip=46.0)
    tip = np.array([[12.0, -3.0, 8.0]])
    fw = tip_to_firmware(ch, tip, [angles])
    assert np.allclose(firmware_to_tip(ch, fw, [angles]), tip)
    # centre commun : le point commande est la pointe + L * axe du corps - L * z
    d = body_axis(ch, [angles])
    assert np.allclose(fw, tip + 46.0 * d - 46.0 * np.array([[0, 0, 1.0]]))


def test_order_matters_and_both_orders_round_trip():
    a_first = HeadChain(axes=[(1, 0, 0), (0, 1, 0)], letters=["A", "B"])
    b_first = HeadChain(axes=[(0, 1, 0), (1, 0, 0)], letters=["B", "A"])
    # meme angle physique A=30, B=40 : l'ordre de montage change la direction
    assert not np.allclose(body_axis(a_first, [[30.0, 40.0]]), body_axis(b_first, [[40.0, 30.0]]), atol=1e-6)
    tip = np.array([[1.0, 2.0, 3.0]])
    for ch, ang in ((a_first, [[30.0, 40.0]]), (b_first, [[40.0, 30.0]])):
        assert np.allclose(firmware_to_tip(ch, tip_to_firmware(ch, tip, ang), ang), tip)


def test_distinct_offsets_are_honoured():
    ch = HeadChain(pivot_to_tip=40.0, offsets=[(0, 0, 0), (0, 0, -10.0)])
    tip = np.array([[0.0, 0.0, 20.0]])
    ang = [[90.0, 90.0]]
    fw = tip_to_firmware(ch, tip, ang)
    assert np.allclose(firmware_to_tip(ch, fw, ang), tip)
    assert not np.allclose(fw, tip_to_firmware(HeadChain(pivot_to_tip=40.0), tip, ang))


def test_description_is_valid_and_lists_ab_on_tool_side():
    m = _machine()
    desc = kin.apex_gantry_ab(GcodeParamsXYZAB(), m)
    kin.validate(desc)
    assert desc["chain"]["tool"] == ["X", "Y", "Z", "A", "B"]
    assert _machine("BA") and kin.apex_gantry_ab(GcodeParamsXYZAB(), _machine("BA"))["chain"]["tool"][-2:] == ["B", "A"]


def _toolpath():
    # trajet a vide, ligne extrudee dont la tete bascule de 0 a 40 deg (B) puis retour
    pos = np.array([[0, 0, 0.3], [10, 0, 0.3], [20, 0, 0.3], [30, 0, 0.3]], dtype=float)
    ang = np.array([[0, 0], [0, 0], [0, 40.0], [0, 40.0]], dtype=float)
    return {"position": pos, "axis_angles": ang, "extrusion": np.array([0, 1, 1, 1.0])}


def test_export_writes_six_axes_header_and_tip_follows_line():
    lines, rep = export_gcode_xyzab(_toolpath(), _machine(), GcodeParamsXYZAB())
    moves = [l for l in lines if l.startswith("G00")]
    assert all(" A" in l and " B" in l and " C0" in l for l in moves)
    assert lines[0].startswith(";") and any("apex-kinematics/1" in l for l in lines)
    assert rep.moves == len(moves) and rep.subdivided_moves >= 1     # le saut de B de 40 deg est decoupe
    assert rep.max_tip_deviation_before > 0.05
    # chaque mouvement relu : la pointe a mi-chemin reste sur la droite (tolerance du decoupage)
    ch = chain_for(_machine(), GcodeParamsXYZAB())
    rows = []
    for l in moves:
        if " E" not in l.split(";")[0]:
            continue
        d = {w[0]: float(w[1:]) for w in l.split(";")[0].split()[1:] if w[0] in "XYZAB"}
        rows.append((d["X"], d["Y"], d["Z"], d["A"], d["B"]))
    rows = np.array(rows)
    tips = firmware_to_tip(ch, rows[:, :3], rows[:, 3:5])
    # dernier mouvement extrude : derniere position + h/2 le long de l'axe du corps
    expected = _toolpath()["position"][-1] + body_axis(ch, [[0, 40.0]])[0] * 0.1
    assert np.allclose(tips[-1], expected, atol=1e-3)


def test_strict_refuses_tip_below_plate(tmp_path):
    tp = _toolpath()
    tp["position"][:, 2] = -1.0
    out = tmp_path / "x.gcode"
    _lines, rep = export_gcode_xyzab(tp, _machine(), GcodeParamsXYZAB(strict=True), str(out))
    assert rep.tip_below_plate > 0 and not rep.written and not out.exists()


def test_main_dispatches_head_ab_machine():
    import main
    assert main.machine_profile_for("aetheris-4").kinematics_kind == "polar-bed"
    g = main.machine_profile_for("aetheris-4-gen1")
    assert g.kinematics_kind == "head-ab" and g.position_coupled_axes == []
    with pytest.raises(ValueError):
        main.machine_profile_for("inconnue")


def test_gen1_pipeline_reaches_every_orientation_with_full_ranges():
    """l_shape sur Gen 1 (A, B +-180 deg) : le filtre median par axe est ignore, aucune orientation
    inatteignable (avec lui : 1 360, mesure du 2026-10-31) ; export sans pointe sous le plateau."""
    import os
    import main
    path = os.path.join(os.path.dirname(__file__), "..", "res", "test_models", "l_shape_overhang.stl")
    if not os.path.exists(path):
        pytest.skip("modele de test absent")
    p = main.PipelineParams()
    p.machine_profile = main.aetheris4_gen1_profile()      # defauts : A +-180, B +-90 (sphere sans ambiguite)
    p.machine_profile.nozzle_model = None
    assert p.machine_toolpath.free_axis_median_window > 1      # le defaut reste actif pour la polaire
    r = main.run_pipeline(path, p)
    assert r.step8_diagnostics == []
    ang = np.array([q.axis_angles for q in r.adjusted_toolpath.points])
    ext = np.array([q.extrusion for q in r.adjusted_toolpath.points]) > 0
    jumps = np.abs(np.diff(np.unwrap(np.radians(ang), axis=0, period=2 * np.pi), axis=0))[ext[1:]] * 180 / np.pi
    # Pas de saut de BRANCHE sur B (B sur +-180 : 186 deg au 99e centile, mesure). Les sauts de A pres de
    # B = +-90 (tete horizontale vers +-X, verrou de cardan : A indetermine) sont reels et ne sont pas testes ici.
    assert np.percentile(jumps[:, 1], 99) < 60.0
    lines, rep = main.export_gcode_checked(r.adjusted_toolpath, p, os.devnull)
    assert rep.tip_below_plate == 0 and rep.moves > 1000
