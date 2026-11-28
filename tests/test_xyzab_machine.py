"""Machine XYZAB (limite #97) : texte du G-code identique en parallele, enveloppe de tete, liste des machines, borne de
hauteur a deux axes, appui le long de l'axe de la buse."""
from __future__ import annotations

import json

import numpy as np
import pytest

from config.machine_profiles import load_machine_profile, machines_main
from core.mesh_types import ToolpathMachine, ToolpathPoint
from element_ardor.step9_primitives import _positioned_toolhead_vertices
from file_io import gcode_xyzab as X


@pytest.fixture(scope="module")
def machine():
    return load_machine_profile("apex-xyzab")


def test_parallel_text_is_the_same_as_on_one_core(monkeypatch):
    n = 3000
    rng = np.random.default_rng(0)
    rows = {"fw": rng.uniform(-50, 50, (n, 3)), "ang": rng.uniform(-180, 180, (n, 2)), "feed": rng.uniform(1, 9e3, n),
            "e": np.cumsum(rng.uniform(0, 0.1, n)), "extrude": rng.random(n) > 0.3}
    comment = rng.integers(0, 4, n)
    one = X._text(rows, comment, 4)
    monkeypatch.setattr(X, "_PARALLEL_MIN_LINES", 10)
    monkeypatch.setattr(X, "_workers", lambda: 3)
    assert X._text(rows, comment, 4) == one
    assert len(one) == n


def test_the_envelope_stands_in_for_the_missing_head_model(machine):
    lo, hi = machine.nozzle_model.bounds
    assert lo[2] == pytest.approx(0.0) and hi[2] == pytest.approx(75.0)
    assert hi[0] == pytest.approx(25.0)


def test_machines_command_lists_the_profiles(capsys):
    assert machines_main([]) == 0
    out = json.loads(capsys.readouterr().out)
    kinds = {m["id"]: m["kinematics"] for m in out["machines"]}
    assert kinds["apex-xyzab"] == "head-ab" and kinds[out["default"]] == "polar"


def test_lowest_head_point_matches_the_placed_head(machine):
    from element_ardor.step10_reorient import _lowest_head_offset
    v = np.asarray(machine.nozzle_model.vertices)
    angles = np.random.default_rng(1).uniform(-180, 180, (20, 2))
    expected = [_positioned_toolhead_vertices(v, np.zeros(3), list(a), machine.rotation_axes)[:, 2].min() for a in angles]
    np.testing.assert_allclose(_lowest_head_offset(v, angles, machine.rotation_axes), expected, atol=1e-9)


def test_a_tilted_head_too_close_to_the_plate_is_straightened(machine):
    from element_ardor.step10_reorient import clamp_head_height
    from element_ardor.step9_collision_sim import CollisionSimParams
    p = ToolpathPoint(position=np.array([0.0, 0.0, 2.0]), axis_angles=[0.0, 80.0], extrusion=1.0, feedrate=1.0)
    sim = CollisionSimParams(path_sampling_step=2.0, safety_margin=1.0)
    out, moved, _dev = clamp_head_height(ToolpathMachine(points=[p]), machine, sim)
    assert moved == 1 and abs(out.points[0].axis_angles[1]) < 80.0


def test_support_is_found_along_the_tool_axis(machine):
    """Buse penchee de 60 deg : le cordon d'appui est a 1,2 mm le long de son axe, donc a 1,04 mm a l'horizontale,
    hors du rayon de recherche vertical (1 mm). Le test vertical le declare en l'air, le test le long de l'axe non."""
    from element_ardor.step9b_support_check import SupportCheckParams, run
    d = np.array([-np.sin(np.radians(60.0)), 0.0, np.cos(np.radians(60.0))])     # R_B(-60) (0, 0, 1)
    below = -1.2 * d
    def pt(x, y, z, e):
        return ToolpathPoint(position=np.array([x, y, z]), axis_angles=[0.0, -60.0], extrusion=e, feedrate=1.0)
    # un cordon n'est de la matiere que si ses deux extremites extrudent ; le trajet a vide (indice 2) les separe
    points = [pt(below[0], -1.0, 6.0 + below[2], 1.0), pt(below[0], 1.0, 6.0 + below[2], 1.0),
              pt(0.0, -0.8, 6.0, 0.0), pt(0.0, -0.5, 6.0, 1.0), pt(0.0, 0.5, 6.0, 1.0)]
    rep = run(ToolpathMachine(points=points), machine, 0.2, 0.0, SupportCheckParams())
    assert {3, 4} <= {e.point_index for e in rep.events}
    assert not {3, 4} & set(rep.tool_axis_unsupported) and 1 in rep.tool_axis_unsupported
