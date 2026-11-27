"""Machine XYZAB (docs/CINEMATIQUE_XYZAB.md) : profil, cinematique directe de la tete, G-code commandant le pivot."""
from __future__ import annotations

import numpy as np
import pytest

from config.machine_profiles import load_machine_profile
from file_io import gcode_xyzab as X
from file_io import kinematics_description as K
from file_io.gcode_export import GcodeParams, export_gcode


@pytest.fixture(scope="module")
def machine():
    return load_machine_profile("apex-xyzab", with_toolhead=False)


def test_profile_is_a_head_ab_machine(machine):
    assert machine.kinematics == "head-ab"
    assert machine.platform_size == (220.0, 220.0)
    assert machine.pivot_to_tip == 40.0
    assert machine.position_coupled_axes == []
    assert [tuple(a.axis) for a in machine.rotation_axes] == [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]


def test_head_direction_follows_a_then_b(machine):
    """A tourne autour de X et porte B (autour de Y au repos) : d = R_A R_B (0, 0, 1)."""
    d = X.head_directions(np.array([[0.0, 0.0], [90.0, 0.0], [0.0, 90.0], [90.0, 90.0]]), machine.rotation_axes)
    np.testing.assert_allclose(d, [[0, 0, 1], [0, -1, 0], [1, 0, 0], [1, 0, 0]], atol=1e-12)


def test_the_description_is_valid(machine):
    desc = X.description(machine, 40.0)
    K.validate(desc)
    assert desc["chain"] == {"bed": ["Y"], "tool": ["X", "Z", "A", "B"]}


def _words(line: str) -> dict:
    return {w[0]: float(w[1:]) for w in line.split(";")[0].split()[1:]}


def test_gcode_commands_the_pivot(machine):
    """Pointe a (1, 2, 0,2) au repos puis buse penchee de 90 deg vers +X : le G-code donne le pivot, a 40 mm de la
    pointe le long de la buse (pointe au sommet du cordon : + h/2 le long de la buse)."""
    toolpath = {"position": np.array([[1.0, 2.0, 0.2], [5.0, 2.0, 0.2]]),
                "axis_angles": np.array([[0.0, 0.0], [0.0, 90.0]]),
                "extrusion": np.array([0.0, 1.0])}
    params = GcodeParams(max_tip_deviation=None)
    lines, rep = export_gcode(toolpath, machine, params)
    moves = [_words(l) for l in lines if l.startswith("G01 X")]
    first = moves[1]                                  # apres l'approche
    h = params.layer_height / 2
    assert first["X"] == pytest.approx(1.0) and first["Y"] == pytest.approx(2.0)
    assert first["Z"] == pytest.approx(0.2 + h + 40.0)
    last_print = [m for m in moves if "E" in m][-1]
    assert last_print["B"] == pytest.approx(90.0) and last_print["A"] == pytest.approx(0.0)
    assert last_print["X"] == pytest.approx(5.0 + h + 40.0) and last_print["Z"] == pytest.approx(0.2)
    assert rep.tip_below_plate == 0 and rep.filament_mm > 0


def test_a_tilt_change_is_split_to_keep_the_tip_on_its_line(machine):
    toolpath = {"position": np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 1.0]]),
                "axis_angles": np.array([[0.0, 0.0], [0.0, 60.0]]),
                "extrusion": np.array([0.0, 1.0])}
    params = GcodeParams(max_tip_deviation=0.05, tip_at_bead_top=False)
    _lines, rep = export_gcode(toolpath, machine, params)
    assert rep.subdivided_moves >= 1
