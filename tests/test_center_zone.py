"""Zone centrale de retractation a l'export (file_io/center_zone.py)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from core.mesh_types import ToolpathMachine, ToolpathPoint
from file_io.center_zone import apply_center_retract_zone, auto_center_radius
from file_io.gcode_export import GcodeParams, export_gcode, toolpath_arrays
from main import _default_machine_profile
from pipeline.step8_machine_toolpath import coupled_axis_angle


@pytest.fixture(scope="module")
def machine():
    return _default_machine_profile()


def _point(machine, x, y, z=1.0, extrusion=1.0):
    c = coupled_axis_angle(np.array([x, y, z]), np.asarray(machine.platform_origin, dtype=float),
                           np.asarray(machine.rotation_axes[0].axis, dtype=float),
                           np.asarray(machine.radial_direction, dtype=float), 1e-6)
    return ToolpathPoint(position=np.array([x, y, z], dtype=float), axis_angles=[c or 0.0, 0.0],
                         extrusion=extrusion, feedrate=1500.0)


def _chord(machine, x0=-10.0, x1=10.0, y=0.0):
    return ToolpathMachine(points=[_point(machine, x0, y, extrusion=0.0), _point(machine, x0, y),
                                   _point(machine, x1, y)])


def test_disabled_or_missing_radius_returns_the_toolpath_untouched(machine):
    tp = _chord(machine)
    for radius in (None, 0.0, -1.0):
        out, rep = apply_center_retract_zone(tp, machine, radius)
        assert out is tp and not rep.applied


def test_chord_through_the_centre_is_split_and_the_inside_becomes_a_travel(machine):
    tp = _chord(machine)
    out, rep = apply_center_retract_zone(tp, machine, 2.0)
    xs = [round(float(p.position[0]), 6) for p in out.points]
    ext = [p.extrusion > 0 for p in out.points]
    assert rep.applied and rep.crossings == 1
    assert rep.removed_length == pytest.approx(4.0)
    # depart (trajet), debut, entree, sortie (trajet), sortie (extrudee), fin
    assert xs == [-10.0, -10.0, -2.0, 2.0, 2.0, 10.0]
    assert ext == [False, True, True, False, True, True]
    assert len(tp.points) == 3, "la trajectoire d'entree n'est jamais mutee"


def test_chord_missing_the_zone_is_unchanged(machine):
    tp = _chord(machine, y=5.0)
    out, rep = apply_center_retract_zone(tp, machine, 2.0)
    assert rep.crossings == 0 and rep.removed_length == 0.0
    assert [tuple(p.position) for p in out.points] == [tuple(p.position) for p in tp.points]


def test_zone_follows_the_platform_axis_not_the_world_origin(machine):
    shifted = replace(machine, platform_origin=np.array([50.0, 0.0, 0.0]))
    out, rep = apply_center_retract_zone(_chord(shifted, x0=40.0, x1=60.0), shifted, 2.0)
    assert rep.crossings == 1
    assert sorted({round(float(p.position[0]), 6) for p in out.points}) == [40.0, 48.0, 52.0, 60.0]
    far, rep_far = apply_center_retract_zone(_chord(shifted, x0=-10.0, x1=10.0), shifted, 2.0)
    assert rep_far.crossings == 0


def test_coupled_angle_of_new_points_comes_from_their_position(machine):
    out, _ = apply_center_retract_zone(_chord(machine), machine, 2.0)
    for p in out.points[2:5]:
        expected = coupled_axis_angle(p.position, np.asarray(machine.platform_origin, dtype=float),
                                      np.asarray(machine.rotation_axes[0].axis, dtype=float),
                                      np.asarray(machine.radial_direction, dtype=float), 1e-6)
        assert p.axis_angles[0] == pytest.approx(expected, abs=1e-6)


def test_segment_ending_inside_the_zone_keeps_only_its_outer_part(machine):
    tp = ToolpathMachine(points=[_point(machine, 10.0, 0.0, extrusion=0.0), _point(machine, 10.0, 0.0),
                                 _point(machine, 0.0, 0.0)])
    out, rep = apply_center_retract_zone(tp, machine, 3.0)
    assert rep.removed_length == pytest.approx(3.0)
    assert out.points[-1].extrusion == 0.0 and out.points[-2].extrusion > 0.0
    assert out.points[-2].position[0] == pytest.approx(3.0)


def test_gcode_retracts_once_and_never_extrudes_inside_the_zone(machine):
    out, _ = apply_center_retract_zone(_chord(machine), machine, 2.0)
    arrays = toolpath_arrays(out)
    lines, report = export_gcode(arrays, machine, GcodeParams(pivot_to_tip=46.0), None)
    before, rep0 = export_gcode(toolpath_arrays(_chord(machine)), machine, GcodeParams(pivot_to_tip=46.0), None)
    assert report.retractions == rep0.retractions + 1
    import re
    radii = []
    for line in lines:
        if line.startswith("G00") and " E" in line:
            radii.append(float(re.search(r" X(-?[\d.]+)", line).group(1)))
    assert min(radii) >= 2.0 - 1e-3, "aucun cordon ne doit etre depose a moins de 2 mm de l'axe"


def test_auto_radius_is_print_speed_over_angular_speed():
    assert auto_center_radius(25.0, 600.0) == pytest.approx(25.0 / np.radians(600.0))
