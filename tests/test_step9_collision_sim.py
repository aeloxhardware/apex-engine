"""
tests/test_step9_collision_sim.py

Tests pour pipeline/step9_collision_sim.py. Fixtures principalement
synthetiques (trimesh.creation.*, convention du projet), sauf un test dedie
qui charge le vrai fichier res/3d/toolheads/aetheris-4/aetheris-4.glb via
io/mesh_loader.py - seule exception justifiee : le loader n'a de sens
qu'applique au fichier reel qu'il cible.

Execution : pytest tests/test_step9_collision_sim.py -v
"""

from __future__ import annotations

import os

import numpy as np
import pytest
import trimesh

from core.mesh_types import AABB, MachineProfile, RotationAxisSpec, ToolpathMachine, ToolpathPoint
from file_io.mesh_loader import load_toolhead_mesh
from pipeline.step9_collision_sim import (
    CollisionSimParams,
    _grow_printed_index_to,
    _point_to_segment_distance,
    _PrintedGeometryIndex,
    _positioned_toolhead_vertices,
    _sample_point_indices,
    run as step9_run,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _small_toolhead() -> trimesh.Trimesh:
    """Petite tete synthetique asymetrique (boite), origine locale a la pointe (bas de la boite)."""
    box = trimesh.creation.box(extents=[4.0, 6.0, 10.0])
    box.apply_translation([0.0, 0.0, 5.0])  # bas de la boite a z=0 -> pointe a l'origine
    return box


def _machine_no_rotation() -> MachineProfile:
    return MachineProfile(rotation_axes=[], platform_origin=None, nozzle_diameter=0.4)


def _straight_up_points(n: int, extrusion: float = 1.0) -> list[ToolpathPoint]:
    return [
        ToolpathPoint(
            position=np.array([0.0, 0.0, float(i) * 2.0]),
            axis_angles=[],
            extrusion=extrusion,
            feedrate=1000.0,
        )
        for i in range(n)
    ]


# ---------------------------------------------------------------------------
# Echantillonnage
# ---------------------------------------------------------------------------

def test_sample_point_indices_subsamples_with_coarse_step():
    points = _straight_up_points(20)
    dense = _sample_point_indices(points, path_sampling_step=0.001)
    coarse = _sample_point_indices(points, path_sampling_step=100.0)
    assert len(coarse) < len(dense)
    assert coarse[0] == 0
    assert coarse[-1] == len(points) - 1


def test_sample_point_indices_empty_points():
    assert _sample_point_indices([], path_sampling_step=1.0) == []


# ---------------------------------------------------------------------------
# Positionnement (N axes generique)
# ---------------------------------------------------------------------------

def test_positioned_toolhead_zero_axes_translates_only():
    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    position = np.array([5.0, 5.0, 5.0])
    result = _positioned_toolhead_vertices(vertices, position, [], [])
    assert np.allclose(result, vertices + position)


def test_positioned_toolhead_single_axis_rotation():
    vertices = np.array([[1.0, 0.0, 0.0]])
    axes = [RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0))]
    result = _positioned_toolhead_vertices(vertices, np.zeros(3), [90.0], axes)
    assert np.allclose(result, [[0.0, 1.0, 0.0]], atol=1e-9)


def test_positioned_toolhead_two_axes_matches_reversed_composition_order():
    """Meme convention que composed_rotation_angles_and_directions : axis_angles[-1] applique en premier."""
    vertices = np.array([[1.0, 0.0, 0.0]])
    axes = [
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 180.0)),
    ]
    # rotation B (axe X, 0deg = no-op ici car le vecteur est colineaire a X) puis C (axe Z, 90deg)
    result = _positioned_toolhead_vertices(vertices, np.zeros(3), [90.0, 0.0], axes)
    assert np.allclose(result, [[0.0, 1.0, 0.0]], atol=1e-9)


def test_positioned_toolhead_three_axes_stays_finite():
    vertices = np.random.default_rng(0).normal(size=(10, 3))
    axes = [
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 180.0)),
        RotationAxisSpec(axis=np.array([0.0, 1.0, 0.0]), angle_range_deg=(-30.0, 30.0)),
    ]
    result = _positioned_toolhead_vertices(vertices, np.array([1.0, 2.0, 3.0]), [10.0, 20.0, 5.0], axes)
    assert np.all(np.isfinite(result))


# ---------------------------------------------------------------------------
# Distance point-segment
# ---------------------------------------------------------------------------

def test_point_to_segment_distance_perpendicular():
    start = np.array([0.0, 0.0, 0.0])
    end = np.array([10.0, 0.0, 0.0])
    point = np.array([5.0, 3.0, 0.0])
    assert _point_to_segment_distance(point, start, end) == pytest.approx(3.0)


def test_point_to_segment_distance_beyond_endpoint_clamped():
    start = np.array([0.0, 0.0, 0.0])
    end = np.array([10.0, 0.0, 0.0])
    point = np.array([15.0, 0.0, 0.0])
    assert _point_to_segment_distance(point, start, end) == pytest.approx(5.0)


def test_point_to_segment_distance_degenerate_segment():
    start = np.array([1.0, 1.0, 1.0])
    end = np.array([1.0, 1.0, 1.0])
    point = np.array([1.0, 1.0, 4.0])
    assert _point_to_segment_distance(point, start, end) == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# Index spatial incremental - jamais en avance
# ---------------------------------------------------------------------------

def test_grow_printed_index_never_looks_ahead():
    """
    _grow_printed_index_to exclut la capsule qui se termine AU point teste
    (la buse est confondue avec le cordon qu'elle vient de deposer, cf.
    step9_collision_sim._grow_printed_index_to) - grown s'arrete donc un
    cran avant up_to_point_idx, pas exactement dessus.
    """
    points = _straight_up_points(10)
    index = _PrintedGeometryIndex(nozzle_diameter=1.0)

    last = _grow_printed_index_to(index, points, up_to_point_idx=3, last_grown_idx=0)
    assert last == 2
    # une capsule tres proche du point 8 (pas encore imprime) ne doit pas etre trouvee
    far_future_aabb_candidates = index.candidates_near(
        AABB(min_corner=points[8].position - 0.5, max_corner=points[8].position + 0.5),
        margin=0.5,
    )
    assert far_future_aabb_candidates == []


# ---------------------------------------------------------------------------
# run() - scenarios de bout en bout
# ---------------------------------------------------------------------------

def test_collision_against_printed_wall_detected():
    toolhead = _small_toolhead()
    machine = _machine_no_rotation()
    machine.nozzle_model = toolhead

    # trajectoire qui imprime un mur vertical (segment x=0,y=0, z:0->20), puis
    # revient au-dessus du mur imprime : le toolhead (boite locale x:+-2,
    # y:+-3, tip a l'origine) chevauche alors directement le mur en XY.
    points = [
        ToolpathPoint(position=np.array([0.0, 0.0, 0.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([0.0, 0.0, 20.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([0.0, 0.0, 10.0]), axis_angles=[], extrusion=0.0, feedrate=1000.0),
    ]
    toolpath = ToolpathMachine(points=points)
    params = CollisionSimParams(path_sampling_step=0.001, safety_margin=1.0)

    report = step9_run(toolpath, machine, params)

    assert report.nozzle_model_checked is True
    assert any(e.obstacle_type == "printed_part" for e in report.events)


def test_trajectory_clear_of_printed_wall_has_no_collision():
    toolhead = _small_toolhead()
    machine = _machine_no_rotation()
    machine.nozzle_model = toolhead

    points = [
        ToolpathPoint(position=np.array([0.0, 0.0, 0.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([0.0, 0.0, 20.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
        # loin du mur imprime
        ToolpathPoint(position=np.array([500.0, 500.0, 10.0]), axis_angles=[], extrusion=0.0, feedrate=1000.0),
    ]
    toolpath = ToolpathMachine(points=points)
    params = CollisionSimParams(path_sampling_step=0.001, safety_margin=1.0)

    report = step9_run(toolpath, machine, params)

    assert not report.has_collisions


def test_safety_margin_boundary_distinguishes_just_inside_vs_outside():
    toolhead = trimesh.creation.box(extents=[0.1, 0.1, 0.1])  # tete minuscule pour un controle precis de la distance
    machine = _machine_no_rotation()
    machine.nozzle_model = toolhead

    wall_points = [
        ToolpathPoint(position=np.array([0.0, 0.0, 0.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([0.0, 0.0, 10.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
    ]
    wall_radius = 0.05  # nozzle_diameter/2 avec nozzle_diameter=0.1

    safety_margin = 2.0
    just_inside = ToolpathPoint(
        position=np.array([wall_radius + safety_margin - 0.5, 0.0, 5.0]), axis_angles=[], extrusion=0.0, feedrate=1000.0
    )
    just_outside = ToolpathPoint(
        position=np.array([wall_radius + safety_margin + 0.5, 0.0, 5.0]), axis_angles=[], extrusion=0.0, feedrate=1000.0
    )

    machine.nozzle_diameter = 0.1
    params = CollisionSimParams(path_sampling_step=0.001, safety_margin=safety_margin)

    report_inside = step9_run(
        ToolpathMachine(points=[*wall_points, just_inside]), machine, params
    )
    report_outside = step9_run(
        ToolpathMachine(points=[*wall_points, just_outside]), machine, params
    )

    assert any(e.obstacle_type == "printed_part" for e in report_inside.events)
    assert not any(e.obstacle_type == "printed_part" for e in report_outside.events)


def test_missing_nozzle_model_returns_unchecked_report_no_exception():
    machine = _machine_no_rotation()
    machine.nozzle_model = None
    toolpath = ToolpathMachine(points=_straight_up_points(5))
    params = CollisionSimParams(path_sampling_step=1.0, safety_margin=1.0)

    report = step9_run(toolpath, machine, params)

    assert report.nozzle_model_checked is False
    assert report.events == []


def test_chuck_checked_is_always_false():
    machine = _machine_no_rotation()
    machine.nozzle_model = _small_toolhead()
    toolpath = ToolpathMachine(points=_straight_up_points(5))
    params = CollisionSimParams(path_sampling_step=1.0, safety_margin=1.0)

    report = step9_run(toolpath, machine, params)

    assert report.chuck_checked is False


def test_zero_rotation_axes_machine_runs_without_crashing():
    machine = _machine_no_rotation()
    machine.nozzle_model = _small_toolhead()
    toolpath = ToolpathMachine(points=_straight_up_points(8))
    params = CollisionSimParams(path_sampling_step=1.0, safety_margin=1.0)

    report = step9_run(toolpath, machine, params)
    assert report.nozzle_model_checked is True


def test_multi_axis_machine_runs_without_crashing():
    machine = MachineProfile(
        rotation_axes=[
            RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
            RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 180.0)),
        ],
        platform_origin=np.array([0.0, 0.0, 0.0]),
        nozzle_diameter=0.4,
    )
    machine.nozzle_model = _small_toolhead()
    points = [
        ToolpathPoint(
            position=np.array([0.0, 0.0, float(i) * 2.0]),
            axis_angles=[float(i * 10 % 360), float(i * 5 % 180)],
            extrusion=1.0,
            feedrate=1000.0,
        )
        for i in range(6)
    ]
    toolpath = ToolpathMachine(points=points)
    params = CollisionSimParams(path_sampling_step=0.001, safety_margin=1.0)

    report = step9_run(toolpath, machine, params)
    assert report.nozzle_model_checked is True


# ---------------------------------------------------------------------------
# Chargement du vrai GLB (seul test touchant un fichier externe)
# ---------------------------------------------------------------------------

def test_load_real_toolhead_glb_has_expected_extents():
    path = os.path.join(
        os.path.dirname(__file__), "..", "res", "3d", "toolheads", "aetheris-4", "aetheris-4.glb"
    )
    mesh = load_toolhead_mesh(path)

    extents = mesh.extents
    assert extents[0] == pytest.approx(76.0, abs=1.0)
    assert extents[1] == pytest.approx(77.15, abs=1.0)
    assert extents[2] == pytest.approx(35.0, abs=1.0)

    distances_to_origin = np.linalg.norm(mesh.vertices, axis=1)
    assert distances_to_origin.min() < 1.0  # la pointe (origine) doit etre quasi sur le maillage
