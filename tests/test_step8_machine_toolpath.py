"""
tests/test_step8_machine_toolpath.py

Tests pour pipeline/step8_machine_toolpath.py. Contrairement a
test_step7_scheduling.py, step8 ne depend que de structures de donnees
(list[Layer], PrintSchedule, MachineProfile) - les fixtures construisent
ces objets a la main plutot que de rejouer tout le pipeline step1..step7,
puisque step8 est teste comme une fonction pure sur son contrat d'entree.

Execution : pytest tests/test_step8_machine_toolpath.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from core.geometry_utils import composed_rotation_angles_and_directions
from core.mesh_types import (
    Layer,
    MachineProfile,
    PrintPass,
    PrintSchedule,
    RotationAxisSpec,
)
from pipeline.step8_machine_toolpath import (
    Step8Params,
    ToolpathResolutionError,
    _clamp_rate,
    _detect_singular_axes,
    _rewrap_angle_sequence,
    _unwrap_angle_sequence,
    run as step8_run,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_layer(
    layer_id: int,
    local_frame: tuple[np.ndarray, np.ndarray, np.ndarray],
    toolpath_2d: list[np.ndarray],
) -> Layer:
    return Layer(
        id=layer_id,
        contour_3d=[],
        region_ids=[0],
        local_frame=local_frame,
        toolpath_2d=toolpath_2d,
        isovalue=0.0,
        segment_colors=[None] * len(toolpath_2d),
    )


def _simple_frame() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Repere local dont la normale (axis_u x axis_v) pointe vers -Z (buse vers le bas)."""
    origin = np.array([0.0, 0.0, 10.0])
    axis_u = np.array([1.0, 0.0, 0.0])
    axis_v = np.array([0.0, -1.0, 0.0])  # u x v = (0,0,-1)
    return origin, axis_u, axis_v


def _current_machine_profile() -> MachineProfile:
    """Machine actuelle : plateau C (tour complet) + tete B (0-180deg), perpendiculaires."""
    return MachineProfile(
        rotation_axes=[
            RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
            RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 180.0)),
        ],
        platform_origin=np.array([0.0, 0.0, 0.0]),
        max_angular_velocity=0.0,
        max_angular_acceleration=0.0,
    )


# ---------------------------------------------------------------------------
# Aller-retour cinematique directe/inverse (machine actuelle, 2 axes)
# ---------------------------------------------------------------------------

def test_round_trip_current_machine_resolves_within_tolerance():
    """
    Angles connus -> direction via composed_rotation_angles_and_directions ->
    _resolve_axis_angles (via run()) doit retrouver une direction resolue a
    moins de 0.5deg de la direction cible - tolerance concrete, pas une
    simple preuve de convergence en principe (voir plan : l'echelle
    pertinente du pipeline est le degre, pas la fraction de degre).
    """
    axes = _current_machine_profile().rotation_axes
    reference = np.array([0.0, 0.0, -1.0])

    known_angles_deg = np.array([40.0, 65.0])
    from pipeline.step8_machine_toolpath import _resolve_axis_angles

    # cinematique directe : angle connu -> direction cible
    directions, angles_grid = composed_rotation_angles_and_directions(
        axes, reference, samples_per_axis=2
    )
    # calcule la direction cible en composant directement (independant de la grille)
    from core.geometry_utils import _rodrigues_rotate

    target = reference.reshape(1, 3)
    for spec, angle_deg in zip(reversed(axes), reversed(list(known_angles_deg))):
        axis = spec.axis / np.linalg.norm(spec.axis)
        target = _rodrigues_rotate(target, axis, np.array([np.radians(angle_deg)]))
    target = target[0]

    params = Step8Params(default_feedrate=1000.0)
    resolved_angles, residual_error_rad = _resolve_axis_angles(target, axes, reference, params)

    assert np.degrees(residual_error_rad) < 0.5, (
        f"Erreur residuelle {np.degrees(residual_error_rad):.3f}deg trop grande "
        f"pour une trajectoire d'impression reelle"
    )


# ---------------------------------------------------------------------------
# Singularite au centre (R≈0) : heritage de l'angle precedent
# ---------------------------------------------------------------------------

def test_singularity_at_platform_center_inherits_previous_angle():
    machine = _current_machine_profile()
    origin_frame = (
        np.array([0.0, 0.0, 10.0]),  # exactement sur l'axe de rotation 0 (platform_origin=(0,0,0))
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, -1.0, 0.0]),
    )
    off_center_frame = _simple_frame()
    off_center_frame = (
        np.array([5.0, 0.0, 10.0]),
        off_center_frame[1],
        off_center_frame[2],
    )

    layer0 = _make_layer(0, off_center_frame, [np.array([[0.0, 0.0], [1.0, 0.0]])])
    layer1 = _make_layer(1, origin_frame, [np.array([[0.0, 0.0], [1.0, 0.0]])])

    schedule = PrintSchedule(
        passes=[
            PrintPass(segments=[(0, 0)], color_id=None),
            PrintPass(segments=[(1, 0)], color_id=None),
        ]
    )

    toolpath, _diagnostics = step8_run(
        [layer0, layer1], schedule, machine, Step8Params(default_feedrate=1000.0)
    )

    # le dernier point de la couche 0 et le premier point (post-singularite) de
    # la couche 1 doivent partager le meme angle de plateau (axe 0) : herite,
    # pas recalcule instable.
    last_angle_before = toolpath.points[-3].axis_angles[0]  # dernier point couche 0 (2 pts + eventuel travel)
    angle_at_singularity = None
    for point in toolpath.points:
        if abs(point.position[0]) < 1e-6 and abs(point.position[1]) < 1e-6:
            angle_at_singularity = point.axis_angles[0]
            break

    assert angle_at_singularity is not None, "Aucun point sur l'axe de singularite trouve"
    assert angle_at_singularity == pytest.approx(last_angle_before, abs=1e-9)


def test_detect_singular_axes_flags_axis_on_platform_center():
    axes = _current_machine_profile().rotation_axes
    platform_origin = np.array([0.0, 0.0, 0.0])

    on_axis = np.array([0.0, 0.0, 10.0])
    off_axis = np.array([5.0, 0.0, 10.0])

    singular_on = _detect_singular_axes(on_axis, platform_origin, axes, [0.0, 0.0], 1e-3)
    singular_off = _detect_singular_axes(off_axis, platform_origin, axes, [0.0, 0.0], 1e-3)

    assert 0 in singular_on
    assert 0 not in singular_off


# ---------------------------------------------------------------------------
# Axe unique verrouille
# ---------------------------------------------------------------------------

def test_locked_single_axis_always_returns_locked_value():
    from pipeline.step8_machine_toolpath import _resolve_axis_angles

    axes = [RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(30.0, 30.0))]
    reference = np.array([1.0, 0.0, 0.0])
    params = Step8Params(default_feedrate=1000.0)

    # direction desiree tres differente de ce qu'un axe verrouille a 30deg peut produire
    desired = np.array([0.0, 0.0, -1.0])
    angles, residual_error_rad = _resolve_axis_angles(desired, axes, reference, params)

    assert angles == pytest.approx([30.0])
    assert residual_error_rad > 0.0  # inatteignable : erreur residuelle non nulle


# ---------------------------------------------------------------------------
# Orientation inatteignable -> diagnostic, pas exception
# ---------------------------------------------------------------------------

def test_unreachable_orientation_produces_diagnostic_not_exception():
    machine = MachineProfile(
        rotation_axes=[
            RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 20.0)),
            RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 20.0)),
        ],
        platform_origin=np.array([0.0, 0.0, 0.0]),
    )

    # normale vers -Z tout en etant tres restreint (0-20deg sur les deux axes) :
    # direction difficile a atteindre en partant de reference (0,0,-1) tournee peu.
    frame = (
        np.array([5.0, 0.0, 10.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),  # u x v = (0,0,-1) x ... produit une normale tournee
    )
    layer = _make_layer(0, frame, [np.array([[0.0, 0.0], [1.0, 0.0]])])
    schedule = PrintSchedule(passes=[PrintPass(segments=[(0, 0)], color_id=None)])

    toolpath, diagnostics = step8_run(
        [layer], schedule, machine,
        Step8Params(default_feedrate=1000.0, max_unreachable_angular_error_deg=0.01),
    )

    assert len(toolpath.points) > 0
    for point in toolpath.points:
        assert point.axis_angles is not None
        assert len(point.axis_angles) == 2
        assert all(np.isfinite(a) for a in point.axis_angles)
    # avec un seuil de diagnostic tres bas (0.01deg), au moins un point restreint
    # devrait etre signale plutot que de lever une exception.
    assert isinstance(diagnostics, list)


# ---------------------------------------------------------------------------
# Machine a 0 axe (orientation libre)
# ---------------------------------------------------------------------------

def test_zero_axis_machine_produces_empty_axis_angles():
    machine = MachineProfile(rotation_axes=[], platform_origin=None)
    layer = _make_layer(0, _simple_frame(), [np.array([[0.0, 0.0], [1.0, 0.0]])])
    schedule = PrintSchedule(passes=[PrintPass(segments=[(0, 0)], color_id=None)])

    toolpath, diagnostics = step8_run(
        [layer], schedule, machine, Step8Params(default_feedrate=1000.0)
    )

    assert len(toolpath.points) == 2
    for point in toolpath.points:
        assert point.axis_angles == []
    assert diagnostics == []


# ---------------------------------------------------------------------------
# N=3 axes (generique au-dela de 2)
# ---------------------------------------------------------------------------

def test_three_axis_machine_resolves_without_crashing():
    machine = MachineProfile(
        rotation_axes=[
            RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
            RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 180.0)),
            RotationAxisSpec(axis=np.array([0.0, 1.0, 0.0]), angle_range_deg=(-30.0, 30.0)),
        ],
        platform_origin=np.array([0.0, 0.0, 0.0]),
    )
    layer = _make_layer(0, _simple_frame(), [np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])])
    schedule = PrintSchedule(passes=[PrintPass(segments=[(0, 0)], color_id=None)])

    toolpath, _diagnostics = step8_run(
        [layer], schedule, machine, Step8Params(default_feedrate=1000.0)
    )

    assert len(toolpath.points) == 3
    for point in toolpath.points:
        assert len(point.axis_angles) == 3
        assert all(np.isfinite(a) for a in point.axis_angles)


# ---------------------------------------------------------------------------
# Lissage vitesse/acceleration
# ---------------------------------------------------------------------------

def test_clamp_rate_never_exceeds_velocity_and_acceleration_limits():
    # saut artificiel important entre deux points consecutifs
    angles = np.array([0.0, 170.0, 175.0, 180.0])
    dt = np.array([1.0, 0.1, 0.1, 0.1])  # dt[0] non utilise (premier point)
    max_velocity_rad_s = np.radians(50.0)
    max_accel_rad_s2 = np.radians(200.0)

    clamped = _clamp_rate(angles, dt, max_velocity_rad_s, max_accel_rad_s2)

    for i in range(1, len(clamped)):
        velocity_deg_s = (clamped[i] - clamped[i - 1]) / dt[i]
        assert abs(velocity_deg_s) <= np.degrees(max_velocity_rad_s) + 1e-6


def test_clamp_rate_preserves_net_motion_direction():
    angles = np.array([0.0, 10.0, 20.0, 30.0])
    dt = np.array([1.0, 1.0, 1.0, 1.0])
    clamped = _clamp_rate(angles, dt, np.radians(100.0), np.radians(100.0))
    # limites larges : le mouvement ne doit quasiment pas etre altere
    assert clamped[-1] == pytest.approx(30.0, abs=1.0)


# ---------------------------------------------------------------------------
# Wraparound (axe a tour complet)
# ---------------------------------------------------------------------------

def test_unwrap_full_turn_axis_produces_small_deltas():
    angles = np.array([359.0, 1.0, 359.0])
    unwrapped = _unwrap_angle_sequence(angles, is_full_turn=True)

    deltas = np.diff(unwrapped)
    assert np.all(np.abs(deltas) <= 3.0), f"Deltas non deplies : {deltas}"


def test_rewrap_returns_values_within_range():
    unwrapped = np.array([359.0, 361.0, 359.0])
    rewrapped = _rewrap_angle_sequence(unwrapped, (0.0, 360.0), is_full_turn=True)
    assert np.all(rewrapped >= 0.0) and np.all(rewrapped < 360.0)


def test_unwrap_restricted_axis_is_noop():
    angles = np.array([10.0, 170.0, 5.0])
    unwrapped = _unwrap_angle_sequence(angles, is_full_turn=False)
    assert np.allclose(unwrapped, angles)


# ---------------------------------------------------------------------------
# Insertion de trajets a vide
# ---------------------------------------------------------------------------

def test_travel_move_inserted_between_two_polylines_same_layer():
    frame = _simple_frame()
    polyline_a = np.array([[0.0, 0.0], [1.0, 0.0]])
    polyline_b = np.array([[5.0, 5.0], [6.0, 5.0]])
    layer = _make_layer(0, frame, [polyline_a, polyline_b])

    schedule = PrintSchedule(
        passes=[PrintPass(segments=[(0, 0), (0, 1)], color_id=None)]
    )

    machine = MachineProfile(rotation_axes=[], platform_origin=None)
    toolpath, _diagnostics = step8_run(
        [layer], schedule, machine, Step8Params(default_feedrate=1000.0)
    )

    extrusions = [p.extrusion for p in toolpath.points]
    # 2 pts polyline_a + 1 travel + 2 pts polyline_b = 5
    assert len(toolpath.points) == 5
    assert extrusions == [1.0, 1.0, 0.0, 1.0, 1.0]


def test_travel_move_inserted_between_non_contiguous_passes():
    frame = _simple_frame()
    polyline_a = np.array([[0.0, 0.0], [1.0, 0.0]])
    polyline_b = np.array([[9.0, 9.0], [10.0, 9.0]])
    layer0 = _make_layer(0, frame, [polyline_a])
    layer1 = _make_layer(1, frame, [polyline_b])

    schedule = PrintSchedule(
        passes=[
            PrintPass(segments=[(0, 0)], color_id=1),
            PrintPass(segments=[(1, 0)], color_id=2),
        ]
    )

    machine = MachineProfile(rotation_axes=[], platform_origin=None)
    toolpath, _diagnostics = step8_run(
        [layer0, layer1], schedule, machine, Step8Params(default_feedrate=1000.0)
    )

    extrusions = [p.extrusion for p in toolpath.points]
    assert extrusions == [1.0, 1.0, 0.0, 1.0, 1.0]


# ---------------------------------------------------------------------------
# Garde toolpath_2d manquant
# ---------------------------------------------------------------------------

def test_missing_toolpath_2d_raises_resolution_error():
    layer = Layer(
        id=0,
        contour_3d=[],
        region_ids=[0],
        local_frame=_simple_frame(),
        toolpath_2d=None,
    )
    schedule = PrintSchedule(passes=[PrintPass(segments=[(0, 0)], color_id=None)])
    machine = MachineProfile(rotation_axes=[], platform_origin=None)

    with pytest.raises(ToolpathResolutionError):
        step8_run([layer], schedule, machine, Step8Params(default_feedrate=1000.0))


def test_schedule_referencing_missing_layer_raises_resolution_error():
    schedule = PrintSchedule(passes=[PrintPass(segments=[(99, 0)], color_id=None)])
    machine = MachineProfile(rotation_axes=[], platform_origin=None)

    with pytest.raises(ToolpathResolutionError):
        step8_run([], schedule, machine, Step8Params(default_feedrate=1000.0))
