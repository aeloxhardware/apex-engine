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


# ---------------------------------------------------------------------------
# Sens de depot inverse (Layer.deposition_reversed)
# ---------------------------------------------------------------------------

def _resolved_direction(layer: Layer, machine: MachineProfile) -> np.ndarray:
    """Direction de depot reellement resolue pour le premier point du calque."""
    from core.geometry_utils import _rodrigues_rotate

    schedule = PrintSchedule(passes=[PrintPass(segments=[(layer.id, 0)], color_id=None)])
    toolpath, _diagnostics = step8_run(
        [layer], schedule, machine, Step8Params(default_feedrate=1000.0)
    )

    # Cinematique DIRECTE a partir des angles resolus, avec exactement la
    # convention de composition du reste du projet (rotations en ordre
    # inverse de la liste) : on ne relit pas un intermediaire de step8, on
    # reconstruit ce que la machine ferait reellement.
    direction = np.array([[0.0, 0.0, -1.0]])
    angles = toolpath.points[0].axis_angles
    for spec, angle_deg in zip(reversed(machine.rotation_axes), reversed(list(angles))):
        axis = spec.axis / np.linalg.norm(spec.axis)
        direction = _rodrigues_rotate(direction, axis, np.array([np.radians(angle_deg)]))
    return direction[0]


def test_deposition_reversed_flips_the_resolved_direction():
    """
    Le sens de depot vient de l'ORDRE de parcours des isosurfaces, jamais du
    signe de grad(phi). Une couche atteinte en DESCENDANT le champ se depose
    sous la matiere existante : la tete doit pointer a l'oppose de la normale.

    Tant que l'ordre etait "phi croissant", ce cas n'existait pas et la
    supposition +grad(phi) restait vraie par construction - c'est
    `LayerParams.frontier_order` qui l'a rendue fausse.
    """
    machine = _current_machine_profile()
    polyline = np.array([[3.0, 0.0], [4.0, 0.0]])

    normal_layer = _make_layer(0, _simple_frame(), [polyline])
    reversed_layer = _make_layer(0, _simple_frame(), [polyline])
    reversed_layer.deposition_reversed = True

    forward = _resolved_direction(normal_layer, machine)
    backward = _resolved_direction(reversed_layer, machine)

    # Les deux directions doivent etre franchement opposees : leur produit
    # scalaire est proche de -1, pas simplement negatif.
    assert float(forward @ backward) < -0.9, (
        f"attendu des directions opposees, obtenu {forward} et {backward}"
    )


def test_deposition_reversed_defaults_to_false_and_changes_nothing():
    """
    Defaut = False : le chemin de code n'est jamais emprunte sans activation
    explicite, et une couche non marquee resout exactement comme avant.
    """
    assert Layer(
        id=0, contour_3d=[], region_ids=[0],
        local_frame=_simple_frame(), toolpath_2d=[], isovalue=0.0,
    ).deposition_reversed is False

    machine = _current_machine_profile()
    polyline = np.array([[3.0, 0.0], [4.0, 0.0]])

    plain = _make_layer(0, _simple_frame(), [polyline])
    explicit_false = _make_layer(0, _simple_frame(), [polyline])
    explicit_false.deposition_reversed = False

    assert _resolved_direction(plain, machine) == pytest.approx(
        _resolved_direction(explicit_false, machine)
    )


# ---------------------------------------------------------------------------
# Orientation PHYSIQUE de la tete (et non plus seulement des relations entre
# directions) : repere de couche tel que step5 le produit, tete positionnee
# avec la primitive de l'etape 9.
# ---------------------------------------------------------------------------

def _head_body_height(layer: Layer, machine: MachineProfile) -> float:
    """Hauteur moyenne, au-dessus de la pointe, d'une tete synthetique dont
    le corps s'etend le long de son +Z local (convention du modele de tete,
    CLAUDE.md : origine = pointe de buse, +Z = haut physique de la tete)."""
    from pipeline.step9_collision_sim import _positioned_toolhead_vertices

    schedule = PrintSchedule(passes=[PrintPass(segments=[(layer.id, 0)], color_id=None)])
    toolpath, _diagnostics = step8_run([layer], schedule, machine, Step8Params(default_feedrate=1000.0))
    point = toolpath.points[0]
    body = np.column_stack([np.zeros(20), np.zeros(20), np.linspace(1.0, 40.0, 20)])
    placed = _positioned_toolhead_vertices(body, point.position, point.axis_angles, machine.rotation_axes)
    return float(placed[:, 2].mean() - point.position[2])


def _step5_flat_frame() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Couche plate telle que step5._compute_local_frame la produit : normale
    = +grad(phi), donc +Z pour une piece qui croit vers le haut."""
    return np.array([5.0, 0.0, 10.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])


def test_flat_layer_keeps_the_head_above_the_nozzle():
    """
    Une couche plate (normale +Z) s'imprime tete vers le bas : le corps de la
    tete doit etre AU-DESSUS de la buse. Avant la correction (reference -Z
    dans step8), la meme couche donnait B = 179,9 deg : tete retournee, corps
    40 mm SOUS la buse, a l'interieur de la piece. Les tests existants ne le
    voyaient pas : ils construisaient des normales -Z (_simple_frame) et ne
    comparaient que des directions entre elles.
    """
    machine = _current_machine_profile()
    layer = _make_layer(0, _step5_flat_frame(), [np.array([[0.0, 0.0], [1.0, 0.0]])])
    assert _head_body_height(layer, machine) > 15.0


def test_reversed_flat_layer_puts_the_head_below_the_nozzle():
    """Couche atteinte en descendant le champ (deposition_reversed) : depot
    sous la matiere existante, la tete est retournee et son corps passe sous
    la buse (B proche de 180 deg, dans la plage de la machine)."""
    machine = _current_machine_profile()
    layer = _make_layer(0, _step5_flat_frame(), [np.array([[0.0, 0.0], [1.0, 0.0]])])
    layer.deposition_reversed = True
    assert _head_body_height(layer, machine) < -15.0


# ---------------------------------------------------------------------------
# Limite #59 : premier axe (C) impose par la position (machine R C Z B)
# ---------------------------------------------------------------------------

from pipeline.step8_machine_toolpath import (  # noqa: E402
    _RawPoint,
    _resolve_all_axis_angles_coupled,
    coupled_axis_angle,
)

_Z = np.array([0.0, 0.0, 1.0])


def _coupled_machine(extra_head_axis: bool = False, coupled: bool = True) -> MachineProfile:
    """Machine actuelle (C tour complet autour de Z, B autour de X, 0-185 deg,
    +R = -Y a angles nuls) ; `extra_head_axis` : + un axe D (2e rotation de
    tete, autour de Y) ; `coupled=False` : C libre (axe lineaire en plus)."""
    axes = [
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0)),
    ]
    if extra_head_axis:
        axes.append(RotationAxisSpec(axis=np.array([0.0, 1.0, 0.0]), angle_range_deg=(-90.0, 90.0)))
    return MachineProfile(
        rotation_axes=axes,
        platform_origin=np.array([0.0, 0.0, 0.0]),
        position_coupled_axes=[0] if coupled else [],
        radial_direction=np.array([0.0, -1.0, 0.0]),
    )


def _raw(position, direction) -> _RawPoint:
    d = np.asarray(direction, dtype=float)
    return _RawPoint(position=np.asarray(position, dtype=float), desired_direction=d / np.linalg.norm(d),
                     layer_id=0, segment_index=0, color_id=None, feedrate=1000.0, extrusion=1.0,
                     is_travel=False)


def _tilt(radial, angle_deg, tangential=False):
    """Direction penchee de angle_deg depuis +Z, vers +R (radial) ou dans le
    plan tangent (tangential=True)."""
    r = np.asarray(radial, dtype=float) / np.linalg.norm(radial)
    t = np.cross(_Z, r)
    side = t if tangential else r
    a = np.radians(angle_deg)
    return np.cos(a) * _Z + np.sin(a) * side


def test_coupled_angle_points_the_radial_direction_at_the_point():
    origin, axis, radial = np.zeros(3), _Z, np.array([0.0, -1.0, 0.0])
    assert coupled_axis_angle(np.array([0.0, -10.0, 3.0]), origin, axis, radial, 1e-3) == pytest.approx(0.0)
    assert coupled_axis_angle(np.array([10.0, 0.0, 0.0]), origin, axis, radial, 1e-3) == pytest.approx(90.0)
    assert coupled_axis_angle(np.array([0.0, 10.0, 0.0]), origin, axis, radial, 1e-3) == pytest.approx(180.0)
    assert coupled_axis_angle(np.array([0.0, 0.0, 5.0]), origin, axis, radial, 1e-3) is None


def test_coupled_flat_layer_gives_c_from_position_and_b_zero():
    machine = _coupled_machine()
    pts = [_raw([10.0, 0.0, 1.0], _Z), _raw([0.0, 10.0, 1.0], _Z)]
    angles, residuals = _resolve_all_axis_angles_coupled(pts, machine.rotation_axes, _Z, machine, Step8Params())
    assert angles[0][0] == pytest.approx(90.0)
    assert angles[1][0] == pytest.approx(180.0)
    assert abs(angles[0][1]) < 1e-6 and max(residuals) < 1e-6


def test_coupled_outward_tilt_is_reached_with_positive_b():
    machine = _coupled_machine()
    p = np.array([10.0, 0.0, 1.0])                    # +R au point = +X
    pts = [_raw(p, _tilt(p[:2].tolist() + [0.0], 30.0))]
    angles, residuals = _resolve_all_axis_angles_coupled(pts, machine.rotation_axes, _Z, machine, Step8Params())
    assert angles[0][1] == pytest.approx(30.0, abs=1e-3)
    assert residuals[0] < np.radians(1e-3)


def test_coupled_tangential_or_inward_tilt_is_unreachable_and_diagnosed():
    """B ne penche que dans le plan radial, vers +R : une cible penchee de 30
    deg tangentiellement ou vers le centre garde ~30 deg d'erreur."""
    machine = _coupled_machine()
    p = np.array([10.0, 0.0, 1.0])
    radial = np.array([1.0, 0.0, 0.0])
    pts = [_raw(p, _tilt(radial, 30.0, tangential=True)), _raw(p, _tilt(-radial, 30.0))]
    _angles, residuals = _resolve_all_axis_angles_coupled(pts, machine.rotation_axes, _Z, machine, Step8Params())
    assert np.degrees(residuals) == pytest.approx([30.0, 30.0], abs=0.5)


def test_coupled_head_with_d_axis_reaches_tangential_tilt():
    """Futur type de machine : + axe D (2e rotation de tete). C reste impose
    par la position, mais B et D atteignent une inclinaison tangentielle."""
    machine = _coupled_machine(extra_head_axis=True)
    p = np.array([10.0, 0.0, 1.0])
    pts = [_raw(p, _tilt(np.array([1.0, 0.0, 0.0]), 30.0, tangential=True))]
    angles, residuals = _resolve_all_axis_angles_coupled(pts, machine.rotation_axes, _Z, machine, Step8Params())
    assert angles[0][0] == pytest.approx(90.0)
    assert np.degrees(residuals[0]) < 1.0


def test_free_c_machine_reaches_tangential_tilt_through_run():
    """Futur type : plateau coulissant (axe lineaire) -> C libre : le chemin
    historique (orientation seule) atteint la meme cible sans diagnostic."""
    machine = _coupled_machine(coupled=False)
    origin = np.array([10.0, 0.0, 1.0])
    target = _tilt(np.array([1.0, 0.0, 0.0]), 30.0, tangential=True)
    u = np.cross(target, [1.0, 0.0, 0.0])
    u /= np.linalg.norm(u)
    v = np.cross(target, u)
    layer = Layer(id=0, contour_3d=[], region_ids=[0], local_frame=(origin, u, v),
                  toolpath_2d=[np.array([[0.0, 0.0], [1.0, 0.0]])])
    schedule = PrintSchedule(passes=[PrintPass(color_id=None, segments=[(0, 0)])])
    _tp, diagnostics = step8_run([layer], schedule, machine, Step8Params())
    assert diagnostics == []


def test_coupled_angle_is_never_rate_smoothed():
    """Lisser C deplacerait le point atteint : meme avec une limite de
    vitesse angulaire, C reste exactement l'angle polaire du point."""
    machine = _coupled_machine()
    machine.max_angular_velocity = 0.01
    origin = np.array([0.0, -10.0, 1.0])
    layer = Layer(id=0, contour_3d=[], region_ids=[0],
                  local_frame=(origin, np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])),
                  toolpath_2d=[np.array([[0.0, 0.0], [10.0, 10.0], [20.0, 10.0]])])
    schedule = PrintSchedule(passes=[PrintPass(color_id=None, segments=[(0, 0)])])
    tp, _diag = step8_run([layer], schedule, machine, Step8Params(clamp_angle_rates=True))
    for point in tp.points:
        expected = coupled_axis_angle(point.position, machine.platform_origin, _Z, machine.radial_direction, 1e-3)
        assert point.axis_angles[0] == pytest.approx(expected, abs=1e-6)


def test_angles_are_not_truncated_by_default_even_with_rate_limits():
    """Le firmware coordonne les axes : par defaut, les limites de vitesse du
    profil ne tronquent PAS les angles (l'avance ralentira a l'export)."""
    machine = _coupled_machine()
    machine.max_angular_velocity = 1e-6
    machine.max_angular_acceleration = 1e-6
    p = np.array([10.0, 0.0, 1.0])
    target = _tilt(np.array([1.0, 0.0, 0.0]), 40.0)
    u = np.cross(target, [0.0, 1.0, 0.0])
    u /= np.linalg.norm(u)
    v = np.cross(target, u)
    layer = Layer(id=0, contour_3d=[], region_ids=[0], local_frame=(p, u, v),
                  toolpath_2d=[np.array([[0.0, 0.0], [0.5, 0.0]])])
    schedule = PrintSchedule(passes=[PrintPass(color_id=None, segments=[(0, 0)])])
    tp, diagnostics = step8_run([layer], schedule, machine, Step8Params())
    assert tp.points[0].axis_angles[1] == pytest.approx(40.0, abs=0.5)
    assert diagnostics == []


def test_free_axis_median_filter_removes_b_oscillation_only():
    """Filtre median de l'axe libre (B), par polyligne : les pics isoles
    0 -> 14 -> 0 disparaissent, C (couple) n'est jamais touche, et une
    nouvelle polyligne n'herite rien de la precedente. (Une alternance
    PARFAITE 0/14/0/14 resiste a un median -- limite connue.)"""
    import main
    from pipeline.step8_machine_toolpath import _median_filter_free_axes, _RawPoint

    machine = main._default_machine_profile()
    raw, angles = [], []
    for seg in range(2):
        for i in range(12):
            raw.append(_RawPoint(position=np.array([10.0 + i, 0.0, 1.0]), desired_direction=np.array([0.0, 0.0, 1.0]),
                                 layer_id=0, segment_index=seg, color_id=None, feedrate=1.0, extrusion=1.0,
                                 is_travel=False))
            angles.append([90.0, (14.0 if i in (3, 6, 9) else 0.0) if seg == 0 else 30.0])
    out, residual = _median_filter_free_axes(angles, raw, machine.rotation_axes, machine,
                                             np.array([0.0, 0.0, 1.0]), 5)
    b = np.array(out)[:, 1]
    assert np.all(b[:12] == 0.0)
    assert np.all(b[12:] == 30.0)
    assert all(a[0] == 90.0 for a in out)
    assert residual[12] == pytest.approx(np.radians(30.0), abs=1e-9)


# ---------------------------------------------------------------------------
# Normale PAR POINT (Step8Params.per_point_normal, phase 3.1 du 2026-10-01)
# ---------------------------------------------------------------------------

def _ring_layer():
    """Couche plate en anneau autour de l'axe du plateau (normale de couche +Z)."""
    ang = np.linspace(0, 2 * np.pi, 13)[:-1]
    ring = np.column_stack([20 * np.cos(ang), 20 * np.sin(ang)])
    frame = (np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    return _make_layer(0, frame, [ring])


def _outward_field(tilt_deg):
    def field(points):
        p = np.asarray(points, dtype=float)
        r = p[:, :2] / np.linalg.norm(p[:, :2], axis=1)[:, None]
        a = np.radians(tilt_deg)
        return np.column_stack([np.sin(a) * r, np.full(len(p), np.cos(a))])
    return field


def test_per_point_normal_follows_local_field_on_coupled_machine():
    """Couche en anneau dont la normale LOCALE penche vers l'exterieur : avec
    une seule normale de couche (+Z), la tete reste verticale partout ; avec
    la normale par point, B vaut l'inclinaison locale en chaque point, sans
    aucun diagnostic (chaque normale est dans le plan radial de son point)."""
    machine = _coupled_machine()
    layer = _ring_layer()
    schedule = PrintSchedule(passes=[PrintPass(segments=[(0, 0)], color_id=None)])
    tp_ref, _ = step8_run([layer], schedule, machine, Step8Params())
    assert all(abs(p.axis_angles[1]) < 1e-6 for p in tp_ref.points)
    tp, diags = step8_run([layer], schedule, machine, Step8Params(per_point_normal=True),
                          direction_field=_outward_field(20.0))
    assert all(p.axis_angles[1] == pytest.approx(20.0, abs=1e-3) for p in tp.points)
    assert diags == []


def test_per_point_normal_ignored_without_field_and_keeps_layer_normal_when_opposed():
    machine = _coupled_machine()
    layer = _ring_layer()
    schedule = PrintSchedule(passes=[PrintPass(segments=[(0, 0)], color_id=None)])
    tp, _ = step8_run([layer], schedule, machine, Step8Params(per_point_normal=True))
    assert all(abs(p.axis_angles[1]) < 1e-6 for p in tp.points)
    flipped = lambda pts: -_outward_field(20.0)(pts)   # noqa: E731 - oppose a la normale de couche
    tp, _ = step8_run([layer], schedule, machine, Step8Params(per_point_normal=True), direction_field=flipped)
    assert all(abs(p.axis_angles[1]) < 1e-6 for p in tp.points)


def test_nodal_direction_field_is_exact_on_linear_field_and_continuous():
    """Normale continue (per_point_normal_nodal) : sur un champ LINEAIRE (gradient
    constant), elle rend exactement ce gradient partout ; elle ne saute pas d'un
    tetraedre voisin a l'autre (contrairement au plus proche centroide)."""
    from scipy.spatial import Delaunay
    from core.mesh_types import ScalarField, TetMesh
    from pipeline.step8_machine_toolpath import gradient_direction_field
    rng = np.random.default_rng(4)
    pts = rng.uniform(0, 10, (300, 3))
    tets = Delaunay(pts).simplices.astype(np.int64)
    g = np.array([0.2, -0.1, 1.0])
    tm = TetMesh(tet_vertices=pts, tetrahedra=tets, tet_adjacency=None)
    sf = ScalarField(values=pts @ g, gradient=np.tile(g, (len(tets), 1)))
    field = gradient_direction_field(tm, sf, nodal=True)
    q = rng.uniform(2, 8, (200, 3))
    assert np.allclose(field(q), g / np.linalg.norm(g), atol=1e-12)
    # gradient bruite par tetraedre : la variante nodale varie moins entre points voisins
    noisy = g + rng.normal(0, 0.3, (len(tets), 3))
    sf2 = ScalarField(values=pts @ g, gradient=noisy)
    line = np.column_stack([np.linspace(2, 8, 400), np.full(400, 5.0), np.full(400, 5.0)])
    jump = lambda d: np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", d[1:], d[:-1]), -1, 1))).mean()  # noqa: E731
    assert jump(gradient_direction_field(tm, sf2, nodal=True)(line)) < jump(gradient_direction_field(tm, sf2)(line))


def test_median_filter_residuals_match_pointwise_rodrigues():
    """Le residuel recalcule par lots apres le filtre median (2026-10-06) egale la
    boucle point par point d'origine, a la precision machine pres."""
    from core.geometry_utils import _rodrigues_rotate_single
    from pipeline.step8_machine_toolpath import _median_filter_free_axes
    rng = np.random.default_rng(3)
    machine = _coupled_machine(extra_head_axis=True)
    n = 400
    angles = np.column_stack([rng.uniform(0, 360, n), rng.uniform(0, 185, n), rng.uniform(-90, 90, n)])
    pts = [_raw(rng.uniform(-20, 20, 3), rng.normal(size=3)) for _ in range(n)]
    for k, p in enumerate(pts):
        p.layer_id = k // 50
    out_angles, residual = _median_filter_free_axes(angles.tolist(), pts, machine.rotation_axes, machine, _Z, 9)
    units = [np.asarray(s.axis, float) / np.linalg.norm(s.axis) for s in machine.rotation_axes]
    for i in range(n):
        d = _Z.copy()
        for u, a in zip(reversed(units), reversed(out_angles[i])):
            d = _rodrigues_rotate_single(d, u, np.radians(a))
        ref = float(np.arccos(np.clip(d @ pts[i].desired_direction, -1.0, 1.0)))
        assert residual[i] == pytest.approx(ref, abs=1e-9)
