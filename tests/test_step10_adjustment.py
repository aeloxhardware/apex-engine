"""
tests/test_step10_adjustment.py

Tests pour pipeline/step10_adjustment.py, sur des trajectoires et des
modeles de tete synthetiques (aucun fichier externe). Les CollisionReport
d'entree sont produits par la VRAIE etape 9 plutot que fabriques a la main,
pour que l'etape 10 soit testee sur les structures qu'elle recevra en
production.

Execution : pytest tests/test_step10_adjustment.py -v
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from core.mesh_types import (
    MachineProfile,
    ResegmentationRequest,
    RotationAxisSpec,
    Step8Diagnostic,
    ToolpathMachine,
    ToolpathPoint,
)
from pipeline.step9_collision_sim import CollisionSimParams, run as collision_run
from pipeline.step10_adjustment import (
    Step10Params,
    _build_reorientation_context,
    _collisions_at_points,
    _direction_for_angles,
    _try_local_reorientation,
    _violates_rate_limits,
    run as adjustment_run,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _z_axis_machine(nozzle_model: trimesh.Trimesh) -> MachineProfile:
    """
    Machine a 1 axe de rotation (plateau Z, tour complet). Choisie pour les
    tests de reorientation : une rotation autour de Z ne change PAS la
    direction de depot (0,0,-1), donc le seuil de surplomb est trivialement
    respecte et le test isole la seule question qui nous interesse ici -
    la reorientation elimine-t-elle la collision.
    """
    return MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0))],
        platform_origin=None,  # pas de disque plateau : isole la collision piece imprimee
        nozzle_model=nozzle_model,
        nozzle_diameter=0.4,
    )


def _offset_toolhead() -> trimesh.Trimesh:
    """Tete ASYMETRIQUE (decalee en +X) : une rotation plateau la deplace reellement."""
    box = trimesh.creation.box(extents=[4.0, 4.0, 4.0])
    box.apply_translation([5.0, 0.0, 2.0])  # occupe x 3..7
    return box


def _axis_straddling_toolhead() -> trimesh.Trimesh:
    """
    Tete FINE centree sur l'axe de rotation (pointe a l'origine). Combinee a
    un mur porte par ce meme axe (_wall_through_axis_toolpath), la collision
    devient INVARIANTE par rotation autour de Z : la distance de chaque point
    de test a l'axe est conservee par la rotation, donc aucune reorientation
    ne peut l'eliminer. C'est le cas qui force les niveaux 2 puis 3.

    Une tete large decalee ne conviendrait pas : le mur traverserait
    l'interieur d'une face sans passer pres d'une arete, et une rotation
    suffirait a l'ecarter (cf. echantillonnage aretes/sommets, limite #21).
    """
    box = trimesh.creation.box(extents=[1.0, 1.0, 4.0])
    box.apply_translation([0.0, 0.0, 2.0])
    return box


def _wall_then_probe_toolpath(probe_position: np.ndarray) -> ToolpathMachine:
    """
    Imprime un mur vertical en (5,0,z) puis amene la tete au point sonde
    (non extrudant). Seul ce dernier point peut collisionner : au point 1 la
    capsule qui s'y termine est exclue (cordon en cours de depot, cf. etape 9).
    """
    return ToolpathMachine(
        points=[
            ToolpathPoint(
                position=np.array([5.0, 0.0, 0.0]), axis_angles=[0.0], extrusion=1.0, feedrate=1000.0
            ),
            ToolpathPoint(
                position=np.array([5.0, 0.0, 10.0]), axis_angles=[0.0], extrusion=1.0, feedrate=1000.0
            ),
            ToolpathPoint(
                position=probe_position, axis_angles=[0.0], extrusion=0.0, feedrate=1000.0
            ),
        ]
    )


def _wall_through_axis_toolpath(probe_position: np.ndarray) -> ToolpathMachine:
    """
    Imprime un mur PORTE PAR L'AXE DE ROTATION (x=0, y=0), puis amene la tete
    au point sonde. Avec _axis_straddling_toolhead, la collision resultante ne
    peut etre resolue par aucune rotation autour de Z.
    """
    return ToolpathMachine(
        points=[
            ToolpathPoint(
                position=np.array([0.0, 0.0, 0.0]), axis_angles=[0.0], extrusion=1.0, feedrate=1000.0
            ),
            ToolpathPoint(
                position=np.array([0.0, 0.0, 10.0]), axis_angles=[0.0], extrusion=1.0, feedrate=1000.0
            ),
            ToolpathPoint(
                position=probe_position, axis_angles=[0.0], extrusion=0.0, feedrate=1000.0
            ),
        ]
    )


def _sim_params() -> CollisionSimParams:
    return CollisionSimParams(path_sampling_step=0.001, safety_margin=1.0)


# ---------------------------------------------------------------------------
# Helpers unitaires
# ---------------------------------------------------------------------------

def test_direction_for_angles_matches_reversed_composition_convention():
    """Meme convention que step8/step9 : axis_angles[-1] applique en premier."""
    axes = [
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 180.0)),
    ]
    # B=90deg bascule (0,0,-1) vers (0,1,0) ; C=0 ne change rien ensuite.
    direction = _direction_for_angles(axes, [0.0, 90.0], np.array([0.0, 0.0, -1.0]))
    assert np.allclose(direction, [0.0, 1.0, 0.0], atol=1e-9)


def test_direction_for_angles_zero_axes_returns_reference():
    direction = _direction_for_angles([], [], np.array([0.0, 0.0, -1.0]))
    assert np.allclose(direction, [0.0, 0.0, -1.0])


def test_rate_limits_ignored_when_max_angular_velocity_unset():
    machine = _z_axis_machine(_offset_toolhead())
    machine.max_angular_velocity = 0.0  # non specifiee
    toolpath = _wall_then_probe_toolpath(np.array([0.0, 0.0, 5.0]))
    assert not _violates_rate_limits(toolpath.points, 2, [180.0], machine, neighbor_margin=1)


def test_rate_limits_reject_large_jump_when_velocity_constrained():
    machine = _z_axis_machine(_offset_toolhead())
    machine.max_angular_velocity = 0.001  # rad/s, tres restrictif
    toolpath = _wall_then_probe_toolpath(np.array([0.0, 0.0, 5.0]))
    assert _violates_rate_limits(toolpath.points, 2, [180.0], machine, neighbor_margin=1)


def test_collisions_at_points_only_evaluates_requested_indices():
    machine = _z_axis_machine(_offset_toolhead())
    toolpath = _wall_then_probe_toolpath(np.array([0.0, 0.0, 5.0]))

    events_all = _collisions_at_points(toolpath, machine, _sim_params(), [0, 1, 2])
    events_early_only = _collisions_at_points(toolpath, machine, _sim_params(), [0, 1])

    assert any(e.point_index == 2 for e in events_all)
    assert events_early_only == []


# ---------------------------------------------------------------------------
# Non-regression : trajectoire deja valide
# ---------------------------------------------------------------------------

def test_valid_trajectory_needs_no_adjustment():
    machine = _z_axis_machine(_offset_toolhead())
    toolpath = _wall_then_probe_toolpath(np.array([500.0, 500.0, 5.0]))  # sonde tres loin du mur
    report = collision_run(toolpath, machine, _sim_params())
    assert not report.has_collisions, "fixture invalide : la trajectoire devait etre saine"

    adjusted, adjustment = adjustment_run(
        toolpath, report, machine, _sim_params(), Step10Params()
    )

    assert adjustment.entries == []
    assert adjustment.converged is True
    assert adjustment.unresolved_collisions == []
    # trajectoire inchangee (angles compris)
    for before, after in zip(toolpath.points, adjusted.points):
        assert after.axis_angles == before.axis_angles
        assert np.allclose(after.position, before.position)


# ---------------------------------------------------------------------------
# Niveau 1 - reorientation locale reussie
# ---------------------------------------------------------------------------

def test_local_reorientation_resolves_collision():
    machine = _z_axis_machine(_offset_toolhead())
    toolpath = _wall_then_probe_toolpath(np.array([0.0, 0.0, 5.0]))
    report = collision_run(toolpath, machine, _sim_params())
    assert report.has_collisions, "fixture invalide : une collision etait attendue"

    adjusted, adjustment = adjustment_run(
        toolpath, report, machine, _sim_params(), Step10Params()
    )

    assert adjustment.converged is True
    assert [entry.method for entry in adjustment.entries] == ["local_reorientation"]
    assert adjustment.unresolved_collisions == []
    # l'angle du point fautif a reellement change, et la collision a disparu
    assert adjusted.points[2].axis_angles != toolpath.points[2].axis_angles
    assert collision_run(adjusted, machine, _sim_params()).events == []


def test_try_local_reorientation_returns_none_when_no_angle_helps():
    machine = _z_axis_machine(_axis_straddling_toolhead())
    toolpath = _wall_through_axis_toolpath(np.array([0.0, 0.0, 5.0]))
    assert collision_run(toolpath, machine, _sim_params()).has_collisions

    context = _build_reorientation_context(toolpath, machine, _sim_params())
    assert _try_local_reorientation(toolpath, 2, machine, _sim_params(), Step10Params(), context) is None


def test_reorientation_rejected_when_it_would_break_overhang_threshold():
    """
    Machine a 2 axes : l'axe B change reellement la direction de depot. Avec
    un seuil de surplomb a 0, aucun candidat autre que l'orientation
    d'origine n'est admissible - la reorientation doit echouer plutot que de
    produire une trajectoire qui reintroduit un surplomb.
    """
    machine = MachineProfile(
        rotation_axes=[
            RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
            RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 180.0)),
        ],
        platform_origin=None,
        nozzle_model=_axis_straddling_toolhead(),
        nozzle_diameter=0.4,
    )
    toolpath = ToolpathMachine(
        points=[
            ToolpathPoint(np.array([0.0, 0.0, 0.0]), [0.0, 0.0], 1.0, 1000.0),
            ToolpathPoint(np.array([0.0, 0.0, 10.0]), [0.0, 0.0], 1.0, 1000.0),
            ToolpathPoint(np.array([0.0, 0.0, 5.0]), [0.0, 0.0], 0.0, 1000.0),
        ]
    )
    params = Step10Params(max_direction_deviation_deg=0.0)

    context = _build_reorientation_context(toolpath, machine, _sim_params())
    assert _try_local_reorientation(toolpath, 2, machine, _sim_params(), params, context) is None


# ---------------------------------------------------------------------------
# Niveau 2 - re-segmentation ciblee via hook
# ---------------------------------------------------------------------------

def test_resegmentation_hook_used_when_reorientation_fails():
    machine = _z_axis_machine(_axis_straddling_toolhead())
    toolpath = _wall_through_axis_toolpath(np.array([0.0, 0.0, 5.0]))
    report = collision_run(toolpath, machine, _sim_params())
    assert report.has_collisions

    received: list[ResegmentationRequest] = []

    def hook(request: ResegmentationRequest) -> ToolpathMachine:
        received.append(request)
        return _wall_through_axis_toolpath(np.array([500.0, 500.0, 5.0]))  # zone degagee

    adjusted, adjustment = adjustment_run(
        toolpath, report, machine, _sim_params(), Step10Params(), resegmentation_hook=hook
    )

    assert adjustment.converged is True
    assert [entry.method for entry in adjustment.entries] == ["region_resegmentation"]
    assert len(received) == 1
    assert received[0].point_index == 2
    # la plage exclue est bien centree sur l'angle fautif, une entree par axe
    assert len(received[0].excluded_axis_ranges_deg) == 1
    low, high = received[0].excluded_axis_ranges_deg[0]
    assert low < 0.0 < high
    assert collision_run(adjusted, machine, _sim_params()).events == []


# ---------------------------------------------------------------------------
# Niveau 3 - support de secours
# ---------------------------------------------------------------------------

def test_fallback_support_when_no_hook_and_no_reorientation():
    machine = _z_axis_machine(_axis_straddling_toolhead())
    toolpath = _wall_through_axis_toolpath(np.array([0.0, 0.0, 5.0]))
    report = collision_run(toolpath, machine, _sim_params())
    assert report.has_collisions

    adjusted, adjustment = adjustment_run(
        toolpath, report, machine, _sim_params(), Step10Params()
    )

    assert [entry.method for entry in adjustment.entries] == ["fallback_support"]
    assert adjustment.entries[0].point_index == 2
    assert adjustment.unresolved_collisions == []  # traite (par support), pas laisse en suspens


def test_unconverged_when_hook_never_resolves_collision():
    """
    Hook qui rend une trajectoire toujours fautive : a l'epuisement des
    iterations, la collision restante bascule en support de secours et
    `converged` passe a False - jamais une exception.
    """
    machine = _z_axis_machine(_axis_straddling_toolhead())
    toolpath = _wall_through_axis_toolpath(np.array([0.0, 0.0, 5.0]))
    report = collision_run(toolpath, machine, _sim_params())

    def useless_hook(request: ResegmentationRequest) -> ToolpathMachine:
        return _wall_through_axis_toolpath(np.array([0.0, 0.0, 5.0]))  # toujours en collision

    adjusted, adjustment = adjustment_run(
        toolpath,
        report,
        machine,
        _sim_params(),
        Step10Params(max_adjustment_iterations=2),
        resegmentation_hook=useless_hook,
    )

    assert adjustment.converged is False
    assert adjustment.iterations_used == 2
    assert any(entry.method == "fallback_support" for entry in adjustment.entries)
    assert any(entry.method == "region_resegmentation" for entry in adjustment.entries)


# ---------------------------------------------------------------------------
# Step8Diagnostic - saute le niveau 1
# ---------------------------------------------------------------------------

def test_step8_diagnostic_skips_reorientation_and_goes_to_resegmentation():
    """
    Une orientation inatteignable a deja epuise la recherche IK a l'etape 8 :
    relancer le niveau 1 a l'identique ne peut rien donner, le niveau 2 doit
    etre sollicite directement.
    """
    machine = _z_axis_machine(_offset_toolhead())
    toolpath = _wall_then_probe_toolpath(np.array([500.0, 500.0, 5.0]))  # aucune collision
    report = collision_run(toolpath, machine, _sim_params())
    assert not report.has_collisions

    diagnostics = [
        Step8Diagnostic(
            point_index=2,
            layer_id=0,
            segment_index=0,
            residual_angle_error_deg=42.0,
            position=toolpath.points[2].position,
        )
    ]
    received: list[ResegmentationRequest] = []

    def hook(request: ResegmentationRequest) -> ToolpathMachine:
        received.append(request)
        return toolpath

    _adjusted, adjustment = adjustment_run(
        toolpath,
        report,
        machine,
        _sim_params(),
        Step10Params(),
        step8_diagnostics=diagnostics,
        resegmentation_hook=hook,
    )

    assert [entry.method for entry in adjustment.entries] == ["region_resegmentation"]
    assert received[0].reason == "orientation inatteignable (etape 8)"
    assert adjustment.converged is True


def test_step8_diagnostic_falls_back_to_support_without_hook():
    machine = _z_axis_machine(_offset_toolhead())
    toolpath = _wall_then_probe_toolpath(np.array([500.0, 500.0, 5.0]))
    report = collision_run(toolpath, machine, _sim_params())

    diagnostics = [
        Step8Diagnostic(
            point_index=2,
            layer_id=0,
            segment_index=0,
            residual_angle_error_deg=42.0,
            position=toolpath.points[2].position,
        )
    ]

    _adjusted, adjustment = adjustment_run(
        toolpath, report, machine, _sim_params(), Step10Params(), step8_diagnostics=diagnostics
    )

    assert [entry.method for entry in adjustment.entries] == ["fallback_support"]


# ---------------------------------------------------------------------------
# Purete : l'entree n'est jamais mutee
# ---------------------------------------------------------------------------

def test_input_toolpath_is_never_mutated():
    machine = _z_axis_machine(_offset_toolhead())
    toolpath = _wall_then_probe_toolpath(np.array([0.0, 0.0, 5.0]))
    original_angles = [list(point.axis_angles) for point in toolpath.points]
    report = collision_run(toolpath, machine, _sim_params())

    adjustment_run(toolpath, report, machine, _sim_params(), Step10Params())

    assert [list(point.axis_angles) for point in toolpath.points] == original_angles


# ---------------------------------------------------------------------------
# clamp_head_height (borne de l'angle de tete selon la hauteur)
# ---------------------------------------------------------------------------

def _clamp_machine():
    from core.mesh_types import MachineProfile, RotationAxisSpec
    head = trimesh.creation.box(extents=[2.0, 2.0, 40.0])
    head.apply_translation([0.0, 0.0, 20.0])            # corps de 40 mm au-dessus de la pointe
    return MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
                       RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))],
        platform_origin=np.array([0.0, 0.0, 0.0]), nozzle_model=head, nozzle_diameter=0.4,
        position_coupled_axes=[0])


def test_clamp_head_height_moves_only_points_whose_head_would_dip():
    from core.mesh_types import ToolpathMachine, ToolpathPoint
    from pipeline.step10_adjustment import clamp_head_height
    machine = _clamp_machine()
    sim = CollisionSimParams(path_sampling_step=1.0, safety_margin=1.0, nozzle_contact_clearance=0.0)
    pts = [ToolpathPoint(position=np.array([10.0, 0.0, 2.0]), axis_angles=[90.0, 170.0], extrusion=1.0,
                         feedrate=1.0),                                      # tete retournee pres du plateau
           ToolpathPoint(position=np.array([10.0, 0.0, 60.0]), axis_angles=[90.0, 170.0], extrusion=1.0,
                         feedrate=1.0),                                      # haut : rien a faire
           ToolpathPoint(position=np.array([10.0, 0.0, 2.0]), axis_angles=[90.0, 10.0], extrusion=1.0,
                         feedrate=1.0)]                                      # tete vers le bas : rien a faire
    tp = ToolpathMachine(points=pts)
    out, moved, dev = clamp_head_height(tp, machine, sim)
    assert moved == 1 and len(dev) == 1
    b = out.points[0].axis_angles[1]
    # corps de 40 mm (+ demi-largeur 1 mm) : z pointe + 40 cos(B) - ... >= marge
    assert 2.0 + 40.0 * np.cos(np.radians(b)) - 1.0 * abs(np.sin(np.radians(b))) >= 1.0 - 0.6
    assert b < 170.0
    assert out.points[1].axis_angles == [90.0, 170.0] and out.points[2].axis_angles == [90.0, 10.0]
    assert tp.points[0].axis_angles == [90.0, 170.0]                    # entree non mutee


def test_clamp_head_height_refuses_non_vertical_other_axes():
    from core.mesh_types import RotationAxisSpec, ToolpathMachine
    from pipeline.step10_adjustment import clamp_head_height
    machine = _clamp_machine()
    machine.rotation_axes[0] = RotationAxisSpec(axis=np.array([0.0, 1.0, 0.0]), angle_range_deg=(0.0, 360.0))
    with pytest.raises(ValueError):
        clamp_head_height(ToolpathMachine(points=[]), machine,
                          CollisionSimParams(path_sampling_step=1.0, safety_margin=1.0))


def test_reorient_coupled_tiers_clears_contact_and_tapers_neighbours():
    from core.mesh_types import MachineProfile, RotationAxisSpec, ToolpathMachine, ToolpathPoint
    from pipeline.step9_collision_sim import FastCollisionChecker, _filter_nozzle_contact_zone, _surface_sample_points
    from pipeline.step10_adjustment import reorient_coupled_tiers
    rod = trimesh.creation.box(extents=[1.0, 1.0, 20.0])
    rod.apply_translation([0.0, 0.0, 10.0])
    machine = MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
                       RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))],
        platform_origin=None, nozzle_model=rod, nozzle_diameter=0.4, position_coupled_axes=[0],
        radial_direction=np.array([0.0, -1.0, 0.0]))
    sim = CollisionSimParams(path_sampling_step=0.5, safety_margin=1.0, nozzle_contact_clearance=1.5,
                             toolhead_sample_spacing=1.0)
    pts = [ToolpathPoint(position=np.array([9.5, 0.0, 5.0]), axis_angles=[90.0, 0.0], extrusion=1.0, feedrate=1.0),
           ToolpathPoint(position=np.array([10.5, 0.0, 5.0]), axis_angles=[90.0, 0.0], extrusion=1.0, feedrate=1.0)]
    for i, x in enumerate(np.arange(6.0, 14.01, 0.5)):
        pts.append(ToolpathPoint(position=np.array([x, 0.0, 0.0]), axis_angles=[90.0, 0.0],
                                 extrusion=0.0 if i == 0 else 1.0, feedrate=1.0))
    tp = ToolpathMachine(points=pts)
    idx = 2 + 8                                            # x = 10 : la tige traverse le cordon
    samples = _filter_nozzle_contact_zone(_surface_sample_points(rod, 1.0), 1.5)
    checker = FastCollisionChecker(tp, machine, sim, samples)
    assert checker.penetrations(idx, pts[idx].position, pts[idx].axis_angles)[1] is not None
    out, log = reorient_coupled_tiers(tp, [idx], machine, sim)
    (i, b, tier), = log
    assert i == idx and b is not None and tier == 0 and 0.0 < b <= 45.0
    assert checker.penetrations(idx, out.points[idx].position, out.points[idx].axis_angles) == (None, None)
    bs = [p.axis_angles[1] for p in out.points]
    assert bs[idx - 2] == pytest.approx(b)                 # +-1 mm : correction entiere
    assert 0.0 < bs[idx - 4] < b                           # 2 mm : attenuee
    assert bs[idx - 6] == 0.0                              # 3 mm : intacte
    assert all(p.axis_angles[0] == 90.0 for p in out.points)
    assert tp.points[idx].axis_angles == [90.0, 0.0]


# ---------------------------------------------------------------------------
# Machine couplee (limite #59) et mode rapide de step9 (session du 2026-10-01)
# ---------------------------------------------------------------------------

def test_coupled_candidates_never_change_the_position_imposed_axis():
    import numpy as np
    from core.mesh_types import RotationAxisSpec
    from pipeline.step10_adjustment import _coupled_candidates
    axes = [RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
            RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))]
    dirs, angles = _coupled_candidates(axes, [123.4, 20.0], np.array([0.0, 0.0, 1.0]), 48)
    assert np.all(angles[:, 0] == 123.4)
    assert angles[:, 1].min() == 0.0 and angles[:, 1].max() == 185.0 and len(angles) == 48
    assert np.allclose(np.linalg.norm(dirs, axis=1), 1.0)


def test_fast_recheck_matches_step9_fast_events():
    """En mode rapide, `_collisions_at_points` rend exactement les evenements
    de step9 rapide aux memes indices (meme controleur, memes echantillons)."""
    import numpy as np
    import trimesh
    from core.mesh_types import MachineProfile, RotationAxisSpec, ToolpathMachine, ToolpathPoint
    from pipeline import step9_collision_sim as s9
    from pipeline.step10_adjustment import _collisions_at_points
    head = trimesh.creation.box(extents=[6, 6, 20])
    head.apply_translation([0, 0, 13])           # corps au-dessus de la pointe
    machine = MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
                       RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))],
        platform_origin=np.array([0.0, 0.0, -50.0]), nozzle_model=head, nozzle_diameter=0.4)
    pts = []
    for k in range(60):                           # une colonne haute, puis un passage dessous a cote
        pts.append(ToolpathPoint(position=np.array([10.0, 0.0, 0.2 * k]), axis_angles=[90.0, 0.0],
                                 extrusion=1.0, feedrate=1.0, color_id=None))
    for k in range(40):
        pts.append(ToolpathPoint(position=np.array([10.0 + 0.25 * k, 2.0, 0.2]), axis_angles=[90.0, 0.0],
                                 extrusion=1.0, feedrate=1.0, color_id=None))
    tp = ToolpathMachine(points=pts)
    params = s9.CollisionSimParams(path_sampling_step=0.1, safety_margin=1.0, fast_kdtree=True)
    ref = {(e.point_index, e.obstacle_type, round(e.penetration_depth, 12)) for e in s9.run(tp, machine, params).events}
    idx = sorted({i for i, _, _ in ref})
    got = {(e.point_index, e.obstacle_type, round(e.penetration_depth, 12))
           for e in _collisions_at_points(tp, machine, params, idx)}
    assert ref and got == ref


def test_exact_recheck_matches_step9_exact_events():
    """En mode EXACT, `_collisions_at_points` doit rendre les evenements de
    step9 exact aux memes indices. Bug reel (corrige le 2026-10-01) : step10
    oubliait `_filter_nozzle_contact_zone`, si bien que la pointe touchait
    toujours le cordon qu'elle depose et que tout point re-teste etait declare
    en collision - aucune reorientation ne pouvait reussir en mode exact."""
    from pipeline import step9_collision_sim as s9
    head = trimesh.creation.box(extents=[6, 6, 20])
    head.apply_translation([0, 0, 10.05])        # corps au-dessus de la pointe, pointe a l'origine
    machine = MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
                       RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))],
        platform_origin=np.array([0.0, 0.0, -50.0]), nozzle_model=head, nozzle_diameter=0.4)
    pts = []
    for k in range(60):
        pts.append(ToolpathPoint(position=np.array([10.0, 0.0, 0.2 * k]), axis_angles=[90.0, 0.0],
                                 extrusion=1.0, feedrate=1.0, color_id=None))
    for k in range(40):
        pts.append(ToolpathPoint(position=np.array([10.0 + 0.25 * k, 2.0, 0.2]), axis_angles=[90.0, 0.0],
                                 extrusion=1.0, feedrate=1.0, color_id=None))
    tp = ToolpathMachine(points=pts)
    params = s9.CollisionSimParams(path_sampling_step=0.1, safety_margin=1.0, fast_kdtree=False)
    assert params.nozzle_contact_clearance > 0
    ref = {(e.point_index, e.obstacle_type, round(e.penetration_depth, 12)) for e in s9.run(tp, machine, params).events}
    sampled = s9._sample_point_indices(tp.points, params.path_sampling_step)
    got = {(e.point_index, e.obstacle_type, round(e.penetration_depth, 12))
           for e in _collisions_at_points(tp, machine, params, sampled)}
    assert ref, "la fixture doit produire des collisions, sinon le test ne verifie rien"
    assert len({i for i, _, _ in ref}) < len(sampled), "et des points sans collision"
    assert got == ref
