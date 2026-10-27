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
    _filter_nozzle_contact_zone,
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
# Zone de garde autour de la pointe de buse (nozzle_contact_clearance)
# ---------------------------------------------------------------------------

def test_filter_nozzle_contact_zone_removes_only_points_within_clearance():
    samples = np.array(
        [
            [0.0, 0.0, 0.0],   # a l'origine (la pointe elle-meme) - exclu
            [1.0, 0.0, 0.0],   # distance 1.0 < clearance=3.0 - exclu
            [2.9, 0.0, 0.0],   # juste sous le seuil - exclu
            [3.0, 0.0, 0.0],   # exactement au seuil - conserve (>=)
            [10.0, 0.0, 0.0],  # loin - conserve
        ]
    )
    filtered = _filter_nozzle_contact_zone(samples, clearance=3.0)
    assert np.allclose(filtered, [[3.0, 0.0, 0.0], [10.0, 0.0, 0.0]])


def test_filter_nozzle_contact_zone_disabled_at_zero_clearance():
    """clearance<=0 doit rendre le tableau strictement inchange - comportement historique."""
    samples = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    filtered = _filter_nozzle_contact_zone(samples, clearance=0.0)
    assert np.array_equal(filtered, samples)


def test_filter_nozzle_contact_zone_falls_back_when_everything_excluded():
    """
    Une tete plus petite que la zone de garde ne doit jamais produire un
    tableau vide (planterait _toolhead_aabb en aval et desactiverait toute
    detection pour cette tete) - repli sur le tableau d'origine intact.
    """
    samples = np.array([[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]])
    filtered = _filter_nozzle_contact_zone(samples, clearance=3.0)
    assert np.array_equal(filtered, samples)


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


def test_nozzle_tip_contact_with_own_cordon_no_longer_false_positive():
    """
    Reproduit le motif mesure sur l_shape_overhang.stl (limite "PRIORITE
    IMMEDIATE") : la pointe de buse (origine locale du modele de tete) est
    par construction au contact du cordon qu'elle vient de deposer - avant
    nozzle_contact_clearance, ce seul point de contact declenchait un
    printed_part a chaque point de trajectoire en extrusion continue, sans
    aucun obstacle reel autre que le cordon legitimement depose juste avant
    (limite #22 n'exclut que LE segment qui se termine au point teste, pas
    ceux, tres proches, deposes juste avant sur la meme trajectoire quasi
    verticale). Avec le defaut (3.0mm), une trajectoire d'extrusion continue
    sans obstacle etranger ne doit plus produire ce faux positif.

    Trajectoire elevee (z de base 20mm) pour isoler cette question du
    plateau, deja couvert par le fixture historique
    (test_collision_against_printed_wall_detected) - pas l'objet de ce test.
    """
    toolhead = _small_toolhead()
    machine = _machine_no_rotation()
    machine.nozzle_model = toolhead

    # Extrusion continue quasi verticale : chaque point est a ~0.3mm du
    # precedent, largement a l'interieur de safety_margin - sans zone de
    # garde, la pointe touche le cordon tout juste depose a chaque pas.
    points = [
        ToolpathPoint(
            position=np.array([0.0, 0.0, 20.0 + float(i) * 0.3]),
            axis_angles=[],
            extrusion=1.0,
            feedrate=1000.0,
        )
        for i in range(8)
    ]
    toolpath = ToolpathMachine(points=points)
    params = CollisionSimParams(path_sampling_step=0.001, safety_margin=1.0)

    report = step9_run(toolpath, machine, params)

    assert not report.has_collisions


def test_obstacle_beyond_clearance_still_detected():
    """
    La zone de garde ne doit exclure QUE le voisinage immediat de la pointe
    - un obstacle touchant le CORPS de la tete (au-dela de clearance) doit
    toujours declencher une collision, sinon le fix aveuglerait toute
    detection au lieu de corriger seulement le faux positif de contact.
    """
    toolhead = _small_toolhead()  # boite x:+-2,y:+-3,z:0->10, tip a l'origine
    machine = _machine_no_rotation()
    machine.nozzle_model = toolhead

    # imprime un mur vertical, puis repositionne la tete de sorte que son
    # CORPS (autour de z=8, largement au-dela de clearance=3.0 de la pointe
    # a z=0) chevauche le mur - la pointe elle-meme (a z=8 dans ce cas) reste
    # loin du mur en XY tant qu'on ne translate qu'en Z, donc translater la
    # base directement pour placer le corps sur le mur.
    points = [
        ToolpathPoint(position=np.array([0.0, 0.0, 0.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([0.0, 0.0, 20.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([0.0, 0.0, 10.0]), axis_angles=[], extrusion=0.0, feedrate=1000.0),
    ]
    toolpath = ToolpathMachine(points=points)
    params = CollisionSimParams(path_sampling_step=0.001, safety_margin=1.0)

    report = step9_run(toolpath, machine, params)

    assert any(e.obstacle_type == "printed_part" for e in report.events)


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

    # Le loader ramene le corps de la tete selon +Z machine (rotation +90 deg
    # autour de X, voir mesh_loader.TOOLHEAD_BODY_TO_MACHINE_UP) : l'axe long
    # de 77.15 mm est donc en Z apres chargement, pas en Y comme dans le
    # fichier brut.
    extents = mesh.extents
    assert extents[0] == pytest.approx(76.0, abs=1.0)
    assert extents[1] == pytest.approx(35.0, abs=1.0)
    assert extents[2] == pytest.approx(76.55, abs=1.0)   # tete v4 de l'interface (R-009) ; v5 du moteur : 77,15

    distances_to_origin = np.linalg.norm(mesh.vertices, axis=1)
    # la pointe (origine) doit etre quasi sur le maillage : la tete v4 a une face de buse plate de 1,1 mm de rayon
    # centree sur l'origine (le sommet le plus proche est a 1,10 mm ; la v5 en avait un a 0,1 mm)
    assert distances_to_origin.min() < 1.5
    assert abs(mesh.bounds[0][2]) < 0.01   # la face de buse est dans le plan z = 0

    # Propriete PHYSIQUE qui compte reellement, et que l'ancien ordre d'axes
    # laissait passer : la tete doit se trouver entierement AU-DESSUS de sa
    # pointe de buse. Sans la rotation du loader, le mesh descendait jusqu'a
    # z=-18.9 (18.9 mm de tete sous le point de depot), ce qui faisait
    # collisionner ~100% des points de trajectoire contre le plateau et la
    # piece - voir mesh_loader.TOOLHEAD_BODY_TO_MACHINE_UP.
    assert mesh.bounds[0][2] > -1.0, "le corps de la tete ne doit jamais pendre sous la buse"
    assert mesh.bounds[1][2] == pytest.approx(76.55, abs=1.0)


# ---------------------------------------------------------------------------
# Controle rapide (fast_kdtree, opt-in)
# ---------------------------------------------------------------------------

def _layered_toolpath(seed: int) -> ToolpathMachine:
    """Polylignes empilees et croisees, dont certaines repassent pres de la
    matiere deja deposee : des collisions ET des passages libres."""
    rng = np.random.default_rng(seed)
    points = []
    for layer in range(6):
        z = 0.3 * layer
        for _line in range(4):
            a = rng.uniform(-8, 8, 2)
            b = rng.uniform(-8, 8, 2)
            for t_i, t in enumerate(np.linspace(0, 1, 25)):
                xy = a + t * (b - a)
                points.append(ToolpathPoint(position=np.array([xy[0], xy[1], z]), axis_angles=[],
                                            extrusion=0.0 if t_i == 0 else 1.0, feedrate=1000.0))
    return ToolpathMachine(points=points)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_fast_kdtree_matches_capsule_loop(seed):
    toolhead = _small_toolhead()
    machine = _machine_no_rotation()
    machine.nozzle_model = toolhead
    toolpath = _layered_toolpath(seed)
    slow = step9_run(toolpath, machine, CollisionSimParams(path_sampling_step=0.5, safety_margin=1.0,
                                                           nozzle_contact_clearance=1.5))
    fast = step9_run(toolpath, machine, CollisionSimParams(path_sampling_step=0.5, safety_margin=1.0,
                                                           nozzle_contact_clearance=1.5, fast_kdtree=True))
    key = lambda r: {(e.point_index, e.obstacle_type): e.penetration_depth for e in r.events}  # noqa: E731
    ks, kf = key(slow), key(fast)
    assert ks, "le jeu d'essai doit produire des collisions"
    # densification a 0,25 mm : la distance a l'axe du cordon est surestimee
    # d'au plus un demi-pas (0,125 mm, point de tete sur l'axe entre deux points)
    tol = 0.25 / 2 + 1e-9
    only_slow = {k for k in ks.keys() - kf.keys() if ks[k] > tol}
    only_fast = {k for k in kf.keys() - ks.keys()}
    assert not only_slow and not only_fast
    for k in ks.keys() & kf.keys():
        assert kf[k] <= ks[k] + 1e-9 and kf[k] == pytest.approx(ks[k], abs=tol)


def test_fast_kdtree_ignores_future_material_and_sees_past():
    """La tete survole une zone imprimee PLUS TARD (beaucoup de points
    futurs proches) : aucune collision ; puis un seul cordon passe proche
    doit etre vu malgre les points futurs qui saturent les K voisins."""
    toolhead = trimesh.creation.box(extents=[1.0, 1.0, 1.0])
    toolhead.apply_translation([0.0, 0.0, 0.5])      # pointe au centre du bas de la boite
    machine = _machine_no_rotation()
    machine.nozzle_model = toolhead
    pts = []

    def add(p, e):
        pts.append(ToolpathPoint(position=np.array(p, dtype=float), axis_angles=[], extrusion=e, feedrate=1.0))

    # tete testee au point 2 : boite x 0..1, y 4.5..5.5, z 3..4
    add([0.0, 5.0, 4.2], 1.0)            # cordon PASSE (deux points en extrusion), a ~0,54 mm du coin (0, 4.5, 4)
    add([1.0, 5.0, 4.2], 1.0)
    add([0.5, 5.0, 3.0], 0.0)            # point teste
    add([20.0, 0.0, 0.0], 0.0)
    for i in range(200):                 # matiere FUTURE, dense, collee a ce coin : sature les K voisins
        add([-0.1 + 0.01 * (i % 20), 4.4 + 0.02 * (i // 20), 4.15], 1.0 if i else 0.0)
    toolpath = ToolpathMachine(points=pts)
    params = CollisionSimParams(path_sampling_step=0.001, safety_margin=1.0, nozzle_contact_clearance=0.0,
                                fast_kdtree=True)
    slow = step9_run(toolpath, machine, CollisionSimParams(path_sampling_step=0.001, safety_margin=1.0,
                                                           nozzle_contact_clearance=0.0))
    fast = step9_run(toolpath, machine, params)
    hits_slow = {e.point_index for e in slow.events if e.obstacle_type == "printed_part"}
    hits_fast = {e.point_index for e in fast.events if e.obstacle_type == "printed_part"}
    assert 2 in hits_slow and 2 in hits_fast
    assert hits_fast == hits_slow


def test_platform_surface_is_at_origin_z():
    """La surface du plateau est a z = platform_origin.z (et non +0,5 mm : disque centre
    sur l'origine avant le 2026-10-06, limite #70). Marge 1 mm : un point a z = 0,9 est en
    contact (penetration 0,1), un point a z = 1,1 ne l'est pas, un point sous le plateau l'est."""
    from pipeline.step9_collision_sim import _check_platform_disc
    origin = np.array([0.0, 0.0, 0.0])
    pen = _check_platform_disc(np.array([[10.0, 0.0, 0.9]]), origin, 1.0)
    assert pen == pytest.approx(0.1, abs=1e-9)
    assert _check_platform_disc(np.array([[10.0, 0.0, 1.1]]), origin, 1.0) is None
    assert _check_platform_disc(np.array([[10.0, 0.0, -0.4]]), origin, 1.0) == pytest.approx(1.0)   # dans le plateau : penetration = marge
    # une origine decalee en z decale la surface d'autant
    assert _check_platform_disc(np.array([[10.0, 0.0, 5.9]]), np.array([0.0, 0.0, 5.0]), 1.0) == pytest.approx(0.1, abs=1e-9)


def test_fast_checker_is_reused_while_material_is_unchanged():
    """Le controleur rapide est partage tant que positions et extrusions sont les memes
    (changer un angle ne change pas la matiere) ; une position modifiee le reconstruit."""
    from core.mesh_types import MachineProfile, ToolpathMachine, ToolpathPoint
    from pipeline.step9_collision_sim import CollisionSimParams, fast_checker_for
    machine = MachineProfile(rotation_axes=[], platform_origin=None, nozzle_diameter=0.4)
    sim = CollisionSimParams(path_sampling_step=1.0, safety_margin=1.0)
    samples = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 2.0]])
    pts = [ToolpathPoint(position=np.array([float(i), 0.0, 0.0]), axis_angles=[], extrusion=1.0, feedrate=1.0)
           for i in range(5)]
    a = fast_checker_for(ToolpathMachine(points=pts), machine, sim, samples)
    turned = [ToolpathPoint(position=p.position.copy(), axis_angles=[], extrusion=1.0, feedrate=2.0) for p in pts]
    assert fast_checker_for(ToolpathMachine(points=turned), machine, sim, samples) is a
    moved = [ToolpathPoint(position=p.position + 0.1, axis_angles=[], extrusion=1.0, feedrate=1.0) for p in pts]
    b = fast_checker_for(ToolpathMachine(points=moved), machine, sim, samples)
    assert b is not a and np.allclose(b.dense.min(axis=0), a.dense.min(axis=0) + 0.1)


def test_printed_part_margin_defaults_to_safety_margin_and_can_be_negative():
    """Marge propre a la matiere deposee : None = safety_margin ; negative =
    recouvrement tolere (la tete est exacte, effleurer une couche n'est pas un choc)."""
    from pipeline.step9_collision_sim import CollisionSimParams, printed_margin_of
    base = CollisionSimParams(path_sampling_step=2.0, safety_margin=1.0)
    assert printed_margin_of(base) == 1.0
    tolerant = CollisionSimParams(path_sampling_step=2.0, safety_margin=1.0, printed_part_margin=-0.15)
    assert printed_margin_of(tolerant) == -0.15 and tolerant.safety_margin == 1.0
