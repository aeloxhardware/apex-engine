"""
tests/test_step8b_travel_routing.py

Tests pour pipeline/step8b_travel_routing.py (CLAUDE.md limite #46).
"""

from __future__ import annotations

import numpy as np
import trimesh

from core.mesh_types import MachineProfile, ToolpathMachine, ToolpathPoint, TravelRoutingDiagnostic
from pipeline.step8_machine_toolpath import Step8Params
from pipeline.step8b_travel_routing import (
    TravelRoutingParams,
    _travel_boundaries,
    route_long_travels,
)
from pipeline.step9_collision_sim import CollisionSimParams


def _small_toolhead() -> trimesh.Trimesh:
    box = trimesh.creation.box(extents=[4.0, 4.0, 8.0])
    box.apply_translation([0.0, 0.0, 4.0])  # pointe (bas) a l'origine locale
    return box


def _machine_no_rotation(nozzle_model=None) -> MachineProfile:
    return MachineProfile(rotation_axes=[], platform_origin=None, nozzle_diameter=0.4, nozzle_model=nozzle_model)


def _point(position, extrusion=1.0) -> ToolpathPoint:
    return ToolpathPoint(position=np.array(position, dtype=float), axis_angles=[], extrusion=extrusion, feedrate=1000.0)


def _default_sim_params() -> CollisionSimParams:
    return CollisionSimParams(path_sampling_step=1.0, safety_margin=0.5, toolhead_sample_spacing=1.0)


# ---------------------------------------------------------------------------
# _travel_boundaries
# ---------------------------------------------------------------------------

def test_travel_boundaries_detects_extrusion_zero_points():
    points = [_point([0, 0, 0]), _point([1, 0, 0]), _point([10, 0, 0], extrusion=0.0), _point([11, 0, 0])]
    assert _travel_boundaries(points) == [2]


def test_travel_boundaries_empty_list():
    assert _travel_boundaries([]) == []


# ---------------------------------------------------------------------------
# route_long_travels - comportement de base
# ---------------------------------------------------------------------------

def test_short_travel_is_left_unchanged():
    points = [_point([0, 0, 0]), _point([1, 0, 0], extrusion=0.0), _point([2, 0, 0])]
    toolpath = ToolpathMachine(points=points)
    params = TravelRoutingParams(length_threshold=5.0, lift_height=3.0)

    routed, diagnostics = route_long_travels(toolpath, _machine_no_rotation(), _default_sim_params(), params)

    assert len(routed.points) == len(points)
    assert diagnostics == []


def test_long_travel_without_toolhead_gets_lifted_without_diagnostic():
    """
    Sans modele de tete charge, aucune validation de collision n'est
    possible : le lift standard est applique sans essai d'augmentation ni
    diagnostic (rien a juger) - voir docstring de route_long_travels.
    """
    points = [_point([0, 0, 0]), _point([50, 0, 0], extrusion=0.0), _point([51, 0, 0])]
    toolpath = ToolpathMachine(points=points)
    params = TravelRoutingParams(length_threshold=5.0, lift_height=3.0, traverse_samples=2)

    routed, diagnostics = route_long_travels(toolpath, _machine_no_rotation(), _default_sim_params(), params)

    assert len(routed.points) > len(points)
    assert diagnostics == []
    # Les points intermediaires inseres doivent tous porter extrusion=0.0
    inserted = routed.points[1:-1]
    assert all(p.extrusion == 0.0 for p in inserted)
    # Un point au moins doit etre leve en Z au-dessus des deux extremites
    max_z = max(p.position[2] for p in inserted)
    assert max_z >= params.lift_height - 1e-9


def test_long_travel_with_toolhead_avoids_printed_capsule_via_lift():
    """
    Une trajectoire imprimee traverse le plan z=0 entre les deux points
    d'extrusion qui encadrent le trajet a vide - une ligne droite non
    routee passerait dedans, mais le chemin lift-and-traverse (a
    lift_height au-dessus) doit rester hors de la capsule de securite.
    """
    toolhead = _small_toolhead()
    machine = _machine_no_rotation(nozzle_model=toolhead)

    printed_a = _point([0, 0, 0])
    printed_b = _point([0, 0, -0.1])  # capsule d'un cordon deja imprime, proche de z=0
    travel_start = _point([0, 5, 0])
    travel_end = _point([0, -5, 0], extrusion=0.0)  # traversee juste au-dessus du cordon si non routee
    after = _point([0, -6, 0])

    points = [printed_a, printed_b, travel_start, travel_end, after]
    toolpath = ToolpathMachine(points=points)
    params = TravelRoutingParams(length_threshold=1.0, lift_height=5.0, max_lift_multiplier=2.0, traverse_samples=3)

    routed, diagnostics = route_long_travels(toolpath, machine, _default_sim_params(), params)

    assert len(routed.points) > len(points)
    assert diagnostics == []  # le lift initial doit suffire a eviter la capsule


def test_diagnostic_emitted_when_even_max_lift_still_collides():
    """
    Une "toiture" de matiere imprimee, en ligne continue le long de tout le
    trajet de traverse (y de -5 a 5) et a une hauteur superieure au lift
    plafond teste (max = lift_height * max_lift_multiplier = 6.0), force
    l'echec du routage a toutes les hauteurs testees (3.0 puis 6.0) : un
    TravelRoutingDiagnostic doit etre emis plutot qu'un succes silencieux
    faux.
    """
    toolhead = _small_toolhead()
    machine = _machine_no_rotation(nozzle_model=toolhead)

    # Ligne continue de cordons imprimes a z=8 (au-dessus du lift plafond
    # de 6.0), couvrant tout le segment y en [-5,5] traverse par le trajet
    # a vide - AUCUNE hauteur de lift testee (3.0 ou 6.0) ne peut l'eviter.
    wall_points = [_point([0, y, 8.0]) for y in np.linspace(-5, 5, 11)]
    travel_start = _point([0, 5, 0])
    travel_end = _point([0, -5, 0], extrusion=0.0)
    after = _point([0, -6, 0])

    points = wall_points + [travel_start, travel_end, after]
    toolpath = ToolpathMachine(points=points)
    params = TravelRoutingParams(
        length_threshold=1.0, lift_height=3.0, max_lift_multiplier=2.0, traverse_samples=5,
        routing_toolhead_sample_spacing=1.0,
    )

    routed, diagnostics = route_long_travels(toolpath, machine, _default_sim_params(), params)

    assert len(diagnostics) == 1
    diagnostic = diagnostics[0]
    assert isinstance(diagnostic, TravelRoutingDiagnostic)
    assert diagnostic.still_colliding is True
    assert diagnostic.max_lift_tested >= params.lift_height


def test_extrusion_points_are_preserved_unchanged():
    """Les points en extrusion (>0) ne sont jamais modifies ni reordonnes par le routage."""
    points = [_point([0, 0, 0]), _point([1, 0, 0]), _point([50, 0, 0], extrusion=0.0), _point([51, 0, 0]), _point([52, 0, 0])]
    toolpath = ToolpathMachine(points=points)
    params = TravelRoutingParams(length_threshold=5.0, lift_height=2.0)

    routed, _ = route_long_travels(toolpath, _machine_no_rotation(), _default_sim_params(), params)

    extrusion_points_before = [p for p in points if p.extrusion > 0.0]
    extrusion_points_after = [p for p in routed.points if p.extrusion > 0.0]
    assert len(extrusion_points_before) == len(extrusion_points_after)
    for before, after in zip(extrusion_points_before, extrusion_points_after):
        assert np.allclose(before.position, after.position)
