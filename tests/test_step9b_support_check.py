"""
tests/test_step9b_support_check.py

Tests pour pipeline/step9b_support_check.py, sur des `ToolpathMachine`
synthetiques construites a la main - meme approche que tests/test_step9_*.py
et tests/test_step10_*.py, qui ne rejouent pas non plus tout le pipeline pour
tester une etape qui ne consomme qu'une trajectoire.

Execution : pytest tests/test_step9b_support_check.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from core.mesh_types import MachineProfile, RotationAxisSpec, ToolpathMachine, ToolpathPoint
from pipeline.step3_segmentation import SegmentationParams
from pipeline.step9_collision_sim import _PrintedGeometryIndex, _point_to_segment_distance
from pipeline.step9b_support_check import (
    SupportCheckParams,
    _has_support_below,
    run as support_check_run,
)


def _point(x: float, y: float, z: float, extrusion: float = 1.0) -> ToolpathPoint:
    return ToolpathPoint(
        position=np.array([x, y, z], dtype=float),
        axis_angles=[],
        extrusion=extrusion,
        feedrate=100.0,
    )


def _machine(nozzle_diameter: float = 0.4) -> MachineProfile:
    return MachineProfile(rotation_axes=[], nozzle_diameter=nozzle_diameter)


def _tilting_machine(nozzle_diameter: float = 0.4) -> MachineProfile:
    """Machine avec un seul axe de rotation reel (tete, X, 0-180deg) - meme
    convention que _default_machine_profile() (main.py)."""
    return MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 180.0))],
        nozzle_diameter=nozzle_diameter,
    )


# Defauts du projet, rappeles ici pour que les fixtures restent lisibles :
# nozzle_diameter=0.4 -> rayon horizontal 1.0 mm ; layer_height_base=0.2 ->
# hauteur verticale 1.5 mm ; tolerance plateau 1.0 mm.
LAYER_HEIGHT = 0.2
PLATFORM_Z = 0.0


def test_stacked_material_is_fully_supported():
    """Empilement vertical trivial : chaque couche repose sur la precedente."""
    toolpath = ToolpathMachine(
        points=[
            _point(0.0, 0.0, 0.0), _point(1.0, 0.0, 0.0),
            _point(1.0, 0.0, 0.2), _point(0.0, 0.0, 0.2),
            _point(0.0, 0.0, 0.4), _point(1.0, 0.0, 0.4),
        ]
    )
    report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    assert report.checked is True
    assert report.events == []
    assert report.tested_extruded_points == 6
    assert report.unsupported_fraction == 0.0


def test_platform_contact_needs_no_support_below():
    """
    Un point pose sur le plateau n'a structurellement rien en dessous, et ne
    doit jamais compter comme un depot dans le vide.
    """
    toolpath = ToolpathMachine(points=[_point(0.0, 0.0, 0.0), _point(1.0, 0.0, 0.0)])
    report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    assert report.events == []


def test_unknown_platform_never_exempts_silently():
    """
    `platform_z=None` (plateau inconnu) ne doit exempter AUCUN point : mieux
    vaut surestimer les evenements que masquer un vrai defaut derriere une
    hypothese de plateau non fournie.
    """
    toolpath = ToolpathMachine(points=[_point(0.0, 0.0, 0.0), _point(1.0, 0.0, 0.0)])
    report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, None, SupportCheckParams()
    )

    assert len(report.events) == 2


def test_ordering_defect_is_distinguished_from_real_overhang():
    """
    Le coeur de l'etape : un point depose en l'air dont le support est
    imprime PLUS TARD est un defaut d'ORDRE (corrigible par reordonnancement),
    pas un surplomb geometrique.
    """
    toolpath = ToolpathMachine(
        points=[
            _point(0.0, 0.0, 5.0), _point(1.0, 0.0, 5.0),   # depose en l'air
            _point(0.0, 0.0, 4.9), _point(1.0, 0.0, 4.9),   # son support, APRES
        ]
    )
    report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    assert report.ordering_defect_count == 2
    # Les points du support sont eux-memes en l'air (rien ne les porte) :
    # geometriquement correct, ce sont de vrais surplombs.
    assert report.true_overhang_count == 2


def test_real_overhang_has_no_support_at_any_time():
    """
    Un point qu'AUCUNE matiere de la trajectoire ne supporte, ni avant ni
    apres, est un surplomb reel - aucun ordre d'impression ne le resout.
    """
    toolpath = ToolpathMachine(
        points=[
            _point(0.0, 0.0, 0.0), _point(1.0, 0.0, 0.0),
            _point(50.0, 50.0, 9.0), _point(51.0, 50.0, 9.0),
        ]
    )
    report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    assert report.ordering_defect_count == 0
    assert report.true_overhang_count == 2
    assert all(event.supporting_point_index is None for event in report.events)


def test_supporting_point_index_is_opt_in():
    """
    La localisation du point porteur coute O(defauts x segments) : elle est
    desactivee par defaut, mais la distinction ordre/surplomb - la donnee qui
    decide quoi corriger - reste TOUJOURS calculee.
    """
    points = [
        _point(0.0, 0.0, 5.0), _point(1.0, 0.0, 5.0),
        _point(0.0, 0.0, 4.9), _point(1.0, 0.0, 4.9),
    ]
    toolpath = ToolpathMachine(points=points)

    default_report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )
    located_report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, PLATFORM_Z,
        SupportCheckParams(locate_supporting_point=True),
    )

    assert SupportCheckParams().locate_supporting_point is False
    assert default_report.ordering_defect_count == located_report.ordering_defect_count

    ordering_default = [e for e in default_report.events if e.has_future_support]
    ordering_located = [e for e in located_report.events if e.has_future_support]
    assert all(e.supporting_point_index is None for e in ordering_default)
    assert all(e.supporting_point_index is not None for e in ordering_located)


def test_search_radius_scales_with_nozzle_diameter():
    """
    Les seuils sont derives de `nozzle_diameter`, jamais des constantes
    absolues : un support horizontalement decale de 1.5 mm est hors de portee
    d'une buse de 0.4 mm (rayon 1.0 mm) mais dans celle d'une buse de 0.8 mm
    (rayon 2.0 mm).
    """
    points = [
        _point(0.0, 0.0, 1.0), _point(0.5, 0.0, 1.0),      # support
        _point(1.5, 0.0, 1.15), _point(2.0, 0.0, 1.15),    # decale de 1.5 mm
    ]
    toolpath = ToolpathMachine(points=points)

    narrow = support_check_run(
        toolpath, _machine(0.4), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )
    wide = support_check_run(
        toolpath, _machine(0.8), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    assert len(wide.events) < len(narrow.events)


def test_material_above_never_counts_as_support():
    """
    Une capsule AU-DESSUS du point ne le supporte pas : la gravite n'a qu'un
    sens. Sans ce test, une inversion de signe passerait inapercue.
    """
    # Le point teste (z=5.0) a de la matiere 0.15 mm AU-DESSUS de lui, bien
    # dans le rayon de recherche - mais au-dessus, donc sans effet portant.
    toolpath = ToolpathMachine(
        points=[
            _point(0.0, 0.0, 5.15), _point(1.0, 0.0, 5.15),  # matiere AU-DESSUS
            _point(40.0, 0.0, 5.15, extrusion=0.0),           # trajet a vide : coupe le cordon
            _point(0.0, 0.0, 5.0), _point(1.0, 0.0, 5.0),     # le point teste
        ]
    )
    report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    # Les 2 points du bas restent non supportes malgre la matiere juste
    # au-dessus d'eux : la gravite n'a qu'un sens.
    assert {event.point_index for event in report.events} == {0, 1, 3, 4}


def test_empty_toolpath_reports_not_checked():
    """
    `checked=False` est un signal EXPLICITE : une liste `events` vide serait
    sinon ambigue entre "verifie, rien a signaler" et "jamais verifie".
    """
    report = support_check_run(
        ToolpathMachine(points=[]), _machine(), LAYER_HEIGHT, PLATFORM_Z,
        SupportCheckParams(),
    )

    assert report.checked is False
    assert report.events == []
    assert report.unsupported_fraction == 0.0


def test_travel_moves_are_not_counted_as_deposits():
    """
    Un trajet a vide ne depose rien : il ne doit ni etre teste, ni diluer la
    fraction rapportee (dont le denominateur est le nombre de points
    REELLEMENT testes, pas la longueur de la trajectoire).
    """
    toolpath = ToolpathMachine(
        points=[
            _point(0.0, 0.0, 0.0), _point(1.0, 0.0, 0.0),
            _point(40.0, 40.0, 9.0, extrusion=0.0),   # trajet a vide en l'air
            _point(41.0, 40.0, 9.0, extrusion=0.0),
        ]
    )
    report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    assert report.tested_extruded_points == 2
    assert report.events == []


def test_tilt_resolves_previously_classified_overhang():
    """
    Un point sans support VERTICAL (ni avant, ni apres, ni dans la piece
    finie) mais dont un cordon existe le long d'une direction atteignable
    par la cinematique de la machine (ici : rotation de tete autour de X,
    0-180deg) n'est ni un defaut d'ordre, ni un surplomb reel - la machine
    sait le deposer en inclinant la tete.

    Reproduit en miniature le cas mesure sur death.stl (194/194 "surplombs
    reels" resolus des qu'on teste une direction inclinee au lieu de la
    seule verticale) : le point teste et son support sont a plus de 1.0mm
    (rayon horizontal) l'un de l'autre au sens vertical, mais a 1.3mm le
    long de la direction inclinee [0,-0.7705,0.6374] (~50deg de la
    verticale, a l'interieur de la plage 0-180deg de l'axe X) - a moins de
    1.0mm horizontalement et 1.5mm verticalement dans le repere LOCAL de
    cette direction.
    """
    # Le segment support est lui-meme ancre sur le plateau (chaque bout pose
    # a z=PLATFORM_Z) - sans quoi il serait teste independamment comme son
    # propre surplomb, ce qui n'est pas ce que ce test verifie.
    support_center = np.array([0.0, 1.002, 4.171])
    tested_point = _point(0.0, 0.0, 5.0)
    support_column_a = _point(support_center[0] - 0.3, support_center[1], PLATFORM_Z)
    support_top_a = _point(support_center[0] - 0.3, support_center[1], support_center[2])
    support_top_b = _point(support_center[0] + 0.3, support_center[1], support_center[2])
    support_column_b = _point(support_center[0] + 0.3, support_center[1], PLATFORM_Z)

    toolpath = ToolpathMachine(
        points=[
            support_column_a, support_top_a, support_top_b, support_column_b,
            tested_point, _point(1.0, 0.0, 5.0),
        ]
    )
    report = support_check_run(
        toolpath, _tilting_machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    events = {e.point_index: e for e in report.events}
    assert 4 in events
    assert events[4].has_future_support is False
    assert events[4].supported_by_tilt is True
    assert report.true_overhang_count == 0
    assert report.tilt_resolved_count == 1


def test_true_overhang_stays_true_when_no_tilt_helps():
    """
    Un point vraiment isole (aucun segment dans le rayon de recherche, a
    AUCUNE distance ni direction) reste un surplomb reel meme avec une
    machine capable d'incliner la tete - la passe multi-directionnelle ne
    doit jamais inventer un support qui n'existe pas.
    """
    toolpath = ToolpathMachine(
        points=[
            _point(0.0, 0.0, 0.0), _point(1.0, 0.0, 0.0),
            _point(50.0, 50.0, 9.0), _point(51.0, 50.0, 9.0),
        ]
    )
    report = support_check_run(
        toolpath, _tilting_machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    assert report.true_overhang_count == 2
    assert report.tilt_resolved_count == 0
    assert all(not e.supported_by_tilt for e in report.events)


def test_empty_rotation_axes_keeps_vertical_only_behavior():
    """
    `rotation_axes=[]` (profil sans cinematique de rotation declaree) doit
    retomber EXACTEMENT sur le comportement historique (verticale
    uniquement) - non-regression explicite : la fixture de
    test_tilt_resolves_previously_classified_overhang, rejouee avec
    `_machine()` (rotation_axes=[]) au lieu de `_tilting_machine()`, ne doit
    PAS resoudre le point par inclinaison (aucune direction candidate autre
    que la verticale n'existe).
    """
    support_center = np.array([0.0, 1.002, 4.171])
    tested_point = _point(0.0, 0.0, 5.0)
    support_column_a = _point(support_center[0] - 0.3, support_center[1], PLATFORM_Z)
    support_top_a = _point(support_center[0] - 0.3, support_center[1], support_center[2])
    support_top_b = _point(support_center[0] + 0.3, support_center[1], support_center[2])
    support_column_b = _point(support_center[0] + 0.3, support_center[1], PLATFORM_Z)

    toolpath = ToolpathMachine(
        points=[
            support_column_a, support_top_a, support_top_b, support_column_b,
            tested_point, _point(1.0, 0.0, 5.0),
        ]
    )
    report = support_check_run(
        toolpath, _machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    events = {e.point_index: e for e in report.events}
    assert 4 in events
    assert events[4].supported_by_tilt is False
    assert report.true_overhang_count == 1
    assert report.tilt_resolved_count == 0


def test_direction_density_matches_segmentation_default():
    """
    Coherence inter-fichiers : la densite de directions testees par defaut
    doit rester alignee avec celle de la segmentation (step3), pas une
    nouvelle constante inventee independamment.
    """
    assert SupportCheckParams().direction_candidate_count == SegmentationParams().direction_candidate_count


def test_ordering_defect_unaffected_by_tilt_check():
    """
    Un point deja classe defaut d'ordre (support vertical futur) n'est
    jamais reevalue par la passe multi-directionnelle - `has_future_support`
    reste la seule source de verite pour ce cas, `supported_by_tilt` reste
    False (categories mutuellement exclusives).
    """
    toolpath = ToolpathMachine(
        points=[
            _point(0.0, 0.0, 5.0), _point(1.0, 0.0, 5.0),   # depose en l'air
            _point(0.0, 0.0, 4.9), _point(1.0, 0.0, 4.9),   # son support, APRES
        ]
    )
    report = support_check_run(
        toolpath, _tilting_machine(), LAYER_HEIGHT, PLATFORM_Z, SupportCheckParams()
    )

    ordering_events = [e for e in report.events if e.has_future_support]
    assert len(ordering_events) == 2
    assert all(not e.supported_by_tilt for e in ordering_events)
    assert report.ordering_defect_count == 2


def test_support_test_matches_step9_distance_primitive():
    """
    Le test de support decompose la distance en composantes verticale et
    horizontale, mais la projection sur le segment doit rester EXACTEMENT
    celle de `step9._point_to_segment_distance` - non-regression contre une
    seconde geometrie qui divergerait silencieusement.
    """
    index = _PrintedGeometryIndex(0.4, 1.0)
    start = np.array([0.0, 0.0, 0.0])
    end = np.array([10.0, 0.0, 0.0])
    index.append_segment(start, end)

    # Point a 0.5 mm au-dessus du MILIEU du segment : la projection clampee
    # tombe a l'interieur, la distance vaut donc exactement 0.5.
    probe = np.array([5.0, 0.0, 0.5])
    assert _point_to_segment_distance(probe, start, end) == pytest.approx(0.5)
    assert _has_support_below(index, probe, horizontal_radius=1.0, vertical_height=1.5)

    # Meme point, mais une hauteur de recherche trop courte pour l'atteindre.
    assert not _has_support_below(index, probe, horizontal_radius=1.0, vertical_height=0.25)
