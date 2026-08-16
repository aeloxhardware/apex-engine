"""
tests/test_step7_scheduling.py

Tests pour pipeline/step7_scheduling.py, sur des mesh synthetiques generes a
la volee (aucun fichier externe requis).

Execution : pytest tests/test_step7_scheduling.py -v

Les fixtures enchainent le pipeline reel step1 -> ... -> step6, pour que
step7 recoive exactement les structures qu'il recevra en production.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from core.mesh_types import ColorMap
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import FieldSolverParams, run as field_run
from pipeline.step5_layers import LayerParams, run as layers_run
from pipeline.step6_toolpath2d import Step6Params, run as toolpath_run
from pipeline.step7_scheduling import (
    HoleCostModel,
    PurgeCostModel,
    SchedulingError,
    Step7Params,
    _color_runs,
    _smooth_edge_colors,
    run as scheduling_run,
)


# ---------------------------------------------------------------------------
# Fixtures : mesh synthetiques -> pipeline complet step1..step6
# ---------------------------------------------------------------------------

def _pipeline_through_step6(
    mesh: trimesh.Trimesh,
    resolution: float = 1.5,
    platform_z: float = -5.0,
    layer_height_base: float = 2.0,
):
    import tempfile
    import os

    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, tet_mesh, _ = ingest_run(path, resolution=resolution)
    finally:
        os.unlink(path)

    overhang_map = overhang_run(mesh_data, OverhangParams(overhang_angle_threshold_deg=45.0))
    regions = segmentation_run(
        mesh_data,
        overhang_map,
        SegmentationParams(platform_z=platform_z, platform_contact_tolerance=0.5),
    )
    field = field_run(tet_mesh, mesh_data, overhang_map, regions, [], FieldSolverParams())
    layers = layers_run(
        tet_mesh, mesh_data, field, regions, LayerParams(layer_height_base=layer_height_base)
    )
    layers = toolpath_run(
        layers, mesh_data, overhang_map,
        Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.3),
    )
    return mesh_data, tet_mesh, regions, field, layers


@pytest.fixture
def flat_box_pipeline():
    """
    Cube pose a plat, 10x10x10, base a z=-5, SUBDIVISE.

    La subdivision est necessaire, pas cosmetique : un cube trimesh brut n'a
    que 12 triangles, dont les faces laterales s'etendent sur toute la
    hauteur (centres a z ~ +/-1.67). Comme `ColorMap.per_face_color` associe
    une couleur ENTIERE a chaque face, une peinture definie par une frontiere
    spatiale (z<0 vs z>0) ne peut tout simplement pas etre representee a
    cette resolution : un point a z=-3.99 herite de la couleur du triangle
    qui le contient, ce qui donne une decoupe diagonale et non horizontale.
    Subdiviser rend les faces assez petites pour que la peinture par face
    approxime reellement la frontiere geometrique.
    """
    mesh = trimesh.creation.box(extents=[10, 10, 10])
    mesh = mesh.subdivide().subdivide()
    return _pipeline_through_step6(mesh, platform_z=-5.0)


def _two_tone_color_map(mesh_data, split_z: float = 0.0, merge_tolerance: float = 0.0) -> ColorMap:
    """
    Peinture bicolore VERTICALE : couleur 0 sous `split_z`, couleur 1
    au-dessus. Chaque couche horizontale ne voit donc qu'une seule couleur -
    utile pour tester la projection, PAS les trous temporaires (voir
    `_lateral_color_map`).
    """
    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    per_face = {
        int(face_id): (0 if center_z < split_z else 1)
        for face_id, center_z in enumerate(face_centers[:, 2])
    }
    return ColorMap(per_face_color=per_face, merge_tolerance=merge_tolerance)


def _lateral_color_map(mesh_data, split_x: float = 0.0, merge_tolerance: float = 0.0) -> ColorMap:
    """
    Peinture bicolore LATERALE (selon X) : les deux couleurs coexistent dans
    CHAQUE couche horizontale.

    C'est la seule configuration qui peut declencher un trou temporaire : un
    trou n'a de sens que si deux couleurs se disputent la meme couche (on en
    imprime une, on differe l'autre). Une decoupe alignee sur l'axe de
    tranchage (`_two_tone_color_map`) ne produit jamais ce conflit - mesure
    directe : une seule couleur par couche, donc zero candidat au report.
    """
    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    per_face = {
        int(face_id): (0 if center_x < split_x else 1)
        for face_id, center_x in enumerate(face_centers[:, 0])
    }
    return ColorMap(per_face_color=per_face, merge_tolerance=merge_tolerance)


# ---------------------------------------------------------------------------
# Cas 1 - Sans ColorMap : le pipeline nu doit rester pleinement fonctionnel
# ---------------------------------------------------------------------------

def test_no_color_map_produces_one_pass_per_layer(flat_box_pipeline):
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline

    colored, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, None, regions, field, Step7Params()
    )

    assert len(schedule.passes) == len(layers), (
        "Sans peinture, chaque couche doit donner exactement une passe"
    )
    assert all(not p.is_hole_fill for p in schedule.passes), "Aucun trou sans couleur"
    assert all(p.color_id is None for p in schedule.passes)
    assert all(layer.segment_colors is not None for layer in colored)


def test_missing_toolpath_raises_explicit_error(flat_box_pipeline):
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline

    from dataclasses import replace

    broken = [replace(layers[0], toolpath_2d=None), *layers[1:]]

    with pytest.raises(SchedulingError):
        scheduling_run(broken, mesh_data, tet_mesh, None, regions, field, Step7Params())


# ---------------------------------------------------------------------------
# Cas 2 - Projection de couleur : surface coloree, interieur libre (Question 1)
# ---------------------------------------------------------------------------

def test_interior_segments_have_no_imposed_color(flat_box_pipeline):
    """
    Question 1 : la couleur est une propriete de surface VISIBLE. Les
    segments internes (infill, parois internes) doivent recevoir None -
    matiere imprimable avec n'importe quelle couleur deja chargee.
    """
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline
    color_map = _two_tone_color_map(mesh_data)

    colored, _schedule = scheduling_run(
        layers, mesh_data, tet_mesh, color_map, regions, field,
        Step7Params(surface_color_distance=0.5),
    )

    all_colors = [c for layer in colored for c in (layer.segment_colors or [])]
    assert all_colors, "Le test necessite des segments"
    assert any(c is None for c in all_colors), (
        "Au moins un segment interne doit etre sans couleur imposee"
    )
    assert any(c is not None for c in all_colors), (
        "Au moins un segment de surface doit recevoir une couleur peinte"
    )


def test_painted_colors_match_geometry(flat_box_pipeline):
    """Les segments colores doivent recevoir la couleur de leur zone geometrique."""
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline
    color_map = _two_tone_color_map(mesh_data, split_z=0.0)

    colored, _schedule = scheduling_run(
        layers, mesh_data, tet_mesh, color_map, regions, field,
        Step7Params(surface_color_distance=0.5),
    )

    from pipeline.step6_toolpath2d import _lift_to_3d

    mismatches = 0
    checked = 0
    for layer in colored:
        for polyline, color in zip(layer.toolpath_2d, layer.segment_colors):
            if color is None:
                continue
            centroid_3d = _lift_to_3d(
                polyline.mean(axis=0).reshape(1, 2), layer.local_frame
            )[0]
            expected = 0 if centroid_3d[2] < 0.0 else 1
            checked += 1
            if color != expected:
                mismatches += 1

    assert checked > 0, "Aucun segment colore trouve"
    # tolerance : les segments a cheval sur la frontiere z=0 sont ambigus
    assert mismatches / checked < 0.25, (
        f"Trop de segments mal colores : {mismatches}/{checked}"
    )


def test_large_surface_distance_colors_everything(flat_box_pipeline):
    """
    Avec un seuil de distance tres grand, meme l'infill interne se voit
    attribuer une couleur - verifie que c'est bien la DISTANCE qui pilote la
    distinction visible/interne, et pas un autre critere implicite.
    """
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline
    color_map = _two_tone_color_map(mesh_data)

    colored, _schedule = scheduling_run(
        layers, mesh_data, tet_mesh, color_map, regions, field,
        Step7Params(surface_color_distance=1e6),
    )

    all_colors = [c for layer in colored for c in (layer.segment_colors or [])]
    assert all(c is not None for c in all_colors), (
        "Avec un seuil infini, aucun segment ne doit rester sans couleur"
    )


# ---------------------------------------------------------------------------
# Cas 3 - merge_tolerance (Decision 4)
# ---------------------------------------------------------------------------

def test_smooth_edge_colors_absorbs_short_runs():
    """Un run plus court que la tolerance est absorbe par son voisin le plus long."""
    colors = [0, 0, 0, 0, 1, 0, 0, 0, 0]
    lengths = np.ones(len(colors))

    smoothed = _smooth_edge_colors(colors, lengths, merge_tolerance=2.0)

    assert smoothed == [0] * 9, f"Le run isole de couleur 1 devait etre absorbe : {smoothed}"


def test_smooth_edge_colors_keeps_long_runs():
    """Un run plus long que la tolerance survit."""
    colors = [0, 0, 0, 0, 1, 1, 1, 1]
    lengths = np.ones(len(colors))

    smoothed = _smooth_edge_colors(colors, lengths, merge_tolerance=2.0)

    assert smoothed == colors, f"Les deux runs longs devaient survivre : {smoothed}"


def test_smooth_edge_colors_disabled_by_zero_tolerance():
    colors = [0, 1, 0, 1]
    lengths = np.ones(len(colors))

    assert _smooth_edge_colors(colors, lengths, merge_tolerance=0.0) == colors


def test_color_runs_segmentation():
    assert _color_runs([0, 0, 1, 1, 1, None]) == [(0, 2, 0), (2, 5, 1), (5, 6, None)]


def test_merge_tolerance_reduces_color_fragmentation(flat_box_pipeline):
    """Une tolerance elevee doit produire moins de segments colores distincts."""
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline

    def count_color_changes(merge_tolerance: float) -> int:
        color_map = _two_tone_color_map(mesh_data, merge_tolerance=merge_tolerance)
        colored, _ = scheduling_run(
            layers, mesh_data, tet_mesh, color_map, regions, field,
            Step7Params(surface_color_distance=0.5),
        )
        return sum(len(layer.toolpath_2d) for layer in colored)

    without = count_color_changes(0.0)
    with_tolerance = count_color_changes(50.0)

    assert with_tolerance <= without, (
        f"Le lissage devrait reduire (ou au pire conserver) la fragmentation : "
        f"{with_tolerance} vs {without}"
    )


# ---------------------------------------------------------------------------
# Cas 4 - Bascule trou vs purge (Decision 6)
# ---------------------------------------------------------------------------

def test_cheap_purge_disables_all_holes(flat_box_pipeline):
    """
    Avec une purge quasi gratuite, aucun trou temporaire ne doit etre retenu.

    Peinture laterale : sans elle le test passerait a vide (aucun trou n'est
    candidat de toute facon), et ne prouverait donc rien sur l'arbitrage.
    """
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline
    color_map = _lateral_color_map(mesh_data)

    _colored, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, color_map, regions, field,
        Step7Params(
            hole_cost_model=HoleCostModel(reheat_penalty=1000.0, adhesion_risk_weight=1000.0),
            purge_cost_model=PurgeCostModel(volume_per_purge=0.0, time_per_purge=0.0),
        ),
    )

    assert all(not p.is_hole_fill for p in schedule.passes), (
        "Une purge gratuite doit toujours battre un trou couteux"
    )


def test_max_open_holes_zero_disables_all_holes(flat_box_pipeline):
    """max_open_holes=0 interdit structurellement tout trou en attente."""
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline
    color_map = _lateral_color_map(mesh_data)

    _colored, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, color_map, regions, field,
        Step7Params(max_open_holes=0),
    )

    assert all(not p.is_hole_fill for p in schedule.passes)


def test_cheap_holes_produce_hole_fill_passes(flat_box_pipeline):
    """
    A l'inverse : trou quasi gratuit + purge tres couteuse doit produire au
    moins une passe de comblement - verifie que la branche "trou" du
    mecanisme est bien atteignable, pas seulement la branche purge.

    Peinture LATERALE obligatoire ici : il faut deux couleurs dans la MEME
    couche pour qu'un report ait lieu.
    """
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline
    color_map = _lateral_color_map(mesh_data)

    _colored, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, color_map, regions, field,
        Step7Params(
            hole_cost_model=HoleCostModel(reheat_penalty=0.0, adhesion_risk_weight=0.0),
            purge_cost_model=PurgeCostModel(volume_per_purge=1000.0, time_per_purge=1000.0),
            max_open_holes=100,
        ),
    )

    assert any(p.is_hole_fill for p in schedule.passes), (
        "Avec un trou gratuit et une purge tres couteuse, des trous temporaires "
        "doivent etre planifies"
    )


# ---------------------------------------------------------------------------
# Cas 5 - Invariant global : chaque segment planifie exactement une fois
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("use_color", [False, True])
def test_every_segment_scheduled_exactly_once(flat_box_pipeline, use_color):
    """
    Invariant le plus important pour step8 : aucun segment perdu, aucun
    duplique. Verifie avec et sans peinture.
    """
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline
    color_map = _two_tone_color_map(mesh_data) if use_color else None

    colored, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, color_map, regions, field, Step7Params()
    )

    expected = {
        (layer.id, index)
        for layer in colored
        for index in range(len(layer.toolpath_2d or []))
    }
    scheduled = [seg for p in schedule.passes for seg in p.segments]

    assert len(scheduled) == len(set(scheduled)), "Un segment est planifie plusieurs fois"
    assert set(scheduled) == expected, (
        f"Segments manquants: {expected - set(scheduled)}, "
        f"segments inconnus: {set(scheduled) - expected}"
    )


def test_pass_layer_range_matches_its_segments(flat_box_pipeline):
    """`layer_range` doit rester coherent avec `segments` (donnee derivee)."""
    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline
    color_map = _two_tone_color_map(mesh_data)

    _colored, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, color_map, regions, field, Step7Params()
    )

    for print_pass in schedule.passes:
        layer_ids = [layer_id for layer_id, _ in print_pass.segments]
        assert print_pass.layer_range == (min(layer_ids), max(layer_ids))


# ---------------------------------------------------------------------------
# Cas 6 - Deadline d'accessibilite (Question 2)
# ---------------------------------------------------------------------------

def test_deadline_is_within_layer_bounds(flat_box_pipeline):
    """
    Une deadline doit toujours rester dans [couche du trou, derniere couche] -
    borne verifiee directement sur la fonction, independamment du planning.
    """
    from pipeline.step7_scheduling import (
        _escape_directions,
        _group_position_3d,
        _group_segments_by_color,
        compute_hole_deadline,
    )
    import scipy.spatial

    mesh_data, tet_mesh, regions, field, layers = flat_box_pipeline
    color_map = _two_tone_color_map(mesh_data)

    colored, _schedule = scheduling_run(
        layers, mesh_data, tet_mesh, color_map, regions, field, Step7Params()
    )

    params = Step7Params()
    groups_per_layer = _group_segments_by_color(colored)
    field_tree = scipy.spatial.cKDTree(tet_mesh.tet_vertices)
    sorted_isovalues = np.sort(np.asarray([l.isovalue for l in colored], dtype=float))
    directions = _escape_directions(params)

    checked = 0
    for layer_index, groups in enumerate(groups_per_layer):
        for group in groups:
            position = _group_position_3d(colored[layer_index], group)
            deadline = compute_hole_deadline(
                position, layer_index, mesh_data, field_tree,
                np.asarray(field.values, dtype=float), sorted_isovalues,
                len(colored), directions, params,
            )
            assert layer_index <= deadline <= len(colored) - 1
            checked += 1

    assert checked > 0


def test_layers_carry_isovalue(flat_box_pipeline):
    """
    Regression : step5 doit renseigner Layer.isovalue (requis par step7), et
    step6 doit le PRESERVER - la fixture va jusqu'a step6, donc cette
    assertion couvre les deux etapes.
    """
    _mesh_data, _tet_mesh, _regions, _field, layers = flat_box_pipeline

    assert all(layer.isovalue is not None for layer in layers)
    isovalues = [layer.isovalue for layer in layers]
    assert len(set(isovalues)) == len(isovalues), "Les isovaleurs doivent etre distinctes"


def test_segment_colors_unset_before_step7(flat_box_pipeline):
    """
    Contrat inverse : `segment_colors` est un champ de l'etape 7 SEULEMENT.
    Il doit rester a None apres step5 et step6 (la fixture va jusqu'a step6),
    sinon c'est que l'une de ces etapes le renseigne par accident.
    """
    _mesh_data, _tet_mesh, _regions, _field, layers = flat_box_pipeline

    assert all(layer.segment_colors is None for layer in layers)
    assert all(layer.toolpath_2d is not None for layer in layers), (
        "step6 doit avoir rempli toolpath_2d"
    )
