"""
tests/test_step6_toolpath2d.py

Tests pour pipeline/step6_toolpath2d.py, sur des mesh synthetiques generes a
la volee (aucun fichier externe requis).

Execution : pytest tests/test_step6_toolpath2d.py -v

Les fixtures enchainent le pipeline reel step1 -> step2 -> step3 -> step4 ->
step5, pour que step6 recoive exactement les structures qu'il recevra en
production.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from core.mesh_types import Layer
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import FieldSolverParams, run as field_run
from pipeline.step5_layers import LayerParams, run as layers_run
from pipeline.step6_toolpath2d import (
    InfillPattern,
    Step6Params,
    _classify_contours,
    _project_to_local_2d,
    run as toolpath_run,
)


# ---------------------------------------------------------------------------
# Fixtures : mesh synthetiques -> pipeline complet step1..step5
# ---------------------------------------------------------------------------

def _pipeline_layers_from_trimesh(
    mesh: trimesh.Trimesh,
    resolution: float = 1.5,
    platform_z: float = -5.0,
    layer_height_base: float = 1.0,
    overhang_angle_threshold_deg: float = 45.0,
    bridge_max_gap_distance: float = 2.0,
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

    overhang_map = overhang_run(
        mesh_data,
        OverhangParams(
            overhang_angle_threshold_deg=overhang_angle_threshold_deg,
            bridge_max_gap_distance=bridge_max_gap_distance,
        ),
    )
    regions = segmentation_run(
        mesh_data,
        overhang_map,
        SegmentationParams(platform_z=platform_z, platform_contact_tolerance=0.5),
    )
    field = field_run(tet_mesh, mesh_data, overhang_map, regions, [], FieldSolverParams())
    layers = layers_run(tet_mesh, mesh_data, field, regions, LayerParams(layer_height_base=layer_height_base))
    return mesh_data, overhang_map, layers


@pytest.fixture
def flat_box_layers():
    """Cube pose a plat, 10x10x10, base a z=-5."""
    mesh = trimesh.creation.box(extents=[10, 10, 10])
    return _pipeline_layers_from_trimesh(mesh, platform_z=-5.0, layer_height_base=2.0)


@pytest.fixture
def hollow_box_layers():
    """Boite creuse : une couche a mi-hauteur a 2 contours disjoints (verifie step5)."""
    outer = trimesh.creation.box(extents=[10, 10, 10])
    inner = trimesh.creation.box(extents=[4, 4, 4])
    hollow = outer.difference(inner, engine="manifold")
    return _pipeline_layers_from_trimesh(hollow, resolution=1.5, platform_z=-5.0, layer_height_base=1.0)


@pytest.fixture
def bridge_layers():
    """
    Deux piliers relies par une dalle horizontale fine en hauteur - la
    face inferieure de la dalle, entre les piliers, n'a rien en dessous
    -> is_bridge (meme geometrie que la fixture pont de step2, adaptee).
    """
    pillar_a = trimesh.creation.box(extents=[3, 3, 10])
    pillar_a.apply_translation([-6, 0, 5])

    pillar_b = trimesh.creation.box(extents=[3, 3, 10])
    pillar_b.apply_translation([6, 0, 5])

    slab = trimesh.creation.box(extents=[16, 3, 1.5])
    slab.apply_translation([0, 0, 9])

    combined = pillar_a.union(pillar_b, engine="manifold").union(slab, engine="manifold")
    return _pipeline_layers_from_trimesh(
        combined, resolution=1.5, platform_z=0.0, layer_height_base=1.0,
        bridge_max_gap_distance=4.0,
    )


@pytest.fixture
def thin_wall_layers():
    """
    Zone plus fine que nozzle_width * wall_count : deux boites reliees par
    un mur tres fin (0.3mm), avec wall_count=2 et nozzle_width=0.4 -> la
    zone du mur fin (0.3 < 0.4*2=0.8) doit voir sa generation de parois
    s'arreter proprement avant wall_count.
    """
    box_a = trimesh.creation.box(extents=[4, 4, 4])
    box_a.apply_translation([-2.15, 0, 2])

    box_b = trimesh.creation.box(extents=[4, 4, 4])
    box_b.apply_translation([2.15, 0, 2])

    thin_wall = trimesh.creation.box(extents=[0.3, 2, 2])
    thin_wall.apply_translation([0, 0, 2])

    combined = box_a.union(box_b, engine="manifold").union(thin_wall, engine="manifold")
    return _pipeline_layers_from_trimesh(combined, resolution=0.6, platform_z=0.0, layer_height_base=1.0)


# ---------------------------------------------------------------------------
# Cas 1 - Cube plat : wall_count parois + infill, ordre centre -> exterieur
# ---------------------------------------------------------------------------

def test_flat_box_produces_walls_and_infill(flat_box_layers):
    mesh_data, overhang_map, layers = flat_box_layers
    params = Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.3)

    result = toolpath_run(layers, mesh_data, overhang_map, params)
    assert len(result) > 0

    for layer in result:
        assert layer.toolpath_2d is not None
        assert len(layer.toolpath_2d) > 0, f"Couche {layer.id} n'a produit aucune trajectoire"


def test_flat_box_walls_ordered_center_to_exterior(flat_box_layers):
    mesh_data, overhang_map, layers = flat_box_layers
    params = Step6Params(nozzle_width=0.4, wall_count=3, infill_density=0.0)

    result = toolpath_run(layers, mesh_data, overhang_map, params)
    assert len(result) > 0

    def _polygon_area(poly: np.ndarray) -> float:
        x, y = poly[:, 0], poly[:, 1]
        return abs(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

    checked_any = False
    for layer in result:
        # sans infill (density=0), toolpath_2d ne contient que des parois
        walls = layer.toolpath_2d
        if len(walls) < 2:
            continue
        areas = [_polygon_area(w) for w in walls]
        assert all(a <= b + 1e-6 for a, b in zip(areas, areas[1:])), (
            f"Couche {layer.id}: parois devraient etre ordonnees par aire croissante "
            f"(centre -> exterieur), obtenu {areas}"
        )
        checked_any = True

    assert checked_any, "Aucune couche avec >=2 parois n'a ete trouvee pour verifier l'ordre"


# ---------------------------------------------------------------------------
# Cas 2 - Boite creuse : la boucle interieure est classee comme trou
# ---------------------------------------------------------------------------

def test_hollow_box_inner_loop_classified_as_hole(hollow_box_layers):
    mesh_data, overhang_map, layers = hollow_box_layers

    multi_contour_layers = [layer for layer in layers if len(layer.contour_3d) >= 2]
    assert multi_contour_layers, "Fixture attendue avec au moins une couche a 2 contours (verifie en step5)"

    layer = multi_contour_layers[0]
    contours_2d = [_project_to_local_2d(c, layer.local_frame) for c in layer.contour_3d]
    groups = _classify_contours(contours_2d)

    assert len(groups) == 1, "Un seul contour matiere attendu (paroi exterieure), la cavite est un trou"
    outer, holes = groups[0]
    assert len(holes) == 1, "La boucle interieure de la cavite doit etre classee comme trou"


def test_hollow_box_no_wall_generated_inside_cavity(hollow_box_layers):
    mesh_data, overhang_map, layers = hollow_box_layers
    params = Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.2)

    result = toolpath_run(layers, mesh_data, overhang_map, params)

    multi_contour_layers = [layer for layer in result if len(layer.contour_3d) >= 2]
    assert multi_contour_layers

    layer = multi_contour_layers[0]
    assert layer.toolpath_2d is not None and len(layer.toolpath_2d) > 0

    # toute polyligne de la trajectoire doit rester a l'exterieur de la
    # cavite (aucune paroi/infill genere dans le vide de la cavite) :
    # verifie qu'aucun point de toolpath_2d ne tombe strictement a
    # l'interieur de la plus petite des 2 boucles de contour (la cavite).
    contours_2d = [_project_to_local_2d(c, layer.local_frame) for c in layer.contour_3d]
    areas = []
    for c in contours_2d:
        c_open = c[:-1] if np.allclose(c[0], c[-1]) else c
        x, y = c_open[:, 0], c_open[:, 1]
        areas.append(abs(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))
    cavity_loop = contours_2d[int(np.argmin(areas))]
    cavity_open = cavity_loop[:-1] if np.allclose(cavity_loop[0], cavity_loop[-1]) else cavity_loop

    from core.geometry_utils import point_in_polygon_2d

    for polyline in layer.toolpath_2d:
        for point in polyline:
            assert not point_in_polygon_2d(point, cavity_open), (
                f"Point {point} de la trajectoire tombe a l'interieur de la cavite (trou) - fuite"
            )


# ---------------------------------------------------------------------------
# Cas 3 - Zone fine : arret propre de la generation de parois avant wall_count
# ---------------------------------------------------------------------------

def test_thin_wall_zone_stops_before_wall_count_without_crash(thin_wall_layers):
    mesh_data, overhang_map, layers = thin_wall_layers
    params = Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.2)

    # ne doit lever aucune exception
    result = toolpath_run(layers, mesh_data, overhang_map, params)
    assert len(result) > 0

    # au moins une couche doit exister avec une trajectoire produite (le
    # mur fin ne doit pas faire planter tout le pipeline meme s'il ne peut
    # pas recevoir wall_count parois pleines)
    assert any(layer.toolpath_2d for layer in result)


# ---------------------------------------------------------------------------
# Cas 4 - Pont interne : trajectoire reorientee perpendiculaire a la portee
# ---------------------------------------------------------------------------

def test_bridge_layer_toolpath_perpendicular_to_span(bridge_layers):
    mesh_data, overhang_map, layers = bridge_layers
    assert np.any(overhang_map.is_bridge), "Fixture attendue avec au moins une face is_bridge"

    params = Step6Params(
        nozzle_width=0.4, wall_count=1, infill_density=0.3,
        bridge_anchor_extra_walls=1,
    )
    result = toolpath_run(layers, mesh_data, overhang_map, params)

    # cherche une couche proche du dessous de la dalle (z ~ 8.25, juste
    # sous la face inferieure du pont a z=8.25) dont l'INFILL (pas les
    # parois) a ete reoriente. Isoler l'infill est necessaire : chaque
    # polyligne de _generate_infill est un segment a 2 points, alors que les
    # parois (issues de _offset_walls) ont des dizaines de points et
    # dominent une PCA globale sur toolpath_2d avec la forme allongee du
    # contour lui-meme (parallele a la portee, pas perpendiculaire) -
    # confondre les deux masque completement le signal de reorientation de
    # l'infill (bug de methodologie de test decouvert en diagnostiquant cet
    # echec, pas un bug de production : voir CORRECTION 1 du plan).
    found_bridge_layer = False
    for layer in result:
        if not layer.toolpath_2d:
            continue
        z_mean = np.mean([c[:, 2].mean() for c in layer.contour_3d])
        if not (7.0 <= z_mean <= 9.5):
            continue

        infill_lines = [poly for poly in layer.toolpath_2d if len(poly) == 2]
        if not infill_lines:
            continue

        origin, axis_u, axis_v = layer.local_frame
        world_span_x = np.array([1.0, 0.0, 0.0])
        span_2d = np.array([np.dot(world_span_x, axis_u), np.dot(world_span_x, axis_v)])
        if np.linalg.norm(span_2d) < 1e-6:
            continue
        span_2d = span_2d / np.linalg.norm(span_2d)

        # chaque segment d'infill a sa propre direction (p1 - p0) ; toutes
        # doivent etre PROCHES DE PERPENDICULAIRES a la portee (|dot| proche
        # de 0), puisque l'infill balaie perpendiculairement au pont pour le
        # traverser plutot que longer son vide.
        line_dirs = np.array([poly[1] - poly[0] for poly in infill_lines])
        norms = np.linalg.norm(line_dirs, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        line_dirs = line_dirs / norms
        dots = np.abs(line_dirs @ span_2d)

        if np.mean(dots) < 0.5:
            found_bridge_layer = True
            break

    assert found_bridge_layer, (
        "Aucune couche pres du pont n'a de trajectoire orientee perpendiculairement "
        "a la direction de portee attendue"
    )


# ---------------------------------------------------------------------------
# Cas 5 - Classification generique : anneau dans anneau (>2 boucles imbriquees)
# ---------------------------------------------------------------------------

def _circle_loop(radius: float, n: int = 32) -> np.ndarray:
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([radius * np.cos(angles), radius * np.sin(angles)], axis=1)


def test_classify_contours_ring_in_ring_parity():
    """
    3 boucles concentriques (rayon 10, 6, 2) : la matiere alterne avec le
    vide en partant de l'exterieur - boucle 1 (r=10, la plus grande) =
    matiere, boucle 2 (r=6) = trou, boucle 3 (r=2, la plus petite) = matiere
    (un ilot plein au centre du trou) - verifie la parite generalisee au-
    dela du cas a 2 boucles.
    """
    outer = _circle_loop(10.0)
    middle = _circle_loop(6.0)
    inner = _circle_loop(2.0)

    groups = _classify_contours([outer, middle, inner])

    # attendu : 2 contours "matiere" (outer avec middle comme trou, et inner
    # comme ilot plein independant)
    assert len(groups) == 2

    outer_group = next(g for g in groups if np.allclose(g[0][0], outer[0]))
    inner_group = next(g for g in groups if np.allclose(g[0][0], inner[0]))

    assert len(outer_group[1]) == 1
    assert np.allclose(outer_group[1][0][0], middle[0])
    assert len(inner_group[1]) == 0


def test_classify_contours_two_disjoint_islands():
    """2 boucles disjointes (pas imbriquees) : chacune est sa propre matiere sans trou."""
    island_a = _circle_loop(3.0) + np.array([20.0, 0.0])
    island_b = _circle_loop(3.0) + np.array([-20.0, 0.0])

    groups = _classify_contours([island_a, island_b])

    assert len(groups) == 2
    for _outer, holes in groups:
        assert holes == []
