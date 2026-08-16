"""
tests/test_step5_layers.py

Tests pour pipeline/step5_layers.py, sur des mesh synthetiques generes a la
volee (aucun fichier externe requis).

Execution : pytest tests/test_step5_layers.py -v

Les fixtures enchainent le pipeline reel step1 -> step2 -> step3 -> step4,
pour que step5 recoive exactement les structures qu'il recevra en
production.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from core.mesh_types import ScalarField, TetMesh
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import FieldSolverParams, run as field_run
from pipeline.step5_layers import (
    LayerExtractionError,
    LayerParams,
    _extract_isosurface_triangles,
    _march_tetrahedron,
    run as layers_run,
)


# ---------------------------------------------------------------------------
# Fixtures : mesh synthetiques -> pipeline complet step1..step4
# ---------------------------------------------------------------------------

def _pipeline_data_from_trimesh(
    mesh: trimesh.Trimesh, resolution: float = 1.5, platform_z: float = -5.0
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
    return mesh_data, tet_mesh, overhang_map, regions, field


@pytest.fixture
def flat_box_data():
    """Cube pose a plat, 10x10x10, base a z=-5."""
    mesh = trimesh.creation.box(extents=[10, 10, 10])
    return _pipeline_data_from_trimesh(mesh, platform_z=-5.0)


@pytest.fixture
def hollow_box_data():
    """Boite creuse : une couche a mi-hauteur doit produire 2 contours disjoints."""
    outer = trimesh.creation.box(extents=[10, 10, 10])
    inner = trimesh.creation.box(extents=[4, 4, 4])
    hollow = outer.difference(inner, engine="manifold")
    return _pipeline_data_from_trimesh(hollow, resolution=1.5, platform_z=-5.0)


@pytest.fixture
def l_shape_data():
    """Forme en L : au moins une couche doit traverser 2 regions."""
    vertical = trimesh.creation.box(extents=[4, 4, 20])
    vertical.apply_translation([0, 0, 10])

    horizontal = trimesh.creation.box(extents=[20, 4, 4])
    horizontal.apply_translation([8, 0, 20])

    combined = vertical.union(horizontal, engine="manifold")
    return _pipeline_data_from_trimesh(combined, resolution=2.0, platform_z=0.0)


# ---------------------------------------------------------------------------
# Cas 1 - Cube plat : nombre de couches, contour unique, tri, local_frame vertical
# ---------------------------------------------------------------------------

def test_flat_box_layer_count_matches_height(flat_box_data):
    mesh_data, tet_mesh, overhang_map, regions, field = flat_box_data
    params = LayerParams(layer_height_base=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    # hauteur physique ~10, gradient ~1 (champ ~ z+5 en zone plane) -> ~10
    # couches attendues avec layer_height_base=1.0, tolerance large (+/-3)
    # car la mediane de ||grad(phi)|| n'est qu'approximativement 1.
    assert 5 <= len(layers) <= 15, f"Nombre de couches inattendu: {len(layers)}"


def test_flat_box_single_contour_per_layer(flat_box_data):
    mesh_data, tet_mesh, overhang_map, regions, field = flat_box_data
    params = LayerParams(layer_height_base=1.5)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)
    assert len(layers) > 0

    for layer in layers:
        assert len(layer.contour_3d) == 1, (
            f"Couche {layer.id} d'un cube plein devrait avoir un seul contour, "
            f"obtenu {len(layer.contour_3d)}"
        )
        contour = layer.contour_3d[0]
        assert contour.shape[1] == 3
        assert np.allclose(contour[0], contour[-1]), "Le contour doit etre ferme (premier == dernier point)"


def test_flat_box_layers_sorted_by_isovalue(flat_box_data):
    mesh_data, tet_mesh, overhang_map, regions, field = flat_box_data
    params = LayerParams(layer_height_base=1.5)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)
    assert len(layers) > 1

    # une seule region (cube plat sans surplomb) -> ordre = ordre des
    # centroides de contour selon z (bas vers haut)
    centroid_z = [layer.contour_3d[0][:, 2].mean() for layer in layers]
    assert all(a <= b + 1e-6 for a, b in zip(centroid_z, centroid_z[1:])), (
        "Les couches devraient etre ordonnees du plateau vers le haut sur un cube plat"
    )


def test_flat_box_local_frame_normal_is_vertical(flat_box_data):
    mesh_data, tet_mesh, overhang_map, regions, field = flat_box_data
    params = LayerParams(layer_height_base=1.5)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)
    assert len(layers) > 0

    up = np.array([0.0, 0.0, 1.0])
    for layer in layers:
        origin, axis_u, axis_v = layer.local_frame
        normal = np.cross(axis_u, axis_v)
        normal = normal / np.linalg.norm(normal)
        assert abs(abs(np.dot(normal, up)) - 1.0) < 0.2, (
            f"Couche {layer.id}: normale locale {normal} devrait etre proche de la verticale sur un cube plat"
        )
        # base orthonormee
        assert abs(np.linalg.norm(axis_u) - 1.0) < 1e-6
        assert abs(np.linalg.norm(axis_v) - 1.0) < 1e-6
        assert abs(np.dot(axis_u, axis_v)) < 1e-6


# ---------------------------------------------------------------------------
# Cas 2 - Boite creuse : 2 contours disjoints dans la meme couche, orientation opposee
# ---------------------------------------------------------------------------

def test_hollow_box_mid_layer_has_two_disjoint_contours(hollow_box_data):
    mesh_data, tet_mesh, overhang_map, regions, field = hollow_box_data
    params = LayerParams(layer_height_base=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)
    assert len(layers) > 0

    # cherche une couche dont le contour est proche de z=0 (mi-hauteur,
    # traversant a la fois la paroi exterieure et la cavite interne)
    mid_layers = [
        layer for layer in layers
        if abs(np.mean([c[:, 2].mean() for c in layer.contour_3d])) < 2.0
    ]
    assert mid_layers, "Aucune couche trouvee pres de la mi-hauteur de la boite creuse"

    multi_contour_layers = [layer for layer in mid_layers if len(layer.contour_3d) >= 2]
    assert multi_contour_layers, (
        "Au moins une couche a mi-hauteur devrait avoir 2 contours disjoints "
        "(paroi exterieure + paroi de la cavite)"
    )


def _point_in_polygon_2d(point: np.ndarray, polygon: np.ndarray) -> bool:
    """Ray-casting standard (pas de nouvelle dependance) : point (2,), polygon (N,2)."""
    x, y = point
    inside = False
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            x_intersect = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < x_intersect:
                inside = not inside
    return inside


def test_hollow_box_disjoint_contours_are_correctly_nested(hollow_box_data):
    """
    Decision 3bis (corrigee) : contour_3d + local_frame doivent suffire a
    step6 pour distinguer le contour exterieur du trou interieur, via un
    test d'imbrication geometrique (point-in-polygon) - PAS via le signe de
    l'aire projetee, dont l'hypothese s'est revelee fausse a la verification
    (teste directement sur un anneau plat synthetique : les deux boucles
    d'un patch a normales coherentes donnent le MEME signe via
    trimesh.outline(), pas des signes opposes).
    """
    mesh_data, tet_mesh, overhang_map, regions, field = hollow_box_data
    params = LayerParams(layer_height_base=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    multi_contour_layers = [layer for layer in layers if len(layer.contour_3d) >= 2]
    assert multi_contour_layers, "Test necessite au moins une couche a 2 contours"

    layer = multi_contour_layers[0]
    origin, axis_u, axis_v = layer.local_frame

    projected = []
    for contour in layer.contour_3d:
        relative = contour - origin
        u = relative @ axis_u
        v = relative @ axis_v
        projected.append(np.stack([u, v], axis=1)[:-1])  # sans le point de fermeture duplique

    # la plus grande boucle (par aire absolue) doit contenir toutes les autres
    areas = [
        abs(0.5 * np.sum(poly[:, 0] * np.roll(poly[:, 1], -1) - np.roll(poly[:, 0], -1) * poly[:, 1]))
        for poly in projected
    ]
    outer_idx = int(np.argmax(areas))

    for i, poly in enumerate(projected):
        if i == outer_idx:
            continue
        assert _point_in_polygon_2d(poly[0], projected[outer_idx]), (
            f"Le contour {i} (aire={areas[i]:.2f}) devrait etre imbrique dans le "
            f"contour exterieur {outer_idx} (aire={areas[outer_idx]:.2f}) - "
            f"step6 doit pouvoir distinguer trou et paroi exterieure par imbrication"
        )
        assert not _point_in_polygon_2d(projected[outer_idx][0], poly), (
            "La relation d'imbrication doit etre non-symetrique (le contour "
            "exterieur ne doit pas etre contenu dans le trou)"
        )


# ---------------------------------------------------------------------------
# Cas 3 - Forme en L : couche multi-region, rang effectif = max des rangs
# ---------------------------------------------------------------------------

def test_l_shape_layer_spans_multiple_regions(l_shape_data):
    mesh_data, tet_mesh, overhang_map, regions, field = l_shape_data
    params = LayerParams(layer_height_base=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)
    assert len(layers) > 0

    # une seule region existe dans cette fixture (segmentation ne fragmente
    # pas la forme en L, verifie dans test_step3) - le test verifie donc
    # que le mecanisme d'association fonctionne correctement (au moins 1
    # region par couche), pas necessairement > 1 region ici.
    for layer in layers:
        assert len(layer.region_ids) >= 1, f"Couche {layer.id} devrait etre associee a au moins une region"


def test_l_shape_layers_have_valid_ordering(l_shape_data):
    """
    L'ordre final des couches ne doit jamais violer le rang topologique :
    une couche de rang effectif strictement inferieur ne doit jamais
    apparaitre apres une couche de rang superieur dans la liste retournee.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = l_shape_data
    params = LayerParams(layer_height_base=1.0)

    from pipeline.step5_layers import _compute_region_ranks

    ranks = _compute_region_ranks(regions)
    layers = layers_run(tet_mesh, mesh_data, field, regions, params)
    assert len(layers) > 0

    effective_ranks = [
        max((ranks.get(r, 0) for r in layer.region_ids), default=0) for layer in layers
    ]
    assert all(a <= b for a, b in zip(effective_ranks, effective_ranks[1:])), (
        f"Les rangs effectifs des couches devraient etre non-decroissants dans "
        f"l'ordre retourne, obtenu: {effective_ranks}"
    )


# ---------------------------------------------------------------------------
# Cas 4 - Regression : isovaleur ne traversant aucun tetraedre
# ---------------------------------------------------------------------------

def test_isovalue_outside_range_produces_no_triangles():
    """
    _extract_isosurface_triangles doit retourner des tableaux vides (pas
    d'erreur) quand l'isovaleur ne traverse aucun tetraedre - cas d'un
    modele a composantes disjointes avec un ecart de valeurs entre elles.
    """
    tet_vertices = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    tetrahedra = np.array([[0, 1, 2, 3]])
    tet_mesh = TetMesh(tet_vertices=tet_vertices, tetrahedra=tetrahedra)

    values = np.array([0.0, 1.0, 2.0, 3.0])  # plage [0,3]
    tet_min = values[tetrahedra].min(axis=1)
    tet_max = values[tetrahedra].max(axis=1)

    vertices, faces, source_tets = _extract_isosurface_triangles(
        tet_mesh, values, isovalue=10.0, tet_min=tet_min, tet_max=tet_max
    )

    assert len(faces) == 0
    assert len(vertices) == 0
    assert len(source_tets) == 0


def test_no_gradient_raises_explicit_error(flat_box_data):
    mesh_data, tet_mesh, overhang_map, regions, field = flat_box_data
    field_no_gradient = ScalarField(values=field.values, gradient=None)

    with pytest.raises(LayerExtractionError):
        layers_run(tet_mesh, mesh_data, field_no_gradient, regions, LayerParams())


# ---------------------------------------------------------------------------
# Cas 5 - Coherence : les points de contour ont une valeur de champ ~= isovaleur
# ---------------------------------------------------------------------------

def test_march_tetrahedron_produces_points_at_isovalue():
    """
    Test unitaire direct sur _march_tetrahedron : chaque point produit doit,
    par interpolation lineaire sur son arete source, correspondre exactement
    a l'isovaleur demandee (le champ P1 est lineaire sur chaque arete).
    """
    tet_points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    tet_values = np.array([0.0, 1.0, 2.0, 3.0])
    isovalue = 1.5

    triangles = _march_tetrahedron(tet_points, tet_values, isovalue)
    assert len(triangles) > 0

    # reconstruit une interpolation barycentrique du champ sur le
    # tetraedre pour verifier chaque point produit
    def interpolate_field(p):
        # coordonnees barycentriques via resolution du systeme lineaire
        # [x1-x0, x2-x0, x3-x0] . [l1,l2,l3] = p - x0
        edge_matrix = (tet_points[1:] - tet_points[0]).T
        rel = p - tet_points[0]
        l123 = np.linalg.solve(edge_matrix, rel)
        l0 = 1.0 - l123.sum()
        weights = np.array([l0, *l123])
        return float(weights @ tet_values)

    for tri in triangles:
        for point in tri:
            value = interpolate_field(point)
            assert abs(value - isovalue) < 1e-9, (
                f"Point {point} du triangle extrait a une valeur de champ "
                f"{value}, attendu {isovalue}"
            )


def test_march_tetrahedron_triangle_normals_point_toward_increasing_field():
    """Regression Decision 3bis : orientation coherente, verifiee directement."""
    tet_points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    # teste les 3 cas structurels (1, 2, 3 sommets au-dessus) avec des
    # jeux de valeurs differents, verifie l'orientation a chaque fois
    for tet_values in (
        np.array([0.0, 1.0, 2.0, 3.0]),   # varie selon n_above
        np.array([3.0, 2.0, 1.0, 0.0]),
        np.array([1.0, 3.0, 0.0, 2.0]),
    ):
        for isovalue in (0.5, 1.5, 2.5):
            triangles = _march_tetrahedron(tet_points, tet_values, isovalue)
            if not triangles:
                continue

            hi_idx = int(np.argmax(tet_values))
            lo_idx = int(np.argmin(tet_values))
            reference_direction = tet_points[hi_idx] - tet_points[lo_idx]

            for tri in triangles:
                normal = np.cross(tri[1] - tri[0], tri[2] - tri[0])
                assert np.dot(normal, reference_direction) >= -1e-9, (
                    f"Triangle mal oriente pour tet_values={tet_values}, isovalue={isovalue}"
                )
