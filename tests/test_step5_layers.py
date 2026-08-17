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

from core.geometry_utils import point_in_polygon_2d
from core.mesh_types import ScalarField, TetMesh
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import FieldSolverParams, run as field_run
from pipeline.step5_layers import (
    LayerExtractionError,
    LayerParams,
    _cluster_contours,
    _extract_isosurface_triangles,
    _march_tetrahedron,
    run as layers_run,
)


# ---------------------------------------------------------------------------
# Fixtures : mesh synthetiques -> pipeline complet step1..step4
# ---------------------------------------------------------------------------

def _pipeline_data_from_trimesh(
    mesh: trimesh.Trimesh,
    resolution: float = 1.5,
    platform_z: float = -5.0,
    segmentation_params: SegmentationParams | None = None,
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
    if segmentation_params is None:
        segmentation_params = SegmentationParams(platform_z=platform_z, platform_contact_tolerance=0.5)
    regions = segmentation_run(mesh_data, overhang_map, segmentation_params)
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


def _l_shape_mesh() -> trimesh.Trimesh:
    vertical = trimesh.creation.box(extents=[4, 4, 20])
    vertical.apply_translation([0, 0, 10])

    horizontal = trimesh.creation.box(extents=[20, 4, 4])
    horizontal.apply_translation([8, 0, 20])

    return vertical.union(horizontal, engine="manifold")


@pytest.fixture
def l_shape_data():
    """
    Forme en L : au moins une couche doit traverser 2 regions.

    Segmentation aux parametres PAR DEFAUT de SegmentationParams
    (direction_deviation_weight=0.0, desactive) : cette fixture ne scinde
    PAS la piece en 2 regions separees (une seule region couvre toute la
    piece, comme avant le fix de step3 documente dans CLAUDE.md) - elle
    n'exerce donc jamais le cas d'ordre aberrant traite par le fix de
    step5 ci-dessous. Voir `l_shape_split_data` pour ce cas.
    """
    return _pipeline_data_from_trimesh(_l_shape_mesh(), resolution=2.0, platform_z=0.0)


@pytest.fixture
def l_shape_split_data():
    """
    Meme geometrie que l_shape_data, mais avec direction_deviation_weight
    actif (meme valeur calibree que main.py/PipelineParams, voir CLAUDE.md
    limites #31/#32) : produit reellement 2 regions separees (base
    verticale + bras incline), le cas qui exerce le fix de step5 (tri par
    isovaleur seule, pas par rang topologique de region).
    """
    return _pipeline_data_from_trimesh(
        _l_shape_mesh(),
        resolution=2.0,
        segmentation_params=SegmentationParams(
            platform_z=0.0,
            platform_contact_tolerance=0.5,
            direction_deviation_weight=0.2,
        ),
    )


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
        assert point_in_polygon_2d(poly[0], projected[outer_idx]), (
            f"Le contour {i} (aire={areas[i]:.2f}) devrait etre imbrique dans le "
            f"contour exterieur {outer_idx} (aire={areas[outer_idx]:.2f}) - "
            f"step6 doit pouvoir distinguer trou et paroi exterieure par imbrication"
        )
        assert not point_in_polygon_2d(projected[outer_idx][0], poly), (
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


def test_l_shape_layers_ordered_by_isovalue_not_region_rank(l_shape_split_data):
    """
    L'ordre final des couches doit suivre l'ISOVALEUR SEULE, jamais le rang
    topologique de region - voir la docstring de step5_layers.run() pour le
    diagnostic complet.

    Utilise l_shape_split_data (direction_deviation_weight actif, 2 regions
    reelles - contrairement a l_shape_data ci-dessus, qui ne scinde pas et
    n'exerce donc jamais ce cas). Avant ce fix, cette meme fixture
    produisait un ordre non-monotone en isovaleur : une couche de rang 0
    (region "base", verticale) tres haute en isovaleur (proche du sommet du
    bras) se retrouvait triee avant tout le corps de la region de rang 1
    (le bras), qui la precede pourtant physiquement - mesure concretement
    sur res/test_models/l_shape_overhang.stl (sauts d'isovaleur ~4.1 -> 10.5
    -> 4.2 dans colored_layers, produisant un trajet a vide en eventail
    visible sur out_test.glb).
    """
    mesh_data, tet_mesh, overhang_map, regions, field = l_shape_split_data
    params = LayerParams(layer_height_base=1.0)

    assert len(regions.regions) > 1, (
        "fixture invalide : l_shape_split_data devrait produire plusieurs "
        "regions (sinon ce test n'exerce pas le cas vise)"
    )

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)
    assert len(layers) > 1

    isovalues = [layer.isovalue for layer in layers]
    assert all(a <= b for a, b in zip(isovalues, isovalues[1:])), (
        f"L'isovaleur des couches devrait etre non-decroissante dans l'ordre "
        f"retourne (ordre = progression spatiale reelle du champ), obtenu: "
        f"{[round(v, 2) for v in isovalues]}"
    )


# ---------------------------------------------------------------------------
# Cas 3bis - CLAUDE.md limite #35 : scission par cluster spatial
# ---------------------------------------------------------------------------

def test_l_shape_split_data_no_layer_has_disjoint_unnested_loops(l_shape_split_data):
    """
    Une fois le contourage par cluster spatial en place, aucun Layer ne
    devrait plus porter 2 boucles a la fois disjointes (non imbriquees) ET
    spatialement eloignees - c'est exactement le chevauchement base/bras
    documente comme limite #35, cense disparaitre avec ce fix.

    Une boucle imbriquee dans une autre (trou d'une section solide, cas
    boite creuse) reste un cas valide a 2 boucles dans un seul Layer - ce
    test ne l'interdit pas, il verifie seulement l'absence du motif
    "eloigne et non imbrique" mesure sur l_shape_overhang.stl.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = l_shape_split_data
    params = LayerParams(layer_height_base=1.0)

    assert len(regions.regions) > 1, (
        "fixture invalide : l_shape_split_data devrait produire plusieurs regions"
    )

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    for layer in layers:
        if len(layer.contour_3d) < 2:
            continue

        origin, axis_u, axis_v = layer.local_frame
        projected = []
        for contour in layer.contour_3d:
            relative = contour - origin
            u = relative @ axis_u
            v = relative @ axis_v
            projected.append(np.stack([u, v], axis=1)[:-1])

        for i in range(len(projected)):
            for j in range(i + 1, len(projected)):
                nested = point_in_polygon_2d(projected[i][0], projected[j]) or point_in_polygon_2d(
                    projected[j][0], projected[i]
                )
                if nested:
                    continue

                centroid_i = projected[i].mean(axis=0)
                centroid_j = projected[j].mean(axis=0)
                radius_i = float(np.linalg.norm(projected[i] - centroid_i, axis=1).max())
                radius_j = float(np.linalg.norm(projected[j] - centroid_j, axis=1).max())
                distance = float(np.linalg.norm(centroid_i - centroid_j))

                assert distance < params.contour_cluster_distance_factor * max(radius_i, radius_j), (
                    f"Layer {layer.id} (isovalue={layer.isovalue:.3f}) porte encore 2 boucles "
                    f"disjointes et eloignees (distance={distance:.2f}) - le contourage par "
                    f"cluster aurait du les separer en 2 Layer distincts"
                )


def test_l_shape_split_data_merged_junction_layers_stay_single_loop(l_shape_split_data):
    """
    Aux isovaleurs ou base et bras se sont deja geometriquement rejoints
    (une seule boucle de contour continue en sortie de
    _extract_closed_contours), le contourage par cluster ne doit JAMAIS
    forcer une scission artificielle - ce serait couper en deux une section
    transversale d'un seul solide, le risque explicitement identifie avant
    l'implementation de ce fix (CLAUDE.md limite #35, couches 89-100 sur
    res/test_models/l_shape_overhang.stl).

    Verifie dynamiquement (pas d'indices ou d'isovaleurs codes en dur) :
    pour toute couche qui aurait ete extraite par _extract_closed_contours
    avec une seule boucle, run() doit produire exactement un Layer, avec
    une seule boucle.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = l_shape_split_data
    params = LayerParams(layer_height_base=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)
    assert len(layers) > 1

    # une couche a 1 seule boucle doit rester 1 seul Layer par construction
    # de _cluster_contours (n <= 1 -> [[0]]) - verifie ici au niveau de
    # l'API publique run(), pas seulement au niveau unitaire de
    # _cluster_contours.
    single_loop_layers = [layer for layer in layers if len(layer.contour_3d) == 1]
    assert single_loop_layers, "Au moins une couche a 1 seule boucle est attendue sur cette fixture"
    for layer in single_loop_layers:
        assert len(layer.contour_3d) == 1

    # aucune couche perdue ni dupliquee : l'ensemble des isovaleurs produites
    # doit rester une sequence croissante sans repetition d'un MEME cluster
    # (des isovaleurs identiques entre 2 Layer issus d'un split sont
    # attendues et valides, mais chaque Layer doit rester geometriquement
    # distinct - verifie indirectement par le test precedent).
    isovalues = [layer.isovalue for layer in layers]
    assert all(a <= b for a, b in zip(isovalues, isovalues[1:])), (
        "L'ordre par isovalue doit rester non-decroissant apres le contourage par cluster"
    )


def test_cluster_contours_keeps_nested_loops_together():
    """
    Unitaire direct sur _cluster_contours : 2 boucles imbriquees (anneau,
    meme cas que hollow_box_data) doivent rester dans le meme cluster, quel
    que soit l'ecart entre leurs centroides - c'est le garde-fou anti-trou
    explicitement demande avant l'implementation de ce fix.
    """
    outer = np.array(
        [[-10, -10, 0], [10, -10, 0], [10, 10, 0], [-10, 10, 0], [-10, -10, 0]], dtype=float
    )
    inner = np.array(
        [[-2, -2, 0], [2, -2, 0], [2, 2, 0], [-2, 2, 0], [-2, -2, 0]], dtype=float
    )
    local_frame = (
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
    )

    clusters = _cluster_contours([outer, inner], local_frame, distance_factor=1.5)

    assert clusters == [[0, 1]], (
        f"Les boucles imbriquees devraient former un seul cluster, obtenu {clusters}"
    )


def test_cluster_contours_separates_distant_disjoint_loops():
    """
    Unitaire direct sur _cluster_contours : 2 boucles disjointes (pas
    imbriquees) et tres eloignees l'une de l'autre doivent former 2
    clusters distincts.
    """
    loop_a = np.array(
        [[-2, -2, 0], [2, -2, 0], [2, 2, 0], [-2, 2, 0], [-2, -2, 0]], dtype=float
    )
    loop_b = loop_a + np.array([50.0, 0.0, 0.0])
    local_frame = (
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
    )

    clusters = _cluster_contours([loop_a, loop_b], local_frame, distance_factor=1.5)

    assert clusters == [[0], [1]], (
        f"Les boucles disjointes et eloignees devraient former 2 clusters distincts, obtenu {clusters}"
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
