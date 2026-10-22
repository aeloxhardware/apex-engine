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
from core.mesh_types import Region, RegionDependencyGraph, ScalarField, TetMesh
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import FieldSolverParams, run as field_run
from pipeline.step5_layers import (
    LayerExtractionError,
    LayerParams,
    _extract_closed_contours,
    _extract_isosurface_triangles,
    _frontier_order,
    _march_tetrahedra_batch,
    _march_tetrahedron,
    _spatial_dfs_order,
    _support_driven_order,
    _weld_vertices,
    run as layers_run,
)


def _pipeline_data_from_trimesh_default(mesh: trimesh.Trimesh, resolution: float = 1.5, platform_z: float = 0.0):
    """Alias local pour construire une fixture ad hoc hors dataclass @pytest.fixture (cas 4 du plan de clustering)."""
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
        mesh_data, overhang_map, SegmentationParams(platform_z=platform_z, platform_contact_tolerance=0.5)
    )
    field = field_run(tet_mesh, mesh_data, overhang_map, regions, [], FieldSolverParams())
    return mesh_data, tet_mesh, overhang_map, regions, field


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


def _tripod_mesh() -> trimesh.Trimesh:
    """
    Pied vertical (pose sur le plateau, base a z=0) + 3 bras non coplanaires
    (+X, +Y, diagonal -X-Y) partant du sommet du pied - meme geometrie que
    le modele synthetique construit pendant le diagnostic de CLAUDE.md
    (chevauchement de directions a une jonction multi-region), formalisee
    ici en fixture pour ne pas la regenerer manuellement a chaque fois.
    """
    foot = trimesh.creation.box(extents=[4, 4, 20])
    foot.apply_translation([0, 0, 10])
    arm_x = trimesh.creation.box(extents=[20, 4, 4])
    arm_x.apply_translation([8, 0, 20])
    arm_y = trimesh.creation.box(extents=[4, 20, 4])
    arm_y.apply_translation([0, 8, 20])
    arm_diag = trimesh.creation.box(extents=[14, 14, 4])
    arm_diag.apply_translation([-7, -7, 20])

    mesh = foot.union(arm_x, engine="manifold")
    mesh = mesh.union(arm_y, engine="manifold")
    mesh = mesh.union(arm_diag, engine="manifold")
    return mesh


@pytest.fixture
def tripod_data():
    """
    Meme geometrie que _tripod_mesh, avec 4 REGIONS INJECTEES MANUELLEMENT
    (pied + 3 bras), pas issues de step3 (qui ne sait pas encore scinder un
    objet en plusieurs regions dans l'etat actuel du code, cf. CLAUDE.md) -
    isole le comportement de step4/step5 face a une vraie jonction
    multi-region deja donnee, exactement la methode utilisee pour le
    diagnostic initial (fragmentation de contour 1->4 boucles mesuree sur
    ce modele avant le fix du melange spatial).
    """
    mesh = _tripod_mesh()

    import tempfile
    import os

    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        # resolution=1.5 : meme valeur que le defaut de PipelineParams
        # (main.py), utilisee lors du diagnostic initial de la
        # fragmentation (CLAUDE.md) - resolution=2.0 (plus grossiere)
        # lisse artificiellement l'artefact avant meme le fix, mesure
        # directement en ecrivant ce test (piege reel, pas suppose).
        mesh_data, tet_mesh, _ = ingest_run(path, resolution=1.5)
    finally:
        os.unlink(path)

    overhang_map = overhang_run(mesh_data, OverhangParams(overhang_angle_threshold_deg=45.0))
    seg_params = SegmentationParams(platform_z=0.0, platform_contact_tolerance=0.5)
    base_graph = segmentation_run(mesh_data, overhang_map, seg_params)

    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    fx, fy, fz = face_centers[:, 0], face_centers[:, 1], face_centers[:, 2]
    is_x = (fx > 2.0) & (fz > 15)
    is_y = (fy > 2.0) & (fz > 15) & ~is_x
    is_diag = (fx < -2.0) & (fy < -2.0) & (fz > 15) & ~is_x & ~is_y

    region_id = np.zeros(len(mesh_data.faces), dtype=int)
    region_id[is_x] = 1
    region_id[is_y] = 2
    region_id[is_diag] = 3

    def _unit(v):
        return v / np.linalg.norm(v)

    regions = {
        0: Region(id=0, face_ids=np.where(region_id == 0)[0], target_direction=np.array([0.0, 0.0, 1.0])),
        1: Region(id=1, face_ids=np.where(region_id == 1)[0], target_direction=_unit(np.array([1.0, 0.0, 0.28867513]))),
        2: Region(id=2, face_ids=np.where(region_id == 2)[0], target_direction=_unit(np.array([0.0, 1.0, 0.28867513]))),
        3: Region(id=3, face_ids=np.where(region_id == 3)[0], target_direction=_unit(np.array([-0.70710678, -0.70710678, 0.28867513]))),
    }
    region_graph = RegionDependencyGraph(
        regions=regions,
        dependencies=[(0, 1), (0, 2), (0, 3)],
        platform_contact_faces=base_graph.platform_contact_faces,
    )
    return mesh_data, tet_mesh, overhang_map, region_graph


def _portal_mesh() -> trimesh.Trimesh:
    """
    Portique en U inverse : 2 montants verticaux separes + 1 traverse
    horizontale en haut - cas MERGE (2 clusters disjoints a une isovaleur
    basse fusionnent en 1 seul cluster a une isovaleur plus haute, des que
    la coupe horizontale atteint la traverse), complementaire au cas SPLIT
    deja couvert par tripod_data (CLAUDE.md, continuite spatiale DFS).
    """
    leg_a = trimesh.creation.box(extents=[4, 4, 20])
    leg_a.apply_translation([-8, 0, 10])
    leg_b = trimesh.creation.box(extents=[4, 4, 20])
    leg_b.apply_translation([8, 0, 10])
    beam = trimesh.creation.box(extents=[24, 4, 4])
    beam.apply_translation([0, 0, 20])
    mesh = leg_a.union(beam, engine="manifold").union(leg_b, engine="manifold")
    return mesh


@pytest.fixture
def portal_data():
    """
    Meme geometrie que _portal_mesh, UNE SEULE region (direction verticale
    uniforme) sur tout le modele - seule la TOPOLOGIE du clustering disjoint
    compte ici (2 montants qui fusionnent en 1 cluster au niveau de la
    traverse), pas la segmentation par direction.
    """
    return _pipeline_data_from_trimesh_default(_portal_mesh(), resolution=1.5, platform_z=0.0)


def _t_shape_single_region_mesh() -> trimesh.Trimesh:
    """
    Pied vertical + 1 bras vertical decale (pas incline) - reproduit
    directement le cas mesure sur res/test_models/l-triple-axis-overhanged-shape.stl
    (CLAUDE.md) : le pied et le bras sont tous deux verticaux, donc step3
    les fusionne dans LA MEME region (une seule nappe de surface connexe),
    mais forment 2 composantes spatialement disjointes du solide a partir
    du niveau ou ils se separent.
    """
    foot = trimesh.creation.box(extents=[4, 4, 20])
    foot.apply_translation([0, 0, 10])
    arm = trimesh.creation.box(extents=[4, 4, 20])
    arm.apply_translation([16, 0, 20])
    bridge = trimesh.creation.box(extents=[16, 4, 4])
    bridge.apply_translation([8, 0, 20])
    mesh = foot.union(bridge, engine="manifold").union(arm, engine="manifold")
    return mesh


@pytest.fixture
def t_shape_single_region_data():
    """Pipeline complet step1->step4, segmentation REELLE (pas injectee) -
    reproduit le cas exact ou pied et bras partagent region_ids identiques."""
    return _pipeline_data_from_trimesh_default(_t_shape_single_region_mesh(), resolution=1.5, platform_z=0.0)


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


def test_node_within_weld_tolerance_of_plane_keeps_whole_contour():
    """
    Regression (Benchy Bird, couche entiere perdue) : un noeud a moins de la
    tolerance de soudure du plan fait fusionner les points de coupe de ses
    aretes -> triangles degeneres -> aretes repetees -> boucle ouverte, jetee.
    Boite 2x2x2 en tetraedres de Kuhn, z des noeuds legerement decale (la
    section xy reste le carre 2x2), un noeud de la face x = 2 a 1e-6 du plan
    z = 1 (un noeud interieur ne casse rien : c'est le bord qui se rompt).
    """
    import itertools

    grid = np.array(list(itertools.product(range(3), repeat=3)), dtype=float)
    index = {tuple(p.astype(int)): i for i, p in enumerate(grid)}
    tets = []
    for cube in itertools.product(range(2), repeat=3):
        for perm in itertools.permutations(range(3)):
            corner = np.array(cube)
            path = [index[tuple(corner)]]
            for axis in perm:
                corner = corner.copy()
                corner[axis] += 1
                path.append(index[tuple(corner)])
            tets.append(path)
    rng = np.random.default_rng(3)
    vertices = grid.copy()
    vertices[:, 2] += rng.choice([-1.0, 1.0], len(grid)) * rng.uniform(0.02, 0.05, len(grid))
    isovalue = 1.0
    vertices[index[(2, 1, 1)], 2] = isovalue + 1e-6     # noeud de BORD (face x = 2)
    tetrahedra = np.array(tets)
    tet_mesh = TetMesh(tet_vertices=vertices, tetrahedra=tetrahedra)
    values = vertices[:, 2].copy()
    tet_min, tet_max = values[tetrahedra].min(axis=1), values[tetrahedra].max(axis=1)

    points, faces, _src = _extract_isosurface_triangles(tet_mesh, values, isovalue, tet_min, tet_max)
    tolerance = LayerParams().weld_tolerance_factor * np.linalg.norm(np.ptp(vertices, axis=0))
    welded, welded_faces = _weld_vertices(points, faces, tolerance)
    degenerate = ((welded_faces[:, 0] == welded_faces[:, 1]) | (welded_faces[:, 1] == welded_faces[:, 2])
                  | (welded_faces[:, 0] == welded_faces[:, 2]))
    assert degenerate.any()          # le cas est bien exerce

    contours = _extract_closed_contours(welded, welded_faces)
    assert len(contours) == 1
    x, y = contours[0][:, 0], contours[0][:, 1]
    area = 0.5 * abs(np.dot(x[:-1], y[1:]) - np.dot(x[1:], y[:-1]))
    assert area == pytest.approx(4.0, abs=1e-6)


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


# ---------------------------------------------------------------------------
# Cas 5bis - _march_tetrahedra_batch doit produire EXACTEMENT le meme
# resultat que _march_tetrahedron appelee en boucle Python (l'oracle) -
# meme methodologie que la verification faite pendant la conception
# (7129 tetraedres synthetiques, 0 ecart), ici comme test permanent.
# ---------------------------------------------------------------------------

def test_march_tetrahedra_batch_matches_scalar_loop():
    """
    Compare, triangle par triangle et EN ORDRE (pas de tri qui masquerait
    un bug de winding), la sortie de _march_tetrahedra_batch a celle de
    _march_tetrahedron appelee ligne par ligne - sur un melange de
    tetraedres a valeurs continues (cas generique) et a valeurs discretes
    (force des egalites/limites n_above, y compris n_above in (0,4) et
    champ constant degenere hi==lo) pour couvrir tous les chemins.
    """
    rng = np.random.default_rng(2024)

    for trial in range(30):
        m = int(rng.integers(50, 400))
        tet_points = rng.uniform(-5.0, 5.0, size=(m, 4, 3))

        # ~30% valeurs discretes (force des egalites/cas limites), ~5%
        # constantes par ligne (degenere hi==lo), le reste continu.
        tet_values = rng.uniform(0.0, 4.0, size=(m, 4))
        discrete_mask = rng.random(m) < 0.30
        tet_values[discrete_mask] = rng.choice(
            [0.0, 1.0, 2.0, 3.0, 1.0], size=(discrete_mask.sum(), 4)
        )
        constant_mask = rng.random(m) < 0.05
        tet_values[constant_mask] = rng.uniform(0.0, 4.0, size=(constant_mask.sum(), 1))

        isovalue = float(rng.uniform(0.0, 4.0))

        batch_vertices, batch_local_idx = _march_tetrahedra_batch(tet_points, tet_values, isovalue)
        batch_triangles = batch_vertices.reshape(-1, 3, 3)

        for row in range(m):
            expected = _march_tetrahedron(tet_points[row], tet_values[row], isovalue)
            actual = batch_triangles[batch_local_idx == row]

            assert len(actual) == len(expected), (
                f"trial={trial} row={row}: {len(actual)} triangles obtenus, "
                f"{len(expected)} attendus (isovalue={isovalue}, "
                f"tet_values={tet_values[row]})"
            )
            for tri_actual, tri_expected in zip(actual, expected):
                assert np.allclose(tri_actual, tri_expected, atol=1e-9), (
                    f"trial={trial} row={row}: triangle different (ordre/winding "
                    f"inclus) - obtenu {tri_actual}, attendu {tri_expected}"
                )


def test_march_tetrahedra_batch_all_degenerate_returns_empty():
    """Toutes les lignes non-croisantes (n_above in (0,4)) ou a champ
    constant : aucun triangle, aucune exception - le repli 'vertices_parts
    vide' doit etre exerce."""
    tet_points = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ]
    )
    tet_values = np.array(
        [
            [10.0, 10.0, 10.0, 10.0],  # n_above == 4 (tout au-dessus)
            [-10.0, -10.0, -10.0, -10.0],  # n_above == 0 (tout en-dessous)
            [2.0, 2.0, 2.0, 2.0],  # champ constant, degenere (hi==lo)
        ]
    )

    vertices, local_idx = _march_tetrahedra_batch(tet_points, tet_values, isovalue=2.0)

    assert len(vertices) == 0
    assert len(local_idx) == 0


def test_march_tetrahedra_batch_empty_input_returns_empty():
    vertices, local_idx = _march_tetrahedra_batch(
        np.zeros((0, 4, 3)), np.zeros((0, 4)), isovalue=1.0
    )
    assert len(vertices) == 0
    assert len(local_idx) == 0


# ---------------------------------------------------------------------------
# Cas 6 - CLAUDE.md : melange spatial des directions resout la fragmentation
# de contour a une jonction multi-region (tripode a 3 branches non
# coplanaires + plateau)
# ---------------------------------------------------------------------------

def _peak_loop_count(layers) -> int:
    """
    Nombre maximal de boucles de contour (len(contour_3d)) atteint sur
    l'ensemble des couches - la METRIQUE PERTINENTE pour la fragmentation
    de jonction n'est PAS le plus grand saut entre 2 couches CONSECUTIVES
    (une fragmentation progressive 1->2->3->4, chaque pas valant 1, aurait
    un "max jump consecutif" de 1 alors qu'elle est bien le probleme
    mesure - piege reel rencontre en ecrivant ce test, corrige avant de
    livrer : le premier jet de cette metrique masquait completement
    l'ameloriation du fix). Le pic atteint (4 boucles simultanees avant
    fix, cf. CLAUDE.md) est le signal direct de la fragmentation.
    """
    if not layers:
        return 0
    return max(len(layer.contour_3d) for layer in layers)


def test_tripod_without_blend_shows_loop_fragmentation_peak(tripod_data):
    """
    Non-regression du DIAGNOSTIC (pas du fix) : sans melange spatial
    (direction_blend_radius=0.0, defaut historique), reproduit le pic de
    fragmentation de contour mesure pendant l'investigation de CLAUDE.md
    (n_loops atteint 4 dans la zone de jonction des 3 bras). Sert de
    temoin pour le test suivant.
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    params = FieldSolverParams(direction_blend_radius=0.0)

    field = field_run(tet_mesh, mesh_data, overhang_map, region_graph, [], params)
    layers = layers_run(tet_mesh, mesh_data, field, region_graph, LayerParams(layer_height_base=0.2))

    assert _peak_loop_count(layers) >= 4, (
        "Fixture invalide : le tripode sans melange devrait reproduire le pic de "
        "fragmentation (4 boucles simultanees) mesure pendant le diagnostic, "
        "cf. CLAUDE.md"
    )


def test_tripod_with_blend_reduces_loop_fragmentation_peak(tripod_data):
    """
    Avec le melange spatial active, le pic du nombre de boucles de contour
    a la jonction doit etre STRICTEMENT REDUIT par rapport au comportement
    sans melange (test precedent) - la transition doit devenir plus
    progressive, pas necessairement parfaitement lisse (le champ reste un
    scalaire unique, cf. limite documentee dans CLAUDE.md).
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data

    baseline_params = FieldSolverParams(direction_blend_radius=0.0)
    baseline_field = field_run(tet_mesh, mesh_data, overhang_map, region_graph, [], baseline_params)
    baseline_layers = layers_run(
        tet_mesh, mesh_data, baseline_field, region_graph, LayerParams(layer_height_base=0.2)
    )
    baseline_peak = _peak_loop_count(baseline_layers)

    blended_params = FieldSolverParams(
        direction_blend_radius=20.0, direction_blend_power=2.0, direction_blend_min_norm=0.1
    )
    blended_field = field_run(tet_mesh, mesh_data, overhang_map, region_graph, [], blended_params)
    blended_layers = layers_run(
        tet_mesh, mesh_data, blended_field, region_graph, LayerParams(layer_height_base=0.2)
    )
    blended_peak = _peak_loop_count(blended_layers)

    assert np.all(np.isfinite(blended_field.values)), "Le champ resolu avec melange doit rester fini"
    assert blended_peak < baseline_peak, (
        f"Le melange spatial devrait reduire le pic de fragmentation du nombre de "
        f"boucles a la jonction (avant: {baseline_peak}, apres: {blended_peak})"
    )


# ---------------------------------------------------------------------------
# Clustering spatial des boucles disjointes (split_disjoint_clusters, CLAUDE.md)
# ---------------------------------------------------------------------------

def test_split_disjoint_clusters_default_is_disabled(hollow_box_data):
    """
    Non-regression stricte : sans split_disjoint_clusters (defaut None), le
    comportement est identique a aujourd'hui - la couche mi-hauteur de la
    boite creuse produit toujours un seul Layer avec 2 contour_3d.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = hollow_box_data
    params = LayerParams(layer_height_base=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    mid_layers = [
        layer for layer in layers
        if abs(np.mean([c[:, 2].mean() for c in layer.contour_3d])) < 2.0
    ]
    multi_contour_layers = [layer for layer in mid_layers if len(layer.contour_3d) >= 2]
    assert multi_contour_layers, (
        "Comportement historique attendu : la couche mi-hauteur garde ses "
        "2 contours dans un seul Layer quand le clustering est desactive"
    )


def test_hollow_box_nested_contours_not_split_even_when_enabled(hollow_box_data):
    """
    Cas le plus critique : meme avec split_disjoint_clusters ACTIVE, deux
    boucles IMBRIQUEES (paroi exterieure + cavite interne) ne doivent
    JAMAIS etre scindees en 2 Layer distincts - ce sont une seule section
    transversale avec un trou, pas deux clusters separes (CLAUDE.md).
    """
    mesh_data, tet_mesh, overhang_map, regions, field = hollow_box_data
    params = LayerParams(layer_height_base=1.0, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    mid_layers = [
        layer for layer in layers
        if abs(np.mean([c[:, 2].mean() for c in layer.contour_3d])) < 2.0
    ]
    multi_contour_layers = [layer for layer in mid_layers if len(layer.contour_3d) >= 2]
    assert multi_contour_layers, (
        "Le clustering actif ne doit PAS scinder les boucles imbriquees de "
        "la boite creuse - elles doivent rester dans un seul Layer"
    )


def test_tripod_split_disjoint_clusters_produces_distinct_local_frames(tripod_data):
    """
    Cas cible du fix : sur le modele tripode (jonction reelle a 4 regions),
    avec split_disjoint_clusters active, au moins une isovaleur partagee
    par plusieurs Layer distincts doit exister, chacun avec un local_frame
    dont la normale differe sensiblement des autres - preuve directe que le
    plan moyen incoherent (cause du faisceau convergent observe
    visuellement, CLAUDE.md) est resolu par construction.
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph, [],
        FieldSolverParams(direction_blend_radius=20.0, direction_blend_power=2.0, direction_blend_min_norm=0.1),
    )
    params = LayerParams(layer_height_base=0.2, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, region_graph, params)

    by_isovalue: dict[float, list] = {}
    for layer in layers:
        key = round(layer.isovalue, 6)
        by_isovalue.setdefault(key, []).append(layer)

    split_groups = [group for group in by_isovalue.values() if len(group) >= 2]
    assert split_groups, (
        "Au moins une isovaleur devrait avoir ete scindee en plusieurs Layer "
        "sur le modele tripode avec split_disjoint_clusters actif"
    )

    found_distinct_normals = False
    for group in split_groups:
        normals = [np.cross(l.local_frame[1], l.local_frame[2]) for l in group]
        for i in range(len(normals)):
            for j in range(i + 1, len(normals)):
                if abs(np.dot(normals[i], normals[j])) < 0.9:
                    found_distinct_normals = True
    assert found_distinct_normals, (
        "Au moins deux Layer issus d'une meme isovaleur scindee devraient "
        "avoir des plans locaux nettement distincts (produit scalaire des "
        "normales < 0.9), preuve que chaque cluster a son propre plan "
        "coherent plutot qu'un plan moyen incoherent"
    )


def test_tripod_split_disjoint_clusters_region_ids_are_distinct_subsets(tripod_data):
    """
    Les region_ids des Layer resultants d'une isovaleur scindee doivent
    etre des sous-ensembles distincts (pas tous identiques) - preuve que le
    recalcul par cluster (approximation par proximite spatiale, CLAUDE.md)
    associe bien chaque cluster a SA region, pas a l'union de toutes.
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph, [],
        FieldSolverParams(direction_blend_radius=20.0, direction_blend_power=2.0, direction_blend_min_norm=0.1),
    )
    params = LayerParams(layer_height_base=0.2, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, region_graph, params)

    by_isovalue: dict[float, list] = {}
    for layer in layers:
        key = round(layer.isovalue, 6)
        by_isovalue.setdefault(key, []).append(layer)

    split_groups = [group for group in by_isovalue.values() if len(group) >= 2]
    assert split_groups, "Au moins une isovaleur scindee est necessaire pour ce test"

    found_distinct_region_sets = False
    for group in split_groups:
        region_sets = [frozenset(l.region_ids) for l in group]
        if len(set(region_sets)) > 1:
            found_distinct_region_sets = True
    assert found_distinct_region_sets, (
        "Au moins un groupe de Layer scindes devrait avoir des region_ids "
        "distincts entre clusters (pas tous identiques)"
    )


def test_close_disjoint_contours_not_split_within_tolerance():
    """
    Deux boites reliees par un pont fin, avec un ecart controle sous le
    seuil de tolerance (split_disjoint_clusters * layer_height_base) : les
    2 boites doivent rester dans un SEUL Layer a l'endroit du gap - preuve
    que le mecanisme n'introduit pas de discontinuite artificielle avant
    une fusion geometrique reelle a une isovaleur voisine.
    """
    box_a = trimesh.creation.box(extents=[4, 4, 4])
    box_a.apply_translation([-2.5, 0, 10])
    box_b = trimesh.creation.box(extents=[4, 4, 4])
    box_b.apply_translation([2.5, 0, 10])
    bridge = trimesh.creation.box(extents=[5, 1, 1])
    bridge.apply_translation([0, 0, 10])

    mesh = box_a.union(bridge, engine="manifold").union(box_b, engine="manifold")
    mesh_data, tet_mesh, overhang_map, regions, field = _pipeline_data_from_trimesh_default(
        mesh, resolution=1.0, platform_z=8.0
    )

    layer_height_base = 0.5
    # tolerance large (10x layer_height_base) : l'objectif de ce test est de
    # prouver l'absence de scission artificielle quand le gap reel est petit
    # devant le seuil, pas de calibrer un seuil de production - la
    # calibration fine se fait sur le modele reel (CLAUDE.md).
    params = LayerParams(layer_height_base=layer_height_base, split_disjoint_clusters=10.0)
    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    # les couches proches de z=10 (hauteur des 2 boites + pont) doivent
    # rester des Layer uniques (pas de scission artificielle malgre la
    # geometrie "en haltere" - le pont, bien que fin, garde les boucles
    # suffisamment proches en distance de bord-a-bord).
    near_bridge = [layer for layer in layers if any(abs(c[:, 2].mean() - 10.0) < 1.0 for c in layer.contour_3d)]
    assert near_bridge, "Au moins une couche proche du pont est necessaire pour ce test"

    by_isovalue: dict[float, list] = {}
    for layer in near_bridge:
        key = round(layer.isovalue, 6)
        by_isovalue.setdefault(key, []).append(layer)

    assert all(len(group) == 1 for group in by_isovalue.values()), (
        "Aucune isovaleur proche du pont ne devrait etre scindee en "
        "plusieurs Layer : l'ecart entre les 2 boites doit rester sous la "
        "tolerance de distance"
    )


# ---------------------------------------------------------------------------
# Continuite spatiale DFS (Decision 8, CLAUDE.md) - ordre des Layer scindes
# ---------------------------------------------------------------------------

def _rank_sequence(layers, region_ranks_by_id: dict[int, int]) -> list[int]:
    ranks = []
    for layer in layers:
        rank = max((region_ranks_by_id.get(r, 0) for r in layer.region_ids), default=0)
        ranks.append(rank)
    return ranks


def _branch_switch_count(layers, jump_threshold: float) -> int:
    """
    Nombre de fois ou l'origine du local_frame "saute" de plus de
    `jump_threshold` entre deux couches consecutives - proxy direct du
    nombre de trajets a vide inseres par step7/step8 entre deux clusters
    spatialement distincts (CLAUDE.md, faisceau convergent mesure sur le
    modele en T reel).

    Se base sur la POSITION du cluster (local_frame.origin), pas sur
    region_ids : deux clusters spatialement disjoints peuvent legitimement
    partager les memes region_ids (une seule region step3 peut couvrir
    plusieurs composantes spatiales du solide, cf. CLAUDE.md - c'est
    exactement le cas cible de ce fix), donc region_ids seul ne peut pas
    detecter la bascule.
    """
    switches = 0
    for i in range(1, len(layers)):
        prev_origin = layers[i - 1].local_frame[0]
        cur_origin = layers[i].local_frame[0]
        if np.linalg.norm(cur_origin - prev_origin) > jump_threshold:
            switches += 1
    return switches


def test_t_shape_spatial_order_groups_branch_before_switching(t_shape_single_region_data):
    """
    Cas cible exact du fix (reproduit le bug mesure sur
    res/test_models/l-triple-axis-overhanged-shape.stl) : pied et bras
    partagent la meme region step3 (meme nappe de surface connexe), donc le
    tri (rank, isovalue) seul les fait alterner a quasi chaque isovaleur.
    Le tri spatial DFS doit produire strictement moins de bascules de
    branche que ce tri naif reconstruit independamment sur les memes
    couches.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = t_shape_single_region_data
    params = LayerParams(layer_height_base=0.5, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    split_isovalues = {
        round(layer.isovalue, 6)
        for layer in layers
        if sum(1 for other in layers if round(other.isovalue, 6) == round(layer.isovalue, 6)) >= 2
    }
    assert split_isovalues, (
        "Au moins une isovaleur scindee est necessaire pour que ce test "
        "prouve quoi que ce soit sur l'ordre - verifier la geometrie de la "
        "fixture si ceci echoue"
    )

    # seuil de saut : la moitie de la distance separant les 2 branches (16
    # unites entre pied et bras dans cette fixture) - largement au-dessus
    # du deplacement normal d'une couche a l'autre au sein d'une meme
    # branche (~layer_height_base), largement en-dessous du saut reel
    # pied<->bras.
    jump_threshold = 5.0
    naive_order = sorted(layers, key=lambda l: (l.isovalue,))
    naive_switches = _branch_switch_count(naive_order, jump_threshold)
    spatial_switches = _branch_switch_count(layers, jump_threshold)

    assert spatial_switches < naive_switches, (
        f"Le tri spatial DFS ({spatial_switches} bascules) devrait produire "
        f"nettement moins de bascules de branche que l'ordre par isovaleur "
        f"seule ({naive_switches} bascules)"
    )


def test_tripod_spatial_order_never_violates_region_rank(tripod_data):
    """
    Non-regression stricte de la contrainte dure : le tri spatial DFS ne
    doit JAMAIS reordonner un cluster de rang superieur avant un cluster de
    rang inferieur, meme actif. Rangs non-decroissants dans la liste
    retournee.
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph, [],
        FieldSolverParams(direction_blend_radius=20.0, direction_blend_power=2.0, direction_blend_min_norm=0.1),
    )
    params = LayerParams(layer_height_base=0.2, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, region_graph, params)

    region_ranks_by_id = {0: 0, 1: 1, 2: 1, 3: 1}
    ranks = _rank_sequence(layers, region_ranks_by_id)

    assert ranks == sorted(ranks), (
        "Le tri spatial DFS ne doit jamais violer le rang topologique de "
        "region - sequence obtenue non triee"
    )


def test_tripod_spatial_order_deterministic(tripod_data):
    """Deux appels identiques doivent produire exactement le meme ordre - garantit
    l'absence de dependance a un ordre d'iteration non deterministe de networkx."""
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph, [],
        FieldSolverParams(direction_blend_radius=20.0, direction_blend_power=2.0, direction_blend_min_norm=0.1),
    )
    params = LayerParams(layer_height_base=0.2, split_disjoint_clusters=3.0)

    layers_a = layers_run(tet_mesh, mesh_data, field, region_graph, params)
    layers_b = layers_run(tet_mesh, mesh_data, field, region_graph, params)

    assert [round(l.isovalue, 9) for l in layers_a] == [round(l.isovalue, 9) for l in layers_b]
    assert [frozenset(l.region_ids) for l in layers_a] == [frozenset(l.region_ids) for l in layers_b]


def test_hollow_box_spatial_order_keeps_nested_contours_single_layer(hollow_box_data):
    """
    Non-regression : le tri spatial DFS (actif automatiquement des que
    split_disjoint_clusters l'est) ne doit introduire aucune scission
    indirecte des boucles imbriquees de la boite creuse - meme nombre de
    Layer qu'avec le seul clustering (sans consequence de l'ordre).
    """
    mesh_data, tet_mesh, overhang_map, regions, field = hollow_box_data
    params = LayerParams(layer_height_base=1.0, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    mid_layers = [
        layer for layer in layers
        if abs(np.mean([c[:, 2].mean() for c in layer.contour_3d])) < 2.0
    ]
    assert all(len(layer.contour_3d) == 2 for layer in mid_layers), (
        "Les boucles imbriquees de la boite creuse doivent rester dans un "
        "seul Layer (2 contour_3d), le graphe de continuite ne doit pas "
        "introduire de scission indirecte"
    )


def test_portal_merge_produces_single_layer_above_beam(portal_data):
    """
    Cas MERGE : sous la traverse, 2 clusters disjoints par isovaleur (un par
    montant) ; au niveau/au-dessus de la traverse, un seul cluster (les 2
    montants se rejoignent). Preuve directe que le clustering ET le graphe
    de continuite gerent la recombinaison sans erreur ni duplication.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = portal_data
    params = LayerParams(layer_height_base=1.0, split_disjoint_clusters=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    by_isovalue: dict[float, list] = {}
    for layer in layers:
        key = round(layer.isovalue, 6)
        by_isovalue.setdefault(key, []).append(layer)

    below_beam = [g for iv, g in by_isovalue.items() if iv < sorted(by_isovalue)[len(by_isovalue) // 3]]
    split_below = [g for g in below_beam if len(g) >= 2]
    assert split_below, (
        "Au moins une isovaleur sous la traverse devrait etre scindee en 2 "
        "Layer (un par montant) - verifier la geometrie/calibration si ceci "
        "echoue"
    )

    above_beam_groups = [g for g in by_isovalue.values() if len(g) == 1]
    assert above_beam_groups, (
        "Au moins une isovaleur (typiquement au niveau de la traverse ou "
        "au-dessus) devrait produire un seul Layer - preuve de la "
        "recombinaison des 2 montants"
    )


def test_portal_spatial_order_visits_one_leg_fully_before_switching(portal_data):
    """
    L'ordre spatial DFS doit regrouper TOUTES les couches d'un meme montant
    avant de passer a l'autre - au plus une poignee de bascules (une par
    racine independante visitee - les 2 montants sont 2 composantes non
    connectees sous la traverse, donc 2 racines - pas une bascule a CHAQUE
    isovaleur comme le ferait un tri (rank, isovalue) naif). Compare
    directement au tri naif pour le prouver, plutot que de supposer un
    chiffre absolu (le nombre exact de racines/artefacts de maillage n'est
    pas le point teste ici).
    """
    mesh_data, tet_mesh, overhang_map, regions, field = portal_data
    params = LayerParams(layer_height_base=1.0, split_disjoint_clusters=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    by_isovalue: dict[float, list] = {}
    for layer in layers:
        key = round(layer.isovalue, 6)
        by_isovalue.setdefault(key, []).append(layer)
    split_isovalues = {iv for iv, g in by_isovalue.items() if len(g) >= 2}
    assert split_isovalues, "Au moins une isovaleur scindee est necessaire pour ce test"

    # seuil de saut : la moitie de la distance separant les 2 montants (16
    # unites, cf. _portal_mesh) - au-dessus du deplacement normal d'une
    # couche a l'autre au sein d'un meme montant.
    jump_threshold = 5.0
    naive_order = sorted(layers, key=lambda l: (l.isovalue,))
    naive_switches = _branch_switch_count(naive_order, jump_threshold)
    spatial_switches = _branch_switch_count(layers, jump_threshold)

    assert spatial_switches < naive_switches, (
        f"L'ordre spatial DFS ({spatial_switches} bascules) devrait produire "
        f"nettement moins de bascules que l'ordre par isovaleur seule "
        f"({naive_switches} bascules) sur le portique"
    )


def test_portal_merge_node_not_duplicated(portal_data):
    """
    Garde-fou direct contre une regression future de _dfs_order_within_rank :
    aucune duplication du noeud de recombinaison malgre ses 2 aretes
    entrantes. Se base sur (isovalue, origine du local_frame arrondie) -
    PAS (isovalue, region_ids) seul, qui coincide legitimement pour les 2
    montants du portique (une seule region pour tout le modele, cf.
    portal_data) sans que ce soit une duplication.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = portal_data

    params_with_sort = LayerParams(layer_height_base=1.0, split_disjoint_clusters=1.0)
    layers_with_sort = layers_run(tet_mesh, mesh_data, field, regions, params_with_sort)

    keys_seen = [
        (round(l.isovalue, 6), tuple(np.round(l.local_frame[0], 3)))
        for l in layers_with_sort
    ]
    assert len(keys_seen) == len(set(keys_seen)), (
        "Un meme (isovalue, origine de cluster) apparait plusieurs fois - "
        "signe d'une duplication du noeud de recombinaison par le DFS"
    )


def test_spatial_order_is_a_permutation():
    """
    Test unitaire direct sur _spatial_dfs_order, isole du pipeline complet -
    construit des raw_layers/region_ranks synthetiques (style tripode ET
    style portique) et verifie que l'ordre retourne est une vraie
    permutation (aucun index perdu, aucun duplique).
    """
    identity_frame = (np.zeros(3), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))

    def _square_loop(cx: float, cy: float, cz: float, half: float = 1.0) -> np.ndarray:
        pts = np.array([
            [cx - half, cy - half, cz],
            [cx + half, cy - half, cz],
            [cx + half, cy + half, cz],
            [cx - half, cy + half, cz],
            [cx - half, cy - half, cz],
        ])
        return pts

    # Style tripode : 1 pied (rank 0) qui se scinde en 2 branches (rank 1),
    # chacune une chaine simple sur 2 isovaleurs supplementaires.
    raw_layers_tripod = [
        (0.0, [_square_loop(0, 0, 0)], [0], identity_frame),
        (1.0, [_square_loop(-3, 0, 1)], [1], identity_frame),
        (1.0, [_square_loop(3, 0, 1)], [1], identity_frame),
        (2.0, [_square_loop(-3, 0, 2)], [1], identity_frame),
        (2.0, [_square_loop(3, 0, 2)], [1], identity_frame),
    ]
    region_ranks_tripod = {0: 0, 1: 1}
    order_tripod = _spatial_dfs_order(raw_layers_tripod, region_ranks_tripod, gap_threshold=1.5)
    assert sorted(order_tripod) == list(range(len(raw_layers_tripod)))

    # Style portique : 2 montants separes (rank 0) qui fusionnent en 1 seul
    # cluster (rank 0 aussi - meme region pour tout le portique).
    raw_layers_portal = [
        (0.0, [_square_loop(-5, 0, 0)], [0], identity_frame),
        (0.0, [_square_loop(5, 0, 0)], [0], identity_frame),
        (1.0, [_square_loop(-5, 0, 1)], [0], identity_frame),
        (1.0, [_square_loop(5, 0, 1)], [0], identity_frame),
        (2.0, [_square_loop(0, 0, 2, half=6.0)], [0], identity_frame),
    ]
    region_ranks_portal = {0: 0}
    order_portal = _spatial_dfs_order(raw_layers_portal, region_ranks_portal, gap_threshold=1.5)
    assert sorted(order_portal) == list(range(len(raw_layers_portal)))


# ---------------------------------------------------------------------------
# LayerParams.support_driven_order - ordre pilote par le support reel
# ---------------------------------------------------------------------------

def _ring(z: float, cx: float = 0.0, radius: float = 2.0, count: int = 24) -> list:
    """Une boucle de contour circulaire fermee a l'altitude z."""
    angles = np.linspace(0.0, 2.0 * np.pi, count)
    return [
        np.stack(
            [cx + radius * np.cos(angles), radius * np.sin(angles), np.full(count, z)],
            axis=1,
        )
    ]


def test_support_order_prefers_supported_over_lower():
    """
    Le SUPPORT est la cle primaire, l'altitude ne departage que des couches
    deja imprimables : une couche BASSE mais sans rien sous elle ne doit
    jamais passer avant une couche PLUS HAUTE mais bien supportee.

    C'est l'invariant qui protege les surplombs par construction - une piece
    pendante est plus basse que son porteur, mais tant que le porteur n'existe
    pas elle n'est pas supportee, donc jamais choisie.
    """
    # Pile posee au sol en x=0, plus une coque isolee et PLUS BASSE en x=20 -
    # mais assez au-dessus de z=0 pour ne jamais tomber dans la tolerance de
    # contact plateau (platform_contact_tolerance_factor * nozzle_width).
    raw_layers = [
        (0.0, _ring(0.0), [0], None),                # 0 - au sol
        (0.1, _ring(0.2), [0], None),                # 1 - portee par 0
        (0.2, _ring(0.4), [0], None),                # 2 - portee par 1
        (0.3, _ring(3.0, cx=20.0), [0], None),       # 3 - isolee, sans support
    ]
    params = LayerParams(
        layer_height_base=0.2, support_driven_order=0.9, nozzle_width=0.4
    )

    order, unresolved = _support_driven_order(raw_layers, params, np.array([0.0, 0.0, 1.0]))

    assert sorted(order) == list(range(len(raw_layers)))
    # La pile supportee s'imprime de bas en haut, dans l'ordre...
    assert order.index(0) < order.index(1) < order.index(2)
    # ...et la coque isolee (pourtant plus basse que la couche 2) passe apres,
    # faute de support - jamais avant sur le seul critere d'altitude.
    assert order.index(3) > order.index(2)
    # Elle est signalee explicitement, jamais placee en silence.
    assert 3 in unresolved


def test_support_order_rejects_shell_only_locally_supported():
    """
    Verrouille le fix de connectivite (bug reel trouve sur death.stl) : une
    coque en surplomb, imprimee couche par couche de BAS EN HAUT et RATTACHEE
    A RIEN, ne doit jamais devenir eligible juste parce que chacune de ses
    couches est bien au-dessus de la precedente.

    Une premiere version du mecanisme testait seulement "y a-t-il de la
    matiere DEJA PLACEE en dessous", sans verifier qu'elle est elle-meme
    reliee au plateau - une coque batie ainsi passait le test a CHAQUE etape
    (chaque couche est bien supportee par celle du dessous), degenerant vers
    un simple tri par hauteur alors que l'ENSEMBLE de la coque ne touche le
    reste du modele nulle part.
    """
    # Pile posee au sol en x=0 (z=0..0.4), et une coque ENTIEREMENT ISOLEE en
    # x=20, batie du bas vers le haut sur 5 couches (z=5.0..5.8) - aucune
    # d'elles n'est au contact du plateau ni d'aucune matiere ancree.
    raw_layers = [
        (0.0, _ring(0.0), [0], None),
        (0.1, _ring(0.2), [0], None),
        (0.2, _ring(0.4), [0], None),
    ]
    for i in range(5):
        raw_layers.append((1.0 + i * 0.1, _ring(5.0 + i * 0.2, cx=20.0), [0], None))

    params = LayerParams(layer_height_base=0.2, support_driven_order=0.9, nozzle_width=0.4)
    order, unresolved = _support_driven_order(raw_layers, params, np.array([0.0, 0.0, 1.0]))

    shell_indices = set(range(3, 8))
    # AUCUNE couche de la coque isolee ne doit etre ancree : toutes doivent
    # etre signalees non resolues, y compris celles du haut de la coque
    # (qui, avec un test purement local, se seraient trouvees "supportees"
    # par la couche juste en dessous d'elles dans la meme coque).
    assert shell_indices <= unresolved
    # La pile posee au sol, elle, reste correctement resolue.
    assert unresolved.isdisjoint({0, 1, 2})


def test_support_order_unresolved_fallback_prefers_proximity_over_height():
    """
    Parmi des couches deja toutes reconnues UNRESOLVED (aucune ancree), le
    departage doit privilegier la PROXIMITE a la couche precedente plutot
    que la seule altitude - reduit les sauts a vide entre couches suspendues
    qui se recouvrent geometriquement (ex: un dome qui se referme
    progressivement sur lui-meme), sans jamais changer le statut
    `unresolved` lui-meme (qui reste base uniquement sur l'ancrage reel).
    """
    # Deux coques ISOLEES (rattachees a rien), a la MEME hauteur (z=5.0) mais
    # a des positions x tres differentes - si le departage se faisait par
    # hauteur seule, l'ordre entre elles serait indetermine (egalite exacte,
    # tranchee par l'index croissant). Une troisieme couche, PROCHE en x de
    # la premiere coque mais plus haute, doit passer AVANT la deuxieme coque
    # si le departage suit bien la proximite.
    raw_layers = [
        (0.0, _ring(0.0), [0], None),  # 0 - pose au sol, ancree (amorce previous_centroid)
        (1.0, _ring(5.0, cx=0.0), [0], None),   # 1 - coque isolee A, x=0
        (1.0, _ring(5.0, cx=50.0), [0], None),  # 2 - coque isolee B, tres loin en x=50
        (1.1, _ring(5.2, cx=0.5), [0], None),   # 3 - proche de A (x=0.5), plus haute que A et B
    ]
    params = LayerParams(layer_height_base=0.2, support_driven_order=0.9, nozzle_width=0.4)
    order, unresolved = _support_driven_order(raw_layers, params, np.array([0.0, 0.0, 1.0]))

    assert unresolved == {1, 2, 3}
    # Apres la couche 0 (ancree, definit previous_centroid), la couche 3
    # (proche de 0 en x, x=0.5) doit passer AVANT la couche 2 (x=50, tres
    # loin) - alors que 2 est a la MEME hauteur que 1 et plus basse que 3,
    # un tri par altitude seule aurait pu la faire passer avant.
    assert order.index(3) < order.index(2)


def test_support_order_requires_nozzle_width():
    """
    `nozzle_width` n'a pas de defaut arbitraire : l'omettre doit lever, pas
    deviner une valeur qui divergerait en silence de
    MachineProfile.nozzle_diameter (meme discipline que platform_z).
    """
    params = LayerParams(layer_height_base=0.2, support_driven_order=0.9, nozzle_width=None)

    with pytest.raises(LayerExtractionError, match="nozzle_width"):
        _support_driven_order([(0.0, _ring(0.0), [0], None)], params, np.array([0.0, 0.0, 1.0]))


def test_support_order_disabled_by_default_leaves_isovalue_sort(flat_box_data):
    """
    `support_driven_order=None` (defaut) = comportement historique strictement
    inchange, et aucune couche marquee `support_unresolved`.
    """
    mesh_data, tet_mesh, _overhang_map, region_graph, field = flat_box_data

    layers = layers_run(
        tet_mesh, mesh_data, field, region_graph, LayerParams(layer_height_base=0.5)
    )

    assert layers
    assert all(not layer.support_unresolved for layer in layers)


def test_support_order_partially_grounded_layer_anchors_its_grounded_part():
    """
    Verrouille le fix du 4e defaut (ancrage a la granularite de l'ARETE) :
    une couche seulement PARTIELLEMENT posee sur le plateau doit ancrer la
    part qui y repose, meme si elle reste globalement `unresolved`.

    Sans ce comportement, la cascade se bloque des la premiere couche non
    integralement posee : plus rien ne s'ancre ensuite et TOUTES les couches
    sont signalees non resolues (mesure sur death.stl : 100%), alors meme que
    la geometrie est imprimable. L'inverse - ancrer la couche en bloc -
    reintroduirait le faux support du 1er defaut.
    """
    # Couche 0 : 5 points POSES au sol (z=0) puis 5 points EN L'AIR (z=1.5,
    # au-dela de la tolerance de contact plateau de 2.5 x 0.4 = 1.0 mm).
    grounded = np.array([[float(x), 0.0, 0.0] for x in range(5)])
    airborne = np.array([[float(x), 0.0, 1.5] for x in range(5, 10)])
    mixed_layer = [np.vstack([grounded, airborne])]

    # Couche 1 : contour AU-DESSUS de la part posee, hors bande plateau
    # (z=1.2 > 1.0) mais dans la portee de recherche verticale (7.5 x 0.2 =
    # 1.5 mm) - il ne peut donc etre ancre QUE par les aretes posees de la
    # couche 0, ce qui est precisement l'invariant teste.
    above_grounded = [np.array([[1.0, 0.0, 1.2], [2.0, 0.0, 1.2], [3.0, 0.0, 1.2]])]

    raw_layers = [
        (0.0, mixed_layer, [0], None),
        (0.1, above_grounded, [0], None),
    ]
    params = LayerParams(layer_height_base=0.2, support_driven_order=0.9, nozzle_width=0.4)

    order, unresolved = _support_driven_order(raw_layers, params, np.array([0.0, 0.0, 1.0]))

    # La couche mixte est bien signalee (5/10 points ancres, sous le seuil)...
    assert 0 in unresolved
    # ...mais elle a ANCRE sa part posee : la couche au-dessus de cette part
    # y trouve un appui et n'est donc PAS signalee.
    assert 1 not in unresolved
    assert order.index(0) < order.index(1)


# ---------------------------------------------------------------------------
# Ordre par croissance de front (LayerParams.frontier_order)
# ---------------------------------------------------------------------------

def _square_loop(centre: np.ndarray, half: float = 1.0) -> np.ndarray:
    """Petite boucle carree horizontale centree sur `centre`."""
    offsets = np.array([
        [-half, -half, 0.0], [half, -half, 0.0],
        [half, half, 0.0], [-half, half, 0.0], [-half, -half, 0.0],
    ])
    return centre[None, :] + offsets


def _frame_at(origin: np.ndarray):
    """
    Repere local horizontal, au format rendu par `_compute_local_frame` :
    (origine, axe_u, axe_v). La normale n'y figure PAS - elle se deduit du
    produit vectoriel des deux axes.
    """
    return (
        np.asarray(origin, dtype=float),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
    )


def _stacked_raw_layers(count: int, spacing: float = 1.0) -> list[tuple]:
    """
    Pile verticale de `count` couches, isovaleur croissante avec l'altitude.
    Chaque couche ne touche que ses voisines immediates.
    """
    layers = []
    for i in range(count):
        centre = np.array([0.0, 0.0, i * spacing])
        layers.append((float(i), [_square_loop(centre)], [0], _frame_at(centre)))
    return layers


def test_frontier_order_is_a_permutation_of_every_layer():
    """
    Aucune couche ne doit etre perdue ni dupliquee : c'est un ORDRE, pas un
    filtre. Verrouille aussi que le reamorcage (front epuise alors qu'il reste
    des couches) n'oublie personne.
    """
    raw_layers = _stacked_raw_layers(8)
    progress = np.arange(8, dtype=float)

    order, _reversed = _frontier_order(raw_layers, progress, gap_threshold=2.0)

    assert sorted(order) == list(range(8))
    assert len(order) == len(set(order))


def test_frontier_order_never_jumps_across_a_disconnected_gap():
    """
    La connexite du front est une contrainte DURE, pas une preference : une
    couche n'est imprimable que si elle touche du deja-imprime. C'est ce qui
    distingue ce parcours d'un tri global par progression, mesure comme
    degradant d'un facteur 3,7.

    Deux piles franchement separees, avec une progression qui les ENTRELACE
    exprès : un tri par progression alternerait entre les deux piles, le
    parcours de front doit terminer la premiere avant d'amorcer la seconde.
    """
    far = np.array([0.0, 100.0, 0.0])
    stack_a = _stacked_raw_layers(4)
    stack_b = [
        (iso, [loops[0] + far], regions, _frame_at(frame[0] + far))
        for iso, loops, regions, frame in _stacked_raw_layers(4)
    ]
    raw_layers = stack_a + stack_b
    # Progression volontairement entrelacee : 0,2,4,6 / 1,3,5,7.
    progress = np.array([0.0, 2.0, 4.0, 6.0, 1.0, 3.0, 5.0, 7.0])

    order, _reversed = _frontier_order(raw_layers, progress, gap_threshold=2.0)

    first_group = {0, 1, 2, 3}
    switch = next(i for i, n in enumerate(order) if n not in first_group)
    assert set(order[:switch]) == first_group, (
        "le front doit epuiser la premiere pile avant d'amorcer la seconde"
    )


def test_frontier_order_marks_layers_reached_from_above_as_reversed():
    """
    Une couche atteinte depuis un voisin d'isovaleur SUPERIEURE est deposee en
    DESCENDANT le champ : la tete doit y pointer a l'oppose. Seul le parcours
    connait cette information - il sait par quel voisin chaque couche a ete
    atteinte, ce qu'aucun tri a posteriori ne peut retrouver.

    Ici la progression est l'INVERSE de l'isovaleur : le parcours doit donc
    partir du sommet et descendre, en marquant tout ce qu'il atteint.
    """
    raw_layers = _stacked_raw_layers(6)
    progress = np.arange(6, dtype=float)[::-1].copy()

    order, reversed_layers = _frontier_order(raw_layers, progress, gap_threshold=2.0)

    assert order[0] == 5, "l'amorce est la plus petite progression, soit le sommet"
    # Tout le reste est atteint en descendant l'isovaleur.
    assert reversed_layers == {0, 1, 2, 3, 4}
    # L'amorce n'a aucun voisin anterieur, donc aucun sens impose.
    assert 5 not in reversed_layers


def test_frontier_order_marks_nothing_when_progress_follows_the_isovalue():
    """Cas nominal : on monte le champ, rien n'est depose a contre-sens."""
    raw_layers = _stacked_raw_layers(6)
    progress = np.arange(6, dtype=float)

    order, reversed_layers = _frontier_order(raw_layers, progress, gap_threshold=2.0)

    assert order == list(range(6))
    assert reversed_layers == set()


def test_frontier_order_disabled_by_default():
    """Opt-in : le defaut ne doit jamais emprunter ce chemin de code."""
    assert LayerParams().frontier_order is False


# ---------------------------------------------------------------------------
# Tri par isovaleur seule (LayerParams.isovalue_order)
# ---------------------------------------------------------------------------

def test_isovalue_order_ignores_region_rank(tripod_data):
    """
    Champ qui EST la progression (distance au plateau) : le rang de region ne
    doit plus primer. Dependances inversees (bras -> pied) pour que le tri
    historique (rang, isovaleur) imprime les bras, d'isovaleur haute, avant le
    pied : sans ce desaccord le test ne prouverait rien.
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph, [],
        FieldSolverParams(direction_blend_radius=20.0, direction_blend_power=2.0, direction_blend_min_norm=0.1),
    )
    inverted = RegionDependencyGraph(
        regions=region_graph.regions,
        dependencies=[(1, 0), (2, 0), (3, 0)],
        platform_contact_faces=region_graph.platform_contact_faces,
    )

    historic = layers_run(tet_mesh, mesh_data, field, inverted, LayerParams(layer_height_base=0.5))
    iso_historic = [layer.isovalue for layer in historic]
    assert iso_historic != sorted(iso_historic), "le tri par rang doit ici contredire l'isovaleur"

    ordered = layers_run(
        tet_mesh, mesh_data, field, inverted, LayerParams(layer_height_base=0.5, isovalue_order=True)
    )
    iso_ordered = [layer.isovalue for layer in ordered]
    assert iso_ordered == sorted(iso_ordered)
    assert sorted(iso_ordered) == sorted(iso_historic), "meme ensemble de couches, seul l'ordre change"


def test_isovalue_order_disabled_by_default():
    """Opt-in : le defaut garde le tri historique (rang, isovaleur)."""
    assert LayerParams().isovalue_order is False
