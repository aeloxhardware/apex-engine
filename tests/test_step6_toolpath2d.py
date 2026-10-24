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

from dataclasses import replace

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
    _generate_infill,
    _generate_spiral_infill,
    _interpolate_out_of_plane,
    _offset_walls,
    _lift_to_3d,
    _oscillating_spiral_polyline,
    _plate_axis_in_local_2d,
    _merge_short_runs,
    _monotonic_runs,
    _out_of_plane_offset,
    _project_to_local_2d,
    _split_ring_into_runs,
    _split_rings_in_toolpath,
    _to_toolpath_3col,
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

    from shapely.geometry import Point, Polygon as ShapelyPolygon

    # STRICTEMENT a l'interieur, a un micrometre pres : depuis que les trous
    # recoivent leur contour (n=0, trace SUR le bord de la cavite), des
    # points sont exactement sur la frontiere, ou un test point-dans-polygone
    # est indecis. Etre sur le bord n'est pas une fuite dans le vide.
    cavity_interior = ShapelyPolygon(cavity_open).buffer(-1e-3)
    for polyline in layer.toolpath_2d:
        for point in polyline:
            # toolpath_2d porte (u, v, offset hors-plan) depuis le fix du
            # debordement (CLAUDE.md) - le test est purement 2D, ne doit
            # voir que (u, v).
            assert not cavity_interior.contains(Point(point[:2])), (
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

    # motif LINES explicite : le test isole l'infill par ses polylignes de 2
    # points (voir ci-dessous), ce que le zigzag par defaut ne produit pas
    params = Step6Params(
        nozzle_width=0.4, wall_count=1, infill_density=0.3,
        bridge_anchor_extra_walls=1, infill_pattern=InfillPattern.LINES,
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
        # traverser plutot que longer son vide. poly[:,:2] : toolpath_2d
        # porte (u, v, offset hors-plan) depuis le fix du debordement
        # (CLAUDE.md), cette comparaison de direction est purement 2D.
        line_dirs = np.array([poly[1, :2] - poly[0, :2] for poly in infill_lines])
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


def test_classify_contours_hole_directly_nested_in_another_hole():
    """
    Reproduit fidelement le vrai bug mesure (KeyError sur
    res/test_models/death.stl, extrait par
    scripts/extract_nested_hole_case.py) : un trou dont l'UNIQUE conteneur
    est lui-meme un trou, avec AUCUN contour de matiere dans sa chaine de
    conteneurs (`containers = {0: [], 1: [0], 2: [1]}`,
    `is_hole = {0: False, 1: True, 2: True}` - exactement le motif observe
    sur la vraie couche fautive de death.stl, boucle 1 = 38 points, boucle
    2 = 69 points, extraits puis reproduits ici avec une geometrie
    synthetique plus simple mais topologiquement IDENTIQUE, verifiee par
    reproduction directe du KeyError avec l'ancienne logique non filtree
    avant d'ecrire ce test).

    Mecanisme geometrique : un anneau non-convexe (`hole_ring`, un trou)
    dont un lobe depasse legerement le contour exterieur (`outer`, non
    convexe lui-meme sur une isosurface reelle) - le point de sonde de
    `hole_ring` reste a l'interieur d'`outer` (1 conteneur, impair, trou),
    mais un second trou (`inner_hole`) place dans le lobe qui depasse a son
    point de sonde a l'interieur de `hole_ring` SEUL, pas d'`outer` (1
    conteneur = `hole_ring`, impair, trou) - sa chaine de conteneurs ne
    contient donc JAMAIS de contour matiere, contrairement a une
    imbrication concentrique propre. La version precedente de
    `_classify_contours` prenait `min(containers[i], key=areas)` SANS
    filtrer par `is_hole`, plantant (`KeyError`) des que ce plus petit
    (et unique) conteneur n'etait pas une cle valide de `groups` (qui n'a
    de cle que pour les contours matiere).
    """
    outer = _circle_loop(10.0, n=64)

    hole_ring = _circle_loop(8.0, n=48)
    # fait deborder un lobe de l'anneau au-dela du rayon 10 (non-convexite),
    # comme observe sur la vraie isosurface fautive
    poke_mask = np.abs(np.arctan2(hole_ring[:, 1], hole_ring[:, 0])) < 0.3
    hole_ring[poke_mask] *= 1.5
    # le point de sonde (index 0) doit rester DANS `outer`, pas dans le
    # lobe qui deborde - decale l'origine du parcours de la boucle
    hole_ring = np.roll(hole_ring, 24, axis=0)

    # place dans le lobe debordant, hors d'`outer` mais dans `hole_ring`
    inner_hole = _circle_loop(0.3, n=16) + np.array([11.0, 0.0])

    # Aucun KeyError : c'est l'assertion principale de ce test (la version
    # buguee levait avant meme d'atteindre les asserts suivants).
    groups = _classify_contours([outer, hole_ring, inner_hole])

    # 1 seul groupe matiere attendu (`outer`), avec `hole_ring` comme son
    # trou - `inner_hole` (chaine de conteneurs sans matiere) est ignore,
    # meme discipline defensive que step5 (`_map_tets_to_regions`) : un
    # trou topologiquement incoherent est ecarte plutot que de planter.
    assert len(groups) == 1
    outer_group = groups[0]
    assert np.allclose(outer_group[0][0], outer[0])
    assert len(outer_group[1]) == 1


def test_classify_contours_two_disjoint_islands():
    """2 boucles disjointes (pas imbriquees) : chacune est sa propre matiere sans trou."""
    island_a = _circle_loop(3.0) + np.array([20.0, 0.0])
    island_b = _circle_loop(3.0) + np.array([-20.0, 0.0])

    groups = _classify_contours([island_a, island_b])

    assert len(groups) == 2
    for _outer, holes in groups:
        assert holes == []


# ---------------------------------------------------------------------------
# Cas 6 - Fix du debordement (CLAUDE.md) : offset hors-plan preserve au
# relevage du contour d'origine (n=0), sur un plan local INCLINE - la
# perte se produit precisement quand normal != axe Z, jamais visible sur
# un plan horizontal simple.
# ---------------------------------------------------------------------------

def test_project_lift_roundtrip_preserves_out_of_plane_points():
    """
    Unitaire direct : un point 3D hors du plan local, projete puis releve
    AVEC son offset hors-plan preserve (_to_toolpath_3col + _lift_to_3d),
    doit redonner EXACTEMENT le point d'origine - contrairement au
    comportement 2 colonnes (sans offset), qui l'aplatit sur le plan
    moyen. Plan local INCLINE (pas seulement horizontal) pour rester
    representatif du cas reel (couche a 55 degres, cf. CLAUDE.md).
    """
    origin = np.array([1.0, -2.0, 0.5])
    axis_u = np.array([0.6, 0.8, 0.0])
    axis_v = np.array([-0.48, 0.36, 0.8])  # orthonormee a axis_u par construction

    assert abs(np.dot(axis_u, axis_v)) < 1e-9
    assert abs(np.linalg.norm(axis_u) - 1.0) < 1e-9
    assert abs(np.linalg.norm(axis_v) - 1.0) < 1e-9

    local_frame = (origin, axis_u, axis_v)

    rng = np.random.default_rng(42)
    points_3d = origin + rng.normal(size=(20, 3)) * 3.0

    p2d = _project_to_local_2d(points_3d, local_frame)
    offset = _out_of_plane_offset(points_3d, local_frame)

    # au moins une partie de l'echantillon aleatoire doit reellement
    # devier du plan (sinon le test ne verifierait rien) :
    assert np.abs(offset).max() > 0.1

    p3col = _to_toolpath_3col(p2d, offset)
    roundtrip = _lift_to_3d(p3col, local_frame)

    assert np.allclose(roundtrip, points_3d, atol=1e-9), (
        "L'aller-retour projection+relevage avec offset hors-plan devrait "
        "redonner exactement les points d'origine"
    )


def test_project_lift_roundtrip_without_offset_flattens_to_plane():
    """
    Non-regression du comportement HISTORIQUE : sans offset (2 colonnes,
    parois offsettees/infill), le relevage retombe sur le plan moyen, pas
    sur les points d'origine - c'est le comportement attendu pour ces
    points 2D nouveaux, sans equivalent 3D d'origine (cf. CLAUDE.md).
    """
    origin = np.array([0.0, 0.0, 0.0])
    axis_u = np.array([1.0, 0.0, 0.0])
    axis_v = np.array([0.0, 1.0, 0.0])
    local_frame = (origin, axis_u, axis_v)

    points_3d = np.array([[1.0, 2.0, 5.0], [3.0, -1.0, -2.0]])
    p2d = _project_to_local_2d(points_3d, local_frame)
    roundtrip = _lift_to_3d(p2d, local_frame)

    assert np.allclose(roundtrip[:, 2], 0.0), "Sans offset, le relevage doit rester sur le plan (z=0 ici)"
    assert not np.allclose(roundtrip, points_3d), (
        "Sans offset, le relevage NE doit PAS retomber sur les points d'origine hors-plan "
        "(comportement attendu, pas une regression)"
    )


def test_run_outer_contour_offset_matches_original_points(flat_box_layers):
    """
    Test d'integration : dans la sortie de run(), la polyligne n=0 (contour
    exterieur, derniere de chaque groupe issu de _offset_walls, cf.
    step6_toolpath2d._offset_walls) doit, une fois relevee en 3D via
    _lift_to_3d, retomber sur les points reels de Layer.contour_3d - pas
    seulement sur le plan moyen approximatif de step5.
    """
    mesh_data, overhang_map, layers = flat_box_layers
    params = Step6Params(nozzle_width=0.4, wall_count=1, infill_density=0.2)

    result = toolpath_run(layers, mesh_data, overhang_map, params)

    checked = 0
    for layer in result:
        if not layer.toolpath_2d:
            continue
        contour_3d = layer.contour_3d[0]  # cube plein, un seul contour par couche
        contour_2d = _project_to_local_2d(contour_3d[:-1], layer.local_frame)

        # la polyligne n=0 est celle dont la projection correspond au
        # contour d'origine (meme technique que step6.run()/_loop_matches)
        for polyline in layer.toolpath_2d:
            if polyline.shape[0] != contour_2d.shape[0]:
                continue
            if not np.allclose(polyline[:, :2], contour_2d, atol=1e-6):
                continue
            lifted = _lift_to_3d(polyline, layer.local_frame)
            assert np.allclose(lifted, contour_3d[:-1], atol=1e-6), (
                f"Couche {layer.id}: le contour n=0 releve en 3D devrait "
                f"retomber exactement sur Layer.contour_3d"
            )
            checked += 1

    assert checked > 0, "Aucune polyligne n=0 trouvee - le test n'a rien verifie"


# ---------------------------------------------------------------------------
# Cas 7 - Offset hors-plan INTERPOLE pour les points generes par shapely
# (parois n>=1, infill). Complete le fix de la limite #31, qui ne couvrait
# que le contour n=0 : mesure sur le modele en T, 74 points d'infill/parois
# internes tombaient hors du solide (max 0.39 mm) parce qu'ils restaient
# colles au plan local MOYEN d'une couche pourtant non plane.
# ---------------------------------------------------------------------------

def test_interpolate_reproduces_source_value_at_coincident_point():
    """Point de requete confondu avec une source : valeur source EXACTE,
    jamais une division par une distance quasi nulle."""
    source_2d = np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0]])
    source_offsets = np.array([1.0, -2.0, 5.0])

    result = _interpolate_out_of_plane(source_2d, source_2d, source_offsets, power=2.0)

    assert np.allclose(result, source_offsets)


def test_interpolate_constant_field_stays_constant():
    """Toutes les sources a la meme valeur : l'interpolation doit rendre
    exactement cette valeur, sans derive numerique."""
    source_2d = np.array([[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]])
    source_offsets = np.full(4, 0.75)
    query = np.array([[2.0, 2.0], [1.0, 3.0], [3.5, 0.5]])

    result = _interpolate_out_of_plane(query, source_2d, source_offsets, power=2.0)

    assert np.allclose(result, 0.75)


def test_interpolate_result_always_bounded_by_sources():
    """
    Propriete d'une moyenne ponderee (et non "exacte sur un champ lineaire",
    que l'IDW ne garantit PAS) : le mecanisme ne peut jamais inventer un
    offset plus extreme que ceux reellement mesures sur le contour.
    """
    rng = np.random.default_rng(7)
    source_2d = rng.normal(size=(24, 2)) * 5.0
    source_offsets = rng.normal(size=24) * 2.0
    query = rng.normal(size=(60, 2)) * 5.0

    result = _interpolate_out_of_plane(query, source_2d, source_offsets, power=2.0)

    assert result.min() >= source_offsets.min() - 1e-12
    assert result.max() <= source_offsets.max() + 1e-12


def test_interpolate_single_source_is_constant():
    """Source unique : constante, pas de division par zero."""
    result = _interpolate_out_of_plane(
        np.array([[1.0, 1.0], [5.0, -3.0]]), np.array([[0.0, 0.0]]), np.array([2.5]), power=2.0
    )
    assert np.allclose(result, 2.5)


def test_interpolate_without_sources_returns_zeros():
    """Aucune source (couche sans boucle appariee) : zeros, pas d'exception -
    meme repli silencieux que le comportement historique."""
    result = _interpolate_out_of_plane(
        np.array([[1.0, 1.0], [2.0, 2.0]]), np.zeros((0, 2)), np.zeros(0), power=2.0
    )
    assert result.shape == (2,)
    assert np.allclose(result, 0.0)


def _wavy_square_layer() -> Layer:
    """
    Couche synthetique franchement NON PLANE : contour carre dont le z ondule
    (amplitude 1.0) autour de son plan moyen z=0. Reproduit le mecanisme reel
    (une couche courbe aplatie sur son plan moyen) sans dependre d'un modele
    dont la non-planarite serait accidentelle.
    """
    n_per_side = 12
    t = np.linspace(0.0, 1.0, n_per_side, endpoint=False)
    side = 20.0
    half = side / 2.0
    xs = np.concatenate([-half + t * side, np.full(n_per_side, half),
                         half - t * side, np.full(n_per_side, -half)])
    ys = np.concatenate([np.full(n_per_side, -half), -half + t * side,
                         np.full(n_per_side, half), half - t * side])
    theta = np.arctan2(ys, xs)
    zs = 1.0 * np.sin(2.0 * theta)

    loop = np.stack([xs, ys, zs], axis=1)
    loop = np.vstack([loop, loop[0]])  # point de fermeture duplique (convention step5)

    local_frame = (np.zeros(3), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    return Layer(id=0, contour_3d=[loop], region_ids=[0], local_frame=local_frame, isovalue=0.0)


def test_generated_points_get_zero_offset_when_disabled(flat_box_layers):
    """
    Comportement historique conserve quand le mecanisme est desactive
    (defaut None) : tout ce qui n'est pas le contour n=0 garde un offset nul.
    """
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _wavy_square_layer()
    params = Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.3)
    assert params.out_of_plane_interpolation_power is None

    result = toolpath_run([layer], mesh_data, overhang_map, params)
    polylines = result[0].toolpath_2d
    assert polylines

    contour_2d = _project_to_local_2d(layer.contour_3d[0][:-1], layer.local_frame)
    generated = [
        p for p in polylines
        if not (p.shape[0] == contour_2d.shape[0] and np.allclose(p[:, :2], contour_2d, atol=1e-9))
    ]
    assert generated, "La fixture doit produire des parois/infill, sinon le test ne verifie rien"
    for poly in generated:
        assert np.allclose(poly[:, 2], 0.0)


def test_generated_points_follow_curved_layer_when_enabled(flat_box_layers):
    """
    Mecanisme actif : les points generes recoivent un offset hors-plan non
    nul, interpole depuis le contour, et jamais hors de l'enveloppe des
    valeurs reellement mesurees sur ce contour.
    """
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _wavy_square_layer()
    params = Step6Params(
        nozzle_width=0.4, wall_count=2, infill_density=0.3,
        out_of_plane_interpolation_power=2.0,
    )

    result = toolpath_run([layer], mesh_data, overhang_map, params)
    polylines = result[0].toolpath_2d

    contour_offsets = _out_of_plane_offset(layer.contour_3d[0][:-1], layer.local_frame)
    contour_2d = _project_to_local_2d(layer.contour_3d[0][:-1], layer.local_frame)
    generated = [
        p for p in polylines
        if not (p.shape[0] == contour_2d.shape[0] and np.allclose(p[:, :2], contour_2d, atol=1e-9))
    ]
    assert generated

    all_offsets = np.concatenate([p[:, 2] for p in generated])
    assert np.abs(all_offsets).max() > 0.05, (
        "Sur une couche franchement non plane, les points generes devraient "
        "recevoir un offset hors-plan non nul"
    )
    assert all_offsets.min() >= contour_offsets.min() - 1e-9
    assert all_offsets.max() <= contour_offsets.max() + 1e-9


def test_interpolation_preserves_exact_outer_contour_offsets(flat_box_layers):
    """
    Non-regression du fix #31 : avec l'interpolation ACTIVE, le contour n=0
    doit toujours retomber EXACTEMENT sur contour_3d - ses offsets sont
    mesures, jamais interpoles.
    """
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _wavy_square_layer()
    params = Step6Params(
        nozzle_width=0.4, wall_count=2, infill_density=0.3,
        out_of_plane_interpolation_power=2.0,
    )

    result = toolpath_run([layer], mesh_data, overhang_map, params)
    contour_3d = layer.contour_3d[0][:-1]
    contour_2d = _project_to_local_2d(contour_3d, layer.local_frame)

    checked = 0
    for poly in result[0].toolpath_2d:
        if poly.shape[0] != contour_2d.shape[0] or not np.allclose(poly[:, :2], contour_2d, atol=1e-9):
            continue
        assert np.allclose(_lift_to_3d(poly, layer.local_frame), contour_3d, atol=1e-9)
        checked += 1

    assert checked > 0, "Contour n=0 introuvable - le test n'a rien verifie"


def test_interpolation_leaves_2d_geometry_untouched(flat_box_layers):
    """
    L'interpolation ne touche QUE la 3e colonne : meme nombre de polylignes,
    memes coordonnees (u,v), meme ordre du centre vers l'exterieur.
    """
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _wavy_square_layer()
    base = Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.3)
    interpolated = Step6Params(
        nozzle_width=0.4, wall_count=2, infill_density=0.3,
        out_of_plane_interpolation_power=2.0,
    )

    without = toolpath_run([layer], mesh_data, overhang_map, base)[0].toolpath_2d
    with_ = toolpath_run([layer], mesh_data, overhang_map, interpolated)[0].toolpath_2d

    assert len(without) == len(with_)
    for a, b in zip(without, with_):
        assert np.allclose(a[:, :2], b[:, :2], atol=1e-12)


# ---------------------------------------------------------------------------
# Decision 5 - reordonnancement intra-anneau pour eviter les depots sans
# support sur une couche courbe (voir CLAUDE.md / step6_toolpath2d.py)
# ---------------------------------------------------------------------------

def _unsupported_indices(ring: np.ndarray, horizontal_radius: float, vertical_height: float) -> list[int]:
    """
    Meme geometrie que le test de support de step9b/step5 (segments
    degeneres a un point), appliquee directement a un anneau (K,3) deja
    depose dans l'ordre - sert a MESURER l'effet du reordonnancement dans
    ces tests, independamment de l'implementation interne.
    """
    unsupported = []
    for i in range(1, len(ring)):
        pos = ring[i]
        prior = ring[:i]
        dz = pos[2] - prior[:, 2]
        horiz = np.hypot(pos[0] - prior[:, 0], pos[1] - prior[:, 1])
        if not np.any((dz > 1e-9) & (dz <= vertical_height) & (horiz <= horizontal_radius)):
            unsupported.append(i)
    return unsupported


def _wavy_ring(n: int, radius: float, amplitude: float, cycles: int) -> np.ndarray:
    """Anneau ferme (K+1,3) dont l'offset (3e colonne) ondule sinusoidalement -
    pas de plateaux (jamais deux offsets consecutifs egaux), point d'espacement
    dense (nozzle-scale) pour que le mecanisme ait une chance geometrique
    reelle de fonctionner (cf. diagnostic : sur un anneau grossier, l'espacement
    entre sommets depasse deja le rayon de recherche independamment de l'ordre)."""
    theta = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    u = radius * np.cos(theta)
    v = radius * np.sin(theta)
    offset = amplitude * np.sin(theta * cycles)
    ring = np.stack([u, v, offset], axis=1)
    return np.vstack([ring, ring[0]])


def test_monotonic_runs_single_run_for_flat_ring():
    """Un anneau plat (comportement historique sans interpolation, offset
    partout nul) ne produit qu'UN SEUL troncon - rien a reordonner."""
    ring = _wavy_ring(n=20, radius=5.0, amplitude=0.0, cycles=3)
    runs = _monotonic_runs(ring[:, 2])
    assert len(runs) == 1
    assert runs[0] == (0, len(ring))


def test_monotonic_runs_partition_without_overlap():
    """Les troncons sont des intervalles [start, stop) DISJOINTS qui
    partitionnent exactement l'anneau - aucun sommet partage ni oublie (bug
    reel trouve et corrige : une version precedente partageait un sommet a
    chaque frontiere, dupliquant des points une fois les troncons reordonnes)."""
    ring = _wavy_ring(n=64, radius=10.0, amplitude=2.0, cycles=6)
    runs = _monotonic_runs(ring[:, 2])
    assert runs[0][0] == 0
    assert runs[-1][1] == len(ring)
    for (_, stop_a), (start_b, _) in zip(runs, runs[1:]):
        assert stop_a == start_b
    total = sum(stop - start for start, stop in runs)
    assert total == len(ring)


def test_split_ring_preserves_every_point():
    """Le decoupage ne peut ni ajouter ni supprimer de points : la
    concatenation des troncons contient exactement le meme multi-ensemble
    d'offsets que l'anneau d'origine (bug reel deja mesure : une version
    precedente partageait un sommet a chaque frontiere et ressortait avec
    PLUS de points que l'entree)."""
    ring = _wavy_ring(n=200, radius=25.0, amplitude=2.0, cycles=6)
    runs = _split_ring_into_runs(ring)

    total = sum(len(r) for r in runs)
    assert total == len(ring)
    concatenated = np.concatenate(runs, axis=0)
    assert np.allclose(sorted(ring[:, 2].tolist()), sorted(concatenated[:, 2].tolist()))


def test_split_ring_runs_are_all_ascending():
    """
    Coeur du fix (Decision 5) : chaque troncon ressort parcouru dans le sens
    MONTANT (un troncon descendant est inverse). C'est ce qui le rend
    auto-portant - chaque point s'appuie sur celui qui vient d'etre depose
    juste avant lui, plus bas et proche horizontalement.
    """
    ring = _wavy_ring(n=400, radius=25.0, amplitude=2.0, cycles=6)
    runs = _split_ring_into_runs(ring)

    assert len(runs) > 1, "la fixture doit vraiment onduler, sinon le test ne verifie rien"
    for run in runs:
        if len(run) >= 2:
            assert run[-1, 2] >= run[0, 2]


def test_split_ring_runs_have_at_least_two_points():
    """
    Un troncon d'UN SEUL point deviendrait une polyligne d'un point : step8 y
    insererait un trajet a vide puis un point isole, et
    step9._grow_printed_index_to (qui exige DEUX extremites en extrusion) n'en
    ferait aucun cordon - ce point cesserait d'etre imprime. `_merge_short_runs`
    doit donc les fusionner dans leur voisin, sans jamais perdre de point.
    """
    # ondulation a haute frequence : force des troncons tres courts
    ring = _wavy_ring(n=40, radius=5.0, amplitude=1.0, cycles=9)
    runs = _split_ring_into_runs(ring)

    assert sum(len(r) for r in runs) == len(ring)
    for run in runs:
        assert len(run) >= 2


def test_split_ring_noop_on_flat_ring():
    """Un anneau plat (un seul troncon) ressort en UNE polyligne identique -
    non-regression stricte du comportement historique (offset nul)."""
    ring = _wavy_ring(n=30, radius=5.0, amplitude=0.0, cycles=4)
    runs = _split_ring_into_runs(ring)
    assert len(runs) == 1
    assert np.allclose(runs[0], ring)


def test_split_rings_in_toolpath_skips_open_polylines():
    """Une ligne d'infill OUVERTE n'est jamais un anneau et ressort
    strictement inchangee : une droite clippee n'a pas de "retour au point de
    depart" a gerer."""
    open_line = np.array([[0.0, 0.0, 2.0], [1.0, 0.0, 1.0], [2.0, 0.0, -2.0], [3.0, 0.0, 0.0]])
    closed_ring = _wavy_ring(n=64, radius=10.0, amplitude=2.0, cycles=6)

    result = _split_rings_in_toolpath([open_line, closed_ring])

    assert np.allclose(result[0], open_line)
    assert len(result) > 2, "l'anneau ondule doit avoir ete decoupe en plusieurs troncons"
    rebuilt = np.concatenate(result[1:], axis=0)
    assert np.allclose(sorted(rebuilt[:, 2].tolist()), sorted(closed_ring[:, 2].tolist()))


def test_split_disabled_by_default_is_bit_identical(flat_box_layers):
    """Comportement historique inchange quand `split_rings_into_runs` reste a
    son defaut (False) - meme discipline opt-in que les autres mecanismes de
    ce fichier (out_of_plane_interpolation_power, etc.)."""
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _wavy_square_layer()
    params = Step6Params(
        nozzle_width=0.4, wall_count=2, infill_density=0.3,
        out_of_plane_interpolation_power=2.0,
    )

    default = toolpath_run([layer], mesh_data, overhang_map, params)[0].toolpath_2d
    explicit_false = toolpath_run(
        [layer], mesh_data, overhang_map, replace(params, split_rings_into_runs=False)
    )[0].toolpath_2d

    assert len(default) == len(explicit_false)
    for a, b in zip(default, explicit_false):
        assert np.allclose(a, b)


def _dense_wavy_square_layer() -> Layer:
    """
    Variante DENSE de `_wavy_square_layer` (48 points, espacement 1.67mm -
    superieur a horizontal_radius=1.0, donc structurellement incapable de
    demontrer l'effet du mecanisme). Meme forme (carre, ondulation
    sinusoidale d'amplitude 1.0 autour du plan moyen), mais espacement proche
    de nozzle_width, pour que le test d'integration porte reellement sur
    l'effet du decoupage et non sur une resolution de fixture insuffisante.
    """
    n_per_side = 60
    t = np.linspace(0.0, 1.0, n_per_side, endpoint=False)
    side = 20.0
    half = side / 2.0
    xs = np.concatenate([-half + t * side, np.full(n_per_side, half),
                         half - t * side, np.full(n_per_side, -half)])
    ys = np.concatenate([np.full(n_per_side, -half), -half + t * side,
                         np.full(n_per_side, half), half - t * side])
    theta = np.arctan2(ys, xs)
    zs = 1.0 * np.sin(2.0 * theta)

    loop = np.stack([xs, ys, zs], axis=1)
    loop = np.vstack([loop, loop[0]])

    local_frame = (np.zeros(3), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    return Layer(id=0, contour_3d=[loop], region_ids=[0], local_frame=local_frame, isovalue=0.0)


def test_split_increases_polyline_count_without_losing_points(flat_box_layers):
    """
    Test d'integration : sur une couche NON PLANE dense, le decoupage doit
    produire STRICTEMENT PLUS de polylignes (chaque frontiere de troncon
    devient une vraie frontiere de segment, donc un vrai trajet a vide insere
    par step8._iter_schedule_points au lieu d'un cordon parasite extrude en
    travers de la piece) sans jamais perdre ni dupliquer un point.
    """
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _dense_wavy_square_layer()

    base_params = Step6Params(
        nozzle_width=0.4, wall_count=2, infill_density=0.3,
        out_of_plane_interpolation_power=2.0, split_rings_into_runs=False,
    )
    split_params = replace(base_params, split_rings_into_runs=True)

    before = toolpath_run([layer], mesh_data, overhang_map, base_params)[0].toolpath_2d
    after = toolpath_run([layer], mesh_data, overhang_map, split_params)[0].toolpath_2d

    assert len(after) > len(before)
    assert sum(len(p) for p in after) == sum(len(p) for p in before)


def test_split_reduces_unsupported_points_end_to_end(flat_box_layers):
    """
    Meme fixture, mesure du VRAI effet : en comptant les points deposes sans
    rien en dessous (avec un trajet a vide a chaque frontiere de polyligne -
    la geometrie que step8 produira reellement), le decoupage en troncons
    montants doit faire baisser le compte.
    """
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _dense_wavy_square_layer()
    horizontal_radius, vertical_height = 1.0, 1.5

    base_params = Step6Params(
        nozzle_width=0.4, wall_count=2, infill_density=0.3,
        out_of_plane_interpolation_power=2.0, split_rings_into_runs=False,
    )
    split_params = replace(base_params, split_rings_into_runs=True)

    def _count_unsupported(toolpath_2d):
        positions, extrusions = [], []
        for seg_idx, poly in enumerate(toolpath_2d):
            pts_3d = _lift_to_3d(poly, layer.local_frame)
            if seg_idx > 0:
                positions.append(pts_3d[0])
                extrusions.append(0.0)   # trajet a vide, comme step8
            for p in pts_3d:
                positions.append(p)
                extrusions.append(1.0)
        positions = np.array(positions)
        extrusions = np.array(extrusions)

        unsupported = 0
        for i in range(1, len(positions)):
            if extrusions[i] <= 0.0:
                continue
            pos = positions[i]
            prior = positions[:i][extrusions[:i] > 0.0]
            if len(prior) == 0:
                unsupported += 1
                continue
            dz = pos[2] - prior[:, 2]
            horiz = np.hypot(pos[0] - prior[:, 0], pos[1] - prior[:, 1])
            if not np.any((dz > 1e-9) & (dz <= vertical_height) & (horiz <= horizontal_radius)):
                unsupported += 1
        return unsupported

    before = _count_unsupported(toolpath_run([layer], mesh_data, overhang_map, base_params)[0].toolpath_2d)
    after = _count_unsupported(toolpath_run([layer], mesh_data, overhang_map, split_params)[0].toolpath_2d)

    assert after < before




# ---------------------------------------------------------------------------
# Direction de fusion des troncons courts (_merge_short_runs)
#
# Cause racine MESUREE de 425 des 1055 defauts d'ordre sur death.stl (40.3%) :
# la fusion arriere systematique recollait un micro-creux (point isole, dz
# median -0.0141mm) a la FIN d'un troncon montant, le faisant deposer APRES un
# point plus haut que lui. Verifie : 425/425 (100%) de ces defauts avaient leur
# predecesseur immediat plus haut, DANS LA MEME polyligne continue.
# ---------------------------------------------------------------------------

def test_merge_short_runs_carries_dip_forward_not_backward():
    """Un point isole PLUS BAS que la fin du troncon precedent doit demarrer le
    troncon SUIVANT (il en devient le point le plus bas), jamais etre colle a la
    fin du precedent (il y serait depose apres plus haut que lui)."""
    #        idx: 0    1    2    3     4    5
    offsets = np.array([0.0, 1.0, 2.0, 1.9, 3.0, 4.0])
    runs = [(0, 3), (3, 4), (4, 6)]     # (3,4) = le micro-creux, 1 seul point

    merged = _merge_short_runs(runs, offsets)

    assert (0, 3) in merged, "le troncon montant ne doit pas absorber le creux"
    assert (3, 6) in merged, "le creux doit demarrer le troncon suivant"


def test_merge_short_runs_still_merges_backward_when_monotony_holds():
    """Si le point isole est PLUS HAUT que la fin du troncon precedent, la
    fusion arriere preserve la montee - comportement conserve."""
    offsets = np.array([0.0, 1.0, 2.0, 2.5, 1.0, 0.5])
    runs = [(0, 3), (3, 4), (4, 6)]

    merged = _merge_short_runs(runs, offsets)

    assert (0, 4) in merged, "2.5 >= 2.0 : la fusion arriere reste monotone"


def test_merge_short_runs_partitions_exactly_with_offsets():
    """Aucun point perdu ni duplique, quelle que soit la direction choisie."""
    offsets = np.array([0.0, 1.0, 0.9, 2.0, 3.0, 2.9, 4.0])
    runs = [(0, 2), (2, 3), (3, 5), (5, 6), (6, 7)]

    merged = _merge_short_runs(runs, offsets)

    assert merged[0][0] == 0
    assert merged[-1][1] == len(offsets)
    for (_, stop_a), (start_b, _) in zip(merged, merged[1:]):
        assert stop_a == start_b
    assert sum(stop - start for start, stop in merged) == len(offsets)


def test_merge_short_runs_legacy_path_unchanged():
    """`offsets=None` conserve strictement l'ancien comportement (fusion
    arriere systematique) - non-regression pour tout appelant sans altitude."""
    runs = [(0, 3), (3, 4), (4, 6)]
    assert _merge_short_runs(runs) == [(0, 4), (4, 6)]


# ---------------------------------------------------------------------------
# Espacement de l'infill : derive de la buse, jamais de la taille de la couche
# ---------------------------------------------------------------------------

def _square(side):
    from shapely.geometry import box
    return box(0.0, 0.0, side, side)


def test_infill_spacing_follows_nozzle_and_density():
    """Carre de 40 mm, buse 0,4, densite 0,2 : une ligne tous les 2 mm, soit
    environ 20 lignes. L'ancienne formule (fondee sur la diagonale de la
    couche) en donnait ~5 : c'est le defaut qui laissait les grandes couches
    presque vides (mesure sur les couches Bird du Benchy)."""
    lines = _generate_infill(_square(40.0), 0.2, InfillPattern.LINES, nozzle_width=0.4)
    assert 18 <= len(lines) <= 22
    ys = sorted(float(line[0][1]) for line in lines)
    assert np.allclose(np.diff(ys), 2.0, atol=1e-6)


def test_infill_spacing_does_not_depend_on_layer_size():
    """Meme densite, couches de tailles differentes : meme espacement, donc
    un nombre de lignes proportionnel a la taille."""
    small = _generate_infill(_square(10.0), 0.2, InfillPattern.LINES, nozzle_width=0.4)
    large = _generate_infill(_square(80.0), 0.2, InfillPattern.LINES, nozzle_width=0.4)
    assert abs(len(large) / len(small) - 8.0) < 1.0


def test_infill_full_density_lines_are_one_nozzle_apart():
    lines = _generate_infill(_square(10.0), 1.0, InfillPattern.LINES, nozzle_width=0.4)
    ys = sorted(float(line[0][1]) for line in lines)
    assert np.allclose(np.diff(ys), 0.4, atol=1e-6)


# ---------------------------------------------------------------------------
# Parois : tous les bords de tous les morceaux (trous et zones coupees en deux)
# ---------------------------------------------------------------------------

def _ring(xmin, ymin, xmax, ymax):
    return np.array([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]], dtype=float)


def _ring_bounds(ring):
    return ring[:, 0].min(), ring[:, 1].min(), ring[:, 0].max(), ring[:, 1].max()


def test_hole_receives_its_contour_and_walls():
    """Carre de 20 mm perce d'un trou de 8 mm, 2 parois : le trou doit etre
    borde par son contour (n=0) ET par une paroi a chaque niveau, placee
    cote matiere (trou elargi de 0,4 puis 0,8 mm). Avant la correction, seul
    le bord exterieur de chaque offset etait trace : le trou n'avait rien."""
    outer, hole = _ring(0, 0, 20, 20), _ring(6, 6, 14, 14)
    walls, _ = _offset_walls(outer, [hole], nozzle_width=0.4, wall_count=2, degenerate_area_ratio=0.01)

    around_hole = [w for w in walls if _ring_bounds(w)[2] - _ring_bounds(w)[0] < 12.0]
    widths = sorted(round(_ring_bounds(w)[2] - _ring_bounds(w)[0], 3) for w in around_hole)
    assert widths == [8.0, 8.8, 9.6], f"bords du trou attendus a 8 / 8,8 / 9,6 mm, obtenus {widths}"


def test_split_offset_keeps_every_piece_with_walls_and_infill():
    """Haltere : deux carres de 10 mm relies par un col de 0,6 mm. Des la
    premiere paroi (0,4 mm), le col disparait et la zone se coupe en deux.
    Chaque moitie doit garder ses parois et son remplissage ; avant la
    correction, seul le plus grand morceau etait garde."""
    outer = np.array([[0, 0], [10, 0], [10, 4.7], [14, 4.7], [14, 0], [24, 0],
                      [24, 10], [14, 10], [14, 5.3], [10, 5.3], [10, 10], [0, 10]], dtype=float)
    walls, innermost = _offset_walls(outer, [], nozzle_width=0.4, wall_count=2, degenerate_area_ratio=0.01)

    inner_walls = walls[:-1]              # le dernier est le contour d'origine
    left = [w for w in inner_walls if _ring_bounds(w)[2] < 12.0]
    right = [w for w in inner_walls if _ring_bounds(w)[0] > 12.0]
    assert len(left) == 2 and len(right) == 2, "2 parois attendues dans chaque moitie"

    infill = _generate_infill(innermost, 0.2, InfillPattern.LINES, nozzle_width=0.4)
    xs = np.concatenate([line[:, 0] for line in infill])
    assert (xs < 12.0).any() and (xs > 12.0).any(), "remplissage attendu dans les deux moities"


# ---------------------------------------------------------------------------
# Motifs de remplissage : zigzag croise-alterne (defaut), grille
# ---------------------------------------------------------------------------

def _segments(polyline):
    return np.diff(np.asarray(polyline)[:, :2], axis=0)


def test_zigzag_links_a_convex_zone_into_one_continuous_path():
    """Carre de 40 mm : un seul zigzag, qui contient toutes les lignes de
    balayage (~20, une tous les 2 mm) reliees par des raccords courts."""
    chains = _generate_infill(_square(40.0), 0.2, InfillPattern.CROSSED_ZIGZAG, nozzle_width=0.4)
    assert len(chains) == 1
    seg = _segments(chains[0])
    sweeps = np.abs(seg[:, 0]) > 10.0          # lignes de balayage (selon u)
    assert 18 <= sweeps.sum() <= 22
    assert np.all(np.linalg.norm(seg[~sweeps], axis=1) <= 2.5 * 2.0 + 1e-9)
    # le zigzag change de sens a chaque ligne
    assert np.all(np.sign(seg[sweeps, 0])[1:] != np.sign(seg[sweeps, 0])[:-1])


def test_zigzag_never_jumps_across_a_hole_or_a_gap():
    """Carre perce d'un grand trou, et haltere coupee en deux par la paroi :
    aucun raccord du zigzag ne traverse le vide - chaque segment reste dans
    la zone a remplir."""
    from shapely.geometry import LineString, box

    holed = box(0, 0, 30, 30).difference(box(8, 8, 22, 22))
    split = box(0, 0, 10, 10).union(box(14, 0, 24, 10))
    for zone in (holed, split):
        chains = _generate_infill(zone, 0.2, InfillPattern.CROSSED_ZIGZAG, nozzle_width=0.4)
        assert len(chains) >= 2, "une zone non convexe doit donner plusieurs zigzags"
        inside = zone.buffer(1e-6)
        for chain in chains:
            for a, b in zip(chain[:-1], chain[1:]):
                assert inside.covers(LineString([a, b]))


def test_grid_crosses_two_sweeps_in_the_same_layer():
    """GRID etait recu et ignore (lignes simples). Il produit maintenant deux
    balayages a 90 deg, chacun a la moitie de la densite."""
    lines = _generate_infill(_square(40.0), 0.2, InfillPattern.GRID, nozzle_width=0.4)
    dirs = np.array([line[-1] - line[0] for line in lines])
    along_u = np.abs(dirs[:, 0]) > np.abs(dirs[:, 1])
    assert 8 <= along_u.sum() <= 11 and 8 <= (~along_u).sum() <= 11


def _flat_square_layer(layer_id, isovalue, side=20.0):
    half = side / 2.0
    loop = np.array([[-half, -half, 0.0], [half, -half, 0.0], [half, half, 0.0],
                     [-half, half, 0.0], [-half, -half, 0.0]])
    frame = (np.zeros(3), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    return Layer(id=layer_id, contour_3d=[loop], region_ids=[0], local_frame=frame, isovalue=isovalue)


def _zigzag_axis(layer):
    """0 si le zigzag de la couche balaie selon u, 1 selon v."""
    chain = max((p for p in layer.toolpath_2d if not np.allclose(p[0], p[-1])), key=len)
    seg = np.abs(_segments(chain))
    return int(seg[:, 1].sum() > seg[:, 0].sum())


def test_zigzag_alternates_between_stacked_layers(flat_box_layers):
    """Lignes croisees-alternees : deux couches EMPILEES (isovaleurs
    consecutives) se croisent a 90 deg. La parite suit le rang de
    l'isovaleur, pas la position dans la liste (l'ordre d'impression en
    arbre fait se succeder des couches de branches differentes), et deux
    morceaux d'une meme isovaleur gardent la meme orientation."""
    mesh_data, overhang_map, _layers = flat_box_layers
    layers = [_flat_square_layer(0, 0.0), _flat_square_layer(1, 0.4),
              _flat_square_layer(2, 0.2), _flat_square_layer(3, 0.2)]
    params = Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.2)
    assert params.infill_pattern is InfillPattern.CROSSED_ZIGZAG, "motif par defaut"
    result = toolpath_run(layers, mesh_data, overhang_map, params)
    axes = [_zigzag_axis(layer) for layer in result]
    assert axes == [0, 0, 1, 1]


def test_zigzag_does_not_retrace_the_innermost_wall():
    """La zone passee a _generate_infill est le polygone de la paroi la plus
    interieure : son bord EST l'axe de cette paroi. Le zigzag doit s'en tenir
    a une demi-buse, sinon ses raccords, qui suivent le bord, repasseraient
    sur la paroi."""
    from shapely.geometry import Point

    zone = _square(20.0)
    chains = _generate_infill(zone, 0.2, InfillPattern.CROSSED_ZIGZAG, nozzle_width=0.4)
    points = np.vstack(chains)
    distances = [zone.exterior.distance(Point(p)) for p in points]
    assert min(distances) >= 0.2 - 1e-9


def test_split_ring_shared_vertices_extrude_every_segment():
    """
    `continuous_rings` (2026-10-01) : avec des troncons a indices DISJOINTS, le
    segment qui relie deux troncons n'est jamais extrude (il devient un trajet
    a vide) - 7,7 % de la longueur perdue sur death.stl. Avec le sommet de
    jonction partage, chaque segment de l'anneau appartient a un troncon, et
    chaque troncon reste montant.
    """
    ring = _wavy_ring(n=400, radius=25.0, amplitude=2.0, cycles=6)
    runs = _split_ring_into_runs(ring, share_vertices=True)
    assert len(runs) > 1

    def edges(poly):
        keys = set()
        for a, b in zip(poly[:-1], poly[1:]):
            keys.add(frozenset((tuple(np.round(a, 9)), tuple(np.round(b, 9)))))
        return keys

    covered = set().union(*(edges(r) for r in runs))
    assert edges(ring) <= covered
    for run in runs:
        assert run[-1, 2] >= run[0, 2]
    # sans partage, des segments manquent : le test verifie bien quelque chose
    disjoint = set().union(*(edges(r) for r in _split_ring_into_runs(ring)))
    assert not edges(ring) <= disjoint


def test_continuous_rings_closes_original_contour(flat_box_layers):
    """Le contour d'origine n=0 sort OUVERT de la classification : sans
    point de fermeture, son dernier segment n'est jamais trace. Avec
    `continuous_rings`, il est ferme (premier point repete)."""
    mesh_data, overhang_map, layers = flat_box_layers
    base = Step6Params(nozzle_width=0.4, wall_count=1, infill_density=0.2)
    closed = Step6Params(nozzle_width=0.4, wall_count=1, infill_density=0.2, continuous_rings=True)

    def original(result):
        found = []
        for layer in result:
            contour_2d = _project_to_local_2d(layer.contour_3d[0][:-1], layer.local_frame)
            for polyline in layer.toolpath_2d or []:
                if len(polyline) >= len(contour_2d) and np.allclose(polyline[: len(contour_2d), :2], contour_2d):
                    found.append((polyline, len(contour_2d)))
        return found

    open_ones = original(toolpath_run(layers, mesh_data, overhang_map, base))
    closed_ones = original(toolpath_run(layers, mesh_data, overhang_map, closed))
    assert open_ones and closed_ones
    assert all(len(p) == n for p, n in open_ones)
    assert all(len(p) == n + 1 and np.allclose(p[0], p[-1]) for p, n in closed_ones)


def test_monotonic_runs_ignores_oscillation_below_min_amplitude():
    # Bruit de 4 um autour d'une couche plate : une seule course, pas une coupe par oscillation.
    offsets = 0.002 * np.sin(np.arange(200) * 1.7)
    assert len(_monotonic_runs(offsets)) > 50                      # defaut historique : coupe partout
    assert _monotonic_runs(offsets, min_amplitude=0.05) == [(0, 200)]


def test_monotonic_runs_still_cuts_a_real_reversal_with_min_amplitude():
    offsets = np.concatenate([np.linspace(0, 1, 50), np.linspace(1, 0, 50)[1:]])
    assert _monotonic_runs(offsets, min_amplitude=0.05) == [(0, 50), (50, 99)]


# ---------------------------------------------------------------------------
# OSCILLATING_SPIRAL
# ---------------------------------------------------------------------------

def _spiral_params(**overrides):
    base = dict(nozzle_width=0.4, wall_count=2, infill_density=0.2,
                infill_pattern=InfillPattern.OSCILLATING_SPIRAL)
    base.update(overrides)
    return Step6Params(**base)


def _disk(radius=13.0):
    from shapely.geometry import Point
    return Point(0.0, 0.0).buffer(radius, resolution=64)


def _raw_spiral(rho, pitch=3.1, amplitude=0.48, seed=0, rank=0):
    rng = np.random.default_rng([seed, rank])
    return _oscillating_spiral_polyline(
        np.zeros(2), 0.0, 13.0, pitch, amplitude, 1.6, rho, 0.4, rng,
        first_sign=1.0, angle_offset=0.0,
    )


@pytest.mark.parametrize("rho", [0.0, 0.5, 1.0])
def test_oscillating_spiral_never_goes_backward_and_never_crosses_itself(rho):
    """Le plateau ne recule jamais (theta croissant) et la trajectoire ne se
    recoupe pas : deux pics purement radiaux consecutifs repasseraient sur la
    meme ligne (defaut trouve a la mise au point)."""
    from shapely.geometry import LineString
    line = _raw_spiral(rho)
    angles = np.unwrap(np.arctan2(line[:, 1], line[:, 0]))
    assert np.all(np.diff(angles) >= -1e-9)
    assert LineString(line).is_simple


def test_oscillating_spiral_has_radial_peaks_that_do_not_advance_only_when_randomized():
    def stationary(line):
        radii = np.hypot(line[:, 0], line[:, 1])
        angles = np.unwrap(np.arctan2(line[:, 1], line[:, 0]))
        return int(np.sum((np.abs(np.diff(angles)) < 1e-9) & (np.abs(np.diff(radii)) > 1e-3)))
    assert stationary(_raw_spiral(0.0)) == 0
    assert stationary(_raw_spiral(1.0)) > 0


def test_oscillating_spiral_amplitude_is_capped_so_neighbouring_turns_never_touch():
    pitch, width = 3.1, 0.4
    line = _raw_spiral(1.0, pitch=pitch, amplitude=50.0)          # seuil demande absurde
    angles = np.unwrap(np.arctan2(line[:, 1], line[:, 0]))
    deviation = np.hypot(line[:, 0], line[:, 1]) - pitch * angles / (2 * np.pi)
    assert np.abs(deviation).max() <= 0.5 * (pitch - width) + 1e-9


@pytest.mark.parametrize("rho", [0.0, 0.5, 1.0])
def test_oscillating_spiral_infill_hits_the_requested_density_inside_the_zone(rho):
    from shapely.geometry import LineString
    params = _spiral_params(spiral_randomization_threshold=rho)
    chains = _generate_spiral_infill(_disk(), params, 0, np.zeros(2))
    zone = _disk().buffer(-0.2, join_style="mitre")
    length = sum(float(np.linalg.norm(np.diff(c, axis=0), axis=1).sum()) for c in chains)
    assert length * 0.4 / zone.area == pytest.approx(0.2, rel=0.12)
    assert all(zone.buffer(1e-3).covers(LineString(c)) for c in chains)


def test_oscillating_spiral_is_reproducible_and_alternates_between_layers():
    params = _spiral_params()
    a = _generate_spiral_infill(_disk(), params, 0, np.zeros(2))
    b = _generate_spiral_infill(_disk(), params, 0, np.zeros(2))
    assert len(a) == len(b) and all(np.array_equal(x, y) for x, y in zip(a, b))
    other_seed = _generate_spiral_infill(_disk(), _spiral_params(spiral_seed=7), 0, np.zeros(2))
    assert not np.array_equal(a[0], other_seed[0])
    odd = _generate_spiral_infill(_disk(), params, 1, np.zeros(2))
    assert not np.array_equal(a[0][: len(odd[0])], odd[0][: len(a[0])])


def test_oscillating_spiral_layers_share_the_same_median_line():
    """Les couches se superposent sur la meme ligne mediane : l'ecart radial
    de chaque couche a la spirale de base reste borne (plus de decalage d'un
    demi-pas entre couches)."""
    pitch = 2.0
    params = _spiral_params(spiral_randomization_threshold=0.7)
    for rank in (0, 1):
        for chain in _generate_spiral_infill(_disk(), params, rank, np.zeros(2)):
            angles = np.unwrap(np.arctan2(chain[:, 1], chain[:, 0]))
            radius = np.hypot(chain[:, 0], chain[:, 1])
            offset = np.mod(radius - pitch * angles / (2 * np.pi) + pitch / 2, pitch) - pitch / 2
            assert np.abs(offset).max() < 0.5 * pitch


def test_oscillating_spiral_connectors_never_make_the_plate_go_backward():
    from pipeline.step6_toolpath2d import _angle_nondecreasing
    for rank in (0, 1):
        for chain in _generate_spiral_infill(_disk(), _spiral_params(spiral_randomization_threshold=1.0),
                                             rank, np.zeros(2)):
            assert _angle_nondecreasing(chain, np.zeros(2))


def test_plate_axis_in_local_2d_follows_the_layer_plane():
    frame = (np.array([5.0, 2.0, 1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    assert np.allclose(_plate_axis_in_local_2d(frame, (0.0, 0.0, 0.0)), [-5.0, -2.0])
    # plan incline de 45 deg autour de x, passant par l'origine : l'axe vertical
    # x = y = 0 le coupe a l'origine, donc (0, 0) en local
    tilt = (np.zeros(3), np.array([1.0, 0.0, 0.0]), np.array([0.0, np.cos(np.pi / 4), np.sin(np.pi / 4)]))
    assert np.allclose(_plate_axis_in_local_2d(tilt, (0.0, 0.0, 3.0)), [0.0, 0.0], atol=1e-9)
    # un axe decale de 1 mm en x se projette a 1 mm sur u
    assert _plate_axis_in_local_2d(tilt, (1.0, 0.0, 0.0))[0] == pytest.approx(1.0)


def test_oscillating_spiral_through_run_is_monotone_around_the_axis_and_alternates(flat_box_layers):
    mesh_data, overhang_map, _layers = flat_box_layers
    layers = [_flat_square_layer(0, 0.0), _flat_square_layer(1, 0.2)]
    params = _spiral_params(spiral_axis_origin=(0.0, 0.0, 0.0), spiral_randomization_threshold=0.5)
    result = toolpath_run(layers, mesh_data, overhang_map, params)
    infills = []
    for layer in result:
        chains = [p for p in layer.toolpath_2d if not np.allclose(p[0, :2], p[-1, :2])]
        assert chains, "le remplissage en spirale doit exister"
        for chain in chains:
            angles = np.unwrap(np.arctan2(chain[:, 1], chain[:, 0]))
            assert np.all(np.diff(angles) >= -1e-9)
        infills.append(max(chains, key=len))
    n = min(len(infills[0]), len(infills[1]))
    assert not np.allclose(infills[0][:n, :2], infills[1][:n, :2])
