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
    field = field_run(tet_mesh, mesh_data, overhang_map, regions, FieldSolverParams())
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


# ---------------------------------------------------------------------------
# Direction de fusion des troncons courts (_merge_short_runs)
#
# Cause racine MESUREE de 425 des 1055 defauts d'ordre sur death.stl (40.3%) :
# la fusion arriere systematique recollait un micro-creux (point isole, dz
# median -0.0141mm) a la FIN d'un troncon montant, le faisant deposer APRES un
# point plus haut que lui. Verifie : 425/425 (100%) de ces defauts avaient leur
# predecesseur immediat plus haut, DANS LA MEME polyligne continue.
# ---------------------------------------------------------------------------
