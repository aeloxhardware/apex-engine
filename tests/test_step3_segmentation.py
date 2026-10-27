"""
tests/test_step3_segmentation.py

Tests pour pipeline/step3_segmentation.py, sur des mesh synthetiques generes
a la volee (aucun fichier externe requis).

Execution : pytest tests/test_step3_segmentation.py -v

Couvre en particulier les limites de conception actees dans CLAUDE.md :
  - #1 : distinction contact-plateau vs surplomb reel (region racine du
    graphe de dependances, aucune arete entrante).
  - #4 : compatibilite face/direction = test de surplomb (pas d'alignement
    direct) - regression sur le bug qui bloquait toutes les faces de murs
    verticaux et empechait toute fusion de region, meme sur un cube trivial.
  - #5 : reference globale (up_vector) pour le graphe de dependances, pas
    une moyenne locale par paire de regions - regression sur les cycles
    artificiels produits par un anneau de regions (sphere tranchee).
  - #6 : bascule fibonacci_sphere_directions (5 axes) vs
    rotation_axis_circle_directions (4 axes reels).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
import trimesh

from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from core.mesh_types import RotationAxisSpec
import scipy.sparse as sp

from core.geometry_utils import build_face_adjacency
from pipeline.step3_segmentation import (
    SegmentationError,
    SegmentationParams,
    _compute_platform_contact,
    _face_direction_compatibility,
    run as segmentation_run,
)


# ---------------------------------------------------------------------------
# Fixtures : mesh synthetiques -> (MeshData, OverhangMap)
# ---------------------------------------------------------------------------

def _pipeline_data_from_trimesh(mesh: trimesh.Trimesh, resolution: float = 1.5):
    """Passe par step1_ingest.run() + step2_overhang.run(), pour tester step3
    exactement dans les conditions du pipeline reel (mêmes fixtures que
    tests/test_step2_overhang.py)."""
    import tempfile
    import os

    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, _, _ = ingest_run(path, resolution=resolution)
    finally:
        os.unlink(path)

    overhang_map = overhang_run(mesh_data, OverhangParams(overhang_angle_threshold_deg=45.0))
    return mesh_data, overhang_map


@pytest.fixture
def flat_box_data():
    """Cube pose a plat, base a z=-5 (extents 10 centres a l'origine)."""
    mesh = trimesh.creation.box(extents=[10, 10, 10])
    return _pipeline_data_from_trimesh(mesh)


@pytest.fixture
def l_shape_data():
    """
    Meme forme en L que test_step2_overhang.py : base verticale (z de 0 a 20)
    + bras horizontal en porte-a-faux (dessous du bras = surplomb franc).
    """
    vertical = trimesh.creation.box(extents=[4, 4, 20])
    vertical.apply_translation([0, 0, 10])

    horizontal = trimesh.creation.box(extents=[20, 4, 4])
    horizontal.apply_translation([8, 0, 20])

    combined = vertical.union(horizontal, engine="manifold")
    return _pipeline_data_from_trimesh(combined, resolution=2.0)


@pytest.fixture
def t_shape_data():
    """
    Pied vertical + 2 bras perpendiculaires en porte-a-faux (+X et -Y) -
    approximation simplifiee de res/test_models/l-triple-axis-overhanged-shape.stl
    (CLAUDE.md, cas reel qui a revele l'absence de scission multi-region).
    """
    foot = trimesh.creation.box(extents=[4, 4, 20])
    foot.apply_translation([0, 0, 10])

    arm_x = trimesh.creation.box(extents=[20, 4, 4])
    arm_x.apply_translation([8, 0, 20])

    arm_y = trimesh.creation.box(extents=[4, 20, 4])
    arm_y.apply_translation([0, -8, 20])

    combined = foot.union(arm_x, engine="manifold").union(arm_y, engine="manifold")
    return _pipeline_data_from_trimesh(combined, resolution=2.0)


@pytest.fixture
def hollow_box_data():
    """Boite creuse, base exterieure a z=-5, cavite interne a mi-hauteur."""
    outer = trimesh.creation.box(extents=[10, 10, 10])
    inner = trimesh.creation.box(extents=[4, 4, 4])
    hollow = outer.difference(inner, engine="manifold")
    return _pipeline_data_from_trimesh(hollow, resolution=1.5)


@pytest.fixture
def sphere_data():
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=5.0)
    return _pipeline_data_from_trimesh(mesh)


@pytest.fixture
def pendant_table_data():
    """
    Table sur pied + excroissance PENDANTE suspendue sous la table, loin du
    pied et s'arretant bien au-dessus du plateau (z=[14,30], plateau a z=0).

    Cas rapporte par l'utilisateur (T retourne avec un element qui pend
    d'une barre). Piege pour un critere de HAUTEUR : le pendant est
    geometriquement plus bas que la table (z=[30,34]) et serait donc
    imprime AVANT elle, alors qu'il y est suspendu et ne peut exister
    qu'apres. Seule l'arborescence de croissance depuis le plateau donne le
    bon ordre : le pendant ne rejoint le plateau qu'a travers la table.
    """
    foot = trimesh.creation.box(extents=[4, 4, 30])
    foot.apply_translation([0, 0, 15])

    table = trimesh.creation.box(extents=[30, 30, 4])
    table.apply_translation([0, 0, 32])

    pendant = trimesh.creation.box(extents=[3, 3, 16])
    pendant.apply_translation([12, 12, 22])

    combined = foot.union(table, engine="manifold").union(pendant, engine="manifold")
    return _pipeline_data_from_trimesh(combined, resolution=2.0)


# ---------------------------------------------------------------------------
# Regression #4 - compatibilite face/direction = surplomb, pas alignement
# ---------------------------------------------------------------------------

def test_wall_faces_compatible_with_up_direction():
    """
    Bug initial corrige (limite #4 du CLAUDE.md) : tester l'alignement direct
    normale/direction bloquait a tort TOUTES les faces de murs verticaux
    (normale perpendiculaire a la verticale), empechant toute fusion de
    region meme sur un cube trivial. Le bon critere est un test de surplomb
    (angle par rapport a "vers le bas local"), sous lequel un mur vertical
    n'est jamais en surplomb, quelle que soit la direction candidate.
    """
    # normales : 4 murs verticaux + un plafond (vers le haut) + un plancher (vers le bas)
    wall_normals = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ]
    )
    ceiling_normal = np.array([[0.0, 0.0, 1.0]])
    floor_normal = np.array([[0.0, 0.0, -1.0]])
    normals = np.vstack([wall_normals, ceiling_normal, floor_normal])

    up_direction = np.array([[0.0, 0.0, 1.0]])
    compatibility = _face_direction_compatibility(normals, up_direction, max_angle_deg=45.0)

    assert np.all(compatibility[:4, 0]), (
        "Les 4 faces de murs verticaux devraient etre compatibles avec la "
        "direction verticale (elles ne sont jamais en surplomb selon 'up'), "
        f"obtenu: {compatibility[:4, 0]}"
    )
    assert compatibility[4, 0], "Le plafond (normale vers le haut) ne doit pas etre en surplomb selon 'up'"
    assert not compatibility[5, 0], "Le plancher (normale vers le bas) doit etre en surplomb selon 'up'"


# ---------------------------------------------------------------------------
# Cas 1 - Cube pose a plat : contact plateau -> region racine sans dependance amont
# ---------------------------------------------------------------------------

def test_flat_box_platform_region_has_no_incoming_dependency(flat_box_data):
    mesh_data, overhang_map = flat_box_data
    params = SegmentationParams(platform_z=-5.0, platform_contact_tolerance=0.5)

    result = segmentation_run(mesh_data, overhang_map, params)

    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    is_bottom_face = face_centers[:, 2] < face_centers[:, 2].min() + 1e-3

    # trouve la/les regions contenant les faces du plateau
    platform_regions = set()
    for region in result.regions.values():
        if np.any(is_bottom_face[region.face_ids]):
            platform_regions.add(region.id)

    assert platform_regions, "Au moins une region devrait contenir les faces du plateau"

    incoming_targets = {b for (_, b) in result.dependencies}
    for region_id in platform_regions:
        assert region_id not in incoming_targets, (
            f"La region {region_id} (contact plateau) ne devrait avoir aucune "
            f"dependance amont - le plateau est la reference la plus basse"
        )


def test_flat_box_produces_acyclic_graph(flat_box_data):
    """Un cube trivial doit segmenter sans lever SegmentationError (regression #4)."""
    mesh_data, overhang_map = flat_box_data
    params = SegmentationParams(platform_z=-5.0, platform_contact_tolerance=0.5)

    result = segmentation_run(mesh_data, overhang_map, params)

    assert len(result.regions) >= 1
    # run() leve SegmentationError si le graphe a un cycle - l'atteindre ici
    # confirme deja l'acyclicite, mais on revalide explicitement le contrat
    # de sortie sur les paires de dependances.
    for a, b in result.dependencies:
        assert a in result.regions and b in result.regions


def test_flat_box_platform_contact_detected_automatically_when_platform_z_none(flat_box_data):
    """
    Detection automatique de platform_z (None = derive min(z) du mesh) :
    sur flat_box_data (base a z=-5, PAS z=0), doit produire EXACTEMENT le
    meme resultat qu'un appel explicite platform_z=-5.0 - preuve directe
    que la detection retrouve la vraie hauteur du plateau, pas 0.0. Corrige
    un bug reel mesure sur res/test_models/death.stl (min(z)=1.245) : avec
    l'ancien defaut fige platform_z=0.0, aucune face de ce cube (base a
    z=-5) n'aurait ete detectee comme contact plateau (CLAUDE.md).
    """
    mesh_data, overhang_map = flat_box_data
    params_auto = SegmentationParams(platform_z=None, platform_contact_tolerance=0.5)
    params_explicit = SegmentationParams(platform_z=-5.0, platform_contact_tolerance=0.5)

    result_auto = segmentation_run(mesh_data, overhang_map, params_auto)
    result_explicit = segmentation_run(mesh_data, overhang_map, params_explicit)

    assert np.array_equal(
        result_auto.platform_contact_faces, result_explicit.platform_contact_faces
    )
    assert np.any(result_auto.platform_contact_faces), (
        "la detection automatique doit trouver au moins une face de contact "
        "plateau sur un cube pose a plat, meme quand sa base n'est pas a z=0"
    )


# ---------------------------------------------------------------------------
# Cas 2 - Forme en L : le bras en surplomb doit obtenir une direction alternative
# ---------------------------------------------------------------------------

def test_l_shape_overhang_region_uses_alternate_direction(l_shape_data):
    mesh_data, overhang_map = l_shape_data
    params = SegmentationParams(platform_z=0.0, platform_contact_tolerance=0.5)

    result = segmentation_run(mesh_data, overhang_map, params)

    up = np.array([0.0, 0.0, 1.0])
    overhang_face_ids = np.where(overhang_map.is_overhang)[0]
    assert len(overhang_face_ids) > 0, "Le bras horizontal devrait avoir des faces en surplomb (cf. test_step2)"

    # la region contenant le plus de faces en surplomb ne devrait pas avoir
    # retenu la verticale pure comme direction d'impression : c'est
    # precisement la raison d'etre de la segmentation en sous-volumes.
    region_of_overhang = {}
    for region in result.regions.values():
        overlap = np.intersect1d(region.face_ids, overhang_face_ids)
        if len(overlap) > 0:
            region_of_overhang[region.id] = len(overlap)

    assert region_of_overhang, "Au moins une region devrait contenir des faces en surplomb"
    dominant_region_id = max(region_of_overhang, key=region_of_overhang.get)
    dominant_direction = result.regions[dominant_region_id].target_direction

    angle_from_up = np.degrees(
        np.arccos(np.clip(np.dot(dominant_direction, up), -1.0, 1.0))
    )
    assert angle_from_up > 1e-3, (
        "La region dominee par le surplomb du bras ne devrait pas avoir "
        "retenu la verticale pure comme direction d'impression"
    )


def test_l_shape_produces_acyclic_graph(l_shape_data):
    mesh_data, overhang_map = l_shape_data
    params = SegmentationParams(platform_z=0.0, platform_contact_tolerance=0.5)

    # ne doit pas lever SegmentationError
    result = segmentation_run(mesh_data, overhang_map, params)
    assert len(result.regions) >= 1


# ---------------------------------------------------------------------------
# Cas 3 - Boite creuse : le fond exterieur (plateau) reste racine du graphe
# ---------------------------------------------------------------------------

def test_hollow_box_platform_region_has_no_incoming_dependency(hollow_box_data):
    mesh_data, overhang_map = hollow_box_data
    params = SegmentationParams(platform_z=-5.0, platform_contact_tolerance=0.5)

    result = segmentation_run(mesh_data, overhang_map, params)

    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    is_bottom_face = face_centers[:, 2] < face_centers[:, 2].min() + 1e-3

    platform_regions = {
        region.id
        for region in result.regions.values()
        if np.any(is_bottom_face[region.face_ids])
    }
    assert platform_regions, "Au moins une region devrait contenir le fond exterieur de la boite"

    incoming_targets = {b for (_, b) in result.dependencies}
    for region_id in platform_regions:
        assert region_id not in incoming_targets, (
            f"La region {region_id} (fond exterieur, contact plateau) ne "
            f"devrait avoir aucune dependance amont"
        )


# ---------------------------------------------------------------------------
# Regression #5 - reference globale up_vector, pas de cycle sur un anneau de regions
# ---------------------------------------------------------------------------

def test_sphere_dependency_graph_has_no_cycle(sphere_data):
    """
    Une sphere segmentee en tranches autour d'un axe de rotation forme un
    anneau de regions adjacentes. Une reference locale par paire produirait
    des cycles artificiels (chaque paire semble coherente localement, mais
    l'ensemble ne l'est pas) - c'est exactement le cas que la reference
    globale (up_vector) doit eviter. platform_z est place tres en dessous du
    mesh pour qu'aucune face ne soit consideree en contact plateau ici : on
    isole la question de l'acyclicite du graphe de dependances.
    """
    mesh_data, overhang_map = sphere_data
    params = SegmentationParams(
        platform_z=-100.0,
        platform_contact_tolerance=0.5,
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0, 360))],
        direction_candidate_count=12,
    )

    # ne doit pas lever SegmentationError (cycle)
    result = segmentation_run(mesh_data, overhang_map, params)
    assert len(result.regions) >= 1


# ---------------------------------------------------------------------------
# Arborescence de croissance depuis le plateau : l'ordre suit "par ou la
# matiere rejoint le plateau", jamais la hauteur (cas de la piece pendante).
# ---------------------------------------------------------------------------

def _regions_containing(result, mesh_data, mask):
    return {r.id for r in result.regions.values() if np.any(mask[r.face_ids])}


def test_pendant_is_printed_after_what_it_hangs_from(pendant_table_data):
    """
    Le pendant est plus BAS que la table mais y est SUSPENDU : il doit venir
    APRES elle. Un critere de hauteur donnerait l'ordre inverse - c'est
    exactement le defaut que l'arborescence corrige.
    """
    mesh_data, overhang_map = pendant_table_data
    params = SegmentationParams(
        platform_z=0.0, platform_contact_tolerance=0.5, max_region_merge_angle_deg=30.0
    )

    result = segmentation_run(mesh_data, overhang_map, params)

    fc = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    pendant_mask = (
        (np.abs(fc[:, 0] - 12) < 2.0) & (np.abs(fc[:, 1] - 12) < 2.0) & (fc[:, 2] < 30.5)
    )
    pendant_regions = _regions_containing(result, mesh_data, pendant_mask)
    assert pendant_regions, "La fixture devrait produire au moins une region de pendant"

    depths = result.region_print_depth
    assert depths is not None, "L'arborescence devrait etre exposee"

    # le pendant descend plus bas que la table, mais doit avoir une
    # profondeur d'arborescence STRICTEMENT superieure a celle des regions
    # de la table (il ne rejoint le plateau qu'a travers elle).
    table_mask = fc[:, 2] >= 29.5
    table_regions = _regions_containing(result, mesh_data, table_mask) - pendant_regions
    assert table_regions, "La fixture devrait produire au moins une region de table"

    shallowest_table = min(depths[r] for r in table_regions)
    for pendant_region in pendant_regions:
        assert depths[pendant_region] > shallowest_table, (
            f"Region pendante {pendant_region} (profondeur "
            f"{depths[pendant_region]}) devrait etre plus profonde que la "
            f"table (profondeur minimale {shallowest_table}) - elle y est "
            f"suspendue. Profondeurs completes : {depths}"
        )

    # ... et le pendant ne doit jamais etre declare AVANT une region non-pendante
    for before, after in result.dependencies:
        if before in pendant_regions and after not in pendant_regions:
            pytest.fail(
                f"Region pendante {before} declaree avant {after} qui la "
                f"supporte (dependances : {sorted(result.dependencies)})"
            )


def test_platform_regions_are_roots_of_the_tree(pendant_table_data):
    """Une region au contact du plateau est une racine : profondeur 0, aucune arete entrante."""
    mesh_data, overhang_map = pendant_table_data
    params = SegmentationParams(platform_z=0.0, platform_contact_tolerance=0.5)

    result = segmentation_run(mesh_data, overhang_map, params)

    platform_regions = {
        r.id
        for r in result.regions.values()
        if np.any(result.platform_contact_faces[r.face_ids])
    }
    assert platform_regions, "La fixture repose sur le plateau"

    incoming = {b for (_, b) in result.dependencies}
    for region_id in platform_regions:
        assert result.region_print_depth[region_id] == 0, (
            f"Region plateau {region_id} devrait etre a profondeur 0"
        )
        assert region_id not in incoming, (
            f"Region plateau {region_id} ne devrait avoir aucune dependance amont"
        )


def test_depths_are_defined_for_every_region_even_without_platform_contact(sphere_data):
    """
    Sans aucune racine plateau (contact tangent : les 3 sommets d'une face
    doivent tenir sous le plateau+tolerance, ce qu'une sphere ne fait pas),
    l'arborescence doit quand meme etre definie pour TOUTE region via une
    racine de repli - jamais de region sans profondeur, jamais d'exception.
    """
    mesh_data, overhang_map = sphere_data
    params = SegmentationParams(
        platform_z=-100.0,
        platform_contact_tolerance=0.5,
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0, 360))],
        direction_candidate_count=12,
    )

    result = segmentation_run(mesh_data, overhang_map, params)

    assert not np.any(result.platform_contact_faces), "Aucune face ne doit toucher le plateau ici"
    assert set(result.region_print_depth) == set(result.regions), (
        "Chaque region doit avoir une profondeur, meme sans racine plateau"
    )
    assert min(result.region_print_depth.values()) == 0, "Une racine de repli doit exister"


def test_internal_cavity_regions_are_reachable_in_the_tree(hollow_box_data):
    """
    Une cavite interne forme une nappe de surface qui ne partage aucune
    arete avec la coque exterieure : elle est INATTEIGNABLE depuis les
    racines plateau. Chaque composante connexe doit donc recevoir sa propre
    racine de repli - aucune region ne doit rester sans profondeur.
    """
    mesh_data, overhang_map = hollow_box_data
    params = SegmentationParams(platform_z=-5.0, platform_contact_tolerance=0.5)

    result = segmentation_run(mesh_data, overhang_map, params)

    assert set(result.region_print_depth) == set(result.regions), (
        f"Toutes les regions doivent avoir une profondeur, y compris celles "
        f"de la cavite interne. Manquantes : "
        f"{set(result.regions) - set(result.region_print_depth)}"
    )


def test_sphere_uses_rotation_axis_circle_when_axis_provided(sphere_data):
    """
    Bascule 4 axes vs 5 axes (limite #6) : quand rotation_axes est fourni
    avec un seul axe a plage complete (0,360) - degenerescence a 1 DOF de
    composed_rotation_directions, equivalente a rotation_axis_circle_directions
    - toutes les directions de region retenues doivent etre soit la
    verticale (up_vector, toujours en premiere position candidate), soit
    perpendiculaires a l'axe de rotation (cercle atteignable par un seul DOF
    de rotation) - jamais une direction hors de ce cercle, qui ne serait pas
    atteignable par une machine 4 axes a 1 DOF.
    """
    mesh_data, overhang_map = sphere_data
    rotation_axis = np.array([0.0, 0.0, 1.0])
    params = SegmentationParams(
        platform_z=-100.0,
        rotation_axes=[RotationAxisSpec(axis=rotation_axis, angle_range_deg=(0, 360))],
        direction_candidate_count=12,
        up_vector=rotation_axis,
    )

    result = segmentation_run(mesh_data, overhang_map, params)

    for region in result.regions.values():
        direction = region.target_direction
        is_up = np.allclose(direction, rotation_axis, atol=1e-6)
        is_on_circle = abs(np.dot(direction, rotation_axis)) < 1e-6
        assert is_up or is_on_circle, (
            f"La direction {direction} de la region {region.id} n'est ni la "
            f"verticale ni sur le cercle atteignable par l'axe de rotation "
            f"{rotation_axis} - inatteignable par une machine 4 axes reelle"
        )


# ---------------------------------------------------------------------------
# Test de coherence supplementaire : chaque face est assignee a une region
# ---------------------------------------------------------------------------

def test_every_face_assigned_to_exactly_one_region(flat_box_data):
    mesh_data, overhang_map = flat_box_data
    params = SegmentationParams(platform_z=-5.0, platform_contact_tolerance=0.5)

    result = segmentation_run(mesh_data, overhang_map, params)

    all_face_ids = np.concatenate([region.face_ids for region in result.regions.values()])
    assert len(all_face_ids) == mesh_data.face_count, (
        "Chaque face du mesh devrait etre assignee a exactement une region"
    )
    assert len(np.unique(all_face_ids)) == mesh_data.face_count, (
        "Aucune face ne devrait etre assignee a plusieurs regions"
    )


# ---------------------------------------------------------------------------
# Cas 4 - Jonction en T (2 bras perpendiculaires) : critere de fusion
# pondere par ecart angulaire (max_region_merge_angle_deg, CLAUDE.md)
# ---------------------------------------------------------------------------

def test_t_shape_default_produces_single_compromise_region(t_shape_data):
    """
    Comportement HISTORIQUE (max_region_merge_angle_deg=None, defaut) :
    non-regression explicite - une seule region de compromis, exactement le
    comportement mesure avant le fix (CLAUDE.md).
    """
    mesh_data, overhang_map = t_shape_data
    params = SegmentationParams(platform_z=0.0, platform_contact_tolerance=0.5)

    result = segmentation_run(mesh_data, overhang_map, params)

    assert len(result.regions) == 1, (
        "Comportement historique attendu sans le critere angulaire : une "
        "seule region de compromis, meme pour 2 bras perpendiculaires"
    )


def test_t_shape_angular_gate_splits_overhang_regions(t_shape_data):
    """
    Avec le critere angulaire active (valeur calibree par mesure directe sur
    res/test_models/l-triple-axis-overhanged-shape.stl, CLAUDE.md), les
    faces en surplomb de chaque bras doivent finir dans des regions
    DISTINCTES du reste (murs verticaux, direction neutre) - la jonction ne
    doit plus etre reduite a une seule direction de compromis appliquee aux
    deux bras.
    """
    mesh_data, overhang_map = t_shape_data
    params = SegmentationParams(
        platform_z=0.0, platform_contact_tolerance=0.5, max_region_merge_angle_deg=30.0
    )

    result = segmentation_run(mesh_data, overhang_map, params)

    assert len(result.regions) >= 2, (
        "Le critere angulaire devrait scinder la jonction en plusieurs "
        "regions plutot que fusionner les 2 bras en une seule direction"
    )

    up = np.array([0.0, 0.0, 1.0])
    non_vertical_regions = [
        r for r in result.regions.values()
        if np.degrees(np.arccos(np.clip(np.dot(r.target_direction, up), -1.0, 1.0))) > 1e-3
    ]
    assert non_vertical_regions, (
        "Au moins une region devrait avoir retenu une direction non-verticale "
        "pour un des surplombs de bras"
    )


def test_t_shape_angular_gate_produces_acyclic_graph(t_shape_data):
    """Meme garde-fou que test_l_shape_produces_acyclic_graph, sur le T-shape avec le critere actif."""
    mesh_data, overhang_map = t_shape_data
    params = SegmentationParams(
        platform_z=0.0, platform_contact_tolerance=0.5, max_region_merge_angle_deg=30.0
    )

    # ne doit pas lever SegmentationError
    result = segmentation_run(mesh_data, overhang_map, params)
    assert len(result.regions) >= 1


def test_l_shape_angular_gate_does_not_over_fragment(l_shape_data):
    """
    Non-regression : la meme valeur calibree sur le T-shape, appliquee au
    L-shape a un seul bras, ne doit pas produire une fragmentation excessive
    - les faces qui restent compatibles avec la verticale doivent rester
    regroupees (seul le vrai surplomb doit se distinguer directionnellement).
    """
    mesh_data, overhang_map = l_shape_data
    params = SegmentationParams(
        platform_z=0.0, platform_contact_tolerance=0.5, max_region_merge_angle_deg=30.0
    )

    result = segmentation_run(mesh_data, overhang_map, params)

    up = np.array([0.0, 0.0, 1.0])
    non_vertical_regions = [
        r for r in result.regions.values()
        if np.degrees(np.arccos(np.clip(np.dot(r.target_direction, up), -1.0, 1.0))) > 1e-3
    ]
    vertical_regions = [r for r in result.regions.values() if r not in non_vertical_regions]

    # au plus 1 direction non-verticale distincte attendue pour un seul bras
    # en surplomb (pas de sur-fragmentation directionnelle) ; les regions a
    # direction verticale (murs, dessus/dessous) peuvent rester nombreuses
    # en nombre de regions sans que ce soit un probleme (meme direction,
    # simple granularite de decoupage - voir CLAUDE.md).
    distinct_non_vertical_directions = {
        tuple(np.round(r.target_direction, 3)) for r in non_vertical_regions
    }
    assert len(distinct_non_vertical_directions) <= 1, (
        f"Sur-fragmentation directionnelle inattendue sur un L a un seul "
        f"bras : {len(distinct_non_vertical_directions)} directions "
        f"non-verticales distinctes trouvees"
    )
