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

import numpy as np
import pytest
import trimesh

from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from core.mesh_types import RotationAxisSpec
from pipeline.step3_segmentation import (
    SegmentationError,
    SegmentationParams,
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


def test_l_shape_with_deviation_weight_splits_base_from_arm(l_shape_data):
    """
    Avec direction_deviation_weight actif, la base verticale (aucun
    surplomb) et le bras horizontal (surplomb franc) doivent former des
    regions SEPAREES, la base retenant une direction proche de la
    verticale - pas seulement "l'angle moyen retenu est plus petit
    qu'avant" (un critere trop faible : une seule region a 30 deg au lieu
    de 54.8 passerait ce test sans resoudre le vrai probleme diagnostique
    sur out_test.glb : la base doit rester PLATE, pas juste moins inclinee).

    Valeurs calibrees empiriquement (balayage direct, pas devinees) :
    - direction_deviation_weight=0.2 : avec region_count_penalty=0.005,
      produit une scission stable sur tout un balayage 0.10-0.50 (marge
      confortable des deux cotes, voir CLAUDE.md).
    - region_count_penalty=0.005 : mesure precise du seuil de rupture -
      la scission se maintient jusqu'a 0.01, se rompt des 0.02 (retombe a
      1 seule region). 0.005 = moitie du seuil de rupture, marge de
      securite. L'ancien defaut (0.1) aurait TOUJOURS empeche la scission.
    """
    mesh_data, overhang_map = l_shape_data
    params = SegmentationParams(
        platform_z=0.0,
        platform_contact_tolerance=0.5,
        direction_deviation_weight=0.2,
        region_count_penalty=0.005,
    )
    result = segmentation_run(mesh_data, overhang_map, params)

    up = np.array([0.0, 0.0, 1.0])

    # Faces du bras identifiees par le surplomb REEL (dessous du bras en
    # porte-a-faux, z=18, cf. l_shape_data), en excluant les coins de la
    # base au contact du plateau (z=0) : ceux-ci sont geometriquement
    # "is_overhang" (normale vers le bas) mais leur risque est force a 0 par
    # _effective_risk_score (ils reposent sur une surface stable, limite #1
    # CLAUDE.md) - les compter comme "surplomb du bras" serait une erreur de
    # test, pas un defaut de step3 (verifie concretement : 2 des 4 faces
    # is_overhang sont a z=0, pas z=18).
    #
    # Faces de la base identifiees loin du coude (z<15, cf. l_shape_data :
    # bras translate a z=20, coude autour de z=18-20.7 mesure), PAS "tout ce
    # qui n'est pas surplomb" - une face proche du coude peut legitimement
    # appartenir a la region du bras sans etre elle-meme la face de surplomb
    # (verifie : region du bras = 5 faces, z de 18.0 a 20.7, dont seulement
    # 2 sont is_overhang). Un seuil large donnerait un faux negatif de test.
    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    is_platform_contact = face_centers[:, 2] < (params.platform_z + params.platform_contact_tolerance)
    overhang_face_ids = np.where(overhang_map.is_overhang & ~is_platform_contact)[0]
    base_face_ids = np.where(face_centers[:, 2] < 15.0)[0]
    assert len(base_face_ids) > 0, "fixture invalide : aucune face de base identifiee sous z=15"
    assert len(overhang_face_ids) > 0, "fixture invalide : aucune face de surplomb du bras (cf. test_step2)"

    def region_ids_for(face_ids: np.ndarray) -> set[int]:
        ids = set()
        for region in result.regions.values():
            if np.intersect1d(region.face_ids, face_ids).size > 0:
                ids.add(region.id)
        return ids

    base_region_ids = region_ids_for(base_face_ids)
    overhang_region_ids = region_ids_for(overhang_face_ids)

    # 1. Ce sont bien des regions DIFFERENTES - le vrai critere, pas un angle moyen.
    assert base_region_ids.isdisjoint(overhang_region_ids), (
        "La base verticale et le bras en surplomb devraient former des "
        "regions separees, pas une region unique inclinee "
        f"(base={base_region_ids}, bras={overhang_region_ids})"
    )

    # 2. La region de la base a retenu une direction proche de la verticale.
    for region_id in base_region_ids:
        angle = np.degrees(np.arccos(np.clip(
            np.dot(result.regions[region_id].target_direction, up), -1.0, 1.0
        )))
        assert angle < 5.0, (
            f"La region de la base (region {region_id}) devrait rester "
            f"quasi-verticale (<5deg), obtenu {angle:.1f}deg"
        )

    # 3. La region du bras garde une direction significativement inclinee -
    #    elle DOIT s'ecarter de la verticale, sinon le surplomb du bras
    #    redevient non supporte (non-regression de
    #    test_l_shape_overhang_region_uses_alternate_direction ci-dessus).
    for region_id in overhang_region_ids:
        angle = np.degrees(np.arccos(np.clip(
            np.dot(result.regions[region_id].target_direction, up), -1.0, 1.0
        )))
        assert angle > 20.0, (
            f"La region du bras (region {region_id}) doit rester inclinee "
            f"pour supporter son surplomb, obtenu {angle:.1f}deg"
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
