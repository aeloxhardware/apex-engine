"""
tests/test_step4_field_solver.py

Tests pour pipeline/step4_field_solver.py, sur des mesh synthetiques generes
a la volee (aucun fichier externe requis).

Execution : pytest tests/test_step4_field_solver.py -v

Les fixtures enchainent le pipeline reel step1 -> step2 -> step3, pour que
step4 recoive exactement les structures qu'il recevra en production.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.spatial
import trimesh

from core.mesh_types import (
    Constraint,
    ConstraintFalloff,
    ConstraintType,
    MeshData,
    Region,
    RegionDependencyGraph,
    TetMesh,
)
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import (
    FieldSolverParams,
    _assemble_stiffness_matrix,
    _blend_tet_directions,
    _branch_indicator,
    _build_gradient_up_operator,
    _clamp_downward_directions,
    _compute_per_tet_gradient,
    _compute_shape_gradients,
    _descending_volume_fraction,
    _diffuse_tet_directions,
    _enforce_upward_directions,
    _identify_platform_nodes,
    _rotate_directions_out_of_overhang,
    _self_overhang_angles,
    _smoothness_modulation_weights,
    build_static_system,
    run as field_run,
    solve_with_constraints,
)


# ---------------------------------------------------------------------------
# Fixtures : mesh synthetiques -> (mesh_data, tet_mesh, overhang_map, regions)
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
    return mesh_data, tet_mesh, overhang_map, regions


@pytest.fixture
def flat_box_data():
    """Cube pose a plat, base a z=-5 (extents 10 centres a l'origine)."""
    mesh = trimesh.creation.box(extents=[10, 10, 10])
    return _pipeline_data_from_trimesh(mesh, platform_z=-5.0)


@pytest.fixture
def l_shape_data():
    """Forme en L : base verticale + bras horizontal en porte-a-faux."""
    vertical = trimesh.creation.box(extents=[4, 4, 20])
    vertical.apply_translation([0, 0, 10])

    horizontal = trimesh.creation.box(extents=[20, 4, 4])
    horizontal.apply_translation([8, 0, 20])

    combined = vertical.union(horizontal, engine="manifold")
    return _pipeline_data_from_trimesh(combined, resolution=2.0, platform_z=0.0)


# ---------------------------------------------------------------------------
# Cas 1 - Cube plat sans contrainte : champ resolu, croissant depuis le plateau
# ---------------------------------------------------------------------------

def test_flat_box_field_increases_away_from_platform(flat_box_data):
    """
    Sans aucune contrainte utilisateur, le champ doit quand meme etre
    resolu (ancrage plateau) et croitre en s'eloignant du plateau : c'est
    ce qui donne a l'etape 5 l'ordre de balayage des isosurfaces.
    """
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data
    params = FieldSolverParams()

    field = field_run(tet_mesh, mesh_data, overhang_map, regions, [], params)

    assert np.all(np.isfinite(field.values)), "Le champ resolu ne doit contenir ni NaN ni inf"

    z = tet_mesh.tet_vertices[:, 2]
    correlation = np.corrcoef(z, field.values)[0, 1]
    assert correlation > 0.5, (
        f"Le champ devrait croitre avec la hauteur z sur un cube pose a plat "
        f"(ancrage a 0 sur le plateau), correlation obtenue: {correlation:.3f}"
    )

    # les noeuds ancres (plateau, z ~= -5) doivent etre proches de 0
    platform_nodes = z < z.min() + 1e-6
    assert np.allclose(field.values[platform_nodes], 0.0, atol=1e-3), (
        "Les noeuds au contact du plateau devraient etre ancres a ~0"
    )


# ---------------------------------------------------------------------------
# Cas 2 - Forme en L : le gradient s'ecarte de la verticale dans le bras
# ---------------------------------------------------------------------------

def test_l_shape_gradient_deviates_from_vertical(l_shape_data):
    """
    Sur une forme en L, la region du bras en porte-a-faux recoit de step3
    une direction d'impression non verticale ; le champ doit en tenir
    compte, donc son gradient ne peut pas rester purement vertical partout.
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data
    params = FieldSolverParams()

    field = field_run(tet_mesh, mesh_data, overhang_map, regions, [], params)

    assert field.gradient is not None
    norms = np.linalg.norm(field.gradient, axis=1)
    significant = norms > 1e-9
    assert np.any(significant), "Le champ ne doit pas etre uniformement constant"

    up = np.array([0.0, 0.0, 1.0])
    cos_to_vertical = np.abs(field.gradient[significant] @ up) / norms[significant]
    angles_deg = np.degrees(np.arccos(np.clip(cos_to_vertical, 0.0, 1.0)))

    assert np.max(angles_deg) > 10.0, (
        f"Le gradient devrait s'ecarter notablement de la verticale quelque part "
        f"dans la forme en L (bras en surplomb), ecart max obtenu: "
        f"{np.max(angles_deg):.2f} deg"
    )


# ---------------------------------------------------------------------------
# Cas 3 - Point d'attache : tire la valeur du champ localement
# ---------------------------------------------------------------------------

def test_attach_point_pulls_field_value_locally(flat_box_data):
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data
    params = FieldSolverParams()

    baseline = field_run(tet_mesh, mesh_data, overhang_map, regions, [], params)

    attach_position = np.array([0.0, 0.0, 3.0])
    target_value = baseline.values.max() * 5.0 + 10.0
    constraint = Constraint(
        type=ConstraintType.ATTACH_POINT,
        position=attach_position,
        target_field_value=target_value,
        influence_radius=3.0,
        falloff=ConstraintFalloff.GAUSSIAN,
        weight=1e4,
    )

    constrained = field_run(tet_mesh, mesh_data, overhang_map, regions, [constraint], params)

    distances = np.linalg.norm(tet_mesh.tet_vertices - attach_position, axis=1)
    near = distances < 3.0
    assert np.any(near), "La fixture doit avoir des noeuds dans le rayon d'influence"

    assert constrained.values[near].mean() > baseline.values[near].mean(), (
        "Un point d'attache avec une valeur cible elevee devrait tirer le champ "
        "vers le haut dans son rayon d'influence"
    )


# ---------------------------------------------------------------------------
# Cas 4 - Regression singularite : aucune contrainte ne doit pas faire echouer
# ---------------------------------------------------------------------------

def test_no_constraints_does_not_raise(flat_box_data):
    """
    Le systeme Laplacien seul est singulier (le champ constant est toujours
    solution). L'ancrage automatique du plateau doit suffire a le rendre
    resoluble meme sans aucune Constraint[] utilisateur - cas courant du
    slicing automatique sans edition manuelle.
    """
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data
    params = FieldSolverParams()

    for constraints in ([], None):
        field = field_run(tet_mesh, mesh_data, overhang_map, regions, constraints, params)
        assert np.all(np.isfinite(field.values)), (
            f"Le champ doit rester fini sans contrainte utilisateur (constraints={constraints!r})"
        )


# ---------------------------------------------------------------------------
# Cas 5 - Coherence des tailles de sortie
# ---------------------------------------------------------------------------

def test_output_shapes_match_mesh_counts(flat_box_data):
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data
    params = FieldSolverParams()

    field = field_run(tet_mesh, mesh_data, overhang_map, regions, [], params)

    assert field.values.shape == (tet_mesh.node_count,)
    assert field.gradient.shape == (tet_mesh.tet_count, 3)


# ---------------------------------------------------------------------------
# Decision 6 - Re-resolution rapide via le systeme statique mis en cache
# ---------------------------------------------------------------------------

def test_cached_static_system_matches_full_run(flat_box_data):
    """
    build_static_system() + solve_with_constraints() doit produire
    exactement le meme resultat que run(), pour les memes contraintes :
    l'enveloppe run() ne doit rien changer au calcul.
    """
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data
    params = FieldSolverParams()

    constraint = Constraint(
        type=ConstraintType.ATTACH_POINT,
        position=np.array([0.0, 0.0, 2.0]),
        target_field_value=7.0,
        influence_radius=3.0,
        falloff=ConstraintFalloff.LINEAR,
        weight=1e3,
    )

    via_run = field_run(tet_mesh, mesh_data, overhang_map, regions, [constraint], params)

    static_system = build_static_system(tet_mesh, mesh_data, overhang_map, regions, params)
    via_cache = solve_with_constraints(static_system, [constraint], params)

    assert np.allclose(via_run.values, via_cache.values, atol=1e-9)


# ---------------------------------------------------------------------------
# CLAUDE.md - Melange spatial des directions de gradient (_blend_tet_directions)
# ---------------------------------------------------------------------------
#
# Fixtures synthetiques minimalistes (pas le pipeline step1->step3 complet) :
# _blend_tet_directions n'a besoin que de MeshData.vertices/faces/kdtree,
# TetMesh.tet_vertices/tetrahedra, et RegionDependencyGraph - construits ici
# a la main pour isoler precisement le comportement de cette fonction pure,
# sans dependre de la segmentation reelle (step3) ni de la tetraedrisation
# reelle (step1/tetgen).

def _single_tet_mesh(centroid: np.ndarray) -> TetMesh:
    """Un seul tetraedre non degenere dont le centroide vaut approximativement `centroid`."""
    offsets = np.array(
        [[-1.0, -1.0, -1.0], [1.0, -1.0, -1.0], [0.0, 1.0, -1.0], [0.0, 0.0, 1.0]]
    )
    verts = centroid + offsets - offsets.mean(axis=0)
    return TetMesh(tet_vertices=verts, tetrahedra=np.array([[0, 1, 2, 3]]))


def _mesh_data_with_faces_at(face_groups: list[np.ndarray]) -> MeshData:
    """
    MeshData minimal : une face triangulaire (degenere, un point x3 - suffit
    car seul le CENTRE de face compte pour le kdtree/les distances) par
    position donnee dans `face_groups` (une liste de positions (3,) par
    groupe, concatenee). Vertices/faces/normales construits pour rester
    coherents en forme avec le contrat MeshData, jamais utilises autrement
    que via face_ids et le kdtree par _blend_tet_directions.
    """
    positions = np.concatenate(face_groups, axis=0)
    n = len(positions)
    vertices = np.repeat(positions, 3, axis=0)
    faces = np.arange(3 * n).reshape(n, 3)
    face_normals = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    vertex_normals = np.repeat(face_normals, 3, axis=0)
    kdtree = scipy.spatial.cKDTree(positions)
    return MeshData(
        vertices=vertices,
        faces=faces,
        face_normals=face_normals,
        vertex_normals=vertex_normals,
        is_watertight=False,
        kdtree=kdtree,
    )


def _region_graph_two_regions(
    face_counts: tuple[int, int], directions: tuple[np.ndarray, np.ndarray]
) -> RegionDependencyGraph:
    offsets = np.cumsum([0, *face_counts])
    regions = {
        i: Region(
            id=i,
            face_ids=np.arange(offsets[i], offsets[i + 1]),
            target_direction=directions[i] / np.linalg.norm(directions[i]),
        )
        for i in range(len(face_counts))
    }
    platform_contact_faces = np.zeros(offsets[-1], dtype=bool)
    return RegionDependencyGraph(regions=regions, dependencies=[], platform_contact_faces=platform_contact_faces)


def test_blend_disabled_by_default_matches_nearest_neighbor():
    """
    direction_blend_radius=0.0 (defaut) doit produire EXACTEMENT le meme
    resultat que le plus-proche-voisin historique - aucune regression de
    comportement sans opt-in explicite.
    """
    region_a_faces = np.array([[10.0, 0.0, 0.0]])
    region_b_faces = np.array([[-10.0, 0.0, 0.0]])
    mesh_data = _mesh_data_with_faces_at([region_a_faces, region_b_faces])
    graph = _region_graph_two_regions(
        (1, 1), (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    )

    tet_mesh = _single_tet_mesh(np.array([10.0, 0.0, 0.0]))  # proche de region_a
    params = FieldSolverParams(direction_blend_radius=0.0)

    result = _blend_tet_directions(tet_mesh, mesh_data, graph, params.up_vector, params)

    assert np.allclose(result[0], [1.0, 0.0, 0.0], atol=1e-6), (
        "Sans melange, un tetraedre proche de region_a doit recevoir EXACTEMENT "
        "sa target_direction (plus-proche-voisin strict, comportement historique)"
    )


def test_blend_near_platform_far_from_region_favors_up_vector():
    """Loin de toute region, proche du plateau : la direction resultante doit etre proche de up_vector."""
    region_faces = np.array([[0.0, 0.0, 100.0]])  # tres loin du tetraedre teste
    # une face de "plateau" a z=0 (necessaire : _compute_platform_level ne
    # traite le plateau comme source que si au moins une face est marquee
    # platform_contact_faces=True, voir test dedie de ce comportement
    # ci-dessous - sans elle, seule la region compte, peu importe la distance).
    platform_faces = np.array([[0.0, 0.0, 0.0]])
    mesh_data = _mesh_data_with_faces_at([platform_faces, region_faces])

    regions = {1: Region(id=1, face_ids=np.array([1]), target_direction=np.array([1.0, 0.0, 0.0]))}
    platform_contact_faces = np.array([True, False])
    graph = RegionDependencyGraph(regions=regions, dependencies=[], platform_contact_faces=platform_contact_faces)

    tet_mesh = _single_tet_mesh(np.array([0.0, 0.0, 0.05]))  # tres proche de z=0 (plateau)
    params = FieldSolverParams(direction_blend_radius=200.0, direction_blend_power=2.0)

    result = _blend_tet_directions(tet_mesh, mesh_data, graph, params.up_vector, params)

    angle_to_up = np.degrees(np.arccos(np.clip(np.dot(result[0], [0.0, 0.0, 1.0]), -1.0, 1.0)))
    assert angle_to_up < 20.0, (
        f"Un tetraedre proche du plateau et loin de toute region devrait recevoir "
        f"une direction proche de up_vector, ecart obtenu: {angle_to_up:.1f} deg"
    )


def test_blend_platform_not_a_source_without_contact_faces():
    """
    Sans aucune face marquee platform_contact_faces=True (segmentation qui
    n'identifie aucun contact plateau), le plateau ne doit PAS etre traite
    comme source d'influence - seule(s) la/les region(s) comptent, quelle
    que soit la distance au plan z=0. Documente explicitement ce
    comportement (source du bug corrige dans le premier jet de
    test_blend_near_platform_far_from_region_favors_up_vector ci-dessus).
    """
    region_faces = np.array([[0.0, 0.0, 100.0]])
    mesh_data = _mesh_data_with_faces_at([region_faces])
    regions = {0: Region(id=0, face_ids=np.array([0]), target_direction=np.array([1.0, 0.0, 0.0]))}
    platform_contact_faces = np.zeros(1, dtype=bool)  # aucun contact plateau declare
    graph = RegionDependencyGraph(regions=regions, dependencies=[], platform_contact_faces=platform_contact_faces)

    tet_mesh = _single_tet_mesh(np.array([0.0, 0.0, 0.05]))
    params = FieldSolverParams(direction_blend_radius=200.0, direction_blend_power=2.0)

    result = _blend_tet_directions(tet_mesh, mesh_data, graph, params.up_vector, params)

    assert np.allclose(result[0], [1.0, 0.0, 0.0], atol=1e-6), (
        "Sans face de contact plateau declaree, la seule region presente doit "
        "l'emporter integralement (le plateau n'existe pas comme source)"
    )


def test_blend_far_from_platform_near_single_region_favors_region_direction():
    """Loin du plateau, proche d'une seule region : quasi winner-take-all en regime etabli."""
    region_faces = np.array([[0.0, 0.0, 10.0]])
    mesh_data = _mesh_data_with_faces_at([region_faces])
    target = np.array([0.6, 0.0, 0.8])
    target = target / np.linalg.norm(target)
    regions = {0: Region(id=0, face_ids=np.array([0]), target_direction=target)}
    platform_contact_faces = np.zeros(1, dtype=bool)
    graph = RegionDependencyGraph(regions=regions, dependencies=[], platform_contact_faces=platform_contact_faces)

    # platform_contact_faces tout a False -> _compute_platform_level retourne
    # None -> le plateau n'est PAS une source d'influence ici (aucune face
    # de contact identifiable) : seule la region compte.
    tet_mesh = _single_tet_mesh(np.array([0.0, 0.0, 10.0]))  # exactement a la region
    params = FieldSolverParams(direction_blend_radius=5.0, direction_blend_power=2.0)

    result = _blend_tet_directions(tet_mesh, mesh_data, graph, params.up_vector, params)

    assert np.allclose(result[0], target, atol=1e-3), (
        "Sans source plateau identifiable et une seule region dans le rayon, "
        "le resultat doit converger vers sa target_direction"
    )


def test_blend_transition_is_monotonic_from_platform_to_region():
    """
    Sur une serie de tetraedres alignes verticalement du plateau vers une
    region, l'angle de la direction resultante par rapport a up_vector doit
    croitre de facon MONOTONE (jamais de saut brutal ni d'oscillation).
    """
    target = np.array([1.0, 0.0, 1.0])
    target = target / np.linalg.norm(target)
    region_faces = np.array([[0.0, 0.0, 20.0]])
    # une face de "plateau" a z=0, coherente avec _compute_platform_level
    # (sommet le plus haut des faces marquees platform_contact_faces)
    mesh_data = _mesh_data_with_faces_at([np.array([[0.0, 0.0, 0.0]]), region_faces])
    platform_contact_faces = np.array([True, False])
    regions = {1: Region(id=1, face_ids=np.array([1]), target_direction=target)}
    graph = RegionDependencyGraph(regions=regions, dependencies=[], platform_contact_faces=platform_contact_faces)

    params = FieldSolverParams(direction_blend_radius=25.0, direction_blend_power=2.0)

    heights = np.linspace(0.5, 19.5, 15)
    angles = []
    for h in heights:
        tet_mesh = _single_tet_mesh(np.array([0.0, 0.0, h]))
        result = _blend_tet_directions(tet_mesh, mesh_data, graph, params.up_vector, params)
        angle = np.degrees(np.arccos(np.clip(np.dot(result[0], [0.0, 0.0, 1.0]), -1.0, 1.0)))
        angles.append(angle)

    diffs = np.diff(angles)
    assert np.all(diffs >= -1e-6), (
        f"L'angle a la verticale devrait croitre de facon monotone du plateau vers la "
        f"region, angles obtenus: {[round(a, 2) for a in angles]}"
    )
    assert angles[-1] > angles[0] + 5.0, "Une vraie transition doit etre mesurable (pas un plat constant)"


def test_blend_opposite_directions_degenerate_case_falls_back_to_one_source():
    """
    Deux sources a directions STRICTEMENT OPPOSEES, poids egaux (tetraedre
    equidistant) : la moyenne ponderee s'annule (norme quasi nulle) - le
    repli doit donner une direction UNITAIRE correspondant a l'une des deux
    sources, jamais un vecteur instable/NaN.
    """
    face_a = np.array([[10.0, 0.0, 0.0]])
    face_b = np.array([[-10.0, 0.0, 0.0]])
    mesh_data = _mesh_data_with_faces_at([face_a, face_b])
    dir_a = np.array([1.0, 0.0, 0.0])
    dir_b = np.array([-1.0, 0.0, 0.0])
    graph = _region_graph_two_regions((1, 1), (dir_a, dir_b))

    tet_mesh = _single_tet_mesh(np.array([0.0, 0.0, 0.0]))  # equidistant des 2 sources
    params = FieldSolverParams(
        direction_blend_radius=50.0, direction_blend_power=2.0, direction_blend_min_norm=0.1
    )

    result = _blend_tet_directions(tet_mesh, mesh_data, graph, params.up_vector, params)

    assert np.all(np.isfinite(result[0])), "Le cas degenere ne doit jamais produire NaN/inf"
    norm = np.linalg.norm(result[0])
    assert abs(norm - 1.0) < 1e-6, f"La direction de repli doit rester unitaire, norme obtenue: {norm:.4f}"
    assert (
        np.allclose(result[0], dir_a, atol=1e-6) or np.allclose(result[0], dir_b, atol=1e-6)
    ), f"Le repli doit correspondre a l'une des 2 sources, obtenu: {result[0]}"


def test_blend_near_opposite_but_not_degenerate_does_not_trigger_fallback():
    """
    Deux sources presque opposees mais avec des poids clairement DIFFERENTS
    (distances tres inegales) : la moyenne ne doit PAS s'annuler - le
    resultat doit rester proche de la source dominante, sans repli force.
    """
    face_near = np.array([[1.0, 0.0, 0.0]])   # tres proche
    face_far = np.array([[-100.0, 0.0, 0.0]])  # tres loin
    mesh_data = _mesh_data_with_faces_at([face_near, face_far])
    dir_near = np.array([1.0, 0.0, 0.0])
    dir_far = np.array([-1.0, 0.0, 0.0])
    graph = _region_graph_two_regions((1, 1), (dir_near, dir_far))

    tet_mesh = _single_tet_mesh(np.array([1.0, 0.0, 0.0]))
    params = FieldSolverParams(
        direction_blend_radius=200.0, direction_blend_power=2.0, direction_blend_min_norm=0.1
    )

    result = _blend_tet_directions(tet_mesh, mesh_data, graph, params.up_vector, params)

    assert np.dot(result[0], dir_near) > 0.9, (
        f"La source dominante (tres proche) devrait l'emporter largement, obtenu: {result[0]}"
    )


def test_same_static_system_reused_with_different_constraints(flat_box_data):
    """
    Le meme systeme statique doit pouvoir servir a plusieurs resolutions
    successives avec des contraintes differentes (edition d'un point
    d'attache), sans etre reassemble et sans etat cache d'une resolution
    a l'autre.
    """
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data
    params = FieldSolverParams()

    static_system = build_static_system(tet_mesh, mesh_data, overhang_map, regions, params)

    low = Constraint(
        type=ConstraintType.ATTACH_POINT,
        position=np.array([0.0, 0.0, 2.0]),
        target_field_value=1.0,
        influence_radius=3.0,
        weight=1e3,
    )
    high = Constraint(
        type=ConstraintType.ATTACH_POINT,
        position=np.array([0.0, 0.0, 2.0]),
        target_field_value=50.0,
        influence_radius=3.0,
        weight=1e3,
    )

    field_low = solve_with_constraints(static_system, [low], params)
    field_high = solve_with_constraints(static_system, [high], params)

    assert not np.allclose(field_low.values, field_high.values), (
        "Deux jeux de contraintes differents doivent produire des champs differents"
    )

    # re-resoudre avec les contraintes initiales doit redonner le meme
    # resultat : aucune contamination d'une resolution a l'autre
    field_low_again = solve_with_constraints(static_system, [low], params)
    assert np.allclose(field_low.values, field_low_again.values, atol=1e-9), (
        "Le systeme statique ne doit pas etre mute par une resolution"
    )


# ---------------------------------------------------------------------------
# Redressement des directions de croissance descendantes (Decision 3bis)
#
# Contexte mesure (death.stl) : 67.6% du volume avait grad(phi).up < 0, ce qui
# empile les couches du haut vers le bas - chaque couche atterrit dans le vide.
# Les fixtures ci-dessous reproduisent la cause en miniature synthetique
# (une region NON-RACINE a direction cible descendante), jamais death.stl.
# ---------------------------------------------------------------------------

def _descending_region_graph(
    regions: RegionDependencyGraph, direction=(0.0, 0.0, -1.0)
) -> RegionDependencyGraph:
    """
    Force la direction cible de TOUTES les regions a pointer vers le bas.

    Reproduit en miniature synthetique la cause mesuree sur death.stl : une
    direction de croissance a composante Z negative, qui empile les couches
    du haut vers le bas. Le forcage est inconditionnel (et non restreint aux
    regions hors plateau) parce que les fixtures synthetiques simples de ce
    fichier ne produisent qu'UNE SEULE region, elle-meme au contact du
    plateau : un forcage conditionnel n'y changerait rien et le test
    passerait sans jamais exercer le mecanisme.
    """
    forced = {
        region_id: Region(
            id=region.id,
            face_ids=region.face_ids,
            target_direction=np.array(direction, dtype=float),
        )
        for region_id, region in regions.regions.items()
    }
    return RegionDependencyGraph(
        regions=forced,
        dependencies=list(regions.dependencies),
        platform_contact_faces=regions.platform_contact_faces,
        region_print_depth=regions.region_print_depth,
    )


def _descending_volume_of(field, tet_mesh) -> float:
    """
    Fraction de volume dont le gradient resolu descend - mesure INDEPENDANTE
    du diagnostic produit par step4 (ne le prend jamais pour argent comptant).
    """
    vertices = tet_mesh.tet_vertices[tet_mesh.tetrahedra]
    volumes = np.abs(
        np.einsum(
            "ij,ij->i",
            vertices[:, 1] - vertices[:, 0],
            np.cross(vertices[:, 2] - vertices[:, 0], vertices[:, 3] - vertices[:, 0]),
        )
    ) / 6.0
    descending = (np.asarray(field.gradient) @ np.array([0.0, 0.0, 1.0])) < -1e-9
    total = float(volumes.sum())
    return float(volumes[descending].sum() / total) if total > 0 else 0.0


def test_upward_growth_operator_matches_per_tet_gradient():
    """
    Test unitaire PUR de l'operateur creux G_up, decouple du solveur : sur un
    tetraedre de coordonnees connues, (G_up @ phi) doit egaler exactement la
    projection sur `up` du gradient calcule par la voie deja en place
    (_compute_per_tet_gradient). Verrouille la construction sparse elle-meme.
    """
    tet_mesh = TetMesh(
        tet_vertices=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        ),
        tetrahedra=np.array([[0, 1, 2, 3]]),
    )
    shape_gradients, _volumes = _compute_shape_gradients(tet_mesh)
    values = np.array([0.0, 2.0, -1.0, 5.0])

    for up in (np.array([0.0, 0.0, 1.0]), np.array([1.0, 2.0, -3.0])):
        operator = _build_gradient_up_operator(tet_mesh, shape_gradients, up)
        normalized_up = up / np.linalg.norm(up)
        expected = (
            _compute_per_tet_gradient(tet_mesh, shape_gradients, values) @ normalized_up
        )
        assert np.allclose(operator @ values, expected)


def test_enforce_upward_directions_flips_only_descending():
    """
    Seules les directions descendantes sont retournees ; une direction
    perpendiculaire a `up` est laissee telle quelle (elle ne descend pas, et
    la retourner rendrait le resultat sensible au bruit numerique).
    """
    directions = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.5, -0.5],
        ]
    )
    result, flipped = _enforce_upward_directions(
        directions, np.array([0.0, 0.0, 1.0]), 1e-9
    )

    assert list(flipped) == [False, True, False, True]
    assert np.allclose(result[0], [0.0, 0.0, 1.0])
    assert np.allclose(result[1], [0.0, 0.0, 1.0])
    assert np.allclose(result[2], [1.0, 0.0, 0.0])
    assert np.allclose(result[3], [-0.5, -0.5, 0.5])
    assert np.allclose(directions[1], [0.0, 0.0, -1.0])
    assert np.all(result @ np.array([0.0, 0.0, 1.0]) >= -1e-9)


def test_upward_growth_disabled_by_default(l_shape_data):
    """
    Non-regression stricte : le defaut (`enforce_upward_growth=False`) ne
    produit aucun diagnostic - comportement historique inchange.
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data
    params = FieldSolverParams()
    assert params.enforce_upward_growth is False

    field = field_run(tet_mesh, mesh_data, overhang_map, regions, [], params)
    assert field.upward_growth_diagnostic is None


def test_upward_growth_collapses_descending_volume(l_shape_data):
    """
    Le coeur du mecanisme : sur une fixture dont les regions non-racines ont
    une direction cible descendante, activer le redressement doit effondrer
    la fraction de volume dont le gradient descend.
    """
    mesh_data, tet_mesh, overhang_map, base_regions = l_shape_data
    regions = _descending_region_graph(base_regions)

    before = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(enforce_upward_growth=False),
    )
    after = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(enforce_upward_growth=True),
    )

    descending_before = _descending_volume_of(before, tet_mesh)
    descending_after = _descending_volume_of(after, tet_mesh)

    assert descending_before > 0.2, (
        f"fixture invalide : seulement {descending_before:.1%} de volume "
        "descendant avant redressement, le test ne discriminerait rien"
    )
    assert descending_after < descending_before / 2.0
    assert after.upward_growth_diagnostic is not None
    assert after.upward_growth_diagnostic.flipped_tet_count > 0


def test_upward_growth_diagnostic_matches_independent_measure(l_shape_data):
    """
    Le residu rapporte par step4 doit correspondre a une mesure recalculee
    independamment depuis le champ retourne - un diagnostic qui derive
    silencieusement de la realite serait pire que pas de diagnostic.
    """
    mesh_data, tet_mesh, overhang_map, base_regions = l_shape_data
    regions = _descending_region_graph(base_regions)

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(enforce_upward_growth=True),
    )
    diagnostic = field.upward_growth_diagnostic
    assert diagnostic is not None
    assert diagnostic.residual_descending_volume_fraction == pytest.approx(
        _descending_volume_of(field, tet_mesh), abs=5e-3
    )


def test_upward_growth_noop_when_all_directions_already_up(flat_box_data):
    """
    Sur un cas deja entierement montant (cube pose a plat), le redressement ne
    retourne AUCUN tetraedre : le mecanisme ne doit jamais s'activer la ou il
    n'y a rien a corriger.
    """
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(enforce_upward_growth=True),
    )
    assert field.upward_growth_diagnostic is not None
    assert field.upward_growth_diagnostic.flipped_tet_count == 0


def test_descending_volume_fraction_weights_by_volume():
    """
    La fraction est ponderee par VOLUME, pas un compte brut de tetraedres :
    un maillage tetraedrique n'est pas uniforme, un compte brut surponderait
    les zones finement maillees.
    """
    grad_dot_up = np.array([-1.0, 1.0, 1.0, 1.0])
    volumes = np.array([9.0, 1.0, 1.0, 1.0])

    fraction, count = _descending_volume_fraction(grad_dot_up, volumes, 1e-9)

    assert count == 1
    assert fraction == pytest.approx(0.75)


# ---------------------------------------------------------------------------
# Boucle de retroaction de surplomb (support_feedback_iterations)
#
# Ferme le circuit ouvert du pipeline historique : step2 evalue le surplomb
# contre la verticale MACHINE avant que le champ n'existe, step9b le mesure
# apres tout sans rien reinjecter. La boucle re-evalue contre `grad(phi)`.
# ---------------------------------------------------------------------------

def test_self_overhang_angles_measure_against_the_field_not_the_machine():
    """
    Le critere DOIT etre evalue contre `-grad(phi)` (la verticale locale du
    champ), jamais contre la verticale machine - c'est toute la raison d'etre
    du mecanisme.

    Meme normale (pointant vers -Z, un plafond) dans deux champs differents :
    surplomb franc quand le depot monte (+Z), aucun surplomb quand le depot
    lui-meme descend (-Z, la matiere est alors posee DANS le sens du plafond).
    """
    normals = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, -1.0]])
    gradients = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])

    angles = np.degrees(_self_overhang_angles(gradients, normals))

    assert angles[0] == pytest.approx(0.0, abs=1e-6), "plafond vs depot montant = surplomb franc"
    assert angles[1] == pytest.approx(180.0, abs=1e-6), "depot descendant : plus un surplomb"


def test_self_overhang_angle_is_pi_where_the_field_is_flat():
    """
    Un gradient quasi nul n'a pas de direction de depot definie : le
    tetraedre doit recevoir pi (jamais un surplomb), pas un angle calcule sur
    du bruit numerique.
    """
    normals = np.array([[0.0, 0.0, -1.0]])
    gradients = np.array([[0.0, 0.0, 1e-18]])

    angles = _self_overhang_angles(gradients, normals)

    assert angles[0] == pytest.approx(np.pi)


def test_rotation_corrects_exactly_to_the_threshold_and_no_further():
    """
    La correction doit etre MINIMALE : pivoter au-dela du seuil eloignerait
    la cible de ce que la geometrie de la region demande, sans rien resoudre
    de plus. On verifie que l'angle apres correction vaut exactement le
    seuil, pas davantage.
    """
    normals = np.array([[0.0, 0.0, -1.0]])
    directions = np.array([[0.0, 0.0, 1.0]])   # surplomb franc (angle 0)
    angles = _self_overhang_angles(directions, normals)
    threshold = np.radians(45.0)

    corrected = _rotate_directions_out_of_overhang(
        directions, normals, angles, np.array([True]), threshold
    )

    new_angle = _self_overhang_angles(corrected, normals)[0]
    assert new_angle == pytest.approx(threshold, abs=1e-6)
    assert np.linalg.norm(corrected[0]) == pytest.approx(1.0, abs=1e-9)


def test_rotation_leaves_unmasked_directions_untouched():
    """Un tetraedre non fautif ne doit jamais voir sa direction modifiee."""
    normals = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, -1.0]])
    directions = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    angles = _self_overhang_angles(directions, normals)

    corrected = _rotate_directions_out_of_overhang(
        directions, normals, angles, np.array([True, False]), np.radians(45.0)
    )

    np.testing.assert_allclose(corrected[1], directions[1])
    assert not np.allclose(corrected[0], directions[0])


def test_support_feedback_disabled_by_default(l_shape_data):
    """
    Opt-in strict : sans activation explicite, aucun diagnostic n'est produit
    et le comportement historique est inchange. `None` (mecanisme desactive)
    se distingue d'un diagnostic present a residu nul.
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data

    field = field_run(tet_mesh, mesh_data, overhang_map, regions, [], FieldSolverParams())

    assert field.support_feedback_diagnostic is None


def test_support_feedback_reduces_self_inflicted_overhang(l_shape_data):
    """
    Sur un vrai porte-a-faux, la boucle doit REDUIRE le volume en surplomb
    auto-infligee entre sa premiere et sa derniere iteration.

    Mesure sur le champ REELLEMENT resolu, jamais sur les cibles corrigees :
    rien ne garantit mathematiquement que le lissage FEM respecte une cible
    pivotee (meme discipline que UpwardGrowthDiagnostic).
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(support_feedback_iterations=3),
    )

    diag = field.support_feedback_diagnostic
    assert diag is not None
    assert diag.iterations_used == 3 or diag.converged
    assert diag.corrected_tet_count > 0, "un porte-a-faux doit declencher des corrections"
    assert (
        diag.residual_self_overhang_volume_fraction
        <= diag.initial_self_overhang_volume_fraction
    )
    assert np.all(np.isfinite(field.values)), "la boucle ne doit jamais casser le systeme"


def test_support_feedback_converges_and_stops_early_on_a_flat_box(flat_box_data):
    """
    Un cube pose a plat n'a aucun surplomb : la boucle doit s'arreter des la
    premiere iteration en signalant la convergence, sans consommer les
    suivantes (chaque iteration coute une resolution complete).
    """
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(support_feedback_iterations=5),
    )

    diag = field.support_feedback_diagnostic
    assert diag is not None
    if diag.initial_self_overhang_tet_count == 0:
        assert diag.converged
        assert diag.iterations_used == 1
        assert diag.corrected_tet_count == 0


def test_support_feedback_iteration_cap_is_hard(l_shape_data):
    """
    Le plafond de 50 iterations est DUR : une valeur demandee superieure est
    bornee, jamais honoree. Chaque iteration coutant une resolution complete,
    un parametrage distrait rendrait sinon le pipeline inutilisable sans
    jamais lever d'erreur.
    """
    from pipeline.step4_field_solver import SUPPORT_FEEDBACK_ITERATION_CAP

    mesh_data, tet_mesh, overhang_map, regions = l_shape_data

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(
            support_feedback_iterations=10_000,
            support_feedback_min_improvement=0.0,  # neutralise l'arret dynamique
        ),
    )

    diag = field.support_feedback_diagnostic
    assert diag is not None
    assert diag.iterations_used <= SUPPORT_FEEDBACK_ITERATION_CAP


def test_support_feedback_stops_early_when_progress_stalls(l_shape_data):
    """
    Arret dynamique : avec une exigence de gain irrealiste (100% par
    iteration), la boucle doit s'arreter des la deuxieme iteration et se
    declarer `stalled` - jamais `converged`, qui est reserve a un surplomb
    reellement resorbe.
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(
            support_feedback_iterations=50,
            support_feedback_min_improvement=1.0,
        ),
    )

    diag = field.support_feedback_diagnostic
    assert diag is not None
    if not diag.converged:
        assert diag.stalled
        assert diag.iterations_used <= 2
    assert len(diag.volume_fraction_history) == diag.iterations_used


def test_support_feedback_history_tracks_every_iteration(l_shape_data):
    """
    L'historique porte une mesure PAR iteration, dans l'ordre : c'est ce qui
    permet de lire la courbe de decroissance (et un eventuel rebond) plutot
    que ses seules extremites.
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(support_feedback_iterations=4),
    )

    diag = field.support_feedback_diagnostic
    assert diag is not None
    assert len(diag.volume_fraction_history) == diag.iterations_used
    assert diag.volume_fraction_history[0] == pytest.approx(
        diag.initial_self_overhang_volume_fraction
    )
    assert diag.volume_fraction_history[-1] == pytest.approx(
        diag.residual_self_overhang_volume_fraction
    )


# ---------------------------------------------------------------------------
# Modulation spatiale du lissage (smoothness_modulation_strength)
# ---------------------------------------------------------------------------

def _straight_bar_data():
    """
    Barreau droit, SANS ramification : l'indicateur doit y rester plat.

    Section volontairement EPAISSE (10mm) et maillage assez fin : sur un
    barreau de 4mm a resolution 2.0, tous les tetraedres touchent la frontiere
    et il ne reste aucun coeur a mesurer - mesure faite, la comparaison
    portait alors sur un ensemble vide.
    """
    mesh = trimesh.creation.box(extents=[10, 10, 30])
    mesh.apply_translation([0, 0, 15])
    return _pipeline_data_from_trimesh(mesh, resolution=1.5, platform_z=0.0)


def _thick_l_data():
    """Meme forme en L que `l_shape_data`, epaissie pour avoir un coeur."""
    vertical = trimesh.creation.box(extents=[10, 10, 24])
    vertical.apply_translation([0, 0, 12])

    horizontal = trimesh.creation.box(extents=[24, 10, 10])
    horizontal.apply_translation([12, 0, 19])

    combined = vertical.union(horizontal, engine="manifold")
    return _pipeline_data_from_trimesh(combined, resolution=1.5, platform_z=0.0)


def _indicator_for(data, params=None):
    mesh_data, tet_mesh, _overhang_map, regions = data
    params = params or FieldSolverParams()
    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)
    platform_nodes = _identify_platform_nodes(tet_mesh, mesh_data, regions)
    return _branch_indicator(
        tet_mesh, shape_gradients, volumes, platform_nodes, params
    )


def test_stiffness_matrix_without_weights_keeps_the_exact_sparsity_pattern():
    """
    Le MOTIF de sparsite ne doit pas dependre des poids : c'est lui qui
    garantit qu'une factorisation directe (`spsolve`) coute exactement la meme
    chose avec et sans modulation. Comparaison stricte des indices, jamais une
    tolerance.
    """
    mesh = trimesh.creation.box(extents=[6, 6, 6])
    _mesh_data, tet_mesh, _om, _rg = _pipeline_data_from_trimesh(
        mesh, resolution=2.0, platform_z=-3.0
    )
    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)

    without = _assemble_stiffness_matrix(tet_mesh, shape_gradients, volumes)
    explicit_ones = _assemble_stiffness_matrix(
        tet_mesh, shape_gradients, volumes, np.ones(tet_mesh.tet_count)
    )

    assert without.shape == explicit_ones.shape
    assert np.array_equal(without.indptr, explicit_ones.indptr)
    assert np.array_equal(without.indices, explicit_ones.indices)
    assert without.data == pytest.approx(explicit_ones.data, rel=1e-12)


def test_modulation_disabled_by_default(l_shape_data):
    """
    Defaut = 0.0 = mecanisme absent. Le diagnostic vaut None, ce qui le rend
    distinguable d'un mecanisme actif dont tous les poids vaudraient 1.0
    (cas normal d'une piece sans ramification).
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [], FieldSolverParams()
    )

    assert FieldSolverParams().smoothness_modulation_strength == 0.0
    assert field.smoothness_modulation_diagnostic is None


def test_modulation_reports_a_diagnostic_when_active(l_shape_data):
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(smoothness_modulation_strength=10.0),
    )

    diag = field.smoothness_modulation_diagnostic
    assert diag is not None
    assert len(diag.indicator_quantiles) == len(diag.indicator_quantile_levels)
    assert diag.min_weight <= diag.median_weight <= 1.0
    assert diag.boundary_tet_count > 0        # tout maillage a une frontiere


def test_modulation_weights_stay_within_floor_and_one(l_shape_data):
    """
    L'indicateur |div(X)| n'est PAS borne (mesure jusqu'a 1.8 en mediane de
    tranche sur real_overhang). Le plancher est donc ce qui borne le contraste
    de raideur, donc le conditionnement du systeme.
    """
    mesh_data, tet_mesh, _overhang_map, regions = l_shape_data
    params = FieldSolverParams(
        smoothness_modulation_strength=1000.0, smoothness_modulation_floor=0.05
    )
    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)
    platform_nodes = _identify_platform_nodes(tet_mesh, mesh_data, regions)

    weights, diagnostic = _smoothness_modulation_weights(
        tet_mesh, shape_gradients, volumes, platform_nodes, params
    )

    assert weights.shape == (tet_mesh.tet_count,)
    assert weights.min() >= params.smoothness_modulation_floor
    assert weights.max() <= 1.0
    # A cette force, le plancher DOIT mordre : sinon le test ne verrouille rien.
    assert diagnostic.floored_tet_count > 0


def test_modulated_stiffness_stays_symmetric_positive_semidefinite(l_shape_data):
    """
    Un poids par element ne doit pas casser les proprietes dont depend tout le
    solveur : `cg` exige un systeme symetrique defini positif, et le noyau du
    Laplacien doit rester les seuls champs constants.
    """
    mesh_data, tet_mesh, _overhang_map, regions = l_shape_data
    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)
    platform_nodes = _identify_platform_nodes(tet_mesh, mesh_data, regions)

    weights, _ = _smoothness_modulation_weights(
        tet_mesh, shape_gradients, volumes, platform_nodes,
        FieldSolverParams(smoothness_modulation_strength=25.0),
    )
    K = _assemble_stiffness_matrix(tet_mesh, shape_gradients, volumes, weights)

    assert abs(K - K.T).max() == pytest.approx(0.0, abs=1e-12)

    rng = np.random.default_rng(0)
    for _ in range(5):
        v = rng.normal(size=K.shape[0])
        assert float(v @ (K @ v)) >= -1e-9

    constant = np.ones(K.shape[0])
    assert float(constant @ (K @ constant)) == pytest.approx(0.0, abs=1e-8)


def test_branch_indicator_is_flat_on_a_straight_bar_and_rises_on_a_branch():
    """
    Le coeur du mecanisme : |div(X)| distingue un CORPS de tube d'un
    BRANCHEMENT. Dans un tube de section constante le front de distance
    n'accelere pas ; a un branchement son aire change d'un coup.

    Comparaison faite sur les CORPS des deux modeles, en excluant les
    tetraedres de frontiere - l'indicateur s'y eleve aussi (source plateau,
    extremites libres), faux positif connu, mesure et documente, qui n'est pas
    ce que ce test cherche a verrouiller.
    """
    bar_indicator, bar_boundary = _indicator_for(_straight_bar_data())
    l_indicator, l_boundary = _indicator_for(_thick_l_data())

    bar_core = bar_indicator[~bar_boundary]
    l_core = l_indicator[~l_boundary]
    assert len(bar_core) > 10 and len(l_core) > 10

    # Le barreau droit n'a AUCUN branchement : son coeur doit rester bas.
    # Le L en a un, et il doit ressortir nettement au-dessus.
    assert np.percentile(l_core, 90) > 3.0 * np.percentile(bar_core, 90)


def test_branch_indicator_ignores_unusable_flow(l_shape_data):
    """
    Un tetraedre dont le flux geodesique n'est pas exploitable ne dit rien :
    lui attribuer 0 le laisse au lissage nominal, plutot que de lui inventer
    une ramification. Sans source plateau, TOUT est inexploitable.
    """
    _mesh_data, tet_mesh, _overhang_map, _regions = l_shape_data
    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)

    indicator, _boundary = _branch_indicator(
        tet_mesh, shape_gradients, volumes,
        np.array([], dtype=np.int64), FieldSolverParams(),
    )

    assert np.all(indicator == 0.0)


# ---------------------------------------------------------------------------
# Borne de plongee (max_downward_growth_angle_deg)
# ---------------------------------------------------------------------------

def test_clamp_downward_preserves_azimuth_and_bounds_the_dive():
    """
    La borne APLATIT une direction descendante, elle ne la retourne pas :
    l'azimut est conserve. C'est ce qui la distingue d'`enforce_upward_growth`,
    dont le flip fait pointer la croissance a l'oppose de ce que la geometrie
    demande.
    """
    up = np.array([0.0, 0.0, 1.0])
    # Plongee a 60 deg sous l'horizontale, azimut a 45 deg dans le plan XY.
    horizontal = np.array([np.cos(np.pi / 4), np.sin(np.pi / 4), 0.0])
    steep = horizontal * np.cos(np.radians(60)) - up * np.sin(np.radians(60))
    directions = np.vstack([steep, np.array([0.0, 0.0, 1.0])])

    clamped, mask = _clamp_downward_directions(directions, up, 30.0, 1e-9)

    assert mask[0] and not mask[1]
    dive = np.degrees(np.arcsin(-clamped[0, 2] / np.linalg.norm(clamped[0])))
    assert dive == pytest.approx(30.0, abs=1e-6)

    azimuth_before = steep[:2] / np.linalg.norm(steep[:2])
    azimuth_after = clamped[0, :2] / np.linalg.norm(clamped[0, :2])
    assert azimuth_after == pytest.approx(azimuth_before, abs=1e-9)

    # Une direction deja montante n'est jamais touchee.
    assert clamped[1] == pytest.approx(directions[1])


def test_clamp_downward_at_ninety_degrees_allows_free_descent():
    """90 deg = descente verticale libre : la borne ne doit plus rien brider."""
    up = np.array([0.0, 0.0, 1.0])
    directions = np.array([[0.0, 0.0, -1.0], [0.3, -0.4, -0.86]])
    directions = directions / np.linalg.norm(directions, axis=1)[:, None]

    clamped, _mask = _clamp_downward_directions(directions, up, 90.0, 1e-9)

    for before, after in zip(directions, clamped):
        assert after == pytest.approx(before, abs=1e-9)


def test_clamp_downward_disabled_by_default_leaves_directions_untouched():
    """None = borne inactive, comportement historique strictement inchange."""
    assert FieldSolverParams().max_downward_growth_angle_deg is None


# ---------------------------------------------------------------------------
# Distance de transition entre branches (transition_distance_mm)
# ---------------------------------------------------------------------------

def _two_direction_field(tet_mesh, split_z: float = 12.0):
    """
    Champ de directions cible synthetique a DEUX zones tranchees : +Z en
    dessous de `split_z`, -Y au-dessus. Reproduit la situation reelle
    (deux branches d'orientation differente separees par une frontiere
    nette) sans dependre d'une segmentation particuliere.
    """
    centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    directions = np.tile(np.array([0.0, 0.0, 1.0]), (tet_mesh.tet_count, 1))
    directions[centers[:, 2] >= split_z] = np.array([0.0, -1.0, 0.0])
    return directions, centers


def test_transition_distance_disabled_by_default(l_shape_data):
    """Opt-in : le defaut ne doit jamais emprunter ce chemin de code."""
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data

    assert FieldSolverParams().transition_distance_mm is None

    field = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [], FieldSolverParams()
    )
    assert field.values.shape == (tet_mesh.node_count,)


def test_diffusion_leaves_a_uniform_direction_field_untouched(l_shape_data):
    """
    Propriete qui rend le mecanisme AUTO-LIMITE aux jonctions : la ou toutes
    les directions cibles sont identiques (l'interieur d'une branche),
    moyenner un tetraedre avec ses voisins ne peut rien changer. C'est ce qui
    dispense d'un masque geometrique ou d'une identification explicite de
    "la branche d'en face".
    """
    _mesh_data, tet_mesh, _overhang_map, _regions = l_shape_data
    uniform = np.tile(np.array([0.0, 0.0, 1.0]), (tet_mesh.tet_count, 1))

    diffused, changed = _diffuse_tet_directions(tet_mesh, uniform, 10.0)

    assert changed == 0
    assert diffused == pytest.approx(uniform, abs=1e-9)


def test_diffusion_widens_the_transition_between_two_directions(l_shape_data):
    """
    Le coeur du mecanisme : plus `transition_distance_mm` est grand, plus la
    zone ou la direction cible est INTERMEDIAIRE (ni +Z ni -Y) s'etend.

    Mesure de la largeur : nombre de tetraedres dont la direction s'ecarte
    notablement des DEUX directions d'origine. Une transition nette n'en
    produit presque aucun ; une transition large en produit beaucoup.
    """
    _mesh_data, tet_mesh, _overhang_map, _regions = l_shape_data
    directions, _centers = _two_direction_field(tet_mesh)

    up = np.array([0.0, 0.0, 1.0])
    side = np.array([0.0, -1.0, 0.0])

    def intermediate_count(field):
        aligned_up = np.abs(field @ up) > 0.99
        aligned_side = np.abs(field @ side) > 0.99
        return int(np.count_nonzero(~(aligned_up | aligned_side)))

    widths = [
        intermediate_count(_diffuse_tet_directions(tet_mesh, directions, d)[0])
        for d in (0.5, 2.0, 8.0)
    ]

    assert widths[0] < widths[1] < widths[2], (
        f"la zone de transition doit s'elargir avec la distance demandee, "
        f"obtenu {widths}"
    )


def test_diffusion_keeps_directions_normalised(l_shape_data):
    """
    Une direction cible est un vecteur UNITAIRE partout ailleurs dans step4
    (`_assemble_gradient_direction_term` renormalise, mais le redressement et
    la borne de plongee supposent deja des vecteurs unitaires en entree) :
    la diffusion ne doit pas rendre la main avec des normes quelconques.
    """
    _mesh_data, tet_mesh, _overhang_map, _regions = l_shape_data
    directions, _centers = _two_direction_field(tet_mesh)

    diffused, changed = _diffuse_tet_directions(tet_mesh, directions, 5.0)

    assert changed > 0
    assert np.linalg.norm(diffused, axis=1) == pytest.approx(
        np.ones(tet_mesh.tet_count), abs=1e-9
    )


def test_transition_distance_changes_the_solved_field():
    """
    Verrouille que le parametre agit REELLEMENT sur le champ resolu, et pas
    seulement sur un intermediaire. Un mecanisme dont le parametre ne change
    pas la sortie est precisement ce qui a ete mesure - puis retire - pour la
    premiere version de ce controle (modulation du poids de lissage : le
    champ changeait, l'orientation des couches non).

    Fixture DEDIEE avec `max_region_merge_angle_deg` : la forme en L par
    defaut ne produit qu'UNE region, donc une seule direction cible, donc
    une diffusion legitimement inerte (verrouille par
    `test_diffusion_leaves_a_uniform_direction_field_untouched`). Il faut
    deux directions distinctes pour que ce controle ait quoi que ce soit a
    etaler - mesure : 1 region sans ce parametre, 4 regions et 2 directions
    distinctes avec.
    """
    vertical = trimesh.creation.box(extents=[4, 4, 20])
    vertical.apply_translation([0, 0, 10])
    horizontal = trimesh.creation.box(extents=[20, 4, 4])
    horizontal.apply_translation([8, 0, 20])
    combined = vertical.union(horizontal, engine="manifold")

    import os
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as handle:
        path = handle.name
    combined.export(path)
    try:
        mesh_data, tet_mesh, _report = ingest_run(path, resolution=2.0)
    finally:
        os.unlink(path)

    overhang_map = overhang_run(
        mesh_data, OverhangParams(overhang_angle_threshold_deg=45.0)
    )
    regions = segmentation_run(
        mesh_data,
        overhang_map,
        SegmentationParams(
            platform_z=0.0,
            platform_contact_tolerance=0.5,
            max_region_merge_angle_deg=30.0,
        ),
    )
    distinct = {tuple(np.round(r.target_direction, 6)) for r in regions.regions.values()}
    assert len(distinct) >= 2, (
        "la fixture doit produire au moins deux directions cibles distinctes, "
        "sinon il n'y a rien a etaler et le test ne verrouille rien"
    )

    sharp = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(transition_distance_mm=0.5),
    )
    wide = field_run(
        tet_mesh, mesh_data, overhang_map, regions, [],
        FieldSolverParams(transition_distance_mm=15.0),
    )

    assert sharp.gradient is not None and wide.gradient is not None
    assert not np.allclose(sharp.gradient, wide.gradient, atol=1e-6)


# ---------------------------------------------------------------------------
# Reutilisation de la sonde du diagnostic de croissance montante
# ---------------------------------------------------------------------------
#
# `build_static_system` resout deja le systeme structurel pour mesurer le
# residu descendant sur le champ REELLEMENT resolu. Quand les contraintes
# utilisateur n'apportent rien - cas par defaut du pipeline -, le systeme que
# `solve_with_constraints` doit resoudre est le MEME. Le resoudre une seconde
# fois coutait 2,0 s sur les 9,45 s de step4 (death.stl @ resolution 0.8).

def test_probe_is_reused_when_constraints_add_nothing(l_shape_data):
    """
    Sans contrainte, `solve_with_constraints` ne doit PAS relancer le
    solveur : il rend la sonde deja calculee, au bit pres.

    L'A/B se fait dans le meme processus, sur le meme systeme statique : on
    lit d'abord le chemin reutilise, puis on efface la sonde pour forcer le
    chemin historique. Comparer deux processus ne prouverait rien sur
    l'egalite bit a bit.
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data
    params = FieldSolverParams(enforce_upward_growth=True)

    static_system = build_static_system(tet_mesh, mesh_data, overhang_map, regions, params)
    assert static_system.probe_values is not None, (
        "le diagnostic de croissance montante etant actif, la sonde doit exister"
    )

    reused = solve_with_constraints(static_system, None, params)

    static_system.probe_values = None  # force la resolution complete
    recomputed = solve_with_constraints(static_system, None, params)

    assert np.array_equal(reused.values, recomputed.values)
    assert np.array_equal(reused.gradient, recomputed.gradient)


def test_probe_reuse_never_aliases_the_cached_array(l_shape_data):
    """
    Le tableau rendu doit etre une COPIE : `solve_with_constraints` est
    appelable plusieurs fois sur le meme systeme statique (son role apres
    edition d'un point d'attache), et une mutation par l'appelant ne doit
    jamais revenir dans le resultat de l'appel suivant.
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data
    params = FieldSolverParams(enforce_upward_growth=True)

    static_system = build_static_system(tet_mesh, mesh_data, overhang_map, regions, params)
    first = solve_with_constraints(static_system, None, params)
    reference = first.values.copy()

    first.values[:] = -12345.0  # mutation hostile par l'appelant

    second = solve_with_constraints(static_system, None, params)
    assert np.array_equal(second.values, reference)


def test_probe_is_not_reused_when_a_constraint_changes_the_system(l_shape_data):
    """
    Des qu'une contrainte contribue reellement, le systeme n'est plus celui
    de la sonde : le solveur DOIT etre relance. Sans cette condition, une
    edition de point d'attache serait silencieusement ignoree - le pire des
    defauts possibles pour ce chemin.
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data
    params = FieldSolverParams(enforce_upward_growth=True)

    static_system = build_static_system(tet_mesh, mesh_data, overhang_map, regions, params)
    without = solve_with_constraints(static_system, None, params)

    constraint = Constraint(
        type=ConstraintType.ATTACH_POINT,
        position=tet_mesh.tet_vertices.mean(axis=0),
        target_field_value=99.0,
        influence_radius=float(np.ptp(tet_mesh.tet_vertices[:, 2])) * 0.5,
        falloff=ConstraintFalloff.LINEAR,
        weight=1e6,
    )
    with_constraint = solve_with_constraints(static_system, [constraint], params)

    assert not np.allclose(without.values, with_constraint.values), (
        "la contrainte doit deplacer le champ, donc interdire la reutilisation"
    )


def test_no_probe_when_upward_growth_diagnostic_is_disabled(l_shape_data):
    """
    Sans diagnostic actif, aucune sonde n'est calculee : le comportement
    historique (une seule resolution, dans solve_with_constraints) est
    strictement conserve.
    """
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data
    params = FieldSolverParams()  # enforce_upward_growth=False par defaut

    static_system = build_static_system(tet_mesh, mesh_data, overhang_map, regions, params)
    assert static_system.probe_values is None
    assert static_system.upward_growth_diagnostic is None

    field = solve_with_constraints(static_system, None, params)
    assert field.values.shape == (tet_mesh.node_count,)


def test_radial_tilt_projection_keeps_only_reachable_component():
    """Phase 3.1 : sur la machine actuelle, la normale ne peut pencher que dans
    le plan (radial, vertical) du point, vers l'exterieur (limite #59)."""
    import numpy as np
    from core.mesh_types import TetMesh
    from pipeline.step4_field_solver import _project_to_reachable_radial_tilt
    base = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    # deux tetraedres : l'un centre vers +X (radiale +X), l'autre sur l'axe
    verts = np.vstack([base + [20, 0, 0], base - 0.25])
    tm = TetMesh(tet_vertices=verts, tetrahedra=np.array([[0, 1, 2, 3], [4, 5, 6, 7]]), tet_adjacency=None)
    d = np.array([[0.3, 0.5, 0.8], [0.6, 0.6, 0.5]])
    d /= np.linalg.norm(d, axis=1)[:, None]
    out, changed = _project_to_reachable_radial_tilt(tm, d, np.zeros(3), np.array([0.0, 0.0, 1.0]))
    rhat = np.array([20.25, 0.25, 0.0]) / np.hypot(20.25, 0.25)
    tang = np.cross([0, 0, 1.0], rhat)
    assert abs(out[0] @ tang) < 1e-12 and out[0] @ rhat >= 0 and abs(np.linalg.norm(out[0]) - 1) < 1e-12
    assert changed[0] and not changed[1]            # sur l'axe : intacte
    inward = -0.8 * rhat + 0.6 * np.array([0, 0, 1.0])
    out, _ = _project_to_reachable_radial_tilt(tm, np.vstack([inward, d[1]]), np.zeros(3), np.array([0.0, 0.0, 1.0]))
    assert np.allclose(out[0], [0, 0, 1])           # penche vers le centre : ramenee a la verticale
