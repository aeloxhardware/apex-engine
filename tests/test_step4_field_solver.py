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
import trimesh

from core.mesh_types import Constraint, ConstraintFalloff, ConstraintType
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import (
    FieldSolverParams,
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
