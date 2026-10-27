"""
tests/test_step4_field_solver.py

Tests pour pipeline/step4_field_solver.py (champ = distance au plateau, methode de la chaleur), sur des mesh synthetiques generes
a la volee (aucun fichier externe requis).

Execution : pytest tests/test_step4_field_solver.py -v

Les fixtures enchainent le pipeline reel step1 -> step2 -> step3, pour que
step4 recoive exactement les structures qu'il recevra en production.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from core.mesh_types import TetMesh
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import (
    FieldSolverParams,
    _compute_shape_gradients,
    _identify_platform_nodes,
    _project_to_reachable_radial_tilt,
    run as field_run,
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
# Champ de la chaleur : distance au plateau a travers la matiere
# ---------------------------------------------------------------------------

def test_flat_box_field_increases_away_from_platform(flat_box_data):
    """Le champ est resolu, nul sur le plateau et croit avec la hauteur."""
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data
    field = field_run(tet_mesh, mesh_data, overhang_map, regions, FieldSolverParams())

    assert np.all(np.isfinite(field.values))
    z = tet_mesh.tet_vertices[:, 2]
    assert np.corrcoef(z, field.values)[0, 1] > 0.9
    platform_nodes = z < z.min() + 1e-6
    assert np.allclose(field.values[platform_nodes], 0.0, atol=1e-3)


def test_output_shapes_match_mesh_counts(flat_box_data):
    mesh_data, tet_mesh, overhang_map, regions = flat_box_data
    field = field_run(tet_mesh, mesh_data, overhang_map, regions, FieldSolverParams())
    assert field.values.shape == (tet_mesh.node_count,)
    assert field.gradient.shape == (tet_mesh.tet_count, 3)


def test_field_is_a_distance_the_gradient_norm_is_close_to_one(l_shape_data):
    """Par construction (flux normalise) ||grad phi|| ~ 1 : c'est une distance, pas un
    compromis lissage/direction. Mediane pesee par le volume."""
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data
    field = field_run(tet_mesh, mesh_data, overhang_map, regions, FieldSolverParams())
    _, volumes = _compute_shape_gradients(tet_mesh)
    norms = np.linalg.norm(field.gradient, axis=1)
    order = np.argsort(norms)
    cumulative = np.cumsum(volumes[order])
    median = norms[order][np.searchsorted(cumulative, cumulative[-1] / 2)]
    assert 0.7 < median < 1.4


def test_overhanging_arm_is_farther_from_the_platform_than_its_support(l_shape_data):
    """Une piece suspendue est plus loin (a travers la matiere) que ce qui la porte :
    le bout du bras a une valeur superieure a celle du pied a la meme hauteur."""
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data
    field = field_run(tet_mesh, mesh_data, overhang_map, regions, FieldSolverParams())
    v = tet_mesh.tet_vertices
    arm_tip = (v[:, 0] > 14) & (v[:, 2] > 17)
    foot_same_height = (np.abs(v[:, 0]) < 2) & (v[:, 2] > 17)
    assert arm_tip.any() and foot_same_height.any()
    assert field.values[arm_tip].mean() > field.values[foot_same_height].mean() + 5.0


def test_platform_nodes_come_from_the_segmentation_contact(flat_box_data):
    mesh_data, tet_mesh, _overhang_map, regions = flat_box_data
    nodes = _identify_platform_nodes(tet_mesh, mesh_data, regions)
    z = tet_mesh.tet_vertices[:, 2]
    assert len(nodes) > 0 and np.all(z[nodes] <= z.min() + 1e-6)


def test_radial_tilt_projection_keeps_only_reachable_component():
    """Sur la machine actuelle, la normale ne peut pencher que dans le plan (radial,
    vertical) du point, vers l'exterieur."""
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


def test_heat_projection_defaults_and_needs_an_origin():
    """La projection sur les inclinaisons atteignables est active par defaut mais sans effet
    tant que `radial_tilt_origin` n'est pas renseigne."""
    assert FieldSolverParams().heat_project_reachable is True
    assert FieldSolverParams().radial_tilt_origin is None


def test_reachable_projection_changes_the_field_when_an_origin_is_given(l_shape_data):
    mesh_data, tet_mesh, overhang_map, regions = l_shape_data
    free = field_run(tet_mesh, mesh_data, overhang_map, regions, FieldSolverParams())
    reachable = field_run(tet_mesh, mesh_data, overhang_map, regions,
                          FieldSolverParams(radial_tilt_origin=np.array([-30.0, 0.0, 0.0])))
    assert np.all(np.isfinite(reachable.values))
    assert not np.allclose(free.values, reachable.values)
