"""
tests/test_step1b_junction_detection.py

Tests pour pipeline/step1b_junction_detection.py, sur un mesh synthetique
genere a la volee (aucun fichier externe requis) - meme geometrie de
jonction multi-region que la fixture `tripod_data` de
tests/test_step5_layers.py (pied + 3 bras non coplanaires), mais n'a besoin
ici que de MeshData + RegionDependencyGraph (pas de tetraedrisation fine ni
de champ scalaire, cf. la nature "passe legere" de la detection).
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest
import trimesh

from core.mesh_types import Region, RegionDependencyGraph
from pipeline.step1_ingest import run as ingest_run
from pipeline.step1b_junction_detection import detect_junction_zones


def _tripod_mesh() -> trimesh.Trimesh:
    """Meme geometrie que tests/test_step5_layers.py::_tripod_mesh."""
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
def tripod_mesh_data_and_graph():
    """
    MeshData + RegionDependencyGraph (4 regions injectees manuellement,
    memes directions divergentes que tripod_data) sur la geometrie tripode
    - pas de tetraedrisation fine necessaire, `detect_junction_zones` ne
    depend que de la surface.
    """
    mesh = _tripod_mesh()

    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, _tet_mesh, _report = ingest_run(path, resolution=1.5)
    finally:
        os.unlink(path)

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
    region_graph = RegionDependencyGraph(regions=regions, dependencies=[(0, 1), (0, 2), (0, 3)])
    return mesh_data, region_graph


def test_detects_zone_at_tripod_junction(tripod_mesh_data_and_graph):
    """
    Les 3 bras du tripode (regions 1, 2, 3) ont des target_direction
    ecartees d'environ 90-120 degres les unes des autres, largement au-dela
    d'un seuil de 30 degres - detect_junction_zones doit trouver au moins
    une zone, centree pres du sommet du pied (z~=20, la jonction physique).
    """
    mesh_data, region_graph = tripod_mesh_data_and_graph

    zones = detect_junction_zones(
        mesh_data, region_graph, angle_threshold_deg=30.0, margin=3.0
    )

    assert len(zones) >= 1, "Le tripode devrait produire au moins une zone de jonction"
    for zone in zones:
        assert zone.radius > 0.0
        # Le sommet du pied (jonction physique des 3 bras) est a z=20 -
        # chaque zone detectee doit etre proche de cette hauteur, pas
        # ailleurs sur le mesh.
        assert 10.0 < zone.centroid[2] < 30.0, (
            f"Zone centree a z={zone.centroid[2]:.2f}, attendu pres de z=20 "
            "(jonction physique des 3 bras)"
        )


def test_no_zone_below_angle_threshold(tripod_mesh_data_and_graph):
    """
    Avec un seuil tres eleve (jamais depasse par la geometrie du tripode,
    angle max mesure < 180 deg par construction), aucune zone ne doit etre
    detectee - comportement par defaut sans faux positif.
    """
    mesh_data, region_graph = tripod_mesh_data_and_graph

    zones = detect_junction_zones(
        mesh_data, region_graph, angle_threshold_deg=179.0, margin=3.0
    )

    assert zones == []


def test_single_region_produces_no_zones():
    """
    Un mesh a une seule region (aucune frontiere inter-region possible) ne
    doit jamais produire de zone, quel que soit le seuil - garde-fou
    trivial mais explicite.
    """
    mesh = trimesh.creation.box(extents=[10, 10, 10])
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, _tet_mesh, _report = ingest_run(path, resolution=2.0)
    finally:
        os.unlink(path)

    region = Region(
        id=0, face_ids=np.arange(len(mesh_data.faces)), target_direction=np.array([0.0, 0.0, 1.0])
    )
    region_graph = RegionDependencyGraph(regions={0: region}, dependencies=[])

    zones = detect_junction_zones(mesh_data, region_graph, angle_threshold_deg=30.0, margin=3.0)

    assert zones == []
