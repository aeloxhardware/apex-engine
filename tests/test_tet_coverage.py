"""
Tests de core/tet_coverage.py : ecart entre la surface et la frontiere du maillage,
et son branchement dans step1 (`check_surface_coverage`).
"""
from __future__ import annotations

import os
import tempfile

import numpy as np
import trimesh

from core.tet_coverage import surface_outside_volume
from pipeline import step1_ingest


def _box_with_tets(extents=(10.0, 10.0, 10.0), resolution=1.5):
    box = trimesh.creation.box(extents=extents)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "box.stl")
        box.export(path)
        mesh_data, tet_mesh, report = step1_ingest.run(path, resolution)
    return box, mesh_data, tet_mesh, report


def test_fully_tetrahedralised_box_has_no_surface_outside_the_volume():
    _, mesh_data, tet_mesh, _ = _box_with_tets()
    rep = surface_outside_volume(mesh_data.vertices, mesh_data.faces, tet_mesh.tet_vertices, tet_mesh.tetrahedra)
    assert rep.far_fraction == 0.0
    assert abs(rep.volume_ratio - 1.0) < 1e-3
    assert rep.clusters == []


def test_removed_tetrahedra_are_reported_as_a_cluster_at_their_location():
    _, mesh_data, tet_mesh, _ = _box_with_tets()
    centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    # on retire tous les tetraedres d'un coin : la peau y reste, le volume non
    corner = np.array([mesh_data.vertices[:, 0].max(), mesh_data.vertices[:, 1].max(), mesh_data.vertices[:, 2].max()])
    removed = np.linalg.norm(centers - corner, axis=1) < 3.0
    assert removed.sum() > 0
    kept = tet_mesh.tetrahedra[~removed]
    rep = surface_outside_volume(mesh_data.vertices, mesh_data.faces, tet_mesh.tet_vertices, kept)
    assert rep.far_fraction > 0.01
    assert len(rep.clusters) >= 1
    # l'amas principal est pres du coin, loin du coin oppose
    assert np.linalg.norm(rep.clusters[0].center - corner) < 6.0
    assert rep.far_area > 0


def test_step1_reports_coverage_only_when_asked():
    box = trimesh.creation.box(extents=(10.0, 10.0, 10.0))
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "box.stl")
        box.export(path)
        _, _, off = step1_ingest.run(path, 1.5)
        _, _, on = step1_ingest.run(path, 1.5, check_surface_coverage=True)
    assert off.surface_outside_volume_fraction is None
    assert on.surface_outside_volume_fraction == 0.0
    assert on.surface_outside_volume_clusters == []
