"""Couches planes partout ou l'impression plane suffit (element_ardor/step4_planar.py), sur des maillages synthetiques
passes par le vrai pipeline (etapes 1 a 3)."""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest
import trimesh

from element_ardor import step4_planar as P
from element_ardor.step1_ingest import run as ingest_run
from element_ardor.step2_overhang import OverhangParams, run as overhang_run
from element_ardor.step3_segmentation import SegmentationParams, run as segmentation_run
from element_ardor.step4_field_solver import FieldSolverParams, _identify_platform_nodes, run as field_run

UP = np.array([0.0, 0.0, 1.0])


def _data(mesh: trimesh.Trimesh, resolution: float = 1.5):
    mesh.apply_translation([0, 0, -mesh.bounds[0][2]])
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, tet_mesh, _ = ingest_run(path, resolution=resolution)
    finally:
        os.unlink(path)
    om = overhang_run(mesh_data, OverhangParams())
    rdg = segmentation_run(mesh_data, om, SegmentationParams(platform_z=0.0, platform_contact_tolerance=0.3))
    return mesh_data, tet_mesh, om, rdg


def _inputs(data) -> P.FieldInputs:
    md, tm, om, rdg = data
    return P.FieldInputs(tm, md, None, None, _identify_platform_nodes(tm, md, rdg), rdg.platform_contact_faces)


def _box(ext, at):
    b = trimesh.creation.box(extents=ext)
    b.apply_translation(at)
    return b


def test_a_box_needs_no_tilt_and_its_field_is_the_height():
    data = _data(_box([10, 10, 10], [0, 0, 5]))
    assert len(P.tilt_points(_inputs(data), P.PlanarParams(), UP)) == 0
    md, tm, om, rdg = data
    sf = field_run(tm, md, om, rdg, FieldSolverParams(planar=P.PlanarParams()))
    z = tm.tet_vertices[:, 2]
    assert np.max(np.abs(sf.values - (z - z.min()))) < 0.02


def test_a_bridge_between_two_pillars_is_printed_flat():
    """Plafond de 12 mm tendu entre deux piliers : un pont, aucune inclinaison."""
    mesh = trimesh.util.concatenate([_box([3, 6, 10], [-7.5, 0, 5]), _box([3, 6, 10], [7.5, 0, 5])])
    mesh = mesh.union(_box([18, 6, 2], [0, 0, 11]), engine="manifold")
    assert len(P.tilt_points(_inputs(_data(mesh)), P.PlanarParams(), UP)) == 0


def test_a_cantilever_needs_tilt_under_its_free_end():
    """Bras de 12 mm en porte-a-faux (appui d'un seul cote) : a incliner au-dela de la levre."""
    mesh = _box([4, 4, 14], [0, 0, 7]).union(_box([16, 4, 3], [6, 0, 15.5]), engine="manifold")
    pts = P.tilt_points(_inputs(_data(mesh)), P.PlanarParams(), UP)
    assert len(pts) > 0
    assert pts[:, 0].max() > 8.0                       # vers le bout libre
    assert np.all(pts[:, 0] > 2.0 + P.PlanarParams().lip - 0.3)   # pas sur la levre contre le pilier


def test_a_small_lip_and_a_gentle_slope_are_printed_flat():
    """Rebord de 0,8 mm (levre, 1,13 mm aux coins) et pente a 45 deg : imprimables a plat."""
    lip = _box([10, 10, 8], [0, 0, 4]).union(_box([11.6, 11.6, 1], [0, 0, 8.5]), engine="manifold")
    assert len(P.tilt_points(_inputs(_data(lip)), P.PlanarParams(), UP)) == 0
    from shapely.geometry import Polygon
    slope = trimesh.creation.extrude_polygon(Polygon([(0, 0), (4, 0), (10, 6), (0, 6)]), 6.0)   # profil en (x, z)
    slope.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, [1, 0, 0]))
    assert len(P.tilt_points(_inputs(_data(slope)), P.PlanarParams(), UP)) == 0


def test_first_layers_stay_flat_and_weight_is_smooth():
    mesh = _box([4, 4, 14], [0, 0, 7]).union(_box([16, 4, 3], [6, 0, 15.5]), engine="manifold")
    data = _data(mesh)
    inputs = _inputs(data)
    w = P.tilt_weight(inputs, P.PlanarParams(), UP)
    tm = data[1]
    c = tm.tet_vertices[tm.tetrahedra].mean(axis=1)
    assert np.all(w[c[:, 2] < P.PlanarParams().flat_height] == 0.0)
    assert w.max() > 0.99 and 0.0 <= w.min()
