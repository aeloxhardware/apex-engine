"""P3 sur la voie chaleur (pipeline/skin/carve.py) : creusement du corps sous les zones de lissage."""
import dataclasses

import numpy as np
import pytest
import trimesh

import main as M
from pipeline import step1_ingest as s1, step2_overhang as s2, step3_segmentation as s3
from pipeline import step4_field_solver as s4, step5_layers as s5
from pipeline.skin import carve, zones as sz


def test_clip_patch_keeps_exactly_the_positive_part_and_stays_watertight():
    rng = np.random.default_rng(1)
    # grille de 6x6 carres coupes en triangles, psi lineaire en x : la partie conservee est x >= 0,37
    xs = np.linspace(0, 1, 7)
    gx, gy = np.meshgrid(xs, xs, indexing="ij")
    v = np.c_[gx.ravel(), gy.ravel(), np.zeros(gx.size)]
    idx = np.arange(gx.size).reshape(gx.shape)
    faces = []
    for i in range(6):
        for j in range(6):
            a, b, c, d = idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]
            faces += [[a, b, c], [a, c, d]] if rng.random() < 0.5 else [[a, b, d], [b, c, d]]
    faces = np.array(faces)
    psi = v[:, 0] - 0.37
    cv, cf, src = carve.clip_patch(v, faces, psi)
    tri = cv[cf]
    area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1).sum()
    assert area == pytest.approx(1.0 - 0.37, abs=1e-12)
    assert (cv[cf][:, :, 0] >= 0.37 - 1e-12).all()
    assert len(src) == len(cf) and src.max() < len(faces)
    mesh = trimesh.Trimesh(cv, cf, process=False)
    assert len(mesh.outline().entities) == 1 and mesh.outline().entities[0].closed   # bord unique : coupe etanche


def _box_pipeline(extent=12.0, resolution=3.0):
    import os, tempfile
    path = os.path.join(tempfile.mkdtemp(), "box.stl")
    trimesh.creation.box(extents=(extent, extent, extent)).export(path)
    p = M.PipelineParams()
    md, tm, rep = s1.run(path, resolution=resolution, near_platform_resolution=p.near_platform_resolution,
                         near_platform_margin=p.near_platform_margin, place_on_platform=True)
    om = s2.run(md, p.overhang)
    seg = dataclasses.replace(p.segmentation, platform_z=rep.platform_z_used)
    graph = s3.run(md, om, seg)
    field = s4.run(tm, md, om, graph, M._effective_field_solver_params(p))
    return tm, field, s5.run(tm, md, field, graph, p.layers)


def _loop_area(loops):
    total = 0.0
    for loop in loops:
        xy = np.asarray(loop)[:, :2]
        total += 0.5 * np.sum(xy[:-1, 0] * xy[1:, 1] - xy[1:, 0] * xy[:-1, 1])
    return abs(total)


@pytest.fixture(scope="module")
def box():
    return _box_pipeline(12.0)


def _lateral_zone(tm, field, z_max=None):
    solid = sz.PrintedSolid(tm.tet_vertices, tm.tetrahedra, field.values)
    sv, sf, nrm, _areas, centers = sz.surface_mesh(tm, solid, 2.0)
    sel = (np.abs(nrm[:, 2]) < 0.3) & (centers[:, 2] > 1e-3)
    if z_max is not None:
        sel &= centers[:, 2] <= z_max
    return sv, sf[sel]


def test_carved_layer_is_the_original_inset_by_the_thickness(box):
    tm, field, layers = box
    sv, sf = _lateral_zone(tm, field)
    cp = carve.CarveParams(thickness=0.2)
    carved, rep = carve.carve_layers(layers, tm, field, carve.ZoneDistance(sv, sf, cp), cp)
    assert rep.touched_layers == len(layers) and rep.dropped_layers == 0
    mid = [c for c in carved if abs(c.isovalue - 6.0) < 0.15][0]
    original = [l for l in layers if l.id == mid.id][0]
    assert _loop_area(original.contour_3d) == pytest.approx(144.0, rel=0.01)
    assert _loop_area(mid.contour_3d) == pytest.approx((12 - 2 * 0.2) ** 2, rel=0.005)   # parois planes


@pytest.fixture(scope="module")
def fine_box():
    return _box_pipeline(12.0, resolution=1.2)


def test_layers_above_the_zone_are_returned_untouched(fine_box):
    tm, field, layers = fine_box
    sv, sf = _lateral_zone(tm, field, z_max=3.0)        # peau seulement en bas
    cp = carve.CarveParams(thickness=0.2)
    carved, rep = carve.carve_layers(layers, tm, field, carve.ZoneDistance(sv, sf, cp), cp)
    by_id = {c.id: c for c in carved}
    high = [l for l in layers if l.isovalue > 8.0]               # loin de la peau (plusieurs aretes de tetraedre)
    assert high and all(by_id[l.id] is l for l in high)
    low = [l for l in layers if 1.0 < l.isovalue < 2.5]
    assert all(by_id[l.id] is not l for l in low)
    assert rep.touched_layers < len(layers)


def test_empty_zone_changes_nothing(box):
    tm, field, layers = box
    cp = carve.CarveParams()
    zone = carve.ZoneDistance(np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64), cp)
    carved, rep = carve.carve_layers(layers, tm, field, zone, cp)
    assert all(a is b for a, b in zip(layers, carved)) and rep.touched_layers == 0


def test_refine_patch_is_conforming_and_preserves_area():
    # deux triangles voisins, un seul marque : le voisin est coupe par fermeture, sans sommet pendant
    v = np.array([[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0.]])
    f = np.array([[0, 1, 2], [0, 2, 3]])
    nv, nf, parent, red = carve.refine_patch(v, f, np.array([True, False]))
    tri = nv[nf]
    assert 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1).sum() == pytest.approx(4.0)
    assert len(nf) == 4 + 2 and red[parent == 0].all() and not red[parent == 1].any()
    mesh = trimesh.Trimesh(nv, nf, process=False)
    assert mesh.outline().entities[0].closed and len(mesh.outline().entities) == 1   # bord = le carre, rien d'autre


def test_cut_line_follows_the_exact_distance_on_a_small_zone(box):
    """Zone petite devant les tetraedres : la ligne de coupe reste a t de la peau (avant raffinement : jusqu'a 1 mm d'ecart)."""
    tm, field, layers = box
    solid = sz.PrintedSolid(tm.tet_vertices, tm.tetrahedra, field.values)
    sv, sf, nrm, _a, centers = sz.surface_mesh(tm, solid, 2.0)
    sel = (np.abs(nrm[:, 0]) > 0.9) & (nrm[:, 0] > 0) & (np.abs(centers[:, 1]) < 2.0) & (np.abs(centers[:, 2] - 6.0) < 2.0)
    assert sel.any()
    cp = carve.CarveParams(thickness=0.2)
    carved, _rep = carve.carve_layers(layers, tm, field, carve.ZoneDistance(sv, sf[sel], cp), cp)
    zmesh = trimesh.Trimesh(sv, sf[sel], process=False)
    surface = trimesh.Trimesh(tm.tet_vertices, sz.boundary_faces(tm.tetrahedra), process=False)
    errors = []
    for before, after in zip(layers, carved):
        if after is before:
            continue
        pts = np.concatenate([np.asarray(c)[:-1] for c in after.contour_3d])
        _c, ds, _t = trimesh.proximity.closest_point(surface, pts)
        cut = pts[ds > 1e-4]
        if len(cut):
            errors.append(trimesh.proximity.closest_point(zmesh, cut)[1] - 0.2)
    errors = np.abs(np.concatenate(errors))
    assert len(errors) > 20 and np.percentile(errors, 95) < 0.02 and errors.max() < 0.08
