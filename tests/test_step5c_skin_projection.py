"""Etape 5c : projection des contours sur la surface d'origine (R-030 point 8, piste D)."""
import numpy as np
import pytest
import trimesh

from core import _cpp
from core.mesh_types import Layer, MeshData
from element_ardor import step5c_skin_projection as s5c


def _mesh_data(mesh: trimesh.Trimesh) -> MeshData:
    return MeshData(vertices=np.asarray(mesh.vertices), faces=np.asarray(mesh.faces),
                    face_normals=np.asarray(mesh.face_normals), vertex_normals=np.asarray(mesh.vertex_normals),
                    is_watertight=True, bvh=None, kdtree=None)


def _layer(contours, layer_id=0):
    return Layer(id=layer_id, contour_3d=contours, region_ids=[0],
                 local_frame=(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0])))


def _loop_on_box(inset: float, size: float = 10.0, z: float = 5.0, points_per_side: int = 4):
    """Rectangle decale vers l'interieur de `inset` en x et de 1,3 `inset` en y (asymetrique : sans egalite exacte
    de distance entre deux faces, dont le choix serait arbitraire), ferme."""
    a, b = inset, size - inset
    c, d = 1.3 * inset, size - 1.3 * inset
    corners = np.array([[a, c], [b, c], [b, d], [a, d]])
    pts = []
    for k in range(4):
        p, q = corners[k], corners[(k + 1) % 4]
        for t in np.linspace(0, 1, points_per_side, endpoint=False):
            pts.append(p + (q - p) * t)
    pts = np.array(pts)
    loop = np.column_stack([pts, np.full(len(pts), z)])
    return np.vstack([loop, loop[:1]])


@pytest.fixture
def box():
    return trimesh.creation.box(extents=[10, 10, 10], transform=trimesh.transformations.translation_matrix([5, 5, 5]))


def test_contour_off_the_surface_is_brought_back_onto_it(box):
    layers = [_layer([_loop_on_box(0.15)])]
    out, report = s5c.run(layers, _mesh_data(box), s5c.SkinProjectionParams())
    query = trimesh.proximity.ProximityQuery(box)
    dist = query.on_surface(out[0].contour_3d[0])[1]
    assert dist.max() < 1e-9
    assert report.moved > 0 and report.max_shift_applied == pytest.approx(0.195, abs=1e-6)


def test_projection_does_not_mutate_its_input(box):
    loop = _loop_on_box(0.15)
    before = loop.copy()
    layers = [_layer([loop])]
    s5c.run(layers, _mesh_data(box), s5c.SkinProjectionParams())
    assert np.array_equal(layers[0].contour_3d[0], before)


def test_closed_loop_stays_closed_and_midpoints_are_on_the_surface(box):
    layers = [_layer([_loop_on_box(0.15, points_per_side=2)])]        # 8 points : coins coupes
    out, _ = s5c.run(layers, _mesh_data(box), s5c.SkinProjectionParams(chord_tolerance=0.005, max_rounds=4))
    loop = out[0].contour_3d[0]
    assert np.array_equal(loop[0], loop[-1])
    mids = 0.5 * (loop[:-1] + loop[1:])
    assert trimesh.proximity.ProximityQuery(box).on_surface(mids)[1].max() < 0.01       # la corde suit les arretes
    assert len(loop) > 9                                                                 # des points ont ete inseres aux coins


def test_without_refinement_a_corner_is_still_cut(box):
    """Garde : la subdivision adaptative est ce qui rend les coins fideles (sans elle, la corde coupe l'arete)."""
    layers = [_layer([_loop_on_box(0.15, points_per_side=2)])]
    out, _ = s5c.run(layers, _mesh_data(box), s5c.SkinProjectionParams(max_rounds=0))
    loop = out[0].contour_3d[0]
    mids = 0.5 * (loop[:-1] + loop[1:])
    assert trimesh.proximity.ProximityQuery(box).on_surface(mids)[1].max() > 0.05


def test_points_far_from_the_surface_are_left_alone(box):
    """Securite : un point a plus de max_shift de la surface (couche interieure) n'est pas deplace."""
    inner = _loop_on_box(3.0)
    out, report = s5c.run([_layer([inner])], _mesh_data(box), s5c.SkinProjectionParams(max_shift=1.5))
    assert np.allclose(out[0].contour_3d[0][: len(inner) - 1], inner[:-1])
    assert report.refused == len(inner) - 1 and report.moved == 0


def test_open_polyline_keeps_its_endpoints_count_order(box):
    line = np.column_stack([np.linspace(1, 9, 5), np.full(5, 0.1), np.full(5, 5.0)])
    out, _ = s5c.run([_layer([line])], _mesh_data(box), s5c.SkinProjectionParams())
    assert out[0].contour_3d[0].shape[0] >= 5
    assert np.allclose(out[0].contour_3d[0][:, 1], 0.0)
    assert np.allclose(out[0].contour_3d[0][0], [1, 0, 5]) and np.allclose(out[0].contour_3d[0][-1], [9, 0, 5])


def test_native_and_python_distance_agree(box):
    if not _cpp.available() or not hasattr(_cpp.module(), "MeshDistance"):
        pytest.skip("noyau natif absent")
    rng = np.random.default_rng(0)
    pts = rng.uniform(-1.0, 11.0, size=(500, 3))
    previous = _cpp.set_enabled(False)
    try:
        py_pts, py_dist = s5c._Surface(box.vertices, box.faces, box.face_normals, 1000).closest(pts)[:2]
        _cpp.set_enabled(True)
        nat_pts, nat_dist = s5c._Surface(box.vertices, box.faces, box.face_normals, 1000).closest(pts)[:2]
    finally:
        _cpp.set_enabled(previous)
    assert np.allclose(py_dist, nat_dist, atol=1e-9)
    assert np.allclose(py_pts, nat_pts, atol=1e-9)


def test_full_run_is_identical_in_both_native_positions(box):
    if not _cpp.available() or not hasattr(_cpp.module(), "MeshDistance"):
        pytest.skip("noyau natif absent")
    layers = [_layer([_loop_on_box(0.15, points_per_side=3)])]
    previous = _cpp.set_enabled(False)
    try:
        a, ra = s5c.run(layers, _mesh_data(box), s5c.SkinProjectionParams())
        _cpp.set_enabled(True)
        b, rb = s5c.run(layers, _mesh_data(box), s5c.SkinProjectionParams())
    finally:
        _cpp.set_enabled(previous)
    assert a[0].contour_3d[0].shape == b[0].contour_3d[0].shape
    assert np.allclose(a[0].contour_3d[0], b[0].contour_3d[0], atol=1e-9)
    assert (ra.points_out, ra.moved) == (rb.points_out, rb.moved)
