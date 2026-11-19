"""Couches pleines de dessus et de dessous (element_ardor/step6_solid.py) sur des couches synthetiques planes."""
import numpy as np
import pytest
from scipy.spatial import cKDTree

from core import line_types, runtime
from core.mesh_types import Layer
from element_ardor import step6_toolpath2d as s6
from element_ardor.step6_params import Step6Params


def _square(half: float, z: float) -> np.ndarray:
    pts = np.array([[-half, -half], [half, -half], [half, half], [-half, half], [-half, -half]], dtype=float)
    return np.column_stack([pts, np.full(len(pts), z)])


def _stack(halves: list) -> list:
    """Une couche plane par demi-cote, 0,2 mm d'ecart, isovaleur = hauteur."""
    out = []
    for k, half in enumerate(halves):
        z = 0.2 * k
        frame = (np.array([0.0, 0.0, z]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
        out.append(Layer(id=k, contour_3d=[_square(half, z)], region_ids=[0], local_frame=frame, isovalue=z))
    return out


class _Map:
    is_bridge = np.zeros(1, dtype=bool)
    bridge_span_direction = np.zeros((1, 3))


class _Mesh:
    kdtree = cKDTree(np.array([[1e6, 1e6, 1e6]]))     # une seule face, jamais un pont


def _infill_length(layer: Layer) -> float:
    total = 0.0
    for poly, kind in zip(layer.toolpath_2d, layer.toolpath_kinds):
        if kind == line_types.INFILL and len(poly) >= 2:
            total += float(np.hypot(*np.diff(np.asarray(poly)[:, :2], axis=0).T).sum())
    return total


@pytest.mark.parametrize("native", [False, True])
def test_top_and_bottom_layers_are_solid(native):
    layers = _stack([10.0] * 12)
    params = Step6Params(solid_top_layers=3, solid_bottom_layers=2, native_layer=False, native_fill=native)
    with runtime.use(runtime.RuntimeOptions(native=native)):
        out = s6.run(layers, _Mesh(), _Map(), params)
    lengths = [_infill_length(layer) for layer in out]
    sparse = lengths[5]
    for k in (0, 1, 9, 10, 11):
        assert lengths[k] > 3.0 * sparse, (k, lengths)
    for k in (2, 3, 4, 5, 6, 7, 8):
        assert lengths[k] < 1.5 * sparse, (k, lengths)


def test_uncovered_ring_under_a_smaller_layer_is_solid():
    """Un epaulement : la couche 5 (20 mm) porte des couches de 10 mm ; l'anneau decouvert de la couche 5 est plein,
    son centre (couvert) ne l'est pas."""
    layers = _stack([10.0] * 6 + [5.0] * 6)
    params = Step6Params(solid_top_layers=2, solid_bottom_layers=0, native_layer=False, native_fill=False)
    with runtime.use(runtime.RuntimeOptions(native=False)):
        out = s6.run(layers, _Mesh(), _Map(), params)
    layer = out[5]
    pts = np.concatenate([np.asarray(p)[:, :2] for p, k in zip(layer.toolpath_2d, layer.toolpath_kinds)
                          if k == line_types.INFILL and len(p) >= 2])
    ring = np.max(np.abs(pts), axis=1) > 5.5
    centre = np.max(np.abs(pts), axis=1) < 4.5
    assert ring.sum() > 5 * max(centre.sum(), 1)


def test_disabled_by_default_changes_nothing():
    layers = _stack([10.0] * 6)
    with runtime.use(runtime.RuntimeOptions(native=False)):
        base = s6.run(layers, _Mesh(), _Map(), Step6Params(native_layer=False, native_fill=False))
    assert Step6Params().solid_top_layers == 0 and Step6Params().solid_bottom_layers == 0
    lengths = [_infill_length(layer) for layer in base]
    assert max(lengths) < 1.5 * min(lengths)


def test_dense_lines_follow_the_long_axis_of_a_band():
    """Une bande de 1 x 12 mm est remplie dans sa longueur, pas en zigzag a travers ; une zone carree garde l'angle de
    la couche."""
    from shapely.geometry import box
    from element_ardor.step6_solid import _piece_angle
    band = shapely_rotate(box(0, 0, 12, 1), 30.0)
    assert abs(np.degrees(_piece_angle(band, np.pi / 2)) % 180.0 - 30.0) < 1.0
    assert _piece_angle(box(0, 0, 5, 5), np.pi / 2) == np.pi / 2


def shapely_rotate(geom, deg):
    import shapely.affinity
    return shapely.affinity.rotate(geom, deg, origin="centroid")


def test_narrow_ring_is_filled_with_concentric_rings():
    """Anneau de 1,6 mm de large (dessus d'une cheminee) : anneaux fermes concentriques, pas de zigzag en bouts ; une
    bande droite de meme largeur garde des lignes dans sa longueur."""
    from shapely.geometry import Point, box
    from element_ardor.step6_solid import _dense_lines, _narrow_ring
    ring = Point(0, 0).buffer(4.0, quad_segs=32).difference(Point(0, 0).buffer(2.4, quad_segs=32))
    assert _narrow_ring(ring, 0.4)
    lines = _dense_lines(ring, 0.4, 0.0)
    assert len(lines) == 4                                 # 1,6 mm = 4 cordons : 2 anneaux par bord, jointifs
    assert all(np.allclose(np.asarray(q)[0], np.asarray(q)[-1]) for q in lines)
    radii = sorted({round(float(np.hypot(*np.asarray(q)[0])), 1) for q in lines})
    assert radii[0] >= 2.4 + 0.15 and radii[-1] <= 4.0 - 0.15
    assert not _narrow_ring(box(0, 0, 20, 1.6), 0.4)
