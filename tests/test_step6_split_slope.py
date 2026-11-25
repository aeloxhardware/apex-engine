"""Decoupage en troncons montants limite aux anneaux inclines (`Step6Params.run_split_min_slope_deg`, limite #92)."""
import numpy as np
import pytest

from core import line_types, runtime
from core.mesh_types import Layer
from element_ardor import step6_toolpath2d as s6
from element_ardor.step6_bridge_runs import _max_slope, _split_rings_with_kinds
from element_ardor.step6_params import Step6Params


def _tilted_layer(angle_deg: float, half: float = 8.0) -> Layer:
    """Couche carree dont le contour monte en pente le long de x (repere local horizontal) : l'anneau monte d'un cote
    et redescend de l'autre."""
    t = np.linspace(-half, half, 41)
    side = [np.column_stack([t, np.full_like(t, -half)]), np.column_stack([np.full_like(t, half), t]),
            np.column_stack([t[::-1], np.full_like(t, half)]), np.column_stack([np.full_like(t, -half), t[::-1]])]
    xy = np.vstack([s[:-1] for s in side] + [side[0][:1]])
    z = 5.0 + np.tan(np.radians(angle_deg)) * np.abs(xy[:, 0])
    frame = (np.array([0.0, 0.0, 5.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    return Layer(id=0, contour_3d=[np.column_stack([xy, z])], region_ids=[0], local_frame=frame, isovalue=5.0)


class _Map:
    is_bridge = np.zeros(1, dtype=bool)
    bridge_span_direction = np.zeros((1, 3))


class _Mesh:
    """Une seule face, jamais un pont ; kd-tree natif si le noyau est actif (sinon le chemin natif est ecarte)."""
    def __init__(self):
        from core.spatial import KDTree
        self.kdtree = KDTree(np.array([[1e6, 1e6, 1e6]]))


def _params(min_slope: float, native: bool) -> Step6Params:
    return Step6Params(out_of_plane_interpolation_power=2.0, split_rings_into_runs=True, continuous_rings=True,
                       run_split_min_amplitude=0.05, run_split_min_slope_deg=min_slope,
                       native_layer=native, native_fill=native)


def _walls(layer: Layer) -> list:
    return [np.asarray(p) for p, k in zip(layer.toolpath_2d, layer.toolpath_kinds)
            if k in (line_types.WALL_OUTER, line_types.WALL_INNER)]


def test_max_slope_is_measured_over_a_window():
    s = np.linspace(0.0, 10.0, 101)
    ramp = np.column_stack([s, np.zeros_like(s), 0.2 * s])                      # 11,3 deg
    assert _max_slope(ramp) == pytest.approx(0.2, rel=1e-6)
    noise = ramp.copy()
    noise[50, 2] += 0.03                                                          # un point decale, 0,1 mm de long
    assert _max_slope(noise) < 0.25                                               # la fenetre de 1 mm l'amortit


def test_a_gentle_ring_is_kept_whole_a_steep_one_is_split():
    s = np.linspace(0.0, 2 * np.pi, 200)
    ring = lambda a: np.column_stack([10 * np.cos(s), 10 * np.sin(s), np.tan(np.radians(a)) * 10 * np.abs(np.cos(s))])
    gentle, steep = ring(5.0), ring(40.0)
    kinds = [line_types.WALL_OUTER, line_types.WALL_OUTER]
    out, _k = _split_rings_with_kinds([gentle, steep], kinds, share_vertices=True, min_amplitude=0.05,
                                      min_slope_deg=20.0)
    assert out[0] is gentle and len(out) > 2
    all_split, _k = _split_rings_with_kinds([gentle, steep], kinds, share_vertices=True, min_amplitude=0.05)
    assert len(all_split) > len(out)                                              # 0 = historique : tout decoupe


@pytest.mark.parametrize("native", [False, True])
def test_layers_under_the_threshold_keep_whole_walls(native):
    with runtime.use(runtime.RuntimeOptions(native=native)):
        gentle = s6.run([_tilted_layer(5.0)], _Mesh(), _Map(), _params(20.0, native))[0]
        legacy = s6.run([_tilted_layer(5.0)], _Mesh(), _Map(), _params(0.0, native))[0]
        steep = s6.run([_tilted_layer(40.0)], _Mesh(), _Map(), _params(20.0, native))[0]
        steep_legacy = s6.run([_tilted_layer(40.0)], _Mesh(), _Map(), _params(0.0, native))[0]
    assert len(_walls(gentle)) < len(_walls(legacy))                               # douce : anneaux entiers
    assert len(_walls(steep)) == len(_walls(steep_legacy))                         # raide : decoupe comme avant
    count = lambda layer: sum(len(p) for p in _walls(layer))
    assert count(gentle) <= count(legacy)                                         # aucun point ajoute


def test_native_and_python_give_the_same_walls_with_the_threshold():
    from core import _cpp
    if not _cpp.available():
        pytest.skip("noyau natif absent")
    out = {}
    layers = [_tilted_layer(5.0), _tilted_layer(40.0)]
    for native in (False, True):
        with runtime.use(runtime.RuntimeOptions(native=native)):
            mesh = _Mesh()
            taken = s6._fill_layers_native(layers, mesh, _Map(), _params(20.0, native), {5.0: 0})
            assert (taken is not None) == native                                  # le chemin natif est bien pris
            out[native] = s6.run(layers, mesh, _Map(), _params(20.0, native))
    for a, b in zip(out[False], out[True]):
        wa, wb = _walls(a), _walls(b)
        assert len(wa) == len(wb)
        for p, q in zip(sorted(wa, key=len), sorted(wb, key=len)):
            assert len(p) == len(q)
