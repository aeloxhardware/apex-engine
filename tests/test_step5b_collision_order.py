"""step5b : ordre des couches contraint par la tete (opt-in, 2026-10-01)."""
import numpy as np

from pipeline.step5b_collision_order import _densify, _stable_topological_order


def test_no_constraint_keeps_order():
    order, broken = _stable_topological_order(5, set())
    assert order == [0, 1, 2, 3, 4] and broken == 0


def test_constraints_already_satisfied_keep_order():
    order, broken = _stable_topological_order(4, {(0, 2), (1, 3), (0, 3)})
    assert order == [0, 1, 2, 3] and broken == 0


def test_backward_constraint_moves_only_what_it_must():
    # la tete sur la couche 3 toucherait la couche 1 : 3 doit passer avant 1,
    # tout le reste garde son ordre d'origine
    order, broken = _stable_topological_order(5, {(3, 1)})
    assert order == [0, 2, 3, 1, 4] or order == [0, 3, 1, 2, 4]
    assert order.index(3) < order.index(1)
    assert broken == 0


def test_cycle_is_broken_and_reported():
    order, broken = _stable_topological_order(3, {(0, 1), (1, 0)})
    assert sorted(order) == [0, 1, 2]
    assert broken == 1


def test_densify_keeps_endpoints_and_spacing():
    line = np.array([[0.0, 0, 0], [10.0, 0, 0]])
    pts = _densify(line, 1.5)
    assert np.allclose(pts[0], line[0]) and np.allclose(pts[-1], line[1])
    assert np.diff(pts[:, 0]).max() <= 1.5 + 1e-12
