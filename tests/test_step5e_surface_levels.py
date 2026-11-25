"""Etape 5e : contours des couches courbes pris sur la surface du STL (element_ardor/step5e_surface_levels.py)."""
from dataclasses import dataclass

import numpy as np
import pytest
import trimesh
from scipy.spatial import Delaunay

from core.mesh_types import Layer, ScalarField, TetMesh
from element_ardor import step5e_surface_levels as s5e

R, H, SLOPE = 5.0, 10.0, 0.3          # cylindre ; phi = z + SLOPE * x : couches inclinees de 16,7 deg


@dataclass
class _MeshData:
    vertices: np.ndarray
    faces: np.ndarray


@pytest.fixture(scope="module")
def cylinder():
    stl = trimesh.creation.cylinder(radius=R, height=H, sections=64)
    stl.apply_translation([0.0, 0.0, H / 2])
    g = np.mgrid[-R:R:11j, -R:R:11j, 0:H:11j].reshape(3, -1).T
    inside = g[np.hypot(g[:, 0], g[:, 1]) < R - 0.3]
    pts = np.vstack([stl.vertices, inside])
    tets = Delaunay(pts).simplices
    tm = TetMesh(tet_vertices=pts, tetrahedra=tets)
    field = ScalarField(values=pts[:, 2] + SLOPE * pts[:, 0])
    return _MeshData(np.asarray(stl.vertices), np.asarray(stl.faces)), tm, field


def _true_level(c: float, n: int = 200) -> np.ndarray:
    t = np.linspace(0.0, 2 * np.pi, n)
    x, y = R * np.cos(t), R * np.sin(t)
    return np.column_stack([x, y, c - SLOPE * x])


def _layer(contour: np.ndarray, c: float) -> Layer:
    frame = (np.array([0.0, 0.0, c]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    return Layer(id=0, contour_3d=[contour], region_ids=[0], local_frame=frame, isovalue=c)


def test_phi_is_exact_for_a_linear_field(cylinder):
    md, tm, field = cylinder
    phi = s5e.phi_at(md.vertices, tm, field)
    assert np.abs(phi - (md.vertices[:, 2] + SLOPE * md.vertices[:, 0])).max() < 1e-9


def test_a_crenellated_tilted_contour_is_replaced_by_the_level_line_on_the_stl(cylinder):
    md, tm, field = cylinder
    c = 5.0
    true = _true_level(c)
    jag = true.copy()
    jag[:, :2] *= (1.0 + 0.06 * np.sign(np.sin(15 * np.linspace(0, 2 * np.pi, len(jag)))))[:, None]   # +/- 0,3 mm
    out, report = s5e.run([_layer(jag, c)], md, tm, field, set())
    assert report.loops == 1 and report.replaced == 1
    loop = out[0].contour_3d[0]
    # sur la surface (polygone a 64 cotes : apotheme R cos(pi/64)) et sur la couche phi = c
    assert np.hypot(loop[:, 0], loop[:, 1]).min() > R * np.cos(np.pi / 64) - 1e-6
    assert np.hypot(loop[:, 0], loop[:, 1]).max() < R + 1e-6
    assert np.abs(loop[:, 2] + SLOPE * loop[:, 0] - c).max() < 1e-6
    assert out[0].isovalue == c and out[0].local_frame is not None


def test_loops_already_resliced_are_left_alone(cylinder):
    md, tm, field = cylinder
    layer = _layer(_true_level(5.0), 5.0)
    out, report = s5e.run([layer], md, tm, field, {(0, 0)})
    assert report.loops == 0 and out[0] is layer


def test_a_loop_far_from_any_level_line_is_kept(cylinder):
    md, tm, field = cylinder
    far = _true_level(5.0)
    far[:, :2] *= 0.4                                 # a 3 mm de la surface : aucune correspondante
    out, report = s5e.run([_layer(far, 5.0)], md, tm, field, set())
    assert report.replaced == 0 and out[0].contour_3d[0] is far


def test_simplify_keeps_corners_and_drops_collinear_points():
    t = np.linspace(0.0, 1.0, 11)[:, None]
    a, b, c = np.array([0.0, 0, 0]), np.array([4.0, 0, 1]), np.array([4.0, 3, 1])
    loop = np.vstack([a + t * (b - a), b + t[1:] * (c - b), c + t[1:] * (a - c)])
    out = s5e.simplify(loop, 0.01)
    assert len(out) == 4 and np.allclose(out[[0, -1]], [a, a])


def test_a_level_through_vertices_still_closes():
    """Isovaleur pile sur des sommets du STL (phi = z, sommets a z = 5 apres subdivision) : la ligne reste fermee.
    Avant, une face dont un sommet valait exactement c n'etait pas coupee : ligne trouee, aucune boucle."""
    box = trimesh.creation.box(extents=[8.0, 8.0, 10.0])
    box.apply_translation([0.0, 0.0, 5.0])
    v, f = trimesh.remesh.subdivide_to_size(np.asarray(box.vertices), np.asarray(box.faces), 1.0)
    assert np.isclose(v[:, 2], 5.0).any()
    seg, _hit = s5e.level_segments(v, f, v[:, 2].copy(), 5.0)
    loops = s5e.closed_loops(seg, 1e-3)
    assert len(loops) == 1 and np.allclose(loops[0][:, 2], 5.0)


def test_only_curved_loops_are_redensified():
    params = s5e.SurfaceLevelParams()
    frame = (np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
    t = np.linspace(0, 2 * np.pi, 9)
    flat = np.column_stack([10 * np.cos(t), 10 * np.sin(t), np.full(9, 2.0)])
    bent = flat.copy()
    bent[:, 2] += 0.5 * np.cos(t)
    assert len(s5e._finish(flat, frame, params)) == 9                     # plane : pas de point ajoute
    assert len(s5e._finish(bent, frame, params)) > 100                    # courbe : redecoupee a 0,3 mm
