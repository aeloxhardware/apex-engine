"""Suite découpée de `test_cpp_equivalence.py` (R-030 point 1, limite de 800 lignes par fichier) : tests déplacés tels quels."""

from __future__ import annotations
from test_cpp_equivalence import restore_switch  # noqa: F401  (fixtures)

import dataclasses
import numpy as np
import pytest
from core import _cpp, runtime
import trimesh
import trimesh.ray.ray_triangle
from test_cpp_equivalence import (  # noqa: F401
    _assert_ik_equivalent,
    _current_machine_axes,
    _ik_fixture,
    pytestmark,
)
from test_cpp_equivalence_spatial import (  # noqa: F401
    _step4_inputs,
)


def _star_region(rng):
    # contour lisse (pas d'angle aigu) : rayon modulé par une sinusoïde
    n = int(rng.integers(40, 120))
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    rad = 10.0 + rng.uniform(0.5, 2.0) * np.sin(int(rng.integers(2, 6)) * ang + rng.uniform(0, 6))
    outer = np.column_stack([rad * np.cos(ang), rad * np.sin(ang)])
    holes = []
    if rng.random() < 0.5:
        t = np.linspace(0, 2 * np.pi, 24, endpoint=False)
        holes.append(np.column_stack([2.5 * np.cos(t), 2.5 * np.sin(t)]))
    return outer, holes


def test_fill_region_native_matches_shapely(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "fill_region"):
        pytest.skip("module natif anterieur au pas 15")
    from shapely.geometry import Polygon
    from element_ardor import step6_toolpath2d as s6

    def length(polys):
        return sum(float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum()) for p in polys)

    rng = np.random.default_rng(15)
    codes = {s6.InfillPattern.LINES: 0, s6.InfillPattern.GRID: 1, s6.InfillPattern.CROSSED_ZIGZAG: 2}
    checked = 0
    for k in range(40):
        outer, holes = _star_region(rng)
        pattern = list(codes)[k % 3]
        angle = float(rng.uniform(0, np.pi))
        walls, innermost = s6._offset_walls(outer, holes, 0.4, 2, 1e-6)
        infill = s6._generate_infill(innermost, 0.2, pattern, 0.4, angle)
        valid, n_infill, n_walls = core.fill_region(outer, holes, 0.4, 2, 1e-6, 0.2, codes[pattern], angle)
        assert valid
        py_walls = walls[: len(walls) - 1 - len(holes)]
        assert len(n_walls) == len(py_walls)
        area_py = sum(abs(Polygon(w).area) for w in py_walls)
        area_c = sum(abs(Polygon(w).area) for w in n_walls)
        assert abs(area_c - area_py) <= 1e-3 * area_py
        for w in n_walls:
            assert np.array_equal(w[0], w[-1])          # anneaux fermes, comme shapely
        assert abs(length(n_infill) - length(infill)) <= 1e-3 * max(length(infill), 1.0)
        checked += length(infill) > 0
    assert checked >= 30, "le jeu d'essai doit produire du remplissage"


def test_fill_region_orientation_and_degenerate():
    core = _cpp.module()
    if not hasattr(core, "fill_region"):
        pytest.skip("module natif anterieur au pas 15")
    from shapely.geometry import LinearRing
    square = np.array([[0, 0], [10, 0], [10, 10], [0, 10]], float)
    hole = np.array([[4, 4], [6, 4], [6, 6], [4, 6]], float)
    valid, _infill, walls = core.fill_region(square, [hole], 0.4, 1, 0.0, 0.2, 2, 0.0)
    assert valid and len(walls) == 2
    # orientation de GEOS : exterieur horaire, trou anti-horaire
    assert not LinearRing(walls[0]).is_ccw and LinearRing(walls[1]).is_ccw
    # zone d'aire nulle : ignoree entierement (comme `_offset_walls` -> [])
    flat = np.array([[0, 0], [5, 0], [10, 0]], float)
    valid, infill, walls = core.fill_region(flat, [], 0.4, 2, 0.0, 0.2, 2, 0.0)
    assert not valid and not infill and not walls


@pytest.mark.parametrize("shape", ["hollow", "l_shape", "sphere"])
def test_overhang_map_native_matches_python(restore_switch, shape):
    core = _cpp.module()
    if not hasattr(core, "overhang_map"):
        pytest.skip("module natif anterieur au pas 21a")
    import os
    import tempfile
    import trimesh
    from element_ardor import step2_overhang as s2
    from element_ardor.step1_ingest import run as ingest_run

    if shape == "hollow":
        mesh = trimesh.creation.box(extents=[10, 10, 10]).difference(
            trimesh.creation.box(extents=[4, 4, 4]), engine="manifold")
    elif shape == "l_shape":
        v = trimesh.creation.box(extents=[4, 4, 20]); v.apply_translation([0, 0, 10])
        h = trimesh.creation.box(extents=[20, 4, 4]); h.apply_translation([8, 0, 20])
        mesh = v.union(h, engine="manifold")
    else:
        mesh = trimesh.creation.icosphere(subdivisions=2, radius=5.0)
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, _, _ = ingest_run(path, resolution=2.0)
    finally:
        os.unlink(path)
    params = s2.OverhangParams()
    _cpp.set_enabled(False)
    ref = s2.run(mesh_data, params)
    _cpp.set_enabled(True)
    assert s2._native_overhang_map(mesh_data, params) is not None, "le natif doit avoir tourne"
    nat = s2.run(mesh_data, params)
    assert np.array_equal(ref.is_overhang, nat.is_overhang)
    assert np.array_equal(ref.is_bridge, nat.is_bridge)
    assert np.allclose(ref.risk_score, nat.risk_score, rtol=0, atol=1e-12)
    assert np.array_equal(np.isnan(ref.bridge_span_direction), np.isnan(nat.bridge_span_direction))
    # La portee est un axe principal : unique (au signe canonique pres) sauf
    # quand deux variances sont egales (plafond carre). On verifie donc ce qui
    # est garanti : meme variance le long des deux directions, par composante.
    bridges = np.flatnonzero(ref.is_bridge)
    for row in np.unique(ref.bridge_span_direction[bridges], axis=0):
        faces = bridges[np.all(ref.bridge_span_direction[bridges] == row, axis=1)]
        pts = mesh_data.vertices[np.unique(mesh_data.faces[faces])]
        pts = pts - pts.mean(axis=0)
        native_rows = np.unique(nat.bridge_span_direction[faces], axis=0)
        assert len(native_rows) == 1
        var_ref, var_nat = np.var(pts @ row), np.var(pts @ native_rows[0])
        assert abs(var_ref - var_nat) <= 1e-9 * max(var_ref, 1e-12)
    if shape == "hollow":
        assert ref.is_bridge.any(), "la cavite doit produire un pont, sinon le test ne verifie rien"


@pytest.mark.parametrize("merge", [None, 30.0])
def test_segmentation_native_matches_python(restore_switch, merge):
    core = _cpp.module()
    if not hasattr(core, "segment_default"):
        pytest.skip("module natif anterieur au pas 21b")
    import os
    import tempfile
    from dataclasses import replace
    import trimesh
    from element_ardor import step2_overhang as s2, step3_segmentation as s3
    from element_ardor.step1_ingest import run as ingest_run

    # T couche + bras pendant : plusieurs regions et dependances
    v = trimesh.creation.box(extents=[4, 4, 20]); v.apply_translation([0, 0, 10])
    h = trimesh.creation.box(extents=[24, 4, 4]); h.apply_translation([0, 0, 22])
    d = trimesh.creation.box(extents=[4, 4, 10]); d.apply_translation([10, 0, 15])
    mesh = v.union(h, engine="manifold").union(d, engine="manifold")
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, tet_mesh, report = ingest_run(path, resolution=2.0)
    finally:
        os.unlink(path)
    om = s2.run(mesh_data, s2.OverhangParams())
    params = replace(s3.SegmentationParams(), platform_z=report.platform_z_used, max_region_merge_angle_deg=merge)
    _cpp.set_enabled(False)
    a = s3.run(mesh_data, om, params)
    _cpp.set_enabled(True)
    b = s3.run(mesh_data, om, params)
    if merge is not None:
        assert len(a.regions) > 1, "la fixture doit produire plusieurs regions"
    assert a.regions.keys() == b.regions.keys()
    assert a.dependencies == b.dependencies
    assert a.region_print_depth == b.region_print_depth
    for k in a.regions:
        assert np.array_equal(a.regions[k].face_ids, b.regions[k].face_ids)
        assert np.array_equal(a.regions[k].target_direction, b.regions[k].target_direction)


def _schedule_signature(colored, schedule):
    layers = [([np.asarray(p).tolist() for p in (l.toolpath_2d or [])], list(l.segment_colors or []), l.color_id)
              for l in colored]
    passes = [(list(p.segments), p.color_id, p.is_hole_fill, tuple(p.layer_range), list(p.deferred_regions))
              for p in schedule.passes]
    return layers, passes


@pytest.mark.parametrize("case", ["none", "two_tone", "lateral", "lateral_merge", "cheap_holes"])
def test_schedule_native_matches_python(restore_switch, case):
    core = _cpp.module()
    if not hasattr(core, "schedule_layers"):
        pytest.skip("module natif anterieur au pas 21c")
    import trimesh
    from element_ardor import step7_scheduling as s7
    from tests.test_step7_scheduling import (
        _lateral_color_map, _pipeline_through_step6, _two_tone_color_map,
    )
    from element_ardor import step6_toolpath2d as s6
    _cpp.set_enabled(True)       # fixtures (kd-tree natif, etc.) construites noyau actif, quel que soit l'environnement
    mesh = trimesh.creation.box(extents=[10, 10, 10]).subdivide().subdivide()
    # Couches de reference par shapely : avec le remplissage natif (Clipper2, arrondi a 1e-6),
    # le milieu d'une arete de la boite symetrique tombe a 7e-8 de la ligne x=0 de la grille de
    # faces, a egale distance de deux triangles de couleurs differentes -- une egalite que les deux
    # implementations departagent a leur facon (tolerance de frontiere declaree, apex/scheduling.hpp).
    mesh_data, tet_mesh, regions, field, layers = _pipeline_through_step6(
        mesh, platform_z=-5.0,
        step6_params=s6.Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.3, native_fill=False))
    params = s7.Step7Params()
    color_map = None
    if case == "two_tone":
        color_map = _two_tone_color_map(mesh_data)
    elif case == "lateral":
        color_map = _lateral_color_map(mesh_data, split_x=0.7)
    elif case == "lateral_merge":
        color_map = _lateral_color_map(mesh_data, split_x=0.7, merge_tolerance=3.0)
    elif case == "cheap_holes":
        color_map = _lateral_color_map(mesh_data, split_x=0.7)
        params = s7.Step7Params(
            hole_cost_model=s7.HoleCostModel(reheat_penalty=0.0, adhesion_risk_weight=0.0),
            purge_cost_model=s7.PurgeCostModel(volume_per_purge=1000.0, time_per_purge=1000.0),
            max_open_holes=100,
        )
    _cpp.set_enabled(False)
    ref = s7.run(layers, mesh_data, tet_mesh, color_map, regions, field, params)
    _cpp.set_enabled(True)
    assert s7._native_schedule(layers, mesh_data, tet_mesh, color_map, field, params) is not None
    nat = s7.run(layers, mesh_data, tet_mesh, color_map, regions, field, params)
    sig_ref, sig_nat = _schedule_signature(*ref), _schedule_signature(*nat)
    assert sig_ref[1] == sig_nat[1]
    assert len(sig_ref[0]) == len(sig_nat[0])
    for (pa, ca, da), (pb, cb, db) in zip(sig_ref[0], sig_nat[0]):
        assert ca == cb and da == db and len(pa) == len(pb)
        for a, b in zip(pa, pb):
            assert np.array_equal(np.asarray(a), np.asarray(b))
    if case == "cheap_holes":
        assert any(p.is_hole_fill for p in nat[1].passes), "la branche trou doit etre exercee"


def test_past_material_batch_matches_single_calls():
    """Requete par lots, en parallele : identique AU BIT aux appels unitaires."""
    core = _cpp.module()
    if not hasattr(core.PastMaterialIndex, "min_past_distance_batch"):
        pytest.skip("module natif anterieur au lot parallele")
    rng = np.random.default_rng(19)
    pts = rng.uniform(0, 20, (5000, 3))
    owners = rng.integers(0, 400, 5000).astype(np.int64)
    index = core.PastMaterialIndex(pts, owners)
    poses = rng.uniform(0, 20, (300, 40, 3))
    limits = rng.integers(0, 400, 300).astype(np.int64)
    batch = index.min_past_distance_batch(poses, limits, 1.5, 4)
    single = np.array([index.min_past_distance(np.ascontiguousarray(poses[k]), int(limits[k]), 1.5)
                       for k in range(300)])
    assert np.array_equal(batch, single)
    assert np.isfinite(batch).sum() > 30, "des poses doivent toucher de la matiere"


@pytest.mark.parametrize("solver", ["direct", "cg"])
@pytest.mark.parametrize("project", [False, True])
def test_heat_method_matches_python(solver, project, restore_switch):
    import dataclasses
    from element_ardor.step4_field_solver import (
        FieldSolverParams, _compute_shape_gradients, _heat_method_distance,
    )

    _mesh_data, tet_mesh, _overhang, _regions = _step4_inputs(resolution=2.0)
    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)
    sources = np.flatnonzero(tet_mesh.tet_vertices[:, 2] < tet_mesh.tet_vertices[:, 2].min() + 1e-6)
    assert len(sources) > 3 and tet_mesh.tet_count > 100

    params = FieldSolverParams(
        convergence_tolerance=1e-12, max_solver_iterations=20000,
        direct_solve_node_threshold=10**9 if solver == "direct" else 1,
        heat_project_reachable=project,
        radial_tilt_origin=np.array([0.5, -0.3, 0.0]) if project else None,
    )
    _cpp.set_enabled(False)
    py = _heat_method_distance(tet_mesh, shape_gradients, volumes, sources, params)
    _cpp.set_enabled(True)
    native = _heat_method_distance(tet_mesh, shape_gradients, volumes, sources, params)

    assert np.isfinite(py).all() and py.max() > 1.0   # garde anti-faux-vert
    assert np.abs(py - native).max() / py.max() < 1e-6


@pytest.mark.parametrize("pivot,retract,speed,temp", [
    (46.0, 0.8, 35.0, 195.0), (45.5, 1.0, 20.0, 210.5), (46.123456789, 0.25, 7.0, 0.0),
    (100.0, 0.001, 1234.5, 250.0),
])
def test_gcode_header_footer_match_python_text(pivot, retract, speed, temp):
    mod = _cpp.module()
    if mod is None or not hasattr(mod, "gcode_header"):
        pytest.skip("module natif absent")
    expected = ["; Apex -- G-code machine R C Z B (file_io/gcode_export.py)",
                f"; pivot->pointe {pivot} mm compense ; retraction {retract} mm a {speed} mm/s",
                "G21", "G90", "M82", "G92 E0", f"M109 S{temp:g}", "; apex_sliced = 1"]
    assert mod.gcode_header(pivot, retract, speed, temp).split("\n") == expected
    assert mod.gcode_footer().split("\n") == ["M104 S0", "M84"]
    assert mod.gcode_header(1e-9, retract, speed, temp) is None   # format diverge : repli Python


def test_reorient_coupled_tiers_native_matches_python(restore_switch):
    import trimesh
    from core.mesh_types import MachineProfile, RotationAxisSpec, ToolpathMachine, ToolpathPoint
    from element_ardor.step9_collision_sim import CollisionSimParams
    from element_ardor.step10_adjustment import reorient_coupled_tiers

    rod = trimesh.creation.box(extents=[1.0, 1.0, 20.0])
    rod.apply_translation([0.0, 0.0, 10.0])
    machine = MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
                       RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))],
        platform_origin=np.array([0.0, 0.0, 0.0]), nozzle_model=rod, nozzle_diameter=0.4,
        position_coupled_axes=[0], radial_direction=np.array([0.0, -1.0, 0.0]))
    sim = CollisionSimParams(path_sampling_step=0.5, safety_margin=1.0, nozzle_contact_clearance=1.5,
                             toolhead_sample_spacing=1.0)
    pts = [ToolpathPoint(position=np.array([9.5, 0.0, 5.0]), axis_angles=[90.0, 0.0], extrusion=1.0, feedrate=1.0),
           ToolpathPoint(position=np.array([10.5, 0.0, 5.0]), axis_angles=[90.0, 0.0], extrusion=1.0, feedrate=1.0)]
    for i, x in enumerate(np.arange(6.0, 14.01, 0.5)):
        pts.append(ToolpathPoint(position=np.array([x, 0.0, 0.2]), axis_angles=[90.0, 0.0],
                                 extrusion=0.0 if i == 0 else 1.0, feedrate=1.0))
    tp = ToolpathMachine(points=pts)
    contacts = [8, 10, 12]

    _cpp.set_enabled(False)
    py_out, py_log = reorient_coupled_tiers(tp, contacts, machine, sim)
    _cpp.set_enabled(True)
    nat_out, nat_log = reorient_coupled_tiers(tp, contacts, machine, sim)

    assert any(b is not None for _, b, _ in py_log)          # garde anti-faux-vert
    assert [(i, t) for i, _, t in py_log] == [(i, t) for i, _, t in nat_log]
    for (_, b0, _), (_, b1, _) in zip(py_log, nat_log):
        assert (b0 is None) == (b1 is None) and (b0 is None or b0 == pytest.approx(b1))
    a_py = np.array([p.axis_angles for p in py_out.points])
    a_nat = np.array([p.axis_angles for p in nat_out.points])
    assert np.abs(a_py - a_nat).max() < 1e-9


def test_clamp_head_height_native_matches_python(restore_switch):
    from tests.test_step10_adjustment import _clamp_machine
    from core.mesh_types import ToolpathMachine, ToolpathPoint
    from element_ardor.step9_collision_sim import CollisionSimParams
    from element_ardor.step10_adjustment import clamp_head_height

    machine = _clamp_machine()
    sim = CollisionSimParams(path_sampling_step=1.0, safety_margin=1.0, nozzle_contact_clearance=0.0)
    rng = np.random.default_rng(7)
    pts = [ToolpathPoint(position=np.array([10.0, 0.0, float(z)]), axis_angles=[90.0, float(b)],
                         extrusion=1.0, feedrate=1.0)
           for z, b in zip(rng.uniform(0.0, 40.0, 300), rng.uniform(0.0, 185.0, 300))]
    tp = ToolpathMachine(points=pts)
    _cpp.set_enabled(False)
    py_out, py_moved, py_dev = clamp_head_height(tp, machine, sim)
    _cpp.set_enabled(True)
    nat_out, nat_moved, nat_dev = clamp_head_height(tp, machine, sim)
    assert py_moved > 5 and py_moved == nat_moved           # garde anti-faux-vert
    assert np.allclose(py_dev, nat_dev, atol=1e-12)
    assert [p.axis_angles for p in py_out.points] == [p.axis_angles for p in nat_out.points]


def test_step10_run_native_matches_python_on_coupled_machine(restore_switch):
    import trimesh
    from core.mesh_types import MachineProfile, RotationAxisSpec, ToolpathMachine, ToolpathPoint
    from element_ardor import step9_collision_sim as s9
    from element_ardor.step10_adjustment import Step10Params, run as step10_run

    head = trimesh.creation.box(extents=[6, 6, 20])
    head.apply_translation([0, 0, 13])
    machine = MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
                       RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))],
        platform_origin=np.array([0.0, 0.0, -50.0]), nozzle_model=head, nozzle_diameter=0.4,
        position_coupled_axes=[0], max_angular_velocity=float(np.radians(1667.0)))
    pts = []
    for k in range(60):
        pts.append(ToolpathPoint(position=np.array([10.0, 0.0, 0.2 * k]), axis_angles=[0.0, 0.0],
                                 extrusion=1.0, feedrate=25.0, color_id=None))
    for k in range(40):
        pts.append(ToolpathPoint(position=np.array([10.0 + 0.25 * k, 2.0, 0.2]), axis_angles=[0.0, 0.0],
                                 extrusion=1.0, feedrate=25.0, color_id=None))
    tp = ToolpathMachine(points=pts)
    sim = s9.CollisionSimParams(path_sampling_step=0.1, safety_margin=1.0, fast_kdtree=True)
    _cpp.set_enabled(False)
    report = s9.run(tp, machine, sim)
    py_out, py_rep = step10_run(tp, report, machine, sim, Step10Params())
    _cpp.set_enabled(True)
    nat_out, nat_rep = step10_run(tp, report, machine, sim, Step10Params())

    assert report.events                                                   # garde anti-faux-vert
    assert any(e.method == "local_reorientation" for e in py_rep.entries)
    assert [(e.point_index, e.method) for e in py_rep.entries] == [(e.point_index, e.method) for e in nat_rep.entries]
    assert py_rep.converged == nat_rep.converged
    a_py = np.array([p.axis_angles for p in py_out.points])
    a_nat = np.array([p.axis_angles for p in nat_out.points])
    assert np.abs(a_py - a_nat).max() < 1e-9


def test_stable_topological_order_matches_python_with_cycles(restore_switch):
    from element_ardor.step5b_collision_order import _stable_topological_order

    rng = np.random.default_rng(3)
    for n, m in [(30, 40), (60, 400), (200, 3000)]:       # le dernier est dense : cycles garantis
        edges = {(int(a), int(b)) for a, b in rng.integers(0, n, size=(m, 2)) if a != b}
        _cpp.set_enabled(False)
        py = _stable_topological_order(n, edges)
        _cpp.set_enabled(True)
        nat = _stable_topological_order(n, edges)
        assert nat == py
    assert py[1] > 0                                       # garde anti-faux-vert : des contraintes violees


def test_step5_default_layers_native_match_python(restore_switch):
    from core.mesh_types import ScalarField
    from element_ardor.step4_field_solver import _compute_shape_gradients, _compute_per_tet_gradient
    from element_ardor.step5_layers import LayerParams, run as layers_run

    mesh_data, tet_mesh, _overhang, regions = _step4_inputs(resolution=2.0)
    shape_gradients, _volumes = _compute_shape_gradients(tet_mesh)
    values = tet_mesh.tet_vertices[:, 2] + 0.3 * tet_mesh.tet_vertices[:, 0]
    field = ScalarField(values=values, gradient=_compute_per_tet_gradient(tet_mesh, shape_gradients, values))
    params = LayerParams(layer_height_base=1.0, isovalue_order=True)

    _cpp.set_enabled(False)
    py = layers_run(tet_mesh, mesh_data, field, regions, params)
    _cpp.set_enabled(True)
    nat = layers_run(tet_mesh, mesh_data, field, regions, params)

    assert len(py) > 5 and len(py) == len(nat)                      # garde anti-faux-vert
    for a, b in zip(py, nat):
        assert a.isovalue == b.isovalue and a.region_ids == b.region_ids
        assert len(a.contour_3d) == len(b.contour_3d)
        for ca, cb in zip(a.contour_3d, b.contour_3d):
            assert np.array_equal(ca, cb)
        for fa, fb in zip(a.local_frame, b.local_frame):
            assert np.allclose(fa, fb, atol=1e-12, rtol=0)


def test_step8_median_filter_and_per_point_normals_match_python(restore_switch):
    from types import SimpleNamespace
    from core.mesh_types import MachineProfile, RotationAxisSpec
    from element_ardor import step8_machine_toolpath as s8

    rng = np.random.default_rng(11)
    n = 400
    axes = [RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
            RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))]
    machine = MachineProfile(rotation_axes=axes, platform_origin=np.zeros(3), nozzle_model=None,
                             nozzle_diameter=0.4, position_coupled_axes=[0])
    desired = rng.normal(size=(n, 3))
    desired /= np.linalg.norm(desired, axis=1)[:, None]
    layer_id = np.repeat(np.arange(8), 50).astype(np.int64)
    segment = (np.arange(n) // 17).astype(np.int64)
    travel = np.zeros(n, bool)
    travel[[0, 60, 61, 200, 399]] = True
    angles = np.column_stack([rng.uniform(0, 360, n), rng.uniform(0, 185, n)])

    def raw():
        return s8._RawPoints(position=rng.normal(size=(n, 3)) * 0 + np.arange(n)[:, None],
                             desired_direction=desired.copy(), layer_id=layer_id, segment_index=segment,
                             color_id=[None] * n, feedrate=np.ones(n), extrusion=np.ones(n), is_travel=travel)

    _cpp.set_enabled(False)
    py_angles, py_res = s8._median_filter_free_axes(angles.tolist(), raw(), axes, machine, np.array([0.0, 0.0, 1.0]), 9)
    _cpp.set_enabled(True)
    nat_angles, nat_res = s8._median_filter_free_axes(angles.tolist(), raw(), axes, machine, np.array([0.0, 0.0, 1.0]), 9)
    assert np.array_equal(np.array(py_angles), np.array(nat_angles))        # une mediane est exacte
    assert np.allclose(py_res, nat_res, atol=1e-12)
    assert not np.array_equal(np.array(py_angles)[:, 1], angles[:, 1])      # le filtre a bien agi

    layers = {i: SimpleNamespace(deposition_reversed=bool(i % 3 == 0)) for i in range(8)}
    field = lambda pos: np.where(np.arange(len(pos))[:, None] % 7 == 0, np.nan, np.roll(desired, 3, axis=0))
    results = []
    for enabled in (False, True):
        _cpp.set_enabled(enabled)
        rp = raw()
        s8._apply_per_point_normals(rp, layers, field)
        results.append(rp.desired_direction.copy())
    assert np.array_equal(results[0], results[1])
    assert not np.array_equal(results[0], desired)                           # des normales ont ete remplacees


@pytest.mark.parametrize("columns", [2, 3])
def test_step8_raw_points_native_identical_to_python(columns, restore_switch):
    from core.mesh_types import Layer, PrintPass, PrintSchedule
    from element_ardor import step8_machine_toolpath as s8

    rng = np.random.default_rng(5)
    layers = []
    for lid in range(6):
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        v = np.cross(u, rng.normal(size=3)); v /= np.linalg.norm(v)
        polys = [rng.normal(size=(int(rng.integers(0, 9)), columns)) for _ in range(5)]
        layer = Layer(id=lid * 3 + 1, contour_3d=[], region_ids=[], local_frame=(rng.normal(size=3), u, v),
                      toolpath_2d=polys)
        layer.deposition_reversed = bool(lid % 2)
        layers.append(layer)
    by_id = {l.id: l for l in layers}
    passes = []
    for color in (None, 0, 1):
        segs = [(layers[int(rng.integers(0, 6))].id, int(rng.integers(0, 5))) for _ in range(12)]
        passes.append(PrintPass(segments=segs, color_id=color, is_hole_fill=False, layer_range=(0, 0), deferred_regions=[]))
    schedule = PrintSchedule(passes=passes)
    params = s8.Step8Params()

    _cpp.set_enabled(False)
    ref = s8._iter_schedule_points(schedule, by_id, params)
    _cpp.set_enabled(True)
    assert s8._iter_schedule_points_native(schedule, by_id, params) is not None
    nat = s8._iter_schedule_points(schedule, by_id, params)

    assert len(ref) == len(nat) > 20 and ref.is_travel.any()                # garde anti-faux-vert
    for name in ("position", "layer_id", "segment_index", "feedrate", "extrusion", "is_travel"):
        assert np.array_equal(getattr(ref, name), getattr(nat, name)), name
    # La normale est normalisee par np.linalg.norm (ddot de BLAS) : quelques ULP d'ecart possibles.
    assert np.allclose(ref.desired_direction, nat.desired_direction, rtol=0, atol=1e-15)
    assert ref.color_id == nat.color_id


def test_step6_native_layers_match_per_layer_path(restore_switch):
    """Couches entieres en natif contre le chemin couche par couche. Sans pont : memes polylignes
    (a l'ULP de la projection BLAS pres). Avec ponts : le chemin couche par couche passe par shapely
    pour eux, le natif par Clipper2 -> meme longueur par couche, decoupage different."""
    from dataclasses import replace
    from element_ardor import step6_toolpath2d as s6
    from element_ardor.step4_field_solver import FieldSolverParams, run as field_run
    from element_ardor.step5_layers import LayerParams, run as layers_run

    _cpp.set_enabled(True)       # le kd-tree de mesh_data doit etre natif pour le chemin par couches entieres
    mesh_data, tet_mesh, overhang, regions = _step4_inputs(resolution=2.0)
    field = field_run(tet_mesh, mesh_data, overhang, regions, FieldSolverParams())
    layers = layers_run(tet_mesh, mesh_data, field, regions, LayerParams(layer_height_base=1.0))
    base = s6.Step6Params(out_of_plane_interpolation_power=2.0, split_rings_into_runs=True,
                          continuous_rings=True, run_split_min_amplitude=0.05)

    _cpp.set_enabled(True)
    assert s6._fill_layers_native(layers, mesh_data, overhang, base, {}) is not None

    def both(params):
        ref = s6.run(layers, mesh_data, overhang, dataclasses.replace(params, native_layer=False))
        return ref, s6.run(layers, mesh_data, overhang, dataclasses.replace(params, native_layer=True))

    def length(polys):
        return sum(float(np.linalg.norm(np.diff(p[:, :2], axis=0), axis=1).sum()) for p in polys)

    ref, nat = both(replace(base, bridge_contour_ratio_threshold=2.0))        # ponts desactives
    assert len(layers) > 5 and sum(len(l.toolpath_2d) for l in ref) > 50      # garde anti-faux-vert
    for a, b in zip(ref, nat):
        assert len(a.toolpath_2d) == len(b.toolpath_2d)
        for pa, pb in zip(a.toolpath_2d, b.toolpath_2d):
            assert pa.shape == pb.shape and np.allclose(pa, pb, atol=1e-9, rtol=0)

    ref, nat = both(base)                                                      # ponts actifs
    differing = 0
    for a, b in zip(ref, nat):
        assert length(a.toolpath_2d) == pytest.approx(length(b.toolpath_2d), rel=2e-3)
        differing += len(a.toolpath_2d) != len(b.toolpath_2d)
    assert differing > 0        # garde : la fixture contient bien des couches en pont


def test_fast_batch_events_match_python_checker(restore_switch):
    import trimesh
    from core.mesh_types import MachineProfile, RotationAxisSpec, ToolpathMachine, ToolpathPoint
    from element_ardor import step9_collision_sim as s9

    head = trimesh.creation.box(extents=[6, 6, 20])
    head.apply_translation([0, 0, 13])
    machine = MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
                       RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))],
        platform_origin=np.array([0.0, 0.0, -0.3]), nozzle_model=head, nozzle_diameter=0.4,
        platform_radius=150.0, platform_thickness=1.0)
    rng = np.random.default_rng(2)
    pts = []
    for k in range(60):
        pts.append(ToolpathPoint(position=np.array([10.0, 0.0, 0.2 * k]), axis_angles=[90.0, float(rng.uniform(0, 40))],
                                 extrusion=1.0, feedrate=1.0, color_id=None))
    for k in range(40):
        pts.append(ToolpathPoint(position=np.array([10.0 + 0.25 * k, 2.0, 0.2]),
                                 axis_angles=[90.0, float(rng.uniform(0, 185))], extrusion=1.0, feedrate=1.0,
                                 color_id=None))
    tp = ToolpathMachine(points=pts)
    params = s9.CollisionSimParams(path_sampling_step=0.1, safety_margin=1.0, fast_kdtree=True)
    samples = s9._filter_nozzle_contact_zone(s9._surface_sample_points(head, params.toolhead_sample_spacing),
                                             params.nozzle_contact_clearance)
    indices = list(range(len(pts)))
    _cpp.set_enabled(False)
    ref = s9._run_fast_kdtree(tp, machine, params, samples, indices)
    _cpp.set_enabled(True)
    checker = s9.FastCollisionChecker(tp, machine, params, samples)
    nat = s9.fast_batch_events(checker, tp, indices)
    assert nat is not None and len(ref) > 5                                  # garde anti-faux-vert
    key = lambda e: (e.point_index, e.obstacle_type)
    assert [key(e) for e in ref] == [key(e) for e in nat]
    assert np.allclose([e.penetration_depth for e in ref], [e.penetration_depth for e in nat], atol=1e-9)
    assert {e.obstacle_type for e in ref} == {"printed_part", "platform"}    # les deux voies sont exercees


def test_native_tetgen_worker_matches_binding():
    """Le maillage de l'executable C++ est identique a celui du binding."""
    import numpy as np
    import trimesh

    from element_ardor import step1_ingest as s1

    if s1._native_tetgen_worker_path() is None:
        pytest.skip("apex_tetgen_worker non construit")
    mesh = trimesh.load("res/test_models/l_shape_overhang.stl")
    mesh.process()
    z = float(mesh.bounds[0][2])
    for near in (None, 0.5):  # sans / avec raffinement plateau (2 regions)
        ref = s1._tetrahedralize_once(mesh, 1.0, near, 3.0, z)
        got = s1._tetrahedralize_once(
            mesh, 1.0, near, 3.0, z, _tetgen_factory=s1._NativeTetgen
        )
        assert len(ref.tetrahedra) > 100  # garde anti-faux-vert
        assert np.array_equal(ref.tet_vertices, got.tet_vertices)
        assert np.array_equal(ref.tetrahedra, got.tetrahedra)


def test_native_grid_fallback_matches_python_grid_search(restore_switch, monkeypatch):
    """
    Repli par grille PORTE en natif (GridSolver) contre la grille Python, sur beaucoup de points
    colineaires a l'axe externe. Seul +Z est une tangence qui met l'analytique en echec (B = 0) ;
    -Z (B = 180, dans la plage [0, 185]) se resout analytiquement : mesure sur la machine reelle,
    43 repli sur 43 pour +Z, 0 sur 43 pour -Z. Les points -Z servent donc de controle : ils ne doivent
    PAS se replier. `Step8Params.native_grid_fallback=False` rend le repli Python, point par point.
    """
    import element_ardor.step8_machine_toolpath as step8

    raw_points, axes, reference, origin, params = _ik_fixture(_current_machine_axes(), seed=303, count=300)
    plus_z = 0
    for index in range(0, 300, 7):
        raw_points[index].desired_direction = np.array([0.0, 0.0, 1.0 if index % 14 == 0 else -1.0])
        plus_z += index % 14 == 0

    _cpp.set_enabled(True)
    python_grid = step8._resolve_all_axis_angles(
        raw_points, axes, reference, origin, dataclasses.replace(params, native_grid_fallback=False))
    python_count = runtime.state().native_fallbacks
    native_grid = step8._resolve_all_axis_angles(raw_points, axes, reference, origin, params)
    native_count = runtime.state().native_fallbacks

    # garde anti-faux-vert : exactement les points +Z se replient (et le repli est bien exerce)
    assert plus_z >= 20 and python_count == plus_z and native_count == python_count, (plus_z, python_count, native_count)
    _assert_ik_equivalent(python_grid, native_grid)


def test_native_grid_fallback_with_three_axes(restore_switch, monkeypatch):
    """Meme controle sur une cinematique a TROIS axes (le repli doit rester generique en N)."""
    from core.mesh_types import RotationAxisSpec
    import element_ardor.step8_machine_toolpath as step8

    axes = _current_machine_axes() + [RotationAxisSpec(axis=np.array([0.0, 1.0, 0.0]), angle_range_deg=(-30.0, 30.0))]
    raw_points, _axes, reference, origin, params = _ik_fixture(axes, seed=404, count=200)
    for index in range(0, 200, 9):
        raw_points[index].desired_direction = np.array([0.0, 0.0, 1.0])
    _cpp.set_enabled(True)
    python_grid = step8._resolve_all_axis_angles(
        raw_points, axes, reference, origin, dataclasses.replace(params, native_grid_fallback=False))
    native_grid = step8._resolve_all_axis_angles(raw_points, axes, reference, origin, params)
    _assert_ik_equivalent(python_grid, native_grid)
