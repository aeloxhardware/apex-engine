"""Suite découpée de `test_cpp_equivalence.py` (R-030 point 1, limite de 800 lignes par fichier) : tests déplacés tels quels."""

from __future__ import annotations
from test_cpp_equivalence import restore_switch  # noqa: F401  (fixtures)

import numpy as np
import pytest
from core import _cpp
import trimesh
import trimesh.ray.ray_triangle
from pipeline._self_intersection_resolve import _intersection_segments as _py_intersection_segments
from pipeline._self_intersection_resolve import generalized_winding_number as _py_generalized_winding_number
from pipeline._self_intersection_resolve import triangle_triangle_segment as _py_triangle_triangle_segment
from pipeline.step6_toolpath2d import _interpolate_out_of_plane as _py_interpolate_out_of_plane
from test_cpp_equivalence import (  # noqa: F401
    pytestmark,
)


_SEGMENT_TOLERANCE = 1e-9  # voir _NATIVE_SEGMENT_TOLERANCE, pipeline/_self_intersection_resolve.py


def test_triangle_triangle_segment_matches_python_random_pairs(restore_switch):
    """
    Balayage aleatoire : la primitive isolee (pas le pipeline complet) doit
    concorder a tolerance serree, avec une garde anti-faux-vert sur le
    nombre de paires reellement croisees tirees.
    """
    rng = np.random.default_rng(777)
    crossing_count = 0
    max_err = 0.0
    for _ in range(500):
        t0 = rng.uniform(-5.0, 5.0, size=(3, 3))
        t1 = rng.uniform(-5.0, 5.0, size=(3, 3))

        _cpp.set_enabled(False)
        expected = _py_triangle_triangle_segment(t0, t1)
        _cpp.set_enabled(True)
        actual = _cpp.module().triangle_triangle_segment(
            np.ascontiguousarray(t0), np.ascontiguousarray(t1)
        )

        assert (expected is None) == (actual is None)
        if expected is None:
            continue
        crossing_count += 1
        err = float(np.max(np.abs(np.asarray(expected) - actual)))
        max_err = max(max_err, err)
        assert err < _SEGMENT_TOLERANCE

    assert crossing_count > 0, "aucune paire croisee tiree : jeu d'essai peu discriminant"


def test_triangle_triangle_segment_degenerate_and_parallel_cases(restore_switch):
    """Cas limites explicites : triangle degenere, plans paralleles, contact ponctuel."""
    degenerate = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    normal_tri = np.array([[0.0, 0.0, -1.0], [1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]])

    _cpp.set_enabled(False)
    expected = _py_triangle_triangle_segment(degenerate, normal_tri)
    assert expected is None, "fixture invalide : triangle degenere attendu"
    _cpp.set_enabled(True)
    actual = _cpp.module().triangle_triangle_segment(
        np.ascontiguousarray(degenerate), np.ascontiguousarray(normal_tri)
    )
    assert actual is None

    coplanar_a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    coplanar_b = np.array([[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]])
    _cpp.set_enabled(False)
    expected = _py_triangle_triangle_segment(coplanar_a, coplanar_b)
    assert expected is None, "fixture invalide : plans paralleles attendus"
    _cpp.set_enabled(True)
    actual = _cpp.module().triangle_triangle_segment(
        np.ascontiguousarray(coplanar_a), np.ascontiguousarray(coplanar_b)
    )
    assert actual is None


def test_generalized_winding_number_matches_python(restore_switch):
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=5.0)
    triangles = np.ascontiguousarray(mesh.triangles, dtype=np.float64)

    rng = np.random.default_rng(42)
    points = rng.uniform(-8.0, 8.0, size=(300, 3))

    _cpp.set_enabled(False)
    expected = _py_generalized_winding_number(points, triangles)
    _cpp.set_enabled(True)
    actual = _cpp.module().generalized_winding_number(
        np.ascontiguousarray(points, dtype=np.float64), triangles
    )

    max_err = float(np.max(np.abs(expected - actual)))
    assert max_err < 1e-10, f"ecart inattendu : {max_err}"
    # Au moins un point clairement interieur (winding ~1) et un exterieur
    # (winding ~0) : garde anti-faux-vert contre un jeu d'essai degenere.
    assert np.any(expected > 0.9)
    assert np.any(expected < 0.1)


def _self_intersecting_turkish_style_mesh() -> "trimesh.Trimesh":
    """Deux nervures qui se traversent - meme motif que turkish.stl (docstring du module)."""
    box = trimesh.creation.box(extents=(4.0, 4.0, 4.0))
    spike = trimesh.Trimesh(
        vertices=np.array([
            [-3.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 3.0],
        ]),
        faces=np.array([[0, 1, 2]]),
        process=False,
    )
    return trimesh.util.concatenate([box, spike])


def test_intersection_segments_same_faces_and_counts(restore_switch):
    """
    `find_intersection_segments` (broad+narrow phase complet) doit affecter
    exactement les MEMES faces, avec le MEME nombre de segments par face,
    que la version Python - seul l'ORDRE a l'interieur d'une face peut
    differer (limite documentee, cpp/include/apex/mesh_intersect.hpp).
    """
    mesh = _self_intersecting_turkish_style_mesh()
    scale = float(mesh.scale) or 1.0
    min_length = 1e-7 * scale

    _cpp.set_enabled(False)
    expected = _py_intersection_segments(mesh, min_length)
    assert expected, "fixture invalide : doit produire des segments d'intersection"

    _cpp.set_enabled(True)
    actual = _py_intersection_segments(mesh, min_length)

    assert set(expected.keys()) == set(actual.keys())
    for face_index in expected:
        assert len(expected[face_index]) == len(actual[face_index])


def test_resolve_self_intersections_agrees_within_documented_tolerance(restore_switch):
    """
    Oracle bout-en-bout sur le vrai cas d'usage (motif turkish.stl) :
    resultat topologiquement PROCHE mais pas garanti identique - la
    divergence assumee est documentee dans mesh_intersect.hpp. On verifie
    ici que l'ecart mesure reste dans l'ordre de grandeur observe (volume
    relatif < 1e-3), pas une exactitude bit-a-bit.
    """
    from pipeline._self_intersection_resolve import resolve_self_intersections

    mesh = _self_intersecting_turkish_style_mesh()

    _cpp.set_enabled(False)
    expected = resolve_self_intersections(mesh)
    _cpp.set_enabled(True)
    actual = resolve_self_intersections(mesh)

    assert (expected is None) == (actual is None)
    if expected is None:
        return
    relative_volume_diff = abs(float(expected.volume) - float(actual.volume)) / max(
        abs(float(expected.volume)), 1e-9
    )
    assert relative_volume_diff < 1e-3
    assert expected.is_watertight == actual.is_watertight


def test_interpolate_out_of_plane_matches_python(restore_switch):
    """Balayage aleatoire, tailles variees, dont des points coincidents forces."""
    rng = np.random.default_rng(11)
    max_err = 0.0
    coincident_trials = 0
    for trial in range(120):
        k = int(rng.integers(1, 60))
        m = int(rng.integers(1, 40))
        query = rng.uniform(-10.0, 10.0, size=(k, 2))
        sources = rng.uniform(-10.0, 10.0, size=(m, 2))
        offsets = rng.uniform(-2.0, 2.0, size=m)
        if trial % 7 == 0:
            # Coincidence EXACTE : exerce la branche `_COINCIDENT_EPS`, qui
            # doit rendre la valeur source telle quelle (jamais 0**-power).
            query[0] = sources[int(rng.integers(0, m))]
            coincident_trials += 1

        _cpp.set_enabled(False)
        expected = _py_interpolate_out_of_plane(query, sources, offsets, 2.0)
        _cpp.set_enabled(True)
        actual = _py_interpolate_out_of_plane(query, sources, offsets, 2.0)

        max_err = max(max_err, float(np.max(np.abs(expected - actual))))

    assert coincident_trials > 0, "aucun cas coincident tire : branche non exercee"
    assert max_err < 1e-12, f"ecart inattendu : {max_err}"


def test_interpolate_out_of_plane_empty_inputs(restore_switch):
    """Les deux cas vides du Python (aucune requete / aucune source)."""
    sources = np.array([[0.0, 0.0], [1.0, 1.0]])
    offsets = np.array([0.5, -0.5])

    for native in (False, True):
        _cpp.set_enabled(native)
        assert len(_py_interpolate_out_of_plane(np.zeros((0, 2)), sources, offsets, 2.0)) == 0
        empty_source = _py_interpolate_out_of_plane(
            np.array([[3.0, 3.0]]), np.zeros((0, 2)), np.zeros(0), 2.0
        )
        assert np.array_equal(empty_source, np.zeros(1))


def test_interpolate_out_of_plane_returns_source_value_on_exact_hit(restore_switch):
    """Une requete posee EXACTEMENT sur une source rend son offset, pas une moyenne."""
    sources = np.array([[0.0, 0.0], [10.0, 0.0]])
    offsets = np.array([1.0, -1.0])
    query = np.array([[0.0, 0.0]])

    for native in (False, True):
        _cpp.set_enabled(native)
        result = _py_interpolate_out_of_plane(query, sources, offsets, 2.0)
        assert result[0] == pytest.approx(1.0, abs=1e-12)


def _mesh_distance_or_skip():
    from core import _cpp
    mod = _cpp.module() if _cpp.available() else None
    if mod is None or not hasattr(mod, "MeshDistance"):
        pytest.skip("apex_core sans MeshDistance (module non reconstruit)")
    return mod


@pytest.mark.parametrize("backend", [0, 1])
def test_mesh_distance_matches_trimesh(backend):
    mod = _mesh_distance_or_skip()
    mesh = trimesh.creation.icosphere(subdivisions=3, radius=10.0)
    mesh.apply_transform(trimesh.transformations.scale_and_translate([1.0, 0.6, 1.3], [1.0, -2.0, 0.5]))
    rng = np.random.default_rng(3)
    pts = rng.uniform(-20, 20, (5000, 3))
    pts[:200] = mesh.vertices[:200] + 1e-3 * rng.normal(size=(200, 3))     # tres pres de la surface
    ref_q, ref_d, _ = trimesh.proximity.closest_point(mesh, pts)
    index = mod.MeshDistance(np.asarray(mesh.vertices, dtype=float), np.asarray(mesh.faces, dtype=np.int64))
    d, q = index.query(pts, backend=backend, threads=4)
    # trimesh n'est PAS exact tout pres des sommets (mesure : jusqu'a 1,9e-5 mm
    # d'erreur sur 15 points sur 5000, toujours par exces). La ou les deux
    # different, l'arbitre est la force brute sur tous les triangles.
    assert np.max(np.abs(d - ref_d)) < 1e-4
    assert np.all(d <= ref_d + 1e-12)
    bad = np.flatnonzero(np.abs(d - ref_d) > 1e-9)
    tri = mesh.triangles
    for i in bad:
        brute = np.linalg.norm(trimesh.triangles.closest_point(tri, np.repeat(pts[i][None], len(tri), 0)) - pts[i],
                               axis=1).min()
        assert abs(d[i] - brute) < 1e-12
    np.testing.assert_allclose(np.linalg.norm(q - pts, axis=1), d, atol=1e-12)
    good = np.abs(d - ref_d) <= 1e-9
    assert np.max(np.linalg.norm(q - ref_q, axis=1)[good & (ref_d > 1e-6)]) < 1e-6


def test_mesh_distance_vulkan_falls_back_to_cpu():
    mod = _mesh_distance_or_skip()
    backend, reason = mod.mesh_distance_backend(2)
    assert backend in (1, 2)
    if backend == 1:
        assert reason


def _random_tet_field(seed: int, n: int = 600):
    from scipy.spatial import Delaunay
    rng = np.random.default_rng(seed)
    pts = rng.uniform(-10.0, 10.0, size=(n, 3))
    tets = Delaunay(pts).simplices.astype(np.int64)
    values = 0.3 * pts[:, 0] + np.sin(pts[:, 1] * 0.4) + 0.05 * pts[:, 2] ** 2 + rng.normal(0, 0.05, n)
    return pts, tets, values


def _python_extract(pts, tets, values, iso, tol):
    from core.mesh_types import TetMesh
    from pipeline.step5_layers import _extract_isosurface_triangles, _weld_vertices
    tm = TetMesh(tet_vertices=pts, tetrahedra=tets, tet_adjacency=None)
    per = values[tets]
    v, f, src = _extract_isosurface_triangles(tm, values, iso, per.min(axis=1), per.max(axis=1))
    wv, wf = _weld_vertices(v, f, tol)
    return wv, wf, src


def test_isosurface_extractor_bit_exact_vs_python():
    core = _cpp.module()
    if not hasattr(core, "IsosurfaceExtractor"):
        pytest.skip("module natif anterieur au pas 13")
    total_caps = total_quads = merged = 0
    for seed in range(4):
        pts, tets, values = _random_tet_field(seed)
        ex = core.IsosurfaceExtractor(pts, tets, values)
        isos = list(np.quantile(values, [0.1, 0.3, 0.5, 0.7, 0.9])) + [float(values[7])]  # iso sur un noeud exact
        for iso in isos:
            for tol in (0.0, 1e-6, 0.05):
                wv, wf, src = _python_extract(pts, tets, values, float(iso), tol)
                nv, nf, nsrc = ex.extract(float(iso), tol)
                assert np.array_equal(src, nsrc)
                assert np.array_equal(wf, nf)
                assert np.array_equal(wv, nv)          # au bit : aucune reduction
                if tol > 0:
                    merged += 3 * len(wf) - len(wv)
            per = values[tets]
            above = (per >= iso).sum(axis=1)
            cand = (per.min(axis=1) <= iso) & (iso <= per.max(axis=1))
            total_caps += int(((above == 1) | (above == 3))[cand].sum())
            total_quads += int((above == 2)[cand].sum())
    # garde anti-faux-vert : les deux cas geometriques et la soudure ont bien ete exerces
    assert total_caps > 100 and total_quads > 100 and merged > 100


def test_isosurface_extractor_empty_and_validation():
    core = _cpp.module()
    if not hasattr(core, "IsosurfaceExtractor"):
        pytest.skip("module natif anterieur au pas 13")
    pts, tets, values = _random_tet_field(9, 200)
    ex = core.IsosurfaceExtractor(pts, tets, values)
    v, f, s = ex.extract(float(values.max()) + 1.0, 1e-6)
    assert v.shape == (0, 3) and f.shape == (0, 3) and s.shape == (0,)
    with pytest.raises(ValueError):
        core.IsosurfaceExtractor(pts, np.array([[0, 1, 2, 10**9]]), values)


def _outline_python(vertices, faces):
    previous = _cpp.set_enabled(False)
    try:
        from pipeline.step5_layers import _extract_closed_contours
        return _extract_closed_contours(vertices, faces)
    finally:
        _cpp.set_enabled(previous)


def test_boundary_loops_match_trimesh_outline(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "boundary_loops") or not hasattr(core, "IsosurfaceExtractor"):
        pytest.skip("module natif anterieur au pas 14")
    import pipeline.step5_layers as s5
    loops_seen = native_used = 0
    for seed in range(6):
        pts, tets, values = _random_tet_field(seed, 800)
        ex = core.IsosurfaceExtractor(pts, tets, values)
        for iso in np.quantile(values, [0.15, 0.35, 0.5, 0.65, 0.85]):
            v, f, _ = ex.extract(float(iso), 1e-6)
            expected = _outline_python(v, f)
            _cpp.set_enabled(True)
            before = s5._NATIVE_OUTLINE_FALLBACKS
            got = s5._extract_closed_contours(v, f)
            native_used += int(s5._NATIVE_OUTLINE_FALLBACKS == before)
            assert len(got) == len(expected)
            for a, b in zip(got, expected):
                assert np.array_equal(a, b)
            loops_seen += len(expected)
    assert loops_seen > 20 and native_used > 20   # garde anti-faux-vert


def test_boundary_loops_fall_back_on_non_manifold_border():
    core = _cpp.module()
    if not hasattr(core, "boundary_loops"):
        pytest.skip("module natif anterieur au pas 14")
    bowtie = np.array([[0, 1, 2], [0, 3, 4]], dtype=np.int64)   # sommet 0 de degre 4
    assert core.boundary_loops(bowtie, 5) is None
    square = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    loops = core.boundary_loops(square, 4)
    assert [list(l) for l in loops] == [[0, 1, 2, 3, 0]]
    degenerate = np.array([[0, 0, 1]], dtype=np.int64)
    assert core.boundary_loops(degenerate, 2) == []


def test_coupled_ik_native_matches_python(restore_switch):
    import pipeline.step8_machine_toolpath as s8
    from core.mesh_types import MachineProfile, RotationAxisSpec
    machine = MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
                       RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))],
        platform_origin=np.zeros(3), position_coupled_axes=[0], radial_direction=np.array([0.0, -1.0, 0.0]))
    rng = np.random.default_rng(17)
    pos = rng.uniform(-40, 40, (600, 3))
    pos[::50, :2] = 0.0                                     # points SUR l'axe : heritage de C
    dirs = rng.normal(size=(600, 3))
    dirs[:, 2] = np.abs(dirs[:, 2]) + 0.3
    dirs[::7] = [0.0, 0.0, 1.0]                             # depots verticaux : repli grille du natif
    raw = [s8._RawPoint(position=p, desired_direction=d / np.linalg.norm(d), layer_id=0, segment_index=0,
                        color_id=None, feedrate=1.0, extrusion=1.0, is_travel=False) for p, d in zip(pos, dirs)]
    ref = np.array([0.0, 0.0, 1.0])
    params = s8.Step8Params()
    _cpp.set_enabled(False)
    a_py, r_py = s8._resolve_all_axis_angles_coupled(raw, machine.rotation_axes, ref, machine, params)
    _cpp.set_enabled(True)
    a_cc, r_cc = s8._resolve_all_axis_angles_coupled(raw, machine.rotation_axes, ref, machine, params)
    a_py, a_cc = np.array(a_py), np.array(a_cc)
    assert np.abs(a_py - a_cc).max() < 1e-7            # degres
    assert np.abs(np.array(r_py) - np.array(r_cc)).max() < 1e-6   # radians (acos mal conditionnee pres de 0)
    # memes diagnostics a 5 deg : l'assertion qui compte reellement (pas 2)
    lim = np.radians(params.max_unreachable_angular_error_deg)
    assert np.array_equal(np.array(r_py) > lim, np.array(r_cc) > lim)
    assert (np.array(r_py) > lim).sum() > 50          # garde : des cibles inatteignables existent


def test_point_in_polygon_2d_native_matches_python(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "point_in_polygon_2d"):
        pytest.skip("module natif anterieur au pas 12")
    from core.geometry_utils import point_in_polygon_2d
    rng = np.random.default_rng(12)
    inside_count = total = 0
    for _ in range(60):
        n = int(rng.integers(6, 40))
        ang = np.sort(rng.uniform(0, 2 * np.pi, n))
        rad = rng.uniform(1.0, 5.0, n)
        poly = np.column_stack([rad * np.cos(ang), rad * np.sin(ang)])
        pts = rng.uniform(-6, 6, (80, 2))
        pts[:5] = poly[:5]                          # sommets exacts : cas limites
        pts[5:8, 1] = poly[:3, 1]                   # meme ordonnee qu'un sommet
        for p in pts:
            _cpp.set_enabled(False)
            a = point_in_polygon_2d(p, poly)
            _cpp.set_enabled(True)
            b = point_in_polygon_2d(p, poly)
            assert a == b
            inside_count += a
            total += 1
    assert 0.1 * total < inside_count < 0.9 * total   # garde : les deux issues exercees


def test_past_material_index_matches_python_fast_checker(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "PastMaterialIndex"):
        pytest.skip("module natif anterieur au pas 19")
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(19)
    dense = rng.uniform(0, 30, (20000, 3))
    owner = rng.integers(0, 5000, 20000).astype(np.int64)
    idx = core.PastMaterialIndex(dense, owner)
    tree = cKDTree(dense)
    owner_pad = np.r_[owner, np.iinfo(np.int64).max]
    reach = 1.2
    hits = 0
    for trial in range(300):
        verts = rng.uniform(0, 30, (60, 3))
        limit = int(rng.integers(0, 5000))
        # reference : algorithme Python de FastCollisionChecker
        dist, ii = tree.query(verts, k=16, distance_upper_bound=reach)
        past = owner_pad[ii] < limit
        d_past = np.where(past, dist, np.inf).min(axis=1)
        sat = np.isfinite(dist[:, -1]) & ~np.isfinite(d_past)
        for vi in np.flatnonzero(sat):
            near = np.asarray(tree.query_ball_point(verts[vi], reach), dtype=int)
            near = near[owner[near] < limit]
            if len(near):
                d_past[vi] = float(np.linalg.norm(dense[near] - verts[vi], axis=1).min())
        ref = float(d_past.min())
        got = idx.min_past_distance(verts, limit, reach)
        assert (ref == got) or (not np.isfinite(ref) and not np.isfinite(got))
        hits += np.isfinite(ref)
    assert 50 < hits < 300          # garde : cas avec et sans matiere passee


def _random_machine_toolpath(seed, n=3000):
    rng = np.random.default_rng(seed)
    steps = rng.normal(0, 0.4, (n, 3))
    steps[:, 2] = np.abs(steps[:, 2]) * 0.05
    pos = np.cumsum(steps, axis=0) + [20.0, 5.0, 0.1]
    pos[100:110, :2] = 0.0                                  # passages au centre (Y herite)
    jumps = rng.choice(n, 60, replace=False)
    pos[jumps] += rng.uniform(-15, 15, (60, 3)) * [1, 1, 0.1]   # trajets longs
    pos[:, 2] = np.maximum(pos[:, 2], 0.05)
    c = (np.degrees(np.arctan2(pos[:, 1], pos[:, 0])) + 90.0) % 360.0
    b = np.clip(rng.normal(20, 30, n), 0, 185)
    ext = (rng.uniform(size=n) > 0.08).astype(float)
    ext[0] = 0
    skin = rng.uniform(size=n) > 0.7
    return {"position": pos, "axis_angles": np.column_stack([c, b]), "extrusion": ext, "is_skin": skin}


@pytest.mark.parametrize("variant", ["defaut", "pivot", "sans_table", "decoupage"])
def test_gcode_export_native_identical_text(restore_switch, variant):
    core = _cpp.module()
    if not hasattr(core, "export_gcode_body"):
        pytest.skip("module natif anterieur au pas 24")
    import main
    from file_io.gcode_export import GcodeParams, export_gcode, head_zlow_from_table
    machine = main._default_machine_profile()
    machine.nozzle_model = None
    if variant == "decoupage":       # bisection max_tip_deviation en natif (2026-10-07)
        params = GcodeParams()
    elif variant == "pivot":
        params = GcodeParams(z_reference="pivot", tilt_sign=-1.0, max_tip_deviation=None)
    else:
        params = GcodeParams(max_tip_deviation=None)
    hz = None if variant == "sans_table" else head_zlow_from_table(np.arange(0, 185.5, 0.5), -np.linspace(1, 30, 371))
    for seed in range(3):
        tp = _random_machine_toolpath(seed)
        _cpp.set_enabled(False)
        a, ra = export_gcode(tp, machine, params, head_zlow=hz)
        _cpp.set_enabled(True)
        b, rb = export_gcode(tp, machine, params, head_zlow=hz)
        assert a == b
        for k in ("moves", "travels", "retractions", "clearance_axes", "filament_mm", "c_mismatch_deg",
                  "tip_z_min", "head_z_min", "tip_below_plate", "head_below_margin", "fw_x_range",
                  "fw_y_range", "fw_z_range", "print_time_s", "subdivided_moves", "max_tip_deviation_before"):
            assert getattr(ra, k) == pytest.approx(getattr(rb, k), rel=1e-13, abs=0.0), k
        if variant == "decoupage":
            assert ra.subdivided_moves > 50                                  # garde : le decoupage agit
        assert ra.retractions > 5 and ra.clearance_axes["R"] + ra.clearance_axes["Z"] > 5   # garde


def _overhang_travel_toolpath():
    """Cordons en porte-a-faux a z = 5 puis trajets a vide SOUS eux a z = 1-1.5,
    tete verticale (B = 0) : l'axe de la tete traverse la matiere deja deposee,
    les degagements sont donc reellement bloques (un cordon court laisse R
    s'echapper, un long bloque tout)."""
    pts, ext = [], []

    def line(a, b, step=0.5):
        a, b = np.asarray(a, float), np.asarray(b, float)
        m = max(1, int(np.ceil(np.linalg.norm(b - a) / step)))
        for k in range(m + 1):
            pts.append(a + (b - a) * k / m)
            ext.append(0.0 if k == 0 else 1.0)

    line([20, 0, 5], [30, 0, 5])            # long porte-a-faux le long de +R
    line([25, 2, 1], [25, 4, 1])            # sous le porte-a-faux mais a cote
    line([24, 0, 1.5], [24, -2, 1.5])       # trajet vers un point SOUS le cordon
    line([12, 0, 5], [12.6, 0, 5])          # porte-a-faux court
    line([12, 0, 1], [12, 3, 1])            # sous lui : R s'echappe
    line([27, 0, 1], [27, -3, 1])           # de nouveau sous le long cordon
    pos = np.array(pts)
    c = (np.degrees(np.arctan2(pos[:, 1], pos[:, 0])) + 90.0) % 360.0
    return {"position": pos, "axis_angles": np.column_stack([c, np.zeros(len(pos))]),
            "extrusion": np.array(ext), "is_skin": np.zeros(len(pos), dtype=bool)}


def test_gcode_export_native_travel_check_identical_text(restore_switch):
    """Controle des degagements (`material_travel_check`) porte en DONNEES dans
    l'export natif (2026-10-13) : texte identique ligne a ligne et memes contacts
    restants que le chemin Python (rappel scipy de la chaine de recherche)."""
    core = _cpp.module()
    if not getattr(core, "EXPORT_TRAVEL_CHECK", False):
        pytest.skip("module natif anterieur au controle des degagements")
    import main
    from file_io.gcode_export import GcodeParams, export_gcode, head_zlow_from_table, material_travel_check
    machine = main._default_machine_profile()
    machine.nozzle_model = None
    hz = head_zlow_from_table(np.arange(0, 185.5, 0.5), -np.linspace(1, 30, 371))
    changed = hits = 0
    for tp in [_overhang_travel_toolpath()] + [_random_machine_toolpath(seed) for seed in range(2)]:
        for params in (GcodeParams(max_tip_deviation=None), GcodeParams()):
            check = material_travel_check(tp, machine, params)
            _cpp.set_enabled(False)
            a, ra = export_gcode(tp, machine, params, head_zlow=hz, travel_check=check)
            plain, rp = export_gcode(tp, machine, params, head_zlow=hz)
            _cpp.set_enabled(True)
            b, rb = export_gcode(tp, machine, params, head_zlow=hz, travel_check=check)
            assert a == b
            for k in ("moves", "travels", "retractions", "clearance_axes", "travel_collisions", "filament_mm",
                      "tip_z_min", "head_z_min", "fw_x_range", "fw_y_range", "fw_z_range", "subdivided_moves"):
                assert getattr(ra, k) == pytest.approx(getattr(rb, k), rel=1e-13, abs=0.0), k
            changed += int(a != plain)
            hits += len(check.hits)
    assert hits > 0 and changed >= 2             # garde : contacts vus, degagements reellement deplaces


def test_feed_plan_native_bit_exact(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "feed_plan_caps"):
        pytest.skip("module natif anterieur au pas 24")
    from file_io.gcode_export import feed_plan
    tp = _random_machine_toolpath(7, 20000)
    _cpp.set_enabled(False)
    a = feed_plan(tp["position"], tp["axis_angles"], 750.0, 300.0, [0], 25.0)
    _cpp.set_enabled(True)
    b = feed_plan(tp["position"], tp["axis_angles"], 750.0, 300.0, [0], 25.0)
    assert a[0] == b[0] and a[1] == b[1] and np.array_equal(a[2], b[2])
    assert (a[2] < 1.0).sum() > 100     # garde : le ralentissement est bien exerce
