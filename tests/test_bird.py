"""Voie Bird (pipeline/bird) : topologie sans pyvista, deformation, juge, tranchage deforme.

Equivalence avec la recherche mesuree hors suite (2026-10-07, Benchy incline, 24 345
tetraedres) : topologie identique au bit, attributs a 3e-8 rad, deformation a 0,3 um
en mediane (12,6 um max) pour 0,99 mm de deplacement ; tranchage deforme identique au bit.
"""
import numpy as np
import pytest
import trimesh

from pipeline.bird import deformation as bd
from pipeline.bird.deformed_slicing import DeformedSlicingParams, barycentric, boundary_faces, run as deformed_run
from pipeline.bird.fold_judge import crossing_zones


def _two_tets():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float)
    cells = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
    return bd.BirdTet(points=points, cells=cells)


def test_neighbour_pairs_follow_shared_node_counts():
    tet = _two_tets()
    pairs = bd.neighbour_pairs(tet.cells, tet.number_of_points)
    for kind in ("point", "edge", "face"):
        assert pairs[kind].tolist() == [[0, 1]]          # les deux cellules partagent une face
    cells = np.array([[0, 1, 2, 3], [3, 4, 5, 6], [2, 3, 7, 8]])
    pairs = bd.neighbour_pairs(cells, 9)
    assert pairs["point"].tolist() == [[0, 1], [0, 2], [1, 2]]
    assert pairs["edge"].tolist() == [[0, 2]]
    assert pairs["face"].tolist() == []


def test_surface_faces_point_outward():
    tet = _two_tets()
    faces, owners = bd.surface_faces_and_owners(tet)
    assert len(faces) == 6 and sorted(set(owners.tolist())) == [0, 1]
    centroid = tet.points.mean(axis=0)
    tri = tet.points[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    assert (np.einsum("ij,ij->i", normals, tri.mean(axis=1) - centroid) > 0).all()


@pytest.fixture(scope="module")
def overhang_tet():
    """Pilier + bras en porte-a-faux (vrai surplomb a 90 deg)."""
    pillar = trimesh.creation.box(extents=[6, 6, 12])
    pillar.apply_translation([0, 0, 6])
    arm = trimesh.creation.box(extents=[10, 6, 3])
    arm.apply_translation([8, 0, 10.5])
    mesh = pillar.union(arm, engine="manifold")
    return bd.tetrahedralize(mesh.vertices, mesh.faces, max_tet_volume=2.0)


def test_one_pass_keeps_base_pinned_and_moves_the_overhang(overhang_tet):
    tet = overhang_tet
    topo, attrs = bd.calculate_tet_attributes(tet)
    cell_mask, node_mask, z_first = bd.planar_base_mask(tet, topo, 30.0, 0.2)
    assert node_mask.any() and z_first == pytest.approx(9.0, abs=1e-6)      # dessous du bras
    params = {**bd.DEFAULT_PARAMS, "MAX_OVERHANG": 30.0, "ROTATION_MULTIPLIER": 1.0,
              "ROTATE_ITERATIONS": 20, "DEFORM_ITERATIONS": 50}
    new, new_attrs, info = bd.run_one_iteration(tet, topo, attrs, params,
                                                planar_base=(cell_mask, node_mask, tet.points.copy()))
    assert np.isfinite(new.points).all()
    assert np.array_equal(new.points[node_mask], tet.points[node_mask])          # base figee
    moved = np.linalg.norm(new.points - tet.points, axis=1)
    assert moved[~node_mask].max() > 0.1                                         # la deformation agit
    assert 0.0 <= info["overhang_fraction"] <= 1.0


def test_fold_judge_sees_a_self_crossing_surface():
    a = trimesh.creation.box(extents=[4, 4, 4])
    b = trimesh.creation.box(extents=[4, 4, 4])
    b.apply_translation([2, 2, 2])
    crossing = trimesh.util.concatenate([a, b])                  # deux boites qui se traversent
    zones = crossing_zones(crossing)
    assert zones and zones[0][0] > 0
    assert crossing_zones(a) == []


def test_barycentric_round_trip():
    rng = np.random.default_rng(0)
    tets = rng.normal(size=(50, 4, 3))
    w = rng.dirichlet(np.ones(4), size=50)
    points = np.einsum("nv,nvk->nk", w, tets)
    assert np.allclose(barycentric(points, tets), w, atol=1e-9)


def test_identity_deformation_maps_toolpaths_onto_themselves(overhang_tet):
    """Sans deformation, le retour dans la piece est l'identite : les trajectoires
    ramenees sont celles de l'etape 6 relevees en 3D, au decalage pres."""
    tet = overhang_tet
    offset = np.array([1.0, 2.0, 3.0])
    res = deformed_run(tet.cells, tet.points, tet.points.copy(), offset, DeformedSlicingParams())
    assert res.stats["layers"] > 10 and res.stats["cut"] == 0 and res.stats["jumps"] == 0
    zs = [np.vstack(L["contour_3d"])[:, 2] - offset[2] for L in res.layers if L["contour_3d"]]
    assert all(np.ptp(z) < 1e-6 for z in zs)                                     # couches planes
    with_paths = [L for L in res.layers if L["toolpath_3d"]]
    assert with_paths
    for L in with_paths:                       # trajectoires ramenees dans le plan de leur couche
        z = np.vstack(L["toolpath_3d"])[:, 2] - offset[2]
        assert np.abs(z - L["isovalue"]).max() < 1e-6


def test_boundary_faces_match_topology_surface(overhang_tet):
    tet = overhang_tet
    faces_a = {tuple(sorted(f)) for f in boundary_faces(tet.cells, tet.points)}
    faces_b = {tuple(sorted(f)) for f in bd.surface_faces_and_owners(tet)[0]}
    assert faces_a == faces_b


def test_native_least_squares_matches_scipy(overhang_tet):
    """Solveur natif (TRF/LSMR transcrit de scipy, jacobienne float32) contre scipy :
    tolerance declaree, pas le bit (QR et racines par LAPACK cote scipy). Benchy
    (hors suite, 2026-10-07) : ecart median 0,15 um, max 25 um pour 0,99 mm de
    deplacement ; passe 62,6 s -> 9,9 s."""
    from core import _cpp
    if not _cpp.available() or not hasattr(_cpp.module(), "bird_lsq_deform"):
        pytest.skip("noyau natif sans bird_lsq_deform (reconstruire cpp/lib)")
    tet = overhang_tet
    topo, attrs = bd.calculate_tet_attributes(tet)
    cell_mask, node_mask, _z = bd.planar_base_mask(tet, topo, 30.0, 0.2)
    params = {**bd.DEFAULT_PARAMS, "MAX_OVERHANG": 30.0, "ROTATION_MULTIPLIER": 1.0,
              "ROTATE_ITERATIONS": 30, "DEFORM_ITERATIONS": 80}
    base = (cell_mask, node_mask, tet.points.copy())
    results = {}
    for native in (False, True):
        previous = _cpp.set_enabled(native)
        try:
            results[native] = bd.run_one_iteration(tet, topo, attrs, params, planar_base=base)
        finally:
            _cpp.set_enabled(previous)
    (py_tet, _a, py_info), (nat_tet, _b, nat_info) = results[False], results[True]
    assert np.allclose(nat_info["rotation_field"], py_info["rotation_field"], atol=1e-5)
    moved = np.linalg.norm(py_tet.points - tet.points, axis=1)
    gap = np.linalg.norm(nat_tet.points - py_tet.points, axis=1)
    assert moved.max() > 0.1
    assert gap.max() < 0.02 * moved.max()
    assert nat_info["deform_cost"] == pytest.approx(py_info["deform_cost"], rel=1e-3)


def test_self_intersecting_surface_is_repaired_before_tetgen(tmp_path):
    """Deux boites qui se traversent, concatenees sans union : TetGen refuse la
    surface brute (auto-intersections, cas de death.stl et turkish.stl) ; la voie
    Bird doit alors passer par la cascade de reparation de step1."""
    a = trimesh.creation.box(extents=(20, 6, 6))
    b = trimesh.creation.box(extents=(6, 20, 6))
    b.apply_translation((0, 0, 2))
    raw = trimesh.util.concatenate([a, b])
    with pytest.raises(RuntimeError):
        bd.tetrahedralize(raw.vertices, raw.faces)
    path = tmp_path / "cross.stl"
    raw.export(path)
    tet = bd.load_and_tetrahedralize(str(path))
    assert tet.number_of_cells > 0
    p = tet.points[tet.cells]
    volume = np.abs(np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                              np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))).sum() / 6.0
    union = trimesh.boolean.union([a, b], engine="manifold").volume
    assert abs(volume - union) < 0.02 * union


def test_native_least_squares_threads_bit_identical(overhang_tet):
    """Solveur LSQ natif parallele (2026-10-13) : J x par lignes, J^T u par colonnes dans
    l'ordre croissant des lignes, residus et jacobienne par tetraedre -- deformation,
    rotations et cout IDENTIQUES AU BIT a 1 fil et a tous les coeurs."""
    from core import _cpp
    if not _cpp.available() or not hasattr(_cpp.module(), "bird_lsq_set_threads"):
        pytest.skip("noyau natif sans solveur LSQ parallele")
    previous = _cpp.set_enabled(True)
    try:
        tet = overhang_tet
        topo, attrs = bd.calculate_tet_attributes(tet)
        cell_mask, node_mask, _z = bd.planar_base_mask(tet, topo, 30.0, 0.2)
        params = {**bd.DEFAULT_PARAMS, "MAX_OVERHANG": 30.0, "ROTATION_MULTIPLIER": 1.0,
                  "ROTATE_ITERATIONS": 30, "DEFORM_ITERATIONS": 80}
        base = (cell_mask, node_mask, tet.points.copy())
        results = []
        for threads in (1, 0):
            _cpp.module().bird_lsq_set_threads(threads, 1)        # seuil 1 : le chemin parallele est exerce
            results.append(bd.run_one_iteration(tet, topo, attrs, params, planar_base=base))
        (tet_a, _a, info_a), (tet_b, _b, info_b) = results
        assert np.array_equal(info_a["rotation_field"], info_b["rotation_field"])
        assert np.array_equal(tet_a.points, tet_b.points)
        assert info_a["deform_cost"] == info_b["deform_cost"]
        assert np.abs(tet_a.points - tet.points).max() > 0.1          # garde : la passe deforme
    finally:
        _cpp.module().bird_lsq_set_threads(0, 4096)
        _cpp.set_enabled(previous)


def test_native_barycentric_candidates_match_lapack(overhang_tet):
    """Retour deforme -> origine natif (apex_core.tet_barycentric_candidates) contre la
    resolution LAPACK : tolerance declaree, pas le bit (LU transcrite ; Benchy passe 4 :
    score <= 4,5e-9, point ramene <= 1,2e-10 mm). Les trois appelants sont compares :
    candidats du tranchage deforme, `machine_toolpath.locate`, `PrintedSolid.phi_at`."""
    from core import _cpp
    from pipeline.bird import deformed_slicing as ds
    from pipeline.bird.machine_toolpath import locate
    from pipeline.bird.skin.zones import PrintedSolid
    from scipy.spatial import cKDTree
    if not _cpp.available() or not hasattr(_cpp.module(), "tet_barycentric_candidates"):
        pytest.skip("noyau natif sans tet_barycentric_candidates (reconstruire cpp/lib)")
    tet = overhang_tet
    original = tet.points
    shear = np.array([[1.0, 0.0, 0.15], [0.0, 1.0, 0.0], [-0.1, 0.0, 1.0]])
    deformed = original @ shear.T + 0.05 * np.sin(original)
    rng = np.random.default_rng(3)
    t = rng.integers(0, len(tet.cells), 3000)
    w = rng.dirichlet(np.ones(4), len(t))
    points_def = np.einsum("nv,nvk->nk", w, deformed[tet.cells][t]) + rng.normal(0, 0.4, (len(t), 3))
    points_org = np.einsum("nv,nvk->nk", w, original[tet.cells][t]) + rng.normal(0, 0.4, (len(t), 3))
    out = {}
    for native in (False, True):
        previous = _cpp.set_enabled(native)
        try:
            assert (ds._barycentric_native() is not None) == native          # chemin demande reellement pris
            d = ds.DeformedToOriginal(tet.cells, deformed, original, k=12)
            tree = cKDTree(original[tet.cells].mean(axis=1))
            out[native] = (d._candidates(points_def), locate(points_org, original, tet.cells, tree, k=12),
                           PrintedSolid(tet.cells, original, deformed, k=12).phi_at(points_org))
        finally:
            _cpp.set_enabled(previous)
    (c_py, loc_py, phi_py), (c_nat, loc_nat, phi_nat) = out[False], out[True]
    for a, b in zip(c_py, c_nat):
        if a.dtype.kind == "f":
            fin = np.isfinite(a)
            assert np.array_equal(fin, np.isfinite(b))
            assert np.abs(a[fin] - b[fin]).max() < 1e-9
        else:
            assert np.array_equal(a, b)                                      # candidats, validite, comptes
    valid = c_py[3]
    assert valid.all(axis=1).sum() < len(valid) and valid.any()              # cas dedans ET dehors tires
    # locate : ex aequo tranches par la regle explicite (LOCATE_TIE), plus par l'arrondi.
    # Sans elle, 2 points sur 3000 differaient ici (scores egaux a 2e-16 pres).
    assert np.array_equal(loc_py, loc_nat)
    inside = np.isfinite(phi_py)
    assert 0 < inside.sum() < len(phi_py)                                    # garde anti-faux-vert
    assert np.array_equal(inside, np.isfinite(phi_nat))
    assert np.abs(phi_py[inside] - phi_nat[inside]).max() < 1e-9
