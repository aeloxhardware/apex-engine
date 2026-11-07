"""Suite découpée de `test_cpp_equivalence.py` (R-030 point 1, limite de 800 lignes par fichier) : tests déplacés tels quels."""

from __future__ import annotations
from test_cpp_equivalence import restore_switch  # noqa: F401  (fixtures)

import numpy as np
import pytest
from core import _cpp
import trimesh
import trimesh.ray.ray_triangle
from core.spatial import RayIntersector as _WrappedRayIntersector
from element_ardor.step1_ingest import _find_self_intersecting_faces as _py_find_self_intersecting_faces
from element_ardor.step1_ingest import _triangles_cross as _py_triangles_cross
from test_cpp_equivalence import (  # noqa: F401
    pytestmark,
)


def _hit_rows(locations: np.ndarray, index_ray: np.ndarray) -> np.ndarray:
    """Impacts sous forme canonique (H,4) triee — comparable sans ordre."""
    if len(locations) == 0:
        return np.zeros((0, 4), dtype=float)
    rows = np.column_stack([np.asarray(locations, dtype=float),
                            np.asarray(index_ray, dtype=float)])
    return rows[np.lexsort(rows.T[::-1])]


def _ray_test_mesh() -> trimesh.Trimesh:
    """Maillage synthetique non trivial : une sphere moins une boite, donc des
    faces internes ET externes, et de vraies aretes partagees (le cas que le
    dedoublonnage de `unique_rows` doit collapser)."""
    sphere = trimesh.creation.icosphere(subdivisions=2, radius=10.0)
    box = trimesh.creation.box(extents=(6.0, 6.0, 30.0))
    return sphere.difference(box, engine="manifold")


@pytest.mark.parametrize(
    "direction_mode", ["vertical", "random_unit", "unnormalized"]
)
def test_bvh_intersects_location_matches_trimesh(direction_mode, restore_switch):
    """
    Memes impacts que `RayMeshIntersector.intersects_location`, AU BIT.

    Les trois modes de direction ne sont pas cosmetiques : `vertical` exerce
    le cas `1/0 = inf` du test de tranches, et `unnormalized` verifie que le
    filtre "forward" utilise bien le produit scalaire avec la direction
    BRUTE (comme `util.diagonal_dot`) et non une direction renormalisee.
    """
    mesh = _ray_test_mesh()
    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)

    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )

    rng = np.random.default_rng(1301)
    lo, hi = mesh.bounds
    count = 800
    origins = rng.uniform(lo, hi, size=(count, 3))
    if direction_mode == "vertical":
        directions = np.tile([0.0, 0.0, -1.0], (count, 1))
    elif direction_mode == "random_unit":
        directions = rng.normal(size=(count, 3))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    else:
        directions = rng.normal(size=(count, 3)) * rng.uniform(0.1, 5.0, size=(count, 1))

    ref_loc, ref_ray, _ref_tri = reference.intersects_location(origins, directions)
    nat_loc, nat_ray, _nat_tri = native.intersects_location(
        np.ascontiguousarray(origins), np.ascontiguousarray(directions)
    )

    assert np.array_equal(_hit_rows(ref_loc, ref_ray), _hit_rows(nat_loc, nat_ray))

    # Garde anti-faux-vert : sans elle, deux resultats VIDES seraient egaux.
    assert len(ref_loc) > 100, f"jeu de test trop pauvre ({len(ref_loc)} impacts)"


def test_bvh_handles_misses_and_empty_input(restore_switch):
    """Rayons qui ratent tout, et lot de rayons vide : ni exception, ni forme
    inattendue, des deux cotes."""
    mesh = _ray_test_mesh()
    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )

    far = np.tile(mesh.bounds[1] + 1000.0, (10, 1))
    away = np.tile([0.0, 0.0, 1.0], (10, 1))
    for loc, ray, _tri in (
        reference.intersects_location(far, away),
        native.intersects_location(np.ascontiguousarray(far), np.ascontiguousarray(away)),
    ):
        assert len(loc) == 0 and len(ray) == 0

    empty = np.zeros((0, 3))
    loc, ray, _tri = native.intersects_location(empty, empty)
    assert len(loc) == 0 and len(ray) == 0


def test_bvh_deduplicates_edge_hits_like_trimesh(restore_switch):
    """
    Un impact PILE sur l'arete partagee par deux triangles ne compte qu'UNE
    fois — c'est le role de `grouping.unique_rows` cote trimesh, et la raison
    pour laquelle un comptage par triangle differe d'un comptage par point.

    Cas construit explicitement : deux triangles coplanaires partageant
    l'arete x=0, et un rayon vertical qui tombe exactement dessus. Un
    dedoublonnage manquant donnerait 2 impacts au lieu d'1, et aucune donnee
    aleatoire ne tomberait sur ce cas.
    """
    vertices = np.array(
        [[-1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 1.0, 0.0],
         [1.0, -1.0, 0.0], [1.0, 1.0, 0.0]]
    )
    faces = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]])
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    origin = np.array([[0.0, 0.0, 5.0]])          # pile au-dessus de l'arete x=0
    direction = np.array([[0.0, 0.0, -1.0]])

    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    ref_loc, ref_ray, _ = reference.intersects_location(origin, direction)

    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )
    nat_loc, nat_ray, _ = native.intersects_location(origin, direction)

    assert np.array_equal(_hit_rows(ref_loc, ref_ray), _hit_rows(nat_loc, nat_ray))
    # Le cas doit VRAIMENT toucher l'arete, sinon le test ne verifie rien.
    assert len(ref_loc) == 1, f"le rayon devrait donner 1 impact dedoublonne, pas {len(ref_loc)}"


def test_wrapped_ray_intersector_dispatches_and_agrees(restore_switch):
    """`core.spatial.RayIntersector` : dispatch correct et memes impacts."""
    mesh = _ray_test_mesh()
    rng = np.random.default_rng(1307)
    lo, hi = mesh.bounds
    origins = rng.uniform(lo, hi, size=(400, 3))
    directions = np.tile([0.0, 0.0, -1.0], (400, 1))

    _cpp.set_enabled(False)
    python_side = _WrappedRayIntersector(mesh)
    assert python_side._native is None
    py_loc, py_ray, _ = python_side.intersects_location(origins, directions)

    _cpp.set_enabled(True)
    native_side = _WrappedRayIntersector(mesh)
    assert native_side._native is not None
    nat_loc, nat_ray, _ = native_side.intersects_location(origins, directions)

    assert np.array_equal(_hit_rows(py_loc, py_ray), _hit_rows(nat_loc, nat_ray))
    assert len(py_loc) > 50


def test_wrapped_ray_intersector_rejects_unsupported_options(restore_switch):
    """
    Une option passee au chemin natif doit LEVER, jamais etre ignoree en
    silence. Aucun appelant du pipeline n'en passe aujourd'hui ; si l'un s'y
    met, l'ecart de comportement doit etre bruyant.
    """
    _cpp.set_enabled(True)
    intersector = _WrappedRayIntersector(_ray_test_mesh())
    with pytest.raises(TypeError, match="non supportees"):
        intersector.intersects_location(
            np.zeros((1, 3)), np.tile([0.0, 0.0, -1.0], (1, 1)), multiple_hits=False
        )


def test_bvh_excludes_the_originating_face(restore_switch):
    """
    NON-REGRESSION D'UN BUG REEL, trouve par la mesure sur le pipeline et non
    par des rayons aleatoires.

    `step2` part d'un point SUR une face, decale de 1e-6 vers le bas, et tire
    vers le bas. L'impact avec sa PROPRE face tombe a `dot = -1e-6`, soit
    exactement le seuil du filtre "vers l'avant" — et le calcul flottant reel
    donne -9.999999994736442e-07, strictement SUPERIEUR au seuil. La narrow
    phase de trimesh l'accepterait donc ; ce qui l'ecarte pour une face
    HORIZONTALE, c'est son BROAD PHASE, qui clippe `t` a
    `buffer_dist = 1e-5` : la boite du rayon s'arrete alors juste sous la
    face, dont l'AABB plate ne la recoupe plus.

    Une premiere version du BVH descendait par un test rayon/AABB classique,
    plus permissif : elle rendait cet impact parasite a distance ~1e-6, et
    `step2` tombait de 49 ponts detectes a 0 sur death.stl. Aucun des 3 000
    rayons aleatoires des tests ci-dessus ne l'avait vu, parce qu'aucun ne
    partait exactement d'une face.

    La fixture est une boite creuse : son plafond de cavite est plat et
    horizontal, exactement la geometrie ou le broad phase tranche.
    """
    outer = trimesh.creation.box(extents=(10.0, 10.0, 10.0))
    inner = trimesh.creation.box(extents=(4.0, 4.0, 4.0))
    mesh = outer.difference(inner, engine="manifold")

    # Point sur le plafond de la cavite (z = +2), decale comme le fait step2.
    origins = np.array([[0.0, 0.0, 2.0 - 1e-6]])
    directions = np.array([[0.0, 0.0, -1.0]])

    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    ref_loc, ref_ray, _ = reference.intersects_location(origins, directions)

    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )
    nat_loc, nat_ray, _ = native.intersects_location(
        np.ascontiguousarray(origins), np.ascontiguousarray(directions)
    )

    assert np.array_equal(_hit_rows(ref_loc, ref_ray), _hit_rows(nat_loc, nat_ray))

    # L'assertion qui verrouille le bug : aucun impact a distance quasi nulle.
    # Sans elle, deux implementations egalement fautives passeraient.
    for loc, ray in ((ref_loc, ref_ray), (nat_loc, nat_ray)):
        distances = np.linalg.norm(loc - origins[np.asarray(ray, dtype=np.int64)], axis=1)
        assert (distances > 1e-4).all(), "impact parasite sur la face d'origine"
    # Garde : le rayon doit VRAIMENT traverser la cavite puis le fond, sinon
    # le test ne verifie rien (il passerait sur un resultat vide).
    assert len(ref_loc) == 2, f"attendu 2 impacts (fond de cavite + fond de boite), vu {len(ref_loc)}"


def test_bvh_matches_trimesh_for_rays_starting_on_faces(restore_switch):
    """
    Cas general du precedent : des centaines d'origines posees SUR les faces,
    decalees comme step2. Ici trimesh accepte parfois l'auto-impact (face non
    plate, dont l'AABB recoupe encore la boite du rayon) — l'oracle exige
    seulement que le natif fasse EXACTEMENT le meme choix, pas qu'il applique
    une regle plus propre.
    """
    mesh = _ray_test_mesh()
    origins = mesh.triangles.mean(axis=1)[:300] + np.array([0.0, 0.0, -1e-6])
    directions = np.tile([0.0, 0.0, -1.0], (len(origins), 1))

    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    ref_loc, ref_ray, _ = reference.intersects_location(origins, directions)

    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )
    nat_loc, nat_ray, _ = native.intersects_location(
        np.ascontiguousarray(origins), np.ascontiguousarray(directions)
    )

    assert np.array_equal(_hit_rows(ref_loc, ref_ray), _hit_rows(nat_loc, nat_ray))
    assert len(ref_loc) > 50, f"jeu de test trop pauvre ({len(ref_loc)} impacts)"


def _step4_inputs(resolution: float = 3.0):
    """Pipeline step1->step3 sur une forme en L, assez fine pour exercer
    plusieurs regions (donc plusieurs sources de melange)."""
    import tempfile
    import os

    import trimesh

    from element_ardor.step1_ingest import run as ingest_run
    from element_ardor.step2_overhang import OverhangParams, run as overhang_run
    from element_ardor.step3_segmentation import SegmentationParams, run as segmentation_run

    vertical = trimesh.creation.box(extents=[6, 6, 24])
    vertical.apply_translation([0, 0, 12])
    horizontal = trimesh.creation.box(extents=[24, 6, 6])
    horizontal.apply_translation([9, 0, 24])
    mesh = vertical.union(horizontal, engine="manifold")

    handle, path = tempfile.mkstemp(suffix=".stl")
    os.close(handle)
    try:
        mesh.export(path)
        mesh_data, tet_mesh, *_ = ingest_run(path, resolution=resolution)
    finally:
        os.unlink(path)

    overhang = overhang_run(mesh_data, OverhangParams())
    regions = segmentation_run(mesh_data, overhang, SegmentationParams(platform_z=0.0))
    return mesh_data, tet_mesh, overhang, regions


def test_shape_gradients_match_numpy_to_machine_precision(restore_switch):
    """
    Cofacteurs (C++) contre `np.linalg.inv` (LAPACK) : meme inverse 3x3 par
    deux algorithmes differents. La tolerance porte sur l'ecart RELATIF - les
    gradients atteignent des magnitudes de plusieurs milliers sur un maillage
    fin, et un seuil absolu y serait ininterpretable.
    """
    from element_ardor.step4_field_solver import _compute_shape_gradients

    _mesh_data, tet_mesh, _overhang, _regions = _step4_inputs()

    _cpp.set_enabled(False)
    gradients_py, volumes_py = _compute_shape_gradients(tet_mesh)
    _cpp.set_enabled(True)
    gradients_native, volumes_native = _compute_shape_gradients(tet_mesh)

    scale = max(float(np.abs(gradients_py).max()), 1e-30)
    assert np.abs(gradients_py - gradients_native).max() / scale < 1e-12
    assert np.allclose(volumes_py, volumes_native, rtol=1e-12, atol=0.0)

    # Garde anti-faux-vert : un maillage vide comparerait deux tableaux vides.
    assert tet_mesh.tet_count >= 90
    assert volumes_py.sum() > 0.0


def test_fem_assembly_matches_scipy_coo_path(restore_switch):
    """
    Assemblage DIRECT en CSR contre le chemin `coo_matrix(...).tocsr()`.

    Deux choses doivent coincider, et la seconde est la plus fragile :
      - les COEFFICIENTS, a la precision machine (l'ordre d'accumulation
        differe de `sum_duplicates`, donc pas au bit) ;
      - le MOTIF de sparsite, exactement. Un motif plus large passerait
        inapercu sur les valeurs (les entrees surnumeraires valant zero) tout
        en changeant le cout de toute factorisation en aval.
    """
    from element_ardor.step4_field_solver import (
        _assemble_stiffness_matrix,
        _compute_shape_gradients,
    )

    _mesh_data, tet_mesh, _overhang, _regions = _step4_inputs()

    _cpp.set_enabled(False)
    gradients, volumes = _compute_shape_gradients(tet_mesh)
    expected = _assemble_stiffness_matrix(tet_mesh, gradients, volumes)

    _cpp.set_enabled(True)
    actual = _assemble_stiffness_matrix(tet_mesh, gradients, volumes)

    assert expected.shape == actual.shape
    assert expected.nnz == actual.nnz, "le motif de sparsite doit etre identique"
    assert np.array_equal(expected.indptr, actual.indptr)
    assert np.array_equal(expected.indices, actual.indices)

    scale = max(float(abs(expected).max()), 1e-30)
    assert abs(expected - actual).max() / scale < 1e-13

    assert expected.nnz > 400, f"jeu de test trop pauvre ({expected.nnz} entrees)"


def test_triangles_cross_matches_python_on_crossing_pair(restore_switch):
    """Deux triangles qui se traversent reellement (pas juste un contact)."""
    t0 = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    t1 = np.array([[1.0, 1.0, -1.0], [1.0, -1.0, -1.0], [1.0, 0.0, 1.0]])

    _cpp.set_enabled(False)
    expected = _py_triangles_cross(t0, t1)
    assert expected is True, "fixture invalide : ces deux triangles doivent se croiser"

    _cpp.set_enabled(True)
    actual = bool(_cpp.module().triangles_cross(
        np.ascontiguousarray(t0), np.ascontiguousarray(t1)
    ))
    assert actual == expected


def test_triangles_cross_matches_python_on_edge_contact(restore_switch):
    """Deux triangles qui se touchent SANS se traverser (bord commun, coplanaires)."""
    t0 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    t1 = np.array([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])

    _cpp.set_enabled(False)
    expected = _py_triangles_cross(t0, t1)
    assert expected is False, "fixture invalide : contact simple, pas une intersection"

    _cpp.set_enabled(True)
    actual = bool(_cpp.module().triangles_cross(
        np.ascontiguousarray(t0), np.ascontiguousarray(t1)
    ))
    assert actual == expected


def test_triangles_cross_matches_python_random_pairs(restore_switch):
    """Balayage aleatoire, dont une fraction reelle de paires qui se croisent."""
    rng = np.random.default_rng(4242)
    mismatches = 0
    crossing_count = 0
    for _ in range(400):
        t0 = rng.uniform(-2.0, 2.0, size=(3, 3))
        t1 = rng.uniform(-2.0, 2.0, size=(3, 3))

        _cpp.set_enabled(False)
        expected = _py_triangles_cross(t0, t1)
        _cpp.set_enabled(True)
        actual = bool(_cpp.module().triangles_cross(
            np.ascontiguousarray(t0), np.ascontiguousarray(t1)
        ))

        if expected:
            crossing_count += 1
        if actual != expected:
            mismatches += 1

    assert crossing_count > 0, "aucune paire croisee tiree : jeu d'essai peu discriminant"
    assert mismatches == 0


def _self_intersecting_test_mesh() -> "trimesh.Trimesh":
    """
    Maillage synthetique watertight portant de VRAIES auto-intersections :
    une boite dont une face a ete remplacee par deux triangles fins tordus
    de facon a traverser une paroi voisine (motif similaire au recouvrement
    CSG mesure sur turkish.stl/death.stl - voir CLAUDE.md).
    """
    import trimesh as _trimesh

    box = _trimesh.creation.box(extents=(4.0, 4.0, 4.0))
    spike = _trimesh.Trimesh(
        vertices=np.array([
            [-3.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 3.0],
        ]),
        faces=np.array([[0, 1, 2]]),
        process=False,
    )
    return _trimesh.util.concatenate([box, spike])


def test_find_self_intersecting_faces_matches_python(restore_switch):
    mesh = _self_intersecting_test_mesh()
    triangles = np.ascontiguousarray(mesh.triangles, dtype=np.float64)
    faces = np.ascontiguousarray(mesh.faces, dtype=np.int64)

    _cpp.set_enabled(False)
    expected = _py_find_self_intersecting_faces(mesh)
    assert expected, "fixture invalide : doit contenir de vraies auto-intersections"

    _cpp.set_enabled(True)
    actual_array = _cpp.module().find_self_intersecting_faces(triangles, faces)
    actual = {int(i) for i in actual_array}

    assert actual == expected


def test_find_self_intersecting_faces_matches_python_on_clean_mesh(restore_switch):
    """Un maillage propre (sans auto-intersection) doit rendre un ensemble vide des deux cotes."""
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=5.0)
    triangles = np.ascontiguousarray(mesh.triangles, dtype=np.float64)
    faces = np.ascontiguousarray(mesh.faces, dtype=np.int64)

    _cpp.set_enabled(False)
    expected = _py_find_self_intersecting_faces(mesh)

    _cpp.set_enabled(True)
    actual_array = _cpp.module().find_self_intersecting_faces(triangles, faces)
    actual = {int(i) for i in actual_array}

    assert expected == set()
    assert actual == expected


def test_step1_dispatches_to_native_when_enabled(restore_switch):
    """
    `step1_ingest._find_self_intersecting_faces` (la fonction du pipeline,
    pas la primitive testee ci-dessus) doit rendre le meme resultat que le
    Python pur, que l'interrupteur soit actif ou non - non-regression du
    branchement dans le fichier du pipeline lui-meme.
    """
    from element_ardor import step1_ingest

    mesh = _self_intersecting_test_mesh()

    _cpp.set_enabled(False)
    expected = step1_ingest._find_self_intersecting_faces(mesh)
    assert expected

    _cpp.set_enabled(True)
    actual = step1_ingest._find_self_intersecting_faces(mesh)

    assert actual == expected
