"""TetraNet : protection des arêtes vives (piste C de R-030 point 8) — caractérisation, DESACTIVEE par defaut.

Mesure du 2026-11-07 : sur une boite, la protection ramene le p99 de l'ecart de peau de 0,13 a 0,015 mm a maille
1,5 mm, mais l'angle dièdre minimal chute de 23 a 5 degres a maille 1,0 mm ; sur death et turkish aucun gain.
Ces tests fixent ce qui est vrai (et ce qui ne l'est pas), ils ne declarent pas la piste resolue."""
import numpy as np
import pytest
import trimesh

from core import _cpp


def _supports_features() -> bool:
    if not _cpp.available() or not hasattr(_cpp.module(), "tetranet_mesh"):
        return False
    try:
        box = trimesh.creation.box(extents=[4, 4, 4])
        _cpp.module().tetranet_mesh(np.ascontiguousarray(box.vertices, dtype=np.float64),
                                    np.ascontiguousarray(box.faces, dtype=np.int64), 2.0, 0.5, 0.05, 0.05, 0, 30.0)
        return True
    except TypeError:
        return False


def _boundary_skin_error(box, out) -> tuple[float, float]:
    tets = np.asarray(out["tetrahedra"])
    verts = np.asarray(out["tet_vertices"])
    faces = np.concatenate([tets[:, [0, 1, 2]], tets[:, [0, 1, 3]], tets[:, [0, 2, 3]], tets[:, [1, 2, 3]]])
    _u, inverse, counts = np.unique(np.sort(faces, axis=1), axis=0, return_inverse=True, return_counts=True)
    boundary = faces[counts[inverse] == 1]
    weights = np.random.default_rng(0).dirichlet([1, 1, 1], size=(len(boundary), 3))
    points = np.einsum("fkj,fjd->fkd", weights, verts[boundary]).reshape(-1, 3)
    dist = trimesh.proximity.ProximityQuery(box).on_surface(points)[1]
    return float(np.percentile(dist, 99)), float(dist.max())


@pytest.fixture(scope="module")
def box():
    mesh = trimesh.creation.box(extents=[10.3, 10.3, 10.3])
    mesh.apply_translation([5.15, 5.15, 5.15])
    return mesh


def _mesh(box, cell, angle):
    return _cpp.module().tetranet_mesh(np.ascontiguousarray(box.vertices, dtype=np.float64),
                                       np.ascontiguousarray(box.faces, dtype=np.int64), cell, 0.5, 0.05, 0.05, 0, angle)


def test_feature_protection_is_off_by_default(box):
    if not _supports_features():
        pytest.skip("noyau natif sans la protection des aretes vives")
    default = _cpp.module().tetranet_mesh(np.ascontiguousarray(box.vertices, dtype=np.float64),
                                          np.ascontiguousarray(box.faces, dtype=np.int64), 1.5)
    explicit_off = _mesh(box, 1.5, 0.0)
    assert np.array_equal(np.asarray(default["tet_vertices"]), np.asarray(explicit_off["tet_vertices"]))
    assert np.array_equal(np.asarray(default["tetrahedra"]), np.asarray(explicit_off["tetrahedra"]))


def test_feature_protection_pulls_the_skin_onto_box_edges(box):
    """Garde : sans la correction l'ecart p99 resterait ~0,13 mm. Echoue si l'arete n'est plus protegee."""
    if not _supports_features():
        pytest.skip("noyau natif sans la protection des aretes vives")
    off = _mesh(box, 1.5, 0.0)
    on = _mesh(box, 1.5, 30.0)
    p99_off, _ = _boundary_skin_error(box, off)
    p99_on, _ = _boundary_skin_error(box, on)
    assert p99_off > 0.1
    assert p99_on < 0.05
    assert on["volume_ratio"] > off["volume_ratio"]
    assert on["min_dihedral_deg"] > 20.0          # qualite des tetraedres preservee a cette maille


@pytest.mark.xfail(strict=True, reason="LIMITE CONNUE (2026-11-07) : a maille 1,0 mm la protection fait chuter l'angle "
                                       "dièdre minimal de 23 a 5 degres ; piste C non fiable, a reprendre")
def test_feature_protection_keeps_quality_at_finer_cell(box):
    if not _supports_features():
        pytest.skip("noyau natif sans la protection des aretes vives")
    off = _mesh(box, 1.0, 0.0)
    on = _mesh(box, 1.0, 30.0)
    assert on["min_dihedral_deg"] >= 0.8 * off["min_dihedral_deg"]
