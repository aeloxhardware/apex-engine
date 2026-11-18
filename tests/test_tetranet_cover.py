"""TetraNet : ecart surface -> bord du maillage (matiere manquante, coins comblés), mesure le 2026-11-17.

Avant, TetraNet ne mesurait que l'ecart du bord du maillage a la surface : une paroi plus mince que la maille, absente du
maillage, ne se voyait pas dans son statut (Benchy : jusqu'a 1,3 mm de surface loin du bord). La mesure dans l'autre sens
ne compte que la peau EXTERIEURE du solide : les faces internes d'un STL a coques qui se chevauchent (turkish : 17 corps)
n'ont pas a etre au bord (enroulement > 1/2 au-dessus d'elles)."""
import numpy as np
import pytest
import trimesh

from core import _cpp


def _native():
    if not _cpp.available() or not hasattr(_cpp.module(), "tetranet_mesh"):
        pytest.skip("noyau natif TetraNet absent")
    module = _cpp.module()
    if "check_cover" not in (module.tetranet_mesh.__doc__ or ""):
        pytest.skip("noyau natif plus ancien que la mesure de couverture : scripts/build_cpp.py --dev")
    return module


def _mesh(mesh: trimesh.Trimesh, cell: float, options: dict | None = None) -> dict:
    return _native().tetranet_mesh(np.ascontiguousarray(mesh.vertices, dtype=np.float64),
                                   np.ascontiguousarray(mesh.faces, dtype=np.int64), cell, options or {})


def _two_overlapping_boxes() -> trimesh.Trimesh:
    """Deux coques fermees qui se chevauchent SANS union booleenne : chacune a des faces dans l'autre."""
    a = trimesh.creation.box(extents=(10, 10, 10))
    b = trimesh.creation.box(extents=(10, 10, 10))
    b.apply_translation((5.0, 0.0, 0.0))
    return trimesh.util.concatenate([a, b])


def test_internal_faces_of_overlapping_shells_are_not_counted_as_missing_material():
    out = _mesh(_two_overlapping_boxes(), 1.0)
    assert out["cover_internal_samples"] > 0
    # les faces internes sont jusqu'a 5 mm a l'interieur : comptees, elles feraient un ecart de plusieurs mm
    assert out["skin_cover_max_error"] < 0.5


def test_a_wall_thinner_than_the_cell_is_reported():
    plate = trimesh.creation.box(extents=(20, 20, 0.3))     # paroi de 0,3 mm, maille de 1 mm
    base = trimesh.creation.box(extents=(20, 20, 4))
    base.apply_translation((0, 0, -5))
    out = _mesh(trimesh.util.concatenate([plate, base]), 1.0)
    assert out["skin_cover_max_error"] > 0.1
    assert out["status"] == 1                                # hors tolerance, dit dans le rapport


def test_cover_check_can_be_disabled():
    out = _mesh(trimesh.creation.box(extents=(10, 10, 10)), 1.0, {"check_cover": False})
    assert out["cover_samples"] == 0 and out["skin_cover_max_error"] == 0.0


def test_unknown_option_is_an_error_not_ignored():
    with pytest.raises(ValueError, match="option inconnue"):
        _mesh(trimesh.creation.box(extents=(10, 10, 10)), 1.0, {"refine_rounds": 3})
