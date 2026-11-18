"""Progression REELLE de TetraNet : trois phases (preparation, tetraedralisation, verification), chacune de 0 a 1,
croissante, issue de compteurs dans les boucles du noyau (jamais d'une duree estimee)."""
import numpy as np
import pytest
import trimesh

from core import _cpp, progress
from element_ardor import step1_build


def _native():
    if not _cpp.available() or not hasattr(_cpp.module(), "tetranet_mesh"):
        pytest.skip("noyau natif TetraNet absent")
    module = _cpp.module()
    if "progress" not in (module.tetranet_mesh.__doc__ or ""):
        pytest.skip("noyau natif plus ancien que la progression : scripts/build_cpp.py --dev")
    return module


def _mesh(progress_cb):
    sphere = trimesh.creation.icosphere(subdivisions=4, radius=20)
    return _native().tetranet_mesh(np.ascontiguousarray(sphere.vertices, dtype=np.float64),
                                   np.ascontiguousarray(sphere.faces, dtype=np.int64), 1.0, progress=progress_cb)


def test_trois_phases_croissantes_de_0_a_1():
    events = []
    out = _mesh(lambda phase, fraction: events.append((phase, fraction)))
    assert out["tets"] > 0
    assert [p for p, _ in events][0] == 0 and [p for p, _ in events][-1] == 2
    for phase in (0, 1, 2):
        fractions = [f for p, f in events if p == phase]
        assert fractions[0] == 0.0 and fractions[-1] == 1.0
        assert fractions == sorted(fractions)
        assert len(set(fractions)) > 3          # de vraies etapes intermediaires, pas 0 puis 1


def test_une_erreur_du_rappel_ne_fait_pas_echouer_le_maillage():
    def boom(phase, fraction):
        raise RuntimeError("barre cassee")
    assert _mesh(boom)["tets"] > 0


def test_sans_abonne_aucun_rappel_est_passe_au_noyau():
    progress.set_substep_sink(None)
    assert step1_build._tetranet_progress_kwargs() == {}


def test_les_phases_du_noyau_deviennent_les_sous_etapes():
    seen = []
    progress.set_substep_sink(lambda name, f: seen.append((name, f)))
    try:
        step1_build._tetranet_progress_kwargs()["progress"](1, 0.5)
    finally:
        progress.set_substep_sink(None)
    assert seen == [("tetrahedralize", 0.5)]


def test_soudure_des_sommets_rapide_sur_coordonnees_float32():
    """Regression : un STL est en float32, donc 29 bits de poids faible nuls une fois en double. L'ancien hachage de la
    soudure les propageait et tombait en collisions (678 000 sommets : 149 s). Maille grossiere : le temps mesure est
    celui de la soudure, pas du maillage."""
    import time
    sphere = trimesh.creation.icosphere(subdivisions=6, radius=20)
    soup = sphere.vertices.astype(np.float32).astype(np.float64)[sphere.faces].reshape(-1, 3)   # sommets non fusionnes
    faces = np.arange(len(soup), dtype=np.int64).reshape(-1, 3)
    start = time.perf_counter()
    out = _native().tetranet_mesh(np.ascontiguousarray(soup), faces, 10.0)
    elapsed = time.perf_counter() - start
    assert out["surface_faces_count"] == len(faces)
    assert elapsed < 6.0, f"soudure de {len(soup)} sommets : {elapsed:.1f} s (collisions de hachage ?)"
