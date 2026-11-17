"""TetraNet : le soudage des sommets reste lineaire sur des coordonnees issues d'un STL (float32).

Defaut mesure le 2026-11-17 (Benchy officiel, 112 569 sommets) : la fonction de hachage du soudage ne melangeait pas les
bits ; or un double issu d'un float32 a ses 29 bits de poids faible nuls, donc presque toutes les cles tombaient dans les
memes buckets de `unordered_map` (puissance de 2 sous MSVC) : soudage 45,9 s sur 48,6 s de maillage, 2,4 s apres
correction (splitmix64). Meme maillage en sortie (204 438 tetraedres). Un STL de grande taille faisait 149 s (678 000 sommets).

Le test ne verifie pas une valeur de temps fine : il borne le temps de tout l'appel avec une marge d'un ordre de grandeur,
sur un maillage dont le soudage seul depassait 45 s avec l'ancien hachage."""
import time

import numpy as np
import pytest
import trimesh

from core import _cpp

pytestmark = pytest.mark.skipif(
    not _cpp.available() or not hasattr(_cpp.module(), "tetranet_mesh"),
    reason="noyau natif avec TetraNet absent",
)

#: Marge large : l'appel corrige dure ~1-2 s ici ; l'ancien hachage depassait 45 s a ce nombre de sommets.
MAX_SECONDS = 15.0


def _float32_sphere():
    sphere = trimesh.creation.icosphere(subdivisions=7, radius=10.0)      # 163 842 sommets
    sphere.vertices = np.asarray(sphere.vertices, dtype=np.float32).astype(np.float64)   # comme un STL
    return sphere


def test_weld_of_float32_coordinates_is_fast():
    sphere = _float32_sphere()
    vertices = np.ascontiguousarray(sphere.vertices, dtype=np.float64)
    assert not np.any(vertices.view(np.uint64) & np.uint64((1 << 29) - 1)), "les coordonnees doivent venir de float32"
    start = time.perf_counter()
    out = _cpp.module().tetranet_mesh(vertices, np.ascontiguousarray(sphere.faces, dtype=np.int64),
                                      6.0, 0.5, 0.05, 0.05, 0, 0.0)
    elapsed = time.perf_counter() - start
    assert out["tets"] > 0
    assert elapsed < MAX_SECONDS, f"maillage de {len(vertices)} sommets float32 : {elapsed:.1f} s (soudage non lineaire ?)"
