"""
core/spatial.py — structures d'acceleration spatiale, chemin Python OU natif.

Deux enveloppes minces, meme discipline : le choix scipy/trimesh vs
`apex_core` est fige A LA CONSTRUCTION selon `core._cpp.enabled()` (comme
`_PrintedGeometryIndex._native`), et un objet cree dans une position de
l'interrupteur ne bascule jamais tout seul.

- `KDTree` — enveloppe de `scipy.spatial.cKDTree` (pas 5). Portage
  ARCHITECTURAL : `cKDTree` est deja natif (C), le gain de vitesse attendu
  est nul ; le but est de retirer la derniere primitive spatiale de scipy du
  coeur.
- `RayIntersector` — enveloppe de `trimesh.ray.ray_triangle.RayMeshIntersector`
  (pas 6). Ici la reference Python n'est PAS native (`pyembree` absent de ce
  venv, verifie) : le gain est reel, mesure ×32 sur death.stl (2418 faces) et
  ×139 sur death_hi-res.stl (226 156 faces), a impacts identiques au BIT.

--- KDTree ---------------------------------------------------------------

Restreinte au sous-ensemble REELLEMENT utilise par le pipeline (verifie site
par site) :

- construction depuis un tableau `(N, 3)` ;
- `.query(points)` — plus proche voisin unique (`k=1`), la seule variante
  appelee dans tout `pipeline/` ;
- `.query_ball_point(point, radius)` — indices dans un rayon autour d'UN
  point (`core/constraints.py`, seul appelant) ;
- l'attribut `.data` — les points, relus par `step3_segmentation`.

--- RayIntersector -------------------------------------------------------

UN SEUL point d'entree est utilise par le pipeline : `intersects_location`,
appele par `step2_overhang` (detection de ponts) et `step7_scheduling`
(deadline d'accessibilite). `intersects_first`, `intersects_any` et
`contains_points` ne sont appeles nulle part — non portes.

Les deux appelants IGNORENT `index_tri` et sont insensibles a l'ordre (ils
reduisent par `.min()` sur les distances). Voir `cpp/include/apex/bvh.hpp`
pour le detail de l'equivalence et la verification que le broad phase de
trimesh ne rate aucun impact.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import scipy.spatial
import trimesh.ray.ray_triangle

from core import _cpp


class KDTree:
    """Voir le docstring du module. Compatible canard avec le sous-ensemble
    de `scipy.spatial.cKDTree` liste ci-dessus, rien de plus."""

    def __init__(self, points: Any):
        self._data = np.ascontiguousarray(points, dtype=np.float64)
        if self._data.ndim != 2 or self._data.shape[1] != 3:
            raise ValueError(f"points doit etre de forme (N, 3), recu {self._data.shape}")

        # `None` => chemin scipy historique, strictement inchange. Meme
        # discipline que `_PrintedGeometryIndex._native`.
        self._native = None
        self._scipy: scipy.spatial.cKDTree | None = None
        if _cpp.enabled():
            self._native = _cpp.module().KDTree(self._data)
        else:
            self._scipy = scipy.spatial.cKDTree(self._data)

    @property
    def data(self) -> np.ndarray:
        """Les points de l'arbre, comme `cKDTree.data`."""
        return self._data

    def query(self, x: Any):
        """Plus proche voisin (`k=1`). Rend `(distances, indices)` pour une
        entree `(M, 3)`, ou `(distance, index)` scalaires pour une entree
        `(3,)` — meme contrat que `cKDTree.query` a `k=1`."""
        pts = np.ascontiguousarray(x, dtype=np.float64)
        if self._native is None:
            return self._scipy.query(pts)

        one_d = pts.ndim == 1
        queries = pts.reshape(1, 3) if one_d else pts
        distances, indices = self._native.query_nn(queries)
        if one_d:
            return float(distances[0]), int(indices[0])
        return distances, indices

    def query_ball_point(self, point: Any, radius: float):
        """Indices des points a distance <= `radius` de `point`. Ordre non
        specifie, comme `cKDTree.query_ball_point` (le seul appelant
        re-filtre derriere)."""
        p = np.ascontiguousarray(point, dtype=np.float64).ravel()
        if self._native is None:
            return self._scipy.query_ball_point(p, radius)
        return self._native.query_radius(p, float(radius))

    def __len__(self) -> int:
        return int(self._data.shape[0])


class RayIntersector:
    """Voir le docstring du module. Compatible canard avec
    `RayMeshIntersector.intersects_location`, rien de plus."""

    def __init__(self, mesh: Any):
        self._native = None
        self._trimesh: Any = None
        if _cpp.enabled():
            # `face_normals` est PASSEE et non recalculee cote C++ : trimesh
            # utilise `mesh.face_normals` quand elles existent, et une normale
            # recalculee donnerait d'autres flottants donc d'autres resultats
            # aux bords de tolerance.
            self._native = _cpp.module().TriangleBVH(
                np.ascontiguousarray(mesh.triangles, dtype=np.float64),
                np.ascontiguousarray(mesh.face_normals, dtype=np.float64),
            )
        else:
            self._trimesh = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)

    def intersects_location(self, ray_origins: Any, ray_directions: Any, **kwargs):
        """Impacts dedoublonnes sur (location, index_ray). Rend le triplet de
        trimesh : `(locations, index_ray, index_tri)`."""
        if self._native is None:
            return self._trimesh.intersects_location(ray_origins, ray_directions, **kwargs)
        if kwargs:
            # Aucun appelant du pipeline ne passe d'option ; une option
            # silencieusement ignoree serait un ecart de comportement invisible.
            raise TypeError(f"options non supportees par le chemin natif : {sorted(kwargs)}")
        return self._native.intersects_location(
            np.ascontiguousarray(ray_origins, dtype=np.float64).reshape(-1, 3),
            np.ascontiguousarray(ray_directions, dtype=np.float64).reshape(-1, 3),
        )
