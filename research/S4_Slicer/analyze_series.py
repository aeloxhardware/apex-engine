"""
Les deux criteres de la boucle automatique, calcules passe par passe sur une
serie produite par detect_folds.py :

  1. CONVERGENCE : deplacement median des noeuds depuis la passe precedente
     (l'utilisateur arretait le bareme quand une passe ne changeait plus rien).
  2. REPLI : auto-intersection de la SURFACE, jugee par l'enfoncement MEDIAN
     des croisements. Seuil propose : une hauteur de couche (0,2 mm) -- un
     croisement moins profond disparait a l'impression. Mesure : maillage
     juge sain par l'utilisateur 0,07 mm ; replis rejetes 0,52 a 1,76 mm.

A lancer avec le venv d'APEX (pas celui de la recherche) : le test
d'auto-intersection est celui de pipeline/step1_ingest.py.
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe analyze_series.py <prefixe>
  ex : analyze_series.py folds_real_overhang_a0.5_m1
"""
import glob
import os
import re
import sys

import numpy as np
import trimesh

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from pipeline.step1_ingest import _find_self_intersecting_faces, _triangles_cross  # noqa: E402

LAYER_HEIGHT = 0.2

prefix = os.path.join(os.path.dirname(__file__), "output_models", sys.argv[1])


def crossing_depth(a, b):
    """Enfoncement de a dans le plan de b (0 si a ne traverse pas ce plan)."""
    n = np.cross(b[1] - b[0], b[2] - b[0])
    n /= np.linalg.norm(n)
    d = (a - b[0]) @ n
    return min(-d.min(), d.max()) if d.min() < 0 < d.max() else 0.0


def surface_crossings(mesh):
    bad = np.array(sorted(_find_self_intersecting_faces(mesh)), dtype=np.int64)
    if len(bad) == 0:
        return 0, 0.0, 0.0
    T, F, C = mesh.triangles, mesh.faces, mesh.triangles_center
    reach = 2.0 * np.linalg.norm(T - C[:, None, :], axis=2).max()
    depths = []
    for k, i in enumerate(bad):
        for j in bad[k + 1:]:
            if set(F[i]) & set(F[j]) or np.linalg.norm(C[i] - C[j]) > reach:
                continue
            if _triangles_cross(T[i], T[j]):
                depths.append(min(crossing_depth(T[i], T[j]), crossing_depth(T[j], T[i])))
    depths = np.asarray(depths) if depths else np.zeros(1)
    return len(bad), float(np.median(depths)), float(depths.max())


def main():
    stls = {int(re.search(r"_pass(\d+)_ovh", f).group(1)): f for f in glob.glob(f"{prefix}_pass*_ovh*.stl")}
    points = {k: np.load(f"{prefix}_pass{k}.npy") for k in [0] + sorted(stls)}
    diag = np.linalg.norm(np.ptp(points[0], axis=0))
    print(f"{os.path.basename(prefix)} : {len(points[0])} noeuds, diagonale {diag:.1f} mm")
    print(f"{'passe':>5} {'MAX_OVH':>7} {'depl. median':>13} {'faces croisees':>15} "
          f"{'enfonc. median':>15} {'max':>7}  verdict")
    for k in sorted(stls):
        ovh = re.search(r"_ovh([\d.]+)\.stl$", stls[k]).group(1)
        move = float(np.median(np.linalg.norm(points[k] - points[k - 1], axis=1)))
        n_bad, med, mx = surface_crossings(trimesh.load(stls[k]))
        verdict = "REPLI" if med > LAYER_HEIGHT else "sain"
        print(f"{k:>5} {ovh:>7} {move:>11.4f} mm {n_bad:>15} {med:>12.3f} mm {mx:>6.3f}  {verdict}")


if __name__ == "__main__":
    main()
