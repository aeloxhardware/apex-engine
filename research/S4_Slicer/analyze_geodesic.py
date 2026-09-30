"""
Repli = deux morceaux de surface qui se croisent alors qu'ils etaient LOIN
l'un de l'autre EN SUIVANT LA SURFACE de la piece d'origine.

Pourquoi (mesure le 2026-09-26) : l'enfoncement des croisements ne separe
pas un repli d'un plissement. Le repli fabrique par l'utilisateur sur le
Benchy (le toit qui fusionne avec la cabine) s'enfonce en mediane de
0,089 mm, a peine plus que des plissements sains (0,07 mm) : les deux
parties glissent l'une dans l'autre. Ce qui les distingue, c'est l'origine
des deux triangles qui se croisent :
  - plissement : la meme paroi qui se froisse -> triangles voisins sur la
    surface d'origine (distance geodesique de quelques mm) ;
  - repli : deux parois distinctes qui se rejoignent -> triangles eloignes
    en suivant la surface, meme s'ils etaient proches a vol d'oiseau.

Pour chaque paire de triangles qui se croisent, on retrouve leurs sommets
dans la piece d'origine (positions pass0) et on calcule la distance
geodesique le long de la surface d'origine (Dijkstra sur le graphe des
aretes de surface).

A lancer avec le venv d'APEX :
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe analyze_geodesic.py <prefixe>
"""
import glob
import os
import re
import sys

import numpy as np
import scipy.sparse as sp
import trimesh
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from pipeline.step1_ingest import _find_self_intersecting_faces, _triangles_cross  # noqa: E402

prefix = os.path.join(os.path.dirname(__file__), "output_models", sys.argv[1])
# Au-dela de cette distance geodesique, l'enfoncement cesse d'etre un
# plissement local. Plafond de la recherche Dijkstra : au-dela, "loin".
FAR_MM = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
SEARCH_LIMIT_MM = 60.0


def crossing_pairs(mesh):
    bad = np.array(sorted(_find_self_intersecting_faces(mesh)), dtype=np.int64)
    T, F, C = mesh.triangles, mesh.faces, mesh.triangles_center
    reach = 2.0 * np.linalg.norm(T - C[:, None, :], axis=2).max()
    pairs = []
    for k, i in enumerate(bad):
        for j in bad[k + 1:]:
            if set(F[i]) & set(F[j]) or np.linalg.norm(C[i] - C[j]) > reach:
                continue
            if _triangles_cross(T[i], T[j]):
                pairs.append((i, j))
    return np.array(pairs, dtype=np.int64).reshape(-1, 2)


def main():
    stls = {int(re.search(r"_pass(\d+)_ovh", f).group(1)): f for f in glob.glob(f"{prefix}_pass*_ovh*.stl")}
    p0 = np.load(f"{prefix}_pass0.npy")
    print(f"{os.path.basename(prefix)} -- seuil 'loin sur la surface' : {FAR_MM} mm")
    print(f"{'passe':>5} {'paires':>7} {'dist. geo. mediane':>19} {'paires LOIN':>12} {'faces LOIN':>11} {'plus longue':>12}")
    for k in sorted(stls):
        mesh = trimesh.load(stls[k])
        pk = np.load(f"{prefix}_pass{k}.npy")
        # bird_core.export_surface_stl recentre avant d'ecrire (milieu en XY,
        # pose a z=0). detect_folds sauve les positions APRES ce recentrage,
        # server.py AVANT : on l'applique ici (sans effet si deja fait).
        lo, hi = pk.min(axis=0), pk.max(axis=0)
        pk = pk - np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, lo[2]])
        # sommet STL (float32) -> noeud du tetraedre -> position d'origine
        dist, node = cKDTree(pk).query(mesh.vertices)
        assert dist.max() < 1e-3, f"correspondance STL/noeuds rompue ({dist.max():.2e} mm)"
        orig = p0[node]

        e = mesh.edges_unique
        w = np.linalg.norm(orig[e[:, 0]] - orig[e[:, 1]], axis=1)
        n = len(orig)
        graph = sp.coo_matrix((np.r_[w, w], (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(n, n)).tocsr()

        pairs = crossing_pairs(mesh)
        if len(pairs) == 0:
            print(f"{k:>5} {0:>7}")
            continue
        va = mesh.faces[pairs[:, 0], 0]
        vb = mesh.faces[pairs[:, 1], 0]
        sources, inv = np.unique(va, return_inverse=True)
        d = dijkstra(graph, directed=False, indices=sources, limit=SEARCH_LIMIT_MM)
        geo = d[inv, vb]
        geo = np.where(np.isfinite(geo), geo, SEARCH_LIMIT_MM)
        far = geo > FAR_MM
        far_faces = len(np.unique(pairs[far].ravel()))
        print(f"{k:>5} {len(pairs):>7} {np.median(geo):>16.1f} mm {far.sum():>12} {far_faces:>11} {geo.max():>9.1f} mm")


if __name__ == "__main__":
    main()
