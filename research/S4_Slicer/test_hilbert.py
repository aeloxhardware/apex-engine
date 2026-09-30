"""
Etape 4d : Hilbert sur depliage ARAP -- MESURE de l'ecart (aucun remede :
decision de l'utilisateur). Pour chaque forme du banc :
  - etirement du depliage (valeurs singulieres s1 >= s2 par triangle, 1 =
    isometrie) : centiles ponderes par l'aire, part de surface a plus de 5 % ;
  - espacement reel sur la surface entre passes voisines (point -> portion de
    courbe la plus proche a plus de 1,5 w le long de la courbe), a l'interieur
    de la zone ; cible w +- 5 % ;
  - couverture (part de surface a moins de w/2), durees.
Ajoute ("forme", "hilbert") a output_models/patterns_bench.pkl.

Usage (venv d'Apex) : ../../.venv/Scripts/python.exe -u test_hilbert.py [--real]
"""
import os
import pickle
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

import surface_paths as sp
from test_patterns import HERE, OVERLAP, W, real_zones, surface_samples
from test_surface_paths import disk, half_cylinder, hemisphere


def weighted_pct(values, weights, pct):
    order = np.argsort(values)
    cum = np.cumsum(weights[order]) / weights.sum()
    return np.interp(np.asarray(pct) / 100.0, cum, values[order])


def run(name, zm):
    d_edge = zm.distance_to_boundary()
    t0 = time.time()
    uv = sp.arap_unfold(zm)
    t_unfold = time.time() - t0
    sv, area = sp.arap_stretch(zm, uv)
    dev = np.abs(sv - 1).max(axis=1)
    s1 = weighted_pct(sv[:, 0], area, [50, 95])
    s2 = weighted_pct(sv[:, 1], area, [5, 50])
    over5 = area[dev > 0.05].sum() / area.sum()
    t0 = time.time()
    res = sp.hilbert(zm, W, OVERLAP, d_edge=d_edge, uv=uv)
    t_path = time.time() - t0
    gaps, pts = sp.path_neighbor_gaps(res.paths, W)
    vtree = cKDTree(zm.V)
    inner = d_edge[vtree.query(pts)[1]] >= sp.first_level(W, OVERLAP) + W
    g = gaps[inner & np.isfinite(gaps)]
    g = g[g < 1.5 * W]
    p5, p50, p95 = np.percentile(g, [5, 50, 95])
    within = ((g >= 0.95 * W) & (g <= 1.05 * W)).mean()
    tree = cKDTree(np.vstack([sp.densify(p) for p in res.paths]))
    cover = (tree.query(surface_samples(zm))[0] <= W / 2 + 1e-6).mean()
    print(f"{name:<26} etirement s1 {s1[0]:.3f}/{s1[1]:.3f} (med/95e) s2 {s2[1]:.3f}/{s2[0]:.3f} (med/5e), "
          f"surface > 5 % : {over5:6.1%} | espacement {p50:.3f} ({p5:.3f}-{p95:.3f}), dans +-5 % : {within:6.1%} "
          f"| couverture {cover:.1%} | {len(res.paths)} passes | ARAP {t_unfold:.1f} s, trace {t_path:.1f} s",
          flush=True)
    return (zm.V, zm.F, res.paths)


if __name__ == "__main__":
    shapes = [(n, zm) for n, zm, *_ in (disk(), half_cylinder(), hemisphere())]
    if "--real" in sys.argv:
        shapes += real_zones()
    saved = {(name, "hilbert"): run(name, zm) for name, zm in shapes}
    path = os.path.join(HERE, "output_models", "patterns_bench.pkl")
    old = pickle.load(open(path, "rb")) if os.path.exists(path) else {}
    old.update(saved)
    pickle.dump(old, open(path, "wb"))
