"""
Banc de P4 (etape 4b) : concentrique conforme, cercles, rectiligne sur les
formes de test (disque, demi-cylindre, demi-sphere) et, avec --real, sur des
zones reelles de P2-bis. Mesures : espacement entre lignes voisines (3D),
couverture (part de la surface a moins de w/2 d'une trajectoire), points hors
zone (distance au bord < premier niveau), duree. Sauve les trajectoires pour
render_patterns.py.

Usage (venv d'Apex) : ../../.venv/Scripts/python.exe -u test_patterns_4b.py [--real]
"""
import os
import pickle
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

import surface_paths as sp
from test_surface_paths import disk, half_cylinder, hemisphere

W = sp.NOZZLE_WIDTH
OVERLAP = 0.2
HERE = os.path.dirname(os.path.abspath(__file__))


def surface_samples(zm, n=20000, seed=0):
    tri = zm.V[zm.F]
    area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    rng = np.random.default_rng(seed)
    f = rng.choice(len(tri), n, p=area / area.sum())
    a, b = rng.random(n), rng.random(n)
    swap = a + b > 1
    a[swap], b[swap] = 1 - a[swap], 1 - b[swap]
    return tri[f, 0] + a[:, None] * (tri[f, 1] - tri[f, 0]) + b[:, None] * (tri[f, 2] - tri[f, 0])


def neighbor_gaps(lines):
    """Pour chaque point, distance a la ligne VOISINE la plus proche (autre ligne)."""
    dense = [sp.densify(pl) for pl in lines if len(pl) >= 2]
    if len(dense) < 2:
        return np.empty(0)
    owner = np.concatenate([np.full(len(d), k) for k, d in enumerate(dense)])
    pts = np.vstack(dense)
    tree = cKDTree(pts)
    gaps = []
    for k, pl in enumerate(lines):
        if len(pl) < 2:
            continue
        dist, idx = tree.query(pl, k=64, distance_upper_bound=3 * W)
        other = (owner[np.minimum(idx, len(owner) - 1)] != k) & np.isfinite(dist)
        best = np.where(other, dist, np.inf).min(axis=1)
        gaps.append(best[np.isfinite(best)])
    return np.concatenate(gaps) if gaps else np.empty(0)


def measure(name, zm, pattern):
    t0 = time.time()
    if pattern == "concentric_conform":
        lines, d_edge = sp.concentric_conform(zm, W, OVERLAP)
        lines = [pl.points for pl in lines]
        spacing_lines = lines
    elif pattern == "concentric":
        lines, d_edge = sp.concentric_circles(zm, W, OVERLAP)
        lines = [pl.points for pl in lines]
        spacing_lines = lines
    else:
        spacing_lines, d_edge = sp.rectilinear(zm, W, OVERLAP, zigzag=False)
        lines, _ = sp.rectilinear(zm, W, OVERLAP, d_edge=d_edge)
    dt = time.time() - t0
    gaps = neighbor_gaps(spacing_lines)
    gaps = gaps[gaps < 1.5 * W]                       # voisin reel (hors bouts de lignes)
    p5, p50, p95 = np.percentile(gaps, [5, 50, 95]) if len(gaps) else (np.nan,) * 3
    samples = surface_samples(zm)
    tree = cKDTree(np.vstack([sp.densify(pl) for pl in lines if len(pl) >= 2]))
    cover = (tree.query(samples)[0] <= W / 2 + 1e-6).mean()
    # hors zone : distance au bord (champ interpole au sommet le plus proche) sous le premier niveau
    vtree = cKDTree(zm.V)
    allp = np.vstack(lines)
    out = (d_edge[vtree.query(allp)[1]] < sp.first_level(W, OVERLAP) - 0.05).mean()
    length = sum(np.linalg.norm(np.diff(pl, axis=0), axis=1).sum() for pl in lines)
    ok = (0.95 * W <= p5) and (p95 <= 1.05 * W) and cover >= 0.95
    print(f"{name:<22} {pattern:<19} {len(lines):>4} passes, {length:7.0f} mm | espacement {p50:.3f} "
          f"({p5:.3f}-{p95:.3f}) | couverture {cover:.1%} | hors zone {out:.2%} | {dt:.1f} s  {'OK' if ok else '<<'}")
    return ok, lines


def real_zones():
    """Trois zones reelles : la plus grande du tilted benchy (coque), la plus
    grande du dino, la plus plane du tilted benchy (normale la plus constante)."""
    out = []
    for model, k in (("tilted_benchy_undertriangulated_repaired_pin", 4), ("dino_pin", 5)):
        skin = np.load(os.path.join(HERE, "output_models", f"loop_{model}_skin_pass{k}.npz"))
        z, a, n = skin["zone"], skin["areas"], skin["normals"]
        areas = np.bincount(z[z >= 0], weights=a[z >= 0])
        big = int(np.argmax(areas))
        out.append((f"{model.split('_')[0]} z{big} ({areas[big]:.0f} mm2)", skin, big))
        if model.startswith("tilted"):
            spread = [np.linalg.norm((n[z == q] * a[z == q, None]).sum(0)) / a[z == q].sum() if areas[q] > 20 else 0
                      for q in range(len(areas))]
            flat = int(np.argmax(spread))
            out.append((f"benchy plane z{flat} ({areas[flat]:.0f} mm2)", skin, flat))
    return [(name, sp.ZoneMesh(skin["vertices"], skin["faces"][skin["zone"] == q])) for name, skin, q in out]


if __name__ == "__main__":
    shapes = [(n, zm) for n, zm, *_ in (disk(), half_cylinder(), hemisphere())]
    if "--real" in sys.argv:
        shapes += real_zones()
    results, saved = [], {}
    for name, zm in shapes:
        for pattern in ("concentric_conform", "concentric", "rectilinear"):
            ok, lines = measure(name, zm, pattern)
            results.append(ok)
            saved[(name, pattern)] = (zm.V, zm.F, lines)
    pickle.dump(saved, open(os.path.join(HERE, "output_models", "patterns_4b.pkl"), "wb"))
    print(f"\n{sum(results)}/{len(results)} conformes")
