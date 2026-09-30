"""
Banc des motifs de P4 (surface_paths.PATTERNS) sur les formes de test
(disque, demi-cylindre, demi-sphere) et, avec --real, sur trois zones reelles
de P2-bis (plus grande zone du tilted benchy = coque, plus grande du dino,
zone la plus plane du tilted benchy).

Mesures, independantes des champs qui ont servi a tracer :
  - espacement : distance 3D de chaque point d'une ligne de niveau k a la
    ligne k+1 (niveaux VOISINS ; les motifs a une seule courbe -- spirales --
    fournissent leurs tours comme niveaux) ; cible buse +-5 % (5e-95e) ;
  - couverture : part de la surface a moins de w/2 d'une trajectoire (>= 95 %) ;
  - hors zone : points dont la distance au bord est sous le premier niveau ;
  - duree du motif.
Sauve output_models/patterns_bench.pkl pour render_patterns.py.

Usage (venv d'Apex) :
  ../../.venv/Scripts/python.exe -u test_patterns.py [--real] [--only-real] [--patterns=a,b]
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


def corner_mask(pts, window=W, angle_deg=25.0):
    """Points a moins d'une buse (le long de la ligne) d'un virage de plus de
    `angle_deg` sur une fenetre d'une buse : dans un coin vif, deux lignes
    decalees sont a w/sin(theta/2) l'une de l'autre (0,57 mm pour un coin
    droit), ecart geometrique de tout decalage de contour, pas une erreur."""
    if len(pts) < 3:
        return np.zeros(len(pts), bool)
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    closed = np.linalg.norm(pts[0] - pts[-1]) < 1e-6
    ia = np.searchsorted(s, s - window / 2).clip(0, len(pts) - 1)
    ib = np.searchsorted(s, s + window / 2).clip(0, len(pts) - 1)
    a, b = pts - pts[ia], pts[ib] - pts
    na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
    ok = (na > 1e-9) & (nb > 1e-9)
    cos = np.ones(len(pts))
    cos[ok] = (a[ok] * b[ok]).sum(1) / (na[ok] * nb[ok])
    sharp = cos < np.cos(np.radians(angle_deg))
    if closed:
        sharp[0] = sharp[-1] = sharp[0] | sharp[-1]
    near = np.zeros(len(pts), bool)
    for q in np.flatnonzero(sharp):
        near |= np.abs(s - s[q]) <= window
    return near


def level_gaps(level_lines, zm, d_edge, skip_corners=False):
    """Ecart des points du niveau k au niveau k+1, mesure a l'INTERIEUR de la
    zone (a plus d'une buse du premier niveau) : au bord, une ligne coupee en
    biais a sa voisine plus courte, ecart geometrique et non d'espacement."""
    vtree = cKDTree(zm.V)
    inner = sp.first_level(W, OVERLAP) + W
    inside, full = {}, {}
    for k, pts in level_lines:
        full.setdefault(k, []).append(pts)
        keep = d_edge[vtree.query(pts)[1]] >= inner
        if skip_corners:
            keep &= ~corner_mask(pts)
        if keep.sum() >= 1:
            inside.setdefault(k, []).append(pts[keep])
    # points INTERIEURS du niveau k contre la ligne k+1 ENTIERE : comparer a
    # une ligne k+1 elle aussi rognee faisait manquer le voisin le plus proche
    # au bord du rognage (faux 0,43 mm sur la petite zone plane)
    gaps = [sp.gaps_between(inside[k], full[k + 1]) for k in sorted(inside) if k + 1 in full]
    g = np.concatenate(gaps) if gaps else np.empty(0)
    return g[g < 1.5 * W]              # au-dela : fin de motif (centre, axe median) sans voisin en face


def measure(name, zm, pattern, **kw):
    t0 = time.time()
    res = sp.PATTERNS[pattern](zm, W, OVERLAP, **kw)
    dt = time.time() - t0
    gaps_all = level_gaps(res.level_lines, zm, res.d_edge)
    gaps = level_gaps(res.level_lines, zm, res.d_edge, skip_corners=True)
    p5, p50, p95 = np.percentile(gaps, [5, 50, 95]) if len(gaps) else (np.nan,) * 3
    p95_all = np.percentile(gaps_all, 95) if len(gaps_all) else np.nan
    paths = [p for p in res.paths if len(p) >= 2]
    if not paths:
        print(f"{name:<24} {pattern:<20} AUCUNE trajectoire")
        return False, res
    tree = cKDTree(np.vstack([sp.densify(p) for p in paths]))
    cover = (tree.query(surface_samples(zm))[0] <= W / 2 + 1e-6).mean()
    vtree = cKDTree(zm.V)
    allp = np.vstack(paths)
    out = (res.d_edge[vtree.query(allp)[1]] < sp.first_level(W, OVERLAP) - 0.05).mean()
    length = sum(np.linalg.norm(np.diff(p, axis=0), axis=1).sum() for p in paths)
    travel = sum(np.linalg.norm(b[0] - a[-1]) for a, b in zip(paths[:-1], paths[1:]))
    ok = bool(len(gaps)) and (0.95 * W <= p5) and (p95 <= 1.05 * W) and cover >= 0.95
    print(f"{name:<24} {pattern:<20} {len(paths):>4} passes, {length:7.0f} mm, vide {travel:6.0f} mm | espacement {p50:.3f} "
          f"({p5:.3f}-{p95:.3f}, coins compris {p95_all:.3f}) | couverture {cover:.1%} | hors zone {out:.2%} | {dt:5.1f} s  {'OK' if ok else '<<'}")
    return ok, res


def real_zones():
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
    flag = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if "=" in a}
    patterns = flag["--patterns"].split(",") if "--patterns" in flag else [p for p in sp.PATTERNS if p != "hilbert"]   # Hilbert : test_hilbert.py
    shapes = [] if "--only-real" in sys.argv else [(n, zm) for n, zm, *_ in (disk(), half_cylinder(), hemisphere())]
    if "--real" in sys.argv or "--only-real" in sys.argv:
        shapes += real_zones()
    results, saved = [], {}
    for name, zm in shapes:
        d_edge = zm.distance_to_boundary()
        for pattern in patterns:
            ok, res = measure(name, zm, pattern, d_edge=d_edge)
            results.append(ok)
            saved[(name, pattern)] = (zm.V, zm.F, res.paths)
    path = os.path.join(HERE, "output_models", "patterns_bench.pkl")
    old = pickle.load(open(path, "rb")) if os.path.exists(path) else {}
    old.update(saved)
    pickle.dump(old, open(path, "wb"))
    print(f"\n{sum(results)}/{len(results)} conformes")
