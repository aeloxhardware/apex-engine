"""
Couverture de la SURFACE d'une piece par les cordons reellement deposes,
relus DEPUIS LE G-CODE (pivot compris) : un point de la surface du STL
d'origine est "decouvert" si aucun cordon extrude ne passe a moins de
`--tol` (defaut : demi-largeur 0,2 + 0,3 mm). Rend le taux de surface
decouverte, les amas de trous (composantes connexes de faces), et une image
de la piece avec les zones decouvertes en rouge.

Usage : .venv/Scripts/python.exe scripts/gcode_surface_coverage.py <gcode> <stl d'origine> <png> [--tol=0.5] [--pivot=46]
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import trimesh  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from file_io.gcode_export import firmware_to_machine, parse_gcode  # noqa: E402


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    gcode, stl, png = [a for a in sys.argv[1:] if not a.startswith("--")][:3]
    tol = float(flags.get("--tol", 0.5))
    pivot = float(flags.get("--pivot", 46.0))
    moves = parse_gcode(open(gcode).read().splitlines())
    tips = np.array([firmware_to_machine(x, y, z, b, pivot) for x, y, z, b, _ in moves])
    ext = np.array([m[4] for m in moves], dtype=bool)
    # cordons densifies a 0,2 mm (segments extrudes uniquement)
    seg = np.flatnonzero(ext[1:]) + 1
    a, b = tips[seg - 1], tips[seg]
    n = np.maximum(1, np.ceil(np.linalg.norm(b - a, axis=1) / 0.2).astype(int))
    rep = np.repeat(np.arange(len(seg)), n)
    off = np.arange(len(rep)) - np.repeat(np.cumsum(n) - n, n)
    t = ((off + 1) / np.repeat(n, n))[:, None]
    beads = a[rep] + t * (b[rep] - a[rep])
    tree = cKDTree(beads)
    part = trimesh.load(stl)
    lo, hi = part.bounds
    part.apply_translation([-(lo[0] + hi[0]) / 2, -(lo[1] + hi[1]) / 2, -lo[2]])
    # echantillonnage de surface ~4 points / mm2, avec la face de chaque point
    count = int(min(400_000, max(20_000, part.area * 4)))
    pts, face = trimesh.sample.sample_surface_even(part, count, seed=0)
    dist, _ = tree.query(pts, workers=-1)
    # la face du dessous (plateau) n'est pas une surface visible a couvrir
    bottom = pts[:, 2] < 0.3
    holes = (dist > tol) & ~bottom
    area_per_pt = part.area / len(pts)
    print(f"{os.path.basename(gcode)} : surface {part.area:.0f} mm2, decouverte (> {tol} mm d'un cordon) "
          f"{holes.mean():.2%} soit ~{holes.sum() * area_per_pt:.0f} mm2 ; ecart max {dist[~bottom].max():.2f} mm")
    # amas : faces dont un point est decouvert, regroupees par adjacence
    bad_faces = np.unique(face[holes])
    if len(bad_faces):
        adj = part.face_adjacency
        keep = np.isin(adj[:, 0], bad_faces) & np.isin(adj[:, 1], bad_faces)
        import scipy.sparse as sp
        index = {f: k for k, f in enumerate(bad_faces)}
        rows = [index[f] for f in adj[keep, 0]]
        cols = [index[f] for f in adj[keep, 1]]
        g = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(bad_faces),) * 2)
        ncomp, lab = connected_components(g, directed=False)
        sizes = np.bincount(lab, weights=part.area_faces[bad_faces])
        big = np.argsort(-sizes)[:8]
        cent = np.array([part.triangles_center[bad_faces[lab == c]].mean(axis=0) for c in big])
        print(f"  {ncomp} amas ; plus grands (aire des faces touchees, mm2 ; centre x y z) :")
        for c, ctr in zip(big, cent):
            print(f"    {sizes[c]:7.1f}  {np.round(ctr, 1)}")
    fig = plt.figure(figsize=(16, 8))
    for k, (el, az) in enumerate(((15, -60), (15, 120))):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        ok = ~holes & ~bottom
        s = np.flatnonzero(ok)[:: max(1, ok.sum() // 60000)]
        ax.scatter(*pts[s].T, s=0.2, c="seagreen", alpha=0.25)
        ax.scatter(*pts[holes].T, s=1.5, c="red")
        ax.view_init(el, az)
        ax.set_box_aspect(np.ptp(pts, axis=0))
        ax.set_title(f"{os.path.basename(gcode)} : rouge = surface sans cordon a moins de {tol} mm")
    plt.tight_layout()
    plt.savefig(png, dpi=100)


if __name__ == "__main__":
    main()
