"""
Aperçu d'un G-code produit par file_io/gcode_export.py : relit le fichier,
rejoue la cinematique inverse (pivot -> pointe, puis R C Z -> x y z) et rend
la trajectoire de la POINTE dans le repere de la piece (PNG, matplotlib).

Usage :
  .venv/Scripts/python.exe scripts/preview_gcode_file.py <fichier.gcode> [--pivot=45]
      [--toolpath=<_machine_toolpath_pass*.pkl>] [--out=<png>]

--toolpath : controle ALLER-RETOUR. Chaque point de la trajectoire d'origine
doit etre retrouve dans le fichier relu a 1e-3 mm pres (valide le pivot, les
signes et le deroulement de Y) ; l'ecart max est affiche.
"""
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from file_io.gcode_export import firmware_to_machine, parse_gcode  # noqa: E402


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    path = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    pivot = float(flags.get("--pivot", 45.0))
    moves = parse_gcode(open(path).read().splitlines())
    tips = np.array([firmware_to_machine(x, y, z, b, pivot) for x, y, z, b, _ in moves])
    extrude = np.array([m[4] for m in moves])
    print(f"{len(moves)} mouvements, dont {extrude.sum()} en extrusion ; boite de la pointe "
          f"{tips.min(axis=0).round(2)} -> {tips.max(axis=0).round(2)}")
    b = np.array([m[3] for m in moves])
    print(f"B : {b.min():.1f} -> {b.max():.1f} deg ; Z pointe min {tips[:, 2].min():.3f} mm")
    if "--toolpath" in flags:
        from scipy.spatial import cKDTree
        ref = np.asarray(pickle.load(open(flags["--toolpath"], "rb"))["position"], dtype=float)
        d, _ = cKDTree(tips).query(ref)
        print(f"aller-retour : ecart max {d.max():.2e} mm, points a plus de 1e-3 mm : {(d > 1e-3).sum()}/{len(ref)}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    segs = np.stack([tips[:-1], tips[1:]], axis=1)
    ext = extrude[1:]
    fig = plt.figure(figsize=(12, 6))
    for n, (elev, azim) in enumerate(((25, -60), (0, -90))):
        ax = fig.add_subplot(1, 2, n + 1, projection="3d")
        order = np.linspace(0, 1, len(segs))
        lc = Line3DCollection(segs[ext], colors=plt.cm.viridis(order[ext]), linewidths=0.3)
        ax.add_collection3d(lc)
        if (~ext).any():
            ax.add_collection3d(Line3DCollection(segs[~ext], colors=(0.6, 0.6, 0.6, 0.25), linewidths=0.2))
        lo, hi = tips.min(axis=0), tips.max(axis=0)
        c, r = (lo + hi) / 2, (hi - lo).max() / 2
        ax.set_xlim(c[0] - r, c[0] + r)
        ax.set_ylim(c[1] - r, c[1] + r)
        ax.set_zlim(c[2] - r, c[2] + r)
        ax.view_init(elev, azim)
        ax.set_title("extrusion (ordre : violet -> jaune), trajets a vide en gris", fontsize=8)
    out = flags.get("--out", os.path.splitext(path)[0] + "_preview.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(out)


if __name__ == "__main__":
    main()
