"""
Rendu de controle d'un G-code produit par file_io/gcode_export.py, relu
DEPUIS LE FICHIER (pas depuis la trajectoire interne), sous DEUX lectures du
modele de tete :
  - "apex"  : le firmware commande le PIVOT de B, la pointe est a 45 mm
              (modele de l'export, limite #61) ;
  - "brute" : X et Z lus comme la position de la POINTE (aucune compensation),
              lecture d'un previsualiseur qui ne connait pas le pivot.

Mesures, par lecture :
  - points de pointe sous le plateau (z < 0) ;
  - points EXTRUDES hors de la piece (distance a la surface du STL d'origine,
    ramene dans le repere machine, et hors du volume) ;
  - segments extrudes longs (> 2 mm : un cordon qui traverse du vide) ;
  - trajets a vide dont le milieu est dans la piece.
Image : vues de face et de dessus, extrusion coloree par ordre d'impression,
trajets a vide en rouge, points sous le plateau en noir.

Usage : .venv/Scripts/python.exe scripts/render_gcode_check.py <gcode> <stl d'origine> <png> [--pivot=45]
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

from file_io.gcode_export import firmware_to_machine, parse_gcode  # noqa: E402


def tips_of(moves, pivot):
    fw = np.array([m[:4] for m in moves], dtype=float)
    return np.array([firmware_to_machine(x, y, z, b, pivot) for x, y, z, b in fw])


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    gcode, stl, png = [a for a in sys.argv[1:] if not a.startswith("--")][:3]
    pivot = float(flags.get("--pivot", 46.0))
    moves = parse_gcode(open(gcode).read().splitlines())
    ext = np.array([m[4] for m in moves], dtype=bool)
    part = trimesh.load(stl)
    lo, hi = part.bounds
    part.apply_translation([-(lo[0] + hi[0]) / 2, -(lo[1] + hi[1]) / 2, -lo[2]])   # repere machine
    rng = np.random.default_rng(0)
    readings = {"apex (pivot 45 mm)": tips_of(moves, pivot), "brute (X,Z = pointe)": tips_of(moves, 0.0)}
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    report = []
    for row, (name, tips) in enumerate(readings.items()):
        seg_ext = ext[1:]
        seg_len = np.linalg.norm(np.diff(tips, axis=0), axis=1)
        long_ext = np.flatnonzero(seg_ext & (seg_len > 2.0))
        below = tips[:, 2] < 0.0
        # hors piece : echantillon de points extrudes (distance signee trimesh : >0 dedans)
        idx = np.flatnonzero(ext)
        sample = rng.choice(idx, min(40000, len(idx)), replace=False)
        signed = trimesh.proximity.signed_distance(part, tips[sample])
        outside = -signed                                  # >0 : hors de la piece
        travel_mid = 0.5 * (tips[:-1] + tips[1:])[~seg_ext]
        tsample = travel_mid[rng.choice(len(travel_mid), min(20000, len(travel_mid)), replace=False)] \
            if len(travel_mid) else np.zeros((0, 3))
        travel_inside = (trimesh.proximity.signed_distance(part, tsample) > 0.5).mean() if len(tsample) else 0.0
        report.append(
            f"[{name}] pointe z min {tips[:, 2].min():.2f} mm ; sous le plateau {below.sum()} pts "
            f"({below.mean():.1%}) ; extrudes hors piece > 1 mm : {(outside > 1.0).mean():.1%} "
            f"(> 3 mm : {(outside > 3.0).mean():.1%}, max {outside.max():.1f} mm) ; "
            f"segments extrudes > 2 mm : {len(long_ext)} (max {seg_len[seg_ext].max():.1f} mm) ; "
            f"trajets a vide traversant la piece : {travel_inside:.1%}")
        step = max(1, len(tips) // 120000)
        sel = np.arange(0, len(tips), step)
        color = np.linspace(0, 1, len(tips))
        for col, (i, j, title) in enumerate(((0, 2, "face (x, z)"), (0, 1, "dessus (x, y)"))):
            ax = axes[row, col]
            m = sel[ext[sel]]
            ax.scatter(tips[m, i], tips[m, j], s=0.2, c=color[m], cmap="viridis")
            t = sel[~ext[sel]]
            ax.scatter(tips[t, i], tips[t, j], s=0.4, c="red", alpha=0.4)
            b = np.flatnonzero(below)
            if len(b):
                ax.scatter(tips[b, i], tips[b, j], s=1.5, c="black")
            for k in long_ext[:400]:
                ax.plot(tips[k:k + 2, i], tips[k:k + 2, j], c="magenta", lw=0.6)
            if j == 2:
                ax.axhline(0.0, c="k", lw=1)
            ax.set_aspect("equal")
            ax.set_title(f"{os.path.basename(gcode)} - lecture {name} - {title}", fontsize=9)
    fig.text(0.01, 0.005, "\n".join(report) + "\nviolet->jaune : ordre d'impression ; rouge : trajets a vide ; "
             "noir : pointe sous le plateau ; magenta : segment extrude > 2 mm", fontsize=8, va="bottom")
    plt.tight_layout(rect=(0, 0.06, 1, 1))
    plt.savefig(png, dpi=110)
    print("\n".join(report))


if __name__ == "__main__":
    main()
