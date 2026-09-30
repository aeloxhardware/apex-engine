"""
Image de controle de P3 : contours de reference (gris) et creuses (rouge) de
quelques couches, avec un zoom la ou l'ecart (epaisseur de lissage) est
visible, et la premiere couche (motif choisi).

Usage (venv de la recherche) :
  .venv/Scripts/python.exe render_carve_compare.py <modele> <passe> <variante>
Sortie : output_models/loop_<...>_carve_compare_pass<k>.png
"""
import os
import pickle
import sys

import matplotlib
import numpy as np
from scipy.spatial import cKDTree

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
args = [a for a in sys.argv[1:] if not a.startswith("--")]
model, k, variant = args[0], int(args[1]), args[2]
prefix = os.path.join(HERE, "output_models", "loop_" + model.replace(" ", "_") + variant)
ref = {L["isovalue"]: L for L in pickle.load(open(f"{prefix}_deformed_slice_pass{k}.pkl", "rb"))}
carved = sorted(pickle.load(open(f"{prefix}_deformed_slice_skin_pass{k}.pkl", "rb")), key=lambda L: L["isovalue"])


def dense(loops, step=0.05):
    out = []
    for c in loops:
        c = np.asarray(c)[:, :2]
        for a, b in zip(c[:-1], c[1:]):
            n = max(1, int(np.linalg.norm(b - a) / step))
            out.append(a + np.linspace(0, 1, n, endpoint=False)[:, None] * (b - a))
    return np.vstack(out) if out else np.empty((0, 2))


layer_arg = [a for a in sys.argv if a.startswith("--layers=")]
picks = ([int(x) - 1 for x in layer_arg[0].split("=", 1)[1].split(",")] if layer_arg
         else [len(carved) // 5, len(carved) // 2, (4 * len(carved)) // 5])
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
for col, i in enumerate(picks):
    L = carved[i]
    R = ref.get(L["isovalue"])
    top, zoom = axes[0, col], axes[1, col]
    for ax in (top, zoom):
        if R:
            for c in R["contour_3d"]:
                c = np.asarray(c)
                ax.plot(c[:, 0], c[:, 1], color="0.55", lw=1.6)
        for c in L["contour_3d"]:
            c = np.asarray(c)
            ax.plot(c[:, 0], c[:, 1], color="red", lw=0.8)
        ax.set_aspect("equal")
    # zoom : point creuse a ~epaisseur du contour de reference
    if R:
        r = dense(R["contour_3d"])
        pts = np.vstack([np.asarray(c)[:, :2] for c in L["contour_3d"]])
        d, _ = cKDTree(r).query(pts)
        cand = np.flatnonzero((d > 0.15) & (d < 0.3))
        center = pts[cand[len(cand) // 2]] if len(cand) else pts[0]
        zoom.set_xlim(center[0] - 2.5, center[0] + 2.5)
        zoom.set_ylim(center[1] - 2.5, center[1] + 2.5)
        zoom.set_title(f"zoom 5 x 5 mm : ecart median {np.median(d[cand]) if len(cand) else float('nan'):.2f} mm")
        top.plot(*center, "o", mfc="none", mec="blue", ms=14)
    top.set_title(f"couche {i + 1}/{len(carved)} (vue de dessus)")

first = carved[0]
ax = axes[0, 3]
for p in first["toolpath_3d"]:
    p = np.asarray(p)
    ax.plot(p[:, 0], p[:, 1], lw=0.6)
ax.set_aspect("equal")
ax.set_title("1re couche : concentrique conforme")
# zoom optionnel : --zoom=<couche>:<x>:<y> (ex. un trou a verifier)
zoom_arg = [a for a in sys.argv if a.startswith("--zoom=")]
ax = axes[1, 3]
if zoom_arg:
    n, zx, zy = zoom_arg[0].split("=", 1)[1].split(":")
    L = carved[int(n) - 1]
    R = ref.get(L["isovalue"])
    for c in (R["contour_3d"] if R else []):
        c = np.asarray(c)
        ax.plot(c[:, 0], c[:, 1], color="0.55", lw=2.5)
    for c in L["contour_3d"]:
        c = np.asarray(c)
        ax.plot(c[:, 0], c[:, 1], color="red", lw=1.0)
    ax.set_xlim(float(zx) - 2.5, float(zx) + 2.5)
    ax.set_ylim(float(zy) - 2.5, float(zy) + 2.5)
    ax.set_aspect("equal")
    ax.set_title(f"couche {n} : zoom 5 x 5 mm sur ({float(zx):.1f}, {float(zy):.1f})")
else:
    ax.axis("off")
fig.text(0.5, 0.005, "gris : contour de reference (sans lissage)   |   rouge : corps creuse   |   "
         "ecart gris-rouge = place de la couche de lissage (0,2 mm) ; pas d'ecart = surface non lissee",
         ha="center", fontsize=12)
fig.suptitle(f"P3 : corps creuse sous les zones de lissage -- {model}", fontsize=15)
fig.tight_layout(rect=(0, 0.03, 1, 1))
out = f"{prefix}_carve_compare_pass{k}.png"
fig.savefig(out, dpi=90)
print(out)
