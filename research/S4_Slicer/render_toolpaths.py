"""
Images de la trajectoire de l'etape 6 (parois + remplissage) sur les couches
sauvees par bird_to_apex.py.

  1. une planche de 6 couches vues dans leur repere local (u, v) : boucles
     parois et contours en bleu, remplissage (zigzags) en orange ;
  2. la piece en 3D a plusieurs stades d'avancement, trajectoire complete
     de toutes les couches deja deposees.

Usage : .venv/Scripts/python.exe render_toolpaths.py <couches.pkl>
"""
import pickle
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection  # noqa: E402

FIRST_COLOR = np.array([0x00, 0x00, 0xFF]) / 255.0
LAST_COLOR = np.array([0xFF, 0xF7, 0x00]) / 255.0

path = sys.argv[1]
layers = pickle.load(open(path, "rb"))
total = len(layers)


def is_closed(poly):
    return len(poly) > 2 and np.linalg.norm(poly[0, :2] - poly[-1, :2]) < 1e-6


def is_infill(poly):
    """Remplissage = polyligne ouverte qui porte au moins une ligne de
    balayage (segment de plus de 0,5 mm exactement selon u ou v). Le contour
    d'origine est lui aussi ouvert : l'ouverture seule ne suffit pas."""
    if is_closed(poly) or len(poly) < 2:
        return False
    d = np.diff(poly[:, :2], axis=0)
    sweep = (np.abs(d[:, 0]) < 1e-9) | (np.abs(d[:, 1]) < 1e-9)
    return bool((sweep & (np.hypot(d[:, 0], d[:, 1]) > 0.5)).any())


# 1. planche de couches en repere local
fracs = [0.10, 0.30, 0.50, 0.65, 0.80, 0.95]
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for ax, frac in zip(axes.ravel(), fracs):
    i = min(int(frac * total), total - 1)
    walls = infill = 0
    for poly in layers[i]["toolpath_2d"]:
        fill = is_infill(poly)
        walls += not fill
        infill += fill
        ax.plot(poly[:, 0], poly[:, 1], color="tab:orange" if fill else "tab:blue", lw=0.6)
    ax.set_aspect("equal")
    ax.set_title(f"couche {i + 1}/{total} ({frac:.0%}) : {walls} parois, {infill} zigzags", fontsize=9)
    ax.tick_params(labelsize=7)
fig.suptitle("Trajectoire de l'etape 6, repere local de chaque couche (mm) -- bleu : parois, orange : remplissage")
fig.tight_layout()
out1 = path.replace(".pkl", "_toolpaths_layers.png")
fig.savefig(out1, dpi=110)
plt.close(fig)

# 2. la piece en 3D a plusieurs stades
pts = np.vstack([t[:, :3] for layer in layers for t in layer["toolpath_3d"] if len(t)])
center = (pts.min(0) + pts.max(0)) / 2
radius = float(np.ptp(pts, axis=0).max()) / 2 * 0.75
fig = plt.figure(figsize=(16, 5.5))
fig.patch.set_facecolor("black")
for k, frac in enumerate((0.35, 0.6, 0.8, 1.0)):
    ax = fig.add_subplot(1, 4, k + 1, projection="3d")
    ax.set_facecolor("black")
    upto = max(1, int(frac * total))
    segs, cols = [], []
    for i in range(upto):
        color = FIRST_COLOR * (1 - i / (total - 1)) + LAST_COLOR * (i / (total - 1))
        for t in layers[i]["toolpath_3d"]:
            if len(t) >= 2:
                segs.extend(np.stack([t[:-1, :3], t[1:, :3]], axis=1))
                cols.extend([color] * (len(t) - 1))
    ax.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=0.25))
    for setter, j in ((ax.set_xlim, 0), (ax.set_ylim, 1), (ax.set_zlim, 2)):
        setter(center[j] - radius, center[j] + radius)
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=15, azim=-60)
    ax.set_axis_off()
    ax.set_title(f"{frac:.0%} (couche {upto}/{total})", color="white", fontsize=10)
fig.tight_layout()
out2 = path.replace(".pkl", "_toolpaths_3d.png")
fig.savefig(out2, dpi=110, facecolor="black")
plt.close(fig)
print(out1)
print(out2)
