"""
GIF d'ordre d'impression a partir des couches sauvees par bird_to_apex.py.

Rendu autonome (venv de la recherche) : scripts/export_layer_order_gif.py
d'Apex importe le pipeline des son chargement, et le venv d'Apex n'a ni
matplotlib ni Pillow. Memes couleurs que lui : premiere couche deposee en
bleu, derniere en jaune.

Usage : .venv/Scripts/python.exe render_layers_gif.py <couches.pkl> [elev] [azim] [--toolpaths]

--toolpaths : trace les trajectoires completes (parois + remplissage,
cle toolpath_3d) au lieu des seuls contours ; la couche en cours de depot
est surlignee en blanc. Sortie suffixee _toolpaths.
"""
import pickle
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection  # noqa: E402
from PIL import Image  # noqa: E402

FIRST_COLOR = np.array([0x00, 0x00, 0xFF]) / 255.0
LAST_COLOR = np.array([0xFF, 0xF7, 0x00]) / 255.0

args = [a for a in sys.argv[1:] if a != "--toolpaths"]
use_toolpaths = "--toolpaths" in sys.argv
key = "toolpath_3d" if use_toolpaths else "contour_3d"
suffix = "_toolpaths" if use_toolpaths else ""
path = args[0]
elev = float(args[1]) if len(args) > 1 else 20.0
azim = float(args[2]) if len(args) > 2 else -60.0
layers = pickle.load(open(path, "rb"))
stride = max(1, len(layers) // 120)
# pas IMPAIR en mode trajectoires : avec un pas pair, la couche surlignee
# serait toujours de la meme parite, et le remplissage croise-alterne
# semblerait toujours balayer dans le meme sens
if use_toolpaths and stride % 2 == 0:
    stride += 1

pts = np.vstack([c[:, :3] for layer in layers for c in layer[key] if len(c)])
center = (pts.min(0) + pts.max(0)) / 2
radius = float(np.ptp(pts, axis=0).max()) / 2 * 1.1

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(111, projection="3d")
fig.patch.set_facecolor("black")
segments, colors, frames = [], [], []
total = len(layers)
for i, layer in enumerate(layers):
    color = FIRST_COLOR * (1 - i / max(total - 1, 1)) + LAST_COLOR * (i / max(total - 1, 1))
    current = []
    for c in layer[key]:
        if len(c) >= 2:
            c = c[:, :3]
            current.extend(np.stack([c[:-1], c[1:]], axis=1))
    segments.extend(current)
    colors.extend([color] * len(current))
    if i % stride and i != total - 1:
        continue
    ax.clear()
    ax.set_facecolor("black")
    for setter, k in ((ax.set_xlim, 0), (ax.set_ylim, 1), (ax.set_zlim, 2)):
        setter(center[k] - radius, center[k] + radius)
    ax.set_box_aspect([1, 1, 1], zoom=1.7)
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    ax.set_title(f"Couche {i + 1}/{total}", color="white", fontsize=12)
    ax.add_collection3d(Line3DCollection(segments, colors=colors,
                                         linewidths=0.3 if use_toolpaths else 0.6))
    if use_toolpaths and current and i != total - 1:
        ax.add_collection3d(Line3DCollection(current, colors="white", linewidths=0.8))
    fig.canvas.draw()
    frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()))
plt.close(fig)

out = path.replace(".pkl", f"{suffix}.gif")
frames[0].save(out, save_all=True, append_images=frames[1:], duration=80, loop=0)
frames[-1].save(path.replace(".pkl", f"{suffix}_final.png"))
for frac in (0.25, 0.5, 0.75):
    frames[int(frac * (len(frames) - 1))].save(path.replace(".pkl", f"{suffix}_{int(frac * 100)}pct.png"))
print(f"GIF : {out}  ({len(frames)} images)")
