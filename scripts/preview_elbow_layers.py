"""Apercu des couches du tube coude, en couleurs alternees (une couleur sur deux).

Usage : python scripts/preview_elbow_layers.py sortie.png [options de gen_elbow_tube_gcode.py]
Rejoue build() du generateur (meme options), sans relire le G-code. Vue de cote
(X,Z) complete + agrandissement de la zone de transition, et vue de dessus.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_elbow_tube_gcode as G


def main() -> None:
    args = sys.argv[1:]
    out = args.pop(0)
    sys.argv = ["gen"] + args
    # memes options que le generateur
    for flag, name in (("--planar-cut", "PLANAR_CUT"), ("--full-infill", "FULL_INFILL")):
        if flag in args:
            setattr(G, name, True)
    if "--pivot" in args:
        G.PIVOT_ABOVE_TIP_MM = float(args[args.index("--pivot") + 1])
    if "--cut-offset" in args:
        G.CUT_OFFSET = float(args[args.index("--cut-offset") + 1])
    orig = G.GCode.__init__

    def init(self, *a, **k):
        orig(self, *a, **k)
        self.trace = []
    G.GCode.__init__ = init
    g = G.build()
    segs = {}
    prev = None
    for x, y, z, ext, layer in g.trace:
        if ext and prev is not None:
            segs.setdefault(layer, []).append(((prev[0], prev[2]), (x, z), (prev[0], prev[1]), (x, y)))
        prev = (x, y, z)
    from matplotlib.collections import LineCollection
    cols = ("#d9822b", "#2b6fd9")
    fig, axs = plt.subplots(1, 3, figsize=(18, 6.5))
    for ax, title in zip(axs, ("vue de cote (X, Z)", "zoom sur la transition plane -> inclinee", "vue de dessus (X, Y)")):
        ax.set_title(title)
        ax.set_aspect("equal")
    side = ([], [])
    top = ([], [])
    for layer, ss in sorted(segs.items()):
        for (a, b, ta, tb) in ss:
            side[layer % 2].append((a, b))
            top[layer % 2].append((ta, tb))
    for par in (0, 1):
        for ax in axs[:2]:
            ax.add_collection(LineCollection(side[par], colors=cols[par], linewidths=0.6))
        axs[2].add_collection(LineCollection(top[par], colors=cols[par], linewidths=0.3))
    axs[0].autoscale()
    axs[2].autoscale()
    axs[1].set_xlim(-6, 12)
    axs[1].set_ylim(17, 28)
    for ln in list(axs[1].collections):
        ln.set_linewidth(2.2)
    plt.tight_layout()
    plt.savefig(out, dpi=130)
    print(out, len(segs), "couches")


if __name__ == "__main__":
    main()
