"""
Profil d'inclinaison d'un G-code Apex : B, pointe et point le plus bas de la
VRAIE tete, par tranche de B, sous plusieurs conventions pivot/pointe.

Sert a confronter le modele de tete de l'export (file_io/gcode_export.py) a
la previsualisation UX de l'utilisateur, qui fait reference : pour chaque
convention, combien de mouvements la previsualisation verrait sous le plateau
et de combien la pointe serait deportee.

Usage (venv d'Apex) :
  .venv/Scripts/python.exe scripts/gcode_tilt_profile.py <fichier.gcode> [--png=<sortie.png>]
      [--pivot=45] [--no-head]

Conventions testees :
  apex         : X = R + L sin B, Z = z + L (cos B - 1)  (export actuel)
  pivot-z      : Z firmware = hauteur du pivot, Z = z + L cos B
  sign-         : deport radial du pivot de signe oppose
  raw          : X/Z lus comme la pointe (aucune compensation)
"""
import os
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

from file_io.gcode_export import firmware_to_machine, parse_gcode  # noqa: E402

BANDS = [(0, 5), (5, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90), (90, 105), (105, 186)]


def tips_for(moves, pivot, convention):
    X, Y, Z, B = moves.T
    if convention == "raw":
        a = np.radians(Y)
        return np.column_stack([X * np.cos(a), X * np.sin(a), Z])
    kw = {"apex": {}, "pivot-z": {"z_reference": "pivot"}, "sign-": {"tilt_sign": -1.0}}[convention]
    return np.array([firmware_to_machine(x, y, z, b, pivot, **kw) for x, y, z, b in moves])


def head_zlow_table(pivot_grid=0.5):
    import main as apex_main
    from pipeline import step9_collision_sim as s9
    machine = apex_main._default_machine_profile()
    if machine.nozzle_model is None:
        return None
    sim = s9.CollisionSimParams(path_sampling_step=2.0, safety_margin=1.0)
    head = s9._filter_nozzle_contact_zone(s9._surface_sample_points(machine.nozzle_model, sim.toolhead_sample_spacing),
                                          sim.nozzle_contact_clearance)
    grid = np.arange(0.0, 185.0 + 1e-9, pivot_grid)
    zl = np.array([s9._positioned_toolhead_vertices(head, np.zeros(3), [0.0, b], machine.rotation_axes)[:, 2].min()
                   for b in grid])
    return grid, zl


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] if "=" in a else "" for a in sys.argv[1:] if a.startswith("--")}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = args[0]
    pivot = float(flags.get("--pivot", 45.0))
    t0 = time.time()
    raw = parse_gcode(open(path).read().splitlines())
    moves = np.array([m[:4] for m in raw], dtype=float)
    ext = np.array([m[4] for m in raw])
    B = moves[:, 3]
    print(f"{path}\n  {len(moves)} mouvements ({ext.sum()} en extrusion), lus en {time.time() - t0:.0f} s")
    print(f"  firmware : X [{moves[:, 0].min():.2f}, {moves[:, 0].max():.2f}]  "
          f"Y [{moves[:, 1].min():.1f}, {moves[:, 1].max():.1f}] deg "
          f"({(moves[:, 1].max() - moves[:, 1].min()) / 360:.1f} tours)  "
          f"Z [{moves[:, 2].min():.2f}, {moves[:, 2].max():.2f}]  B [{B.min():.1f}, {B.max():.1f}]")
    zl = None if "--no-head" in flags else head_zlow_table()
    results = {}
    for conv in ("apex", "pivot-z", "sign-", "raw"):
        tips = tips_for(moves, pivot, conv)
        results[conv] = tips
        low = None
        if zl is not None:
            grid, table = zl
            low = tips[:, 2] + table[np.clip(np.round(B / (grid[1] - grid[0])).astype(int), 0, len(grid) - 1)]
        print(f"\n== convention {conv} ==")
        print(f"{'B':>9} {'mvts':>8} {'pointe<0':>9} {'pointe<0,1':>10} {'min pointe':>10}"
              + (f" {'tete<0':>7} {'min tete':>9}" if low is not None else ""))
        for lo, hi in BANDS:
            m = (B >= lo) & (B < hi)
            if not m.any():
                continue
            row = (f"{lo:>3}-{hi:<4}  {m.sum():>8} {(tips[m, 2] < 0).sum():>9} {(tips[m, 2] < 0.1).sum():>10} "
                   f"{tips[m, 2].min():>10.3f}")
            if low is not None:
                row += f" {(low[m] < 0).sum():>7} {low[m].min():>9.2f}"
            print(row)
        print(f"  total pointe < 0 : {(tips[:, 2] < 0).sum()} ({100 * (tips[:, 2] < 0).mean():.2f} %)")
    png = flags.get("--png")
    if png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = len(B)
        step = max(1, n // 200000)
        x = np.arange(n)[::step]
        fig, ax = plt.subplots(3, 1, figsize=(16, 10), sharex=True)
        ax[0].scatter(x, B[::step], s=0.3, c=np.where(ext, "tab:blue", "0.7")[::step])
        ax[0].set_ylabel("B (deg)")
        ax[0].set_title(os.path.basename(path) + " : bleu = extrusion, gris = trajets")
        for conv, col in (("apex", "tab:blue"), ("raw", "tab:red")):
            ax[1].scatter(x, results[conv][::step, 2], s=0.3, c=col, label=conv)
        ax[1].axhline(0, color="k", lw=0.8)
        ax[1].set_ylabel("z pointe (mm)")
        ax[1].legend(markerscale=20)
        ax[2].scatter(x, moves[::step, 2], s=0.3)
        ax[2].axhline(0, color="k", lw=0.8)
        ax[2].set_ylabel("Z firmware")
        ax[2].set_xlabel("mouvement")
        plt.tight_layout()
        plt.savefig(png, dpi=110)
        print(f"\ngraphe : {png}")
    print(f"duree {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
