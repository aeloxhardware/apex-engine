"""
Essai BORNE de P3 (skin_carve.carved_region, distance exacte adaptative) sur
quelques couches seulement, avant tout tranchage complet : temps, pic memoire
du processus, et controle de la distance du nouveau contour a la peau.

Usage (venv d'Apex) :
  ../../.venv/Scripts/python.exe test_carve_layer.py <modele> <passe> <variante> [n_couches] [epaisseur]
Les couches testees sont reparties parmi les couches TOUCHEES (la mediane
d'abord). Entrees : _cells/_pass0/_pass<k>.npy, _skin_pass<k>.npz,
_deformed_slice_pass<k>.pkl (isovaleurs).
"""
import os
import pickle
import sys
import time

import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import skin_carve  # noqa: E402


def main():
    model, k, variant = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    n_layers = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    thickness = float(sys.argv[5]) if len(sys.argv) > 5 else 0.2
    prefix = os.path.join(HERE, "output_models", "loop_" + model.replace(" ", "_") + variant)
    cells = np.load(f"{prefix}_cells.npy")
    p0 = np.load(f"{prefix}_pass0.npy")
    pk = np.load(f"{prefix}_pass{k}.npy")
    skin = np.load(f"{prefix}_skin_pass{k}.npz")
    isovalues = sorted({L["isovalue"] for L in pickle.load(open(f"{prefix}_deformed_slice_pass{k}.pkl", "rb"))})

    t0 = time.time()
    psi = skin_carve.skin_psi(p0, skin, thickness)
    exact = skin_carve.ExactDistance(skin, cells, pk, p0, thickness)
    touched = [h for h in isovalues if skin_carve.layer_is_touched(cells, pk, psi, h)]
    print(f"psi + preparation : {time.time() - t0:.1f} s ; {len(touched)}/{len(isovalues)} couches touchees ; "
          f"memoire {skin_carve.process_bytes() / 2**20:.0f} Mo")
    order = [len(touched) // 2] + list(np.linspace(0, len(touched) - 1, n_layers).astype(int))
    picks = list(dict.fromkeys(order))[:n_layers]

    from bird_to_apex import barycentric
    for i in picks:
        h = touched[i]
        t0 = time.time()
        region = skin_carve.carved_region(cells, pk, psi, h, exact=exact, step=thickness)
        dt = time.time() - t0
        loops = skin_carve.region_loops(region, h)
        pts = np.vstack([l[::max(1, len(l) // 400)] for l in loops]) if loops else np.empty((0, 3))
        # retour dans la piece reelle par le tetraedre qui contient le point
        z = pk[cells][:, :, 2]
        cut = np.flatnonzero((z.min(axis=1) <= h) & (z.max(axis=1) >= h))
        real = []
        for p in pts:
            with np.errstate(all="ignore"):
                b = barycentric(np.repeat(p[None], len(cut), 0), pk[cells[cut]])
            inside = np.flatnonzero((b >= -1e-6).all(axis=1))
            if len(inside):
                real.append(b[inside[0]] @ p0[cells[cut[inside[0]]]])
        real = np.array(real)
        _c, dist, _t = trimesh.proximity.closest_point(exact.mesh, real) if len(real) else (0, np.array([]), 0)
        near = dist[np.abs(dist - thickness) < 0.5]
        print(f"couche {isovalues.index(h) + 1} (z={h:.2f}) : {dt:.1f} s, {len(loops)} boucles, "
              f"aire {region.area:.1f} mm2, memoire {skin_carve.process_bytes() / 2**20:.0f} Mo ; "
              f"points pres de la peau : {len(near)}/{len(real)}, distance mediane "
              f"{np.median(near) if len(near) else float('nan'):.3f} mm")


if __name__ == "__main__":
    main()
