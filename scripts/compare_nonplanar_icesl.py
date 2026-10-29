"""Compare nos surfaces de couches à celles d'un G-code d'IceSL (même modèle, même position).

    python scripts/compare_nonplanar_icesl.py modele.stl icesl.gcode [--center 155,155]

Pour chaque couche IceSL j, mesure z_IceSL - S_j(x,y) aux points extrudés. Sert de repère,
pas de critère de copie : IceSL est propriétaire, seule sa sortie est lue.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import trimesh

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from nonplanar import from_mesh, build  # noqa: E402


def read_layers(path: str):
    layers, cur = [], []
    x = y = z = e = 0.0
    for line in open(path):
        if line.startswith("; <layer>"):
            if cur:
                layers.append(np.array(cur))
            cur = []
            continue
        if not line.startswith("G1"):
            continue
        d = {t[0]: float(t[1:]) for t in line.split(";")[0].split()[1:] if t[0] in "XYZE"}
        nx, ny, nz, ne = d.get("X", x), d.get("Y", y), d.get("Z", z), d.get("E", e)
        if "E" in d and ne > e and ("X" in d or "Y" in d):
            cur.append((nx, ny, nz))
        x, y, z, e = nx, ny, nz, ne
    if cur:
        layers.append(np.array(cur))
    return layers


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stl")
    ap.add_argument("gcode")
    ap.add_argument("--center", default="155,155")
    ap.add_argument("--layer-height", type=float, default=0.2)
    ap.add_argument("--first-layer", type=float, default=0.3)
    a = ap.parse_args()
    mesh = trimesh.load(a.stl)
    mesh.apply_translation([0, 0, -mesh.bounds[0, 2]])
    cx, cy = (float(v) for v in a.center.split(","))
    c = mesh.bounds.mean(axis=0)
    mesh.apply_translation([cx - c[0], cy - c[1], 0])
    stack = build(from_mesh(mesh), a.layer_height, a.first_layer)
    layers = [l for l in read_layers(a.gcode) if len(l) > 100]
    print(f"IceSL : {len(layers)} couches ; nous : {stack.count}")
    n = min(len(layers), stack.count)
    for j in sorted({1, 2, 5, 10, n // 2, n - 1, n}):
        L = layers[j - 1]
        L = L[(L[:, 0] > mesh.bounds[0, 0]) & (L[:, 0] < mesh.bounds[1, 0])]
        ours = stack.z(j, L[:, 0], L[:, 1])
        d = L[:, 2] - ours
        print(f"couche {j:2d} : z IceSL - nôtre  moyenne {d.mean():+.3f}  écart-type {d.std():.3f}  "
              f"min {d.min():+.3f}  max {d.max():+.3f}")


if __name__ == "__main__":
    main()
