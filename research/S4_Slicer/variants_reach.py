"""Part des normales du corps Bird atteignables selon la variante de machine
(Benchy passe 4, 1 polyligne sur 10). Ensembles de directions atteignables,
dans le repere ou le rayon du point est -Y (C impose) : base = R_X(b) z,
b in [0,185] ; +D = R_X(b) R_Y(d) z, d in [-90,90] ; +Y (C libre) : sphere
entiere avec C tour complet et B in [0,185] (limite #18), donc ecart nul."""
import os, pickle, sys
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from full_step89 import locate, tet_gradients

# usage : variants_reach.py [<modele> <passe> <variante>] ; sans print_order (P5 non fait),
# les points testes sont les centres des tetraedres (normale = gradient de phi).
_a = [a for a in sys.argv[1:] if not a.startswith("--")]
model, k, variant = (_a[0], _a[1], _a[2]) if len(_a) > 2 else (
    "tilted benchy undertriangulated_repaired", "4", "_pin")
pre = "output_models/loop_" + model.replace(" ", "_") + variant
cells, p0, pk = (np.load(f"{pre}_{s}.npy") for s in ("cells", "pass0", f"pass{k}"))
g = tet_gradients(p0, cells, pk[:, 2]); g /= np.linalg.norm(g, axis=1, keepdims=True)
tree = cKDTree(p0[cells].mean(1))
if os.path.exists(f"{pre}_print_order_pass{k}.pkl") and "--centres" not in sys.argv:
    order = pickle.load(open(f"{pre}_print_order_pass{k}.pkl", "rb"))["passes"]
    polys = [p for ps in order if ps["kind"] == "body" for p in ps["points"]][::10]
    P = np.vstack(polys); N = g[locate(P, p0, cells, tree)]
else:
    P = p0[cells].mean(1); N = g
    print("(centres des tetraedres : pas d'ordre d'impression pour cette variante)")
a = -np.pi / 2 - np.arctan2(P[:, 1], P[:, 0])
c, s = np.cos(a), np.sin(a)
L = np.column_stack([c * N[:, 0] - s * N[:, 1], s * N[:, 0] + c * N[:, 1], N[:, 2]])

def rx(b):
    return np.stack([np.zeros_like(b), -np.sin(b), np.cos(b)], -1)
b = np.radians(np.arange(0, 185.01, 0.5))
base = rx(b)
B, D = np.meshgrid(b, np.radians(np.arange(-180, 180.01, 1.0)))
B, D = B.ravel(), D.ravel()
v = np.stack([np.sin(D), np.zeros_like(D), np.cos(D)], -1)          # R_Y(d) z
plusD = np.stack([v[:, 0], np.cos(B) * v[:, 1] - np.sin(B) * v[:, 2],
                  np.sin(B) * v[:, 1] + np.cos(B) * v[:, 2]], -1)    # puis R_X(b)
r = P[:, :2] / np.maximum(np.linalg.norm(P[:, :2], axis=1, keepdims=True), 1e-9)
print(f"normales vers le centre (radial < -0,05) : {np.mean((N[:, :2] * r).sum(1) < -0.05):.1%}")
for name, S in (("actuelle (R C Z B)", base), ("+ D (sans limite)", plusD)):
    dot = cKDTree(S).query(L)[0]
    err = np.degrees(2 * np.arcsin(np.clip(dot / 2, 0, 1)))
    print(f"{name:20s}: ecart median {np.median(err):5.1f} deg, > 5 deg {np.mean(err > 5):6.1%}, "
          f"> 45 deg {np.mean(err > 45):5.1%}, max {err.max():5.1f}")
print(f"{'+ Y, + Y + D':20s}: ecart 0 partout (sphere complete) ; {len(P)} points testes")
