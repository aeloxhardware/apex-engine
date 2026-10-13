"""Part de surface dont la normale est inatteignable avec C impose par la position, selon le
decalage du modele sur le plateau (point 7 du DT). Pas de pipeline : geometrique pur.
Le modele est pose sur le plateau (z min = 0), centre sur (dx, dy). Au point p (rayon r, azimut a),
la tete ne balaie que le DEMI-plan (R>=0, Z) (B in [0,185] : sin B >= 0). Normale n = (n_r, n_t, n_z) :
si n_r >= 0 l'erreur est asin(|n_t|) ; si n_r < 0 (penche vers le centre) la direction atteignable la
plus proche est +-Z, erreur = acos(|n_z|).
Usage : measure_offset_reach.py <stl> ..."""
import sys
import numpy as np
import trimesh

OFFSETS = (0.0, 30.0, 60.0, 90.0, 120.0)
DIRS = [(1, 0), (-1, 0), (0, 1), (0, -1)]
for path in sys.argv[1:]:
    m = trimesh.load(path, force="mesh")
    c = m.triangles_center
    c = c - np.array([(m.bounds[0][0] + m.bounds[1][0]) / 2, (m.bounds[0][1] + m.bounds[1][1]) / 2, m.bounds[0][2]])
    n, area = m.face_normals, m.area_faces
    steep = np.abs(n[:, 2]) < 0.7     # flancs (hors faces quasi horizontales)
    print(f"{path}  ({len(area)} faces, flancs {100*area[steep].sum()/area.sum():.0f} % de l'aire)")
    for off in OFFSETS:
        rows = []
        for dx, dy in DIRS if off else [(0, 0)]:
            p = c[:, :2] + np.array([dx * off, dy * off])
            r = np.hypot(p[:, 0], p[:, 1])
            tx, ty = np.where(r > 1e-9, -p[:, 1] / np.maximum(r, 1e-9), 0), np.where(r > 1e-9, p[:, 0] / np.maximum(r, 1e-9), 0)
            rx, ry = np.where(r > 1e-9, p[:, 0] / np.maximum(r, 1e-9), 1.0), np.where(r > 1e-9, p[:, 1] / np.maximum(r, 1e-9), 0.0)
            nr, nt = n[:, 0] * rx + n[:, 1] * ry, n[:, 0] * tx + n[:, 1] * ty
            err = np.degrees(np.where(nr >= 0, np.arcsin(np.clip(np.abs(nt), 0, 1)), np.arccos(np.clip(np.abs(n[:, 2]), 0, 1))))
            w = area[steep]
            rows.append([100 * w[err[steep] > t].sum() / w.sum() for t in (15, 45)])
        a = np.array(rows)
        print(f"  decalage {off:5.0f} mm : flancs > 15 deg : moy {a[:,0].mean():5.1f} % (min {a[:,0].min():5.1f}, max {a[:,0].max():5.1f}) ; > 45 deg : moy {a[:,1].mean():5.1f} % (min {a[:,1].min():5.1f}, max {a[:,1].max():5.1f})")
