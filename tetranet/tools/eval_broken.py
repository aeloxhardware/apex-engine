"""Mesure TetraNet sur des STL cassés, contre le maillage de la pièce saine.

Usage : python eval_broken.py <ref_prefix> <dossier_stl> <cli> [h]
  <ref_prefix> : sortie du CLI sur la pièce saine (<prefix>.nodes / .tets)
"""
import glob, json, os, subprocess, sys, tempfile
import numpy as np

ref_prefix, folder, cli = sys.argv[1:4]
h = sys.argv[4] if len(sys.argv) > 4 else "1.0"
VOX = 2.0

def load(prefix):
    n = np.fromfile(prefix + ".nodes", "<f8").reshape(-1, 3); t = np.fromfile(prefix + ".tets", "<i4").reshape(-1, 4)
    p = n[t]; vol = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]) / 6
    return p.mean(1), vol

def field(c, vol):
    """Volume de matière par voxel de 2 mm (insensible à un décalage du réseau de quelques 0,1 mm)."""
    k = np.floor(c / VOX).astype(np.int64)
    d = {}
    for key, v in zip(map(tuple, k), vol): d[key] = d.get(key, 0.0) + v
    return d

def overlap(a, b):
    keys = set(a) | set(b)
    num = sum(min(a.get(k, 0.0), b.get(k, 0.0)) for k in keys); den = sum(max(a.get(k, 0.0), b.get(k, 0.0)) for k in keys)
    miss = sum(max(a.get(k, 0.0) - b.get(k, 0.0), 0.0) for k in keys); extra = sum(max(b.get(k, 0.0) - a.get(k, 0.0), 0.0) for k in keys)
    return num / den, miss, extra

rc, rv = load(ref_prefix); ref = field(rc, rv); refvol = rv.sum()
print(f"référence : {len(rv)} tétraèdres, volume {refvol:.0f} mm3", flush=True)
rows = []
for stl in sorted(glob.glob(os.path.join(folder, "*.stl"))):
    name = os.path.basename(stl)[:-4]; pref = os.path.join(tempfile.gettempdir(), "tn_eval_" + name)
    out = subprocess.run([cli, stl, "--cell", h, "--alpha", "0.5", "--out", pref], capture_output=True, text=True)
    try: r = json.loads(out.stdout.strip().splitlines()[-1])
    except Exception: rows.append(dict(name=name, error=(out.stderr or out.stdout)[-200:])); print(rows[-1]); continue
    if r["tets"] == 0: rows.append(dict(name=name, status=r["status"], reason=r["reason"])); print(rows[-1]); continue
    c, v = load(pref); occ = field(c, v)
    iou, miss, extra = overlap(ref, occ)
    row = dict(name=name, status=r["status"], seconds=round(r["seconds"], 2), tets=r["tets"], sign_conflicts=r["sign_conflicts"],
               vol_vs_ref=round(float(v.sum() / refvol), 4), iou=round(iou, 4),
               missing_pct=round(100 * miss / refvol, 2), extra_pct=round(100 * extra / refvol, 2),
               nonmanifold=r["nonmanifold_edges"], min_dihedral=round(r["min_dihedral_deg"], 2), skin_max=round(r["skin_max_error"], 3),
               dropped=r["dropped_degenerate"], neg_vol=int((v <= 0).sum()))
    rows.append(row); print(json.dumps(row), flush=True)
json.dump(rows, open(os.path.join(folder, "eval.json"), "w"), indent=1)
