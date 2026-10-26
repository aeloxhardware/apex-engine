"""Un seul trou de rayon croissant sur un site fin ; mesure ce que TetraNet garde du volume.
Usage : python sweep_holes.py <stl> <report.json> <ref_prefix> <cli> [h]"""
import json, subprocess, sys, tempfile, os
import numpy as np, trimesh
sys.path.insert(0, os.path.dirname(__file__))
stl, rep, ref_prefix, cli = sys.argv[1:5]; h = sys.argv[5] if len(sys.argv) > 5 else "1.0"
site = np.array(json.load(open(rep))["sites"][int(os.environ.get("SITE", "0"))]["p"])
m = trimesh.load(stl, process=False); m.merge_vertices()
V, F = np.asarray(m.vertices), np.asarray(m.faces); cen = V[F].mean(1)
def vol(prefix):
    n = np.fromfile(prefix + ".nodes", "<f8").reshape(-1, 3); t = np.fromfile(prefix + ".tets", "<i4").reshape(-1, 4); p = n[t]
    return np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]).sum() / 6
ref = vol(ref_prefix)
import struct
for r in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0):
    keep = np.linalg.norm(cen - site, axis=1) >= r
    T = V[F[keep]].astype("<f4"); path = os.path.join(tempfile.gettempdir(), f"hole_{r}.stl")
    with open(path, "wb") as f:
        f.write(b"x".ljust(80, b" ")); f.write(struct.pack("<I", int(keep.sum())))
        rec = np.zeros(int(keep.sum()), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]); rec["v"] = T; f.write(rec.tobytes())
    pref = path[:-4]
    out = subprocess.run([cli, path, "--cell", h, "--alpha", "0.5", "--out", pref], capture_output=True, text=True)
    j = json.loads(out.stdout.strip().splitlines()[-1])
    print(f"rayon {r:3.1f} mm ({int((~keep).sum()):5d} faces retirées) : volume gardé {vol(pref)/ref:6.3f}  conflits de signe {j['sign_conflicts']}", flush=True)
