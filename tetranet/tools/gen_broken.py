"""Fabrique des versions volontairement cassées d'un STL (un type de défaut par fichier, plus un
fichier global), aux endroits les plus fins de la pièce, pour éprouver un mailleur.

Usage : python gen_broken.py <modele.stl> <dossier_sortie>
Types : holes, flipped, selfint, nonmanifold, cracks, degenerate, global.
"""
import json, os, struct, sys
import numpy as np
import trimesh

src, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
rng = np.random.default_rng(7)

m = trimesh.load(src, process=False)
m.merge_vertices()
V0, F0 = np.asarray(m.vertices, float), np.asarray(m.faces, int)
print("modèle :", len(V0), "sommets,", len(F0), "faces,", "étanche" if m.is_watertight else "non étanche")

def write_stl(path, V, F):
    T = V[F].astype("<f4")
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]).astype("<f4")
    with open(path, "wb") as f:
        f.write(b"tetranet broken".ljust(80, b" ")); f.write(struct.pack("<I", len(F)))
        rec = np.zeros(len(F), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        rec["n"] = n; rec["v"] = T
        f.write(rec.tobytes())

# ---- endroits les plus fins : épaisseur le long de la normale entrante
mesh = trimesh.Trimesh(V0, F0, process=False)
cent = mesh.triangles_center; nrm = mesh.face_normals
samp = rng.choice(len(F0), size=min(8000, len(F0)), replace=False)
orig = cent[samp] - nrm[samp] * 1e-3
loc, ray, _ = mesh.ray.intersects_location(orig, -nrm[samp], multiple_hits=False)
thick = np.full(len(samp), np.inf)
thick[ray] = np.linalg.norm(loc - orig[ray], axis=1)
order = np.argsort(thick)
sites = []
for k in order:
    if thick[k] < 0.3: continue            # artefact (auto-impact)
    p = cent[samp[k]]
    if all(np.linalg.norm(p - s["p"]) > 9 for s in sites):
        sites.append(dict(p=p, n=nrm[samp[k]], thickness=float(thick[k])))
    if len(sites) == 8: break
print("endroits les plus fins :", [(np.round(s["p"], 1).tolist(), round(s["thickness"], 2)) for s in sites])

def ball(p, r): return np.linalg.norm(cent - p, axis=1) < r

def holes(V, F, sizes=(1, 2, 4, 6)):
    keep = np.ones(len(F), bool)
    for s, r in zip(sites, sizes):
        keep &= ~(np.linalg.norm(V[F].mean(1) - s["p"], axis=1) < r)
    return V, F[keep]

def flipped(V, F, radii=(3, 5, 8, 12)):
    F = F.copy()
    for s, r in zip(sites, radii):
        sel = np.linalg.norm(V[F].mean(1) - s["p"], axis=1) < r
        F[sel] = F[sel][:, ::-1]
    return V, F

def selfint(V, F, depth=2.5, radii=(4, 5, 6, 7)):
    V = V.copy()
    vn = mesh.vertex_normals
    for s, r in zip(sites, radii):
        sel = np.linalg.norm(V - s["p"], axis=1) < r
        w = np.clip(1 - np.linalg.norm(V[sel] - s["p"], axis=1) / r, 0, 1)[:, None]
        V[sel] -= vn[sel] * depth * (s["thickness"] / max(s["thickness"], 1.0) + 1.0) * w
    return V, F

def nonmanifold(V, F, fins=120, dup=800):
    V = list(V); F = list(F)
    for s in sites[:4]:
        idx = np.where(np.linalg.norm(np.asarray(V)[np.asarray(F)].mean(1) - s["p"], axis=1) < 6)[0]
        pick = rng.choice(idx, size=min(fins // 4, len(idx)), replace=False)
        for fi in pick:
            a, b, c = F[fi]
            n = np.cross(np.asarray(V[b]) - V[a], np.asarray(V[c]) - V[a]); n /= np.linalg.norm(n) + 1e-12
            V.append(((np.asarray(V[a]) + V[b]) / 2 + n * 1.0).tolist()); F.append([a, b, len(V) - 1])  # ailette : 3e face sur l'arête ab
        pick = rng.choice(idx, size=min(dup // 4, len(idx)), replace=False)
        F.extend([F[i] for i in pick])                                   # faces dupliquées
    return np.asarray(V), np.asarray(F)

def cracks(V, F, gaps=(0.1, 0.3, 0.6, 1.0)):
    V = V.copy(); F = F.copy(); Vl = list(V)
    for s, g in zip(sites, gaps):
        # une fente : les faces d'un côté d'un plan passant par le site sont décalées de g
        n = s["n"]; t = np.cross(n, [0, 0, 1.0]); t /= np.linalg.norm(t) + 1e-12
        cen = V[F].mean(1)
        sel = (np.linalg.norm(cen - s["p"], axis=1) < 7) & (((cen - s["p"]) @ t) > 0)
        for fi in np.where(sel)[0]:
            new = []
            for v in F[fi]:
                Vl.append((np.asarray(V[v]) + t * g).tolist()); new.append(len(Vl) - 1)
            F[fi] = new
    return np.asarray(Vl), F

def degenerate(V, F):
    V = V.copy(); F = F.copy()
    V += rng.normal(0, 0.02, V.shape)                                    # bruit de sommets
    pick = rng.choice(len(F), size=2000, replace=False)                  # arêtes effondrées
    for fi in pick: V[F[fi][1]] = V[F[fi][0]]
    return V, F

ALL = {"holes": holes, "flipped": flipped, "selfint": selfint, "nonmanifold": nonmanifold, "cracks": cracks, "degenerate": degenerate}
report = {"sites": [dict(p=s["p"].tolist(), thickness=s["thickness"]) for s in sites]}
for name, fn in ALL.items():
    V, F = fn(V0, F0)
    write_stl(os.path.join(outdir, f"{name}.stl"), np.asarray(V), np.asarray(F)); report[name] = dict(faces=int(len(F)))
V, F = V0, F0
for name in ("selfint", "holes", "flipped", "nonmanifold", "cracks", "degenerate"):
    V, F = ALL[name](V, F)
write_stl(os.path.join(outdir, "global.stl"), np.asarray(V), np.asarray(F)); report["global"] = dict(faces=int(len(F)))
json.dump(report, open(os.path.join(outdir, "report.json"), "w"), indent=1)
print("écrit :", sorted(os.listdir(outdir)))
