"""
P4 (etape 4e) -- trajectoires des couches de lissage non planaire pour TOUTES
les zones d'un modele, a partir de la sortie de P2-bis (_skin_pass<k>.npz).

Par zone :
  - maillage de calcul (surface_paths.ZoneMesh), motif choisi (motif general,
    ou motif par zone), angle des lignes (general ou par zone) ;
  - chaque trajectoire descendue de t/2 sous la surface le long de la normale
    aux sommets (normales des faces de la zone moyennees par l'aire : pas de
    fente ouverte/fermee aux aretes des facettes) -- l'axe du cordon est au
    milieu de la couche d'epaisseur t creusee par P3 ;
  - reechantillonnage (segments <= 0,5 mm, points a moins de 0,05 mm fusionnes),
    normale de surface en chaque point (orientation de la tete pour P5),
    facteur d'extrusion 1 (0 plus tard pour l'ironing) ;
  - controles : distance des points a la surface = t/2 +- 0,05 mm, couverture,
    duree. Une zone en echec est SIGNALEE (erreur gardee), jamais ignoree.
Sortie : loop_<modele>_skin_paths_pass<k>.pkl.

Usage (venv d'Apex) :
  ../../.venv/Scripts/python.exe -u skin_paths.py <modele> <passe> [--pattern=rectilinear]
      [--zone-pattern=12:concentric,40:hilbert] [--line-angle=45] [--zone-angle=12:0]
      [--overlap=0.2] [--thickness=0.2] [--concentric-order=outside_in] [--zones=3,7]
"""
import safe_io  # noqa: F401 -- ecritures robustes aux verrous transitoires
import os
import pickle
import sys
import time

import numpy as np
import trimesh
from scipy.spatial import cKDTree

import surface_paths as sp
from skin_carve import check_memory

HERE = os.path.dirname(os.path.abspath(__file__))
MAX_SEGMENT = 0.5          # mm, longueur maximale d'un segment de trajectoire
MIN_SEGMENT = 0.05         # mm, points plus proches fusionnes
DEPTH_TOLERANCE = 0.05     # mm, ecart admis a t/2 sous la surface


def parse_map(text, cast):
    out = {}
    for item in filter(None, (text or "").split(",")):
        z, v = item.split(":")
        out[int(z)] = cast(v)
    return out


def vertex_normals(zm, face_normals):
    """Normale par sommet du maillage de calcul : normales des faces
    d'origine (via zm.parent) ponderees par l'aire des faces raffinees."""
    tri = zm.V[zm.F]
    area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    fn = face_normals[zm.parent] * area[:, None]
    vn = np.zeros_like(zm.V)
    for c in range(3):
        np.add.at(vn, zm.F[:, c], fn)
    return vn / np.maximum(np.linalg.norm(vn, axis=1), 1e-12)[:, None]


def resample(path):
    """Segments <= MAX_SEGMENT ; points a moins de MIN_SEGMENT fusionnes."""
    keep = [path[0]]
    for p in path[1:-1]:
        if np.linalg.norm(p - keep[-1]) >= MIN_SEGMENT:
            keep.append(p)
    if len(keep) > 1 and np.linalg.norm(path[-1] - keep[-1]) < MIN_SEGMENT:
        keep[-1] = path[-1]                         # l'extremite reelle est conservee
    else:
        keep.append(path[-1])
    keep = np.asarray(keep)
    out = [keep[:1]]
    for a, b in zip(keep[:-1], keep[1:]):
        n = max(1, int(np.ceil(np.linalg.norm(b - a) / MAX_SEGMENT)))
        out.append(a + np.linspace(0, 1, n + 1)[1:, None] * (b - a))
    return np.vstack(out)


def zone_paths(zm, face_normals, pattern, angle, overlap, thickness, order, width, edge_rings=1,
               edge_overlap=sp.EDGE_OVERLAP, remedy=None):
    """Trajectoires d'une zone. Pour Hilbert : part de surface etiree par le
    depliage ; au-dela de sp.HILBERT_WARN_FRACTION sans remede choisi, la zone
    est tracee telle quelle (A) et un AVERTISSEMENT structure est rendu (choix
    A / B / C expliques, pour l'interface) ; remede B = decoupe, C = grille
    resserree."""
    d_edge = zm.distance_to_boundary()
    warning, kw = None, {}
    fn = pattern
    if pattern == "hilbert":
        uv, frac = sp.hilbert_stretch(zm)
        if remedy == "B":
            fn, kw = sp.hilbert_split, {"piece_rings": edge_rings, "piece_overlap": edge_overlap}
            edge_rings = 0                                  # cordons poses morceau par morceau
        elif remedy == "C":
            fn = sp.hilbert_geodesic
        else:
            kw = {"uv": uv}
            if frac > sp.HILBERT_WARN_FRACTION and remedy != "A":
                warning = {"type": "hilbert_stretch", "stretched_area_fraction": frac,
                           "tolerance": sp.HILBERT_STRETCH_TOLERANCE,
                           "message": (f"Surface trop bombee pour Hilbert : {frac:.0%} de la zone s'etire de "
                                       f"plus de {sp.HILBERT_STRETCH_TOLERANCE:.0%} une fois mise a plat."),
                           "options": sp.HILBERT_OPTIONS, "applied": "A"}
    res = sp.with_edge_ring(zm, fn, width, overlap, edge_rings, edge_overlap, d_edge=d_edge,
                            angle_deg=angle, order=order, **kw)
    return (*_finish(zm, face_normals, res, thickness), d_edge, warning)


def _finish(zm, face_normals, res, thickness):
    vn = vertex_normals(zm, face_normals)
    tree = cKDTree(zm.V)
    paths, normals = [], []
    for p in res.paths:
        if len(p) < 2:
            continue
        p = resample(p)
        n = vn[tree.query(p)[1]]
        paths.append(p - n * thickness / 2)
        normals.append(n)
    return paths, normals


def main():
    model, k = sys.argv[1], int(sys.argv[2])
    flag = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[3:] if "=" in a}
    general = flag.get("--pattern", "rectilinear")
    zone_pattern = parse_map(flag.get("--zone-pattern"), str)
    line_angle = float(flag.get("--line-angle", 45.0))
    zone_angle = parse_map(flag.get("--zone-angle"), float)
    overlap = float(flag.get("--overlap", 0.2))
    thickness = float(flag.get("--thickness", 0.2))
    order = flag.get("--concentric-order", "outside_in")
    only = {int(z) for z in flag["--zones"].split(",")} if "--zones" in flag else None
    edge_rings = int(flag.get("--edge-rings", 1))
    edge_overlap = float(flag.get("--edge-overlap", sp.EDGE_OVERLAP))
    remedy = flag.get("--hilbert-remedy")
    zone_remedy = parse_map(flag.get("--zone-hilbert-remedy"), str)
    for name in [general, *zone_pattern.values()]:
        if name not in sp.PATTERNS:
            sys.exit(f"motif inconnu : {name} (disponibles : {', '.join(sp.PATTERNS)})")
    for r in [remedy, *zone_remedy.values()]:
        if r not in (None, "A", "B", "C"):
            sys.exit(f"remede Hilbert inconnu : {r} (A, B ou C)")

    prefix = os.path.join(HERE, "output_models", f"loop_{model}")
    skin = np.load(f"{prefix}_skin_pass{k}.npz")
    normals = skin["normals"]
    width = float(skin["nozzle_width"])
    zone, zone_moment = skin["zone"], {z: int(m) for z, m in enumerate(skin["zone_moment"])}
    merges = [tuple(r) for r in skin["zone_merges"]] if "zone_merges" in skin.files else []
    if "zone_merges" not in skin.files:
        print("ATTENTION : sortie P2-bis anterieure au traitement des micro-zones (relancer surface_skin.py)")
    zones = [z for z in sorted(zone_moment) if (zone == z).any() and (only is None or z in only)]
    print(f"{model} passe {k} : {len(zones)} zones, motif general {general}, angle {line_angle} deg, "
          f"recouvrement {overlap}, epaisseur {thickness}, cordons de bord {edge_rings}", flush=True)

    out, t_all = {}, time.time()
    for z in zones:
        faces = np.flatnonzero(zone == z)
        pattern = zone_pattern.get(z, general)
        angle = zone_angle.get(z, line_angle)
        entry = {"moment": int(zone_moment[z]), "pattern": pattern, "angle": angle, "faces": faces,
                 "area": float(skin["areas"][faces].sum()), "thickness": thickness, "edge_rings": edge_rings,
                 "merged_from": [int(m[1]) for m in merges if m[0] == z]}
        t0 = time.time()
        try:
            zm = sp.ZoneMesh(skin["vertices"], skin["faces"][faces])
            paths, pnormals, _d, warning = zone_paths(zm, normals[faces], pattern, angle, overlap, thickness,
                                                      order, width, edge_rings, edge_overlap,
                                                      zone_remedy.get(z, remedy))
            if warning:
                entry["warning"] = warning
            entry.update(paths=paths, normals=pnormals, extrusion=[1.0] * len(paths))
            if paths:
                allp = np.vstack(paths)
                sub = trimesh.Trimesh(skin["vertices"], skin["faces"][faces], process=False)
                depth = trimesh.proximity.closest_point(sub, allp)[1]
                entry["depth_ok"] = float((np.abs(depth - thickness / 2) <= DEPTH_TOLERANCE).mean())
                tree = cKDTree(np.vstack([sp.densify(p + n * thickness / 2) for p, n in zip(paths, pnormals)]))
                samples = trimesh.sample.sample_surface(sub, max(200, int(entry["area"] * 20)), seed=0)[0]
                entry["coverage"] = float((tree.query(samples)[0] <= width / 2 + 0.02).mean())
            else:
                entry["depth_ok"] = entry["coverage"] = 0.0
        except Exception as exc:                      # zone signalee, jamais ignoree en silence
            entry.update(paths=[], normals=[], extrusion=[], error=f"{type(exc).__name__}: {exc}")
        entry["seconds"] = time.time() - t0
        out[z] = entry
        msg = entry.get("error") or (f"{len(entry['paths']):4d} passes, profondeur ok {entry['depth_ok']:6.1%}, "
                                     f"couverture {entry['coverage']:6.1%}")
        print(f"  zone {z:3d} moment {entry['moment']:3d} {entry['area']:7.1f} mm2 {pattern:<20} {msg} "
              f"| {entry['seconds']:.1f} s", flush=True)
        check_memory(f"zone {z}")
    path = f"{prefix}_skin_paths_pass{k}.pkl"
    pickle.dump({"model": model, "pass": k, "general_pattern": general, "zones": out}, open(path, "wb"))
    ok = [e for e in out.values() if "error" not in e and e["paths"]]
    area = sum(e["area"] for e in out.values())
    print(f"\n{len(ok)}/{len(out)} zones tracees ; couverture ponderee "
          f"{sum(e['coverage'] * e['area'] for e in ok) / max(area, 1e-9):.1%} ; profondeur ok ponderee "
          f"{sum(e['depth_ok'] * e['area'] for e in ok) / max(sum(e['area'] for e in ok), 1e-9):.1%} ; "
          f"{time.time() - t_all:.0f} s au total")
    warned = {z: e["warning"] for z, e in out.items() if "warning" in e}
    if warned:
        print(f"\nAVERTISSEMENT Hilbert ({len(warned)} zone(s) : {sorted(warned)}) -- tracees telles quelles (A).")
        for z, w in sorted(warned.items()):
            print(f"  zone {z} : {w['message']}")
        print("Choix possibles (--hilbert-remedy=A|B|C ou --zone-hilbert-remedy=zone:B,...) :")
        for key, (title, text) in sp.HILBERT_OPTIONS.items():
            print(f"  {key}. {title} : {text}")
    errors = {z: e["error"] for z, e in out.items() if "error" in e}
    if errors:
        print(f"zones en echec ({len(errors)}) : " + "; ".join(f"{z}: {m}" for z, m in errors.items()))
    print(path)


if __name__ == "__main__":
    main()
