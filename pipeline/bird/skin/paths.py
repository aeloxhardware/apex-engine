# Porte de research/S4_Slicer/skin_paths.py (2026-10-07) : calcul copie TEL QUEL, main() devenu
# trace() (entrees/sorties en memoire). Voie Bird d'Apex, derive de S4 Slicer (Joshua Bird, GPL-3.0).
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
import os
import pickle
import sys
import time

import numpy as np
import trimesh
from scipy.spatial import cKDTree

from pipeline.bird.skin import surface_paths as sp
from pipeline.bird.skin.skin_carve import check_memory

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


def dedupe_shared_passes(out, width, radius_factor=0.5, min_cos=0.9):
    """Retire d'une zone les morceaux de passe deja deposes par une zone imprimee AVANT elle
    (moment plus petit, puis numero) : meme endroit (< radius_factor * w) et meme normale.
    Mesure (2026-10-06, Squirtle passe 8) : des zones voisines tracent chacune la passe de bord
    de leur frontiere commune (zones 7/10 : ~75 mm a 0,02 mm l'une de l'autre). Rend le nombre
    de points retires. Les morceaux restants plus courts qu'une largeur de buse sont abandonnes."""
    order = sorted((z for z, e in out.items() if e.get("paths")), key=lambda z: (out[z]["moment"], z))
    done_pts, done_nrm, removed = [], [], 0
    for z in order:
        e = out[z]
        if done_pts:
            tree = cKDTree(np.vstack(done_pts))
            allnrm = np.vstack(done_nrm)
            paths, normals, ext = [], [], []
            for p, n, x in zip(e["paths"], e["normals"], e["extrusion"]):
                p, n = np.asarray(p, float), np.asarray(n, float)
                dist, j = tree.query(p)
                dup = (dist < radius_factor * width) & ((n * allnrm[np.minimum(j, len(allnrm) - 1)]).sum(1) > min_cos)
                removed += int(dup.sum())
                if not dup.any():
                    paths.append(p); normals.append(n); ext.append(x)
                    continue
                edges = np.diff(np.concatenate([[1], dup.astype(int), [1]]))
                for a, b in zip(np.flatnonzero(edges == -1), np.flatnonzero(edges == 1)):
                    q = p[a:b]
                    if len(q) >= 2 and np.linalg.norm(np.diff(q, axis=0), axis=1).sum() >= width:
                        paths.append(q); normals.append(n[a:b]); ext.append(x)
            e["paths"], e["normals"], e["extrusion"] = paths, normals, ext
        for p, n in zip(e["paths"], e["normals"]):
            p, n = np.asarray(p, float), np.asarray(n, float)
            q = sp.densify(p)
            # normale du sommet le plus proche pour chaque point densifie
            nn = cKDTree(p).query(q)[1]
            done_pts.append(q); done_nrm.append(n[nn])
    return removed



def trace(skin, pattern="rectilinear", zone_pattern=None, line_angle=45.0, zone_angle=None, overlap=0.2,
          thickness=0.2, order="outside_in", only=None, edge_rings=1, edge_overlap=sp.EDGE_OVERLAP,
          remedy=None, zone_remedy=None, dedupe=True, log=False, cache=None):
    """Trajectoires de lissage de toutes les zones (P4), en memoire. `skin` :
    sortie de zones.compute. Rend {zone: entree} (ce que la recherche ecrivait
    dans _skin_paths_pass<k>.pkl, cle "zones").

    `cache` (dict, optionnel) : resultats par zone, cles par une empreinte de TOUT
    ce que le calcul d'une zone lit (indices et coordonnees de ses sommets,
    normales, aires, parametres). La boucle de rognage retrace a chaque tour des
    zones inchangees : elles sont reprises du cache, resultat identique (copie
    profonde, `dedupe_shared_passes` modifiant les entrees ensuite)."""
    general = pattern
    zone_pattern = zone_pattern or {}
    zone_angle = zone_angle or {}
    zone_remedy = zone_remedy or {}
    for name in [general, *zone_pattern.values()]:
        if name not in sp.PATTERNS:
            raise ValueError(f"motif inconnu : {name} (disponibles : {', '.join(sp.PATTERNS)})")
    for r in [remedy, *zone_remedy.values()]:
        if r not in (None, "A", "B", "C"):
            raise ValueError(f"remede Hilbert inconnu : {r} (A, B ou C)")
    normals = skin["normals"]
    width = float(skin["nozzle_width"])
    zone, zone_moment = skin["zone"], {z: int(m) for z, m in enumerate(skin["zone_moment"])}
    merges = [tuple(r) for r in skin["zone_merges"]] if "zone_merges" in skin else []
    zones = [z for z in sorted(zone_moment) if (zone == z).any() and (only is None or z in only)]
    out, t_all = {}, time.time()
    for z in zones:
        faces = np.flatnonzero(zone == z)
        pattern = zone_pattern.get(z, general)
        angle = zone_angle.get(z, line_angle)
        entry = {"moment": int(zone_moment[z]), "pattern": pattern, "angle": angle, "faces": faces,
                 "area": float(skin["areas"][faces].sum()), "thickness": thickness, "edge_rings": edge_rings,
                 "merged_from": [int(m[1]) for m in merges if m[0] == z]}
        t0 = time.time()
        key = None
        if cache is not None:
            import copy
            import hashlib
            h = hashlib.sha256()
            zf = skin["faces"][faces]
            for arr in (zf, skin["vertices"][zf], normals[faces], skin["areas"][faces]):
                h.update(np.ascontiguousarray(arr).tobytes())
            h.update(repr((pattern, angle, overlap, thickness, order, width, edge_rings, edge_overlap,
                           zone_remedy.get(z, remedy))).encode())
            key = h.hexdigest()
            if key in cache:
                cache["_hits"] = cache.get("_hits", 0) + 1
                entry.update(copy.deepcopy(cache[key]))
                entry["seconds"] = time.time() - t0
                out[z] = entry
                continue
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
        if key is not None:
            cache["_misses"] = cache.get("_misses", 0) + 1
            import copy
            cache[key] = copy.deepcopy({k: entry[k] for k in ("paths", "normals", "extrusion", "warning", "depth_ok",
                                                              "coverage", "error") if k in entry})
        msg = entry.get("error") or (f"{len(entry['paths']):4d} passes, profondeur ok {entry['depth_ok']:6.1%}, "
                                     f"couverture {entry['coverage']:6.1%}")
        if log: print(f"  zone {z:3d} moment {entry['moment']:3d} {entry['area']:7.1f} mm2 {pattern:<20} {msg} "
              f"| {entry['seconds']:.1f} s", flush=True)
        check_memory(f"zone {z}")
    if dedupe:
        removed = dedupe_shared_passes(out, width)
        if log: print(f"passes partagees entre zones retirees : {removed} points")
    return out
