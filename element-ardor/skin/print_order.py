# Porte de pipeline/bird/skin/print_order.py (branche bird_old, 2026-10-28) pour la voie chaleur : imports seulement ;
# `body` = couches du corps creuse (voir `body_from_layers`). Derive de S4 Slicer (Joshua Bird, GPL-3.0).
# Porte de research/S4_Slicer/skin_print.py (2026-10-07) : calcul copie TEL QUEL, main() devenu
# order_and_check() (entrees/sorties en memoire). Voie Bird d'Apex, derive de S4 Slicer (Joshua Bird, GPL-3.0).
"""
P5 -- ordre d'impression des couches de lissage non planaire et verification
de collision avec la VRAIE tete (plan approuve, phase P5).

Entrees (meme modele, meme passe) :
  - corps creuse : loop_<m>_deformed_slice_skin_pass<k>.pkl (P3), couches dans
    l'ordre d'impression, trajectoires dans le repere du STL ;
  - zones : loop_<m>_skin_pass<k>.npz (P2-bis), repere Bird (plateau z = 0,
    axe du plateau en x = y = 0) ;
  - trajectoires de lissage : loop_<m>_skin_paths_pass<k>.pkl (P4), repere Bird.

Ce que fait P5 :
  1. chaque zone est inseree juste apres la DERNIERE couche du corps (dans
     l'ordre d'impression) dont le rang est <= son moment : la matiere qui la
     porte est alors entierement deposee ;
  2. angles de tete par point = cinematique REELLE de la machine (profil par
     defaut d'Apex : C impose par la position, limite #59), cible = normale
     de surface ; ecart au-dela de --head-tolerance (45 deg) signale ;
  3. verification avec le modele de tete reel (aetheris-4.glb) et les
     primitives de l'etape 9 : points de tete echantillonnes, zone de contact
     de la buse exclue (3 mm), plateau (z < 0) et matiere DEJA imprimee a ce
     moment (capsules : corps + zones precedentes), marge de securite 1 mm ;
  4. une zone dont un point collisionne est RETIREE et signalee (fichier
     _skin_excluded_pass<k>.json) : relancer P3 avec --skin-exclude=<fichier>
     annule son creusement (une seule fois, plan P5), puis relancer P5.

Sortie : loop_<m>_print_order_pass<k>.pkl -- liste ordonnee de passes
{kind: "body"|"skin", layer|zone, points (N,3) repere Bird, angles (N, axes)
pour le lissage, normals}, et le rapport.

Usage (venv d'Apex) :
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe -u skin_print.py <modele STL> <passe> [variante]
      [--head-tolerance=45] [--sampling=0 (0 = tous les points)] [--margin=1.0] [--feed=25 (mm/s)]
"""
import json
import os
import pickle
import sys
import time

import numpy as np

from pipeline.skin._memory import check_memory  # noqa: E402
from pipeline.step8_machine_toolpath import Step8Params, _RawPoint, _resolve_all_axis_angles_coupled  # noqa: E402
from pipeline.step8_machine_toolpath import _resolve_axis_angles  # noqa: E402
from pipeline.step9_collision_sim import (  # noqa: E402
    CollisionSimParams, _filter_nozzle_contact_zone, _positioned_toolhead_vertices, _surface_sample_points,
)

HEAD_TOLERANCE_DEG = 45.0


def head_angles(points, normals, machine, step8):
    """Angles de tete (N, axes) et ecart a la normale (deg) par point."""
    axes = machine.rotation_axes
    if machine.position_coupled_axes:
        raw = [_RawPoint(position=p, desired_direction=n, layer_id=0, segment_index=0, color_id=None,
                         feedrate=1.0, extrusion=1.0, is_travel=False) for p, n in zip(points, normals)]
        a, err = _resolve_all_axis_angles_coupled(raw, axes, np.array([0.0, 0.0, 1.0]), machine, step8)
        return np.asarray(a), np.degrees(err)
    out, err, prev = [], [], None
    for n in normals:
        a, e = _resolve_axis_angles(n, axes, np.array([0.0, 0.0, 1.0]), step8, prev)
        out.append(a)
        err.append(np.degrees(e))
        prev = a
    return np.asarray(out), np.asarray(err)


def feed_plan(points, angles, machine, feed_mm_s):
    """Temps de parcours d'une passe quand le firmware COORDONNE les axes
    (reponse de l'utilisateur du 2026-09-28) : les angles restent EXACTS et
    c'est le mouvement entier qui ralentit la ou un axe d'orientation (non
    couple a la position) atteint sa vitesse ou son acceleration maximale.
    Mouvement 1D de l'angle par segment : vitesse angulaire plafonnee par
    passes avant/arriere (v^2 <= v0^2 + 2 a |dtheta|, remise a zero a chaque
    changement de sens), comme un planificateur de firmware.
    Rend (temps nominal, temps reel, facteur de vitesse par segment)."""
    seg = np.linalg.norm(np.diff(points, axis=0), axis=1)
    t_nom = seg / max(feed_mm_s, 1e-9)
    t_real = t_nom.copy()
    vmax_deg = np.degrees(machine.max_angular_velocity)
    amax_deg = np.degrees(machine.max_angular_acceleration)
    coupled = set(machine.position_coupled_axes or [])
    if len(seg) and vmax_deg > 0 and amax_deg > 0:
        for ax in range(angles.shape[1]):
            if ax in coupled:
                continue                           # C : limites inconnues, jamais lisse
            d = np.diff(np.unwrap(np.radians(angles[:, ax])))
            d = np.degrees(d)
            mag = np.abs(d)
            cap = np.full(len(d), vmax_deg)
            prev = 0.0
            for i in range(len(d)):                # avant
                if i and np.sign(d[i]) != np.sign(d[i - 1]):
                    prev = 0.0
                cap[i] = min(cap[i], np.sqrt(prev ** 2 + 2 * amax_deg * mag[i]))
                prev = cap[i]
            nxt = 0.0
            for i in range(len(d) - 1, -1, -1):    # arriere
                if i < len(d) - 1 and np.sign(d[i]) != np.sign(d[i + 1]):
                    nxt = 0.0
                cap[i] = min(cap[i], np.sqrt(nxt ** 2 + 2 * amax_deg * mag[i]))
                nxt = cap[i]
            t_real = np.maximum(t_real, mag / np.maximum(cap, 1e-9))
    factor = np.where(t_real > 0, t_nom / np.maximum(t_real, 1e-12), 1.0)
    return float(t_nom.sum()), float(t_real.sum()), factor


def sample_indices(points, step):
    """Indices espaces d'environ `step` mm le long de la polyligne (premier et
    dernier compris), comme l'etape 9."""
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    marks = np.arange(0, s[-1] + 1e-9, step)
    idx = np.unique(np.r_[np.searchsorted(s, marks), len(points) - 1])
    return idx[idx < len(points)]



def _skin_check_native():
    """apex_core.SkinCheckIndex si `APEX_USE_CPP=1`."""
    from core import _cpp
    module = _cpp.module() if _cpp.enabled() else None
    return getattr(module, "SkinCheckIndex", None) if module is not None else None


def _lines(polylines):
    lines = [np.asarray(p, dtype=float).reshape(-1, 3) for p in polylines]
    starts = np.r_[0, np.cumsum([len(p) for p in lines])].astype(np.int64)
    pts = np.vstack(lines) if lines else np.zeros((0, 3))
    return np.ascontiguousarray(pts), starts


def _check_native(native_cls, sequence, body, paths, skin_angles, offset, head, machine, sim, plate_r, reach,
                  densify, tolerance, over_tol, report, collision_points):
    """Partie 3 de order_and_check en C++ : remplit `report` et `collision_points`
    exactement comme la boucle Python, rend le nombre de points testes."""
    index = native_cls(reach + 1.0)
    axes = np.array([np.asarray(spec.axis, dtype=float) / np.linalg.norm(np.asarray(spec.axis, dtype=float))
                     for spec in machine.rotation_axes]).reshape(-1, 3)
    head = np.ascontiguousarray(head, dtype=float)
    tested = 0
    for kind, ref in sequence:
        if kind == "body":
            pts, starts = _lines([t - offset for t in body[ref]["toolpath_3d"]])
            index.add_polylines(pts, starts, densify)
            continue
        z = ref
        e = paths[z]
        angles_z, err_z = skin_angles[z]
        hits = {"plateau": 0, "piece": 0, "tolerance": 0}
        for p, err in zip(e["paths"], err_z):
            far = np.flatnonzero(err > tolerance)
            hits["tolerance"] += len(far)
            collision_points.extend(p[far].tolist())
        pts, starts = _lines(e["paths"])
        rad = np.radians(np.vstack([np.asarray(a, dtype=float).reshape(len(p), -1)
                                    for p, a in zip(e["paths"], angles_z)])) if len(pts) else np.zeros((0, len(axes)))
        idx = [np.arange(len(p)) if sim.path_sampling_step <= 0 else sample_indices(p, sim.path_sampling_step)
               for p in e["paths"]]
        tested_idx = np.concatenate([i + s for i, s in zip(idx, starts[:-1])]).astype(np.int64) if idx else             np.zeros(0, np.int64)
        # Version parallele (2026-10-13) quand le module la porte : meme resultat au bit
        # (tests/test_bird_skin.py), points testes repartis sur les coeurs.
        check = getattr(index, "check_zone_parallel", None) or index.check_zone
        plate, dist = check(head, axes, pts, starts, np.ascontiguousarray(np.cos(rad)),
                            np.ascontiguousarray(np.sin(rad)), tested_idx, plate_r, sim.safety_margin,
                            reach, densify)
        tested += len(tested_idx)
        hits["plateau"] += int(plate.sum())
        hits["piece"] += int((dist < reach).sum())
        hit = (dist < reach) | plate
        collision_points.extend(pts[tested_idx[hit]].tolist())
        closest = float(dist.min()) if len(dist) else np.inf
        report[z] = {"moment": e["moment"], "area": e["area"], "hits": hits,
                     "worst_penetration": float(max(0.0, reach - closest)),
                     "over_tolerance": over_tol[z],
                     "max_head_deviation": float(max(x.max() for x in err_z)) if err_z else 0.0}
        if not any(hits.values()):
            index.add_polylines(pts, starts, densify)
        check_memory(f"zone {z}")
    return tested


def order_and_check(body, skin, paths, machine, sim, offset, tolerance=HEAD_TOLERANCE_DEG, plate_r=150.0,
                    exclude=(), feed=25.0, log=False):
    """P5 en memoire : insere chaque zone de lissage apres la derniere couche du
    corps qui la porte, calcule les angles de tete, verifie les collisions avec
    la VRAIE tete, retire les zones en contact. `body` : couches du corps
    (deformed_slicing, toolpath_3d dans le repere de sortie, `offset` = repere
    Bird -> sortie) ; `skin` : zones.compute ; `paths` : paths.trace.
    Rend {passes, report, excluded, collision_points, monotone_ranks}."""
    t_all = time.time()
    isovalues = skin["isovalues"]
    ranks = np.array([int(np.clip(np.searchsorted(isovalues, L["isovalue"] - 1e-9), 0, len(isovalues) - 1))
                      for L in body])
    monotone = bool(np.all(np.diff(ranks) >= 0))

    # --- 1. ordre : chaque zone apres la derniere couche de rang <= moment ------
    skip = set(int(z) for z in (exclude or ()))     # zones deja retirees (creusement annule par P3)
    zones = sorted((z for z, e in paths.items() if e["paths"] and z not in skip), key=lambda z: paths[z]["moment"])
    after = {}
    for z in zones:
        carrying = np.flatnonzero(ranks <= paths[z]["moment"])
        after.setdefault(int(carrying.max()) if len(carrying) else -1, []).append(z)
    sequence = [("skin", z) for z in after.get(-1, [])]
    for i in range(len(body)):
        sequence.append(("body", i))
        sequence.extend(("skin", z) for z in after.get(i, []))
    if log: print(f"P5 : {len(body)} couches de corps (rangs {'croissants' if monotone else 'NON monotones'}), "
          f"{len(zones)} zones de lissage inserees", flush=True)

    # --- 2. angles de tete ------------------------------------------------------
    step8 = Step8Params()
    t0 = time.time()
    skin_angles, over_tol = {}, {}
    for z in zones:
        e = paths[z]
        per_path_a, per_path_e = [], []
        e["feed_factor"], e["time_nominal"], e["time_real"] = [], 0.0, 0.0
        for p, n in zip(e["paths"], e["normals"]):
            a, err = head_angles(p, n, machine, step8)          # angles EXACTS (axes coordonnes)
            tn, tr, factor = feed_plan(p, a, machine, feed)
            e["feed_factor"].append(factor)
            e["time_nominal"] += tn
            e["time_real"] += tr
            per_path_a.append(a)
            per_path_e.append(err)
        skin_angles[z] = (per_path_a, per_path_e)
        over_tol[z] = float(np.concatenate(per_path_e).__gt__(tolerance).mean())
    tn = sum(paths[z]["time_nominal"] for z in zones)
    tr = sum(paths[z]["time_real"] for z in zones)
    if log: print(f"vitesse de depot {feed:.0f} mm/s ; tete B limitee a {np.degrees(machine.max_angular_velocity):.0f} deg/s et "
          f"{np.degrees(machine.max_angular_acceleration):.0f} deg/s2 : temps des passes de lissage "
          f"{tn / 60:.1f} min -> {tr / 60:.1f} min (x{tr / max(tn, 1e-9):.2f})", flush=True)
    if log: print(f"angles de tete ({'C impose par la position' if machine.position_coupled_axes else 'C libre'}) : "
          f"{time.time() - t0:.1f} s ; zones avec des points au-dela de {tolerance:.0f} deg : "
          f"{sum(1 for v in over_tol.values() if v > 0)}", flush=True)

    # --- 3. collisions avec la vraie tete -----------------------------------------
    # Meme critere que l'etape 9 (point de tete a moins de marge + rayon de
    # cordon d'une capsule imprimee), mais vectorise : la matiere imprimee est
    # un nuage de points (trajectoires densifiees a DENSIFY mm, erreur <= DENSIFY/2)
    # dans un arbre k-d. La boucle Python par capsule de l'etape 9 mettait plus
    # de 10 min sur le Benchy (1,16 M de points de corps, tete de 77 mm).
    from scipy.spatial import cKDTree
    DENSIFY = 0.25
    BLOCK = 256
    reach = sim.safety_margin + machine.nozzle_diameter / 2

    def dense(polylines):
        out = []
        for p in polylines:
            if len(p) < 2:
                out.append(p)
                continue
            seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
            n = np.maximum(1, np.ceil(seg / DENSIFY).astype(int))
            t = np.concatenate([np.arange(m) / m for m in n])
            i = np.repeat(np.arange(len(seg)), n)
            out.append(np.vstack([p[i] + t[:, None] * (p[i + 1] - p[i]), p[-1:]]))
        return out

    head = _surface_sample_points(machine.nozzle_model, sim.toolhead_sample_spacing)
    head = _filter_nozzle_contact_zone(head, sim.nozzle_contact_clearance)
    body_pts, body_tree, body_upto = [], None, -1
    skin_pts = []
    report, t0, tested = {}, time.time(), 0
    collision_points = []
    native = _skin_check_native()
    if native is not None:
        # Meme controle en C++ (apex/skin_check.hpp) : memes obstacles, meme distance
        # bornee, meme densification ; la boucle ci-dessous n'est alors pas executee.
        tested = _check_native(native, sequence, body, paths, skin_angles, offset, head, machine, sim, plate_r,
                               reach, DENSIFY, tolerance, over_tol, report, collision_points)
    for kind, ref in (sequence if native is None else ()):
        if kind == "body":
            body_pts.extend(dense([t - offset for t in body[ref]["toolpath_3d"]]))
            continue
        z = ref
        e = paths[z]
        if body_upto != len(body_pts):                 # le corps a avance : arbre reconstruit
            body_tree = cKDTree(np.vstack(body_pts)) if body_pts else None
            body_upto = len(body_pts)
        skin_tree = cKDTree(np.vstack(skin_pts)) if skin_pts else None
        angles_z, err_z = skin_angles[z]
        hits = {"plateau": 0, "piece": 0, "tolerance": 0}
        for p, err in zip(e["paths"], err_z):          # au-dela de la tolerance tete/normale : rogne aussi
            far = np.flatnonzero(err > tolerance)
            hits["tolerance"] += len(far)
            collision_points.extend(p[far].tolist())
        closest = np.inf
        # Correctifs du 2026-09-30 (decision utilisateur) : TOUS les points sont
        # testes (un contact de 0,02 mm tombait entre deux echantillons de 2 mm
        # sur le dino), et la tete est aussi testee contre les lignes deja
        # deposees de la zone en cours et le debut de la ligne en cours, comme
        # full_step89.py --fast.
        zone_done = []
        for p, a in zip(e["paths"], angles_z):
            zone_tree = cKDTree(np.vstack(zone_done)) if zone_done else None
            idx = (np.arange(len(p)) if sim.path_sampling_step <= 0
                   else sample_indices(p, sim.path_sampling_step))
            block_id, block_tree = -1, None
            for i in idx:
                tested += 1
                verts = _positioned_toolhead_vertices(head, p[i], list(a[i]), machine.rotation_axes)
                # debut de la ligne : arbre des blocs complets (reconstruit tous
                # les BLOCK points, sinon cout quadratique -- lignes de 20 000
                # points sur le dino) + petit arbre du bloc en cours
                if i >= 1 and i // BLOCK != block_id:
                    block_id = i // BLOCK
                    start = block_id * BLOCK
                    block_tree = cKDTree(dense([p[:start + 1]])[0]) if start >= 1 else None
                start = block_id * BLOCK if i >= 1 else 0
                partial_tree = cKDTree(dense([p[start:i]])[0]) if i - start >= 1 else None
                # plateau = disque z = 0 de rayon plate_r, meme marge
                on_plate = ((verts[:, 2] < sim.safety_margin)
                            & (np.linalg.norm(verts[:, :2], axis=1) < plate_r + sim.safety_margin)).any()
                if on_plate:
                    hits["plateau"] += 1
                d = np.inf
                for tree in (body_tree, skin_tree, zone_tree, block_tree if i >= 1 else None, partial_tree):
                    if tree is not None:
                        d = min(d, float(tree.query(verts, distance_upper_bound=reach + 1.0)[0].min()))
                closest = min(closest, d)
                if d < reach:
                    hits["piece"] += 1
                if d < reach or on_plate:
                    collision_points.append(p[i].tolist())
            zone_done.extend(dense([p]))
        report[z] = {"moment": e["moment"], "area": e["area"], "hits": hits,
                     "worst_penetration": float(max(0.0, reach - closest)),
                     "over_tolerance": over_tol[z],
                     "max_head_deviation": float(max(x.max() for x in err_z)) if err_z else 0.0}
        if not any(hits.values()):
            skin_pts.extend(dense(e["paths"]))         # zone deposee : obstacle pour la suite
        check_memory(f"zone {z}")
    if log: print(f"collisions : {tested} points de lissage testes en {time.time() - t0:.1f} s", flush=True)

    excluded = sorted(z for z, r in report.items() if any(r["hits"].values()))
    area_all = sum(r["area"] for r in report.values())
    area_ex = sum(report[z]["area"] for z in excluded)
    if log: print(f"\nzones en collision (RETIREES) : {len(excluded)}/{len(report)} ({area_ex:.0f} mm2 sur {area_all:.0f})")
    for z in excluded:
        r = report[z]
        if log: print(f"  zone {z:3d} moment {r['moment']:3d} {r['area']:7.1f} mm2 : plateau {r['hits']['plateau']}, "
              f"piece {r['hits']['piece']}, au-dela de {tolerance:.0f} deg {r['hits']['tolerance']} (penetration max {r['worst_penetration']:.2f} mm), "
              f"ecart tete max {r['max_head_deviation']:.0f} deg")
    kept = [z for z in report if z not in excluded]
    dev = [report[z]["max_head_deviation"] for z in kept]
    if dev:
        if log: print(f"zones gardees : {len(kept)}, ecart tete / normale max par zone : mediane {np.median(dev):.0f} deg, "
              f"max {max(dev):.0f} deg")

    excluded_set = set(excluded)
    passes = []
    for kind, ref in sequence:
        if kind == "body":
            passes.append({"kind": "body", "layer": ref,
                           "points": [t - offset for t in body[ref]["toolpath_3d"]]})
        elif ref not in excluded_set:
            passes.append({"kind": "skin", "zone": ref, "points": paths[ref]["paths"],
                           "normals": paths[ref]["normals"], "angles": skin_angles[ref][0],
                           "feed_factor": paths[ref]["feed_factor"]})
    if log: print(f"P5 : {time.time() - t_all:.0f} s au total", flush=True)
    return {"passes": passes, "report": report, "excluded": excluded, "monotone_ranks": monotone,
            "head_tolerance": tolerance, "collision_points": collision_points}
