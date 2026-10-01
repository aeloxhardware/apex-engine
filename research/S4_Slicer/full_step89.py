"""
Toute la trajectoire (corps Bird + passes de lissage, ordre de P5) dans les
etapes 8 et 9 d'Apex : cinematique REELLE (C impose par la position) et
VRAIE tete (aetheris-4). Jamais fait avant le 2026-09-29 : P5 ne verifiait
que le lissage.

- corps : direction visee = normale de la couche Bird dans la piece reelle,
  soit le gradient de phi (hauteur deformee) sur le tetraedre d'origine du point ;
- lissage : angles deja calcules par P5 ;
- premier point de chaque polyligne = trajet a vide (extrusion 0) ;
- etape 9 : boucle de step9.run reproduite avec ses primitives, pour pouvoir
  n'evaluer qu'un echantillon (--eval-every) ou borner la duree
  (--max-seconds) ; l'index de matiere, lui, croit TOUJOURS en entier.

Usage (venv d'Apex) :
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe -u full_step89.py <modele> <passe> <variante>
      [--eval-every=1] [--max-seconds=0] [--limit-points=0]
"""
import safe_io  # noqa: F401 -- ecritures robustes aux verrous transitoires
import os
import pickle
import sys
import time

import numpy as np
from scipy.ndimage import median_filter
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import main as apex_main  # noqa: E402
from bird_to_apex import barycentric, tet_gradients  # noqa: E402
from core.mesh_types import ToolpathMachine, ToolpathPoint  # noqa: E402
from pipeline import step9_collision_sim as s9  # noqa: E402
from pipeline.step8_machine_toolpath import Step8Params, _RawPoint, _resolve_all_axis_angles_coupled  # noqa: E402


def locate(points, p0, cells, centroid_tree, k=12):
    """Tetraedre d'origine contenant chaque point (sinon le plus proche des k candidats)."""
    _d, cand = centroid_tree.query(points, k=k)
    best = cand[:, 0].copy()
    score = np.full(len(points), -np.inf)
    for j in range(k):
        with np.errstate(all="ignore"):
            b = barycentric(points, p0[cells[cand[:, j]]])
        s = np.nan_to_num(b.min(axis=1), nan=-np.inf)
        better = s > score
        best[better], score[better] = cand[better, j], s[better]
    return best


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    model, k, variant = args[0], int(args[1]), args[2] if len(args) > 2 else ""
    every = int(flags.get("--eval-every", 1))
    max_s = float(flags.get("--max-seconds", 0))
    limit = int(flags.get("--limit-points", 0))
    plate_r = float(flags.get("--plate-radius", 150.0))   # plateau de 300 mm (utilisateur, 2026-09-29)
    smooth_b = int(flags.get("--smooth-b", 0))            # filtre median de B (points), 0 = aucun
    prefix = os.path.join(HERE, "output_models", "loop_" + model.replace(" ", "_") + variant)
    order = pickle.load(open(f"{prefix}_print_order_pass{k}.pkl", "rb"))["passes"]
    cells, p0, pk = (np.load(f"{prefix}_{s}.npy") for s in ("cells", "pass0", f"pass{k}"))
    machine = apex_main._default_machine_profile()
    step8 = Step8Params()

    # --- etape 8 : angles --------------------------------------------------------
    t0 = time.time()
    grad = tet_gradients(p0, cells, pk[:, 2])
    grad /= np.maximum(np.linalg.norm(grad, axis=1, keepdims=True), 1e-12)
    tree = cKDTree(p0[cells].mean(axis=1))
    polys = []          # (points (N,3), angles (N,axes), kind, passe)
    poly_normals = []   # normale visee par point (couche du corps, surface du lissage)
    body_err = []
    for pi, ps in enumerate(order):
        if ps["kind"] == "skin":
            polys.extend((p, np.asarray(a), "skin", pi) for p, a in zip(ps["points"], ps["angles"]))
            poly_normals.extend(np.asarray(nn) for nn in ps["normals"])
            continue
        for p in ps["points"]:
            n = grad[locate(p, p0, cells, tree)]
            raw = [_RawPoint(position=q, desired_direction=d, layer_id=ps["layer"], segment_index=0, color_id=None,
                             feedrate=1.0, extrusion=1.0, is_travel=False) for q, d in zip(p, n)]
            a, err = _resolve_all_axis_angles_coupled(raw, machine.rotation_axes, np.array([0.0, 0.0, 1.0]),
                                                      machine, step8)
            if smooth_b > 1 and len(p) > 2:
                # Filtre median de B le long de la polyligne (2026-09-30) : la
                # normale est constante par tetraedre ; pres de l'axe, B oscille
                # 0 <-> ~14 deg d'un point a l'autre (5,9 % des segments du corps
                # de Squirtle), ce qui ralentit le firmware (x3,7). Applique AVANT
                # la borne en hauteur et la reorientation : le controle de
                # collision valide les angles lisses. Ecart a la normale recalcule.
                a = np.asarray(a, dtype=float)
                a[:, 1] = median_filter(a[:, 1], size=smooth_b, mode="nearest")
                sb, cb = np.sin(np.radians(a[:, 1])), np.cos(np.radians(a[:, 1]))
                sc, cc = np.sin(np.radians(a[:, 0])), np.cos(np.radians(a[:, 0]))
                head = np.column_stack([sb * sc, -sb * cc, cb])
                err = np.arccos(np.clip(np.einsum("ij,ij->i", head, n), -1.0, 1.0))
            body_err.append(np.degrees(err))
            polys.append((p, np.asarray(a), "body", pi))
            poly_normals.append(n)
        if limit and sum(len(q[0]) for q in polys) >= limit:
            break
    body_err = np.concatenate(body_err) if body_err else np.zeros(0)
    pts, pass_of, poly_of = [], [], []
    for j, (p, a, kind, pi) in enumerate(polys):
        for i, (q, ang) in enumerate(zip(p, a)):
            pts.append(ToolpathPoint(position=np.asarray(q), axis_angles=list(ang), extrusion=0.0 if i == 0 else 1.0,
                                     feedrate=1.0, color_id=1 if kind == "skin" else 0))
        pass_of.extend([pi] * len(p))
        poly_of.extend([j] * len(p))
    pass_of, poly_of = np.asarray(pass_of), np.asarray(poly_of)
    print(f"etape 8 : {len(pts)} points ({len(polys)} polylignes) en {time.time() - t0:.0f} s ; corps : ecart tete / "
          f"normale de couche median {np.median(body_err):.1f} deg, 95e {np.percentile(body_err, 95):.1f}, "
          f"max {body_err.max():.1f} ; > 5 deg {np.mean(body_err > 5):.1%}, > 45 deg {np.mean(body_err > 45):.1%}",
          flush=True)

    # --- etape 9 -----------------------------------------------------------------
    sim = s9.CollisionSimParams(path_sampling_step=2.0, safety_margin=1.0)
    head = s9._filter_nozzle_contact_zone(s9._surface_sample_points(machine.nozzle_model, sim.toolhead_sample_spacing),
                                          sim.nozzle_contact_clearance)
    lo_b, hi_b = machine.rotation_axes[1].angle_range_deg
    grid_b = np.arange(lo_b, hi_b + 1e-9, 0.5)
    zlow = np.array([s9._positioned_toolhead_vertices(head, np.zeros(3), [0.0, float(bb)], machine.rotation_axes)
                     [:, 2].min() for bb in grid_b])

    def clamp_all():
        """Borne de B selon la hauteur de buse (redondance avec la
        reorientation, decision de l'utilisateur du 2026-09-29) : point le plus
        bas de la VRAIE tete pour chaque B (C tourne autour de Z : sans effet
        sur z) ; chaque point est ramene au B autorise le plus proche."""
        clamped, devs = 0, []
        for q in pts:
            ok = q.position[2] + zlow >= sim.safety_margin
            if not ok.any() or ok[np.abs(grid_b - q.axis_angles[1]).argmin()]:
                continue
            best = grid_b[ok][np.abs(grid_b[ok] - q.axis_angles[1]).argmin()]
            devs.append(abs(best - q.axis_angles[1]))
            q.axis_angles = [q.axis_angles[0], float(best)] + list(q.axis_angles[2:])
            clamped += 1
        return clamped, devs

    if "--height-clamp" in sys.argv:
        clamped, devs = clamp_all()
        print(f"borne de B selon la hauteur : {clamped} points ramenes ({clamped / len(pts):.2%}) ; |dB| median "
              f"{np.median(devs) if devs else 0:.1f} deg, max {max(devs) if devs else 0:.1f}", flush=True)
    samples = s9._sample_point_indices(pts, sim.path_sampling_step)[::every]
    events, last, t0 = [], 0, time.time()
    if "--fast" in sys.argv:
        # Meme critere que P5 (skin_print.py) : matiere imprimee = points des
        # polylignes deja deposees, densifies a 0,25 mm, dans un arbre k-d ;
        # contact si un point de tete est a moins de marge + rayon de cordon ;
        # plateau = plan z < marge. Arbre de tout ce qui precede la passe
        # courante (reconstruit a chaque passe) + arbre du debut de la passe.
        reach = sim.safety_margin + machine.nozzle_diameter / 2

        def dense(p, step=0.25):
            if len(p) < 2:
                return p
            seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
            m = np.maximum(1, np.ceil(seg / step).astype(int))
            t = np.concatenate([np.arange(c) / c for c in m])
            i = np.repeat(np.arange(len(seg)), m)
            return np.vstack([p[i] + t[:, None] * (p[i + 1] - p[i]), p[-1:]])

        starts = np.r_[0, np.cumsum([len(q[0]) for q in polys])]
        # paliers de tolerance tete / normale visee : un ecart de depot vaut
        # mieux qu'une collision (dernier palier = sans limite)
        tiers = [float(x) for x in flags.get("--reorient-tiers", "45,90,181").split(",")]
        b_grid = np.arange(lo_b, hi_b + 1e-9, 2.5)

        def head_dir(angles):
            return s9._positioned_toolhead_vertices(np.array([[0.0, 0.0, 1.0]]), np.zeros(3), angles,
                                                    machine.rotation_axes)[0]
        dpolys = [dense(np.asarray(q[0])) for q in polys]

        def apply_found(found_list, full=1.0, taper=3.0):
            """Ecrit chaque B trouve dans la trajectoire : plein a +-full mm
            (abscisse curviligne, meme polyligne), degressif jusqu'a +-taper
            mm. Recouvrement : la correction ponderee la plus ample gagne."""
            best = {}
            for idx, _z, b0, f, _k in found_list:
                if f is None:
                    continue
                j = int(poly_of[idx])
                p = np.asarray(polys[j][0])
                s_arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
                dist = np.abs(s_arc - s_arc[idx - starts[j]])
                w = np.clip((taper - dist) / (taper - full), 0, 1)
                for li in np.flatnonzero(w > 0):
                    g = int(starts[j] + li)
                    delta = w[li] * (f[0] - b0)
                    if g not in best or abs(delta) > abs(best[g]):
                        best[g] = delta
            for g, delta in best.items():
                a = pts[g].axis_angles
                pts[g].axis_angles = [a[0], float(np.clip(a[1] + delta, lo_b, hi_b))] + list(a[2:])
            return len(best)

        rounds = int(flags.get("--apply-rounds", 4)) if "--apply" in sys.argv else 1
        reorient = "--reorient" in sys.argv or "--apply" in sys.argv
        for rnd in range(rounds):
          events, reoriented = [], []
          before_tree, before_pass, prefix_tree, prefix_key = None, -1, None, None
          for n, idx in enumerate(samples):
            pi, j = int(pass_of[idx]), int(poly_of[idx])
            if pi != before_pass:
                done = [dpolys[m] for m in range(len(polys)) if polys[m][3] < pi]
                before_tree = cKDTree(np.vstack(done)) if done else None
                before_pass = pi
            if prefix_key != j:
                same = [dpolys[m] for m in range(j) if polys[m][3] == pi]
                prefix_tree = cKDTree(np.vstack(same)) if same else None
                prefix_key = j
            local = idx - starts[j]
            partial = np.asarray(polys[j][0][:local]) if local >= 1 else None
            partial_tree = cKDTree(dense(partial)) if partial is not None and len(partial) else None

            def check(angles):
                """(penetration plateau ou None, penetration piece ou None)."""
                verts = s9._positioned_toolhead_vertices(head, pts[idx].position, angles, machine.rotation_axes)
                over_disc = np.linalg.norm(verts[:, :2], axis=1) < plate_r + sim.safety_margin
                low = verts[over_disc, 2]
                plate = float(sim.safety_margin - low.min()) if len(low) and low.min() < sim.safety_margin else None
                d = np.inf
                for tree in (before_tree, prefix_tree, partial_tree):
                    if tree is not None:
                        d = min(d, float(tree.query(verts, distance_upper_bound=reach + 1.0)[0].min()))
                return plate, (float(reach - d) if d < reach else None)

            plate, part = check(pts[idx].axis_angles)
            if plate is not None:
                events.append(s9.CollisionEvent(idx, "platform", plate))
            if part is not None:
                events.append(s9.CollisionEvent(idx, "printed_part", part))
            if reorient and (plate is not None or part is not None):
                # Niveau 1 de step10 avec le controle rapide : autre B (C reste
                # impose par la position), du plus proche au plus eloigne, dans
                # la tolerance tete / normale visee.
                a0 = list(pts[idx].axis_angles)
                target = poly_normals[j][local]
                target = target / max(np.linalg.norm(target), 1e-12)
                cands = sorted(b_grid, key=lambda b: abs(b - a0[1]))
                devs = {float(b): np.degrees(np.arccos(np.clip(
                    head_dir([a0[0], float(b)] + a0[2:]) @ target, -1, 1))) for b in cands}
                found = None
                for tier, tol in enumerate(tiers):
                    for b in cands:
                        dv = devs[float(b)]
                        if dv > tol or (tier and dv <= tiers[tier - 1]):
                            continue
                        if pts[idx].position[2] + zlow[np.abs(grid_b - b).argmin()] < sim.safety_margin:
                            continue
                        if check([a0[0], float(b)] + a0[2:]) == (None, None):
                            found = (float(b), dv, tier)
                            break
                    if found:
                        break
                reoriented.append((idx, pts[idx].position[2], a0[1], found,
                                   "skin" if pts[idx].color_id else "body"))
            if n and n % 20000 == 0:
                print(f"  {n}/{len(samples)} points testes, {time.time() - t0:.0f} s, {len(events)} evenements",
                      flush=True)
          ok = [x for x in reoriented if x[3] is not None]
          by_tier = [sum(1 for x in ok if x[3][2] == t) for t in range(len(tiers))]
          print(f"tour {rnd + 1} : {len({e.point_index for e in events})} points en contact ; resolus par un autre B : "
                f"{len(ok)}/{len(reoriented)} (paliers {'/'.join(f'{t:g}' for t in tiers)} deg : {by_tier})",
                flush=True)
          if "--apply" not in sys.argv or not ok or rnd == rounds - 1:
              break                  # le dernier tour ne fait que VERIFIER
          changed = apply_found(reoriented)
          if "--height-clamp" in sys.argv:
              clamp_all()
          print(f"  applique a {changed} points (lissage +-1 / 3 mm)", flush=True)
        bad = [x for x in reoriented if x[3] is None]
        if bad:
            z = np.array([x[1] for x in bad])
            print(f"  NON resolus : {len(bad)} ; z median {np.median(z):.1f} mm ; B median "
                  f"{np.median([x[2] for x in bad]):.0f} deg")
        pickle.dump(reoriented, open(f"{prefix}_reorient_pass{k}.pkl", "wb"))
        samples_iter = []
    else:
        samples_iter = samples
    index = s9._PrintedGeometryIndex(machine.nozzle_diameter, sim.grid_cell_size_multiplier)
    for n, idx in enumerate(samples_iter):
        last = s9._grow_printed_index_to(index, pts, idx, last)
        events.extend(s9._evaluate_point(idx, head, pts[idx], machine, index, sim.safety_margin))
        if n and n % 2000 == 0:
            el = time.time() - t0
            print(f"  {n}/{len(samples)} points testes, {el:.0f} s ({n / el:.0f} pts/s), {len(events)} collisions",
                  flush=True)
        if max_s and time.time() - t0 > max_s:
            print(f"  arret a {max_s:.0f} s : {n + 1}/{len(samples)} points testes (jusqu'au point {idx}/{len(pts)})")
            break
    el = time.time() - t0
    kinds = {}
    for e in events:
        key = (e.obstacle_type, "skin" if pts[e.point_index].color_id else "body")
        kinds[key] = kinds.get(key, 0) + 1
    hit_pts = sorted({e.point_index for e in events})
    print(f"etape 9 : {len(samples)} echantillons (1 sur {every}) en {el:.0f} s ; {len(events)} evenements sur "
          f"{len(hit_pts)} points ; par type : {kinds}")
    if events:
        pen = np.array([e.penetration_depth for e in events])
        print(f"  penetration mediane {np.median(pen):.2f} mm, max {pen.max():.2f} mm")
    pickle.dump({"events": [(e.point_index, e.obstacle_type, e.penetration_depth) for e in events],
                 "positions": np.array([pts[i].position for i in hit_pts]) if hit_pts else np.zeros((0, 3)),
                 "body_err": body_err}, open(f"{prefix}_full_step89_pass{k}.pkl", "wb"))
    if "--apply" in sys.argv:
        # trajectoire machine CORRIGEE (borne + reorientation), entree de
        # l'export G-code : positions (repere plateau), angles, extrusion,
        # corps / lissage
        pickle.dump({"position": np.array([q.position for q in pts]),
                     "axis_angles": np.array([q.axis_angles for q in pts]),
                     "extrusion": np.array([q.extrusion for q in pts]),
                     "is_skin": np.array([bool(q.color_id) for q in pts])},
                    open(f"{prefix}_machine_toolpath_pass{k}.pkl", "wb"))


if __name__ == "__main__":
    main()
