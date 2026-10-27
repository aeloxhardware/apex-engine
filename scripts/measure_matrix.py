"""
Matrice de mesure step9b : un modele x une variante de parametres, une ligne
JSON ajoutee a un fichier de resultats.

Sert a mesurer les mecanismes opt-in sur les QUATRE modeles de
res/test_models/ avant toute activation par defaut (regle de CLAUDE.md), et a
refaire la reference step9b devenue obsolete (limites #57, #58).

Usage (venv d'Apex, APEX_USE_CPP=1 conseille) :
  .venv/Scripts/python.exe -u scripts/measure_matrix.py <modele.stl> <variante> [--resolution=0.8]
      [--layer-height=0.2] [--out=out/measure_matrix.jsonl]

Variantes : voir apply_variant (ref, median9, ppn, corder, nocorder, heat_t1, noreach, split3, ...).
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from main import PipelineParams  # noqa: E402
from pipeline.step9b_support_check import SupportCheckParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402


def apply_variant(params: PipelineParams, name: str) -> None:
    if name == "ref":
        return
    if name == "legacy_runs":       # coupure a la moindre oscillation (avant run_split_min_amplitude)
        params.toolpath2d.run_split_min_amplitude = 1e-12
        return
    if name == "placed":            # modele pose a z = 0 a l'ingestion (2026-10-13)
        params.place_on_platform = True
        return
    if name == "ref_noradial":      # reference historique (avant l'activation par defaut, limite #64)
        params.reachable_tilt_from_machine = False
        return
    elif name == "median9":
        params.machine_toolpath.free_axis_median_window = 9
    elif name == "ppn":
        params.machine_toolpath.per_point_normal = True
    elif name == "ppn_nodal":
        params.machine_toolpath.per_point_normal = True
        params.machine_toolpath.per_point_normal_nodal = True
    elif name == "continuous":      # parois continues, Step6Params.continuous_rings (2026-10-01)
        params.toolpath2d.continuous_rings = True
    elif name in ("native_fill", "continuous_native"):   # pas 15 : parois/remplissage Clipper2
        from pipeline import step6_toolpath2d
        step6_toolpath2d._NATIVE_FILL = True
        params.toolpath2d.continuous_rings = True
    elif name == "corder":          # step5b : ordre contraint par la tete
        from pipeline.step5b_collision_order import CollisionOrderParams
        params.collision_order = CollisionOrderParams()
    elif name == "corder_support":  # step5b + contraintes d'appui
        from pipeline.step5b_collision_order import CollisionOrderParams
        params.collision_order = CollisionOrderParams(support_depth=0.3)
    elif name == "nocorder":        # sans l'ordre contraint par la tete
        params.collision_order = None
    elif name == "heat_t1":         # temps de diffusion de Crane (h^2) au lieu de 100 h^2
        params.field_solver.heat_method_time_factor = 1.0
    elif name == "noreach":         # flux de la chaleur non projete sur les inclinaisons atteignables
        params.field_solver.heat_project_reachable = False
    elif name == "split3":          # couches scindees en morceaux disjoints
        params.layers.split_disjoint_clusters = 3.0
    else:
        raise SystemExit(f"variante inconnue : {name}")


def _layer_thickness(result) -> float:
    """Epaisseur reelle mediane entre couches consecutives (distance mediane
    entre centroides de contours d'isovaleurs voisines) - pour egaliser les
    comparaisons quand une variante change l'espacement (limites #52, #55)."""
    iso = {}
    for L in result.layers:
        if L.contour_3d:
            pts = np.vstack([np.asarray(c) for c in L.contour_3d])
            iso.setdefault(round(float(L.isovalue), 9), []).append(pts)
    keys = sorted(iso)
    if len(keys) < 3:
        return float("nan")
    from scipy.spatial import cKDTree
    d = []
    for a, b in zip(keys[:-1], keys[1:]):
        pa = np.vstack(iso[a])
        pb = np.vstack(iso[b])
        if len(pa) > 2000:
            pa = pa[np.linspace(0, len(pa) - 1, 2000).astype(int)]
        dist, _ = cKDTree(pb).query(pa)
        d.append(float(np.median(dist)))
    return float(np.median(d))


def _rodrigues_batch(v, axis, ang):
    axis = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    c, sn = np.cos(ang)[:, None], np.sin(ang)[:, None]
    return v * c + np.cross(axis, v) * sn + np.outer(v @ axis, axis) * (1 - c)


def head_normal_deviation(result, params) -> dict:
    """Angle (deg) entre la direction de la tete (cinematique directe des
    angles resolus) et la normale LOCALE de couche (+grad(phi) du tetraedre le
    plus proche), sur les points en extrusion. Metrique PHYSIQUE, independante
    de la cible que step8 s'est donnee - un compteur de diagnostics peut
    baisser sans que le depot s'ameliore."""
    from pipeline.step8_machine_toolpath import gradient_direction_field
    pts = [p for p in result.toolpath.points if p.extrusion > 0]
    if not pts:
        return {}
    axes = params.machine_profile.rotation_axes
    ang = np.radians(np.array([p.axis_angles for p in pts], dtype=float))
    d = np.tile([0.0, 0.0, 1.0], (len(pts), 1))
    for k in range(len(axes) - 1, -1, -1):
        d = _rodrigues_batch(d, axes[k].axis, ang[:, k])
    n = gradient_direction_field(result.tet_mesh, result.scalar_field)(np.array([p.position for p in pts]))
    ok = np.all(np.isfinite(n), axis=1)
    cosang = np.clip(np.einsum("ij,ij->i", d[ok], n[ok]), -1.0, 1.0)
    dev = np.degrees(np.arccos(cosang))
    return {"dev_median": round(float(np.median(dev)), 2), "dev_p95": round(float(np.percentile(dev, 95)), 2),
            "dev_gt15_pct": round(100 * float((dev > 15).mean()), 2),
            "dev_gt45_pct": round(100 * float((dev > 45).mean()), 2)}


def print_time_estimate(result, params, speed=25.0, travel=80.0) -> dict:
    """Temps d'impression estime (heures) : nominal (longueur / vitesse) et
    reel, avec le ralentissement du firmware aux limites de B
    (file_io.gcode_export.feed_plan, limite #59). Mesure le cout des sauts
    de B, qu'un compteur seul ne chiffre pas."""
    from file_io.gcode_export import feed_plan
    pts = result.toolpath.points
    if len(pts) < 2:
        return {}
    pos = np.array([p.position for p in pts], dtype=float)
    ang = np.array([p.axis_angles for p in pts], dtype=float)
    ext = np.array([p.extrusion > 0 for p in pts])
    m = params.machine_profile
    _, _, slow = feed_plan(pos, ang, float(np.degrees(m.max_angular_velocity or 0.0)),
                           float(np.degrees(m.max_angular_acceleration or 0.0)),
                           list(m.position_coupled_axes or []), speed)
    seg = np.linalg.norm(np.diff(pos, axis=0), axis=1)
    pr = ext[1:]
    nominal = seg[pr].sum() / speed + seg[~pr].sum() / travel
    real = (seg[pr] / (speed * np.maximum(slow[pr], 1e-9))).sum() + seg[~pr].sum() / travel
    return {"time_nominal_h": round(nominal / 3600, 3), "time_real_h": round(real / 3600, 3)}


def main() -> None:
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    model, variant = args[0], args[1]
    out = flags.get("--out", os.path.join(ROOT, "out", "measure_matrix.jsonl"))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    params = PipelineParams(resolution=float(flags.get("--resolution", 0.8)))
    params.layers.layer_height_base = float(flags.get("--layer-height", 0.2))
    params.support_check = SupportCheckParams()
    apply_variant(params, variant)
    t0 = time.time()
    result = run_slicing_only(model, params)
    elapsed = time.time() - t0
    rep = result.support_check_report
    pts = result.toolpath.points
    travel = [float(np.linalg.norm(pts[i].position - pts[i - 1].position))
              for i in range(1, len(pts)) if pts[i].extrusion <= 0.0]
    extruded = sum(float(np.linalg.norm(pts[i].position - pts[i - 1].position))
                   for i in range(1, len(pts)) if pts[i].extrusion > 0.0)
    b = np.array([p.axis_angles[1] for p in pts]) if pts and len(pts[0].axis_angles) > 1 else np.zeros(0)
    b_jumps = int((np.abs(np.diff(b)) > 10.0).sum()) if len(b) > 1 else 0
    row = {
        "model": os.path.basename(model), "variant": variant, "resolution": params.resolution,
        "layer_height_base": params.layers.layer_height_base,
        "apex_use_cpp": os.environ.get("APEX_USE_CPP", "0"),
        "seconds": round(elapsed, 1), "layers": len(result.layers), "points": len(pts),
        "layer_thickness_measured": round(_layer_thickness(result), 4),
        "tested": rep.tested_extruded_points, "unsupported": len(rep.events),
        "unsupported_fraction": rep.unsupported_fraction,
        "ordering_defects": rep.ordering_defect_count, "tilt_resolved": rep.tilt_resolved_count,
        "true_overhang": rep.true_overhang_count,
        "extruded_mm": round(extruded, 1),
        "travel_total_mm": round(sum(travel), 1), "travel_max_mm": round(max(travel, default=0.0), 2),
        "step8_diagnostics": len(result.step8_diagnostics or []),
        "b_jumps_gt10": b_jumps,
        **head_normal_deviation(result, params),
        **print_time_estimate(result, params),
        "time": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    if "--collisions" in sys.argv and params.machine_profile.nozzle_model is not None:
        from pipeline import step9_collision_sim as s9
        t1 = time.time()
        rep9 = s9.run(result.toolpath, params.machine_profile, params.collision)
        sampled = len(s9._sample_point_indices(pts, params.collision.path_sampling_step))
        hit = {e.point_index for e in rep9.events}
        row["collision_points"] = len(hit)
        row["collision_fraction"] = round(len(hit) / max(sampled, 1), 4)
        row["collision_platform"] = sum(e.obstacle_type == "platform" for e in rep9.events)
        row["collision_seconds"] = round(time.time() - t1, 1)
    with open(out, "a", encoding="utf-8") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row, indent=1))


if __name__ == "__main__":
    main()
