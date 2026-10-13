"""
Phase 3 du plan bird-into-apex-cpp : voie Bird contre voie FEM sur un modele.

Une ligne JSON par voie dans --out (defaut out/bird/phase3.jsonl) :
  - step9b (matiere deposee sans appui, defauts d'ordre, surplombs reels), mesure sur
    la trajectoire AVANT l'etape 10 dans les deux voies (comme measure_matrix.py) ;
  - collisions de la vraie tete apres l'etape 9 (Bird : apres les paliers de
    reorientation, que la voie FEM n'active pas par defaut) et non resolues apres l'etape 10 ;
  - orientations inatteignables (Step8Diagnostic), points de trajectoire ;
  - temps de calcul, temps d'impression estime (tous axes, Z recommande), tours de Y.

Usage (venv d'Apex, APEX_USE_CPP=1) :
  .venv/Scripts/python.exe -u scripts/measure_bird_vs_fem.py <modele.stl> [--ways=fem,bird] [--out=...]
      [--keep-gcode=<dossier>]   (G-code conserves : <modele>_<voie>.gcode)
"""
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def gcode_metrics(toolpath, params, keep_path=None):
    import tempfile
    import main as M
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        _lines, rep = M.export_gcode_checked(toolpath, params, keep_path or os.path.join(tmp, "x.gcode"))
    return {"moves": rep.moves, "lines": rep.lines, "retractions": rep.retractions,
            "clearance_axes": rep.clearance_axes, "travel_collisions": rep.travel_collisions,
            "print_time_h": rep.print_time_s / 3600.0, "print_time_fw_h": rep.print_time_fw_s / 3600.0,
            "y_turns": rep.y_turns, "tip_z_min": rep.tip_z_min, "head_z_min": rep.head_z_min,
            "filament_m": rep.filament_mm / 1000.0, "written": rep.written, "export_s": round(time.time() - t0, 1)}


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    model = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    ways = flags.get("--ways", "fem,bird").split(",")
    out = flags.get("--out", os.path.join(ROOT, "out", "bird", "phase3.jsonl"))
    os.makedirs(os.path.dirname(out), exist_ok=True)

    import main as M
    from pipeline.step9b_support_check import SupportCheckParams, run as support_run

    for way in ways:
        params = M.PipelineParams()
        t0 = time.time()
        row = {"model": os.path.basename(model), "way": way, "resolution": params.resolution}
        if way == "fem":
            params.support_check = SupportCheckParams()
            r = M.run_pipeline(model, params)
            sc = r.support_check_report
            toolpath, events, adj_rep, diag = r.toolpath, r.collision_report.events, r.adjustment_report, \
                r.step8_diagnostics
            row["layers"] = len(r.layers)
        else:
            loop_result = None
            if "--from-pass" in flags:      # <prefixe>,<k> : loop_<m>_pin_ de la recherche, passe k
                import numpy as np
                from pipeline.bird import loop as bird_loop
                prefix, k = flags["--from-pass"].rsplit(",", 1)
                cells, p0, pk = (np.load(f"{prefix}{s}.npy") for s in ("cells", "pass0", f"pass{k}"))
                loop_result = bird_loop.from_saved_passes(cells, p0, pk, reason=f"passe {k} enregistree")
                row["from_pass"] = flags["--from-pass"]
            r = M.run_bird_pipeline(model, params, loop_result=loop_result)
            toolpath, events, adj_rep, diag = r.toolpath, r.collision_report.events, r.adjustment_report, \
                r.step8_diagnostics
            sc = support_run(toolpath, params.machine_profile, params.layers.layer_height_base, 0.0,
                             SupportCheckParams())
            row["tetrahedra"] = int(len(r.loop.original.cells))
            row["accepted_passes"] = len(r.loop.passes) - 1
            row["stop_reason"] = r.loop.stop_reason
            body = getattr(r.slicing, "body", None) or getattr(r.slicing, "layers", [])
            row["layers"] = len(body)
            if hasattr(r.slicing, "timings"):
                row["skin_timings"] = r.slicing.timings
                row["skin_excluded_zones"] = r.slicing.excluded_zones
        row["compute_s"] = round(time.time() - t0, 1)
        row.update({
            "points": len(toolpath.points),
            "unsupported_fraction": sc.unsupported_fraction,
            "unsupported_events": len(sc.events),
            "ordering_defects": sc.ordering_defect_count,
            "true_overhangs": sc.true_overhang_count,
            "collision_points": len({e.point_index for e in events}),
            "unresolved_after_step10": len(adj_rep.unresolved_collisions),
            # supports de secours : collisions que l'etape 10 n'a pas su reorienter (comptees
            # comme traitees par `converged`, d'ou ce champ separe)
            "fallback_supports": sum(1 for e in adj_rep.entries if e.method == "fallback_support"),
            "local_reorientations": sum(1 for e in adj_rep.entries if e.method == "local_reorientation"),
            "step8_diagnostics": len(diag),
        })
        keep = flags.get("--keep-gcode")
        keep_path = (os.path.join(keep, f"{os.path.splitext(os.path.basename(model))[0]}_{way}.gcode")
                     if keep else None)
        row.update(gcode_metrics(r.adjusted_toolpath, params, keep_path))
        with open(out, "a", encoding="utf-8") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
