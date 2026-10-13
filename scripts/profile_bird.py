"""
Profil par phase de la voie Bird (`main.run_bird_pipeline` + export), phase 0 du
plan « Bird entier en C++, puis Vulkan » (2026-10-13).

Chronometre en temps REEL (perf_counter) les fonctions de chaque phase, en les
enveloppant a l'execution : pas le surcout x2 de cProfile, et les appels en
sous-processus (pools de P2/P3) sont comptes a leur duree vue par l'appelant.
Une fonction recursive ou rappelee par elle-meme n'est comptee qu'au niveau le
plus externe. Les noms importes ailleurs (`from x import f`) sont remplaces
partout ou ils pointent vers la meme fonction.

Sortie : tableau trie par temps (les postes > 1 s marques), et une ligne JSON
dans --out (defaut out/bird/profile_bird.jsonl).

Usage (garde __main__ obligatoire sous Windows : P2 et P3 lancent des processus) :
  APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/profile_bird.py <modele.stl> [--out=...]
      [--no-skin] [--no-export] [--from-pass=<prefixe>,<k>] [--gcode=<fichier>]
"""
import functools
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# (module, attribut, phase) : attribut "Classe.methode" pour une methode.
TARGETS = [
    ("pipeline.bird.loop", "run", "boucle"),
    ("pipeline.bird.deformation", "load_and_tetrahedralize", "boucle/maillage"),
    ("pipeline.bird.deformation", "build_topology", "boucle/topologie"),
    ("pipeline.bird.deformation", "update_tet_attributes", "boucle/attributs (Dijkstra 1)"),
    ("pipeline.bird.deformation", "calculate_path_length_to_base_gradient", "boucle/gradient (Dijkstra 2, plans)"),
    ("pipeline.bird.deformation", "optimize_rotations", "boucle/LSQ rotations"),
    ("pipeline.bird.deformation", "calculate_deformation", "boucle/LSQ deformation"),
    ("pipeline.bird.deformation", "surface_mesh", "boucle/surface"),
    ("pipeline.bird.fold_judge", "judge", "boucle/juge des replis"),
    ("pipeline.bird.skin.chain", "run", "lissage"),
    ("pipeline.bird.deformed_slicing", "run", "lissage/tranchage deforme"),
    ("pipeline.step5_layers", "run", "lissage/tranchage/etape 5"),
    ("pipeline.step6_toolpath2d", "run", "lissage/tranchage/etape 6"),
    ("pipeline.bird.deformed_slicing", "DeformedToOriginal._candidates", "lissage/tranchage/candidats barycentriques"),
    ("pipeline.bird.deformed_slicing", "DeformedToOriginal._viterbi", "lissage/tranchage/Viterbi"),
    ("pipeline.bird.skin.zones", "compute", "lissage/P2 zones"),
    ("pipeline.bird.skin.zones", "last_accessible_moment_parallel", "lissage/P2/dernier moment accessible"),
    ("pipeline.bird.skin.zones", "split_along_level", "lissage/P2/decoupe par niveau"),
    ("pipeline.bird.skin.surface_paths", "merge_micro_zones", "lissage/P2/fusion des micro-zones"),
    ("pipeline.bird.skin.chain", "carve_hook", "lissage/P3 creusement (construction)"),
    ("pipeline.bird.skin.paths", "trace", "lissage/P4 trajectoires"),
    ("pipeline.bird.skin.surface_paths", "lines_multi", "lissage/P4/lignes de niveau"),
    ("pipeline.bird.skin.paths", "dedupe_shared_passes", "lissage/P4/doublons"),
    ("pipeline.bird.skin.print_order", "order_and_check", "lissage/P5 ordre + collisions"),
    ("pipeline.bird.skin.print_order", "head_angles", "lissage/P5/angles de tete"),
    ("pipeline.bird.skin.print_order", "feed_plan", "lissage/P5/planification d'avance"),
    ("pipeline.bird.skin.print_order", "_check_native", "lissage/P5/controle natif"),
    ("pipeline.bird.machine_toolpath", "run_from_passes", "etape 8 Bird"),
    ("pipeline.bird.machine_toolpath", "run", "etape 8 Bird"),
    ("pipeline.bird.machine_toolpath", "locate", "etape 8 Bird/localisation"),
    ("pipeline.step8_machine_toolpath", "resolve_raw_points", "etape 8 Bird/IK"),
    ("main", "clamp_head_height", "borne de B"),
    ("main", "collision_run", "etape 9"),
    ("main", "reorient_coupled_tiers", "paliers de reorientation"),
    ("main", "adjustment_run", "etape 10"),
    ("main", "export_gcode_checked", "export"),
    ("file_io.gcode_export", "material_travel_check", "export/controle des degagements (index)"),
]

STATS = {}
LSQ = {}      # compteurs internes du solveur natif (bird_lsq_last_profile), cumules par phase


def _lsq_profile():
    from core import _cpp
    if not _cpp.enabled() or not hasattr(_cpp.module(), "bird_lsq_last_profile"):
        return None
    return _cpp.module().bird_lsq_last_profile()


def _wrap(fn, label):
    depth = {"n": 0}

    @functools.wraps(fn)
    def timed(*args, **kwargs):
        depth["n"] += 1
        t0 = time.perf_counter()
        try:
            return fn(*args, **kwargs)
        finally:
            depth["n"] -= 1
            if depth["n"] == 0:
                s = STATS.setdefault(label, {"seconds": 0.0, "calls": 0})
                s["seconds"] += time.perf_counter() - t0
                s["calls"] += 1
                if "LSQ" in label:
                    prof = _lsq_profile()
                    if prof is not None:
                        acc = LSQ.setdefault(label, {})
                        for k, v in prof.items():
                            if k in ("rows", "cols", "nnz"):
                                acc[k] = v
                            else:
                                acc[k] = acc.get(k, 0) + v
    timed.__profile_bird_original__ = fn
    return timed


def install():
    import importlib
    import main  # noqa: F401  (charge l'orchestrateur et ses imports)
    for mod_name, attr, label in TARGETS:
        try:
            mod = importlib.import_module(mod_name)
        except ImportError:
            continue
        owner, name = mod, attr
        if "." in attr:
            cls_name, name = attr.split(".")
            owner = getattr(mod, cls_name, None)
            if owner is None:
                continue
        fn = getattr(owner, name, None)
        if fn is None or hasattr(fn, "__profile_bird_original__"):
            continue
        wrapped = _wrap(fn, label)
        setattr(owner, name, wrapped)
        if owner is mod:   # remplacer aussi les `from mod import fn` deja faits ailleurs
            for other in list(sys.modules.values()):
                d = getattr(other, "__dict__", None)
                if d is None or other is mod:
                    continue
                for key, value in list(d.items()):
                    if value is fn:
                        setattr(other, key, wrapped)


def main_cli():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    model = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    out = flags.get("--out", os.path.join(ROOT, "out", "bird", "profile_bird.jsonl"))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    install()
    import tempfile

    import main as M
    params = M.PipelineParams()
    bird = M.BirdPipelineParams()
    if "--no-skin" in sys.argv:
        bird.skin = None
    loop_result = None
    if "--from-pass" in flags:      # <prefixe>,<k> : passe deja calculee, boucle non relancee
        import numpy as np
        from pipeline.bird import loop as bird_loop
        prefix, k = flags["--from-pass"].rsplit(",", 1)
        cells, p0, pk = (np.load(f"{prefix}{s}.npy") for s in ("cells", "pass0", f"pass{k}"))
        loop_result = bird_loop.from_saved_passes(cells, p0, pk, reason=f"passe {k} enregistree")
    t0 = time.perf_counter()
    r = M.run_bird_pipeline(model, params, bird, loop_result=loop_result)
    if "--no-export" not in sys.argv:
        if "--gcode" in flags:          # G-code conserve (controle d'identite)
            M.export_gcode_checked(r.adjusted_toolpath, params, flags["--gcode"])
        else:
            with tempfile.TemporaryDirectory() as tmp:
                M.export_gcode_checked(r.adjusted_toolpath, params, os.path.join(tmp, "x.gcode"))
    total = time.perf_counter() - t0

    sizes = {"tetrahedra": int(len(r.loop.original.cells)), "nodes": int(len(r.loop.original.points)),
             "accepted_passes": len(r.loop.passes) - 1, "toolpath_points": len(r.adjusted_toolpath.points)}
    print(f"\n{os.path.basename(model)} : {total:.1f} s au total ; {sizes}")
    print(f"{'phase':58s} {'s':>9s} {'appels':>7s}")
    for label, s in sorted(STATS.items(), key=lambda kv: -kv[1]["seconds"]):
        mark = "  <-- > 1 s" if s["seconds"] > 1.0 else ""
        print(f"{label:58s} {s['seconds']:9.2f} {s['calls']:7d}{mark}")
    for label, acc in LSQ.items():
        tot = acc.get("total_s", 0.0) or 1e-12
        print(f"\n{label} (natif) : {acc.get('total_s', 0):.2f} s ; produits J.x {acc.get('matvec_s', 0):.2f} s "
              f"({100 * acc.get('matvec_s', 0) / tot:.0f} %), J^T.u {acc.get('rmatvec_s', 0):.2f} s "
              f"({100 * acc.get('rmatvec_s', 0) / tot:.0f} %), residus {acc.get('fun_s', 0):.2f} s, "
              f"jacobienne {acc.get('jac_s', 0):.2f} s ; nfev {acc.get('nfev', 0)}, iterations LSMR "
              f"{acc.get('lsmr_iterations', 0)} ; J {acc.get('rows', 0)} x {acc.get('cols', 0)}, nnz {acc.get('nnz', 0)}")
    row = {"model": os.path.basename(model), "total_s": round(total, 2), "apex_use_cpp": os.environ.get("APEX_USE_CPP", "0"),
           "skin": bird.skin is not None, **sizes,
           "phases": {k: {"seconds": round(v["seconds"], 3), "calls": v["calls"]} for k, v in STATS.items()},
           "lsq": LSQ, "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(out, "a", encoding="utf-8") as f:
        f.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main_cli()
