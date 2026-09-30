"""
Boucle de rognage local du lissage non planaire (decision utilisateur du
2026-09-28) : P2-bis -> P4 -> P3 -> P5, repetee tant que la VRAIE tete touche
quelque part. A chaque tour, les points de collision de P5 sont AJOUTES au
fichier cumule ; P2-bis retire du lissage les faces a moins de --trim-radius
de ces points (collision objet), les zones sont recalculees et tout est
relance. Arret : zero collision, ou --max-iter tours (les zones encore en
collision sont alors retirees par P5, comme le prevoyait le plan).

Usage (venv d'Apex) :
  ../../.venv/Scripts/python.exe -u skin_loop.py "<modele STL>" <passe> <variante> --first-layer=<motif>
      [--max-iter=3] [--trim-radius=2.0] [--head-tolerance=45] [--pattern=rectilinear]
Journaux : output_models/skin_loop_<etape>_it<n>.log
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def run(step, cmd, it, env):
    log = os.path.join(HERE, "output_models", f"skin_loop_{step}_it{it}.log")
    t0 = time.time()
    with open(log, "w") as f:
        code = subprocess.call(cmd, cwd=HERE, stdout=f, stderr=subprocess.STDOUT, env=env)
    print(f"  {step:<6} {time.time() - t0:6.0f} s  (code {code}, {os.path.basename(log)})", flush=True)
    if code:
        sys.exit(f"echec de {step} : voir {log}")
    return log


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    model, k, variant = args[0], args[1], args[2] if len(args) > 2 else ""
    max_iter = int(flags.get("--max-iter", 3))
    py = sys.executable
    env = dict(os.environ, APEX_USE_CPP="1")
    name = "loop_" + model.replace(" ", "_") + variant
    prefix = os.path.join(HERE, "output_models", name)
    cumulative = f"{prefix}_skin_trim_pass{k}.json"
    json.dump({"points": []}, open(cumulative, "w"))
    p2_extra = [f"--trim-radius={flags.get('--trim-radius', '2.0')}",
                f"--head-tolerance={flags.get('--head-tolerance', '45')}", f"--skin-trim={cumulative}"]
    # Refonte du 2026-09-30 (--p3-each-turn pour l'ancien comportement) : pendant
    # la boucle, P5 teste contre le corps NON creuse (sur-ensemble de la matiere,
    # plus prudent) et P3 ne tourne qu'UNE fois, a la fin, suivi d'un P5 de
    # verification sur le vrai corps creuse.
    each_turn = "--p3-each-turn" in sys.argv
    tol = f"--head-tolerance={flags.get('--head-tolerance', '45')}"
    p3 = [py, "-u", "bird_slice_deformed.py", model, k, variant, "--skin", f"--first-layer={flags['--first-layer']}"]
    p5 = [py, "-u", "skin_print.py", model, k, variant, tol]
    clean = False
    for it in range(1, max_iter + 1):
        print(f"tour {it}", flush=True)
        run("P2", [py, "-u", "surface_skin.py", model, k, variant, *p2_extra], it, env)
        run("P4", [py, "-u", "skin_paths.py", name[len("loop_"):], k,
                   f"--pattern={flags.get('--pattern', 'rectilinear')}"], it, env)
        if each_turn:
            run("P3", p3, it, env)
        # pendant la boucle : un point tous les 2 mm (24 s au lieu de 240 s sur
        # le dino) ; le P5 FINAL, sur le corps creuse, teste tous les points
        log = run("P5", p5 + ["--sampling=2"] + ([] if each_turn else ["--body-uncarved"]), it, env)
        new = json.load(open(f"{prefix}_skin_collisions_pass{k}.json"))["points"]
        summary = [ln.strip() for ln in open(log) if ln.startswith("zones en collision")]
        print(f"  -> {len(new)} points de collision ; {summary[0] if summary else ''}", flush=True)
        if not new:
            clean = True
            break
        acc = json.load(open(cumulative))
        acc["points"].extend(new)
        json.dump(acc, open(cumulative, "w"))
    if clean:
        if each_turn:
            print("zero collision : termine")
            return
        run("P3", p3, "final", env)
        log = run("P5", p5, "final", env)
        new = json.load(open(f"{prefix}_skin_collisions_pass{k}.json"))["points"]
        if not new:
            print("zero collision (verifie sur le corps creuse) : termine")
            return
        print(f"verification sur le corps creuse : {len(new)} points de collision -> zones retirees")
    else:
        print(f"{max_iter} tours sans zero collision : zones restantes retirees, tranchage refait sans leur creusement")
    # Tour final (plan P5) : les zones encore en collision sont retirees, leur
    # creusement est annule (P3 --skin-exclude, UNE fois), puis verification.
    # Exclusion CUMULEE jusqu'a zero contact (correctif du 2026-09-30) : le P5
    # final est exhaustif et voit des contacts que la boucle (2 mm, corps non
    # creuse) n'a pas vus ; sans nouveau P3, les zones qu'il retire laissaient
    # leur creusement SANS couche de lissage (Squirtle : 5 zones, 106 mm2).
    kept_file = f"{prefix}_skin_excluded_final_pass{k}.json"
    excluded_zones = set(json.load(open(f"{prefix}_skin_excluded_pass{k}.json"))["excluded_zones"])
    final_exclusion_rounds(py, model, k, variant, p3, p5, prefix, kept_file, excluded_zones, env)


def final_exclusion_rounds(py, model, k, variant, p3, p5, prefix, kept_file, excluded_zones, env, rounds=3):
    for r in range(1, rounds + 1):
        json.dump({"excluded_zones": sorted(int(z) for z in excluded_zones)}, open(kept_file, "w"))
        run("P3", p3 + [f"--skin-exclude={kept_file}"], f"final{r}", env)
        run("P5", p5 + [f"--exclude={kept_file}"], f"final{r}", env)
        new = json.load(open(f"{prefix}_skin_collisions_pass{k}.json"))["points"]
        more = set(json.load(open(f"{prefix}_skin_excluded_pass{k}.json"))["excluded_zones"])
        print(f"  -> zones retirees : {sorted(excluded_zones)} ; points de collision restants : {len(new)}", flush=True)
        if not new:
            return
        excluded_zones |= more
    print(f"  ATTENTION : {rounds} tours d'exclusion sans zero contact -- zones en contact non creusees "
          f"mais trous possibles, voir {prefix}_skin_excluded_pass{k}.json")


if __name__ == "__main__":
    main()
