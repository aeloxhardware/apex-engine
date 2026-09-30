"""
Boucle automatique de la methode Bird -- remplace la methode manuelle de
l'utilisateur (une passe, on regarde, on decide).

Conception de l'utilisateur (2026-09-26), avec les deux juges mesures contre
ses propres verdicts :

  - MAX_OVERHANG DIMINUE au fil du bareme (30, 20, 10, 5, 0) tant que chaque
    passe change significativement la piece ;
  - si une passe ne change plus grand-chose, ROTATION_MULTIPLIER AUGMENTE
    (+MULT_STEP) a la passe suivante ;
  - UN SEUL parametre change par passe, pour savoir lequel reculer en cas de
    repli : on jette la passe, on repart du maillage precedent, et on essaie
    la valeur a mi-chemin entre la derniere reussie et celle qui a echoue ;
  - bareme epuise mais la piece bouge encore : on REPETE la meme passe
    (si une repetition replie, rien a reculer : fin) ;
  - fin : quand meme le plus petit ecart replie la piece (on garde la
    derniere passe saine), quand le multiplicateur atteint MULT_MAX, ou au
    plafond de passes.

Juges :
  - REPLI : fold_judge, plus grande zone ou la surface se traverse
    > FOLD_ZONE_MM2. Benchy : zones saines <= 13,2 mm2, plus petit repli
    de l'utilisateur 69,7 mm2.
  - PAS DE CHANGEMENT SIGNIFICATIF : deplacement median de la passe
    < STALL_RATIO x le plus grand deplacement deja vu. Benchy : 0,2 % a
    convergence. real_overhang fin : plateau tremblant a 14-22 %, d'ou 25 %
    -- CALE sur ce modele, pas mesure contre le jugement de l'utilisateur.
Tout est calibre sur deux modeles au plus : a verifier sur death.stl.

La condition d'arret de la conception d'origine (surplomb sous un seuil)
n'est PAS implementee : aucune mesure de surplomb fiable n'existe encore.

Pas de recentrage entre les passes (comme server.py, ou le bareme manuel de
l'utilisateur a ete valide) : l'export recentre une COPIE.

Le juge des replis vit dans le venv d'Apex (pipeline/step1_ingest), la
deformation dans celui de la recherche : il est appele en sous-processus.

BASE PLANE (2026-09-27) : par defaut, tout le volume sous le surplomb le
plus bas de la piece d'origine est FIGE (bird_core.planar_base_mask) --
rotation nulle, noeuds a leur position d'origine a chaque passe. Sans cela,
la face posee sur le plateau se creusait en cuvette dans l'espace deforme
(Benchy : 6,75 mm, coupee par 34 couches imprimees du centre vers le bord).
--no-pin reproduit le comportement d'origine (reference de comparaison).
--base-contact : variante "contact + marge" -- seules les BASE_CONTACT_LAYERS
premieres couches sont figees (toujours sous le premier surplomb), Bird reste
libre au-dessus. Sorties suffixees _pinc.
--reachable-tilt : parametre SELECTIONNABLE (decision de l'utilisateur du
2026-09-29) -- les couches ne penchent que du cote que la tete de la machine
ACTUELLE atteint (B d'un seul cote : normale vers +R). Chaque cellule tourne
autour de l'axe tangentiel d'un angle signe ; la normale de couche dans la
piece vaut R^T z, qui penche vers le CENTRE pour un angle positif. On borne
donc MAX_POS_ROTATION_DEG a 0 (le lissage reste entre les cibles et 0).
Sorties suffixees _reach.
--reachable-normal=<poids> : contrainte portant sur la NORMALE FINALE des
couches (termes ajoutes a la deformation, cf. bird_core
._calculate_deformation_pinned) au lieu du signe de la rotation par cellule.
Selectionnable. Sorties suffixees _rn<poids>.
--offset=<dx>,<dy> : decalage de la piece sur le plateau (mm), pour les tests
(dans Apex, l'interface fixe l'origine du modele assemble au centre du
plateau). Les directions radiales -- donc les inclinaisons atteignables --
en dependent. Sauvegarde dans _shift.npy ; sorties suffixees _off<dx>_<dy>.

Usage (venv de la recherche) :
  .venv/Scripts/python.exe bird_loop.py <modele> [max_tet_volume] [--no-pin | --base-contact]
<modele> : input_models/<modele>.stl, sinon res/test_models/<modele>.stl.
Sorties : output_models/loop_<modele>[_a<vol>][_pin]_pass<k>.{stl,npy},
_final.stl, _report.json (une ligne par tentative) et, base figee,
_base_mask.npy (masque des noeuds figes).
"""
import json
import os
import subprocess
import sys
import time

import numpy as np

import bird_core

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
APEX_PYTHON = os.path.join(HERE, "..", "..", ".venv", "Scripts", "python.exe")

OVERHANG_SCHEDULE = [30.0, 20.0, 10.0, 5.0, 0.0]
MULT_START = 1.0
MULT_STEP = 0.5
MULT_MAX = 4.0
FOLD_ZONE_MM2 = 30.0
STALL_RATIO = 0.25
MIN_OVERHANG_STEP_DEG = 1.25
MIN_MULT_STEP = 0.125
MAX_PASSES = 16
# Base plane : surplomb qui en marque la fin (meme critere que Bird, premier
# seuil du bareme) et marge sous ce surplomb (hauteur de couche du tranchage,
# main.PipelineParams d'Apex -- le venv de la recherche ne peut pas l'importer).
BASE_OVERHANG_DEG = OVERHANG_SCHEDULE[0]
LAYER_HEIGHT = 0.2
# variante "contact + marge" : nombre habituel de couches pleines du dessous
BASE_CONTACT_LAYERS = 3

BASE_PARAMS = {
    "NEIGHBOUR_LOSS_WEIGHT": 20,
    "SET_INITIAL_ROTATION_TO_ZERO": False,
    "INITIAL_ROTATION_FIELD_SMOOTHING": 30,
    "MAX_POS_ROTATION_DEG": 360,
    "MAX_NEG_ROTATION_DEG": -360,
    "ROTATE_ITERATIONS": 100,
    "DEFORM_ITERATIONS": 1000,
    "STEEP_OVERHANG_COMPENSATION": True,
}


def judge_fold(stl_path):
    env = dict(os.environ, APEX_USE_CPP="1")
    out = subprocess.run([APEX_PYTHON, os.path.join(HERE, "fold_judge.py"), "--json", stl_path],
                         capture_output=True, text=True, env=env, check=True)
    return json.loads(out.stdout.strip().splitlines()[-1])


def export(tet, path):
    bird_core.export_surface_stl(tet.copy(deep=True), path)


def model_path(model_name):
    for folder in (os.path.join(HERE, "input_models"), os.path.join(ROOT, "res", "test_models")):
        path = os.path.join(folder, f"{model_name}.stl")
        if os.path.exists(path):
            return path
    raise FileNotFoundError(model_name)


def main():
    shift_arg = next((a for a in sys.argv[1:] if a.startswith("--offset=")), None)
    rn_arg = next((a for a in sys.argv[1:] if a.startswith("--reachable-normal=")), None)
    rn_weight = float(rn_arg.split("=", 1)[1]) if rn_arg else 0.0
    shift = np.zeros(3)
    if shift_arg:
        shift[:2] = [float(x) for x in shift_arg.split("=", 1)[1].split(",")]
    args = [a for a in sys.argv[1:] if a not in ("--no-pin", "--base-contact", "--reachable-tilt")
            and not a.startswith("--offset=") and not a.startswith("--reachable-normal=")]
    reachable = "--reachable-tilt" in sys.argv
    pin = "--no-pin" not in sys.argv
    contact = "--base-contact" in sys.argv
    model_name = args[0] if args else "benchy upsidedown tilted"
    max_tet_volume = float(args[1]) if len(args) > 1 else None
    prefix = os.path.join(HERE, "output_models", "loop_" + model_name.replace(" ", "_")
                          + (f"_a{max_tet_volume:g}" if max_tet_volume else "")
                          + (("_pinc" if contact else "_pin") if pin else "")
                          + ("_reach" if reachable else "")
                          + (f"_rn{rn_weight:g}" if rn_arg else "")
                          + (f"_off{shift[0]:g}_{shift[1]:g}" if shift_arg else ""))

    # dessous aplati sur la premiere couche (voir load_and_tetrahedralize)
    tet = bird_core.load_and_tetrahedralize(model_path(model_name), max_tet_volume=max_tet_volume,
                                            flatten_bottom=LAYER_HEIGHT)
    if shift_arg:
        tet.points = np.asarray(tet.points) + shift
    np.save(f"{prefix}_shift.npy", shift)
    tet, state = bird_core.calculate_tet_attributes(tet)
    planar_base, z_first = None, None
    if pin:
        cell_mask, node_mask, z_first = bird_core.planar_base_mask(
            tet, BASE_OVERHANG_DEG, LAYER_HEIGHT,
            max_height=BASE_CONTACT_LAYERS * LAYER_HEIGHT if contact else None)
        planar_base = (cell_mask, node_mask, np.asarray(tet.points).copy())
        np.save(f"{prefix}_base_mask.npy", node_mask)
        z0 = np.asarray(tet.points)[:, 2].min()
        print(f"base figee : {cell_mask.sum()} cellules, {node_mask.sum()} noeuds, "
              f"jusqu'a {np.asarray(tet.points)[node_mask, 2].max() - z0:.2f} mm du plateau (surplomb le plus bas "
              f"a partir de {z_first - z0:.2f} mm)")
    # Connectivite + positions d'origine : bird_to_apex.py en a besoin pour
    # transporter la hauteur deformee vers le maillage d'Apex.
    np.save(f"{prefix}_cells.npy", np.asarray(tet.field_data["cells"]))
    np.save(f"{prefix}_pass0.npy", np.asarray(tet.points))
    print(f"{model_name} : {tet.number_of_cells} tetraedres")
    print(f"juges : repli si plus grande zone > {FOLD_ZONE_MM2} mm2 ; "
          f"pas de changement si deplacement < {STALL_RATIO:.0%} du plus grand deja vu")
    print(f"{'passe':>5} {'MAX_OVH':>7} {'MULT':>6} {'depl. med.':>11} {'zone max':>10} {'duree':>6}  decision")

    # Derniers reglages ACCEPTES, et valeurs visees pour la prochaine passe.
    ok_ovh, ok_mult = None, MULT_START
    next_ovh_index = 0
    try_ovh, try_mult = OVERHANG_SCHEDULE[0], MULT_START
    changed = "overhang"             # parametre modifie par la tentative en cours
    biggest_move = 0.0
    accepted = 0
    report = []
    stop_reason = f"plafond de {MAX_PASSES} passes atteint"

    for _ in range(MAX_PASSES):
        params = {**BASE_PARAMS, "MAX_OVERHANG": try_ovh, "ROTATION_MULTIPLIER": try_mult}
        if reachable:
            params["MAX_POS_ROTATION_DEG"] = 0
        if rn_weight:
            params["REACHABLE_NORMAL_WEIGHT"] = rn_weight
        t0 = time.time()
        new_tet, info = bird_core.run_one_iteration(tet, state, params, planar_base=planar_base)
        dt = time.time() - t0
        move = float(np.median(np.linalg.norm(np.asarray(new_tet.points) - np.asarray(tet.points), axis=1)))
        candidate = f"{prefix}_try{len(report) + 1}_ovh{try_ovh:g}_m{try_mult:g}.stl"
        export(new_tet, candidate)
        fold = judge_fold(candidate)
        entry = {"attempt": len(report) + 1, "max_overhang": try_ovh, "rotation_multiplier": try_mult,
                 "changed": changed, "median_move_mm": move, "seconds": round(dt, 1), **fold}

        used_ovh, used_mult = try_ovh, try_mult

        def show(label):
            print(f"{entry.get('pass', '-'):>5} {used_ovh:>7g} {used_mult:>6g} {move:>9.3f}mm "
                  f"{fold['largest_zone_mm2']:>7.1f}mm2 {dt:>5.0f}s  {label}")

        if fold["largest_zone_mm2"] > FOLD_ZONE_MM2:
            # Repli : passe rejetee, on recule le SEUL parametre qui a change.
            if changed == "repeat":
                entry["decision"] = "REPLI sur une passe repetee -- rien a reculer"
                report.append(entry)
                show(entry["decision"])
                stop_reason = "repli en repetant les derniers reglages"
                break
            if changed == "overhang":
                base = ok_ovh if ok_ovh is not None else OVERHANG_SCHEDULE[0] + 10.0
                mid, step_too_small = (base + try_ovh) / 2.0, abs(base - try_ovh) / 2.0 < MIN_OVERHANG_STEP_DEG
                try_ovh = mid
            else:
                mid, step_too_small = (ok_mult + try_mult) / 2.0, abs(try_mult - ok_mult) / 2.0 < MIN_MULT_STEP
                try_mult = mid
            if step_too_small:
                entry["decision"] = f"REPLI -- ecart minimal atteint sur {changed}"
                report.append(entry)
                show(entry["decision"])
                stop_reason = f"repli meme avec le plus petit ecart sur {changed}"
                break
            entry["decision"] = f"REPLI -- passe rejetee, essai a {changed}={mid:g}"
            report.append(entry)
            show(entry["decision"])
            continue

        # Passe acceptee.
        stalled = accepted > 0 and move < STALL_RATIO * biggest_move
        accepted += 1
        tet = new_tet
        ok_ovh, ok_mult = try_ovh, try_mult
        biggest_move = max(biggest_move, move)
        if changed == "overhang" and next_ovh_index < len(OVERHANG_SCHEDULE) \
                and try_ovh == OVERHANG_SCHEDULE[next_ovh_index]:
            next_ovh_index += 1
        export(tet, f"{prefix}_pass{accepted}.stl")
        np.save(f"{prefix}_pass{accepted}.npy", np.asarray(tet.points))
        entry["pass"] = accepted

        # Prochaine passe : un seul parametre change.
        if stalled:
            if try_mult + MIN_MULT_STEP > MULT_MAX:
                entry["decision"] = "acceptee -- multiplicateur au maximum"
                report.append(entry)
                show(entry["decision"])
                stop_reason = f"multiplicateur au maximum ({MULT_MAX})"
                break
            changed, try_mult = "multiplier", min(try_mult + MULT_STEP, MULT_MAX)
            entry["decision"] = f"acceptee, peu de changement -> multiplicateur {try_mult:g}"
        elif next_ovh_index >= len(OVERHANG_SCHEDULE):
            changed = "repeat"
            entry["decision"] = "acceptee, bareme epuise mais ca bouge encore -> meme passe"
        else:
            changed, try_ovh = "overhang", OVERHANG_SCHEDULE[next_ovh_index]
            entry["decision"] = f"acceptee -> MAX_OVERHANG {try_ovh:g}"
        report.append(entry)
        show(entry["decision"])

    final = f"{prefix}_final.stl"
    export(tet, final)
    with open(f"{prefix}_report.json", "w", encoding="utf-8") as f:
        json.dump({"model": model_name, "planar_base": pin, "lowest_overhang_z": z_first,
                   "stop_reason": stop_reason, "accepted_passes": accepted,
                   "final_max_overhang": ok_ovh, "final_rotation_multiplier": ok_mult,
                   "fold_zone_mm2": FOLD_ZONE_MM2, "stall_ratio": STALL_RATIO,
                   "attempts": report}, f, indent=2, ensure_ascii=False)
    print()
    print(f"arret : {stop_reason} -- {accepted} passe(s) acceptee(s), "
          f"derniers reglages MAX_OVERHANG={ok_ovh:g}, ROTATION_MULTIPLIER={ok_mult:g}")
    print(f"resultat : {final}")


if __name__ == "__main__":
    main()
