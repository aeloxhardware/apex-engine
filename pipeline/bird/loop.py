"""
Boucle automatique de la methode Bird (porte de research/S4_Slicer/bird_loop.py).

Conception de l'utilisateur (2026-09-26) :
  - MAX_OVERHANG diminue au fil du bareme (30, 20, 10, 5, 0) tant que chaque
    passe change significativement la piece ;
  - si une passe ne change plus grand-chose, ROTATION_MULTIPLIER augmente ;
  - UN SEUL parametre change par passe ; un repli rejette la passe et essaie la
    valeur a mi-chemin entre la derniere reussie et celle qui a echoue ;
  - bareme epuise mais la piece bouge encore : on repete la meme passe ;
  - fin : repli meme au plus petit ecart, multiplicateur au maximum, ou plafond.
Resultat valide visuellement par l'utilisateur sur le Benchy (passe 4).
LIMITE : seuils cales sur deux modeles au plus (repli 30 mm2, decrochage 25 %).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from pipeline.bird import deformation as bd
from pipeline.bird.fold_judge import as_judged_mesh, judge


@dataclass
class BirdLoopParams:
    overhang_schedule: tuple[float, ...] = (30.0, 20.0, 10.0, 5.0, 0.0)
    mult_start: float = 1.0
    mult_step: float = 0.5
    mult_max: float = 4.0
    fold_zone_mm2: float = 30.0
    stall_ratio: float = 0.25
    min_overhang_step_deg: float = 1.25
    min_mult_step: float = 0.125
    max_passes: int = 16
    layer_height: float = 0.2              # aplatissement du dessous et base plane
    pin_base: bool = True                  # base imprimee a plat (planar_base_mask)
    base_contact_layers: Optional[int] = None   # variante « contact + marge » (None = jusqu'au 1er surplomb)
    max_tet_volume: Optional[float] = None
    # Maillage trop grossier (2026-10-07) : la deformation se calcule par tetraedre, et
    # TetGen par defaut ne raffine pas un STL de quelques faces (l_shape : 30 tetraedres).
    # Si le volume moyen d'un tetraedre depasse `max_mean_tet_volume` (mm3), on refait la
    # tetraedrisation avec `refine_tet_volume`. References de la recherche : Benchy
    # 0,63 mm3 en moyenne (inchange par cette regle), real_overhang deforme a 0,5 mm3.
    # None = desactive.
    max_mean_tet_volume: Optional[float] = 2.0
    refine_tet_volume: float = 0.5
    reachable_tilt: bool = False           # rotations positives interdites (MAX_POS_ROTATION_DEG = 0)
    reachable_normal_weight: float = 0.0   # contrainte de normale atteignable (0 = desactivee)
    offset_xy: tuple[float, float] = (0.0, 0.0)
    deformation_params: dict = field(default_factory=lambda: dict(bd.DEFAULT_PARAMS))


@dataclass
class BirdLoopResult:
    original: bd.BirdTet                   # maillage d'origine (repere Bird : axe du plateau en x = y = 0)
    final: bd.BirdTet                      # meme connectivite, points deformes de la derniere passe acceptee
    passes: list[np.ndarray]               # points de chaque passe acceptee (passes[0] = origine)
    topology: bd.BirdTopology
    base_node_mask: Optional[np.ndarray]
    lowest_overhang_z: Optional[float]
    stop_reason: str
    report: list[dict]
    final_max_overhang: Optional[float]
    final_rotation_multiplier: float


def run(stl_path: str, params: BirdLoopParams = BirdLoopParams(),
        log: Optional[Callable[[str], None]] = None) -> BirdLoopResult:
    say = log or (lambda _msg: None)
    tet = bd.load_and_tetrahedralize(stl_path, max_tet_volume=params.max_tet_volume,
                                     flatten_bottom=params.layer_height)
    if params.max_tet_volume is None and params.max_mean_tet_volume is not None:
        p = tet.points[tet.cells]
        volume = float(np.abs(np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                                        np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))).sum() / 6.0)
        if volume / max(tet.number_of_cells, 1) > params.max_mean_tet_volume:
            say(f"maillage grossier ({tet.number_of_cells} tetraedres, {volume / tet.number_of_cells:.1f} mm3 "
                f"en moyenne) : retetraedrisation a {params.refine_tet_volume} mm3")
            tet = bd.load_and_tetrahedralize(stl_path, max_tet_volume=params.refine_tet_volume,
                                             flatten_bottom=params.layer_height)
    if any(params.offset_xy):
        tet = tet.with_points(tet.points + np.array([params.offset_xy[0], params.offset_xy[1], 0.0]))
    topo, attrs = bd.calculate_tet_attributes(tet)
    original = tet

    planar_base, z_first, node_mask = None, None, None
    if params.pin_base:
        max_height = (params.base_contact_layers * params.layer_height) if params.base_contact_layers else None
        cell_mask, node_mask, z_first = bd.planar_base_mask(tet, topo, params.overhang_schedule[0],
                                                            params.layer_height, max_height=max_height)
        planar_base = (cell_mask, node_mask, tet.points.copy())
        say(f"base figee : {cell_mask.sum()} cellules, {node_mask.sum()} noeuds")
    say(f"{tet.number_of_cells} tetraedres")

    schedule = params.overhang_schedule
    ok_ovh, ok_mult = None, params.mult_start
    next_ovh_index = 0
    try_ovh, try_mult = schedule[0], params.mult_start
    changed = "overhang"
    biggest_move = 0.0
    accepted = 0
    report: list[dict] = []
    passes = [tet.points.copy()]
    stop_reason = f"plafond de {params.max_passes} passes atteint"

    for _ in range(params.max_passes):
        p = {**params.deformation_params, "MAX_OVERHANG": try_ovh, "ROTATION_MULTIPLIER": try_mult}
        if params.reachable_tilt:
            p["MAX_POS_ROTATION_DEG"] = 0
        if params.reachable_normal_weight:
            p["REACHABLE_NORMAL_WEIGHT"] = params.reachable_normal_weight
        new_tet, new_attrs, info = bd.run_one_iteration(tet, topo, attrs, p, planar_base=planar_base)
        move = float(np.median(np.linalg.norm(new_tet.points - tet.points, axis=1)))
        surface = bd.surface_mesh(new_tet, topo, recentre=True)
        fold = judge(as_judged_mesh(surface.vertices, surface.faces))
        entry = {"attempt": len(report) + 1, "max_overhang": try_ovh, "rotation_multiplier": try_mult,
                 "changed": changed, "median_move_mm": move, **fold}

        if fold["largest_zone_mm2"] > params.fold_zone_mm2:
            if changed == "repeat":
                entry["decision"] = "REPLI sur une passe repetee -- rien a reculer"
                report.append(entry)
                stop_reason = "repli en repetant les derniers reglages"
                break
            if changed == "overhang":
                base = ok_ovh if ok_ovh is not None else schedule[0] + 10.0
                mid, too_small = (base + try_ovh) / 2.0, abs(base - try_ovh) / 2.0 < params.min_overhang_step_deg
                try_ovh = mid
            else:
                mid, too_small = (ok_mult + try_mult) / 2.0, abs(try_mult - ok_mult) / 2.0 < params.min_mult_step
                try_mult = mid
            if too_small:
                entry["decision"] = f"REPLI -- ecart minimal atteint sur {changed}"
                report.append(entry)
                stop_reason = f"repli meme avec le plus petit ecart sur {changed}"
                break
            entry["decision"] = f"REPLI -- passe rejetee, essai a {changed}={mid:g}"
            report.append(entry)
            say(entry["decision"])
            continue

        stalled = accepted > 0 and move < params.stall_ratio * biggest_move
        accepted += 1
        tet, attrs = new_tet, new_attrs
        passes.append(tet.points.copy())
        ok_ovh, ok_mult = try_ovh, try_mult
        biggest_move = max(biggest_move, move)
        if changed == "overhang" and next_ovh_index < len(schedule) and try_ovh == schedule[next_ovh_index]:
            next_ovh_index += 1
        entry["pass"] = accepted

        if stalled:
            if try_mult + params.min_mult_step > params.mult_max:
                entry["decision"] = "acceptee -- multiplicateur au maximum"
                report.append(entry)
                stop_reason = f"multiplicateur au maximum ({params.mult_max})"
                break
            changed, try_mult = "multiplier", min(try_mult + params.mult_step, params.mult_max)
            entry["decision"] = f"acceptee, peu de changement -> multiplicateur {try_mult:g}"
        elif next_ovh_index >= len(schedule):
            changed = "repeat"
            entry["decision"] = "acceptee, bareme epuise mais ca bouge encore -> meme passe"
        else:
            changed, try_ovh = "overhang", schedule[next_ovh_index]
            entry["decision"] = f"acceptee -> MAX_OVERHANG {try_ovh:g}"
        report.append(entry)
        say(f"passe {accepted} : {entry['decision']} (depl. {move:.3f} mm, zone {fold['largest_zone_mm2']:.1f} mm2)")

    return BirdLoopResult(original=original, final=tet, passes=passes, topology=topo, base_node_mask=node_mask,
                          lowest_overhang_z=z_first, stop_reason=stop_reason, report=report,
                          final_max_overhang=ok_ovh, final_rotation_multiplier=ok_mult)
