"""
Export G-code pour la machine ACTUELLE (R C Z B, limite #59 de CLAUDE.md).

Ecrit de zero le 2026-09-30, sans reprendre les generateurs de formes fixes
`scripts/gen_*_gcode.py` (decision de l'utilisateur). Seuls les FAITS machine
releves dans ces scripts, verifies avec l'utilisateur, sont repris :

- commande `G00 Y<C> X<R> Z<Z> B<B> E<e> F<f>` ; en-tete `G21 G90 M82 G92 E0
  M109 S<t>` ;
- firmware Y = angle polaire du point, `atan2(y, x)` en degres. Avec la
  convention d'Apex (`coupled_axis_angle`, +R = -Y a angles nuls), c'est
  exactement `C - 90` : l'export le recalcule depuis la position et controle
  l'accord avec `axis_angles[0]` (`GcodeReport.c_mismatch_deg`) ;
- firmware B = B d'Apex (B > 0 penche le CORPS de la tete vers +R, donc la
  buse vers le centre, meme convention que step9._positioned_toolhead_vertices) ;
- le firmware commande le PIVOT de B, a `pivot_to_tip` (46 mm) de la pointe,
  sans compensation : l'export la fait. Reference : a B = 0 les coordonnees
  firmware sont celles de la pointe (calibration buse verticale), donc
  X = R + L sin B et Z = z + L (cos B - 1) ;
- le firmware traite Y (degres) et B comme des axes lineaires : F est
  recalcule par segment pour garder la vitesse de la POINTE constante ;
- le firmware (RRF) ne fait AUCUNE cinematique, ni polaire ni pivot
  (utilisateur, 2026-09-30) : tout est ici.

Modele pivot/pointe PARAMETRE (`tilt_sign`, `z_reference`) : la
previsualisation UX de l'utilisateur, qui fait reference, place des
trajectoires sous le plateau la ou ce modele n'en voit aucune. Convention en
attente de confirmation ; les valeurs par defaut sont celles d'origine.

Les points de trajectoire sont au milieu de la couche (step5) ; la pointe est
emise au SOMMET du cordon (`tip_at_bead_top`, h/2 le long de l'axe de tete).
`GcodeReport` porte des garde-fous (pointe et tete les plus basses, etendues
X/Y/Z firmware) et `strict` refuse d'ecrire un fichier dont la pointe passe
sous le plateau.

Les limites de vitesse et d'acceleration de B (MachineProfile) ralentissent
tout le mouvement, le firmware coordonnant les axes (limite #59) :
`feed_plan` donne ce facteur par segment.

Trajets a vide (decision de l'utilisateur, 2026-09-30) : retraction, puis
degagement DYNAMIQUE sur un seul axe a la fois, choisi selon l'orientation de
la tete : Z si la tete est plutot verticale (+Z tete vers le bas, -Z tete
retournee), sinon R dans le sens du corps de la tete. Un degagement qui ferait
passer la tete sous la marge du plateau bascule sur l'autre axe.

Entree : le dictionnaire des trajectoires machine (`position` (N,3) repere
machine, `axis_angles` (N,2) degres, `extrusion` (N,) -- 0 au premier point
de chaque polyligne, donc le segment qui y mene est un trajet a vide --,
`is_skin` (N,) optionnel), tel qu'ecrit par la recherche
(`full_step89.py --apply`) ; la conversion depuis `ToolpathMachine` est
`toolpath_arrays`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from core import progress as _progress


@dataclass
class GcodeParams:
    nozzle_temp: float = 195.0
    filament_diameter: float = 1.75
    layer_height: float = 0.2            # epaisseur du cordon du corps (mm)
    skin_thickness: float = 0.2          # epaisseur des couches de lissage (mm)
    line_width: float = 0.4
    flow: float = 1.0
    print_speed: float = 25.0            # mm/s a la pointe
    travel_speed: float = 80.0           # mm/s a la pointe
    rotation_feed: float = 3000.0        # F (unites firmware / min) d'un mouvement sans deplacement de pointe
    max_feed: float = 12000.0            # plafond de F
    retract_length: float = 0.8          # mm de filament (utilisateur, 2026-09-30)
    retract_speed: float = 35.0          # mm/s (utilisateur, 2026-09-30)
    retract_min_travel: float = 2.0      # pas de retraction pour un trajet plus court (mm)
    travel_clearance: float = 1.0        # degagement le long d'un seul axe (mm)
    # Si un controle de contact est fourni (`travel_check`), un degagement en contact
    # est repris avec l'autre axe, puis avec un degagement double, quadruple... jusqu'a
    # `travel_clearance * travel_clearance_max_factor` (2026-10-06).
    travel_clearance_max_factor: float = 4.0
    pivot_to_tip: float = 46.0           # pivot de B -> pointe (mm) : 46, mesure de l'utilisateur (2026-10-01 ; 45 auparavant)
    # Modele pivot/pointe, parametre en attendant la confirmation de
    # l'utilisateur (sa previsualisation UX et ce modele divergent, 2026-09-30).
    # tilt_sign : signe du deport radial du pivot quand B > 0 (+1 = corps vers +R).
    # z_reference : "tip" = a B = 0 le Z firmware est celui de la pointe
    # (Z = z + L (cos B - 1)) ; "pivot" = Z firmware est la hauteur du pivot
    # (Z = z + L cos B).
    tilt_sign: float = 1.0
    z_reference: str = "tip"
    # b_sign : B firmware = b_sign * B d'Apex. +1 : B POSITIF penche le corps de la tete vers +R
    # (la buse vers le centre) ; sur la machine, B NEGATIF envoie la pointe vers l'EXTERIEUR
    # (croquis de l'utilisateur, 2026-10-06). Un essai a -1 du meme jour reposait sur un
    # contresens (« B negatif = vers l'avant » decrit le sens de l'axe, pas une inversion
    # a appliquer) : corrige. Toute la geometrie reste en B d'Apex ; seul le B emis peut changer.
    b_sign: float = 1.0
    # Les points de trajectoire sont au MILIEU de la couche (step5 : premiere
    # isovaleur a v_min + delta/2) : la pointe doit etre au SOMMET du cordon,
    # donc relevee de h/2 le long de l'axe de la tete.
    tip_at_bead_top: bool = True
    strict: bool = False                 # refuse d'ecrire si la pointe passe sous le plateau
    # Cumul de rotation de Y (plateau) tolere, en tours ; None = aucune limite. Utile si des
    # cables ou un joint tournant limitent la rotation continue (Benchy ~32 tours, dino ~621,
    # 2026-10-06). Avec `strict`, un depassement refuse l'ecriture, comme une pointe sous le plateau.
    max_y_turns: Optional[float] = None
    platform_margin: float = 1.0         # marge tete / plateau (mm)
    singularity_radius: float = 1e-6     # R sous lequel Y est herite du point precedent
    decimals: int = 4
    # Vitesses maximales des axes du firmware (config.g actif, M203 / 60) pour l'estimation du
    # temps `print_time_fw_s` : (X mm/s, Y deg/s, Z mm/s, B deg/s). Z : 8000 mm/min = 133 mm/s
    # RECOMMANDE (l'actuel, 400 mm/min = 6,7 mm/s, multiplie le temps du Benchy par ~6 : le
    # pivot de B couple Z a B). Mettre (1666.7, 600.0, 6.67, 1666.7) pour l'etat actuel de config.g.
    # None = pas d'estimation. Vitesses seules, accelerations ignorees : borne BASSE.
    fw_axis_speeds: Optional[tuple] = (1666.7, 600.0, 133.3, 1666.7)
    # Le firmware interpole X, Z et B LINEAIREMENT. Le pivot etant a
    # `pivot_to_tip` de la pointe, un changement de B entre deux points proches
    # fait sortir la POINTE de sa trajectoire PENDANT le mouvement (mesure du
    # 2026-10-01 : 1 953 mouvements extrudes du Benchy a plus de 1 mm, jusqu'a
    # 26 mm ; visibles dans la previsualisation de l'utilisateur). Chaque
    # mouvement est donc decoupe jusqu'a ce que la pointe reste a moins de cette
    # distance (mm) de la droite visee. None = desactive (historique).
    max_tip_deviation: Optional[float] = 0.05


@dataclass
class GcodeReport:
    lines: int = 0
    moves: int = 0
    travels: int = 0
    subdivided_moves: int = 0             # mouvements decoupes (max_tip_deviation)
    max_tip_deviation_before: float = 0.0  # pire ecart de la pointe AVANT decoupage (mm)
    retractions: int = 0
    clearance_axes: dict = field(default_factory=lambda: {"Z": 0, "R": 0, "none": 0})
    filament_mm: float = 0.0
    print_time_nominal_s: float = 0.0
    print_time_s: float = 0.0            # avec le ralentissement aux limites de B
    c_mismatch_deg: float = 0.0          # ecart max |Y - (C - 90)|, controle de repere
    travel_collisions: int = 0           # degagements signales par `travel_check`
    # garde-fous (dans le modele de tete d'Apex)
    tip_z_min: float = float("inf")      # pointe la plus basse emise (mm)
    head_z_min: float = float("inf")     # point le plus bas de la vraie tete (si head_zlow)
    tip_below_plate: int = 0             # mouvements dont la pointe passe sous z = 0
    head_below_margin: int = 0           # mouvements dont la tete passe sous platform_margin
    fw_x_range: tuple = (0.0, 0.0)
    fw_y_range: tuple = (0.0, 0.0)       # Y deroule : etendue en degres (tours cumules)
    fw_z_range: tuple = (0.0, 0.0)
    written: bool = False                # faux si `strict` a refuse l'ecriture
    y_turns: float = 0.0                 # etendue de Y deroule, en tours (max - min) / 360
    print_time_fw_s: float = 0.0         # avec les vitesses max de TOUS les axes (fw_axis_speeds)


def toolpath_arrays(toolpath) -> dict:
    """ToolpathMachine (Apex) -> dictionnaire d'entree de `export_gcode`."""
    pts = toolpath.points
    return {"position": np.array([p.position for p in pts], dtype=float),
            "axis_angles": np.array([p.axis_angles for p in pts], dtype=float),
            "extrusion": np.array([p.extrusion for p in pts], dtype=float),
            "is_skin": np.array([bool(getattr(p, "is_skin", False)) for p in pts], dtype=bool)}


def head_body_direction_rz(b_deg):
    """Direction du corps de la tete (pointe -> pivot) dans le plan (R, Z)."""
    b = np.radians(b_deg)
    return np.sin(b), np.cos(b)


def _z_offset(cz, z_reference):
    if z_reference == "tip":
        return cz - 1.0
    if z_reference == "pivot":
        return cz
    raise ValueError(f"z_reference inconnu : {z_reference!r} (attendu 'tip' ou 'pivot')")


def _encrypted_header(params):
    """Bloc d'en-tete chiffre Apex (`sliced` + `meta`, R-010) ; `[]` si le chiffrement est indisponible."""
    from config.schema import APEX_CONTRACT_VERSION
    from file_io.apex_header import header_block_lines
    meta = {"generator": "Apex engine", "contract": APEX_CONTRACT_VERSION,
            "kinematics": {"pivot_to_tip": float(params.pivot_to_tip), "b_sign": float(params.b_sign),
                           "tilt_sign": float(params.tilt_sign), "z_reference": params.z_reference},
            "layer_height": float(params.layer_height), "line_width": float(params.line_width),
            "print_speed": float(params.print_speed), "travel_speed": float(params.travel_speed),
            "nozzle_temp": float(params.nozzle_temp)}
    return header_block_lines(meta)


KINEMATICS_HEADER_KEY = "apex_kinematics"


def kinematics_header_line(params):
    """Ligne d'en-tete `; apex_kinematics pivot_to_tip=46 b_sign=1 tilt_sign=1 z_reference=tip` (R-007).

    Dit avec quelle cinematique pivot/pointe le fichier a ete ecrit : l'apercu et le noyau d'edition de
    l'interface la relisent au lieu de la supposer. Quand le bloc chiffre (R-010) sera ecrit, ces valeurs
    rejoindront `meta` ; la ligne en clair reste, comme `; apex_sliced = 1`."""
    return (f"; {KINEMATICS_HEADER_KEY} pivot_to_tip={float(params.pivot_to_tip):g} "
            f"b_sign={float(params.b_sign):g} tilt_sign={float(params.tilt_sign):g} "
            f"z_reference={params.z_reference}")


def parse_kinematics_header(lines):
    """Relit `kinematics_header_line` ; `None` si le fichier n'en porte pas (G-code d'une ancienne version)."""
    prefix = f"; {KINEMATICS_HEADER_KEY} "
    for line in lines:
        line = line.strip()
        if line.startswith(prefix):
            kv = dict(item.split("=", 1) for item in line[len(prefix):].split())
            return {"pivot_to_tip": float(kv["pivot_to_tip"]), "b_sign": float(kv["b_sign"]),
                    "tilt_sign": float(kv["tilt_sign"]), "z_reference": kv["z_reference"]}
        if line[:3] in ("G00", "G01"):
            break   # l'en-tete est fini au premier mouvement
    return None


def machine_to_firmware(position, b_deg, pivot_to_tip, previous_y=None, singularity_radius=1e-6,
                        tilt_sign=1.0, z_reference="tip", b_sign=1.0):
    """Pointe (x, y, z) + B d'Apex -> (X, Y, Z, B) firmware (B emis = b_sign * B). Y deroule par rapport a
    `previous_y` (jamais de saut de 360 deg) ; herite sous `singularity_radius`."""
    x, y, z = (float(v) for v in position)
    r = float(np.hypot(x, y))
    if r < singularity_radius and previous_y is not None:
        fy = previous_y
    else:
        fy = float(np.degrees(np.arctan2(y, x)))
        if previous_y is not None:
            fy = previous_y + ((fy - previous_y + 180.0) % 360.0 - 180.0)
    sr, cz = head_body_direction_rz(b_deg)
    return (r + tilt_sign * pivot_to_tip * sr, fy, z + pivot_to_tip * _z_offset(cz, z_reference), float(b_sign * b_deg))


def firmware_to_machine(fw_x, fw_y, fw_z, fw_b, pivot_to_tip, tilt_sign=1.0, z_reference="tip", b_sign=1.0):
    """Inverse de `machine_to_firmware` (aperçu et controle aller-retour)."""
    sr, cz = head_body_direction_rz(fw_b * b_sign)
    r = fw_x - tilt_sign * pivot_to_tip * sr
    z = fw_z - pivot_to_tip * _z_offset(cz, z_reference)
    a = np.radians(fw_y)
    return np.array([r * np.cos(a), r * np.sin(a), z])


def head_axis_vectors(position, b_deg):
    """Direction 3D du corps de la tete (pointe -> pivot) en chaque point :
    sin B selon +R (radiale du point), cos B selon +Z. Au centre (R ~ 0), la
    radiale est indefinie ; sin B y est negligeable en pratique, +X est pris."""
    pos = np.atleast_2d(np.asarray(position, dtype=float))
    b = np.radians(np.atleast_1d(np.asarray(b_deg, dtype=float)))
    r = np.hypot(pos[:, 0], pos[:, 1])
    radial = np.zeros_like(pos)
    ok = r > 1e-9
    radial[ok, 0] = pos[ok, 0] / r[ok]
    radial[ok, 1] = pos[ok, 1] / r[ok]
    radial[~ok, 0] = 1.0
    return radial * np.sin(b)[:, None] + np.array([0.0, 0.0, 1.0]) * np.cos(b)[:, None]


def feed_plan(points, angles, max_velocity_deg, max_accel_deg, coupled_axes, speed):
    """Temps d'une suite de points quand le firmware COORDONNE les axes : les
    angles restent exacts, le mouvement entier ralentit la ou un axe
    d'orientation (non couple a la position) atteint sa vitesse ou son
    acceleration maximale. Mouvement 1D de l'angle par segment, vitesse
    plafonnee par passes avant/arriere (v^2 <= v0^2 + 2 a |dtheta|, remise a
    zero a chaque changement de sens). Rend (temps nominal, temps reel,
    facteur de vitesse par segment, <= 1)."""
    seg = np.linalg.norm(np.diff(points, axis=0), axis=1)
    t_nom = seg / max(speed, 1e-9)
    t_real = t_nom.copy()
    if len(seg) and max_velocity_deg > 0 and max_accel_deg > 0:
        for ax in range(angles.shape[1]):
            if ax in coupled_axes:
                continue
            d = np.degrees(np.diff(np.unwrap(np.radians(angles[:, ax]))))
            mag = np.abs(d)
            from core import _cpp
            if _cpp.enabled() and hasattr(_cpp.module(), "feed_plan_caps"):
                # pas 24 : memes passes avant/arriere en natif, au bit
                cap = _cpp.module().feed_plan_caps(np.ascontiguousarray(d), float(max_velocity_deg),
                                                   float(max_accel_deg))
                t_real = np.maximum(t_real, mag / np.maximum(cap, 1e-9))
                continue
            cap = np.full(len(d), float(max_velocity_deg))
            prev = 0.0
            for i in range(len(d)):
                if i and np.sign(d[i]) != np.sign(d[i - 1]):
                    prev = 0.0
                cap[i] = min(cap[i], np.sqrt(prev ** 2 + 2 * max_accel_deg * mag[i]))
                prev = cap[i]
            nxt = 0.0
            for i in range(len(d) - 1, -1, -1):
                if i < len(d) - 1 and np.sign(d[i]) != np.sign(d[i + 1]):
                    nxt = 0.0
                cap[i] = min(cap[i], np.sqrt(nxt ** 2 + 2 * max_accel_deg * mag[i]))
                nxt = cap[i]
            t_real = np.maximum(t_real, mag / np.maximum(cap, 1e-9))
    factor = np.where(t_real > 0, t_nom / np.maximum(t_real, 1e-12), 1.0)
    return float(t_nom.sum()), float(t_real.sum()), factor


def firmware_time_estimate(fw_moves, params: "GcodeParams") -> float:
    """Temps (s) d'une suite de mouvements firmware (X, Y, Z, B, extrude?) : chaque mouvement dure
    le plus long de sa duree nominale a la vitesse de pointe et de |delta axe| / vitesse max de
    l'axe (le firmware coordonne les axes et ralentit l'ensemble). Accelerations ignorees."""
    if params.fw_axis_speeds is None or len(fw_moves) < 2:
        return 0.0
    m = np.array([mv[:4] for mv in fw_moves], dtype=float)
    ext = np.array([mv[4] for mv in fw_moves], dtype=bool)
    a, bb = np.radians(m[:, 1]), np.radians(m[:, 3] * params.b_sign)
    r = m[:, 0] - params.tilt_sign * params.pivot_to_tip * np.sin(bb)
    z = m[:, 2] - params.pivot_to_tip * (np.cos(bb) - 1.0 if params.z_reference == "tip" else np.cos(bb))
    tip = np.stack([r * np.cos(a), r * np.sin(a), z], axis=1)
    length = np.linalg.norm(np.diff(tip, axis=0), axis=1)
    nominal = np.where(ext[1:], length / params.print_speed, length / params.travel_speed)
    d = np.abs(np.diff(m, axis=0))
    caps = np.asarray(params.fw_axis_speeds, dtype=float)
    return float(np.maximum(nominal, (d / np.maximum(caps, 1e-9)).max(axis=1)).sum())


def clearance_move(position, b_deg, params, head_zlow: Optional[Callable[[float], float]] = None,
                   blocked: Optional[Callable[[np.ndarray], bool]] = None):
    """Degagement sur UN seul axe (Z ou R) selon l'orientation de la tete.
    Rend (axe, point degage) ; axe 'none' si aucun des deux n'est sur.

    `blocked(point) -> bool` (optionnel) signale un point de degagement en contact
    avec la matiere deja deposee : le degagement le plus court non bloque est alors
    retenu (les deux axes au meme pas avant de doubler le pas, jusqu'a
    `travel_clearance_max_factor`). Si tous sont bloques, le premier candidat sur
    pour le plateau est conserve : jamais pire que sans controle."""
    sr, cz = head_body_direction_rz(b_deg)
    x, y, z = (float(v) for v in position)
    r = float(np.hypot(x, y))
    radial = np.array([x, y, 0.0]) / r if r > 1e-9 else np.array([1.0, 0.0, 0.0])
    order = ["Z", "R"] if abs(cz) >= abs(sr) else ["R", "Z"]
    factors = [1.0]
    while blocked is not None and factors[-1] * 2.0 <= params.travel_clearance_max_factor + 1e-9:
        factors.append(factors[-1] * 2.0)
    fallback = None
    for factor in factors:
        for axis in order:
            if axis == "Z":
                step = np.array([0.0, 0.0, params.travel_clearance * factor * (1.0 if cz >= 0 else -1.0)])
            else:
                if r < 1e-9 and sr < 0:
                    continue
                step = radial * params.travel_clearance * factor * (1.0 if sr >= 0 else -1.0)
            cand = np.asarray(position, dtype=float) + step
            lowest = cand[2] + (head_zlow(b_deg) if head_zlow is not None else 0.0)
            if not (lowest >= params.platform_margin or cand[2] >= position[2]):
                continue
            if blocked is None or not blocked(cand):
                return axis, cand
            if fallback is None:
                fallback = (axis, cand)
    if fallback is not None:
        return fallback
    return "none", np.asarray(position, dtype=float)


def material_travel_check(toolpath: dict, machine, params: "GcodeParams", densify_step: float = 0.25,
                          axis_samples: int = 16) -> Callable[[int, np.ndarray, float], bool]:
    """`travel_check(i, point, B)` de la chaine de recherche (export_gcode.py,
    2026-10-06), devenu fonction du moteur : un point de degagement est en contact
    si de la matiere deposee AVANT le point de trajectoire i (cordons densifies a
    `densify_step`, positions de `toolpath` non decalees) est a <= rayon de buse
    d'un des `axis_samples` points de l'axe de la tete, de la pointe au pivot
    (tete orientee par C = atan2(x, -y) et B). Porte ses donnees en attribut
    `.native` : l'export natif fait le meme controle sans rappel Python. `.hits`
    garde les contacts vus par le chemin Python."""
    from scipy.spatial import cKDTree
    from pipeline import step9_collision_sim as s9

    pos = np.asarray(toolpath["position"], dtype=float)
    ext = np.asarray(toolpath["extrusion"], dtype=float)
    # Densification vectorisee (2026-10-13 ; 9,7 s -> moins d'une seconde sur le Benchy) :
    # memes operations que la boucle d'origine (k/m puis a + t*seg), points identiques au bit.
    seg_ids = np.flatnonzero(ext[1:] > 0) + 1
    if len(seg_ids):
        seg = pos[seg_ids] - pos[seg_ids - 1]
        m = np.maximum(1, np.ceil(np.linalg.norm(seg, axis=1) / densify_step).astype(np.int64))
        order = np.repeat(seg_ids, m)
        first = np.cumsum(m) - m
        k = np.arange(len(order)) - np.repeat(first, m) + 1
        t = (k / np.repeat(m, m))[:, None]
        dense = pos[order - 1] + t * seg[np.repeat(np.arange(len(seg_ids)), m)]
    else:
        dense, order = np.zeros((0, 3)), np.zeros(0, dtype=np.int64)
    tree = cKDTree(dense) if len(dense) else None
    radius = machine.nozzle_diameter / 2
    axis_t = np.linspace(0.0, params.pivot_to_tip, axis_samples)
    hits = []

    def travel_check(i, point, b):
        if tree is None:
            return False
        c = float(np.degrees(np.arctan2(point[0], -point[1])))
        axis_pts = s9._positioned_toolhead_vertices(np.column_stack([0 * axis_t, 0 * axis_t, axis_t]), point,
                                                    [c, b], machine.rotation_axes)
        for lst in tree.query_ball_point(axis_pts, radius):
            if lst and (order[lst] < i).any():
                hits.append((i, np.asarray(point).tolist(), b))
                return True
        return False

    axes = [np.asarray(a.axis, dtype=float) / np.linalg.norm(a.axis) for a in machine.rotation_axes]
    travel_check.hits = hits
    if len(axes) == 2:
        travel_check.native = {"points": np.ascontiguousarray(dense), "owners": order.astype(np.int64),
                               "axis_t": axis_t, "radius": float(radius), "axis_c": axes[0], "axis_b": axes[1],
                               "max_factor": float(params.travel_clearance_max_factor)}
    return travel_check


def head_zlow_from_table(grid: np.ndarray, values: np.ndarray) -> Callable[[float], float]:
    """`head_zlow` tabule (valeur au point de grille de B le plus proche,
    premier en cas d'egalite) ; porte sa table en attribut `.table`, ce qui
    permet a l'export natif (pas 24) de l'utiliser sans rappel Python."""
    grid = np.asarray(grid, dtype=float)
    values = np.asarray(values, dtype=float)

    def head_zlow(b):
        return float(values[np.abs(grid - b).argmin()])

    head_zlow.table = (grid, values)
    return head_zlow


def export_gcode(toolpath: dict, machine, params: GcodeParams, path: Optional[str] = None,
                 head_zlow: Optional[Callable[[float], float]] = None,
                 travel_check: Optional[Callable[[int, np.ndarray, float], bool]] = None):
    """Ecrit le G-code (si `path`) et rend (lignes, GcodeReport).

    `head_zlow(B)` : hauteur du point le plus bas de la vraie tete sous la
    pointe pour cet angle (negative), pour refuser un degagement vers le
    plateau. `travel_check(i, point, B) -> bool` : controle de collision d'un
    point de degagement (vrai = contact), compte dans le rapport."""
    pos = np.asarray(toolpath["position"], dtype=float)
    ang = np.asarray(toolpath["axis_angles"], dtype=float)
    ext = np.asarray(toolpath["extrusion"], dtype=float)
    skin = np.asarray(toolpath.get("is_skin", np.zeros(len(pos), bool)), dtype=bool)
    n = len(pos)
    if params.tip_at_bead_top and n:
        half = np.where(skin, params.skin_thickness, params.layer_height) * 0.5
        pos = pos + head_axis_vectors(pos, ang[:, 1]) * half[:, None]
    rep = GcodeReport()
    fil_area = np.pi * (params.filament_diameter / 2) ** 2
    L, dec = params.pivot_to_tip, params.decimals
    coupled = list(getattr(machine, "position_coupled_axes", None) or [])
    vmax = float(np.degrees(getattr(machine, "max_angular_velocity", 0.0) or 0.0))
    amax = float(np.degrees(getattr(machine, "max_angular_acceleration", 0.0) or 0.0))
    _, _, slow = feed_plan(pos, ang, vmax, amax, coupled, params.print_speed) if n > 1 else (0, 0, np.ones(0))

    # En-tete natif (pas 24, suite) quand le format de repr(float) est le meme ;
    # sinon (valeur hors [1e-4, 1e16)) le texte Python.
    from core import _cpp
    native_header = None
    if _cpp.enabled() and hasattr(_cpp.module(), "gcode_header"):
        native_header = _cpp.module().gcode_header(
            float(L), float(params.retract_length), float(params.retract_speed),
            float(params.nozzle_temp))
    if native_header is not None:
        out = native_header.split("\n")
    else:
        out = ["; Apex -- G-code machine R C Z B (file_io/gcode_export.py)",
               f"; pivot->pointe {L} mm compense ; retraction {params.retract_length} mm a {params.retract_speed} mm/s",
               "G21", "G90", "M82", "G92 E0", f"M109 S{params.nozzle_temp:g}",
               # Marqueur « tranche par Apex » lu par l'interface (onglet Edition : modification en
               # place, sinon EXPORT COPY seulement ; cf. slicer_ux/CDC_INTEGRATION_NOEUDS.md, D4).
               # En clair tant que le moteur n'ecrit pas le bloc d'en-tete chiffre Apex
               # (slicer_ux/docs/APEX_HEADER_FORMAT.md), ou il rejoindra les metadonnees.
               "; apex_sliced = 1"]
    out.append(kinematics_header_line(params))   # R-007 : la cinematique d'ecriture voyage avec le fichier
    out += _encrypted_header(params)              # R-010 : meme information dans le bloc chiffre (sans le paquet : rien)
    state = {"fw": None, "e": 0.0, "retracted": False}
    guard = {"x": [], "y": [], "z": []}

    def fmt(v):
        s = f"{v:.{dec}f}".rstrip("0").rstrip(".")
        return "0" if s in ("-0", "") else s

    # Pas 24 du portage C++ : corps du fichier en natif, texte identique ligne
    # a ligne (tests/test_cpp_equivalence.py), decoupage `max_tip_deviation`
    # compris depuis le 2026-10-07, controle des degagements depuis le 2026-10-13
    # (`material_travel_check`, porte en donnees). Hors de portee du natif : un
    # `travel_check` quelconque (rappel Python sans `.native`) ou un `head_zlow` sans table.
    table = getattr(head_zlow, "table", None) if head_zlow is not None else (np.zeros(0), np.zeros(0))
    travel_native = getattr(travel_check, "native", None) if travel_check is not None else None
    if (_cpp.enabled() and hasattr(_cpp.module(), "export_gcode_body")
            and (travel_check is None
                 or (travel_native is not None and getattr(_cpp.module(), "EXPORT_TRAVEL_CHECK", False)))
            and table is not None and n):
        body, st = _cpp.module().export_gcode_body(
            np.ascontiguousarray(pos), np.ascontiguousarray(ang[:, :2]), np.ascontiguousarray(ext),
            np.ascontiguousarray(skin, dtype=np.uint8), np.ascontiguousarray(slow, dtype=float),
            np.ascontiguousarray(table[0], dtype=float), np.ascontiguousarray(table[1], dtype=float),
            {"layer_height": params.layer_height, "skin_thickness": params.skin_thickness,
             "line_width": params.line_width, "flow": params.flow, "print_speed": params.print_speed,
             "travel_speed": params.travel_speed, "rotation_feed": params.rotation_feed,
             "max_feed": params.max_feed, "retract_length": params.retract_length,
             "retract_speed": params.retract_speed, "retract_min_travel": params.retract_min_travel,
             "travel_clearance": params.travel_clearance, "pivot_to_tip": L,
             "platform_margin": params.platform_margin, "singularity_radius": params.singularity_radius,
             "tilt_sign": params.tilt_sign, "b_sign": params.b_sign, "z_reference": params.z_reference,
             "decimals": dec,
             "fil_area": fil_area, "max_tip_deviation": params.max_tip_deviation},
            *(() if travel_native is None else (travel_native,)))
        _z_offset(1.0, params.z_reference)   # meme erreur que le chemin Python si z_reference est inconnu
        if "subdivided_moves" not in st:
            # Module natif PERIME (anterieur au decoupage `max_tip_deviation`, pas 24 etendu) : il
            # ecrirait des mouvements non decoupes, donc une pointe qui quitte sa trajectoire. Son
            # resultat est ecarte et le chemin Python, identique, le remplace (jamais un fichier faux).
            n_loop = n
        else:
            rep.travel_collisions = int(st.get("travel_collisions", 0))
            out += body.split(chr(10))
            rep.moves, rep.travels, rep.retractions = st["moves"], st["travels"], st["retractions"]
            rep.subdivided_moves = st["subdivided_moves"]
            rep.max_tip_deviation_before = st["max_tip_deviation_before"]
            rep.clearance_axes = {"Z": st["clearance_Z"], "R": st["clearance_R"], "none": st["clearance_none"]}
            rep.filament_mm, rep.c_mismatch_deg = st["filament_mm"], st["c_mismatch_deg"]
            rep.tip_z_min, rep.head_z_min = st["tip_z_min"], st["head_z_min"]
            rep.tip_below_plate, rep.head_below_margin = st["tip_below_plate"], st["head_below_margin"]
            if st["any"]:
                guard["x"], guard["y"], guard["z"] = list(st["x"]), list(st["y"]), list(st["z"])
            n_loop = 0
    else:
        n_loop = n

    def move(point, b, tip_speed, e_add=0.0, comment=""):
        """Mouvement vers (point, b), DECOUPE si besoin pour que la pointe suive
        la droite visee a `max_tip_deviation` pres (voir GcodeParams)."""
        tol = params.max_tip_deviation
        prev = state["fw"]
        if tol is None or prev is None:
            emit(point, b, tip_speed, e_add, comment)
            return
        a_tip, a_b = state["tip"], float(prev[3]) * params.b_sign
        pieces = []

        def deviation(p0, b0, fw0, p1, b1):
            fw1 = machine_to_firmware(p1, b1, L, fw0[1], params.singularity_radius,
                                      params.tilt_sign, params.z_reference, params.b_sign)
            mid_fw = [(u + v) / 2.0 for u, v in zip(fw0, fw1)]
            mid_tip = firmware_to_machine(*mid_fw, L, params.tilt_sign, params.z_reference, params.b_sign)
            return float(np.linalg.norm(mid_tip - (np.asarray(p0) + np.asarray(p1)) / 2.0)), fw1

        def split(p0, b0, fw0, p1, b1, depth):
            dev, fw1 = deviation(p0, b0, fw0, p1, b1)
            if depth == 0:
                rep.max_tip_deviation_before = max(rep.max_tip_deviation_before, dev)
            if dev <= tol or depth >= 10:
                pieces.append((np.asarray(p1, dtype=float), b1))
                return fw1
            pm = (np.asarray(p0, dtype=float) + np.asarray(p1, dtype=float)) / 2.0
            bm = (b0 + b1) / 2.0
            fwm = split(p0, b0, fw0, pm, bm, depth + 1)
            return split(pm, bm, fwm, p1, b1, depth + 1)

        split(a_tip, a_b, prev, np.asarray(point, dtype=float), float(b), 0)
        if len(pieces) > 1:
            rep.subdivided_moves += 1
        total = float(np.linalg.norm(np.asarray(point, dtype=float) - a_tip))
        start = a_tip
        for k, (p, bb) in enumerate(pieces):
            share = float(np.linalg.norm(p - start)) / total if total > 1e-12 else 1.0 / len(pieces)
            emit(p, bb, tip_speed, e_add * share, comment if k == len(pieces) - 1 else "")
            start = p

    def emit(point, b, tip_speed, e_add=0.0, comment=""):
        prev = state["fw"]
        fw = machine_to_firmware(point, b, L, None if prev is None else prev[1], params.singularity_radius,
                                 params.tilt_sign, params.z_reference, params.b_sign)
        tz = float(point[2])
        rep.tip_z_min = min(rep.tip_z_min, tz)
        rep.tip_below_plate += int(tz < 0.0)
        if head_zlow is not None:
            hz = tz + head_zlow(b)
            rep.head_z_min = min(rep.head_z_min, hz)
            rep.head_below_margin += int(hz < params.platform_margin)
        guard["x"].append(fw[0])
        guard["y"].append(fw[1])
        guard["z"].append(fw[2])
        tip_d = 0.0 if prev is None else float(np.linalg.norm(point - state["tip"]))
        fw_d = 0.0 if prev is None else float(np.linalg.norm(np.subtract(fw, prev)))
        if tip_d > 1e-9:
            feed = tip_speed * 60.0 * fw_d / tip_d
        else:
            feed = params.rotation_feed
        feed = min(max(feed, 1.0), params.max_feed)
        state["e"] += e_add
        words = [f"X{fmt(fw[0])}", f"Y{fmt(fw[1])}", f"Z{fmt(fw[2])}", f"B{fmt(fw[3])}"]
        if e_add:
            words.append(f"E{fmt(state['e'])}")
        out.append("G00 " + " ".join(words) + f" F{fmt(feed)}" + (f" ; {comment}" if comment else ""))
        state["fw"], state["tip"] = fw, np.asarray(point, dtype=float)
        rep.moves += 1

    def retract(sign):
        state["e"] += sign * params.retract_length
        out.append(f"G01 E{fmt(state['e'])} F{fmt(params.retract_speed * 60.0)}")
        state["retracted"] = sign < 0

    last_kind = None
    for i in range(n_loop):
        _progress.report(i, n_loop, span=(0.4, 1.0))   # R-017 : points ecrits
        kind = "lissage" if skin[i] else "corps"
        b = float(ang[i, 1])
        if i > 0:
            y_fw = machine_to_firmware(pos[i], b, L, None, params.singularity_radius,
                                       params.tilt_sign, params.z_reference, params.b_sign)[1]
            if np.hypot(*pos[i, :2]) > 1e-3:
                rep.c_mismatch_deg = max(rep.c_mismatch_deg,
                                         abs((y_fw - (ang[i, 0] - 90.0) + 180.0) % 360.0 - 180.0))
        if kind != last_kind:
            out.append(f"; --- {kind} ---")
            last_kind = kind
        if i == 0 or ext[i] <= 0:
            # trajet a vide vers pos[i]
            rep.travels += 1
            if i == 0:
                axis, up = clearance_move(pos[0], b, params, head_zlow)
                move(up, b, params.travel_speed, comment="approche")
                rep.clearance_axes[axis] += 1
                move(pos[0], b, params.travel_speed)
                continue
            length = float(np.linalg.norm(pos[i] - pos[i - 1]))
            if length >= params.retract_min_travel:
                retract(-1)
                rep.retractions += 1
                b0 = float(ang[i - 1, 1])
                blk0 = (lambda p, i=i, b0=b0: travel_check(i, p, b0)) if travel_check is not None else None
                blk1 = (lambda p, i=i, b=b: travel_check(i, p, b)) if travel_check is not None else None
                axis0, up0 = clearance_move(pos[i - 1], b0, params, head_zlow, blocked=blk0)
                axis1, up1 = clearance_move(pos[i], b, params, head_zlow, blocked=blk1)
                rep.clearance_axes[axis0] += 1
                if travel_check is not None:
                    rep.travel_collisions += int(travel_check(i, up0, b0)) + int(travel_check(i, up1, b))
                move(up0, b0, params.travel_speed, comment=f"degagement {axis0}")
                move(up1, b, params.travel_speed)
                move(pos[i], b, params.travel_speed)
                retract(+1)
            else:
                move(pos[i], b, params.travel_speed)
            continue
        seg = float(np.linalg.norm(pos[i] - pos[i - 1]))
        h = params.skin_thickness if skin[i] else params.layer_height
        e_add = seg * h * params.line_width * params.flow / fil_area
        rep.filament_mm += e_add
        speed = params.print_speed * (float(slow[i - 1]) if len(slow) >= i else 1.0)
        move(pos[i], b, speed, e_add=e_add)
    if n_loop:
        retract(-1)
        axis, up = clearance_move(pos[-1], float(ang[-1, 1]), params, head_zlow)
        move(up, float(ang[-1, 1]), params.travel_speed, comment=f"degagement final {axis}")
    if _cpp.enabled() and hasattr(_cpp.module(), "gcode_footer"):
        out += _cpp.module().gcode_footer().split("\n")
    else:
        out += ["M104 S0", "M84"]
    printing = ext > 0
    if n > 1:
        seg = np.linalg.norm(np.diff(pos, axis=0), axis=1)
        pr = printing[1:]
        rep.print_time_nominal_s = float(seg[pr].sum() / params.print_speed + seg[~pr].sum() / params.travel_speed)
        rep.print_time_s = float((seg[pr] / (params.print_speed * np.maximum(slow[pr], 1e-9))).sum()
                                 + seg[~pr].sum() / params.travel_speed)
    rep.lines = len(out)
    if params.fw_axis_speeds is not None and n > 1:
        rep.print_time_fw_s = firmware_time_estimate(parse_gcode(out), params)
    if guard["x"]:
        rep.fw_x_range = (float(min(guard["x"])), float(max(guard["x"])))
        rep.fw_y_range = (float(min(guard["y"])), float(max(guard["y"])))
        rep.y_turns = (rep.fw_y_range[1] - rep.fw_y_range[0]) / 360.0
        rep.fw_z_range = (float(min(guard["z"])), float(max(guard["z"])))
    if params.strict and (rep.tip_below_plate
                          or (params.max_y_turns is not None and rep.y_turns > params.max_y_turns)):
        return out, rep
    if path:
        rep.written = True
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(out) + "\n")
    return out, rep


def parse_gcode(lines):
    """Relit les mouvements G00 d'un fichier produit par `export_gcode` :
    liste de (X, Y, Z, B, extrude?) firmware, dans l'ordre."""
    moves, e_prev = [], 0.0
    cur = {}
    for line in lines:
        code = line.split(";", 1)[0].strip()
        if not code:
            continue
        words = code.split()
        if words[0] not in ("G00", "G0", "G01", "G1"):
            continue
        vals = {w[0]: float(w[1:]) for w in words[1:]}
        e = vals.get("E")
        extrude = e is not None and e > e_prev + 1e-12 and any(k in vals for k in "XYZB")
        if e is not None:
            e_prev = e
        if any(k in vals for k in "XYZB"):
            cur.update({k: vals[k] for k in "XYZB" if k in vals})
            moves.append((cur["X"], cur["Y"], cur["Z"], cur["B"], extrude))
    return moves
