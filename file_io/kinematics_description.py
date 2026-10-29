"""Description de la cinématique d'une machine, écrite au début du G-code (`apex-kinematics/1`).

Objectif (demande de l'utilisateur, 2026-10-29) : un format que n'importe quel lecteur du G-code
(aperçu, interface, autre logiciel) peut relire pour connaître EXACTEMENT l'architecture de la
machine et recalculer la position de la pointe depuis les valeurs des axes XYZABC.

Principe : une CHAÎNE SÉRIE de liaisons, comme URDF ou la table d'axes d'un post-processeur
5 axes, de deux côtés d'un repère machine fixe : le côté PIÈCE (porte-plateau : ce qui porte la
pièce) et le côté OUTIL (ce qui porte la buse). Chaque axe, un par lettre de XYZABC (les six sont
TOUJOURS listés, un axe inutilisé a `type: "unused"`), dit :

  letter        lettre du G-code
  type          "linear" | "rotary" | "unused"
  side          "tool" | "bed" (côté de la chaîne où il se trouve)
  axis          vecteur directeur unitaire, dans le repère machine À L'ORIGINE de TOUS les axes
  origin        un point de l'axe (mm) pour un axe rotatif ; absent pour un axe linéaire
  range         [min, max] en mm ou degrés ; null si illimité (rotation déroulée)
  positive      sens positif, en une phrase (jamais « le sens du moteur »)

  Sens des rotations : RÈGLE UNIQUE, main droite autour de `axis` (`conventions.rotation`) ; `positive`
  n'est qu'un commentaire lisible.

`part_origin` : le point de la pièce qui coïncide avec le pivot quand tous les axes sont à 0 (sans lui,
deux lecteurs peuvent placer la pièce différemment sur le plateau).

L'ordre des axes dans chaque côté est l'ordre d'application depuis le repère machine : la pointe
est obtenue en appliquant les liaisons du côté outil du repère machine vers la buse, puis
`tool_offset`. `conventions` porte ce qui n'est pas géométrique (compensation du pivot, signes).

Limite assumée : une chaîne série décrit les machines cartésiennes, polaires, à plateau ou tête
tournants, 4 à 6 axes. Un couplage moteurs → axes (CoreXY, delta) se déclare dans `coupling`
(`null` ici : les axes du G-code sont déjà les axes cartésiens de la machine).
"""
from __future__ import annotations

import json
from typing import Optional

SCHEMA = "apex-kinematics/1"
ROTATION_RULE = "main droite autour de axis, appliquee dans l'ordre de chain (de la machine vers la piece ou l'outil)"
JSON_KEY = "apex_kinematics_json"


def _unused(letter: str) -> dict:
    return {"letter": letter, "type": "unused", "note": "toujours 0"}


def cartesian(travel: Optional[dict] = None) -> dict:
    """Imprimante cartésienne : X, Y, Z linéaires ; A, B, C inutilisés."""
    travel = travel or {}
    lin = lambda l, a: {"letter": l, "type": "linear", "side": "tool" if l != "Y" else "bed",   # noqa: E731
                        "axis": a, "range": travel.get(l), "positive": "coordonnée croissante"}
    axes = [lin("X", [1, 0, 0]), lin("Y", [0, 1, 0]), lin("Z", [0, 0, 1]),
            _unused("A"), _unused("B"), _unused("C")]
    return {"schema": SCHEMA, "machine": "cartesian", "units": {"linear": "mm", "angular": "deg"},
            "word_order": "XYZABC", "frame": "repère machine droit, Z vers le haut",
            "axes": axes, "chain": {"bed": ["Y"], "tool": ["X", "Z"]},
            "tool_offset": [0, 0, 0], "part_origin": [0, 0, 0], "coupling": None,
            "conventions": {"rotation": ROTATION_RULE}}


def apex_polar(params, machine=None) -> dict:
    """Machine Apex actuelle : plateau tournant Y, tête X (radiale) + Z + pivot B.

    Valeurs lues dans `params` (GcodeParams : pivot→pointe, signes) et, si fourni, dans le
    `MachineProfile` (plage de B, centre du plateau). Repère : celui de la pièce quand Y = 0.
    """
    b_range = [0.0, 185.0]
    origin = [0.0, 0.0, 0.0]
    if machine is not None:
        if len(getattr(machine, "rotation_axes", [])) >= 2:
            lo, hi = machine.rotation_axes[1].angle_range_deg
            b_range = [float(lo), float(hi)]
        if getattr(machine, "platform_origin", None) is not None:
            origin = [float(v) for v in machine.platform_origin]
    L = float(params.pivot_to_tip)
    axes = [
        {"letter": "X", "type": "linear", "side": "tool", "axis": [1, 0, 0], "range": None,
         "positive": "la tête s'éloigne du centre du plateau (+R)"},
        {"letter": "Y", "type": "rotary", "side": "bed", "axis": [0, 0, 1], "origin": origin, "range": None,
         "positive": "Y vaut l'angle polaire (degrés, anti-horaire) du point de la pièce amené sous la tête ; "
                     "rotation déroulée, plusieurs tours possibles"},
        {"letter": "Z", "type": "linear", "side": "tool", "axis": [0, 0, 1], "range": None,
         "positive": "la tête monte"},
        _unused("A"),
        {"letter": "B", "type": "rotary", "side": "tool", "axis": [0, 1, 0],
         "origin": "pivot de la tête (à `pivot_to_tip` au-dessus de la pointe à B = 0)", "range": b_range,
         "positive": "le corps de la tête penche vers +X (+R), la buse vers le centre"},
        _unused("C"),
    ]
    return {"schema": SCHEMA, "machine": "apex-polar-4axis", "units": {"linear": "mm", "angular": "deg"},
            "word_order": "XYZABC",
            "frame": "repère de la pièce quand Y = 0 : +X = direction radiale de la tête, Z vers le haut",
            "axes": axes, "chain": {"bed": ["Y"], "tool": ["X", "Z", "B"]},
            "tool_offset": [0, 0, -L], "part_origin": origin, "coupling": None,
            "conventions": {
                "rotation": ROTATION_RULE,
                "tool_offset_from": "pivot de B",
                "pivot_compensation": "le G-code commande le PIVOT : X = R + L·sin(B), Z = z + L·(cos(B) − 1) ; "
                                      "à B = 0, X et Z sont ceux de la pointe",
                "pivot_to_tip": L, "b_sign": float(params.b_sign), "tilt_sign": float(params.tilt_sign),
                "z_reference": params.z_reference,
                "position_coupled_axes": ["Y"],
                "note": "Y est imposé par la position du point (la tête ne balaie qu'un demi-plan)"}}


def validate(desc: dict) -> None:
    """Refuse une description incohérente (`ValueError`) : six axes XYZABC exactement une fois, chaînes
    qui ne citent que des axes utilisés du bon côté, vecteurs d'axe unitaires, pivot des axes rotatifs."""
    if desc.get("schema") != SCHEMA:
        raise ValueError(f"schéma inconnu : {desc.get('schema')!r}")
    letters = [a["letter"] for a in desc["axes"]]
    if letters != list("XYZABC"):
        raise ValueError(f"les six axes XYZABC doivent figurer une fois, dans cet ordre : {letters}")
    by = {a["letter"]: a for a in desc["axes"]}
    for side in ("bed", "tool"):
        for l in desc["chain"][side]:
            if l not in by or by[l]["type"] == "unused" or by[l]["side"] != side:
                raise ValueError(f"axe {l} cité côté {side} mais déclaré autrement")
    used = {l for side in ("bed", "tool") for l in desc["chain"][side]}
    for a in desc["axes"]:
        if a["type"] == "unused":
            if a["letter"] in used:
                raise ValueError(f"axe inutilisé {a['letter']} présent dans une chaîne")
            continue
        if a["letter"] not in used:
            raise ValueError(f"axe {a['letter']} déclaré mais absent des chaînes")
        v = a["axis"]
        if abs(sum(c * c for c in v) - 1.0) > 1e-9:
            raise ValueError(f"axe {a['letter']} : vecteur non unitaire {v}")
        if a["type"] == "rotary" and "origin" not in a:
            raise ValueError(f"axe rotatif {a['letter']} sans `origin`")
    if len(desc["part_origin"]) != 3 or len(desc["tool_offset"]) != 3:
        raise ValueError("part_origin et tool_offset sont des vecteurs 3D")


def header_lines(desc: dict) -> list[str]:
    """Lignes de commentaire : résumé lisible, puis le JSON complet sur UNE ligne."""
    validate(desc)
    parts = []
    for a in desc["axes"]:
        if a["type"] == "unused":
            parts.append(f"{a['letter']}=unused")
        else:
            rng = "unlimited" if a.get("range") is None else f"{a['range'][0]:g}..{a['range'][1]:g}"
            parts.append(f"{a['letter']}={a['type']}/{a['side']}/{rng}")
    return [f"; kinematics {desc['machine']} ({SCHEMA}) : " + " ".join(parts),
            f"; {JSON_KEY} " + json.dumps(desc, ensure_ascii=True, separators=(",", ":"))]


def parse_header(lines) -> Optional[dict]:
    """Relit le JSON écrit par `header_lines` ; `None` si absent (fichier d'une ancienne version)."""
    prefix = f"; {JSON_KEY} "
    for line in lines:
        line = line.strip()
        if line.startswith(prefix):
            return json.loads(line[len(prefix):])
        if line[:2] in ("G0", "G1"):
            break
    return None
