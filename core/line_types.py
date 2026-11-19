"""Types de ligne d'une trajectoire : ce que dépose chaque cordon (parois, remplissage, pont, peau).

Un type est un petit entier (`uint8`) porté par polyligne dans `Layer.toolpath_kinds` (étapes 6 et 7), puis par
point dans `ToolpathPoint.line_type` (étape 8), enfin écrit dans le G-code par des étiquettes `;TYPE:...` (convention de
Cura, lue par la plupart des visionneuses) pour que l'aperçu puisse colorer les trajectoires, comme le fait
Element Aeolis (`element_aeolis/slicer.py`, `LINE_TYPES`).

Les codes sont partagés avec le noyau natif (`cpp/src/core/step6_layer.cpp`, `cpp/src/core/gcode_export.cpp`) : toute
modification ici se reporte là-bas, et inversement.
"""
from __future__ import annotations

import numpy as np

NONE = 0            # inconnu : trajet à vide, ou trajectoire sans information de type
WALL_OUTER = 1      # contour d'origine de la couche (extérieur et bord de chaque cavité)
WALL_INNER = 2      # parois décalées vers l'intérieur
INFILL = 3          # remplissage
BRIDGE = 4          # remplissage d'une zone en pont (balayé en travers de la portée)
SKIN = 5            # lissage non planaire (peau)

# Étiquette `;TYPE:` de chaque code. NONE sur un mouvement imprimé donne `CUSTOM` (étiquette générique de Cura).
LABELS = {NONE: "CUSTOM", WALL_OUTER: "WALL-OUTER", WALL_INNER: "WALL-INNER", INFILL: "FILL", BRIDGE: "BRIDGE",
          SKIN: "SKIN"}


def inherit_kinds(old_polys: list, old_kinds: list[int] | None, new_polys: list) -> list[int]:
    """Types des polylignes `new_polys`, redécoupées depuis `old_polys` (étape 7, découpe par couleur).

    Sans peinture, rien n'est redécoupé : mêmes polylignes, mêmes types. Avec peinture, chaque polyligne d'origine d'au
    moins deux points est remplacée par des tranches CONTIGUËS (`poly[a:b + 1]`) dont les arêtes la recouvrent une fois
    chacune : on les rattache à leur parente en cumulant les arêtes. Si le découpage ne se laisse pas rattacher, tout
    est rendu `NONE` (jamais un type faux).
    """
    if old_kinds is None or len(old_kinds) != len(old_polys):
        return [NONE] * len(new_polys)
    if len(new_polys) == len(old_polys) and all(len(a) == len(b) for a, b in zip(old_polys, new_polys)):
        return list(old_kinds)
    kinds: list[int] = []
    cursor = 0
    for poly, kind in zip(old_polys, old_kinds):
        edges = len(poly) - 1
        if edges < 1:
            continue                       # ni arête ni tranche (rejetée par l'étape 7)
        remaining = edges
        while remaining > 0 and cursor < len(new_polys):
            used = len(new_polys[cursor]) - 1
            if used < 1 or used > remaining:
                return [NONE] * len(new_polys)
            kinds.append(kind)
            remaining -= used
            cursor += 1
        if remaining != 0:
            return [NONE] * len(new_polys)
    if cursor != len(new_polys):
        return [NONE] * len(new_polys)
    return kinds


def kinds_array(values) -> np.ndarray:
    """Tableau `uint8` des types (forme attendue par le noyau natif)."""
    return np.asarray(values, dtype=np.uint8)
