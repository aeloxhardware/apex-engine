"""Les six axes du G-code d'Apex : X Y Z A B C, toujours tous écrits.

Décision de l'utilisateur (2026-10-29) : pour TOUS les algorithmes de tranchage, le
transcripteur écrit une coordonnée pour chacun des six axes XYZABC sur chaque mouvement, même
si l'un d'eux (ou plusieurs) reste à 0 tout au long du fichier. Un consommateur n'a ainsi
jamais à deviner un axe absent. Le C++ (`cpp/src/core/gcode_export.cpp`) reproduit cet ordre.
"""
from __future__ import annotations

from typing import Callable, Mapping

AXES = "XYZABC"


def axes_words(values: Mapping[str, float], fmt: Callable[[float], str]) -> list[str]:
    """Mots `X.. Y.. Z.. A.. B.. C..` dans cet ordre ; un axe absent de `values` vaut 0."""
    unknown = set(values) - set(AXES)
    if unknown:
        raise ValueError(f"axes inconnus : {sorted(unknown)} (attendus : {AXES})")
    return [f"{a}{fmt(float(values.get(a, 0.0)))}" for a in AXES]
