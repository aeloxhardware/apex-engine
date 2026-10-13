"""
Longueur pivot-pointe L (`GcodeParams.pivot_to_tip`) : valeur courante de la machine.

Le firmware ne fait aucune cinematique (CLAUDE.md limite #61) : L est applique
par le slicer, dans le G-code. Tout G-code produit avant un changement de L
(changement de buse, nouvelle mesure) est donc a REGENERER.

La valeur vit dans `config/pivot_to_tip.json` : une liste d'enregistrements
(date, buse, L, note), le dernier fait foi. Sans fichier : 46 mm (mesure de
l'utilisateur, 2026-10-01).
"""
from __future__ import annotations

import datetime
import json
import os
from dataclasses import asdict, dataclass
from typing import Sequence

DEFAULT_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config", "pivot_to_tip.json")
DEFAULT_PIVOT_TO_TIP = 46.0


@dataclass
class PivotLengthRecord:
    date: str
    nozzle: str
    l_mm: float
    note: str = ""


def load_history(path: str = DEFAULT_PATH) -> list[PivotLengthRecord]:
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as f:
        return [PivotLengthRecord(**{k: row[k] for k in ("date", "nozzle", "l_mm", "note") if k in row}) for row in json.load(f)]


def save_history(history: Sequence[PivotLengthRecord], path: str = DEFAULT_PATH) -> None:
    """Ecriture atomique (fichier temporaire puis remplacement)."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump([asdict(r) for r in history], f, indent=2, ensure_ascii=False)
        f.write("\n")
    os.replace(tmp, path)


def record_pivot_length(l_mm: float, nozzle: str = "", note: str = "", path: str = DEFAULT_PATH) -> PivotLengthRecord:
    record = PivotLengthRecord(date=datetime.date.today().isoformat(), nozzle=nozzle, l_mm=float(l_mm), note=note)
    save_history(load_history(path) + [record], path)
    return record


def current_pivot_to_tip(path: str = DEFAULT_PATH) -> tuple[float, str]:
    """(L, origine) : derniere valeur enregistree, sinon la valeur par defaut."""
    history = load_history(path)
    if history:
        last = history[-1]
        return last.l_mm, f"valeur du {last.date} ({last.nozzle or 'buse non precisee'})"
    return DEFAULT_PIVOT_TO_TIP, "valeur par defaut (aucune valeur enregistree)"
