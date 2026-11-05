"""
config/machine_profiles.py — profils de machine lus depuis des fichiers de données (R-030 point 4).

Aucune constante de machine dans le code : rayon du plateau, tête, cinématique et limites vivent dans
`config/machines/<id>.json`. Ajouter une machine = ajouter un fichier. Le modèle 3D de la tête n'est pas
chargé ici (`main.ensure_toolhead` le résout à partir de `toolhead_id`).
"""

from __future__ import annotations

import json
import os
import re

from config import schema
from core.mesh_types import MachineProfile

MACHINES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "machines")


def default_machine_id() -> str:
    """Profil utilisé quand l'appelant n'en désigne aucun (`config/machines/index.json`)."""
    with open(os.path.join(MACHINES_DIR, "index.json"), encoding="utf-8") as f:
        return str(json.load(f)["default_machine"])


def list_machines() -> list[str]:
    return sorted(n[:-5] for n in os.listdir(MACHINES_DIR) if n.endswith(".json") and n != "index.json")


def load_machine_profile(machine_id: str) -> MachineProfile:
    """Profil `machine_id`, sans modèle de tête. Leve `ValueError` si l'identifiant est invalide,
    `FileNotFoundError` s'il est inconnu."""
    if not re.fullmatch(r"[A-Za-z0-9_-]+", machine_id or "") or machine_id == "index":
        raise ValueError(f"identifiant de machine invalide : {machine_id!r}")
    path = os.path.join(MACHINES_DIR, f"{machine_id}.json")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"profil de machine inconnu : {machine_id!r} (connus : {list_machines()})")
    with open(path, encoding="utf-8") as f:
        return schema.from_jsonable(MachineProfile, json.load(f))
