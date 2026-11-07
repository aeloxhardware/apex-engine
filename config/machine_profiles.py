"""
config/machine_profiles.py — profils de machine en DONNEES (R-030 point 4).

Une machine = un fichier `res/machines/<id>.json`. Ajouter une machine ou une tete ne demande aucune modification
du code : le coeur ne connait ni `aetheris-4`, ni le rayon du plateau, ni les limites d'axes. L'identifiant par
defaut est lui aussi une donnee (`res/machines/index.json`).
"""
from __future__ import annotations

import json
import math
import os
import re
from typing import Optional

import numpy as np

from core.mesh_types import MachineProfile, RotationAxisSpec

MACHINES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "res", "machines")
_ID = re.compile(r"[A-Za-z0-9_-]+")


def default_machine_id() -> str:
    with open(os.path.join(MACHINES_DIR, "index.json"), encoding="utf-8") as handle:
        return str(json.load(handle)["default"])


def machine_ids() -> list[str]:
    return sorted(n[:-5] for n in os.listdir(MACHINES_DIR) if n.endswith(".json") and n != "index.json")


def load_machine_profile(machine_id: Optional[str] = None, *, with_toolhead: bool = True) -> MachineProfile:
    """Construit le `MachineProfile` d'une machine. `with_toolhead` charge le GLB de la tete ; un echec de
    chargement laisse `nozzle_model=None` (la simulation de collision se desactive, jamais un plantage)."""
    machine_id = machine_id or default_machine_id()
    if not _ID.fullmatch(machine_id):
        raise ValueError(f"identifiant de machine invalide : {machine_id!r}")
    with open(os.path.join(MACHINES_DIR, f"{machine_id}.json"), encoding="utf-8") as handle:
        data = json.load(handle)
    nozzle_model = None
    if with_toolhead:
        try:
            from file_io.mesh_loader import load_toolhead_by_id
            nozzle_model = load_toolhead_by_id(data["toolhead_id"])
        except Exception as exc:  # pragma: no cover - chemin de repli
            import sys
            print(f"[machine] Avertissement : tete {data['toolhead_id']} non chargee ({exc}) - "
                  f"simulation de collision desactivee.", file=sys.stderr)
    return MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array(a["axis"], dtype=float),
                                        angle_range_deg=tuple(float(v) for v in a["angle_range_deg"]))
                       for a in data["rotation_axes"]],
        platform_origin=np.array(data["platform_origin"], dtype=float),
        nozzle_model=nozzle_model,
        toolhead_id=data["toolhead_id"],
        nozzle_diameter=float(data["nozzle_diameter"]),
        max_angular_velocity=float(math.radians(data["max_angular_velocity_deg_s"])),
        max_angular_acceleration=float(math.radians(data["max_angular_acceleration_deg_s2"])),
        platform_radius=float(data["platform_radius"]),
        platform_thickness=float(data["platform_thickness"]),
        position_coupled_axes=[int(i) for i in data["position_coupled_axes"]],
        radial_direction=np.array(data["radial_direction"], dtype=float),
    )


def machine_data(machine_id: Optional[str] = None) -> dict:
    """Contenu brut du fichier de profil (sans construire de MachineProfile)."""
    machine_id = machine_id or default_machine_id()
    if not _ID.fullmatch(machine_id):
        raise ValueError(f"identifiant de machine invalide : {machine_id!r}")
    with open(os.path.join(MACHINES_DIR, f"{machine_id}.json"), encoding="utf-8") as handle:
        return json.load(handle)


def default_toolhead_id() -> str:
    """Tete de la machine par defaut (donnee, pas une constante du code)."""
    return str(machine_data()["toolhead_id"])
