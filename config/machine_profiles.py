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
    profile = MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array(a["axis"], dtype=float),
                                        angle_range_deg=tuple(float(v) for v in a["angle_range_deg"]))
                       for a in data["rotation_axes"]],
        platform_origin=np.array(data["platform_origin"], dtype=float),
        nozzle_model=None,
        toolhead_id=data["toolhead_id"],
        nozzle_diameter=float(data["nozzle_diameter"]),
        max_angular_velocity=float(math.radians(data["max_angular_velocity_deg_s"])),
        max_angular_acceleration=float(math.radians(data["max_angular_acceleration_deg_s2"])),
        platform_radius=float(data["platform_radius"]),
        platform_thickness=float(data["platform_thickness"]),
        position_coupled_axes=[int(i) for i in data["position_coupled_axes"]],
        radial_direction=None if data.get("radial_direction") is None else np.array(data["radial_direction"], dtype=float),
        kinematics=str(data.get("kinematics", "polar")),
        platform_size=None if data.get("platform_size") is None else tuple(float(v) for v in data["platform_size"]),
        pivot_to_tip=None if data.get("pivot_to_tip") is None else float(data["pivot_to_tip"]),
        toolhead_envelope=data.get("toolhead_envelope"),
    )
    if with_toolhead:
        ensure_toolhead(profile)
    return profile


def ensure_toolhead(machine_profile: MachineProfile) -> None:
    """Charge la tete que designe `toolhead_id` quand le profil n'en porte pas (profil venu d'une requete
    JSON : `nozzle_model` n'est pas serialise, R-020) ; sans `toolhead_id`, l'enveloppe parametrique du profil
    (`toolhead_envelope`, limite #97). Echec : journal et `nozzle_model=None`, jamais un crash
    (step9 degrade en `nozzle_model_checked=False`)."""
    import sys
    if machine_profile.nozzle_model is not None:
        return
    if not machine_profile.toolhead_id:
        if machine_profile.toolhead_envelope:
            from file_io.toolhead_envelope import envelope_mesh
            machine_profile.nozzle_model = envelope_mesh(machine_profile.toolhead_envelope)
        return
    try:
        from file_io.mesh_loader import load_toolhead_by_id

        machine_profile.nozzle_model = load_toolhead_by_id(machine_profile.toolhead_id)
    except Exception as exc:
        print(f"[main] Avertissement : tete {machine_profile.toolhead_id!r} non chargee ({exc}) - "
              f"simulation de collision desactivee (nozzle_model=None).", file=sys.stderr)


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


def machines_main(argv: list[str]) -> int:
    """`main.py machines` : profils disponibles en JSON (identifiant, nom, architecture, profil complet) et machine
    par defaut. L'interface en tire son choix de machine au lieu de l'inventer (R-045)."""
    import argparse
    from config.schema import APEX_CONTRACT_VERSION
    argparse.ArgumentParser(prog="main.py machines", description="Machines connues du moteur.").parse_args(argv)
    machines = []
    for machine_id in machine_ids():
        if not _ID.fullmatch(machine_id):          # fichier d'exemple (`*.example.json`), pas un profil
            continue
        data = machine_data(machine_id)
        machines.append({"id": machine_id, "name": data.get("name", machine_id),
                         "kinematics": data.get("kinematics", "polar"), "profile": data})
    print(json.dumps({"contract": APEX_CONTRACT_VERSION, "default": default_machine_id(), "machines": machines},
                     ensure_ascii=False), flush=True)
    return 0
