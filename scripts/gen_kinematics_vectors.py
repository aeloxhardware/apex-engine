"""
Genere tests/vectors/kinematics_vectors.json (R-008, RULES §6) : vecteurs de test de la cinematique
pivot/pointe, tires de `machine_to_firmware` / `firmware_to_machine` (file_io/gcode_export.py).

Chaque vecteur : parametres (pivot_to_tip, b_sign, tilt_sign, z_reference), pointe (x, y, z) et B d'Apex,
`previous_y` facultatif (Y deroule), sortie firmware (X, Y, Z, B). L'inverse (`firmware_to_machine`) doit
redonner la pointe, sauf au centre (R < singularity_radius) ou Y est herite. Les implementations de
l'interface (`gcodeParse.ts` `nozzleTip`, noyau `native/`) rejouent ce fichier.

Usage : python scripts/gen_kinematics_vectors.py [--check]
--check : n'ecrit rien, sort en erreur si le fichier n'est plus a jour.
"""
from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from file_io.gcode_export import machine_to_firmware  # noqa: E402

PATH = os.path.join(ROOT, "tests", "vectors", "kinematics_vectors.json")
SINGULARITY_RADIUS = 1e-6
TOLERANCE_MM = 1e-6

CONFIGS = [
    {"pivot_to_tip": 46.0, "b_sign": 1.0, "tilt_sign": 1.0, "z_reference": "tip"},    # defaut actuel
    {"pivot_to_tip": 46.0, "b_sign": 1.0, "tilt_sign": 1.0, "z_reference": "pivot"},
    {"pivot_to_tip": 45.0, "b_sign": -1.0, "tilt_sign": -1.0, "z_reference": "tip"},
]
POINTS = [
    ((30.0, 0.0, 5.0), 0.0, None), ((30.0, 0.0, 5.0), 45.0, None), ((30.0, 0.0, 5.0), 90.0, None),
    ((0.0, 25.0, 12.5), 135.0, None), ((-20.0, 0.0, 0.2), 185.0, None), ((-10.0, -10.0, 3.0), 60.0, None),
    ((12.0, 9.0, 7.0), 30.0, 700.0),          # Y deroule : reste dans la continuite de previous_y
    ((12.0, -9.0, 7.0), 30.0, -350.0),
    ((0.0, 0.0, 4.0), 20.0, 123.5),           # centre : Y herite
]


def build() -> str:
    vectors = []
    for cfg in CONFIGS:
        for position, b, prev in POINTS:
            x, y, z, fb = machine_to_firmware(position, b, cfg["pivot_to_tip"], previous_y=prev,
                                              singularity_radius=SINGULARITY_RADIUS, tilt_sign=cfg["tilt_sign"],
                                              z_reference=cfg["z_reference"], b_sign=cfg["b_sign"])
            vectors.append({**cfg, "tip": list(position), "b_apex": b, "previous_y": prev,
                            "firmware": {"X": round(x, 9), "Y": round(y, 9), "Z": round(z, 9), "B": round(fb, 9)},
                            "inverse_exact": not (abs(position[0]) < 1e-6 and abs(position[1]) < 1e-6)})
    doc = {"format": "apex-kinematics-vectors-1",
           "note": "Genere par scripts/gen_kinematics_vectors.py depuis file_io/gcode_export.py. NE PAS EDITER. "
                   "Unites : mm et degres. B d'Apex : 0 = tete vers le bas, jamais negatif. Y firmware = "
                   "atan2(y, x) deroule contre previous_y. Au centre (R < singularity_radius) Y est herite.",
           "singularity_radius": SINGULARITY_RADIUS, "tolerance_mm": TOLERANCE_MM, "vectors": vectors}
    return json.dumps(doc, indent=1, ensure_ascii=False) + "\n"


def main() -> int:
    text = build()
    if "--check" in sys.argv:
        current = open(PATH, encoding="utf-8").read() if os.path.exists(PATH) else ""
        if current != text:
            print("kinematics_vectors.json n'est plus a jour : relancer scripts/gen_kinematics_vectors.py")
            return 1
        return 0
    os.makedirs(os.path.dirname(PATH), exist_ok=True)
    with open(PATH, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
