"""Temps d'un G-code selon les limites par axe du firmware (config.g actif).

Modele : chaque mouvement dure max(duree nominale a la vitesse de pointe, |d_axe| /
vitesse max de l'axe). Les accelerations ne sont pas modelisees ici (borne basse).
Usage : estimate_gcode_time.py <fichier.gcode> [--pivot-to-tip=<mm>]
"""
import sys
import numpy as np
sys.path.insert(0, ".")
from file_io.gcode_export import parse_gcode, firmware_to_machine  # noqa: F401
from file_io.pivot_length import current_pivot_to_tip  # noqa: E402

# L du G-code : derniere valeur de config/pivot_to_tip.json, sinon 46 mm ; --pivot-to-tip=<mm> l'impose.
_flags = {a.split('=')[0]: a.split('=', 1)[1] for a in sys.argv[2:] if a.startswith('--') and '=' in a}
PIVOT = float(_flags['--pivot-to-tip']) if '--pivot-to-tip' in _flags else current_pivot_to_tip()[0]

# config.g : M203 (mm/min -> /s). X 100000, Y 36000, Z 400, B 100000.
LIMITS = {"config.g": (1666.7, 600.0, 400 / 60, 1666.7),
          "config.g sans Z": (1666.7, 600.0, 1e9, 1666.7),
          "4axisinput.g": (333.3, 360.0, 133.3, 360.0)}

moves = np.array([m[:4] for m in parse_gcode(open(sys.argv[1], encoding="utf-8"))], float)
ext = np.array([m[4] for m in parse_gcode(open(sys.argv[1], encoding="utf-8"))])
d = np.abs(np.diff(moves, axis=0))
# Pointe en 3D machine (inverse du pivot a PIVOT mm), longueur reelle de chaque mouvement.
a = np.radians(moves[:, 1])
sr = np.sin(np.radians(moves[:, 3]))
cz = np.cos(np.radians(moves[:, 3]))
rr = moves[:, 0] - PIVOT * sr
tipxyz = np.stack([rr * np.cos(a), rr * np.sin(a), moves[:, 2] - PIVOT * (cz - 1.0)], axis=1)
tip = np.linalg.norm(np.diff(tipxyz, axis=0), axis=1)
pr = ext[1:]
nominal = np.where(pr, tip / 25.0, tip / 80.0)
print(f"{sys.argv[1]}: {len(d)} mouvements, nominal {nominal.sum()/3600:.2f} h")
for name, (vx, vy, vz, vb) in LIMITS.items():
    t = np.maximum.reduce([nominal, d[:, 0] / vx, d[:, 1] / vy, d[:, 2] / vz, d[:, 3] / vb])
    over = (t > nominal * 1.001)
    print(f"  {name:16s}: {t.sum()/3600:6.2f} h  (mouvements ralentis {100*over.mean():.1f} %)")
