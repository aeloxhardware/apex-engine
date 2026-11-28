"""Processus de mise en texte d'un morceau de G-code XYZAB (`file_io/gcode_xyzab.py`, `_text`).

Lance par `python -m file_io.gcode_text_worker` : il ne reimporte jamais le script de l'appelant (un pool de
`multiprocessing` sous Windows le reimporte, et un script sans garde `__main__` relancait alors tout le tranchage dans
chaque processus, 2026-11-28). Entree : les colonnes du morceau (`np.savez`) sur l'entree standard ; sortie : le texte
ASCII des lignes de mouvement."""
from __future__ import annotations

import io
import sys

import numpy as np


def main() -> None:
    from file_io.gcode_xyzab import COLUMNS, _format_chunk
    data = np.load(io.BytesIO(sys.stdin.buffer.read()))
    cols = tuple(data[name] for name in COLUMNS)
    sys.stdout.buffer.write(_format_chunk((cols, int(data["decimals"]))).encode("ascii"))


if __name__ == "__main__":
    main()
