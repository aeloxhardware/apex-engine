#!/usr/bin/env python3
"""Contrôle du registre des palliatifs (`docs/PALLIATIFS.md`, limite dure n° 8).

Chaque ligne PAL-xxx doit citer un test existant. La revue a lieu à chaque release (RULES §11.3 : plus de dates).

    python scripts/check_palliatives.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    problems = []
    rows = [ln for ln in (ROOT / "docs" / "PALLIATIFS.md").read_text(encoding="utf-8").splitlines()
            if ln.startswith("| PAL-")]
    if not rows:
        problems.append("aucun palliatif enregistré")
    for row in rows:
        cells = [c.strip() for c in row.strip("|").split("|")]
        ident, tests = cells[0], cells[3]
        refs = re.findall(r"`([^`]+)`", tests)
        if not refs:
            problems.append(f"{ident}: aucun test cité")
        for ref in refs:
            if ref.startswith("tests/") and not (ROOT / ref.split("::")[0]).is_file():
                problems.append(f"{ident}: test introuvable : {ref}")
    for problem in problems:
        print(problem)
    print(f"check_palliatives : {'ECHEC' if problems else 'ok'} ({len(rows)} palliatif(s))")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
