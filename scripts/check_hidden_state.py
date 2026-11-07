#!/usr/bin/env python3
"""Contrôle de l'état caché (apex/integration/RULES.md §11, limite dure n° 3 et 4 ; R-030 points 3 et 4).

Échoue (code 1) si le code du moteur (hors tests, scripts et recherche) contient :
  - une instruction `global` ;
  - une lecture d'environnement `APEX_*` autre que celles autorisées (EMPLACEMENT du module natif) ;
  - un nom de machine ou de tête écrit en dur (`"aetheris-4"`) ou un rayon de plateau par défaut.

    python scripts/check_hidden_state.py
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCOPE = ["core", "element_ardor", "element_aeolis", "file_io", "config", "main.py", "slice_cli.py"]
ALLOWED_ENV = {"APEX_CPP_LIB_DIR"}          # emplacement du binaire natif, ne change pas le résultat
MACHINE_LITERALS = re.compile(r"""["']aetheris-4["']""")
CONSTANTS = re.compile(r"\bPLATFORM_(RADIUS|HALF_HEIGHT|CENTER_OFFSET)\b")


def _files():
    for entry in SCOPE:
        path = ROOT / entry
        if path.is_file():
            yield path
        else:
            yield from sorted(p for p in path.rglob("*.py") if "__pycache__" not in p.parts)


def main() -> int:
    problems: list[str] = []
    for path in _files():
        rel = path.relative_to(ROOT).as_posix()
        text = path.read_text(encoding="utf-8", errors="replace")
        tree = ast.parse(text)
        for node in ast.walk(tree):
            if isinstance(node, ast.Global):
                problems.append(f"{rel}:{node.lineno}: `global {', '.join(node.names)}`")
            elif isinstance(node, ast.Constant) and isinstance(node.value, str):
                for name in re.findall(r"\bAPEX_[A-Z_]+\b", node.value):
                    # chaînes qui sont des clés d'environnement (pas la prose des docstrings)
                    if re.fullmatch(r"APEX_[A-Z_]+", node.value) and name not in ALLOWED_ENV \
                            and name not in {"APEX_CONTRACT_VERSION", "APEX_VERSION"}:
                        problems.append(f"{rel}:{node.lineno}: variable d'environnement {name}")
        for number, line in enumerate(text.splitlines(), 1):
            code = line.split("#", 1)[0]
            if MACHINE_LITERALS.search(code):
                problems.append(f"{rel}:{number}: nom de machine/tête en dur")
            if CONSTANTS.search(code):
                problems.append(f"{rel}:{number}: constante de plateau en dur")
    for problem in problems:
        print(problem)
    print(f"check_hidden_state : {'ECHEC, ' + str(len(problems)) + ' problème(s)' if problems else 'ok'}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
