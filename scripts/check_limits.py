#!/usr/bin/env python3
"""Contrôle des limites dures (apex/integration/RULES.md §11.1) avec cliquet.

Seuils : fichier source <= 800 lignes, fonction <= 80 lignes, signature <= 8 paramètres.
`scripts/limits_baseline.json` liste la dette existante : une valeur ne peut que diminuer et
aucune nouvelle entrée n'est admise.

    python scripts/check_limits.py            contrôle (code 1 au moindre dépassement nouveau)
    python scripts/check_limits.py --update   abaisse / retire les entrées réduites (n'en ajoute jamais)
    python scripts/check_limits.py --report   rapport de dette (pour apex/RELEASES.md)
    python scripts/check_limits.py --init     écrit la dette de départ (une seule fois)
"""
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BASELINE = Path(__file__).resolve().parent / "limits_baseline.json"
MAX_FILE, MAX_FN, MAX_PARAMS = 800, 80, 8
# Code tiers, sorties et environnements : hors du code d'Apex.
SKIP_PARTS = {".git", ".venv", "venv", "__pycache__", "node_modules", "out", "build", "third_party", "research"}


def _files() -> list[Path]:
    out = []
    for p in sorted(ROOT.rglob("*.py")):
        rel = p.relative_to(ROOT)
        if any(part in SKIP_PARTS for part in rel.parts):
            continue
        out.append(p)
    return out


def _measure(path: Path) -> dict[str, int]:
    rel = path.relative_to(ROOT).as_posix()
    text = path.read_text(encoding="utf-8", errors="replace")
    result = {f"file:{rel}": len(text.splitlines())}
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return result
    seen: dict[str, int] = {}

    def visit(node: ast.AST, scope: str) -> None:
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                label = f"{scope}.{child.name}" if scope else child.name
                n = seen.get(label, 0) + 1
                seen[label] = n
                key = f"{rel}::{label}#{n}"
                length = (child.end_lineno or child.lineno) - child.lineno + 1
                if length > MAX_FN:
                    result[f"fn:{key}"] = length
                a = child.args
                names = [x.arg for x in [*a.posonlyargs, *a.args, *a.kwonlyargs]]
                if a.vararg:
                    names.append(a.vararg.arg)
                if a.kwarg:
                    names.append(a.kwarg.arg)
                count = len([x for x in names if x not in ("self", "cls")])
                if count > MAX_PARAMS:
                    result[f"params:{key}"] = count
                visit(child, label)
            elif isinstance(child, ast.ClassDef):
                visit(child, f"{scope}.{child.name}" if scope else child.name)
            else:
                visit(child, scope)

    visit(tree, "")
    return result


def _limit_of(key: str) -> int:
    if key.startswith("file:"):
        return MAX_FILE
    return MAX_PARAMS if key.startswith("params:") else MAX_FN


def current() -> dict[str, int]:
    out: dict[str, int] = {}
    for p in _files():
        for k, v in _measure(p).items():
            if v > _limit_of(k):
                out[k] = v
    return out


def _dump(d: dict[str, int]) -> None:
    BASELINE.write_text(json.dumps(dict(sorted(d.items())), indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    mode = argv[1] if len(argv) > 1 else ""
    cur = current()
    if mode == "--init":
        _dump(cur)
        print(f"baseline écrite : {len(cur)} dépassements")
        return 0
    if not BASELINE.exists():
        print("limits_baseline.json absent (python scripts/check_limits.py --init)", file=sys.stderr)
        return 2
    base = json.loads(BASELINE.read_text(encoding="utf-8"))
    fresh = [f"{k} = {v}" for k, v in cur.items() if k not in base]
    worse = [f"{k} : {base[k]} -> {v}" for k, v in cur.items() if k in base and v > base[k]]
    better = [k for k, v in cur.items() if k in base and v < base[k]]
    fixed = [k for k in base if k not in cur]
    if mode == "--report":
        by = lambda p: sum(1 for k in cur if k.startswith(p))
        print(f"fichiers > {MAX_FILE} : {by('file:')}\nfonctions > {MAX_FN} : {by('fn:')}\nsignatures > {MAX_PARAMS} : {by('params:')}")
        for k, v in sorted((kv for kv in cur.items() if kv[0].startswith("file:")), key=lambda kv: -kv[1]):
            print(f"  {v}\t{k[5:]}")
        return 0
    if mode == "--update":
        nxt = dict(base)
        for k in fixed:
            del nxt[k]
        for k in better:
            nxt[k] = cur[k]
        _dump(nxt)
        print(f"cliquet resserré : {len(fixed)} retirées, {len(better)} réduites")
        if fresh or worse:
            print("des régressions subsistent (non ajoutées)", file=sys.stderr)
            return 1
        return 0
    if fresh or worse:
        print("check_limits : RÉGRESSION", file=sys.stderr)
        for line in fresh:
            print("  nouveau  " + line, file=sys.stderr)
        for line in worse:
            print("  aggravé  " + line, file=sys.stderr)
        return 1
    print(f"check_limits : ok ({len(cur)} dettes, {len(better) + len(fixed)} à resserrer avec --update)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
