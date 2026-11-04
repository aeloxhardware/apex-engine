#!/usr/bin/env python3
"""Écrit docs/LICENCES_DEPENDANCES.md : licence de chaque dépendance figée dans requirements.lock.

Lit les métadonnées des paquets INSTALLÉS (importlib.metadata) : à lancer dans l'environnement construit depuis le
verrou. Une licence introuvable est écrite « inconnue », jamais devinée (R-030 point 5, RULES §11.1).

    python scripts/gen_licences.py            écrit le fichier
    python scripts/gen_licences.py --check    échoue si le fichier n'est pas à jour (CI)
"""
from __future__ import annotations

import re
import sys
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LOCK = ROOT / "requirements.lock"
OUT = ROOT / "docs" / "LICENCES_DEPENDANCES.md"
# Licences à copyleft fort, à signaler en tête : le projet est libre (RULES, décision de l'utilisateur) et les accepte.
COPYLEFT = re.compile(r"Affero|AGPL|(?<![LA])GPL|(?<!Lesser )General Public License|GPL3", re.I)
# Licence DÉCLARÉE par le paquet différente de ce qu'il embarque : l'enveloppe PyPI de TetGen est sous MIT mais compile
# les sources de TetGen (AGPL-3.0). On l'écrit ici plutôt que de laisser croire « MIT » (voir CLAUDE.md, section licences).
EMBEDDED = {"tetgen": "MIT (enveloppe) ; embarque TetGen, AGPL-3.0"}


def locked() -> list[tuple[str, str]]:
    out = []
    for line in LOCK.read_text(encoding="utf-8").splitlines():
        m = re.match(r"^([A-Za-z0-9_.\-]+)==([^\s\\]+)", line)
        if m:
            out.append((m.group(1), m.group(2)))
    return sorted(out, key=lambda x: x[0].lower())


def licence_of(name: str) -> str:
    if name.lower() in EMBEDDED:
        return EMBEDDED[name.lower()]
    try:
        md = metadata.metadata(name)
    except metadata.PackageNotFoundError:
        return "inconnue (paquet non installé)"
    expr = (md.get("License-Expression") or "").strip()
    if expr:
        return expr
    lic = (md.get("License") or "").strip()
    if lic and len(lic) < 80 and "\n" not in lic:
        return lic
    classes = [c.split("::")[-1].strip() for c in md.get_all("Classifier") or [] if c.startswith("License ::")]
    return "; ".join(classes) if classes else "inconnue"


def render() -> str:
    rows = [(n, v, licence_of(n)) for n, v in locked()]
    lines = [
        "# Licences des dépendances",
        "",
        "Généré par `scripts/gen_licences.py` depuis `requirements.lock` : ne pas éditer à la main.",
        "Le logiciel est libre et gratuit (décision de l'utilisateur) : les licences à copyleft fort ci-dessous sont acceptées.",
        "Ce tableau n'est pas un avis juridique.",
        "",
        "| paquet | version figée | licence déclarée |",
        "|---|---|---|",
    ]
    for n, v, lic in rows:
        flag = " **(copyleft fort)**" if COPYLEFT.search(lic) else ""
        lines.append(f"| {n} | {v} | {lic}{flag} |")
    lines += [
        "",
        "Hors du verrou Python : TetGen est aussi compilé depuis ses sources (`cpp/third_party/tetgen/`, AGPL-3.0), Clipper2",
        "(`cpp/third_party/clipper2/`, BSL-1.0) et Eigen (`cpp/third_party/eigen/`, MPL-2.0) ; leurs licences sont dans leurs dossiers.",
        "",
    ]
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    text = render()
    if "--check" in argv:
        if not OUT.exists() or OUT.read_text(encoding="utf-8") != text:
            print("docs/LICENCES_DEPENDANCES.md n'est pas à jour (python scripts/gen_licences.py)", file=sys.stderr)
            return 1
        print("licences à jour")
        return 0
    OUT.write_text(text, encoding="utf-8")
    print(f"écrit : {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
