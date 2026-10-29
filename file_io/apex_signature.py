"""Signature d'un G-code Apex : une ligne fixe, en clair, qui ne signifie rien.

L'aperçu (interface) la cherche dans les premières lignes pour savoir que le fichier vient du
moteur Apex ; il lit alors le bloc d'en-tête chiffré (`file_io/apex_header.py`, format v1) et y
trouve la version, la machine, la cinématique. La chaîne est volontairement sans sens : rien à
deviner ni à imiter par hasard. Ce n'est PAS une protection (le dépôt est public), seulement un
repère fiable, moins fragile que de classer un fichier d'après les mots (B, C) qu'il contient.

La même constante est écrite par les deux transcripteurs (polaire et non planaire) et copiée côté
interface (`slicer-ux`). La changer casse la reconnaissance des fichiers déjà écrits.
"""
from __future__ import annotations

import re
from pathlib import Path

SIGNATURE = "; 85Ejj-sV6dI-NXTGl-JiodA-0b8aZ"
SCAN_LINES = 40                      # l'interface cherche dans les 40 premières lignes
_VERSION_FILE = Path(__file__).resolve().parent.parent / "config" / "APEX_VERSION"
_VERSION_RE = re.compile(r"^[a-z]+-\d+(\.\d+){7}$")


def has_signature(lines) -> bool:
    for n, line in enumerate(lines):
        if n >= SCAN_LINES:
            break
        if line.rstrip() == SIGNATURE:
            return True
    return False


def apex_version() -> str:
    """Numéro de version de l'application (`apex/VERSION`), copié dans `config/APEX_VERSION` à chaque
    montée de version (RULES §9.1). `unknown` si le fichier manque."""
    try:
        v = _VERSION_FILE.read_text(encoding="utf-8").strip()
    except OSError:
        return "unknown"
    return v if _VERSION_RE.match(v) else "unknown"


def common_meta(generator: str) -> dict:
    """Champs de `meta` communs aux deux transcripteurs."""
    from config.schema import APEX_CONTRACT_VERSION
    return {"generator": generator, "apex_version": apex_version(), "contract": APEX_CONTRACT_VERSION}
