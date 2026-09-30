"""
core/_cpp.py — chargement du noyau natif `apex_core` et interrupteur de bascule.

Role pendant la migration vers C++ (voir le plan de portage) : chaque fonction
portee garde son implementation Python A COTE de la version C++, et le site
d'appel choisit entre les deux. C'est ce qui permet aux 305 tests pytest de
servir d'ORACLE D'EQUIVALENCE — un test peut appeler les deux implementations
dans le meme processus et comparer leurs sorties.

Deux principes non negociables :

1. **Le Python reste la valeur par defaut.** `APEX_USE_CPP` doit etre positionne
   explicitement pour activer le natif. Un portage ne s'active jamais tout seul
   avant d'avoir ete mesure — meme discipline que tous les mecanismes opt-in du
   projet (`support_driven_order`, `anchor_growth_weight`, ...).

2. **L'absence du module n'est jamais une erreur.** Si `cpp/lib/` ne contient
   pas de `.pyd` (personne n'a compile, ou une autre plateforme), tout retombe
   silencieusement sur Python. Lancer les tests ne doit jamais exiger un
   compilateur C++.

Ce fichier disparaitra a la bascule finale, en meme temps que les bindings.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Optional

_LIB_DIR = Path(__file__).resolve().parent.parent / "cpp" / "lib"

_module: Optional[Any] = None
_load_error: Optional[str] = None

if _LIB_DIR.is_dir():
    lib_path = str(_LIB_DIR)
    if lib_path not in sys.path:
        sys.path.insert(0, lib_path)

try:
    import apex_core as _module  # type: ignore[no-redef]
except ImportError as exc:  # pragma: no cover - depend de la presence du build
    _module = None
    _load_error = str(exc)


# Surcharge explicite, prioritaire sur la variable d'environnement.
# `None` = suivre l'environnement. Utilise par les tests d'equivalence, qui
# doivent exercer les DEUX chemins dans un meme processus.
_forced: Optional[bool] = None

_TRUE_VALUES = {"1", "true", "yes", "on"}


def available() -> bool:
    """Le module natif est-il compile et importable ?"""
    return _module is not None


def load_error() -> Optional[str]:
    """Message d'import si le module natif est absent, `None` sinon."""
    return _load_error


def enabled() -> bool:
    """Faut-il emprunter le chemin C++ ? Faux tant que rien n'est compile."""
    if _module is None:
        return False
    if _forced is not None:
        return _forced
    return os.environ.get("APEX_USE_CPP", "").strip().lower() in _TRUE_VALUES


def set_enabled(value: Optional[bool]) -> Optional[bool]:
    """
    Force (`True`/`False`) ou libere (`None`) le choix d'implementation.

    Retourne l'ancienne valeur, pour que l'appelant puisse la restaurer —
    indispensable dans un test, qui ne doit pas laisser fuiter son reglage.
    """
    global _forced
    previous = _forced
    _forced = value
    return previous


#: Modes de requete de `apex_core.PrintedGrid` — voir cpp/include/apex/printed_index.hpp.
GRID_MODE_FAITHFUL = 0
GRID_MODE_STAMPED = 1


def grid_mode() -> int:
    """
    Mode de requete de l'index spatial natif, via `APEX_CPP_GRID_MODE`.

    `faithful` transcrit exactement la structure Python (un ensemble de
    hachage par requete) ; `stamped` (defaut) dedoublonne par estampilles.
    Les deux existent pour MESURER SEPAREMENT le gain du changement de langage
    et celui de la correction algorithmique — un seul chiffre global ne
    permettrait pas de savoir lequel des deux a porte le resultat.
    """
    if os.environ.get("APEX_CPP_GRID_MODE", "").strip().lower() == "faithful":
        return GRID_MODE_FAITHFUL
    return GRID_MODE_STAMPED


def module() -> Any:
    """Le module natif. Leve si absent — a n'appeler qu'apres `enabled()`."""
    if _module is None:
        raise RuntimeError(
            "Le noyau natif apex_core n'est pas disponible "
            f"({_load_error}). Compiler avec : python scripts/build_cpp.py"
        )
    return _module
