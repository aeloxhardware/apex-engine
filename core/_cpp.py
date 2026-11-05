"""
core/_cpp.py — chargement du noyau natif `apex_core` et interrupteur de bascule.

Role pendant la migration vers C++ (voir le plan de portage) : chaque fonction
portee garde son implementation Python A COTE de la version C++, et le site
d'appel choisit entre les deux. C'est ce qui permet aux 305 tests pytest de
servir d'ORACLE D'EQUIVALENCE — un test peut appeler les deux implementations
dans le meme processus et comparer leurs sorties.

Deux principes non negociables :

1. **Le Python reste la valeur par defaut.** `RuntimeOptions.use_cpp` doit etre positionne
   explicitement (par un point d'entree) pour activer le natif. Un portage ne s'active jamais tout seul
   avant d'avoir ete mesure — meme discipline que tous les mecanismes opt-in du
   projet (`support_driven_order`, `anchor_growth_weight`, ...).

2. **L'absence du module n'est jamais une erreur.** Si `cpp/lib/` ne contient
   pas de `.pyd` (personne n'a compile, ou une autre plateforme), tout retombe
   silencieusement sur Python. Lancer les tests ne doit jamais exiger un
   compilateur C++.

Ce fichier disparaitra a la bascule finale, en meme temps que les bindings.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Optional

from core.runtime import GRID_MODE_FAITHFUL as _MODE_FAITHFUL
from core.runtime import RuntimeOptions

_DEFAULT_LIB_DIR = Path(__file__).resolve().parent.parent / "cpp" / "lib"

# Options installees par un POINT D'ENTREE (`configure`) ; ni cette fonction ni le reste de la
# bibliotheque ne lisent `os.environ`. Par defaut : tout eteint (Python seul).
_options: RuntimeOptions = RuntimeOptions()

# Le module natif se charge a la PREMIERE utilisation, pas a l'import : `configure` peut ainsi
# choisir le dossier (cpp/build_dev/lib, voir scripts/build_cpp.py --dev) avant le chargement.
_module: Optional[Any] = None
_load_error: Optional[str] = None
_loaded_from: Optional[Path] = None
_loaded = False


def _lib_dir(options: RuntimeOptions) -> Path:
    return Path(options.lib_dir).resolve() if options.lib_dir else _DEFAULT_LIB_DIR


def _ensure_loaded() -> None:
    global _module, _load_error, _loaded_from, _loaded
    if _loaded:
        return
    _loaded = True
    _loaded_from = _lib_dir(_options)
    if _loaded_from.is_dir() and str(_loaded_from) not in sys.path:
        sys.path.insert(0, str(_loaded_from))
    try:
        import apex_core as module  # type: ignore[import-not-found]
        _module = module
    except ImportError as exc:  # pragma: no cover - depend de la presence du build
        _module = None
        _load_error = str(exc)


def configure(options: RuntimeOptions) -> RuntimeOptions:
    """
    Installe les options d'execution (point d'entree uniquement). Rend les options precedentes.

    Le dossier du module natif ne peut plus changer une fois le module charge : `lib_dir`
    doit etre fixe avant la premiere utilisation (leve sinon, au lieu d'ignorer en silence).
    """
    global _options
    if _loaded and options.lib_dir and _lib_dir(options) != _loaded_from:
        raise RuntimeError(
            f"Le noyau natif est deja charge depuis {_loaded_from} : lib_dir={options.lib_dir!r} "
            "doit etre configure avant sa premiere utilisation."
        )
    previous, _options = _options, options
    return previous


def options() -> RuntimeOptions:
    """Options d'execution installees."""
    return _options


# Surcharge explicite, prioritaire sur les options. `None` = suivre les options.
# Utilise par les tests d'equivalence, qui doivent exercer les DEUX chemins dans un meme processus.
_forced: Optional[bool] = None


def available() -> bool:
    """Le module natif est-il compile et importable ?"""
    _ensure_loaded()
    return _module is not None


def load_error() -> Optional[str]:
    """Message d'import si le module natif est absent, `None` sinon."""
    _ensure_loaded()
    return _load_error


def enabled() -> bool:
    """Faut-il emprunter le chemin C++ ? Faux tant que rien n'est compile."""
    _ensure_loaded()
    if _module is None:
        return False
    if _forced is not None:
        return _forced
    return _options.use_cpp


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


# --- Backend GPU (Vulkan), phase 1 du plan « portage C++, puis Vulkan » ---
# Meme discipline que le natif : DESACTIVE par defaut (`RuntimeOptions.use_gpu`), jamais
# une erreur s'il manque (module construit sans Vulkan, pas de pilote, pas de
# GPU) : `gpu_status()` dit pourquoi, et chaque site d'appel retombe sur le CPU.
_gpu_forced: Optional[bool] = None


def gpu_status() -> dict:
    """{ok, reason, device, fp64, timestamp_period_ns, compiled} ; ok = False avec
    la raison si le module natif est absent ou anterieur au socle Vulkan."""
    _ensure_loaded()
    if _module is None:
        return {"ok": False, "reason": f"module natif absent ({_load_error})", "compiled": False}
    if not hasattr(_module, "gpu_status"):
        return {"ok": False, "reason": "module natif anterieur au socle Vulkan", "compiled": False}
    return dict(_module.gpu_status())


def gpu_enabled() -> bool:
    """Faut-il emprunter le chemin GPU ? Exige le natif actif ET un GPU utilisable."""
    if not enabled():
        return False
    wanted = _gpu_forced if _gpu_forced is not None else _options.use_gpu
    return bool(wanted) and bool(gpu_status().get("ok"))


def set_gpu_enabled(value: Optional[bool]) -> Optional[bool]:
    """Force ou libere le choix GPU (tests) ; rend l'ancienne valeur."""
    global _gpu_forced
    previous = _gpu_forced
    _gpu_forced = value
    return previous


#: Modes de requete de `apex_core.PrintedGrid` — voir cpp/include/apex/printed_index.hpp.
GRID_MODE_FAITHFUL = 0
GRID_MODE_STAMPED = 1


def grid_mode() -> int:
    """
    Mode de requete de l'index spatial natif (`RuntimeOptions.grid_mode`).

    `faithful` transcrit exactement la structure Python (un ensemble de
    hachage par requete) ; `stamped` (defaut) dedoublonne par estampilles.
    Les deux existent pour MESURER SEPAREMENT le gain du changement de langage
    et celui de la correction algorithmique — un seul chiffre global ne
    permettrait pas de savoir lequel des deux a porte le resultat.
    """
    return GRID_MODE_FAITHFUL if _options.grid_mode == _MODE_FAITHFUL else GRID_MODE_STAMPED


def module() -> Any:
    """Le module natif. Leve si absent — a n'appeler qu'apres `enabled()`."""
    _ensure_loaded()
    if _module is None:
        raise RuntimeError(
            "Le noyau natif apex_core n'est pas disponible "
            f"({_load_error}). Compiler avec : python scripts/build_cpp.py"
        )
    return _module
