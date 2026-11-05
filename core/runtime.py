"""
core/runtime.py — options d'exécution du moteur (R-030 point 3).

Aucune variable d'environnement ne change le résultat d'un calcul : la bibliothèque
(`core/`, `element-*/`, `file_io/`) ne lit JAMAIS `os.environ`. Les choix d'exécution
(noyau natif, GPU, mode de requête de l'index, worker TetGen) tiennent dans cet objet
immuable, construit par un POINT D'ENTRÉE (`main.py`, `slice_cli.py`, un script, le
`conftest.py` des tests) puis installé par `core._cpp.configure`.

`options_from_environment` est la seule fonction qui lit l'environnement : c'est une
commodité de ligne de commande des points d'entrée, pas un chemin de la bibliothèque.

Dette ouverte (R-030 point 3, passe 3) : l'objet est encore installé dans un état de module
de `core._cpp` au lieu d'être passé en paramètre à chaque étape.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Optional

_TRUE_VALUES = {"1", "true", "yes", "on"}

GRID_MODE_STAMPED = "stamped"
GRID_MODE_FAITHFUL = "faithful"


@dataclass(frozen=True)
class RuntimeOptions:
    use_cpp: bool = False                 # emprunter le noyau natif (apex_core)
    use_gpu: bool = False                 # emprunter le backend Vulkan (exige use_cpp)
    grid_mode: str = GRID_MODE_STAMPED    # requête de l'index natif : "stamped" ou "faithful"
    lib_dir: Optional[str] = None         # dossier du module natif et du worker TetGen (None : cpp/lib)
    tetgen_worker: bool = True            # worker TetGen C++ s'il est construit (False : worker Python)
    tetgen_log: Optional[str] = None      # fichier où journaliser les temps par phase de TetGen

    def __post_init__(self) -> None:
        if self.grid_mode not in (GRID_MODE_STAMPED, GRID_MODE_FAITHFUL):
            raise ValueError(f"grid_mode inconnu : {self.grid_mode!r}")


def options_from_environment(environ: Optional[Mapping[str, str]] = None) -> RuntimeOptions:
    """Lit APEX_USE_CPP, APEX_USE_GPU, APEX_CPP_GRID_MODE, APEX_CPP_LIB_DIR, APEX_TETGEN_WORKER
    et APEX_TETGEN_LOG. À n'appeler que depuis un point d'entrée."""
    if environ is None:
        import os
        environ = os.environ

    def flag(name: str) -> bool:
        return environ.get(name, "").strip().lower() in _TRUE_VALUES

    grid = environ.get("APEX_CPP_GRID_MODE", "").strip().lower()
    return RuntimeOptions(
        use_cpp=flag("APEX_USE_CPP"),
        use_gpu=flag("APEX_USE_GPU"),
        grid_mode=GRID_MODE_FAITHFUL if grid == GRID_MODE_FAITHFUL else GRID_MODE_STAMPED,
        lib_dir=environ.get("APEX_CPP_LIB_DIR") or None,
        tetgen_worker=environ.get("APEX_TETGEN_WORKER", "1") != "0",
        tetgen_log=environ.get("APEX_TETGEN_LOG") or None,
    )


def configure_from_environment() -> RuntimeOptions:
    """Point d'entrée : lit l'environnement et installe les options. Rend les options installées."""
    from core import _cpp
    options = options_from_environment()
    _cpp.configure(options)
    return options
