"""
core/runtime.py — options d'execution EXPLICITES du moteur et etat de run (R-030 point 3).

Remplace les variables d'environnement (`APEX_USE_CPP`, `APEX_USE_GPU`, `APEX_CPP_GRID_MODE`,
`APEX_TETGEN_WORKER`, `APEX_TETGEN_LOG`) et les `global` mutables (compteurs de repli, cache du
controleur de collision, progression) : plus rien ne se lit dans l'environnement du processus.

  - `RuntimeOptions` : dataclass figee, sérialisable (porte par `PipelineParams.runtime`, donc par le JSON
    de la requete). C'est la SEULE source de ce qui change le calcul (natif ou Python, GPU, mode de grille,
    execution speculative de TetGen, nombre de fils de step6).
  - `RuntimeState` : ce qui s'accumule PENDANT un run (compteurs de repli, cache du controleur rapide,
    progression). Un run = un etat, jamais partage entre deux runs.
  - `use(options)` ouvre un run (options + etat neufs) pour la duree d'un bloc `with`. Le run est porte par
    un `contextvars.ContextVar` : deux runs concurrents (fils, interface) ne se voient pas.
    Les fils d'execution doivent etre lances avec `bind(fn)` pour heriter du run courant.

Hors de tout `use(...)`, un run implicite par defaut existe (options par defaut : Python, pas de GPU) : c'est
ce que voient les tests unitaires et les appels isoles.

Seul reste dans l'environnement `APEX_CPP_LIB_DIR` (core/_cpp.py) : l'EMPLACEMENT du module natif a l'import,
qui ne change pas le resultat mais le binaire charge (voir scripts/build_cpp.py --dev).
"""
from __future__ import annotations

import contextlib
import contextvars
import dataclasses
import threading
from dataclasses import dataclass, field
from typing import Any, Callable, Optional


@dataclass(frozen=True)
class RuntimeOptions:
    """Ce qui change la facon de calculer, jamais le contrat de sortie."""

    native: bool = False                   # noyau C++ (apex_core) quand il est compile
    gpu: bool = False                      # backend Vulkan (exige `native`)
    grid_mode: str = "stamped"             # "stamped" | "faithful" : mode de requete de l'index spatial natif
    tetgen_worker: bool = True             # TetGen dans l'executable separe quand il existe
    tetgen_log: Optional[str] = None       # fichier ou TetGen journalise ses phases (diagnostic seulement)
    speculative_tetgen: bool = True        # tentative finale lancee en parallele de la sonde (step1)
    step6_threads: Optional[int] = None    # fils de step6 ; None = min(8, coeurs)


@dataclass
class RuntimeState:
    """Etat d'un run : compteurs, caches, progression. Jamais partage entre deux runs."""

    options: RuntimeOptions = field(default_factory=RuntimeOptions)
    # step8 : replis de la resolution analytique vers la grille (None = chemin natif non emprunte)
    native_fallbacks: Optional[int] = None
    # step5 : contours rendus par trimesh faute de pouvoir etre traces en natif
    outline_fallbacks: int = 0
    # step9 : un seul controleur rapide, reconstruit si ses entrees changent
    fast_checker_entry: Optional[tuple] = None
    # progression (core/progress.py)
    progress_sink: Optional[Callable[[float], None]] = None
    stage_switcher: Optional[Callable[[str], None]] = None
    progress_last_fraction: float = 0.0
    progress_last_time: float = 0.0
    progress_scale: tuple = (0.0, 1.0)
    lock: threading.Lock = field(default_factory=threading.Lock, repr=False)


_current: contextvars.ContextVar[Optional[RuntimeState]] = contextvars.ContextVar("apex_runtime", default=None)


def state() -> RuntimeState:
    """L'etat du run courant (le run implicite par defaut s'il n'y en a pas d'ouvert)."""
    current = _current.get()
    if current is None:
        current = RuntimeState()
        _current.set(current)
    return current


def options() -> RuntimeOptions:
    return state().options


def replace_options(**changes: Any) -> RuntimeOptions:
    """Change des options du run courant ; rend les anciennes (pour les restaurer). Pour les tests et les
    outils de mesure : le code du pipeline ne l'appelle jamais."""
    current = state()
    previous = current.options
    current.options = dataclasses.replace(previous, **changes)
    return previous


@contextlib.contextmanager
def use(opts: Optional[RuntimeOptions] = None):
    """Ouvre un run (options et etat neufs) pour la duree du bloc."""
    token = _current.set(RuntimeState(options=opts if opts is not None else RuntimeOptions()))
    try:
        yield _current.get()
    finally:
        _current.reset(token)


def bind(fn: Callable) -> Callable:
    """`fn` s'executera dans le run courant meme depuis un autre fil (ThreadPoolExecutor.submit/map)."""
    # Le run implicite doit exister AVANT la copie pour que le fil et l'appelant partagent le meme etat.
    state()
    context = contextvars.copy_context()

    def bound(*args, **kwargs):
        return context.copy().run(fn, *args, **kwargs)

    return bound
