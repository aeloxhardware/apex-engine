"""
Progression REELLE a l'interieur d'une etape (R-017, RULES §5.3).

`run_pipeline` ouvre une etape avec `begin(sink)`; les boucles de l'etape appellent `report(done, total)`
(travail fait / travail a faire) ou `milestone(f)` (une fraction connue apres une phase indivisible). Rien
n'est estime : une fraction n'est emise que pour du travail effectivement termine. Sans `begin`, tout est
un no-op (tests unitaires, appels hors du pipeline).

  - la fraction est CROISSANTE dans une etape et ne vaut jamais 1.0 ici (le 1.0 est celui de `run_pipeline`) ;
  - debit limite : au plus un message par pas de `MIN_STEP` (1 %) et par `MIN_INTERVAL` (100 ms) ;
  - `phase(lo, hi)` : une sous-etape (ex. step5b dans l'etape « layers ») occupe [lo, hi] de la barre ;
  - sur fils d'execution (step6) : verrou, et la fraction monte quel que soit l'ordre d'arrivee.
"""
from __future__ import annotations

import contextlib
import threading
import time
from typing import Callable, Optional

MIN_STEP = 0.01
MIN_INTERVAL = 0.1
MAX_FRACTION = 0.999

_lock = threading.Lock()
_sink: Optional[Callable[[float], None]] = None
_last_fraction = 0.0
_last_time = 0.0
_scale = (0.0, 1.0)
_clock = time.monotonic


def begin(sink: Optional[Callable[[float], None]]) -> None:
    """Ouvre une etape : `sink(fraction)` recevra des fractions croissantes dans ]0, 1[."""
    global _sink, _last_fraction, _last_time, _scale
    with _lock:
        _sink, _last_fraction, _last_time, _scale = sink, 0.0, -1e9, (0.0, 1.0)


def end() -> None:
    global _sink
    with _lock:
        _sink = None


@contextlib.contextmanager
def phase(lo: float, hi: float):
    """Les appels de `report`/`milestone` du bloc occupent [lo, hi] de la barre de l'etape."""
    global _scale
    with _lock:
        previous, _scale = _scale, (lo, hi)
    try:
        yield
    finally:
        with _lock:
            _scale = previous


def milestone(fraction: float) -> None:
    """Fraction (0..1 de la phase courante) atteinte apres un travail indivisible, emise sans limite de debit
    mais toujours croissante."""
    _emit(fraction, throttle=False)


def report(done: float, total: float, span: Optional[tuple] = None) -> None:
    """`done` unites de travail terminees sur `total` ; limite de debit. `span=(lo, hi)` : ce travail occupe
    seulement [lo, hi] de la phase courante (boucle dont le nombre de tours n'est pas connu d'avance)."""
    if total > 0:
        fraction = done / total
        _emit(span[0] + (span[1] - span[0]) * fraction if span else fraction, throttle=True)


def _emit(fraction: float, throttle: bool) -> None:
    global _last_fraction, _last_time
    with _lock:
        if _sink is None:
            return
        lo, hi = _scale
        value = min(MAX_FRACTION, lo + (hi - lo) * min(max(fraction, 0.0), 1.0))
        if value <= _last_fraction:
            return
        now = _clock()
        if throttle and (value - _last_fraction < MIN_STEP or now - _last_time < MIN_INTERVAL):
            return
        _last_fraction, _last_time = value, now
        sink = _sink
    sink(value)
