"""R-030 point 3 : aucune variable d'environnement ne change le resultat dans la bibliotheque."""
import re
from pathlib import Path

import pytest

from core import _cpp
from core.runtime import RuntimeOptions, options_from_environment

ROOT = Path(__file__).resolve().parent.parent
LIBRARY_DIRS = ["core", "element-ardor", "element-aeolis", "file_io", "config"]
ALLOWED = {ROOT / "core" / "runtime.py"}   # seul lecteur de l'environnement


def test_library_never_reads_the_environment():
    pattern = re.compile(r"os\.environ|os\.getenv|\bgetenv\(")
    offenders = []
    for d in LIBRARY_DIRS:
        for path in (ROOT / d).rglob("*.py"):
            if path in ALLOWED:
                continue
            for number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
                if pattern.search(line) and not line.lstrip().startswith("#"):
                    offenders.append(f"{path.relative_to(ROOT)}:{number}")
    assert not offenders, offenders


def test_options_from_environment_parses_every_variable():
    env = {"APEX_USE_CPP": "1", "APEX_USE_GPU": "yes", "APEX_CPP_GRID_MODE": "Faithful",
           "APEX_CPP_LIB_DIR": "/x", "APEX_TETGEN_WORKER": "0", "APEX_TETGEN_LOG": "/log"}
    o = options_from_environment(env)
    assert (o.use_cpp, o.use_gpu, o.grid_mode, o.lib_dir, o.tetgen_worker, o.tetgen_log) == (
        True, True, "faithful", "/x", False, "/log")
    assert options_from_environment({}) == RuntimeOptions()


def test_unknown_grid_mode_is_refused():
    with pytest.raises(ValueError):
        RuntimeOptions(grid_mode="x")


def test_configure_drives_enabled_grid_mode_and_returns_previous():
    previous_forced = _cpp.set_enabled(None)
    previous = _cpp.configure(RuntimeOptions(use_cpp=False, grid_mode="faithful"))
    try:
        assert _cpp.enabled() is False
        assert _cpp.grid_mode() == _cpp.GRID_MODE_FAITHFUL
        _cpp.configure(RuntimeOptions(use_cpp=True))
        assert _cpp.enabled() is _cpp.available()
        assert _cpp.grid_mode() == _cpp.GRID_MODE_STAMPED
    finally:
        _cpp.configure(previous)
        _cpp.set_enabled(previous_forced)


def test_changing_lib_dir_after_load_is_refused():
    _cpp.available()   # force le chargement
    previous = _cpp.options()
    with pytest.raises(RuntimeError):
        _cpp.configure(RuntimeOptions(lib_dir=str(ROOT / "nonexistent_lib_dir")))
    assert _cpp.options() is previous
