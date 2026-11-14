"""R-030 point 3 : la bibliotheque ne lit pas l'environnement et n'a aucun `global`."""
import re

import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LIBRARY_DIRS = ["core", "element_ardor", "element_aeolis", "file_io", "config"]
ALLOWED = {ROOT / "core" / "runtime.py"}   # seul lecteur de l'environnement


@pytest.mark.xfail(strict=True, reason="dette connue : APEX_CPP_LIB_DIR lu dans core/_cpp.py (l.38-39) "
                   "et element_ardor/step1_tetgen_run.py (l.199) ; a retirer quand ces lectures passent par RuntimeOptions")
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



def test_library_has_no_global_statement():
    """R-030 point 3 : aucun `global` mutable dans la bibliotheque."""
    offenders = []
    for d in ["core", "element_ardor", "element_aeolis", "file_io", "config"]:
        for path in (ROOT / d).rglob("*.py"):
            for number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
                if re.match(r"\s*global\s+\w", line):
                    offenders.append(f"{path.relative_to(ROOT)}:{number}")
    assert not offenders, offenders
