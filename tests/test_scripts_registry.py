"""Registre des scripts (limite dure n° 9) : tout script de `scripts/` est listé dans `scripts/README.md` et compile."""
import py_compile
import re
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def _listed() -> set[str]:
    return set(re.findall(r"^\| `([\w.]+)` \|", (SCRIPTS / "README.md").read_text(encoding="utf-8"), re.M))


def test_every_script_is_listed_and_every_listed_script_exists():
    present = {p.name for p in SCRIPTS.iterdir() if p.suffix in (".py", ".json")}
    listed = _listed()
    assert present - listed == set(), f"scripts non listés dans scripts/README.md : {sorted(present - listed)}"
    assert listed - present == set(), f"scripts listés mais absents : {sorted(listed - present)}"


def test_every_script_compiles(tmp_path):
    for path in sorted(SCRIPTS.glob("*.py")):
        py_compile.compile(str(path), cfile=str(tmp_path / (path.name + "c")), doraise=True)
