"""R-030 point 4 : aucune constante de machine dans le code, tout vient d'un profil de donnees."""
import re
from pathlib import Path

import numpy as np
import pytest

from config.machine_profiles import MACHINES_DIR, default_machine_id, list_machines, load_machine_profile
from core.mesh_types import MachineProfile

ROOT = Path(__file__).resolve().parent.parent
CODE = [ROOT / "main.py", ROOT / "slice_cli.py"] + [
    p for d in ("core", "element-ardor", "element-aeolis", "file_io", "config") for p in (ROOT / d).rglob("*.py")]


def test_every_profile_file_loads_and_is_listed():
    machines = list_machines()
    assert machines and default_machine_id() in machines
    for machine_id in machines:
        profile = load_machine_profile(machine_id)
        assert isinstance(profile, MachineProfile) and profile.toolhead_id and profile.platform_radius > 0


def test_aetheris_4_profile_values():
    p = load_machine_profile("aetheris-4")
    assert p.platform_radius == 150.0 and p.toolhead_id == "aetheris-4" and p.position_coupled_axes == [0]
    assert p.max_angular_velocity == float(np.radians(1666.7))
    assert p.max_angular_acceleration == float(np.radians(1300.0))
    assert [a.angle_range_deg for a in p.rotation_axes] == [(0.0, 360.0), (0.0, 185.0)]
    assert np.array_equal(p.radial_direction, [0.0, -1.0, 0.0])


def test_unknown_and_invalid_ids_are_refused():
    with pytest.raises(FileNotFoundError):
        load_machine_profile("n-existe-pas")
    for bad in ("", "../x", "index", "a/b"):
        with pytest.raises(ValueError):
            load_machine_profile(bad)


def test_class_default_hides_no_machine():
    assert MachineProfile().toolhead_id == "" and MachineProfile().platform_radius is None


def test_code_holds_no_machine_constant():
    literal = re.compile(r"""["']aetheris-4["']|\bPLATFORM_RADIUS\b|platform_radius\s*=\s*[0-9]""")
    offenders = []
    for path in CODE:
        if path.name == "machine_profiles.py" and path.parent.name == "config":
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if literal.search(line) and not line.lstrip().startswith("#"):
                offenders.append(f"{path.relative_to(ROOT)}:{number}")
    assert not offenders, offenders


def test_machines_dir_is_data_only():
    assert {p.suffix for p in Path(MACHINES_DIR).iterdir()} <= {".json", ".md"}
