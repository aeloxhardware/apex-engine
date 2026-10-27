"""R-020 : la tete se designe par un identifiant, jamais par un chemin."""
import pytest

import main
from core.mesh_types import MachineProfile
from file_io.mesh_loader import load_toolhead_by_id, toolhead_path


def test_toolhead_path_rejects_anything_but_an_identifier():
    for bad in ("../x", "a/b", "a\\b", "", "a.glb", "a b"):
        with pytest.raises(ValueError):
            toolhead_path(bad)


def test_unknown_toolhead_is_a_file_not_found():
    with pytest.raises(FileNotFoundError):
        load_toolhead_by_id("n-existe-pas")


def test_ensure_toolhead_loads_the_designated_head_and_never_raises(capsys):
    profile = MachineProfile(toolhead_id="aetheris-4")
    main.ensure_toolhead(profile)
    assert profile.nozzle_model is not None and len(profile.nozzle_model.vertices) > 0
    unknown = MachineProfile(toolhead_id="n-existe-pas")
    main.ensure_toolhead(unknown)
    assert unknown.nozzle_model is None and "n-existe-pas" in capsys.readouterr().err


def test_ensure_toolhead_keeps_a_head_already_loaded():
    sentinel = object()
    profile = MachineProfile(nozzle_model=sentinel)
    main.ensure_toolhead(profile)
    assert profile.nozzle_model is sentinel
