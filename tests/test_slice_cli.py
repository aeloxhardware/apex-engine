"""
Protocole d'appel du tranchage (slice_cli, RULES du depot apex §5 ; contrats C2 et C3).
"""
import json
import os
import types

import pytest
import trimesh

import slice_cli
from config import schema


def _request(tmp_path, **overrides):
    data = {"contract": schema.APEX_CONTRACT_VERSION, "mesh": "m.stl", "output_dir": str(tmp_path)}
    data.update(overrides)
    return json.dumps(data)


def test_steps_are_unique_and_end_with_export():
    steps = slice_cli._steps()
    slugs = [s["slug"] for s in steps]
    assert len(slugs) == len(set(slugs)) and slugs[-1] == "export"
    assert all(s["label"] for s in steps)


def test_steps_command_prints_contract_and_steps(capsys):
    assert slice_cli.steps_main([]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["contract"] == schema.APEX_CONTRACT_VERSION and out["steps"] == slice_cli._steps()


def test_schema_carries_the_contract_number():
    from scripts import gen_params
    assert json.loads(gen_params.build_schema())["x-apex-contract"] == schema.APEX_CONTRACT_VERSION


@pytest.mark.parametrize("text, fragment", [
    ("pas du json", "illisible"),
    ("[]", "objet JSON"),
    (json.dumps({"contract": schema.APEX_CONTRACT_VERSION, "mesh": "m", "output_dir": "o", "extra": 1}),
     "inconnu"),
    (json.dumps({"mesh": "m", "output_dir": "o"}), "obligatoire"),
    (json.dumps({"contract": schema.APEX_CONTRACT_VERSION + 1, "mesh": "m", "output_dir": "o"}), "contrat"),
    (json.dumps({"contract": True, "mesh": "m", "output_dir": "o"}), "contrat"),
    (json.dumps({"contract": schema.APEX_CONTRACT_VERSION, "mesh": "m", "output_dir": "o",
                 "output_name": "dossier/x.gcode"}), "output_name"),
    (json.dumps({"contract": schema.APEX_CONTRACT_VERSION, "mesh": "m", "output_dir": "o",
                 "params": {"resolutoin": 1.0}}), "resolutoin"),
    (json.dumps({"contract": schema.APEX_CONTRACT_VERSION, "mesh": "m", "output_dir": "o",
                 "gcode_params": {"nope": 1}}), "nope"),
])
def test_invalid_requests_are_refused(text, fragment):
    with pytest.raises(slice_cli.InvalidRequest, match=fragment):
        slice_cli.parse_request(text)


def test_partial_params_keep_defaults_and_gcode_is_strict_by_default(tmp_path):
    from main import PipelineParams
    req = slice_cli.parse_request(_request(tmp_path, params={"resolution": 2.5}, gcode_params={"print_speed": 20}))
    assert req["params"].resolution == 2.5
    assert req["params"].layers.layer_height_base == PipelineParams().layers.layer_height_base
    assert req["gcode_params"].print_speed == 20 and req["gcode_params"].strict is True


@pytest.mark.parametrize("tip_below, depth, turns, limit, expected", [
    (1, -0.05, 3.0, None, "pour 1 mouvement (jusqu'à 0,05 mm"),
    (4, -0.17, 3.0, None, "pour 4 mouvements (jusqu'à 0,17 mm"),
    (0, 1.0, 12.5, 10.0, "(12,5 tours) dépasse la limite (10,0 tours)"),
])
def test_export_refusal_message_agrees_in_number(tip_below, depth, turns, limit, expected):
    report = types.SimpleNamespace(tip_below_plate=tip_below, tip_z_min=depth, y_turns=turns)
    message = slice_cli.export_refusal_message(report, types.SimpleNamespace(max_y_turns=limit))
    assert expected in message
    assert "1 mouvements" not in message


def _box_stl(tmp_path):
    path = tmp_path / "cube.stl"
    trimesh.creation.box(extents=(8, 8, 8)).export(path)
    return str(path)


def _events(stderr: str):
    out = []
    for line in stderr.splitlines():
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue   # journal libre du pipeline
        if isinstance(event, dict) and event.get("type") in ("step", "log"):
            out.append(event)
    return out


def test_slice_end_to_end(tmp_path, capsys):
    request = tmp_path / "req.json"
    out_dir = tmp_path / "sortie"
    request.write_text(_request(tmp_path, mesh=_box_stl(tmp_path), output_dir=str(out_dir),
                                params={"resolution": 3.0}))
    assert slice_cli.slice_main(["--request", str(request)]) == 0
    captured = capsys.readouterr()

    # sortie standard : une seule ligne, le resultat
    lines = captured.out.strip().splitlines()
    assert len(lines) == 1
    result = json.loads(lines[0])
    assert result["status"] == "ok" and result["contract"] == schema.APEX_CONTRACT_VERSION
    assert os.path.isfile(result["gcode"]) and result["report"]["moves"] > 0
    assert result["toolpath_points"] > 0 and "placement_offset_z" in result

    # progression : chaque etape annoncee, dans l'ordre, de 0 a 1, et seulement celles-la
    steps = [e for e in _events(captured.err) if e["type"] == "step"]
    expected = [s["slug"] for s in slice_cli._steps()]
    assert [e["slug"] for e in steps if e["fraction"] == 0.0] == expected
    assert [e["slug"] for e in steps if e["fraction"] == 1.0] == expected
    assert all(0.0 <= e["fraction"] <= 1.0 for e in steps)

    # aucun fichier temporaire ne reste
    assert sorted(os.listdir(out_dir)) == [os.path.basename(result["gcode"])]


def test_slice_missing_mesh_is_invalid_request(tmp_path, capsys):
    request = tmp_path / "req.json"
    request.write_text(_request(tmp_path, mesh=str(tmp_path / "absent.stl")))
    assert slice_cli.slice_main(["--request", str(request)]) == 2
    captured = capsys.readouterr()
    assert json.loads(captured.out.strip().splitlines()[-1])["status"] == "invalid_request"
    assert any(e["type"] == "log" and e["level"] == "error" for e in _events(captured.err))


def test_slice_unreadable_model_is_refused_with_exit_3(tmp_path, capsys):
    bad = tmp_path / "vide.stl"
    bad.write_text("ceci n'est pas un maillage")
    request = tmp_path / "req.json"
    request.write_text(_request(tmp_path, mesh=str(bad)))
    assert slice_cli.slice_main(["--request", str(request)]) == 3
    result = json.loads(capsys.readouterr().out.strip().splitlines()[-1])
    assert result["status"] in ("model_refused", "surface_inconsistent") and result["message"]
    assert not [f for f in os.listdir(tmp_path) if f.endswith((".gcode", ".part"))]
