"""Contrat UI -> moteur : aller-retour JSON des parametres et schema derive
des annotations de type (config/schema.py)."""
import dataclasses
import json

import numpy as np
import pytest

from config import schema
from config.params import ROOTS, GcodeParams, PipelineParams
from core.mesh_types import RotationAxisSpec


def _eq(a, b):
    if dataclasses.is_dataclass(a):
        return type(a) is type(b) and all(
            _eq(getattr(a, f.name), getattr(b, f.name)) for f in dataclasses.fields(a) if f.name != "nozzle_model")
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return np.array_equal(np.asarray(a, dtype=float), np.asarray(b, dtype=float))
    if isinstance(a, (list, tuple)):
        return type(a) is type(b) and len(a) == len(b) and all(_eq(x, y) for x, y in zip(a, b))
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(_eq(a[k], b[k]) for k in a)
    return a == b


@pytest.mark.parametrize("name", sorted(ROOTS))
def test_round_trip_defaults(name):
    cls = ROOTS[name]
    obj = cls()
    text = schema.dumps(obj)
    back = schema.loads(cls, text)
    assert _eq(obj, back)
    assert json.loads(schema.dumps(back)) == json.loads(text)


def test_round_trip_modified_values():
    p = PipelineParams()
    p.resolution = 0.8
    p.segmentation.max_region_merge_angle_deg = 30.0
    p.toolpath2d.infill_pattern = type(p.toolpath2d.infill_pattern)("lines") \
        if "lines" in [m.value for m in type(p.toolpath2d.infill_pattern)] else p.toolpath2d.infill_pattern
    p.coupled_reorientation_tiers = (45.0, 90.0, 181.0)
    p.layers.support_rotation_axes = [RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0))]
    back = schema.loads(PipelineParams, schema.dumps(p))
    assert _eq(p, back)


def test_partial_update_keeps_defaults_and_rejects_typos():
    back = schema.from_jsonable(PipelineParams, {"layers": {"layer_height_base": 0.3}})
    assert back.layers.layer_height_base == 0.3
    assert back.layers.weld_tolerance_factor == PipelineParams().layers.weld_tolerance_factor
    assert back.resolution == PipelineParams().resolution
    with pytest.raises(ValueError):
        schema.from_jsonable(PipelineParams, {"layers": {"layer_heigth_base": 0.3}})
    with pytest.raises(ValueError):
        schema.from_jsonable(PipelineParams, {"machine_toolpath": {"free_axis_median_window": 2.5}})


def test_nozzle_model_is_runtime_only():
    s = schema.json_schema(PipelineParams)
    mp = s["$defs"]["MachineProfile"]["properties"]["nozzle_model"]
    assert mp.get("x-runtime-only") is True
    assert "nozzle_model" not in json.loads(schema.dumps(PipelineParams()))["machine_profile"]


@pytest.mark.parametrize("name", sorted(ROOTS))
def test_every_field_in_schema(name):
    """Chaque champ de chaque dataclass de l'arbre figure dans le schema, avec un type."""
    s = schema.json_schema(ROOTS[name])
    for path, _t in schema.walk_fields(ROOTS[name]):
        node = s
        for part in path.split("."):
            prop = node["properties"][part]
            node = prop
            ref = prop.get("$ref") or next((a.get("$ref") for a in prop.get("anyOf", []) if "$ref" in a), None)
            if ref is None:
                items = prop.get("items") or next((a.get("items") for a in prop.get("anyOf", []) if "items" in a), None)
                ref = items.get("$ref") if isinstance(items, dict) else None
            if ref:
                node = s["$defs"][ref.split("/")[-1]]
        assert prop, path


def test_schema_defaults_are_json():
    for cls in ROOTS.values():
        json.dumps(schema.json_schema(cls))


def test_generated_files_up_to_date():
    """config/params.schema.json et cpp/include/apex/params.hpp sont regeneres
    depuis les dataclasses : un champ ajoute sans relancer le generateur echoue ici."""
    import subprocess
    import sys
    r = subprocess.run([sys.executable, "scripts/gen_params.py", "--check"], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_native_params_header_matches_python_defaults():
    from core import _cpp
    core = _cpp.module() if _cpp.available() else None
    if core is None or not hasattr(core, "params_defaults_probe"):
        pytest.skip("module natif absent ou anterieur au pas 11")
    d = core.params_defaults_probe()
    p, g = PipelineParams(), GcodeParams()
    assert d["resolution"] == p.resolution
    assert d["near_platform_resolution"] == p.near_platform_resolution
    assert d["layer_height_base"] == p.layers.layer_height_base
    assert d["out_of_plane_interpolation_power"] == p.toolpath2d.out_of_plane_interpolation_power
    assert d["split_rings_into_runs"] == p.toolpath2d.split_rings_into_runs
    assert d["infill_pattern"] == p.toolpath2d.infill_pattern.value
    assert d["heat_method_time_factor"] == p.field_solver.heat_method_time_factor
    assert d["heat_project_reachable"] == p.field_solver.heat_project_reachable
    assert d["rotation_axes_count"] == len(p.machine_profile.rotation_axes)
    assert d["b_max_deg"] == p.machine_profile.rotation_axes[1].angle_range_deg[1]
    assert d["position_coupled_axes"] == list(p.machine_profile.position_coupled_axes)
    assert d["pivot_to_tip"] == g.pivot_to_tip and d["tip_at_bead_top"] == g.tip_at_bead_top
    assert d["z_reference"] == g.z_reference
