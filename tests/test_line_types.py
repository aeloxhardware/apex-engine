"""Types de ligne des trajectoires (`core/line_types.py`) : étape 6 -> 7 -> 8 -> étiquettes `;TYPE:` du G-code.

Un type par polyligne (paroi extérieure, paroi intérieure, remplissage, pont), porté jusqu'au G-code pour que
l'aperçu colore les trajectoires comme Element Aeolis. Chaque test qui touche au natif est joué dans les deux positions
du commutateur (Python, natif) ; le natif absent ou périmé fait sauter ses tests, jamais échouer.
"""
import dataclasses

import numpy as np
import pytest

from core import _cpp, line_types
from core.mesh_types import Layer, PrintPass, PrintSchedule


@pytest.fixture
def restore_switch():
    previous = _cpp.set_enabled(None)
    yield
    _cpp.set_enabled(previous)


def _native_has(name: str) -> bool:
    return _cpp.available() and bool(getattr(_cpp.module(), name, False))


# --- héritage à l'étape 7 -------------------------------------------------------------------------------------


def _line(n: int, y: float = 0.0) -> np.ndarray:
    return np.column_stack([np.arange(n, dtype=float), np.full(n, y), np.zeros(n)])


def test_inherit_kinds_identity_without_split():
    old = [_line(4), _line(3), _line(5)]
    kinds = [line_types.INFILL, line_types.WALL_INNER, line_types.WALL_OUTER]
    assert line_types.inherit_kinds(old, kinds, [p.copy() for p in old]) == kinds


def test_inherit_kinds_follows_contiguous_slices():
    """Découpe par couleur : `poly[a:b + 1]`, tranches contiguës dont les arêtes recouvrent la parente."""
    old = [_line(7), _line(4)]
    kinds = [line_types.WALL_OUTER, line_types.INFILL]
    new = [old[0][0:3], old[0][2:7], old[1][0:4]]                 # 2 + 4 arêtes, puis 3
    assert line_types.inherit_kinds(old, kinds, new) == [line_types.WALL_OUTER, line_types.WALL_OUTER,
                                                         line_types.INFILL]


def test_inherit_kinds_skips_polylines_dropped_by_the_split():
    """Une polyligne d'un point n'a ni arête ni tranche (rejetée par l'étape 7) : elle ne décale rien."""
    old = [_line(1), _line(4)]
    kinds = [line_types.INFILL, line_types.WALL_INNER]
    new = [old[1][0:2], old[1][1:4]]
    assert line_types.inherit_kinds(old, kinds, new) == [line_types.WALL_INNER, line_types.WALL_INNER]


def test_inherit_kinds_never_guesses():
    """Découpage non rattachable ou types absents : tout est `NONE`, jamais un type faux."""
    old = [_line(5)]
    assert line_types.inherit_kinds(old, [line_types.INFILL], [_line(9)]) == [line_types.NONE]
    assert line_types.inherit_kinds(old, None, [_line(5)]) == [line_types.NONE]


# --- étape 6 -----------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def box_layers():
    import trimesh
    from test_step6_toolpath2d import _pipeline_layers_from_trimesh
    return _pipeline_layers_from_trimesh(trimesh.creation.box(extents=[10, 10, 10]), platform_z=-5.0,
                                         layer_height_base=2.0)


def _fill(box_layers, native: bool, **overrides):
    from element_ardor import step6_toolpath2d as s6
    mesh_data, overhang, layers = box_layers
    params = s6.Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.3, native_layer=native, **overrides)
    return s6.run(layers, mesh_data, overhang, params)


@pytest.mark.parametrize("split", [False, True])
def test_step6_python_kinds_parallel_to_polylines(box_layers, restore_switch, split):
    _cpp.set_enabled(False)
    # ponts desactives : la premiere couche de la boite est vue comme un pont (remplissage typé BRIDGE, test dédié plus bas)
    out = _fill(box_layers, native=False, split_rings_into_runs=split, continuous_rings=split,
                bridge_contour_ratio_threshold=2.0)
    assert out
    for layer in out:
        assert layer.toolpath_kinds is not None and len(layer.toolpath_kinds) == len(layer.toolpath_2d)
        kinds = set(layer.toolpath_kinds)
        assert line_types.WALL_OUTER in kinds and line_types.INFILL in kinds     # garde : il y a bien un remplissage
        assert line_types.WALL_INNER in kinds                                     # wall_count = 2
        first_wall = min(i for i, k in enumerate(layer.toolpath_kinds) if k != line_types.INFILL)
        assert all(k == line_types.INFILL for k in layer.toolpath_kinds[:first_wall])   # remplissage d'abord
        assert layer.toolpath_kinds[-1] == line_types.WALL_OUTER                  # contour d'origine en dernier


@pytest.mark.skipif(not _native_has("FILL_LAYERS_KINDS"), reason="module natif sans types de ligne")
@pytest.mark.parametrize("split", [False, True])
def test_step6_native_kinds_match_python(box_layers, restore_switch, split):
    """Mêmes types, dans les deux positions du commutateur (couches entières natives ou couche par couche)."""
    options = dict(split_rings_into_runs=split, continuous_rings=split, bridge_contour_ratio_threshold=2.0)
    _cpp.set_enabled(False)
    ref = _fill(box_layers, native=False, **options)
    _cpp.set_enabled(True)
    nat = _fill(box_layers, native=True, **options)
    for a, b in zip(ref, nat):
        assert len(a.toolpath_2d) == len(b.toolpath_2d)
        assert list(a.toolpath_kinds) == list(b.toolpath_kinds)


def test_step6_bridge_infill_is_typed_bridge(restore_switch):
    import trimesh
    from test_step6_toolpath2d import _pipeline_layers_from_trimesh
    from element_ardor import step6_toolpath2d as s6
    pillar = trimesh.creation.box(extents=[3, 3, 10])
    a, b = pillar.copy(), pillar.copy()
    a.apply_translation([-6, 0, 5])
    b.apply_translation([6, 0, 5])
    slab = trimesh.creation.box(extents=[16, 3, 1.5])
    slab.apply_translation([0, 0, 9])
    mesh = a.union(b, engine="manifold").union(slab, engine="manifold")
    mesh_data, overhang, layers = _pipeline_layers_from_trimesh(mesh, platform_z=0.0, bridge_max_gap_distance=4.0)
    _cpp.set_enabled(False)
    # une paroi seulement, sans renfort d'ancrage : avec 4 parois de 0,4 mm, un pilier de 3 mm n'a plus de remplissage
    params = s6.Step6Params(native_layer=False, wall_count=1, bridge_anchor_extra_walls=0, infill_density=0.3)
    out = s6.run(layers, mesh_data, overhang, params)
    kinds = {k for layer in out for k in layer.toolpath_kinds}
    assert line_types.BRIDGE in kinds, "la dalle doit produire un remplissage de pont"
    assert line_types.INFILL in kinds          # les piliers, eux, restent en remplissage ordinaire


# --- étape 8 ------------------------------------------------------------------------------------------------------


def test_step8_points_carry_the_line_type_of_their_polyline():
    from element_ardor import step8_machine_toolpath as s8
    from core.mesh_types import ToolpathPoint
    frame = (np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
    polys = [_line(3), _line(4, 1.0)]
    layer = Layer(id=7, contour_3d=[], region_ids=[], local_frame=frame, toolpath_2d=polys,
                  toolpath_kinds=[line_types.INFILL, line_types.WALL_OUTER])
    schedule = PrintSchedule(passes=[PrintPass(segments=[(7, 0), (7, 1)], color_id=None, is_hole_fill=False,
                                               layer_range=(0, 0), deferred_regions=[])])
    raw = s8._iter_schedule_points(schedule, {7: layer}, s8.Step8Params())
    assert raw.is_travel.sum() == 1                                      # un trajet à vide entre les deux polylignes
    points = [ToolpathPoint(position=p, axis_angles=[0.0, 0.0], extrusion=float(e), feedrate=1.0)
              for p, e in zip(raw.position, raw.extrusion)]
    s8.tag_line_types(points, raw, {7: layer})
    assert [p.line_type for p in points] == [line_types.INFILL] * 3 + [line_types.NONE] + [line_types.WALL_OUTER] * 4


def test_step8_untyped_layer_stays_untyped():
    from element_ardor import step8_machine_toolpath as s8
    from core.mesh_types import ToolpathPoint
    frame = (np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
    layer = Layer(id=1, contour_3d=[], region_ids=[], local_frame=frame, toolpath_2d=[_line(3)])
    schedule = PrintSchedule(passes=[PrintPass(segments=[(1, 0)], color_id=None, is_hole_fill=False,
                                               layer_range=(0, 0), deferred_regions=[])])
    raw = s8._iter_schedule_points(schedule, {1: layer}, s8.Step8Params())
    points = [ToolpathPoint(position=p, axis_angles=[0.0, 0.0], extrusion=1.0, feedrate=1.0) for p in raw.position]
    s8.tag_line_types(points, raw, {1: layer})
    assert {p.line_type for p in points} == {line_types.NONE}


# --- G-code -------------------------------------------------------------------------------------------------------


def _typed_toolpath(seed: int = 0):
    from test_cpp_equivalence_mesh import _random_machine_toolpath
    tp = _random_machine_toolpath(seed)
    n = len(tp["position"])
    rng = np.random.default_rng(seed)
    runs = np.repeat(rng.integers(0, 6, size=n // 20 + 1), 20)[:n]      # types par séries de 20 points, 0 compris
    tp["line_type"] = runs.astype(np.uint8)
    return tp


def _export(tp, native: bool):
    import main
    from file_io.gcode_export import GcodeParams, export_gcode
    machine = main._default_machine_profile()
    machine.nozzle_model = None
    _cpp.set_enabled(native)
    return export_gcode(tp, machine, GcodeParams(max_tip_deviation=None))[0]


def test_gcode_labels_precede_the_first_printed_move_of_each_type(restore_switch):
    tp = _typed_toolpath()
    out = _export(tp, native=False)
    labels = [line for line in out if line.startswith(";TYPE:")]
    assert labels, "aucune étiquette écrite"
    assert all(a != b for a, b in zip(labels, labels[1:]))              # une étiquette seulement au changement de type
    assert {l[len(";TYPE:"):] for l in labels} <= set(line_types.LABELS.values())
    # l'étiquette précède immédiatement un mouvement (jamais en fin de fichier ni avant un trajet)
    for i, line in enumerate(out):
        if line.startswith(";TYPE:"):
            nxt = next(l for l in out[i + 1:] if not l.startswith(";"))
            assert nxt.startswith(("G00", "G01"))


def test_gcode_without_line_types_is_unchanged(restore_switch):
    tp = _typed_toolpath()
    typed = _export(tp, native=False)
    del tp["line_type"]
    plain = _export(tp, native=False)
    assert not any(line.startswith(";TYPE:") for line in plain)
    assert [l for l in typed if not l.startswith(";TYPE:")] == plain    # les étiquettes sont les seules différences


@pytest.mark.skipif(not _native_has("EXPORT_LINE_TYPES"), reason="module natif sans types de ligne")
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_gcode_labels_native_identical_to_python(restore_switch, seed):
    tp = _typed_toolpath(seed)
    ref = _export(tp, native=False)
    nat = _export(tp, native=True)
    assert ref == nat
    assert sum(l.startswith(";TYPE:") for l in ref) > 5                  # garde anti-faux-vert


def test_parse_gcode_ignores_type_labels(restore_switch):
    from file_io.gcode_export import parse_gcode
    tp = _typed_toolpath()
    typed, plain = _export(tp, native=False), None
    del tp["line_type"]
    plain = _export(tp, native=False)
    assert parse_gcode(typed) == parse_gcode(plain)


def test_toolpath_arrays_exposes_line_types():
    from core.mesh_types import ToolpathMachine, ToolpathPoint
    from file_io.gcode_export import toolpath_arrays
    pts = [ToolpathPoint(position=np.zeros(3), axis_angles=[0.0, 0.0], extrusion=1.0, feedrate=1.0, line_type=t)
           for t in (0, 3, 1)]
    arrays = toolpath_arrays(ToolpathMachine(points=pts))
    assert arrays["line_type"].dtype == np.uint8 and arrays["line_type"].tolist() == [0, 3, 1]


def test_skin_points_are_typed_skin():
    from element_ardor.skin import inject
    from core.mesh_types import ToolpathMachine, ToolpathPoint
    frame = (np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
    layer = Layer(id=0, contour_3d=[], region_ids=[], local_frame=frame)
    body = ToolpathMachine(points=[ToolpathPoint(position=np.zeros(3), axis_angles=[0.0, 0.0], extrusion=1.0,
                                                 feedrate=1.0, line_type=line_types.INFILL)],
                           layer_ids=np.array([0]))
    passes = [{"kind": "body", "layer": 0},
              {"kind": "skin", "points": [np.array([[0.0, 0, 1], [1.0, 0, 1]])], "angles": [np.zeros((2, 2))]}]
    out, _ = inject.inject(body, [layer], passes, feedrate=10.0)
    assert [p.line_type for p in out.points] == [line_types.INFILL, line_types.SKIN, line_types.SKIN]
    assert dataclasses.is_dataclass(out.points[0])
