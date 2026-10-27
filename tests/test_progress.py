"""R-017 : fractions intermediaires reelles (core/progress.py) et leur cablage dans le pipeline."""
import numpy as np
import pytest

from core import progress


@pytest.fixture
def sink():
    got = []
    clock = [0.0]
    old = progress._clock
    progress._clock = lambda: clock[0]
    progress.begin(got.append)
    yield got, clock
    progress.end()
    progress._clock = old


def test_report_is_noop_without_a_stage():
    progress.end()
    progress.report(1, 2)          # ne leve rien
    progress.milestone(0.5)


def test_fractions_are_increasing_throttled_and_below_one(sink):
    got, clock = sink
    for done in range(1001):
        clock[0] += 1.0                      # le temps ne limite pas : seul le pas de 1 % compte
        progress.report(done, 1000)
    assert got == sorted(got) and len(set(got)) == len(got)
    assert all(0.0 < f < 1.0 for f in got)
    assert all(b - a >= progress.MIN_STEP - 1e-12 for a, b in zip(got, got[1:]))
    assert len(got) <= 101


def test_throttle_by_time(sink):
    got, clock = sink
    progress.report(50, 100)
    clock[0] += progress.MIN_INTERVAL / 2
    progress.report(90, 100)                 # trop tot
    assert got == [0.5]
    clock[0] += progress.MIN_INTERVAL
    progress.report(90, 100)
    assert got == [0.5, 0.9]


def test_never_goes_backwards(sink):
    got, _ = sink
    progress.milestone(0.6)
    progress.milestone(0.2)
    assert got == [0.6]


def test_phase_scales_into_its_range(sink):
    got, _ = sink
    with progress.phase(0.7, 1.0):
        progress.milestone(0.5)
    assert got == [pytest.approx(0.85)]
    progress.milestone(0.9)                  # hors phase : echelle pleine
    assert got[-1] == pytest.approx(0.9)


def test_span_for_loops_of_unknown_length(sink):
    got, _ = sink
    progress.report(1, 2, span=(0.0, 0.8))
    assert got == [pytest.approx(0.4)]


def test_native_step6_chunks_match_a_single_call(monkeypatch):
    """Le remplissage natif par tranches rend les memes polylignes qu'un appel unique."""
    from dataclasses import replace
    from core.mesh_types import Layer
    from pipeline import step6_toolpath2d as s6

    layers = []
    for k in range(70):
        ring = np.array([[0, 0, k], [10, 0, k], [10, 10, k], [0, 10, k]], dtype=float) + k * 0.01
        layers.append(Layer(id=k, contour_3d=[ring], region_ids=[0],
                            local_frame=(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0])),
                            isovalue=float(k)))
    calls = []

    class FakeModule:
        @staticmethod
        def fill_layers(frames, loop_begin, point_begin, points, sweep, config, kdtree, is_bridge, span, threads):
            calls.append(len(frames))
            assert loop_begin[0] == 0 and point_begin[0] == 0
            assert loop_begin[-1] == len(point_begin) - 1 and point_begin[-1] == len(points)
            return [[np.array([[sweep[i], float(points[point_begin[loop_begin[i]], 2]), 0.0]])] for i in range(len(frames))]

    monkeypatch.setattr(s6, "_native_fill_module", lambda: FakeModule)
    monkeypatch.setattr(s6, "_NATIVE_LAYER", True)
    mesh_data = type("M", (), {"kdtree": type("K", (), {"_native": object()})()})()
    overhang = type("O", (), {"is_bridge": np.zeros(3, bool), "bridge_span_direction": np.zeros((3, 3))})()
    params = s6.Step6Params()
    rank = {l.isovalue: r for r, l in enumerate(layers)}
    out = s6._fill_layers_native(layers, mesh_data, overhang, params, rank)
    assert sum(calls) == 70 and len(calls) > 1 and all(c >= s6._NATIVE_MIN_CHUNK or c == calls[-1] for c in calls)
    assert [float(l.toolpath_2d[0][0, 1]) for l in out] == [pytest.approx(float(k) + k * 0.01) for k in range(70)]


def test_failed_slice_leaves_no_stage_open(tmp_path, capsys):
    """Un echec en cours d'etape ne laisse pas de rappel vers le flux de journal, qui va etre ferme."""
    import json
    import slice_cli
    from config import schema
    bad = tmp_path / "vide.stl"
    bad.write_text("ceci n'est pas un maillage")
    request = tmp_path / "req.json"
    request.write_text(json.dumps({"contract": schema.APEX_CONTRACT_VERSION, "mesh": str(bad), "output_dir": str(tmp_path)}))
    assert slice_cli.slice_main(["--request", str(request)]) == 3
    assert progress._sink is None
