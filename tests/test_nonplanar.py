"""Tranchage non planaire (paquet `nonplanar/`) : maillages synthétiques générés à la volée."""
from pathlib import Path

import numpy as np
import pytest
import trimesh

from nonplanar import NotAHeightfieldError, SliceParams, build, from_mesh, slice_mesh

MODELS = Path(__file__).resolve().parent.parent / "res" / "test_models" / "non_planar"


def wave_plate(size=20.0, n=41, base=0.0, amp=1.5, mean=3.0) -> trimesh.Trimesh:
    """Plaque dont le dessus est z = mean + amp·sin(2πx/size·2) ; dessous plat ; flancs verticaux."""
    xs = np.linspace(0, size, n)
    X, Y = np.meshgrid(xs, xs)
    Z = base + mean + amp * np.sin(2 * np.pi * X / size * 2)
    top_v = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    bot_v = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, base)])
    idx = np.arange(X.size).reshape(n, n)
    quads = np.array([[idx[j, i], idx[j, i + 1], idx[j + 1, i + 1], idx[j + 1, i]]
                      for j in range(n - 1) for i in range(n - 1)])
    top_f = np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
    verts = np.vstack([top_v, bot_v])
    off = X.size
    bot_f = np.vstack([quads[:, [0, 2, 1]], quads[:, [0, 3, 2]]]) + off
    faces = [top_f, bot_f]
    # flancs : quatre bords
    def wall(a, b):
        out = []
        for p, q in zip(a[:-1], a[1:]):
            out += [[p, q, q + off], [p, q + off, p + off]]
        return np.array(out)
    edge_x0 = idx[:, 0][::-1]
    edge_x1 = idx[:, -1]
    edge_y0 = idx[0, :]
    edge_y1 = idx[-1, :][::-1]
    for e in (edge_y0, edge_x1, edge_y1, edge_x0):
        faces.append(wall(e, None))
    m = trimesh.Trimesh(verts, np.vstack(faces), process=True)
    m.fix_normals()
    return m


def test_heightfield_matches_analytic_surface():
    m = wave_plate()
    f = from_mesh(m, cell=0.25)
    xx, yy = np.meshgrid(f.xs, f.ys)
    expected = 3.0 + 1.5 * np.sin(2 * np.pi * xx / 20 * 2)
    assert np.nanmax(np.abs(f.top - expected)) < 0.05          # facettes planes
    assert np.allclose(f.bottom, 0.0)


@pytest.mark.parametrize("inter", [False, True])
def test_layers_nested_end_on_top_and_bounded_thickness(inter):
    f = from_mesh(wave_plate(), cell=0.25)
    stack = build(f, intermediate=inter)
    xx, yy = np.meshgrid(f.xs[::6], f.ys[::6])
    x, y = xx.ravel(), yy.ravel()
    prev = f.bottom_at(x, y)
    for k, j in stack.sequence():
        z = stack.z(k, x, y) if j == 0 else stack.z_sub(k, j, x, y)
        t = stack.thickness(k, x, y) if j == 0 else stack.thickness_sub(k, j, x, y)
        active = stack.active_mask(k, j)[::6, ::6].ravel() & (t > 1e-6)
        assert (z >= prev - 1e-9)[active].all()                 # couches emboîtées
        assert (t[active] <= 0.25 + 1e-9).all()                 # plafond
        assert (t[active] >= 0.1 - 1e-9).all()                  # plancher
        prev = np.where(active, z, prev)
    assert np.allclose(prev, f.top_at(x, y), atol=1e-6)         # la dernière couche EST le dessus


def test_intermediate_layers_cover_only_thick_columns():
    f = from_mesh(wave_plate(amp=2.5, mean=3.5), cell=0.25)
    stack = build(f, intermediate=True)
    subs = [(k, j) for k, j in stack.sequence() if j > 0]
    assert subs, "colonnes hautes : il faut des inter-couches"
    full = f.top.size - np.isnan(f.top).sum()
    for k, j in subs:
        assert stack.active_mask(k, j).sum() < full             # jamais toute la pièce
    plain = build(f, intermediate=False)
    assert all(j == 0 for _, j in plain.sequence())


def test_overhang_is_refused():
    pillar = trimesh.creation.box([2, 2, 5])
    pillar.apply_translation([0, 0, 2.5])
    cap = trimesh.creation.box([10, 10, 1])
    cap.apply_translation([0, 0, 5.5])
    mushroom = trimesh.util.concatenate([pillar, cap])
    with pytest.raises(NotAHeightfieldError):
        from_mesh(mushroom, cell=0.5)


def test_gcode_extrusion_monotone_and_z_inside_part():
    m = wave_plate()
    gcode, rep = slice_mesh(m, SliceParams(grid_cell=0.5))
    assert rep.layer_count > 5
    f = from_mesh(m, cell=0.5)
    last_e, zmax = 0.0, 0.0
    for line in gcode.splitlines():
        if line.startswith("G1") and " E" in line and " X" in line:
            d = {t[0]: float(t[1:]) for t in line.split(";")[0].split()[1:]}
            assert d["E"] >= last_e - 1e-9
            last_e = d["E"]
            zmax = max(zmax, d["Z"])
            assert f.bottom_at(d["X"], d["Y"]) - 2e-3 <= d["Z"] <= f.top_at(d["X"], d["Y"]) + 2e-3  # G-code à 3 décimales
    assert zmax == pytest.approx(float(np.nanmax(f.top)), abs=0.06)
    assert rep.thickness_max <= 0.25 + 1e-6
    assert rep.thickness_min >= 0.1 - 1e-6


def test_place_on_bed_moves_bottom_to_zero():
    m = wave_plate(base=-2.0)
    _, rep = slice_mesh(m, SliceParams(grid_cell=0.5))
    assert rep.offset[2] == pytest.approx(2.0)


@pytest.mark.skipif(not MODELS.exists(), reason="modèles de test absents")
@pytest.mark.parametrize("name", ["plate_5", "plate_2", "plate_3", "plate_1_1", "plate_6"])
def test_user_models_are_heightfields(name):
    f = from_mesh(trimesh.load(MODELS / f"{name}.stl"), cell=0.5)
    assert np.nanmax(f.thickness) > 0


def _head_envelope():
    from nonplanar.headcheck import load_envelope
    return load_envelope()


def test_head_clears_gentle_surface_but_hits_a_steep_wall():
    """Marche de 14 mm au pied de laquelle la tête (large à partir de 10 mm de hauteur) touche ; loin : rien."""
    from nonplanar.headcheck import check_points
    env = _head_envelope()

    class F:   # champ minimal : xs, ys, top, bottom
        xs = np.linspace(0, 60, 61)
        ys = np.linspace(0, 60, 61)
        top = np.zeros((61, 61))
        bottom = np.zeros((61, 61))

    f = F()
    f.top = np.where(f.xs[None, :] > 30, 15.2, 1.0) * np.ones((61, 1))
    t = lambda x, y: np.full_like(x, 0.2)                # noqa: E731
    # La buse dépose sur la partie basse (z = 1,2), 6 mm avant la marche (matière à 15,0).
    hit, pen = check_points(env, f, lambda x, y: np.where(x > 30, 15.2, 1.2), t,
                            np.array([[24.0, 30.0]]), np.array([1.2]))
    assert hit == 1 and pen > 0.1
    # Loin de la marche : rien.
    hit, _ = check_points(env, f, lambda x, y: np.where(x > 30, 15.2, 1.2), t,
                          np.array([[5.0, 30.0]]), np.array([1.2]))
    assert hit == 0


def test_head_report_in_slice():
    m = wave_plate(amp=0.5, mean=2.0)
    _, rep = slice_mesh(m, SliceParams(grid_cell=0.5))
    assert rep.head.checked and rep.head.points > 0
    assert rep.head.colliding_points == 0                    # ondulation douce
    assert 15 < rep.head.clearance_slope_deg < 45


def test_every_motion_carries_the_six_axes_and_header_describes_kinematics():
    import json
    from file_io import kinematics_description as kin
    gcode, _ = slice_mesh(wave_plate(), SliceParams(grid_cell=0.5, check_head=False))
    lines = gcode.splitlines()
    moves = [l for l in lines if l.startswith(("G0 ", "G1 ")) and any(f" {a}" in l for a in "XYZ")]
    assert moves
    for l in moves:
        words = [w[0] for w in l.split(";")[0].split()[1:]]
        assert "".join(w for w in words if w in "XYZABC") == "XYZABC", l   # tous, dans l'ordre
        assert " A0.000 " in l + " " and " B0.000 " in l + " " and " C0.000 " in l + " "
    desc = kin.parse_header(lines)
    assert desc["machine"] == "cartesian"
    assert [a["letter"] for a in desc["axes"]] == list("XYZABC")
    assert {a["letter"] for a in desc["axes"] if a["type"] == "unused"} == set("ABC")
    json.dumps(desc)


def test_walls_follow_only_the_part_outline_not_internal_region_boundaries():
    from shapely.geometry import box
    from nonplanar import toolpath as tp
    outline = box(0, 0, 20, 20)
    region = box(0, 0, 12, 20)                  # frontière interne en x = 12 (colonne terminée au-delà)
    path = tp.layer_polylines(region, None, 3, tp.ToolpathParams(), outline)
    walls = [pl.xy for pl in path.polylines if pl.kind.startswith("wall")]
    assert walls
    allw = np.vstack(walls)
    assert (allw[:, 0] < 12 + 1e-5).all()                       # rien au-delà de la région
    assert not ((abs(allw[:, 0] - 12) < 0.05) & (allw[:, 1] > 1) & (allw[:, 1] < 19)).any()   # pas de paroi sur x = 12
    infill = np.vstack([pl.xy for pl in path.polylines if pl.kind == "infill"])
    assert infill[:, 0].max() > 11.5                            # le remplissage va jusqu'à la frontière interne


def test_encrypted_header_names_the_machine():
    from file_io.apex_header import cipher_available, read_apex_header
    if not cipher_available():
        pytest.skip("paquet cryptography absent : pas de bloc chiffre")
    gcode, _ = slice_mesh(wave_plate(), SliceParams(grid_cell=0.5, check_head=False))
    block = read_apex_header(gcode.splitlines())
    assert block["sliced"] is True
    assert block["meta"]["machine"] == "cartesian"
    assert block["meta"]["kinematics_description"]["axes"][0]["letter"] == "X"


def test_signature_is_first_line_and_meta_carries_the_version():
    from file_io import apex_header
    from file_io.apex_signature import SIGNATURE, apex_version
    gcode, _ = slice_mesh(wave_plate(), SliceParams(grid_cell=0.5, check_head=False))
    lines = gcode.splitlines()
    assert lines[0] == SIGNATURE
    if apex_header.cipher_available():
        meta = apex_header.read_apex_header(lines)["meta"]
        assert meta["apex_version"] == apex_version() != "unknown"
        assert meta["algorithm"] == "nonplanar-heightfield/1"
        assert meta["algorithm_name"] == "Element Aeolis" and meta["generator"] == "Apex engine (Element Aeolis)"


# -- progression réelle, vectorisations -----------------------------------------------------------
def _quiet_throttle(monkeypatch):
    from core import progress
    monkeypatch.setattr(progress, "MIN_INTERVAL", 0.0)
    monkeypatch.setattr(progress, "MIN_STEP", 0.0)


def test_progress_follows_every_stage_with_real_fractions(monkeypatch):
    """Chaque étape s'ouvre à 0, finit à 1 et ne recule jamais ; les étapes qui bouclent sur du travail
    (surface, remplissage, collision, écriture) émettent des fractions INTERMÉDIAIRES, pas un saut."""
    from nonplanar.slicer import AEOLIS_STEPS
    _quiet_throttle(monkeypatch)
    events = []
    slice_mesh(wave_plate(), SliceParams(), progress=lambda slug, f: events.append((slug, f)))
    slugs = [s for s, _ in AEOLIS_STEPS if s != "import"]
    assert [s for s, f in events if f == 0.0] == slugs
    assert [s for s, f in events if f == 1.0] == slugs
    for slug in slugs:
        fr = [f for s, f in events if s == slug]
        assert fr == sorted(fr) and fr[0] == 0.0 and fr[-1] == 1.0
    for slug in ("heightfield", "toolpath2d", "collision", "export"):
        assert any(0.0 < f < 1.0 for s, f in events if s == slug), slug


def test_slice_path_opens_import_stage(monkeypatch, tmp_path):
    from nonplanar.slicer import slice_path
    _quiet_throttle(monkeypatch)
    wave_plate(n=11).export(tmp_path / "p.stl")
    events = []
    slice_path(str(tmp_path / "p.stl"), SliceParams(check_head=False),
               progress=lambda slug, f: events.append((slug, f)))
    assert events[0] == ("import", 0.0) and ("import", 1.0) in events
    assert events.index(("import", 1.0)) < events.index(("heightfield", 0.0))


def test_progress_is_silent_without_callback():
    from core import progress
    slice_mesh(wave_plate(n=11), SliceParams(check_head=False))
    progress.report(1, 2)          # aucun rappel resté branché


def test_densify_xy_matches_linspace_construction():
    from nonplanar.toolpath import densify_xy
    rng = np.random.default_rng(3)
    xy = rng.random((40, 2)) * 25
    ref = [xy[0:1]]
    for a, b in zip(xy[:-1], xy[1:]):
        m = max(1, int(np.ceil(float(np.linalg.norm(b - a)) / 0.5)))
        ref.append(a + (b - a) * np.linspace(0, 1, m + 1)[1:, None])
    assert np.array_equal(densify_xy(xy, 0.5), np.vstack(ref))


def test_column_matches_stack_calls():
    f = from_mesh(wave_plate())
    stack = build(f)
    rng = np.random.default_rng(4)
    x, y = rng.random(200) * 20, rng.random(200) * 20
    col = stack.sample(x, y)
    for k in (1, 2, 5, stack.count):
        assert np.array_equal(col.z(k), stack.z(k, x, y), equal_nan=True)
        assert np.array_equal(col.thickness(k), stack.thickness(k, x, y), equal_nan=True)
    if stack.count > 1:
        assert np.array_equal(col.z_sub(1, 1), stack.z_sub(1, 1, x, y), equal_nan=True)


def test_rasterize_matches_triangle_loop():
    from nonplanar.heightfield import _rasterize
    m = wave_plate(n=15)
    xs = np.linspace(0, 20, 41)
    ys = np.linspace(0, 20, 41)
    n = m.face_normals
    for sel, keep_max in ((n[:, 2] > 1e-9, True), (n[:, 2] < -1e-9, False)):
        out = np.full((41, 41), np.nan)
        for t in m.triangles[sel]:
            (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) = t
            det = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
            if abs(det) < 1e-14:
                continue
            gx, gy = np.meshgrid(xs, ys)
            a = ((y1 - y2) * (gx - x2) + (x2 - x1) * (gy - y2)) / det
            b = ((y2 - y0) * (gx - x2) + (x0 - x2) * (gy - y2)) / det
            c = 1.0 - a - b
            inside = (a >= -1e-9) & (b >= -1e-9) & (c >= -1e-9)
            cand = np.where(inside, a * z0 + b * z1 + c * z2, np.nan)
            out = np.fmax(out, cand) if keep_max else np.fmin(out, cand)
        assert np.array_equal(_rasterize(m.triangles[sel], xs, ys, keep_max), out, equal_nan=True)


def test_footprint_falls_back_when_top_faces_overlap():
    """L'union « de couverture » exige des triangles qui ne se recouvrent pas ; sinon, union robuste."""
    from nonplanar.heightfield import _footprint
    a = trimesh.creation.box((10, 10, 2))
    b = trimesh.creation.box((10, 10, 2))
    b.apply_translation((5, 0, 0))
    both = trimesh.util.concatenate([a, b])
    assert _footprint(both).area == pytest.approx(150.0, abs=1e-4)
