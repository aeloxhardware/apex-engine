"""Tranchage non planaire (paquet `element_aeolis/`) : maillages synthétiques générés à la volée."""
from pathlib import Path

import numpy as np
import pytest
import trimesh

from element_aeolis import NotAHeightfieldError, SliceParams, build, from_mesh, slice_mesh

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
    # Ce test porte sur le moteur « surface en hauteur » : on lève la garde de la buse (sinon repli général).
    gcode, rep = slice_mesh(m, SliceParams(grid_cell=0.5, enforce_nozzle_limits=False))
    assert rep.mode == "heightfield" and rep.layer_count > 5
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
    from element_aeolis.headcheck import load_envelope
    return load_envelope()


def test_head_clears_gentle_surface_but_hits_a_steep_wall():
    """Marche de 14 mm au pied de laquelle la tête (large à partir de 10 mm de hauteur) touche ; loin : rien."""
    from element_aeolis.headcheck import check_points
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


def test_head_check_on_a_part_lower_than_the_head_does_not_crash():
    """Pièce plus basse que tout point de l'enveloppe : aucun contact possible, pas de min() sur tableau vide."""
    from element_aeolis.headcheck import check_points
    env = _head_envelope()

    class F:
        xs = np.linspace(0, 10, 11)
        ys = np.linspace(0, 10, 11)
        top = np.full((11, 11), 0.4)
        bottom = np.zeros((11, 11))

    t = lambda x, y: np.full_like(x, 0.2)                # noqa: E731
    assert check_points(env, F(), lambda x, y: np.full_like(x, 0.4), t,
                        np.array([[5.0, 5.0]]), np.array([0.4])) == (0, 0.0)


def test_head_report_in_slice():
    m = wave_plate(amp=0.5, mean=2.0)
    _, rep = slice_mesh(m, SliceParams(grid_cell=0.5))
    assert rep.head.checked and rep.head.points > 0
    assert rep.head.colliding_points == 0                    # ondulation douce
    assert rep.head.clearance_slope_deg == pytest.approx(SliceParams().nozzle_angle_deg)   # enveloppe des réglages de buse
    # Une vraie tête (GLB) remplace l'enveloppe des réglages.
    from pathlib import Path

    from config.machine_profiles import default_toolhead_id
    from file_io.mesh_loader import toolhead_path
    head = Path(toolhead_path(default_toolhead_id()))
    if head.exists():
        _, rep = slice_mesh(m, SliceParams(grid_cell=0.5, head_model=str(head)))
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
    from element_aeolis import toolpath as tp
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
    from element_aeolis.slicer import AEOLIS_STEPS
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
    from element_aeolis.slicer import slice_path
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
    from element_aeolis.toolpath import densify_xy
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
    from element_aeolis.heightfield import _rasterize
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
    from element_aeolis.heightfield import _footprint
    a = trimesh.creation.box((10, 10, 2))
    b = trimesh.creation.box((10, 10, 2))
    b.apply_translation((5, 0, 0))
    both = trimesh.util.concatenate([a, b])
    assert _footprint(both).area == pytest.approx(150.0, abs=1e-4)


# --------------------------------------------------------------------------------- tranchage général

def cavity_box(outer=30.0, inner=20.0, height=12.0, cavity_h=6.0) -> trimesh.Trimesh:
    """Boîte fermée creuse : une cavité interne (ni surface en hauteur, ni imprimable par un moteur 2,5D)."""
    box = trimesh.creation.box([outer, outer, height])
    box.apply_translation([0, 0, height / 2])
    hole = trimesh.creation.box([inner, inner, cavity_h])
    hole.apply_translation([0, 0, height / 2])
    return box.difference(hole)


def test_cone_minorant_matches_brute_force():
    from element_aeolis.general import cone_minorant
    rng = np.random.default_rng(3)
    f = np.full((14, 17), np.inf)
    for _ in range(4):
        f[rng.integers(0, 14), rng.integers(0, 17)] = rng.uniform(0, 5)
    step = 0.3
    got = cone_minorant(f, step)
    jj, ii = np.indices(f.shape)
    src = np.argwhere(np.isfinite(f))
    brute = np.min([f[a, b] + step * np.hypot(jj - a, ii - b) for a, b in src], axis=0)
    assert (got >= brute - 1e-9).all()                       # le chanfrein surestime la distance euclidienne...
    src_min = np.min([f[a, b] + 1.0824 * step * np.hypot(jj - a, ii - b) for a, b in src], axis=0)
    assert (got <= src_min + 1e-9).all()                     # ... d'au plus 8,24 % (1/cos 22,5°)
    assert (got <= f + 1e-12).all()


def test_relief_respects_the_nozzle_limits():
    from element_aeolis.general import build_relief
    xs = ys = np.arange(0, 40.001, 0.5)
    X, Y = np.meshgrid(xs, ys)
    top = 10.0 + 3.0 * np.exp(-((X - 20) ** 2 + (Y - 20) ** 2) / 40.0)          # bosse gaussienne
    for angle, reach in ((10.0, 2.0), (25.0, 1.0), (45.0, 4.0)):
        r = build_relief(top, xs, ys, angle, reach, z_bottom=0.0)
        assert r.d.max() <= 1e-9 and r.d.min() >= -reach - 1e-9
        assert r.max_gradient <= np.tan(np.radians(angle)) * 1.01
        # D(q) − D(p) ≤ min(r·tan α, reach) sur des couples tirés au hasard
        rng = np.random.default_rng(0)
        a, b = rng.integers(0, len(xs), (2, 400)), rng.integers(0, len(xs), (2, 400))
        dist = np.hypot(xs[a[0]] - xs[b[0]], ys[a[1]] - ys[b[1]])
        diff = r.d[a[1], a[0]] - r.d[b[1], b[0]]
        assert (diff <= np.minimum(dist * np.tan(np.radians(angle)) * 1.01, reach) + 1e-9).all()


def test_relief_vanishes_without_nozzle_slack():
    from element_aeolis.general import build_relief
    xs = ys = np.arange(0, 10.001, 0.5)
    top = 5.0 + np.outer(np.linspace(0, 1, len(ys)), np.ones(len(xs)))
    assert np.abs(build_relief(top, xs, ys, 0.0, 2.0, 0.0).d).max() == 0.0     # angle nul : planaire
    assert np.abs(build_relief(top, xs, ys, 30.0, 0.0, 0.0).d).max() == 0.0    # buse sans longueur : planaire


def test_section_region_keeps_cavities_as_holes():
    from element_aeolis.general import section_region
    m = cavity_box()
    z = m.vertices[:, 2]
    r = section_region(m.vertices, m.faces, z - 6.0)         # coupe au milieu de la cavité
    assert r.area == pytest.approx(30 * 30 - 20 * 20, rel=1e-6)
    r = section_region(m.vertices, m.faces, z - 2.0)         # sous la cavité : plein
    assert r.area == pytest.approx(30 * 30, rel=1e-6)


def test_cavity_model_is_sliced_by_default_with_the_general_engine():
    m = cavity_box()
    gcode, rep = slice_mesh(m, SliceParams())
    assert rep.mode == "general" and "volume" in rep.fallback_reason
    assert rep.layer_count >= 55
    assert rep.head.checked and rep.head.points > 0 and rep.head.colliding_points == 0
    assert rep.thickness_min >= 0.1 - 1e-6 and rep.thickness_max <= 0.25 + 1e-6
    zs = [float(l.split(" Z")[1].split()[0]) for l in gcode.splitlines() if l.startswith("G1") and " E" in l and " Z" in l]
    assert min(zs) == pytest.approx(0.2, abs=1e-3) and max(zs) == pytest.approx(12.0, abs=0.02)
    assert '"aeolis_mode"' not in gcode                       # le bloc d'en-tête est chiffré, pas en clair


def test_general_engine_fills_the_part_volume():
    from element_aeolis.toolpath import ToolpathParams
    m = cavity_box()
    _, rep = slice_mesh(m, SliceParams(toolpath=ToolpathParams(infill_density=1.0)))
    volume = rep.filament_length * np.pi * (1.75 / 2) ** 2
    assert 0.85 * m.volume < volume < 1.05 * m.volume


def test_mode_heightfield_refuses_and_mode_general_accepts_a_heightfield():
    with pytest.raises(NotAHeightfieldError):
        slice_mesh(cavity_box(), SliceParams(mode="heightfield"))
    _, rep = slice_mesh(wave_plate(), SliceParams(mode="general", grid_cell=0.5))
    assert rep.mode == "general" and rep.fallback_reason == ""
    with pytest.raises(ValueError):
        slice_mesh(wave_plate(), SliceParams(mode="inconnu"))


def test_default_nozzle_limits_send_a_steep_plate_to_the_general_engine_and_raised_limits_keep_it():
    m = wave_plate(amp=1.5, mean=3.0)                          # pente maximale environ 25°
    _, low = slice_mesh(m, SliceParams(grid_cell=0.5))
    assert low.mode == "general" and "limites de la buse" in low.fallback_reason
    assert low.slope_max_deg <= SliceParams().nozzle_angle_deg + 1e-6       # les couches respectent l'angle
    assert low.head.colliding_points == 0
    _, high = slice_mesh(m, SliceParams(grid_cell=0.5, nozzle_angle_deg=80.0, nozzle_reach=50.0))
    assert high.mode == "heightfield" and high.head.colliding_points == 0
    _, off = slice_mesh(m, SliceParams(grid_cell=0.5, enforce_nozzle_limits=False))
    assert off.mode == "heightfield"                                        # l'ancien comportement reste disponible


def wedge_plate(size=20.0, deg=8.0, base=4.0) -> trimesh.Trimesh:
    """Prisme : dessous plat, dessus en pente douce de `deg` degrés, flancs verticaux."""
    rise = np.tan(np.radians(deg)) * size
    verts = np.array([[0, 0, 0], [size, 0, 0], [size, size, 0], [0, size, 0],
                      [0, 0, base], [size, 0, base + rise], [size, size, base + rise], [0, size, base]], float)
    faces = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                      [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    m = trimesh.Trimesh(verts, faces, process=True)
    m.fix_normals()
    return m


def test_refinement_keeps_the_mesh_watertight_and_the_volume():
    from element_aeolis.general import _refine
    m = cavity_box()
    v, f = _refine(m.vertices, m.faces, 0.7)
    r = trimesh.Trimesh(v, f, process=False)
    assert len(f) > 10 * len(m.faces)
    assert r.is_watertight and r.volume == pytest.approx(m.volume, rel=1e-9)   # aucune jonction en T
    assert np.max(np.linalg.norm(v[f][:, :, :2] - np.roll(v[f][:, :, :2], 1, axis=1), axis=2)) <= 0.7 + 1e-9


def test_relief_does_not_lose_material():
    """Le relief déplace les couches, pas la matière : même quantité de filament avec et sans relief. (Le défaut
    d'origine, un raffinement à jonctions en T qui ouvrait les contours, faisait tomber death.stl de 113 000 à
    5 010 mm extrudés ; le test d'étanchéité du raffinement le couvre directement.)"""
    m = wedge_plate()
    flat = slice_mesh(m, SliceParams(mode="general", nozzle_angle_deg=0.0))[1]
    warped = slice_mesh(m, SliceParams(mode="general", nozzle_angle_deg=10.0, nozzle_reach=2.0))[1]
    assert warped.relief_mm > 0.1 and flat.relief_mm == 0.0
    assert warped.filament_length == pytest.approx(flat.filament_length, rel=0.06)


def test_general_layers_follow_a_gentle_top_and_stay_inside_the_limits():
    """Plaque à dessus en pente douce (8°) : la pente des couches ne dépasse jamais l'angle de la buse ni leur dénivelé
    sa longueur, et elles épousent le dessus (écart au dessus bien plus petit que sans relief)."""
    from element_aeolis.general import analyze
    kw = dict(layer_height=None, first_layer=None, min_layer=0.1, max_layer=0.25)
    plan = analyze(wedge_plate(), nozzle_angle_deg=10.0, nozzle_reach=2.0, **kw)
    assert plan.relief.amp > 0.1
    top_layer = plan.stack.count
    for k in (top_layer // 2, top_layer):
        z = plan.stack.z_at(k, plan.relief.d)
        gy, gx = np.gradient(z, plan.relief.ys, plan.relief.xs)
        assert np.hypot(gx, gy).max() <= np.tan(np.radians(10.0)) * 1.01
        assert z.max() - z.min() <= 2.0 + 1e-9
    flat = analyze(wedge_plate(), nozzle_angle_deg=0.0, nozzle_reach=2.0, **kw)
    assert flat.relief.amp == 0.0


def _type_runs(gcode: str):
    """Pour chaque étiquette ;TYPE: le nombre de mouvements extrudés qui la suivent, par type."""
    cur, counts = None, {}
    for line in gcode.splitlines():
        if line.startswith(";TYPE:"):
            cur = line[6:]
        elif line.startswith("; <layer>"):
            cur = None
        elif line.startswith("G1") and " X" in line and " E" in line and cur is not None:
            counts[cur] = counts.get(cur, 0) + 1
        elif line.startswith("G1") and " X" in line and " E" in line:
            counts["<sans type>"] = counts.get("<sans type>", 0) + 1
    return counts


@pytest.mark.parametrize("engine", ["heightfield", "general"])
def test_every_extruded_move_carries_a_line_type(engine):
    m = wave_plate(amp=0.5, mean=2.0) if engine == "heightfield" else cavity_box()
    gcode, rep = slice_mesh(m, SliceParams(grid_cell=0.5, enforce_nozzle_limits=False)
                            if engine == "heightfield" else SliceParams())
    assert rep.mode == engine
    counts = _type_runs(gcode)
    assert "<sans type>" not in counts
    assert set(counts) <= {"WALL-OUTER", "WALL-INNER", "FILL", "SKIN"}
    assert {"WALL-OUTER", "WALL-INNER", "SKIN"} <= set(counts)       # parois et peau toujours présentes
    assert counts["WALL-OUTER"] > 0 and counts.get("FILL", 0) + counts["SKIN"] > counts["WALL-OUTER"] * 0.2


def test_solid_infill_is_marked_apart_from_sparse_infill():
    from shapely.geometry import box
    from element_aeolis import toolpath as tp
    region = box(0, 0, 20, 20)
    path = tp.layer_polylines(region, box(0, 0, 20, 10), 3, tp.ToolpathParams(), region)
    kinds = {pl.kind for pl in path.polylines}
    assert {"infill", "solid_infill", "wall_outer", "wall_inner"} <= kinds
