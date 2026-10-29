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
