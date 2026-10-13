"""Lissage non planaire de la voie Bird (pipeline/bird/skin) : briques elementaires.

L'equivalence avec la recherche a ete mesuree hors suite le 2026-10-07 sur le Benchy
incline, passe 4 (CLAUDE.md limite #75) : P2, P3, P4 et P5 IDENTIQUES AU BIT.
"""
import numpy as np
from shapely.geometry import Polygon

from pipeline.bird.skin import surface_patterns
from pipeline.bird.skin.zones import face_adjacency, layer_rank, minimal_stabbing, split_along_level


def test_minimal_stabbing_touches_every_interval_with_fewest_points():
    lo = np.array([0, 2, 5, 6, 9])
    hi = np.array([3, 4, 7, 8, 10])
    points = minimal_stabbing(lo, hi)
    assert all(((points >= a) & (points <= b)).any() for a, b in zip(lo, hi))
    assert len(points) == 3                                     # [3, 7, 10]


def test_layer_rank_is_the_carrying_layer():
    iso = np.array([0.1, 0.3, 0.5])
    assert layer_rank(np.array([0.05, 0.1, 0.2, 0.5, 9.0]), iso).tolist() == [0, 0, 1, 2, 2]


def test_split_along_level_cuts_exactly_and_keeps_area():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], dtype=float)
    faces = np.array([[0, 1, 2], [1, 3, 2]])
    phi = vertices[:, 0] + vertices[:, 1]                       # lineaire sur chaque face
    v, f, p, parent = split_along_level(vertices, faces, phi, np.ones(2, bool), 0.5)
    area = lambda vv, ff: 0.5 * np.linalg.norm(np.cross(vv[ff[:, 1]] - vv[ff[:, 0]], vv[ff[:, 2]] - vv[ff[:, 0]]),
                                               axis=1).sum()
    assert np.isclose(area(v, f), area(vertices, faces))
    assert np.allclose(p[len(vertices):], 0.5)                  # sommets crees sur la ligne de niveau
    # aucune face ne chevauche plus la ligne de niveau
    assert not ((p[f].min(axis=1) < 0.5 - 1e-9) & (p[f].max(axis=1) > 0.5 + 1e-9)).any()
    assert set(parent.tolist()) <= {0, 1}


def test_face_adjacency_is_symmetric():
    faces = np.array([[0, 1, 2], [1, 3, 2], [3, 4, 2]])
    adj = face_adjacency(faces).toarray()
    assert (adj == adj.T).all() and adj[0, 1] and adj[1, 2] and not adj[0, 2]


def test_every_pattern_stays_inside_its_zone():
    zone = Polygon([(0, 0), (12, 0), (12, 8), (0, 8)])
    for pattern in surface_patterns.PATTERNS:
        paths = surface_patterns.generate(pattern, zone, 0.4)
        assert paths, pattern
        pts = np.vstack(paths)
        assert (pts[:, 0] >= -1e-6).all() and (pts[:, 0] <= 12 + 1e-6).all(), pattern
        assert (pts[:, 1] >= -1e-6).all() and (pts[:, 1] <= 8 + 1e-6).all(), pattern


def test_native_skin_check_matches_brute_force():
    """SkinCheckIndex (apex/skin_check.hpp) contre un calcul brut : obstacles = matiere
    permanente + lignes precedentes de la zone + debut de la ligne en cours, distance
    bornee a reach + 1. Equivalence sur le Benchy (2026-10-07) : memes zones exclues et
    memes points de collision que la boucle Python, penetration a 6e-15 mm."""
    import pytest
    from core import _cpp
    if not _cpp.available() or not hasattr(_cpp.module(), "SkinCheckIndex"):
        pytest.skip("noyau natif sans SkinCheckIndex")
    rng = np.random.default_rng(3)
    step, reach = 0.25, 1.2
    body = [rng.uniform(0, 10, (6, 3)) for _ in range(3)]
    zone = [rng.uniform(0, 10, (40, 3)) * [1, 1, 0.3] for _ in range(2)]
    head = rng.normal(scale=2.0, size=(30, 3))

    def dense(p):
        if len(p) < 2:
            return p
        seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
        m = np.maximum(1, np.ceil(seg / step).astype(int))
        t = np.concatenate([np.arange(c) / c for c in m])
        i = np.repeat(np.arange(len(seg)), m)
        return np.vstack([p[i] + t[:, None] * (p[i + 1] - p[i]), p[-1:]])

    index = _cpp.module().SkinCheckIndex(reach + 1.0)
    pts = np.vstack(body)
    starts = np.r_[0, np.cumsum([len(b) for b in body])].astype(np.int64)
    index.add_polylines(pts, starts, step)
    zpts = np.vstack(zone)
    zstarts = np.r_[0, np.cumsum([len(b) for b in zone])].astype(np.int64)
    tested = np.arange(len(zpts), dtype=np.int64)
    zeros = np.zeros((len(zpts), 0))
    plate, dist = index.check_zone(head, np.zeros((0, 3)), zpts, zstarts, zeros, zeros, tested, 150.0, 1.0, reach, step)
    material = [dense(b) for b in body]
    expected = []
    for li, line in enumerate(zone):
        for i in range(len(line)):
            obstacles = material + [dense(z) for z in zone[:li]] + ([dense(line[:i])] if i >= 1 else [])
            obs = np.vstack(obstacles)
            verts = head + line[i]
            d = np.sqrt(((verts[:, None, :] - obs[None, :, :]) ** 2).sum(axis=2)).min()
            expected.append(d if d < reach + 1.0 else np.inf)
    expected = np.array(expected)
    assert np.isfinite(expected).any()
    finite = np.isfinite(expected)
    assert np.array_equal(np.isfinite(dist), finite)
    assert np.allclose(dist[finite], expected[finite], rtol=0, atol=1e-12)
    assert np.array_equal(plate, ((head[None, :, 2] + zpts[:, None, 2] < 1.0)
                                  & (np.linalg.norm(head[None, :, :2] + zpts[:, None, :2], axis=2) < 151.0)).any(axis=1))


def test_native_skin_check_parallel_identical_to_sequential():
    """check_zone_parallel (2026-10-13) : chaque point teste ne voit que le PREFIXE de la
    zone qui lui revient (bizarrerie des blocs de 256 points comprise) ; resultat
    IDENTIQUE AU BIT a check_zone, quel que soit le nombre de fils. Lignes de plus de
    512 points (deux debuts de bloc), points testes en debut de ligne, ligne d'un point,
    axes de rotation reels."""
    import pytest
    from core import _cpp
    if not _cpp.available() or not hasattr(_cpp.module().SkinCheckIndex(1.0), "check_zone_parallel"):
        pytest.skip("noyau natif sans check_zone_parallel")
    rng = np.random.default_rng(11)
    step, reach = 0.25, 1.2
    index = _cpp.module().SkinCheckIndex(reach + 1.0)
    body = [np.cumsum(rng.normal(scale=0.3, size=(200, 3)), axis=0) + [5, 5, 1] for _ in range(3)]
    pts = np.vstack(body)
    index.add_polylines(pts, np.r_[0, np.cumsum([len(b) for b in body])].astype(np.int64), step)
    zone = [np.cumsum(rng.normal(scale=0.15, size=(n, 3)), axis=0) * [1, 1, 0.3] + [5, 5, 2]
            for n in (600, 1, 300, 45)]
    zone[3] += [60.0, 60.0, 0.0]          # ligne isolee : son premier point ne voit aucun obstacle
    zpts = np.vstack(zone)
    zstarts = np.r_[0, np.cumsum([len(z) for z in zone])].astype(np.int64)
    axes = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    rad = np.radians(np.column_stack([rng.uniform(0, 360, len(zpts)), rng.uniform(0, 90, len(zpts))]))
    head = rng.normal(scale=1.5, size=(60, 3))
    all_idx = np.arange(len(zpts), dtype=np.int64)
    sampled = np.unique(np.r_[zstarts[:-1], 256, 512, 600 + 1 + 256, rng.choice(len(zpts), 200, replace=False)])
    for tested in (all_idx, sampled.astype(np.int64)):
        args = (head, axes, zpts, zstarts, np.cos(rad), np.sin(rad), tested, 150.0, 1.0, reach, step)
        plate0, dist0 = index.check_zone(*args)
        assert np.isfinite(dist0).sum() > 50 and np.isinf(dist0).any()      # garde : cas varies
        for threads in (1, 3, 0):
            plate1, dist1 = index.check_zone_parallel(*args, threads=threads)
            assert np.array_equal(plate1, plate0)
            assert np.array_equal(dist1, dist0)
