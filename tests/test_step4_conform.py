"""Couches conformes sous les grandes surfaces en pente (element_ardor/step4_conform.py)."""
import dataclasses
import os

import numpy as np
import pytest
import trimesh

import main as M
from core import runtime
from element_ardor.step4_conform import ConformParams
from element_ardor.step4_planar import PlanarParams

SLOPE = 0.2                     # dessus en pente : z = 8 + 0,2 x (11,3 deg), x de 0 a 24


def _wedge() -> trimesh.Trimesh:
    xs, ys = (0.0, 24.0), (0.0, 16.0)
    pts = [[x, y, 0.0] for x in xs for y in ys] + [[x, y, 8.0 + SLOPE * x] for x in xs for y in ys]
    return trimesh.convex.convex_hull(np.array(pts))


def _field(tmp_path, conform):
    path = os.path.join(str(tmp_path), "wedge.stl")
    _wedge().export(path)
    p = M.PipelineParams()
    p.field_solver = dataclasses.replace(p.field_solver, planar=PlanarParams(conform=conform))
    with runtime.use(p.runtime):
        stage = lambda *_a, **_k: None
        md, tm, rep = M._ingest_stage(path, p, None, None, stage, [None])
        _om, _graph, field = M._analysis_stages(p, md, tm, rep, stage)
    return tm, field.values


@pytest.fixture(scope="module")
def fields(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("conform")
    tm, planar = _field(tmp, None)
    _tm2, conform = _field(tmp, ConformParams(depth=6.0))
    return tm.tet_vertices, planar, conform


def _top(x: np.ndarray) -> np.ndarray:
    return x[:, 2] >= 8.0 + SLOPE * x[:, 0] - 0.3


def test_the_sloped_top_becomes_a_layer(fields):
    """Plan : phi suit z le long du dessus (variation 4,8 mm) ; conforme : phi est presque constant sur le dessus."""
    x, planar, conform = fields
    top = _top(x) & (x[:, 0] > 1.0) & (x[:, 0] < 23.0)
    assert np.ptp(planar[top]) > 4.0
    assert np.ptp(conform[top]) < 0.25 * np.ptp(planar[top])


def test_layers_deep_below_the_top_are_unchanged(fields):
    x, planar, conform = fields
    deep = x[:, 2] < 8.0 + SLOPE * x[:, 0] - 6.5
    assert deep.sum() > 20
    assert np.abs(conform[deep] - planar[deep]).max() < 0.05


def test_the_field_still_grows_upwards(fields):
    """Le decalage ne retourne jamais les couches : phi croit avec z dans chaque colonne de noeuds."""
    x, _planar, conform = fields
    key = np.round(x[:, :2], 1)
    _u, col = np.unique(key, axis=0, return_inverse=True)
    col = col.ravel()
    for c in np.unique(col):
        sel = np.flatnonzero(col == c)
        if len(sel) < 3:
            continue
        order = sel[np.argsort(x[sel, 2])]
        assert np.all(np.diff(conform[order]) > -1e-6)


def test_narrow_tops_are_left_planar(tmp_path):
    """Un dessus en pente plus etroit que `min_width` n'est pas retenu : le champ est celui du plan."""
    _tm, planar = _field(tmp_path, None)
    _tm2, narrow = _field(tmp_path, ConformParams(min_width=50.0))
    assert np.abs(planar - narrow).max() < 1e-9


def test_a_ceiling_under_a_conformal_top_stays_flat(tmp_path):
    """Cavite fermee sous le dessus en pente (plafond horizontal a z = 6) : le plafond reste une couche plane (il
    s'imprime en pont accroche aux parois), le dessus devient conforme. Sans cela, les couches paralleles au dessus
    refermaient le plafond en biais, couche apres couche, chacune posee dans le vide (Benchy : plafond de la cabine)."""
    cavity = trimesh.creation.box(extents=[16.0, 8.0, 4.0])
    cavity.apply_translation([12.0, 8.0, 4.0])
    solid = _wedge().difference(cavity)
    path = os.path.join(str(tmp_path), "hollow.stl")
    solid.export(path)
    p = M.PipelineParams()
    p.field_solver = dataclasses.replace(p.field_solver, planar=PlanarParams(conform=ConformParams(depth=6.0)))
    with runtime.use(p.runtime):
        stage = lambda *_a, **_k: None
        md, tm, rep = M._ingest_stage(path, p, None, None, stage, [None])
        _om, _graph, field = M._analysis_stages(p, md, tm, rep, stage)
    x, phi = tm.tet_vertices, field.values
    ceiling = (np.abs(x[:, 2] - 6.0) < 0.15) & (x[:, 0] > 5.0) & (x[:, 0] < 19.0) & (np.abs(x[:, 1] - 8.0) < 3.0)
    top = _top(x) & (x[:, 0] > 1.0) & (x[:, 0] < 23.0)
    assert ceiling.sum() > 5 and top.sum() > 5
    assert np.ptp(phi[ceiling]) < 0.3                       # plafond : une couche
    assert np.ptp(phi[top]) < 0.25 * (SLOPE * 22.0)        # dessus : conforme


def test_the_head_keeps_the_printing_direction_under_conformal_layers(tmp_path):
    """Sous le dessus conforme, le gradient s'incline avec les couches, mais la direction d'impression
    (`ScalarField.direction`, axe de la tete a l'etape 8) reste celle du champ plan : verticale. Avant, la tete
    s'inclinait et passait sous le plateau (Benchy : 91 mouvements sous la marge)."""
    path = os.path.join(str(tmp_path), "wedge.stl")
    _wedge().export(path)
    p = M.PipelineParams()
    p.field_solver = dataclasses.replace(p.field_solver, planar=PlanarParams(conform=ConformParams(depth=6.0)))
    with runtime.use(p.runtime):
        stage = lambda *_a, **_k: None
        md, tm, rep = M._ingest_stage(path, p, None, None, stage, [None])
        _om, _graph, field = M._analysis_stages(p, md, tm, rep, stage)
    assert field.direction is not None
    centres = tm.tet_vertices[tm.tetrahedra].mean(axis=1)
    top = (centres[:, 2] > 8.0 + SLOPE * centres[:, 0] - 1.0) & (centres[:, 0] > 2.0) & (centres[:, 0] < 22.0)
    unit = lambda g: g / np.linalg.norm(g, axis=1)[:, None]
    tilt = lambda g: np.degrees(np.arccos(np.clip(unit(g[top])[:, 2], -1.0, 1.0)))
    assert np.median(tilt(field.gradient)) > 5.0          # les couches suivent la pente (11,3 deg)
    assert np.median(tilt(field.direction)) < 1.0         # la tete reste verticale


def test_without_conformal_layers_there_is_no_separate_direction(tmp_path):
    path = os.path.join(str(tmp_path), "wedge.stl")
    _wedge().export(path)
    p = M.PipelineParams()
    p.field_solver = dataclasses.replace(p.field_solver, planar=PlanarParams(conform=None))
    with runtime.use(p.runtime):
        stage = lambda *_a, **_k: None
        md, tm, rep = M._ingest_stage(path, p, None, None, stage, [None])
        _om, _graph, field = M._analysis_stages(p, md, tm, rep, stage)
    assert field.direction is None                        # l'etape 8 suit alors le gradient, comme avant


def _levels(tmp_path, mesh, conform):
    path = os.path.join(str(tmp_path), "piece.stl")
    mesh.export(path)
    p = M.PipelineParams()
    p.field_solver = dataclasses.replace(p.field_solver, planar=PlanarParams(conform=conform))
    with runtime.use(p.runtime):
        stage = lambda *_a, **_k: None
        md, tm, rep = M._ingest_stage(path, p, None, None, stage, [None])
        _om, _graph, field = M._analysis_stages(p, md, tm, rep, stage)
    return tm, field


def test_the_sloped_top_and_a_flat_top_are_levels(tmp_path):
    """Le dessus en pente rendu plan et un dessus plat sont des niveaux du champ (limite #95), le plus grand d'abord :
    leur valeur est celle du champ sur la surface."""
    step = trimesh.creation.box(extents=(6.0, 6.0, 5.3))
    step.apply_translation((-4.0, 8.0, 2.65))                      # dessus plat a z = 5,3, colle au coin
    tm, field = _levels(tmp_path, trimesh.boolean.union([_wedge(), step], engine="manifold"), ConformParams(depth=6.0))
    assert field.levels is not None and len(field.levels) == 2
    x = tm.tet_vertices
    top = _top(x) & (x[:, 0] > 1.0) & (x[:, 0] < 23.0)
    assert abs(field.levels[0] - np.median(field.values[top])) < 0.1
    assert abs(field.levels[1] - 5.3) < 0.1


def test_without_conformal_layers_there_are_no_levels(tmp_path):
    _tm, field = _levels(tmp_path, _wedge(), None)
    assert field.levels is None
