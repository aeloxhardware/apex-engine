"""Etape 5d : contours des couches planes recoupes dans le STL d'origine (element_ardor/step5d_planar_sections.py)."""
import dataclasses
import os

import numpy as np
import pytest
import trimesh

import main as M
from core import runtime
from element_ardor import step5d_planar_sections as s5d


@pytest.fixture(scope="module")
def cylinder(tmp_path_factory):
    """Cylindre de rayon 6 mm (64 facettes), hauteur 6 mm, tranche jusqu'a l'etape 5 (sans l'etape 5d)."""
    mesh = trimesh.creation.cylinder(radius=6.0, height=6.0, sections=64)
    mesh.apply_translation(-mesh.bounds[0])
    path = os.path.join(str(tmp_path_factory.mktemp("s5d")), "cyl.stl")
    mesh.export(path)
    p = dataclasses.replace(M.PipelineParams(), planar_sections=None)
    with runtime.use(p.runtime):
        stage = lambda *_a, **_k: None
        md, tm, rep = M._ingest_stage(path, p, None, None, stage, [None])
        _om, graph, field = M._analysis_stages(p, md, tm, rep, stage)
        layers = M._layers_stage(p, tm, md, field, graph)
    return md, tm, field, layers


def _polygon_distance(points, radius=6.0, sections=64):
    """Distance au bord du polygone regulier (section exacte du cylindre facette) dans le plan XY."""
    c = np.array([radius, radius])
    xy = points[:, :2] - c
    ang = np.arctan2(xy[:, 1], xy[:, 0])
    half = np.pi / sections
    local = ang % (2 * half) - half                              # angle depuis le milieu de la facette (sommet a 0)
    apothem = radius * np.cos(half)
    return np.abs(np.hypot(xy[:, 0], xy[:, 1]) * np.cos(local) - apothem)


def test_planar_layers_get_the_exact_section_of_the_stl(cylinder):
    md, tm, field, layers = cylinder
    out, report = s5d.run(layers, md, tm, field)
    assert report.planar_layers == len(layers) and report.replaced_loops > 0
    inner = [L for L in out if 0.3 < np.median(np.concatenate(L.contour_3d)[:, 2]) < 5.7]
    pts = np.concatenate([np.concatenate(L.contour_3d) for L in inner])
    assert _polygon_distance(pts).max() < 1e-6          # simplification a 0,01 mm : les sommets du polygone restent
    # meme isovaleur, meme repere, meme nombre de boucles
    for a, b in zip(layers, out):
        assert a.isovalue == b.isovalue and len(a.contour_3d) == len(b.contour_3d)
        assert all(np.allclose(u, v) for u, v in zip(a.local_frame, b.local_frame))


def test_a_loop_without_a_closed_section_is_kept(cylinder):
    """Section ouverte (STL non etanche) : la boucle de l'isosurface est gardee telle quelle."""
    md, tm, field, layers = cylinder
    keep = np.ones(len(md.faces), bool)
    keep[np.argmax(md.vertices[md.faces][:, :, 0].mean(axis=1))] = False      # un trou dans la paroi
    open_md = dataclasses.replace(md, faces=md.faces[keep])
    out, report = s5d.run(layers, open_md, tm, field)
    assert report.kept_loops > 0
    hit = [k for k, (a, b) in enumerate(zip(layers, out)) if a.contour_3d[0] is b.contour_3d[0]]
    assert hit


def test_tilted_layers_are_left_alone(cylinder):
    md, tm, field, layers = cylinder
    tilted = dataclasses.replace(field, gradient=np.tile([0.3, 0.0, 1.0], (len(field.gradient), 1)))
    out, report = s5d.run(layers, md, tm, tilted)
    assert report.planar_layers == 0 and all(a is b for a, b in zip(layers, out))


def test_only_the_loops_far_from_tilted_material_are_resliced(tmp_path):
    """Deux cylindres a 20 mm l'un de l'autre ; champ incline dans le premier seulement (comme la cabine d'un Benchy
    incline) : ses boucles gardent l'isosurface, celles du second sont recoupees dans le STL."""
    a = trimesh.creation.cylinder(radius=4.0, height=4.0, sections=48)
    b = a.copy()
    b.apply_translation([20.0, 0.0, 0.0])
    mesh = trimesh.util.concatenate([a, b])
    mesh.apply_translation(-mesh.bounds[0])
    path = str(tmp_path / "deux.stl")
    mesh.export(path)
    p = dataclasses.replace(M.PipelineParams(), planar_sections=None, allow_multi_mesh=True)
    with runtime.use(p.runtime):
        stage = lambda *_a, **_k: None
        md, tm, rep = M._ingest_stage(path, p, None, None, stage, [None])
        _om, graph, field = M._analysis_stages(p, md, tm, rep, stage)
        layers = M._layers_stage(p, tm, md, field, graph)
    centres = tm.tet_vertices[tm.tetrahedra].mean(axis=1)
    grad = field.gradient.copy()
    grad[centres[:, 0] < 10.0] = [0.3, 0.0, 1.0]
    out, report = s5d.run(layers, md, tm, dataclasses.replace(field, gradient=grad))
    assert report.replaced_loops > 0 and report.planar_layers > 0
    for before, after in zip(layers, out):
        for u, v in zip(before.contour_3d, after.contour_3d):
            if np.asarray(u)[:, 0].mean() < 10.0:
                assert u is v                                  # cote incline : isosurface gardee


def _circle(r, cx=0.0, n=128, wobble=0.0):
    t = np.linspace(0, 2 * np.pi, n)
    rr = r + wobble * np.sign(np.sin(9 * t))                    # creneaux de +/- wobble
    return np.column_stack([cx + rr * np.cos(t), rr * np.sin(t), np.full(n, 47.9)])


def _match(loop, candidates):
    p = s5d.PlanarSectionParams()
    dense = [s5d._densify(c, p.sample_step) for c in candidates]
    trees = [s5d.cKDTree(d[:, :2]) for d in dense]
    return s5d._match(loop, candidates, dense, trees, p)


def test_a_crenellated_loop_near_a_top_is_matched_to_its_unrivalled_section():
    """Derniere couche sous un dessus (cheminee du Benchy) : l'isosurface crenele a +/- 0,6 mm ; la section voisine la
    plus proche est a 2 mm. Avant, p99 > 0,5 mm suffisait a garder le contour crenele (« labyrinthe »)."""
    ring, other = _circle(3.5), _circle(1.5)
    assert _match(_circle(3.5, wobble=0.6), [other, ring]) is ring


def test_an_ambiguous_loose_match_is_refused():
    """Deux sections a peu pres aussi proches : aucune n'est sure, la boucle de l'isosurface est gardee."""
    a, b = _circle(3.5), _circle(4.3)
    assert _match(_circle(3.9, wobble=0.1), [a, b]) is None
    assert _match(_circle(3.5, wobble=0.1), [a, b]) is a            # sous match_tol : la plus proche, sans condition
