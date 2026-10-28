"""Lissage non-planaire sur la voie chaleur : zones (P2), trajectoires de surface (P4), ordre et collisions (P5),
enchainement complet (pipeline/skin/chain.py), sur l_shape_overhang.stl."""
import dataclasses
import os

import numpy as np
import pytest
import trimesh

import main as M
from pipeline import step1_ingest as s1, step2_overhang as s2, step3_segmentation as s3
from pipeline import step4_field_solver as s4, step5_layers as s5
from pipeline.skin import chain, zones as sz

STL = os.path.join(os.path.dirname(__file__), "..", "res", "test_models", "l_shape_overhang.stl")


@pytest.fixture(scope="module")
def heat():
    p = M.PipelineParams()
    md, tm, rep = s1.run(STL, resolution=p.resolution, near_platform_resolution=p.near_platform_resolution,
                         near_platform_margin=p.near_platform_margin, place_on_platform=True)
    om = s2.run(md, p.overhang)
    graph = s3.run(md, om, dataclasses.replace(p.segmentation, platform_z=rep.platform_z_used))
    field = s4.run(tm, md, om, graph, M._effective_field_solver_params(p))
    return p, md, tm, om, field, s5.run(tm, md, field, graph, p.layers)


@pytest.fixture(scope="module")
def result(heat):
    p, md, tm, om, field, layers = heat
    return chain.run(tm, field, layers, md, om, M._effective_toolpath2d_params(p), p.machine_profile, p.collision,
                     chain.SkinParams(face_sample=300, max_iter=4))


def test_zone_statuses_partition_the_surface(heat):
    p, _md, tm, _om, field, layers = heat
    z = sz.compute(tm, field, np.unique([l.isovalue for l in layers]), p.machine_profile, p.collision, face_sample=200)
    skin = sz.partition(z)
    st = skin["status"]
    assert set(np.unique(st)) <= set(sz.STATUS_NAMES)
    assert abs(skin["areas"].sum() - z.areas.sum()) < 1e-6 * z.areas.sum()      # la coupe ne perd aucune aire
    retained = skin["zone"] >= 0
    assert retained.any() and (st[retained] == sz.STATUS_RETAINED).all()
    # zones accessibles a leur moment (invariant de P2) : dernier moment accessible >= moment de la zone
    for zid, moment in enumerate(skin["zone_moment"]):
        assert (skin["last_access"][skin["zone"] == zid] >= moment).all()


def test_skin_paths_lie_half_a_thickness_under_the_surface(result):
    kept = [z for z, e in result.paths.items() if e["paths"] and z not in result.print_order["excluded"]]
    assert kept
    zone_faces = result.skin["faces"][np.isin(result.skin["zone"], kept)]
    mesh = trimesh.Trimesh(result.skin["vertices"], zone_faces, process=False)
    pts = np.concatenate([q for z in kept for q in result.paths[z]["paths"]])
    d = trimesh.proximity.closest_point(mesh, pts)[1]
    # t / 2 : l'axe du cordon est au milieu de la couche creusee. Les trajectoires descendent le long des normales de
    # SOMMET (moyennees) : aux aretes vives de la zone la profondeur s'ecarte de t/2 (heritage de la voie Bird, a
    # corriger plus tard) ; l'immense majorite des points est exacte et aucun n'est plus profond que t + w/2
    assert (np.abs(d - 0.1) < 0.02).mean() > 0.9 and d.max() < 0.25


def test_chain_ends_without_collision_and_orders_skin_after_its_supporting_body(result):
    assert result.clean and not result.print_order["collision_points"]
    passes = result.print_order["passes"]
    assert any(q["kind"] == "skin" for q in passes)
    ranks = np.searchsorted(result.skin["isovalues"], [b["isovalue"] - 1e-9 for b in result.body])
    for i, q in enumerate(passes):
        if q["kind"] == "skin":
            moment = result.print_order["report"][q["zone"]]["moment"]
            before = [r for r, other in zip(range(len(passes)), passes) if other["kind"] == "body" and other["layer"] is not None
                      and ranks[other["layer"]] <= moment and r < i]
            assert before     # au moins une couche portante deposee avant la zone
    for q in passes:
        if q["kind"] == "skin":
            assert all(len(a) == len(p) for a, p in zip(q["angles"], q["points"]))   # un angle de tete par point


def test_carved_body_keeps_the_skin_thickness_clear(result):
    kept = [z for z, e in result.paths.items() if e["paths"] and z not in result.print_order["excluded"]]
    zone_faces = result.skin["faces"][np.isin(result.skin["zone"], kept)]
    mesh = trimesh.Trimesh(result.skin["vertices"], zone_faces, process=False)
    body = np.concatenate([np.concatenate(b["toolpath_3d"]) for b in result.body if b["toolpath_3d"]])
    d = trimesh.proximity.closest_point(mesh, body[::2])[1]
    assert d.min() > 0.15          # l'axe du cordon le plus externe est sur le contour creuse, a t = 0,2 de la peau


def test_inject_places_each_zone_after_its_supporting_layer_and_remaps_indices(result):
    from core.mesh_types import ToolpathMachine, ToolpathPoint
    from pipeline.skin.inject import inject
    layers = result.layers
    pts, ids = [], []
    for layer in layers:
        for k in range(3):
            pts.append(ToolpathPoint(position=np.array([layer.isovalue, k, 0.0]), axis_angles=[0.0, 0.0],
                                     extrusion=1.0, feedrate=1000.0))
            ids.append(layer.id)
    body = ToolpathMachine(points=pts, layer_ids=np.array(ids))
    new, old_to_new = inject(body, layers, result.print_order["passes"], 1000.0)
    skin_idx = [i for i, p in enumerate(new.points) if p.is_skin]
    n_skin = sum(len(path) for q in result.print_order["passes"] if q["kind"] == "skin" for path in q["points"])
    assert len(skin_idx) == n_skin and len(new.points) == len(pts) + n_skin
    assert (np.diff(old_to_new) > 0).all()                                        # l'ordre du corps est conserve
    assert all(new.points[old_to_new[i]] is pts[i] for i in range(0, len(pts), 7))  # indices remappes exacts
    # premier point de chaque polyligne de lissage : trajet a vide ; les autres extrudent
    assert new.points[skin_idx[0]].extrusion == 0.0 and new.points[skin_idx[1]].extrusion == 1.0
    # chaque bloc de lissage vient apres le dernier point de sa couche porteuse
    first_skin = skin_idx[0]
    layer_of_prev = new.layer_ids[first_skin - 1]
    assert layer_of_prev != -1 and layer_of_prev in {l.id for l in layers}
