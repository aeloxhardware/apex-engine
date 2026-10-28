"""
Mesure du creusement (P3) sur la voie chaleur.

  - couches non touchees : identiques a l'etape 5 (meme objet) ;
  - lignes de coupe (points des contours creuses qui ne sont pas sur la surface d'origine) : leur distance a la peau
    retenue doit valoir l'epaisseur `t` ; on rapporte l'ecart moyen, p95 et maximal ;
  - aire retiree : somme sur les couches de (aire d'origine - aire creusee) x espacement des couches, comparee a
    l'aire retenue x t ;
  - temps.

Usage : APEX_USE_CPP=0 python scripts/measure_skin_carve.py <modele.stl|--box=20> [--zones=p2|all|lateral]
        [--faces=300] [--thickness=0.2]
  zones=p2 : retenues par P2 (pipeline/skin/zones.py, echantillon de faces) ; all : toutes les faces hors plateau ;
  lateral : faces dont la normale est horizontale (|nz| < 0.3).
"""
from __future__ import annotations

import dataclasses
import sys
import tempfile
import time

import numpy as np
import trimesh
from scipy.spatial import cKDTree

ROOT = __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
sys.path.insert(0, ROOT)

import importlib  # noqa: E402

M = importlib.import_module("main")
from pipeline import step1_ingest as s1, step2_overhang as s2, step3_segmentation as s3  # noqa: E402
from pipeline import step4_field_solver as s4, step5_layers as s5  # noqa: E402
from pipeline.skin import carve, zones as sz  # noqa: E402


def polygon_area(loops):
    """Aire projetee sur le plan moyen d'un jeu de boucles (parite : trous soustraits) ; approximation par la
    formule du lacet sur le plan de meilleur ajustement, signe par imbrication non gere : aire = |somme des aires
    signees|, suffisante pour comparer deux versions d'une meme couche."""
    total = 0.0
    for loop in loops:
        pts = np.asarray(loop)
        c = pts - pts.mean(axis=0)
        u, s, vt = np.linalg.svd(c, full_matrices=False)
        xy = c @ vt[:2].T
        a = 0.5 * np.sum(xy[:-1, 0] * xy[1:, 1] - xy[1:, 0] * xy[:-1, 1])
        total += a
    return abs(total)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = dict(a[2:].split("=") for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    thickness = float(opts.get("thickness", 0.2))
    mode = opts.get("zones", "p2")
    if "box" in opts:
        stl = tempfile.mktemp(suffix=".stl")
        e = float(opts["box"])
        b = trimesh.creation.box(extents=(e, e, e))
        b.export(stl)
    else:
        stl = args[0]
    p = M.PipelineParams()
    t0 = time.time()
    md, tm, rep = s1.run(stl, resolution=p.resolution, near_platform_resolution=p.near_platform_resolution,
                         near_platform_margin=p.near_platform_margin, place_on_platform=p.place_on_platform)
    om = s2.run(md, p.overhang)
    seg = dataclasses.replace(p.segmentation, platform_z=rep.platform_z_used) if p.segmentation.platform_z is None else p.segmentation
    graph = s3.run(md, om, seg)
    field = s4.run(tm, md, om, graph, M._effective_field_solver_params(p))
    layers = s5.run(tm, md, field, graph, p.layers)
    isovalues = np.sort(np.unique([l.isovalue for l in layers]))
    print(f"{stl}: {len(tm.tetrahedra)} tetra, {len(layers)} couches, phi max {field.values.max():.1f} ({time.time()-t0:.0f} s)", flush=True)

    t1 = time.time()
    if mode == "p2":
        z = sz.compute(tm, field, isovalues, p.machine_profile, p.collision, face_sample=int(opts.get("faces", 300)))
        keep = z.retained
    else:
        solid = sz.PrintedSolid(tm.tet_vertices, tm.tetrahedra, field.values)
        sv, sf, nrm, areas, _c = sz.surface_mesh(tm, solid, 2.0)
        free = sv[sf][:, :, 2].max(axis=1) >= tm.tet_vertices[:, 2].min() + sz.PLATE_CONTACT_TOL
        keep = np.flatnonzero(free & ((np.abs(nrm[:, 2]) < 0.3) if mode == "lateral" else True))
        z = sz.SkinZones(sv, sf, nrm, areas, np.zeros(len(sf), np.int8), np.zeros(len(sf), int), np.zeros(len(sf)),
                         np.ones(len(sf), bool))
    zone_faces = z.faces[keep]
    zone_area = z.areas[keep].sum()
    print(f"zones ({mode}) : {len(keep)} faces, {zone_area:.0f} mm2 ({time.time()-t1:.0f} s)", flush=True)

    cp = carve.CarveParams(thickness=thickness)
    t2 = time.time()
    zone = carve.ZoneDistance(z.vertices, zone_faces, cp)
    carved, rep_c = carve.carve_layers(layers, tm, field, zone, cp)
    print(f"creusement {rep_c.seconds['carve']:.1f} s ({rep_c.seconds['distance_calls']} distances exactes) : {rep_c.touched_layers} couches touchees sur "
          f"{len(layers)}, {rep_c.dropped_layers} supprimees, {rep_c.nodes_in_skin} nœuds dans la peau", flush=True)

    by_id = {l.id: l for l in carved}
    same = sum(1 for l in layers if l.id in by_id and by_id[l.id] is l)
    print(f"couches rendues telles quelles : {same}")
    # lignes de coupe : points hors de la surface d'origine
    surface = trimesh.Trimesh(tm.tet_vertices, sz.boundary_faces(tm.tetrahedra), process=False)
    zmesh = trimesh.Trimesh(z.vertices, zone_faces, process=False) if len(zone_faces) else None
    pts = []
    removed = 0.0
    spacing = float(np.median(np.diff(isovalues))) if len(isovalues) > 1 else 0.0
    for l in carved:
        if l.contour_3d is not layers[[x.id for x in layers].index(l.id)].contour_3d:
            orig = layers[[x.id for x in layers].index(l.id)]
            removed += max(0.0, polygon_area(orig.contour_3d) - polygon_area(l.contour_3d)) * spacing
            allp = np.concatenate([np.asarray(c)[:-1] for c in l.contour_3d])
            sample = allp[:: max(1, len(allp) // 400)]
            pts.append(sample)
    if pts and zmesh is not None:
        pts = np.concatenate(pts)
        _c, d_surface, _t = trimesh.proximity.closest_point(surface, pts)
        cutp = pts[d_surface > 1e-4]
        _c, d_zone, _t = trimesh.proximity.closest_point(zmesh, cutp)
        err = d_zone - thickness
        print(f"lignes de coupe : {len(cutp)} points ; distance a la peau retenue - t : moyenne {err.mean():+.3f} mm, "
              f"p95 {np.percentile(np.abs(err),95):.3f} mm, max {np.abs(err).max():.3f} mm (t = {thickness})")
    print(f"volume retire INDICATIF (aires projetees sur le plan moyen, peu fiable sur couches courbes) {removed:.0f} mm3 ; aire retenue x t = {zone_area*thickness:.0f} mm3")


if __name__ == "__main__":
    main()
