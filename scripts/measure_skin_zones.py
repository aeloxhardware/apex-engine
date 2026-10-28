"""
Jalon 1 du lissage non-planaire sur la voie CHALEUR : quelle part de la surface peut recevoir une couche de
surface (P2 de la voie Bird, `pipeline/bird/skin/zones.py` de la branche `bird_old`) ?

Pour chaque face de la surface (re-decoupee a 5 largeurs de buse comme P2) :
  - moment porteur  = rang de la couche qui porte son sommet le plus haut (phi de la methode de la chaleur) ;
  - orientation     = IK de l'etape 8 (C impose par la position), buse vers -n, ecart tete / normale juge au centre
                      et aux trois sommets ;
  - dernier moment accessible = la matiere imprimee est {phi <= phi(t)} (croissante, comme avec Bird) : la vraie
                      tete (marge de securite comprise) posee au centre de la face ne doit pas entrer dans la piece
                      avant que la face soit portee ; plateau = plan relevee de la marge.
  Statuts comme P2 : retenue / inatteignable / collision plateau / collision objet / contact plateau.

Echantillon de faces (--faces=N, graine fixe) : le calcul complet coute des minutes en Python.
Usage : APEX_USE_CPP=0 python scripts/measure_skin_zones.py <modele.stl> [--faces=300] [--tolerance=45]
"""
from __future__ import annotations

import dataclasses
import sys
import time

import numpy as np

ROOT = __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
sys.path.insert(0, ROOT)

import importlib  # noqa: E402

M = importlib.import_module("main")
from pipeline import step1_ingest as s1, step2_overhang as s2, step3_segmentation as s3  # noqa: E402
from pipeline import step4_field_solver as s4, step5_layers as s5  # noqa: E402
from pipeline.skin import zones as sz  # noqa: E402


def main():
    stl = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    opts = dict(a[2:].split("=") for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    n_sample, tolerance = int(opts.get("faces", 300)), float(opts.get("tolerance", 45))
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
    z = sz.compute(tm, field, isovalues, p.machine_profile, p.collision,
                   sz.SkinZoneParams(head_tolerance_deg=tolerance), face_sample=n_sample)
    sel = np.flatnonzero(z.evaluated)
    area = z.areas[sel].sum()
    print(f"{len(z.faces)} faces ({z.areas.sum():.0f} mm2) ; echantillon {len(sel)} faces ({time.time()-t1:.0f} s)")
    for code in (0, 1, 2, 3):
        m = z.status[sel] == code
        print(f"  {sz.STATUS_NAMES[code]:20s} {m.sum():5d} faces  {100*z.areas[sel][m].sum()/area:5.1f} % de l'aire")
    dev = z.deviation_deg[sel]
    print(f"ecart tete/normale : mediane {np.median(dev):.1f} deg, p90 {np.percentile(dev,90):.1f} deg")
    for tol in (5, 15, 30, 45, 60):
        print(f"  ecart <= {tol:2d} deg : {100*z.areas[sel][dev <= tol].sum()/area:5.1f} % de l'aire")


if __name__ == "__main__":
    main()
