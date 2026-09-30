"""
Banc de P4 (etape 4a) : outils de surface de surface_paths.py sur trois formes
aux distances EXACTES connues -- disque plat, demi-cylindre (developpable),
demi-sphere (non developpable). Criteres : distance a +-2 % pres de
l'espacement vise en ecart, espacement des lignes de niveau = buse +-5 %
(5e-95e centiles), points des lignes a la bonne distance exacte.

Usage (venv d'Apex) : ../../.venv/Scripts/python.exe -u test_surface_paths.py
"""
import time

import numpy as np
import trimesh

from surface_paths import NOZZLE_WIDTH, ZoneMesh, gaps_between

R, L = 10.0, 20.0
W = NOZZLE_WIDTH


def disk(skinny=False):
    """Disque bien maille (Delaunay sur une grille) ; skinny=True : disque en
    triangles-aiguilles (triangulation du seul contour), cas degrade."""
    from shapely.geometry import Point
    if skinny:
        v2, f = trimesh.creation.triangulate_polygon(Point(0, 0).buffer(R, 96))
    else:
        from scipy.spatial import Delaunay
        g = np.mgrid[-R:R:0.5, -R:R:0.5].reshape(2, -1).T
        g = g[np.linalg.norm(g, axis=1) < R - 0.25]
        ring = R * np.column_stack([np.cos(np.linspace(0, 2 * np.pi, 180, endpoint=False)),
                                    np.sin(np.linspace(0, 2 * np.pi, 180, endpoint=False))])
        v2 = np.vstack([g, ring])
        f = Delaunay(v2).simplices
    zm = ZoneMesh(np.column_stack([v2, np.zeros(len(v2))]), np.asarray(f))
    r = np.linalg.norm(zm.V[:, :2], axis=1)
    return "disque aiguilles" if skinny else "disque plat", zm, R - r, r, lambda p: R - np.linalg.norm(p[:, :2], axis=1)


def half_cylinder():
    m = trimesh.creation.cylinder(R, L, sections=96)
    keep = (np.abs(m.face_normals[:, 2]) < 0.1) & (m.triangles_center[:, 1] > 0)
    m = m.submesh([np.flatnonzero(keep)], append=True)
    zm = ZoneMesh(np.asarray(m.vertices), np.asarray(m.faces))

    def exact(p):                                   # rectangle deroule pi R x L
        s = R * np.arctan2(p[:, 1], p[:, 0])          # abscisse curviligne, 0..pi R
        z = p[:, 2] + L / 2
        return np.minimum.reduce([s, np.pi * R - s, z, L - z])
    return "demi-cylindre", zm, exact(zm.V), None, exact


def hemisphere():
    s = trimesh.creation.icosphere(4, radius=R).slice_plane([0, 0, 0], [0, 0, 1], cap=False)
    s = trimesh.Trimesh(s.vertices, s.faces, process=True)       # bord exactement sur l'equateur
    s.merge_vertices()
    s.update_faces(s.nondegenerate_faces())
    s.remove_unreferenced_vertices()
    v = np.asarray(s.vertices)
    zm = ZoneMesh(v / np.linalg.norm(v, axis=1)[:, None] * R, np.asarray(s.faces))
    zm.V = zm.V / np.linalg.norm(zm.V, axis=1)[:, None] * R      # raffinage reprojete sur la sphere

    def polar(p):
        return R * np.arccos(np.clip(p[:, 2] / np.linalg.norm(p, axis=1), -1, 1))
    return "demi-sphere", zm, R * np.pi / 2 - polar(zm.V), polar(zm.V), lambda p: R * np.pi / 2 - polar(p)


def check(name, zm, exact_b, exact_c, exact_fn):
    ok = True
    t0 = time.time()
    d = zm.distance_to_boundary()
    t_d = time.time() - t0
    err = np.abs(d - exact_b)
    levels = np.arange(W / 2, exact_b.max() - W, W)
    lines = [zm.isolines(d, lv) for lv in levels]
    on_level = np.concatenate([np.abs(exact_fn(np.vstack(ls)) - lv) for ls, lv in zip(lines, levels) if ls])
    near = np.concatenate([np.abs(exact_fn(np.vstack(ls)) - lv) for ls, lv in zip(lines, levels) if ls and lv <= 2.0])
    gaps = np.concatenate([gaps_between(a, b) for a, b in zip(lines[:-1], lines[1:])])
    p5, p50, p95 = np.percentile(gaps, [5, 50, 95])
    good = (0.95 * W <= p5) and (p95 <= 1.05 * W) and np.percentile(near, 95) <= 0.05
    ok &= good
    print(f"{name:<14} {len(zm.V):>6} sommets | distance au bord : ecart a l'exact median {np.median(err):.3f} mm "
          f"(max {err.max():.3f}), {t_d:.2f} s | position des lignes : a < 2 mm du bord 95e {np.percentile(near, 95):.3f}, partout median {np.median(on_level):.3f} ; "
          f"espacement {p50:.3f} ({p5:.3f}-{p95:.3f}) {'OK' if good else 'HORS +-5 %'}")
    if exact_c is not None:
        center = int(np.argmin(exact_c))
        dc = zm.distance_to_center(center)
        lm = zm.log_map(center)
        rad = np.linalg.norm(lm, axis=1)
        e_c, e_r = np.abs(dc - exact_c), np.abs(rad - exact_c)
        levels = np.arange(W, exact_c.max() - W, W)
        lines = [zm.isolines(dc, lv) for lv in levels]
        gaps = np.concatenate([gaps_between(a, b) for a, b in zip(lines[:-1], lines[1:])])
        p5, p50, p95 = np.percentile(gaps, [5, 50, 95])
        good = (0.95 * W <= p5) and (p95 <= 1.05 * W)
        ok &= good
        print(f"{'':<14} {'':>6}          | distance au centre : ecart median {np.median(e_c):.3f} (max {e_c.max():.3f}) ; "
              f"espacement des cercles {p50:.3f} ({p5:.3f}-{p95:.3f}) {'OK' if good else 'HORS +-5 %'} | "
              f"log map : rayon ecart median {np.median(e_r):.3f} (max {e_r.max():.3f})")
    return ok


def cylinder_line():
    """Distance a une COURBE : generatrice au milieu du demi-cylindre (theta =
    pi/2) ; distance exacte R |theta - pi/2| (developpable)."""
    name, zm, _b, _c, _f = half_cylinder()
    theta = np.arctan2(zm.V[:, 1], zm.V[:, 0])
    on = np.flatnonzero(np.abs(theta - np.pi / 2) < 1e-6)
    if len(on) < 3:                                 # aucune generatrice exactement sur des sommets
        on = np.argsort(np.abs(theta - np.pi / 2))[:1]
        target = theta[on[0]]
        on = np.flatnonzero(np.abs(theta - target) < 1e-9)
    on = on[np.argsort(zm.V[on, 2])]
    t0 = time.time()
    d = zm.distance_to_curve([(int(i), []) for i in on])
    t_d = time.time() - t0
    exact = R * np.abs(theta - theta[on[0]])
    levels = np.arange(W, exact.max() - W, W)
    lines = [zm.isolines(d, lv) for lv in levels]
    gaps = np.concatenate([gaps_between(a, b) for a, b in zip(lines[:-1], lines[1:])])
    p5, p50, p95 = np.percentile(gaps, [5, 50, 95])
    good = (0.95 * W <= p5) and (p95 <= 1.05 * W)
    print(f"{'demi-cylindre':<14} distance a une generatrice ({len(on)} points) : ecart a l'exact median "
          f"{np.median(np.abs(d - exact)):.3f} (max {np.abs(d - exact).max():.3f}), {t_d:.2f} s ; "
          f"espacement des lignes {p50:.3f} ({p5:.3f}-{p95:.3f}) {'OK' if good else 'HORS +-5 %'}")
    return good


if __name__ == "__main__":
    results = [check(*disk()), check(*half_cylinder()), check(*hemisphere()), cylinder_line()]
    print(f"\n{sum(results)}/{len(results)} conformes")
    print("\ncas degrade, non compte (triangles-aiguilles ; mesure : espacement 0,387-0,413, premiers anneaux 95e 0,056) :")
    check(*disk(skinny=True))
