"""
Banc d'essai des 10 motifs de surface_patterns.py sur des zones planes.

Pour chaque motif et chaque forme (carre, disque, L concave, anneau troue) :
  - dedans    : aucun segment ne sort de la zone ;
  - couverture: part de la zone couverte par les cordons (largeur = spacing) ;
  - densite   : longueur tracee x spacing / aire (1 = ni trou ni surepaisseur) ;
  - passes    : nombre de polylignes (proxy des trajets a vide).
Seuils : dedans obligatoire, couverture >= 95 %, densite entre 0,85 et 1,3.

Usage (venv d'Apex) : ../../.venv/Scripts/python.exe test_surface_patterns.py
"""
import math

import numpy as np
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import unary_union

import surface_patterns as sp

SPACING = 0.4
SHAPES = {
    "carre 20": box(0, 0, 20, 20),
    "disque r10": Point(0, 0).buffer(10, quad_segs=64),
    "L concave": Polygon([(0, 0), (20, 0), (20, 6), (6, 6), (6, 20), (0, 20)]),
    "anneau": Point(0, 0).buffer(10, quad_segs=64).difference(Point(0, 0).buffer(4, quad_segs=64)),
}


def check(pattern, zone):
    paths = sp.generate(pattern, zone, SPACING)
    inside_zone = zone.buffer(1e-6)
    inside = all(inside_zone.covers(LineString(p)) for p in paths)
    length = sum(float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum()) for p in paths)
    beads = unary_union([LineString(p).buffer(SPACING / 2, cap_style="flat") for p in paths])
    # la zone utile est l'aire des AXES de cordon : un cordon deborde de
    # spacing/2 de part et d'autre, on mesure donc sur la zone elargie d'autant
    covered = beads.intersection(zone).area / zone.area
    density = length * SPACING / zone.area
    ok = inside and covered >= 0.95 and 0.85 <= density <= 1.3
    return ok, inside, covered, density, len(paths)


def main():
    failures = 0
    print(f"{'motif':<21} {'forme':<11} {'dedans':>6} {'couv.':>6} {'densite':>7} {'passes':>6}")
    for pattern in sp.PATTERNS:
        for name, zone in SHAPES.items():
            ok, inside, covered, density, n = check(pattern, zone)
            failures += not ok
            print(f"{pattern:<21} {name:<11} {'oui' if inside else 'NON':>6} {covered:>6.1%} "
                  f"{density:>7.2f} {n:>6}  {'' if ok else '<-- ECHEC'}")
    print(f"\n{len(sp.PATTERNS) * len(SHAPES) - failures}/{len(sp.PATTERNS) * len(SHAPES)} cas conformes")


if __name__ == "__main__":
    main()
