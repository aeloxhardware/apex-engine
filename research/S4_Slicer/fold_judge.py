"""
Juge des replis pour la boucle automatique de la methode Bird.

Definition de l'utilisateur (dessin du 2026-09-26) : une barre qui se
courbe en U est saine, meme si sa face interieure se comprime ; une barre
dont la SURFACE passe a travers elle-meme est un repli.

Mesure : les triangles de la surface deformee qui se traversent
(step1_ingest._find_self_intersecting_faces d'Apex) sont regroupes en zones
d'un seul tenant (triangles croises partageant une arete). Un plissement,
c'est beaucoup de petites zones eparpillees ; un repli, c'est une grande
zone ou un pan entier traverse un autre. Le critere est donc l'AIRE DE LA
PLUS GRANDE ZONE.

Criteres ecartes, tous mesures contre les verdicts de l'utilisateur :
tetraedres retournes et recouvrements entre tetraedres non voisins (voient
la compression normale des parois), enfoncement median des croisements (le
toit qui fusionne avec la cabine glisse, 0,089 mm, comme un plissement a
0,07 mm), distance geodesique entre triangles croises (ecart x4 seulement).

A lancer avec le venv d'APEX :
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe fold_judge.py <fichier.stl> [...]
"""
import os
import sys

import numpy as np
import scipy.sparse as sp
import trimesh
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from pipeline.step1_ingest import _find_self_intersecting_faces  # noqa: E402


def crossing_zones(mesh):
    """Zones de croisement de la surface, triees par aire decroissante :
    liste de (aire_mm2, etendue_mm, nombre_de_faces)."""
    bad = np.array(sorted(_find_self_intersecting_faces(mesh)), dtype=np.int64)
    if len(bad) == 0:
        return []
    index = -np.ones(len(mesh.faces), dtype=np.int64)
    index[bad] = np.arange(len(bad))
    adj = mesh.face_adjacency
    both = (index[adj] >= 0).all(axis=1)
    graph = sp.coo_matrix(
        (np.ones(both.sum()), (index[adj[both, 0]], index[adj[both, 1]])), shape=(len(bad), len(bad))
    )
    n_zones, label = connected_components(graph, directed=False)
    centers = mesh.triangles_center[bad]
    areas = np.bincount(label, weights=mesh.area_faces[bad], minlength=n_zones)
    counts = np.bincount(label, minlength=n_zones)
    zones = [(float(areas[z]), float(np.ptp(centers[label == z], axis=0).max()), int(counts[z]))
             for z in range(n_zones)]
    return sorted(zones, reverse=True)


if __name__ == "__main__":
    if sys.argv[1] == "--json":
        # Mode machine, appele par bird_loop.py depuis l'autre venv.
        import json
        zones = crossing_zones(trimesh.load(sys.argv[2]))
        print(json.dumps({
            "largest_zone_mm2": zones[0][0] if zones else 0.0,
            "largest_zone_extent_mm": zones[0][1] if zones else 0.0,
            "zones": len(zones),
            "total_mm2": sum(a for a, _, _ in zones),
        }))
        sys.exit(0)
    for path in sys.argv[1:]:
        zones = crossing_zones(trimesh.load(path))
        biggest = ", ".join(f"{a:.1f} mm2 ({e:.1f} mm)" for a, e, _ in zones[:3]) or "-"
        total = sum(a for a, _, _ in zones)
        print(f"{os.path.basename(path):<55} {len(zones):>4} zones | plus grandes : {biggest} | total {total:.0f} mm2")
