"""Mesure : pymeshlab repare-t-il turkish.stl sans detruire la geometrie ?

Critere eliminatoire du plan (cascade mesuree). Protocole ITERATIF avec
ESCALADE :
  selection EXACTE des faces auto-intersectantes (pas de rayon arbitraire)
  -> retrait cible -> bouchage par vraie triangulation 3D -> re-selection,
  jusqu'a convergence (zero auto-intersection).

Deux mecanismes appris par la mesure, pas supposes :
 1. Une seule passe ne suffit pas : retirer un patch expose de nouvelles
    intersections (mesure : 73 puis 12 a la passe suivante).
 2. Le processus peut entrer en CYCLE LIMITE : le bouchage d'un trou
    recree exactement les faces intersectantes qu'on vient de retirer
    (mesure : 12 intersections identiques sur 7 passes consecutives). Quand
    le compte se repete, on DILATE la selection d'un anneau supplementaire
    avant de retirer, pour donner au bouchage la place de produire un patch
    non intersectant. C'est l'escalade prevue au plan : partir du minimum
    strict, n'elargir que sous contrainte mesuree.
"""
import sys
import time

import numpy as np
import trimesh

sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

MODEL = sys.argv[1] if len(sys.argv) > 1 else "res/test_models/turkish.stl"
MAX_PASSES = int(sys.argv[2]) if len(sys.argv) > 2 else 12
OUT = sys.argv[3] if len(sys.argv) > 3 else "out/_repaired_probe.npz"

raw = trimesh.load(MODEL)
cleaned, _report = s1._clean_and_repair(raw)
f0 = len(cleaned.faces)
z0 = (cleaned.vertices[:, 2].min(), cleaned.vertices[:, 2].max())
vol0 = cleaned.volume
print(f"[source] faces={f0} watertight={cleaned.is_watertight} "
      f"z=[{z0[0]:.2f},{z0[1]:.2f}] volume={vol0:.2f}")

import pymeshlab


def _close_holes(ms, strict=True):
    """MeshLab exige un maillage 2-manifold pour boucher : on repare la
    non-manifoldness AVANT chaque bouchage, et on n'echoue jamais dessus."""
    for op in (ms.meshing_remove_duplicate_faces,
               ms.meshing_remove_duplicate_vertices,
               ms.meshing_repair_non_manifold_edges):
        try:
            op()
        except Exception:
            pass
    try:
        ms.meshing_close_holes(maxholesize=300, selfintersection=strict,
                               newfaceselected=False)
        return True
    except Exception:
        return False


ms = pymeshlab.MeshSet()
ms.add_mesh(pymeshlab.Mesh(
    vertex_matrix=np.asarray(cleaned.vertices, dtype=np.float64),
    face_matrix=np.asarray(cleaned.faces, dtype=np.int32),
))

t0 = time.perf_counter()
dilations = 0
previous = None
converged = False
for it in range(MAX_PASSES):
    ms.compute_selection_by_self_intersections_per_face()
    n_sel = ms.current_mesh().selected_face_number()
    if n_sel == 0:
        print(f"[passe {it}] 0 auto-intersection -> CONVERGENCE")
        converged = True
        break
    if previous is not None and n_sel == previous:
        dilations += 1  # cycle limite detecte : elargir d'un anneau
    previous = n_sel
    for _ in range(dilations):
        ms.apply_selection_dilatation()
    n_removed = ms.current_mesh().selected_face_number()
    ms.meshing_remove_selected_faces()
    before = ms.current_mesh().face_number()
    _close_holes(ms, strict=True)
    after = ms.current_mesh().face_number()
    print(f"[passe {it}] {n_sel} intersections (dilat={dilations} -> "
          f"{n_removed} faces retirees), bouchage {before} -> {after}")
if not converged:
    print(f"[!] pas de convergence en {MAX_PASSES} passes")

_close_holes(ms, strict=False)
elapsed = time.perf_counter() - t0

m = ms.current_mesh()
repaired = trimesh.Trimesh(
    vertices=np.asarray(m.vertex_matrix(), dtype=np.float64),
    faces=np.asarray(m.face_matrix(), dtype=np.int64),
    process=False,
)
trimesh.repair.fix_normals(repaired)
trimesh.repair.fix_winding(repaired)

f1 = len(repaired.faces)
z1 = (repaired.vertices[:, 2].min(), repaired.vertices[:, 2].max())
vol1 = repaired.volume
n_open = int((np.bincount(
    trimesh.grouping.group_rows(repaired.edges_sorted, require_count=1).shape[:1]
    or [0]) if False else len(trimesh.grouping.group_rows(repaired.edges_sorted, require_count=1))))
print()
print(f"[resultat] faces {f0} -> {f1} ({100.0*f1/f0:.1f}% preserve) en {elapsed:.2f}s")
print(f"[resultat] z=[{z1[0]:.2f},{z1[1]:.2f}] (attendu [{z0[0]:.2f},{z0[1]:.2f}])")
print(f"[resultat] volume ecart {100.0*abs(vol1-vol0)/abs(vol0):.3f}%")
print(f"[resultat] watertight={repaired.is_watertight} aretes_ouvertes={n_open}")
# JAMAIS un export STL ici : le STL stocke des triangles bruts en
# float32 SANS sommets partages. Le round-trip re-fusionne par proximite
# et perd de la precision, ce qui RECREE les degenerescences de contact
# exact qu'on vient d'eliminer - mesure directe : 0 auto-intersection
# avant export, 51 apres rechargement du STL. Transmission sans perte
# obligatoire (float64, indices de faces conserves).
np.savez(OUT, vertices=np.asarray(repaired.vertices, dtype=np.float64),
         faces=np.asarray(repaired.faces, dtype=np.int64))
print(f"[resultat] exporte (sans perte) -> {OUT}")
