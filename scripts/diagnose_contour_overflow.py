"""
scripts/diagnose_contour_overflow.py

Diagnostic ad hoc : mesure combien de points de Layer.contour_3d (sortie
step5, AVANT tout offset/step6) tombent hors du solide source, sur
res/test_models/death.stl. Reutilise scripts/export_slicing_preview.py::
run_slicing_only telle quelle (pas de reimplementation de l'orchestration)
pour garantir les memes params/comportement que le GLB deja genere.

Contexte : artefact visuel rapporte par l'utilisateur sur out/death_preview.glb -
des contours de couche (vert->rouge, _layers_path3d) forment des pics
verticaux etroits qui depassent nettement de la silhouette du mesh, pres
d'un renflement/surplomb lateral (proche d'un cercle visible cote droit du
modele). CLAUDE.md documente (limite #31) que Layer.contour_3d brut est
"correct a 100%, 0 point hors solide" - mais cette mesure date de modeles
synthetiques ~30 triangles, jamais revalidee sur death.stl (12662 triangles,
reparation d'auto-intersections active, limite #43).

Usage :
    .venv/Scripts/python.exe scripts/diagnose_contour_overflow.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import trimesh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"
OUTSIDE_TOL_MM = 0.005  # abaisse par rapport a la limite #31 (0.05mm) - l'utilisateur
# rapporte un debordement FIN mais visible sur des zigzags serres pres d'un motif
# festonne dense (haut du modele), distinct du renflement lateral du diagnostic
# precedent - besoin de ne pas manquer un debordement sous 0.05mm


def main() -> None:
    params = PipelineParams()

    result = run_slicing_only(MODEL_PATH, params)
    # colored_layers (sortie step7) : c'est EXACTEMENT ce que le GLB affiche
    # (file_io/visualization_export.py ligne 186, _layers_path3d(result.colored_layers, ...)),
    # pas result.layers (step5+step6) - contour_3d n'est modifie ni par step6
    # ni par step7 en theorie (jamais reecrit, verifie par lecture de code),
    # mais on mesure sur la donnee REELLEMENT rendue pour ecarter tout doute.
    layers = result.colored_layers

    all_points: list[np.ndarray] = []
    layer_of_point: list[int] = []
    for layer in layers:
        for loop in layer.contour_3d:
            pts = np.asarray(loop)
            if len(pts) >= 2 and np.allclose(pts[0], pts[-1]):
                pts = pts[:-1]  # exclut le point de fermeture duplique
            all_points.append(pts)
            layer_of_point.extend([layer.id] * len(pts))

    all_points = np.concatenate(all_points, axis=0)
    layer_of_point = np.array(layer_of_point)

    # Mesure contre le mesh SOURCE original (pas le mesh interne
    # potentiellement repare par pymeshfix, limite #43) : on veut savoir si
    # le contour deborde du modele tel que fourni par l'utilisateur.
    # ProximityQuery uniquement (limite #15 : jamais mesh.contains(), qui
    # explose en RAM/temps sur ~165k points sans backend embree - mesure
    # directe : >4.6 Go RSS et >15 min sans terminer, process tue).
    source_mesh = trimesh.load(MODEL_PATH, force="mesh")
    prox = trimesh.proximity.ProximityQuery(source_mesh)

    closest, dist, triangle_id = prox.on_surface(all_points)
    # Signe du debordement via la normale de la face la plus proche : un
    # point est "dehors" si (point - closest) pointe globalement dans le
    # meme sens que la normale de la face la plus proche (produit scalaire
    # positif). Approximation locale suffisante pour ce diagnostic (pas
    # besoin d'un test d'inclusion global couteux).
    face_normals = source_mesh.face_normals[triangle_id]
    to_point = all_points - closest
    signed = np.einsum("ij,ij->i", to_point, face_normals)
    outside_mask = (dist > OUTSIDE_TOL_MM) & (signed > 0)
    n_outside = int(outside_mask.sum())
    n_total = len(all_points)

    # Distribution complete des distances signees (dehors uniquement), meme
    # si n_outside=0 a un seuil donne - un debordement fin mais reel doit
    # apparaitre ici, pas seulement dans le compte binaire.
    outside_only_dist = dist[signed > 0]
    if len(outside_only_dist) > 0:
        percentiles = [50, 90, 99, 99.9, 100]
        pvals = np.percentile(outside_only_dist, percentiles)
        print("Distribution des distances (points du cote EXTERIEUR de la face "
              "la plus proche, quel que soit le seuil) :")
        for p, v in zip(percentiles, pvals):
            print(f"  p{p:<5} = {v:.5f}mm")
        print(f"  max = {outside_only_dist.max():.5f}mm sur {len(outside_only_dist)}/{n_total} points "
              f"du cote exterieur (avant tout seuil)")

    print(f"\n=== RESULTAT : {n_outside}/{n_total} points de contour_3d hors solide "
          f"(tol={OUTSIDE_TOL_MM}mm) ===")

    if n_outside == 0:
        print("Aucun point hors solide sur contour_3d brut a ce seuil : l'artefact "
              "visuel ne vient PAS de la geometrie de step5 (Layer.contour_3d), "
              "meme a une tolerance fine. Voir la distribution ci-dessus pour la "
              "distance reelle max mesuree. A chercher cote export GLB "
              "(file_io/visualization_export.py, ex: Path3D/Line qui relie des "
              "boucles DIFFERENTES entre elles) ou rendu du visualiseur.")
        return

    outside_points = all_points[outside_mask]
    outside_dist = dist[outside_mask]
    outside_layers = layer_of_point[outside_mask]

    order = np.argsort(-outside_dist)
    print("\nTop 30 points les plus eloignes du solide (x, y, z, dist_mm, layer_id):")
    for i in order[:30]:
        p = outside_points[i]
        print(f"  ({p[0]:8.3f}, {p[1]:8.3f}, {p[2]:8.3f})  dist={outside_dist[i]:.4f}mm  layer={outside_layers[i]}")

    centroid = outside_points.mean(axis=0)
    bbox_min = outside_points.min(axis=0)
    bbox_max = outside_points.max(axis=0)
    print(f"\nCentroide des points hors-solide : {centroid}")
    print(f"Bounding box : min={bbox_min}  max={bbox_max}")

    # Combien de couches distinctes sont touchees, et sur quelle plage
    # d'isovaleur/altitude Z - utile pour juger si ca correspond a UNE zone
    # localisee (surplomb) ou un probleme diffus.
    touched_layers = np.unique(outside_layers)
    print(f"\n{len(touched_layers)} couche(s) distincte(s) touchee(s) sur {len(layers)} : "
          f"ids {touched_layers.tolist()[:20]}{'...' if len(touched_layers) > 20 else ''}")

    r_xy = np.linalg.norm(outside_points[:, :2], axis=1)
    print(f"Rayon XY des points hors-solide : min={r_xy.min():.2f} max={r_xy.max():.2f} "
          f"median={np.median(r_xy):.2f}")
    all_r_xy = np.linalg.norm(source_mesh.vertices[:, :2], axis=1)
    print(f"Rayon XY du mesh source        : min={all_r_xy.min():.2f} max={all_r_xy.max():.2f} "
          f"median={np.median(all_r_xy):.2f}")


if __name__ == "__main__":
    main()
