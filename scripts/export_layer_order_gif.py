"""
scripts/export_layer_order_gif.py

Genere un GIF anime montrant les couches (Layer.contour_3d) s'ajouter UNE PAR
UNE dans leur ORDRE D'IMPRESSION REEL (result.colored_layers, sortie de
step5/step7 - la meme liste ordonnee que consomme step8), pour distinguer
visuellement un defaut d'ordre qu'un GLB statique ne montre pas.

Ne fait AUCUN calcul de pipeline propre : reutilise
scripts.export_slicing_preview.run_slicing_only (step1->step8, step9/step10
court-circuites) exactement comme file_io.visualization_export, puis ne fait
que du rendu matplotlib.

Usage :
    python scripts/export_layer_order_gif.py <mesh.stl> [-o sortie.gif]
        [--support-driven-order FRACTION] [--layer-stride N] [--fps N]
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib

matplotlib.use("Agg")  # rendu offscreen : pas de fenetre, requis en environnement non interactif

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection  # noqa: E402

from file_io.visualization_export import export_parts_glb  # noqa: E402
from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

# Meme degrade que file_io.visualization_export (limite deja actee avec
# l'utilisateur) : premier depose en bleu, dernier en jaune.
FIRST_COLOR = np.array([0x00, 0x00, 0xFF]) / 255.0
LAST_COLOR = np.array([0xFF, 0xF7, 0x00]) / 255.0
# Couche support_unresolved (Layer.support_unresolved, step5) : rouge vif,
# pour la distinguer immediatement du degrade normal - c'est le signal qu'on
# cherche a inspecter visuellement (une couche imprimee faute de support
# suffisant, repoussee en dernier recours par le glouton).
UNRESOLVED_COLOR = np.array([1.0, 0.1, 0.1])


def _layer_segments(layer) -> list[np.ndarray]:
    """Toutes les arêtes (paires de points consécutifs) des boucles de contour d'une couche."""
    segments = []
    for loop in layer.contour_3d:
        points = np.asarray(loop, dtype=float)
        if len(points) < 2:
            continue
        segments.extend(np.stack([points[:-1], points[1:]], axis=1))
    return segments


def _color_for_index(index: int, total: int) -> np.ndarray:
    if total <= 1:
        return FIRST_COLOR
    t = index / (total - 1)
    return FIRST_COLOR * (1.0 - t) + LAST_COLOR * t


def build_gif(
    layers: list,
    output_path: str,
    layer_stride: int = 1,
    fps: int = 12,
    elev: float = 20.0,
    azim: float = -60.0,
) -> None:
    selected = layers[::layer_stride] if layer_stride > 1 else layers
    total = len(selected)
    if total == 0:
        raise ValueError("Aucune couche a animer (liste vide).")

    all_points = np.vstack(
        [np.asarray(loop, dtype=float) for layer in selected for loop in layer.contour_3d if len(loop) > 0]
    )
    bbox_min, bbox_max = all_points.min(axis=0), all_points.max(axis=0)
    center = (bbox_min + bbox_max) / 2.0
    radius = max(float(np.max(bbox_max - bbox_min)) / 2.0, 1e-6) * 1.1

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("black")
    fig.patch.set_facecolor("black")

    frames = []
    accumulated_segments: list[np.ndarray] = []
    accumulated_colors: list[np.ndarray] = []

    for i, layer in enumerate(selected):
        color = UNRESOLVED_COLOR if getattr(layer, "support_unresolved", False) else _color_for_index(i, total)
        for segment in _layer_segments(layer):
            accumulated_segments.append(segment)
            accumulated_colors.append(color)

        ax.clear()
        ax.set_facecolor("black")
        ax.set_xlim(center[0] - radius, center[0] + radius)
        ax.set_ylim(center[1] - radius, center[1] + radius)
        ax.set_zlim(center[2] - radius, center[2] + radius)
        ax.set_box_aspect([1, 1, 1])
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(
            f"Couche {i + 1}/{total}", color="white", fontsize=12
        )

        collection = Line3DCollection(accumulated_segments, colors=accumulated_colors, linewidths=1.0)
        ax.add_collection3d(collection)

        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())
        frames.append(buf[:, :, :3].copy())

    plt.close(fig)

    from PIL import Image

    images = [Image.fromarray(frame) for frame in frames]
    duration_ms = int(1000 / max(fps, 1))
    images[0].save(
        output_path,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh_path")
    parser.add_argument("-o", "--output-gif", default="layer_order.gif")
    parser.add_argument("--resolution", type=float, default=None)
    parser.add_argument(
        "--layer-stride", type=int, default=1, help="1 = toutes les couches, N = une sur N"
    )
    parser.add_argument("--fps", type=int, default=12, help="images par seconde du GIF")
    parser.add_argument(
        "--support-driven-order", type=float, default=None, metavar="FRACTION",
        help="Ordonne les couches par support reel au lieu de l'isovaleur croissante "
             "(voir LayerParams.support_driven_order).",
    )
    parser.add_argument(
        "--split-disjoint-clusters", type=float, default=None, metavar="FACTEUR",
        help="Scinde une couche en plusieurs Layer quand ses contours forment des "
             "composantes spatialement disjointes (voir LayerParams). REQUIS pour "
             "que --frontier-order serve a quelque chose.",
    )
    parser.add_argument(
        "--frontier-order", action="store_true",
        help="Ordre par CROISSANCE DE FRONT : une couche n'est imprimable que si "
             "elle touche de la matiere deja imprimee (voir LayerParams.frontier_order). "
             "Sans --split-disjoint-clusters, sans effet.",
    )
    parser.add_argument(
        "--geodesic-order", action="store_true",
        help="Ordre d'impression dont la PROGRESSION vient de la distance "
             "geodesique au plateau (Dijkstra a travers la matiere), la continuite "
             "spatiale restant celle du parcours DFS. Voir LayerParams.geodesic_order.",
    )
    parser.add_argument(
        "--support-feedback", type=int, default=None, metavar="MAX_ITER",
        help="Boucle de retroaction de surplomb dans step4 : re-evalue le surplomb "
             "contre le champ REELLEMENT resolu et corrige les directions cibles "
             "fautives, jusqu'a stagnation du gain (plafond dur a 50). "
             "Voir FieldSolverParams.support_feedback_iterations.",
    )
    parser.add_argument(
        "--feedback-geodesic-target", action="store_true",
        help="La boucle de retroaction corrige les tetraedres fautifs avec leur "
             "direction GEODESIQUE LOCALE au lieu d'une rotation minimale "
             "(voir FieldSolverParams.support_feedback_use_geodesic_target). "
             "Sans --support-feedback, sans effet.",
    )
    parser.add_argument(
        "--geodesic-directions", action="store_true",
        help="Directions cibles issues de la distance geodesique au plateau "
             "(voir FieldSolverParams.use_geodesic_directions, limite #52).",
    )
    parser.add_argument(
        "--smoothness-modulation", type=float, default=None, metavar="STRENGTH",
        help="Module le poids de lissage PAR TETRAEDRE, abaisse la ou la piece se "
             "ramifie (|div(X)| du flux geodesique). Mesure : sans effet quand la "
             "direction cible est CONSTANTE sur le modele (une seule region step3, "
             "ce qui est le cas de real_overhang) - rien a departager entre. "
             "Voir FieldSolverParams.smoothness_modulation_strength.",
    )
    parser.add_argument(
        "--branch-cut-threshold", type=float, default=None, metavar="SEUIL",
        help="Segmentation step3 par STRUCTURE VOLUMIQUE plutot que par "
             "orientation de face : coupe le graphe d'adjacence tetraedre-"
             "tetraedre la ou |div(X)| depasse ce seuil, une composante = une "
             "branche = sa propre region. Corrige le defaut mesure sur "
             "real_overhang.stl (3 branches jamais separees par l'angle, a "
             "aucun seuil) - residu de planeite 1.1%% -> 0.0%%, inclinaison "
             "du pilier 50deg -> 0deg. Voir SegmentationParams.branch_cut_threshold. "
             "0.05 mesure comme separant proprement les 3 branches sur "
             "real_overhang (0.07 les refusionne deja).",
    )
    parser.add_argument(
        "--anchor-growth-weight", type=float, default=None, metavar="POIDS",
        help="Penalise les directions candidates opposees au sens de "
             "croissance depuis l'ancrage reel de la region (plateau ou "
             "region parent). Corrige la degenerescence mesuree du critere "
             "par defaut (risque residuel nul pour plusieurs candidats, "
             "argmin tranche par ordre d'indice - regions differentes "
             "convergent vers la MEME direction). Voir "
             "SegmentationParams.anchor_growth_weight. NON monotone : ne pas "
             "interpoler entre deux valeurs mesurees.",
    )
    parser.add_argument(
        "--transition-distance", type=float, default=None, metavar="MM",
        help="Distance (mm) sur laquelle les couches commencent a se deformer "
             "pour prendre l'orientation de la branche suivante. Diffuse le "
             "champ de DIRECTIONS CIBLES (auto-limite aux jonctions : a "
             "l'interieur d'une branche les directions sont constantes, les "
             "moyenner ne change rien). PETIT = cassure tardive et nette ; "
             "GRAND = transition precoce et progressive. Mesure sur "
             "real_overhang.stl (inclinaison mediane a z=+5) : OFF/2mm -> 0 "
             "deg, 20mm -> 8 deg, 40mm -> 17 deg. Voir "
             "FieldSolverParams.transition_distance_mm.",
    )
    parser.add_argument(
        "--export-parts-glb", metavar="CHEMIN.glb", default=None,
        help="OPTION DESACTIVABLE, INDEPENDANTE du GIF : exporte en plus un "
             "fichier GLB dedie montrant chaque PARTIE (region step3) du "
             "modele avec une fleche de sa propre couleur, partant du "
             "centroide de la partie, de la longueur de son etendue reelle "
             "le long de sa direction, et pointant dans sa direction "
             "d'impression (Region.target_direction). Voir "
             "file_io.visualization_export.export_parts_glb. Ne modifie en "
             "rien la generation du GIF ni du GLB de visualisation existant.",
    )
    parser.add_argument(
        "--layer-height", type=float, default=None, metavar="MM",
        help="Ecrase layer_height_base. A utiliser avec --geodesic-directions "
             "pour comparer A EPAISSEUR DE COUCHE REELLE EGALE : la geodesique "
             "produit des couches ~40%% plus epaisses a valeur nominale identique, "
             "et comparer sans compenser fausserait la lecture.",
    )
    parser.add_argument("--planar-parts", action="store_true",
                        help="parties a couches planes (step4c, 2026-10-19)")
    parser.add_argument("--reeb", type=float, default=None, metavar="MM",
                        help="parties planes decoupees par le squelette (graphe de Reeb geodesique, tranches de MM)")
    parser.add_argument("--merge-angle", type=float, default=None, metavar="DEG",
                        help="max_region_merge_angle_deg de step3 (decoupage plus fin)")
    parser.add_argument("--elev", type=float, default=20.0, help="Elevation de la camera (degres)")
    parser.add_argument("--azim", type=float, default=-60.0, help="Azimut de la camera (degres)")
    args = parser.parse_args(argv)

    params = PipelineParams()
    if args.resolution is not None:
        params.resolution = args.resolution
    if args.support_driven_order is not None:
        params.layers.support_driven_order = args.support_driven_order
        params.layers.nozzle_width = params.machine_profile.nozzle_diameter
        params.layers.support_rotation_axes = params.machine_profile.rotation_axes
    if args.split_disjoint_clusters is not None:
        params.layers.split_disjoint_clusters = args.split_disjoint_clusters
    if args.frontier_order:
        params.layers.frontier_order = True
    if args.geodesic_order:
        params.layers.geodesic_order = True
    if args.support_feedback is not None:
        params.field_solver = replace(
            params.field_solver,
            support_feedback_iterations=args.support_feedback,
            support_feedback_use_geodesic_target=args.feedback_geodesic_target,
        )
    if args.geodesic_directions:
        params.field_solver = replace(params.field_solver, use_geodesic_directions=True)
    if args.smoothness_modulation is not None:
        params.field_solver = replace(
            params.field_solver,
            smoothness_modulation_strength=args.smoothness_modulation,
        )
    if args.branch_cut_threshold is not None:
        params.segmentation = replace(
            params.segmentation, branch_cut_threshold=args.branch_cut_threshold
        )
    if args.anchor_growth_weight is not None:
        params.segmentation = replace(
            params.segmentation, anchor_growth_weight=args.anchor_growth_weight
        )
    if args.transition_distance is not None:
        params.field_solver = replace(
            params.field_solver, transition_distance_mm=args.transition_distance
        )
    if args.merge_angle is not None:
        params.segmentation = replace(params.segmentation, max_region_merge_angle_deg=args.merge_angle)
    if args.planar_parts:
        params.field_solver = replace(params.field_solver, planar_parts=True)
        params.collision_order = None
        if args.reeb is not None:
            params.field_solver = replace(params.field_solver, planar_parts_partition="reeb",
                                          planar_reeb_band_mm=args.reeb)
    if args.layer_height is not None:
        params.layers.layer_height_base = args.layer_height

    result = run_slicing_only(args.mesh_path, params)

    if args.export_parts_glb is not None:
        export_parts_glb(result, args.export_parts_glb)
        n_parts = len(result.region_dependency_graph.regions)
        print(f"[layer-gif] GLB des parties ecrit ({n_parts} partie(s)) : {args.export_parts_glb}")

    reversed_count = sum(
        1 for layer in result.colored_layers
        if getattr(layer, "deposition_reversed", False)
    )
    if reversed_count:
        # Signale explicitement ce que le GIF NE MONTRE PAS : l'orientation de
        # la tete. Deux animations identiques peuvent correspondre a une buse
        # pointant vers le bas ou vers le haut.
        print(f"[layer-gif] {reversed_count}/{len(result.colored_layers)} couche(s) "
              f"deposees a contre-sens (tete inversee) - NON VISIBLE sur le GIF")

    diag = getattr(result.scalar_field, "support_feedback_diagnostic", None)
    if diag is not None:
        curve = " -> ".join(f"{v:.2%}" for v in diag.volume_fraction_history)
        print(f"[layer-gif] retroaction : {diag.iterations_used} iteration(s), "
              f"converged={diag.converged}, stalled={diag.stalled}")
        print(f"[layer-gif]   surplomb auto-inflige : {curve}")

    print(f"[layer-gif] rendu de {len(result.colored_layers)} couche(s)...")
    build_gif(
        result.colored_layers,
        args.output_gif,
        layer_stride=args.layer_stride,
        fps=args.fps,
        elev=args.elev,
        azim=args.azim,
    )
    size = os.path.getsize(args.output_gif)
    print(f"[layer-gif] GIF ecrit : {os.path.abspath(args.output_gif)} ({size} octets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
