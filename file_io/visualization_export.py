"""
file_io/visualization_export.py

Export d'une scene GLB de visualisation a partir d'un PipelineResult
(main.run_pipeline) : mesh source + couches courbes + trajectoire finale
ajustee + marqueurs de collisions non resolues, en un seul fichier
ouvrable dans n'importe quel visualiseur 3D classique (Blender, etc.) -
meme format que le modele de tete (res/3d/toolheads/aetheris-4.glb),
choisi pour cette coherence plutot qu'introduire un second format.

Ce fichier ne fait AUCUN calcul de pipeline : il ne fait que RE-LIRE les
resultats deja produits par les 10 etapes et les assembler en scene.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import trimesh
from trimesh.path import Path3D
from trimesh.path.entities import Line

if TYPE_CHECKING:
    from main import PipelineResult

LAYER_COLORMAP_LOW = np.array([40, 80, 220, 255])   # premiere couche : bleu
LAYER_COLORMAP_HIGH = np.array([220, 60, 40, 255])  # derniere couche : rouge
EXTRUDING_COLOR = np.array([20, 200, 20, 255])      # segments avec extrusion > 0
TRAVEL_COLOR = np.array([200, 200, 200, 120])       # trajets a vide (extrusion == 0)
COLLISION_MARKER_COLOR = [230, 20, 20, 255]
COLLISION_MARKER_RADIUS = 1.5


def _layer_color(layer_index: int, layer_count: int) -> np.ndarray:
    """Degrade lineaire premiere -> derniere couche - rend visible l'ordre
    d'impression, pas une couleur arbitraire par couche."""
    if layer_count <= 1:
        return LAYER_COLORMAP_LOW
    t = layer_index / (layer_count - 1)
    return (LAYER_COLORMAP_LOW * (1.0 - t) + LAYER_COLORMAP_HIGH * t).astype(np.uint8)


def _layers_path3d(layers: list, layer_stride: int) -> Path3D | None:
    """
    Un Path3D regroupant TOUTES les couches (une ou plusieurs entites Line
    fermees par couche, cf. Layer.contour_3d = liste de polylignes fermees
    (K,3)), colore par index de couche. None si aucune couche n'a de
    contour exploitable.
    """
    vertices: list[np.ndarray] = []
    entities: list[Line] = []
    colors: list[list[int]] = []

    selected_layers = layers[::layer_stride] if layer_stride > 1 else layers
    layer_count = len(selected_layers)

    for layer_index, layer in enumerate(selected_layers):
        color = _layer_color(layer_index, layer_count).tolist()
        for loop in layer.contour_3d:
            if len(loop) < 2:
                continue
            start_idx = len(vertices)
            vertices.extend(list(loop))
            point_indices = list(range(start_idx, start_idx + len(loop)))
            entities.append(Line(points=point_indices))
            colors.append(color)

    if not entities:
        return None

    return Path3D(entities=entities, vertices=np.array(vertices), colors=colors)


def _toolpath_path3d(toolpath) -> Path3D | None:
    """
    Un Path3D reliant les points de trajectoire DANS L'ORDRE, colore par
    segment selon `extrusion` (segment i-1 -> i) : distingue visuellement
    les portions extrudees des trajets a vide, utile pour juger
    l'imprimabilite au premier coup d'oeil.
    """
    points = toolpath.points
    if len(points) < 2:
        return None

    vertices = np.array([p.position for p in points])
    entities = []
    colors = []
    for i in range(1, len(points)):
        entities.append(Line(points=[i - 1, i]))
        is_extruding = points[i].extrusion > 0.0 and points[i - 1].extrusion > 0.0
        colors.append((EXTRUDING_COLOR if is_extruding else TRAVEL_COLOR).tolist())

    return Path3D(entities=entities, vertices=vertices, colors=colors)


def _collision_markers(positions: list[np.ndarray]) -> trimesh.Trimesh | None:
    """Petites spheres rouges aux positions des collisions NON resolues -
    signal visuel immediat sans avoir a lire AdjustmentReport separement."""
    if not positions:
        return None

    spheres = []
    for position in positions:
        sphere = trimesh.creation.icosphere(subdivisions=1, radius=COLLISION_MARKER_RADIUS)
        sphere.apply_translation(np.asarray(position, dtype=float))
        sphere.visual.face_colors = COLLISION_MARKER_COLOR
        spheres.append(sphere)

    return trimesh.util.concatenate(spheres)


def export_visualization_glb(
    result: "PipelineResult",
    output_path: str,
    layer_stride: int = 1,
) -> None:
    """
    Assemble une scene GLB : mesh source (semi-transparent) + couches
    courbes (degrade par ordre d'impression) + trajectoire finale ajustee
    (verte = extrusion, grise = trajet a vide) + marqueurs rouges aux
    collisions non resolues (s'il y en a).

    Args:
        result: sortie de main.run_pipeline (contient toutes les etapes
            intermediaires, cette fonction ne recalcule rien).
        output_path: chemin du fichier .glb a ecrire.
        layer_stride: 1 = toutes les couches, N = une couche sur N (utile
            pour un modele a beaucoup de couches, evite une scene trop
            dense pour etre lisible).
    """
    scene = trimesh.Scene()

    source_mesh = trimesh.Trimesh(
        vertices=result.mesh_data.vertices, faces=result.mesh_data.faces, process=False
    )
    source_mesh.visual.face_colors = [200, 200, 200, 60]  # semi-transparent : voir les couches a travers
    scene.add_geometry(source_mesh, node_name="source_mesh", geom_name="source_mesh")

    layers_path = _layers_path3d(result.colored_layers, layer_stride)
    if layers_path is not None:
        scene.add_geometry(layers_path, node_name="layers", geom_name="layers")

    toolpath_path = _toolpath_path3d(result.adjusted_toolpath)
    if toolpath_path is not None:
        scene.add_geometry(toolpath_path, node_name="toolpath", geom_name="toolpath")

    # CollisionEvent ne porte pas de `position` propre (seulement
    # point_index) - on la retrouve via le point correspondant de la
    # trajectoire ajustee.
    unresolved_positions = [
        result.adjusted_toolpath.points[event.point_index].position
        for event in result.adjustment_report.unresolved_collisions
        if event.point_index < len(result.adjusted_toolpath.points)
    ]
    markers = _collision_markers(unresolved_positions)
    if markers is not None:
        scene.add_geometry(markers, node_name="unresolved_collisions", geom_name="unresolved_collisions")

    scene.export(output_path)
