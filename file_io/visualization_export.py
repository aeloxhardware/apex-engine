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
from scipy.spatial import cKDTree
import trimesh
from trimesh.path import Path3D
from trimesh.path.entities import Line

if TYPE_CHECKING:
    from main import PipelineResult

LAYER_COLORMAP_LOW = np.array([40, 80, 220, 255])   # premiere couche : bleu
LAYER_COLORMAP_HIGH = np.array([220, 60, 40, 255])  # derniere couche : rouge
EXTRUDING_COLOR = np.array([20, 200, 20, 255])      # segments avec extrusion > 0
# Trajets a vide (extrusion == 0) : bleu vif OPAQUE. Volontairement plus
# clair/cyan que LAYER_COLORMAP_LOW ci-dessus (le bleu de la premiere
# couche), pour rester distinguable a l'oeil malgre la proximite de teinte.
# Opacite pleine (et non 120 comme a l'origine) : ces trajets sont ce qu'on
# cherche a inspecter (liaisons inter-branches tracees en ligne droite, sans
# planification - voir step8._iter_schedule_points), pas un decor de fond.
TRAVEL_COLOR = np.array([0, 150, 255, 255])
COLLISION_MARKER_COLOR = [230, 20, 20, 255]
COLLISION_MARKER_RADIUS = 1.5
UP_ARROW_COLOR = [0, 220, 0, 255]

# Les 3 premieres parties (l'immense majorite des modeles de ce projet -
# real_overhang.stl en a 3) recoivent rouge/vert/bleu, demandes
# explicitement pour leur lisibilite immediate. Au-dela, des couleurs sont
# generees (voir `_part_color`), ALEATOIRES mais contraintes a rester
# ELOIGNEES de toutes les couleurs deja choisies - jamais une palette fixe
# qui finirait par cycler et faire porter la meme couleur a deux parties
# distinctes sur un modele a beaucoup de branches.
PART_ARROW_COLORS = np.array([
    [230, 30, 30, 255],   # rouge
    [30, 200, 60, 255],   # vert
    [40, 100, 255, 255],  # bleu
])


def _part_color(index: int, previous_colors: list[np.ndarray]) -> np.ndarray:
    """
    Couleur de la partie d'indice `index` (0-based) : rouge/vert/bleu pour
    les 3 premieres, puis un tirage ALEATOIRE dans l'espace HSV (teinte
    uniforme sur tout le cercle, saturation et valeur elevees pour rester
    vive) rejete et retire tant qu'il reste trop proche - en distance
    circulaire de TEINTE, pas RGB brute (RGB confondrait par exemple un
    rouge fonce et un rouge clair comme "eloignes" alors qu'ils se
    ressemblent visuellement) - d'une couleur deja retenue.

    `previous_colors` : les couleurs deja assignees aux parties precedentes
    (RGBA uint8), pour que chaque nouveau tirage les evite toutes, pas
    seulement la derniere.
    """
    if index < len(PART_ARROW_COLORS):
        return PART_ARROW_COLORS[index]

    import colorsys

    def hue_of(rgba: np.ndarray) -> float:
        h, _s, _v = colorsys.rgb_to_hsv(*(rgba[:3].astype(np.float64) / 255.0))
        return h

    previous_hues = [hue_of(c) for c in previous_colors]
    rng = np.random.default_rng(index)  # deterministe : meme modele -> memes couleurs

    min_hue_distance = 0.12  # sur un cercle de teinte unitaire (~43 deg)
    for _attempt in range(50):
        hue = float(rng.uniform(0.0, 1.0))
        far_enough = all(
            min(abs(hue - h), 1.0 - abs(hue - h)) >= min_hue_distance
            for h in previous_hues
        )
        if far_enough or not previous_hues:
            break
    r, g, b = colorsys.hsv_to_rgb(hue, 0.75, 0.95)
    return np.array([int(r * 255), int(g * 255), int(b * 255), 255], dtype=np.uint8)

# Coloration du MESH SOURCE par ordre d'impression (option
# `color_by_print_order`) : bleu = imprime en premier, jaune = en dernier.
# Deliberement distinct du degrade des COUCHES (bleu -> rouge ci-dessus) :
# les deux peuvent etre affiches ensemble, une palette commune les rendrait
# indistinguables.
PRINT_ORDER_FIRST_COLOR = np.array([0, 0, 255, 255])      # premier : #0000ff
PRINT_ORDER_LAST_COLOR = np.array([255, 247, 0, 255])     # dernier : #fff700
# Un mesh colore par ordre d'impression doit rester LISIBLE : le defaut
# semi-transparent (alpha 60, pour voir les couches au travers) delave les
# couleurs au point de les rendre inexploitables. L'option bascule donc
# aussi l'opacite - c'est le but meme de l'option, inspecter le modele.
PRINT_ORDER_ALPHA = 235


def _print_order_face_colors(mesh_data, layers) -> "np.ndarray | None":
    """
    Une couleur par face du mesh source, du JAUNE (matiere deposee en
    premier) au BLEU (deposee en dernier), suivant l'ORDRE D'IMPRESSION
    REEL des couches - `layers` etant la liste deja ordonnee que step8 a
    consommee (sortie step7), pas une reconstruction.

    Chaque face est rattachee a la couche dont un point de contour est le
    plus proche de son centre (cKDTree sur tous les points de contour,
    chacun etiquete par l'index de sa couche dans l'ordre d'impression),
    puis coloree par cet index normalise. Le degrade est donc CONTINU et
    aussi fin que le nombre de couches (~200 sur death.stl), la ou une
    coloration par profondeur de region ne donnerait que quelques aplats.

    Rattachement par plus proche voisin plutot qu'une evaluation du champ
    scalaire : `phi` seul n'est pas l'ordre d'impression final des que le
    modele a plusieurs regions ou branches (step5 reordonne les couches par
    rang topologique puis continuite spatiale). Les couches ordonnees sont
    la seule source de verite sur l'ordre reellement imprime.

    Retourne None si aucune couche n'a de contour exploitable : l'appelant
    retombe alors sur le mesh gris semi-transparent, jamais une exception
    ni une couleur inventee.
    """
    if not layers:
        return None

    points: list[np.ndarray] = []
    point_layer_index: list[int] = []
    for layer_index, layer in enumerate(layers):
        for loop in layer.contour_3d:
            if len(loop) == 0:
                continue
            loop_array = np.asarray(loop, dtype=float)
            points.append(loop_array)
            point_layer_index.extend([layer_index] * len(loop_array))

    if not points:
        return None

    contour_points = np.vstack(points)
    layer_of_point = np.asarray(point_layer_index, dtype=np.int64)

    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    _distances, nearest = cKDTree(contour_points).query(face_centers)
    face_layer = layer_of_point[np.asarray(nearest, dtype=np.int64)]

    span = max(int(face_layer.max()) - int(face_layer.min()), 1)
    ratio = np.clip((face_layer - face_layer.min()) / span, 0.0, 1.0)

    colors = (
        PRINT_ORDER_FIRST_COLOR[None, :] * (1.0 - ratio[:, None])
        + PRINT_ORDER_LAST_COLOR[None, :] * ratio[:, None]
    ).astype(np.uint8)
    colors[:, 3] = PRINT_ORDER_ALPHA
    return colors


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


def _up_arrow(mesh_data) -> trimesh.Trimesh:
    """
    Fleche verte pointant vers +Z (convention deja actee partout dans ce
    fichier : le plateau est a z=platform_z, l'impression progresse en Z
    croissant) - repere visuel pour s'orienter dans la scene GLB, demande
    explicitement par l'utilisateur. Taille proportionnelle a l'etendue du
    mesh source (pas une taille fixe, resterait illisible ou disproportionnee
    selon l'echelle du modele). Positionnee a cote de la bounding box, pas en
    collision visuelle avec le modele.
    """
    bbox_min = mesh_data.vertices.min(axis=0)
    bbox_max = mesh_data.vertices.max(axis=0)
    extent = bbox_max - bbox_min
    scale = float(np.max(extent)) if np.max(extent) > 0 else 1.0

    shaft_radius = scale * 0.015
    shaft_height = scale * 0.6
    head_radius = shaft_radius * 2.5
    head_height = scale * 0.15

    shaft = trimesh.creation.cylinder(radius=shaft_radius, height=shaft_height)
    shaft.apply_translation([0.0, 0.0, shaft_height / 2.0])

    head = trimesh.creation.cone(radius=head_radius, height=head_height)
    head.apply_translation([0.0, 0.0, shaft_height])

    arrow = trimesh.util.concatenate([shaft, head])
    arrow.visual.face_colors = UP_ARROW_COLOR

    origin_offset = np.array(
        [bbox_min[0] - scale * 0.2, bbox_min[1] - scale * 0.2, bbox_min[2]]
    )
    arrow.apply_translation(origin_offset)
    return arrow


def export_visualization_glb(
    result: "PipelineResult",
    output_path: str,
    layer_stride: int = 1,
    color_by_print_order: bool = True,
    show_toolpath: bool = False,
) -> None:
    """
    Assemble une scene GLB : mesh source (degrade par ordre d'impression
    par defaut) + couches courbes (degrade par ordre d'impression) +
    marqueurs rouges aux collisions non resolues (s'il y en a). La
    trajectoire finale ajustee est omise par defaut (voir show_toolpath).

    Args:
        result: sortie de main.run_pipeline (contient toutes les etapes
            intermediaires, cette fonction ne recalcule rien).
        output_path: chemin du fichier .glb a ecrire.
        layer_stride: 1 = toutes les couches, N = une couche sur N (utile
            pour un modele a beaucoup de couches, evite une scene trop
            dense pour etre lisible).
        color_by_print_order: True (defaut) = mesh colore en degrade
            continu suivant l'ordre d'impression reel des couches (BLEU
            #0000ff depose en premier -> JAUNE #fff700 en dernier). False
            = mesh source gris semi-transparent (comportement historique).
            Les COUCHES gardent leur propre degrade (bleu -> rouge) : seul
            le mesh change. Rend le mesh quasi opaque : c'est le but de
            l'option (lire l'ordre sur le modele lui-meme).
        show_toolpath: False (defaut) = trajectoire omise de la scene,
            desactivee temporairement a la demande de l'utilisateur. True
            = comportement historique (Path3D vert/bleu par extrusion).
    """
    scene = trimesh.Scene()

    source_mesh = trimesh.Trimesh(
        vertices=result.mesh_data.vertices, faces=result.mesh_data.faces, process=False
    )
    print_order_colors = (
        _print_order_face_colors(result.mesh_data, result.colored_layers)
        if color_by_print_order
        else None
    )
    if print_order_colors is not None:
        # glTF stocke les couleurs PAR SOMMET, et les sommets d'un mesh
        # ferme sont partages entre faces voisines : assigner des
        # `face_colors` tels quels fait ecraser toutes les couleurs par une
        # seule ou moyennee a l'export - mesure directe sur le modele en T
        # (5 regions) : 21 couleurs delavees dans le .glb relu, aucun jaune
        # ni bleu pur, les frontieres entre regions brouillees par la
        # moyenne aux sommets partages. Dedoubler les sommets
        # (`unmerge_vertices`) donne a chaque face les siens, ce qui
        # preserve exactement la coloration par region. Uniquement sur ce
        # chemin : le mesh gris uniforme par defaut n'en a pas besoin.
        source_mesh.unmerge_vertices()
        source_mesh.visual.face_colors = print_order_colors
    else:
        source_mesh.visual.face_colors = [200, 200, 200, 60]  # semi-transparent : voir les couches a travers
    scene.add_geometry(source_mesh, node_name="source_mesh", geom_name="source_mesh")

    up_arrow = _up_arrow(result.mesh_data)
    scene.add_geometry(up_arrow, node_name="up_arrow", geom_name="up_arrow")

    layers_path = _layers_path3d(result.colored_layers, layer_stride)
    if layers_path is not None:
        scene.add_geometry(layers_path, node_name="layers", geom_name="layers")

    if show_toolpath:
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


def _direction_arrow(
    origin: np.ndarray, direction: np.ndarray, length: float, color: np.ndarray
) -> trimesh.Trimesh:
    """
    Fleche generique pointant depuis `origin` selon `direction` (unitaire),
    de longueur `length`. Meme construction que `_up_arrow` (cylindre + cone),
    mais orientation et longueur PARAMETRIQUES au lieu de fixees a +Z - c'est
    tout ce qui differe, aucune duplication de la geometrie de base.
    """
    shaft_radius = max(length * 0.04, 1e-6)
    shaft_height = length * 0.8
    head_radius = shaft_radius * 2.5
    head_height = length * 0.2

    shaft = trimesh.creation.cylinder(radius=shaft_radius, height=shaft_height)
    shaft.apply_translation([0.0, 0.0, shaft_height / 2.0])

    head = trimesh.creation.cone(radius=head_radius, height=head_height)
    head.apply_translation([0.0, 0.0, shaft_height])

    arrow = trimesh.util.concatenate([shaft, head])
    arrow.visual.face_colors = color

    # Rotation de +Z (repere local de construction ci-dessus) vers
    # `direction` : axe = Z x direction, angle = acos(Z . direction).
    z_axis = np.array([0.0, 0.0, 1.0])
    dot = float(np.clip(np.dot(z_axis, direction), -1.0, 1.0))
    if dot > 1.0 - 1e-9:
        rotation = np.eye(4)
    elif dot < -1.0 + 1e-9:
        # Direction opposee a Z : rotation de 180 deg autour de n'importe
        # quel axe perpendiculaire (X convient toujours ici).
        rotation = trimesh.transformations.rotation_matrix(np.pi, [1.0, 0.0, 0.0])
    else:
        axis = np.cross(z_axis, direction)
        axis = axis / np.linalg.norm(axis)
        angle = float(np.arccos(dot))
        rotation = trimesh.transformations.rotation_matrix(angle, axis)

    arrow.apply_transform(rotation)
    arrow.apply_translation(origin)
    return arrow


def export_parts_glb(result: "PipelineResult", output_path: str) -> None:
    """
    Scene GLB DEDIEE aux PARTIES du modele (regions de
    `RegionDependencyGraph`, telles que produites par step3 - voir
    `SegmentationParams.branch_cut_threshold` pour la segmentation par
    structure volumique plutot que par orientation de face) : le mesh
    source COLORE PAR PARTIE (chaque face prend la couleur de sa region,
    meme palette que les fleches - une face sans region assignee reste gris
    semi-transparent) + UNE FLECHE PAR PARTIE, partant du centroide de la
    partie, pointant dans sa direction d'impression
    (`Region.target_direction`), et LONGUE DE L'ETENDUE REELLE de la partie
    projetee sur cette direction - pas une taille arbitraire, pour que la
    fleche donne une lecture immediate de quelle portion du modele chaque
    direction gouverne.

    Fichier SEPARE de `export_visualization_glb` (mecanisme independant,
    jamais mele a l'export existant) : sert a DIAGNOSTIQUER la segmentation
    elle-meme (combien de parties, ou, dans quel sens), pas a visualiser le
    resultat du slicing. Ne fait AUCUN calcul de pipeline, comme le reste de
    ce fichier - relit uniquement `result.region_dependency_graph` et
    `result.mesh_data`.

    Args:
        result: sortie de main.run_pipeline (ou tout objet portant au moins
            `mesh_data` et `region_dependency_graph`, meme forme que
            PipelineResult).
        output_path: chemin du fichier .glb a ecrire.
    """
    scene = trimesh.Scene()

    source_mesh = trimesh.Trimesh(
        vertices=result.mesh_data.vertices, faces=result.mesh_data.faces, process=False
    )

    # Une couleur par partie, calculee UNE FOIS (voir `_part_color` :
    # rouge/vert/bleu pour les 3 premieres, puis aleatoire mais eloigne de
    # tout ce qui precede) et reutilisee identique pour la face ET la
    # fleche de chaque partie - jamais deux tirages independants qui
    # pourraient diverger.
    regions = result.region_dependency_graph.regions
    sorted_regions = sorted(regions.items())
    part_colors: list[np.ndarray] = []
    for i, _region_id in enumerate([rid for rid, _r in sorted_regions]):
        part_colors.append(_part_color(i, part_colors))

    # Couleur PAR FACE selon sa region - gris semi-transparent par defaut
    # pour toute face qu'aucune region ne revendique (ne devrait pas arriver
    # en pratique, step3 couvre toutes les faces, mais explicite plutot que
    # suppose).
    face_colors = np.tile(
        np.array([200, 200, 200, 60], dtype=np.uint8), (result.mesh_data.face_count, 1)
    )
    for i, (region_id, region) in enumerate(sorted_regions):
        face_ids = np.asarray(region.face_ids)
        if len(face_ids) == 0:
            continue
        face_colors[face_ids] = part_colors[i]

    # Meme precaution que `export_visualization_glb` (voir sa docstring) :
    # glTF stocke les couleurs PAR SOMMET, et des sommets partages entre
    # faces de regions differentes seraient moyennes/ecrases a l'export sans
    # dedoublement prealable.
    source_mesh.unmerge_vertices()
    source_mesh.visual.face_colors = face_colors
    scene.add_geometry(source_mesh, node_name="source_mesh", geom_name="source_mesh")

    for i, (region_id, region) in enumerate(sorted_regions):
        face_ids = np.asarray(region.face_ids)
        if len(face_ids) == 0:
            continue
        vertex_ids = np.unique(result.mesh_data.faces[face_ids])
        points = result.mesh_data.vertices[vertex_ids]

        direction = np.asarray(region.target_direction, dtype=float)
        norm = np.linalg.norm(direction)
        if norm < 1e-9:
            continue
        direction = direction / norm

        centroid = points.mean(axis=0)
        # Longueur = etendue REELLE de la partie le long de sa propre
        # direction (projection des sommets sur `direction`, pas la
        # diagonale de la bounding box - une partie fine et longue
        # perpendiculairement a sa direction ne doit pas gonfler la fleche).
        projections = (points - centroid) @ direction
        length = float(projections.max() - projections.min())
        if length < 1e-6:
            length = float(np.ptp(points, axis=0).max()) or 1.0

        color = part_colors[i]
        # Centre la fleche sur la partie : origine reculee d'une demi-longueur
        # depuis le centroide, dans le sens OPPOSE a la direction, pour que la
        # fleche s'etende de part et d'autre du centroide plutot que de
        # partir de lui vers l'exterieur uniquement.
        arrow_origin = centroid - direction * (length / 2.0)
        arrow = _direction_arrow(arrow_origin, direction, length, color)
        scene.add_geometry(
            arrow, node_name=f"part_{region_id}_arrow", geom_name=f"part_{region_id}_arrow"
        )

    scene.export(output_path)
