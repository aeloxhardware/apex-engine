"""
pipeline/step6_toolpath2d.py

Etape 6 du pipeline : remplissage 2D de chaque couche (parois + infill),
dans le repere local de la couche produit par step5.

Depend de Layer (step5, contour_3d + local_frame), MeshData (step1, kdtree
pour le pont surface/volume), OverhangMap (step2, is_bridge +
bridge_span_direction).

Ce que ce fichier fait, dans l'ordre, par Layer :
  1. Projette chaque boucle de contour_3d dans le plan local (local_frame).
  2. Classe les boucles projetees en contours "matiere" / trous par
     imbrication geometrique (point_in_polygon_2d), PAS par signe d'aire -
     voir CLAUDE.md limite #11 et la Decision 3bis de step5.
  3. Genere les parois par offset polygonal successif vers l'interieur
     (shapely.buffer), wall_count niveaux, en s'arretant proprement si
     l'offset fait disparaitre une zone trop fine avant d'atteindre
     wall_count.
  4. Genere l'infill a l'interieur du dernier offset.
  5. Pour les boucles associees a une zone is_bridge (step2, via le meme
     pont kdtree surface/volume que step4/step5) : ajoute des parois
     d'ancrage supplementaires sur le perimetre puis reoriente l'infill
     pour traverser perpendiculairement a bridge_span_direction.

Sortie : la meme liste de Layer, avec toolpath_2d rempli (repere local 2D,
pas 3D - step8 se chargera de le relever au repere machine).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum

import os

import numpy as np
import shapely

import shapely.affinity
from shapely.geometry import LineString, MultiPolygon, Polygon

from core import _cpp
from core import progress as _progress
from core.geometry_utils import point_in_polygon_2d
from core.mesh_types import Layer, MeshData, OverhangMap

from pipeline.step6_params import (  # noqa: F401  (réexport : imports existants)
    InfillPattern,
    Step6Params,
    ToolpathGenerationError,
)
from pipeline.step6_geometry import (  # noqa: F401  (réexport : imports existants)
    _COINCIDENT_EPS,
    _classify_contour_indices,
    _classify_contours,
    _frame_normal,
    _interpolate_out_of_plane,
    _lift_to_3d,
    _loop_matches,
    _offset_walls,
    _out_of_plane_for_polyline,
    _out_of_plane_offset,
    _polygon_pieces,
    _project_to_local_2d,
    _signed_area_2d,
    _to_shapely_polygon,
    _to_toolpath_3col,
)
from pipeline.step6_infill import (  # noqa: F401  (réexport : imports existants)
    _MAX_INFILL_LINES,
    _SPIRAL_DENSITY_ITERATIONS,
    _SPIRAL_DENSITY_TOLERANCE,
    _SPIRAL_MIN_RADIUS_PITCHES,
    _SPIRAL_SAMPLE_STEP_FACTOR,
    _SPIRAL_STOP_PROBABILITY,
    _ZIGZAG_MAX_CONNECTOR_SPACINGS,
    _angle_nondecreasing,
    _boundary_connector,
    _clipped_scanlines,
    _connect_zigzag,
    _generate_infill,
    _generate_spiral_infill,
    _infill_angle_for_layer,
    _layer_rank,
    _oscillating_spiral_polyline,
    _plate_axis_in_local_2d,
)
from pipeline.step6_bridge_runs import (  # noqa: F401  (réexport : imports existants)
    _bridge_direction_for_loop,
    _is_closed_ring,
    _merge_short_runs,
    _monotonic_runs,
    _reorient_bridge_toolpath,
    _split_ring_into_runs,
    _split_rings_in_toolpath,
)


# ---------------------------------------------------------------------------
# Decision 5 - Decoupage des anneaux en troncons MONOTONES
#
# CAUSE RACINE MESUREE : sur death.stl, 70% des defauts d'ordre residuels ont
# leur support a EXACTEMENT 2 positions plus loin dans la trajectoire, avec un
# ecart horizontal/vertical minuscule (< 0.5mm / < 0.03mm) - donc DEUX SOMMETS
# CONSECUTIFS D'UN MEME ANNEAU. Sur une couche dont le contour ondule hors du
# plan local (cf. out_of_plane_interpolation_power), un anneau ferme DOIT a la
# fois monter et descendre pour revenir a son point de depart, et le sens de
# parcours de shapely est arbitraire vis-a-vis de cette ondulation : un point
# depose pendant une DESCENTE locale n'a rien en dessous de lui au moment T
# (son predecesseur immediat dans le parcours est plus HAUT).
#
# Ce que fait CE fichier, et uniquement ca : decoupe chaque anneau ferme en
# troncons MONOTONES en offset hors-plan, chacun stocke dans le sens MONTANT
# (troncon descendant parcouru en sens inverse - memes segments physiques).
# Un troncon montant est auto-portant : chaque point s'appuie sur celui qui
# vient d'etre depose juste avant lui, plus bas et proche horizontalement.
#
# Chaque troncon devient une POLYLIGNE SEPAREE de Layer.toolpath_2d, jamais
# une concatenation. C'est un bug reel deja mesure : step8._iter_schedule_points
# n'insere un trajet a vide (extrusion=0.0) QUE lorsque la cle
# (layer_id, segment_index) change, jamais a l'interieur d'une polyligne. Une
# version precedente concatenait les troncons : la buse extrudait donc en ligne
# droite entre la fin d'un troncon et le debut du suivant (cordon parasite en
# travers de la piece), et step9_collision_sim._grow_printed_index_to - qui
# n'ajoute une capsule que si les DEUX extremites ont extrusion > 0.0 -
# enregistrait ces cordons parasites comme de la matiere REELLE, faussant a la
# hausse toute mesure de support en aval.
#
# L'ORDRE des troncons entre eux n'est deliberement PAS decide ici : c'est le
# travail de pipeline/step7b_deposition_order.py, qui dispose de l'etat
# incremental 3D reel (toutes couches, toutes parois, tout infill) et des
# directions d'inclinaison reellement atteignables, la ou ce fichier ne voit
# qu'un seul anneau dans le repere local 2D d'une seule couche.
# ---------------------------------------------------------------------------


# Nombre de fils de step6 (None = min(8, coeurs)) ; 1 = serie historique.
_STEP6_THREADS = None

# Pas 15 (limite #56) : parois et remplissage d'une zone par le noyau natif
# (Clipper2 au lieu de shapely/GEOS). PAS au bit pres - meme offset en
# onglet, autre arithmetique et autre point de depart d'anneau -, d'ou
# l'option. ACTIVE PAR DEFAUT quand le noyau natif est actif (2026-10-26) :
# `Step6Params.native_fill=False` rend shapely (limite #82). Les ponts restent en Python.
_NATIVE_PATTERN_CODES = {InfillPattern.LINES: 0, InfillPattern.GRID: 1, InfillPattern.CROSSED_ZIGZAG: 2}


def _native_fill_module(params: "Step6Params"):
    if not params.native_fill:
        return None
    from core import _cpp
    module = _cpp.module() if _cpp.enabled() else None
    return module if module is not None and hasattr(module, "fill_region") else None


def run(
    layers: list[Layer],
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    params: Step6Params,
) -> list[Layer]:
    """
    Genere la trajectoire 2D remplie de chaque couche (parois + infill,
    gestion des ponts), dans le repere local (2D) de chaque Layer.

    Args:
        layers: sortie de pipeline/step5_layers.py.
        mesh_data: sortie de pipeline/step1_ingest.py (kdtree, pont
            surface/volume pour la detection de pont).
        overhang_map: sortie de pipeline/step2_overhang.py (is_bridge,
            bridge_span_direction).
        params: parametres de parois/infill/ponts.

    Returns:
        Nouvelle liste de Layer (pas de mutation en place), avec
        toolpath_2d rempli - liste de polylignes (K,3) dans le repere local
        de la couche (u, v, offset hors-plan), ordonnees du centre vers
        l'exterieur. La 3e colonne (offset) preserve la distance reelle au
        plan approximatif de step5 (Decision 6), pour que step8 releve les
        points sur la couche courbe reelle plutot que sur le plan moyen :
          - contour original (paroi la plus exterieure, n=0) : offset EXACT,
            mesure sur contour_3d, jamais interpole (limite #31) ;
          - parois offsettees (n>=1) et infill : offset INTERPOLE depuis les
            boucles de contour du meme groupe si
            `Step6Params.out_of_plane_interpolation_power` est renseigne,
            sinon nul (comportement historique).

        Si `Step6Params.split_rings_into_runs` est actif, chaque anneau ferme
        (paroi, contour original) est en plus DECOUPE en troncons monotones,
        chacun emis comme une polyligne separee et parcouru dans le sens
        montant - voir Decision 5 plus haut. Le nombre de polylignes augmente,
        le nombre total de POINTS est inchange.
    """
    isovalue_rank = {
        value: rank
        for rank, value in enumerate(sorted({l.isovalue for l in layers if l.isovalue is not None}))
    }

    native_layers = _fill_layers_native(layers, mesh_data, overhang_map, params, isovalue_rank)
    if native_layers is not None:
        return native_layers

    finished = [0]

    def one_layer(layer: Layer) -> Layer:
        out = _fill_layer(layer, mesh_data, overhang_map, params, isovalue_rank)
        finished[0] += 1                       # R-017 : couches remplies (le verrou est dans `_progress`)
        _progress.report(finished[0], len(layers))
        return out

    # Couches INDEPENDANTES (fonction pure par couche) : fils d'execution,
    # ordre de sortie conserve (`map`), resultat identique par construction.
    # shapely et le noyau natif relachent le GIL sur leurs calculs.
    workers = _STEP6_THREADS if _STEP6_THREADS is not None else min(8, os.cpu_count() or 1)
    if workers > 1 and len(layers) > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=workers) as pool:
            return list(pool.map(one_layer, layers))
    return [one_layer(layer) for layer in layers]


# Remplissage de couches ENTIERES en natif (cpp/src/core/step6_layer.cpp) : classification,
# pont, parois/remplissage Clipper2, offset hors-plan, decoupage en troncons, en parallele.
# `Step6Params.native_layer=False` rend le chemin couche par couche de Python.
_NATIVE_CHUNKS = 20      # tranches de couches du chemin natif (une fraction de progression par tranche)
_NATIVE_MIN_CHUNK = 32   # jamais moins : un appel natif trop court ne parallelise plus


def _fill_layers_native(layers, mesh_data, overhang_map, params, isovalue_rank):
    """Liste de Layer avec `toolpath_2d` rempli, ou None si le chemin natif ne s'applique pas
    (noyau inactif, motif non natif, kd-tree sans version native)."""
    module = _native_fill_module(params)
    if (not params.native_layer or module is None or not hasattr(module, "fill_layers")
            or params.infill_pattern not in _NATIVE_PATTERN_CODES or not layers):
        return None
    kdtree = getattr(mesh_data.kdtree, "_native", None)
    if kdtree is None:
        return None
    frames = np.array([np.concatenate([np.asarray(f, dtype=float) for f in layer.local_frame])
                       for layer in layers], dtype=float)
    loop_begin = np.zeros(len(layers) + 1, dtype=np.int64)
    sizes: list[int] = []
    chunks: list[np.ndarray] = []
    for k, layer in enumerate(layers):
        loop_begin[k + 1] = loop_begin[k] + len(layer.contour_3d)
        for c in layer.contour_3d:
            arr = np.asarray(c, dtype=float).reshape(-1, 3)
            sizes.append(len(arr))
            chunks.append(arr)
    point_begin = np.zeros(len(sizes) + 1, dtype=np.int64)
    np.cumsum(sizes, out=point_begin[1:])
    points = np.ascontiguousarray(np.concatenate(chunks, axis=0)) if chunks else np.zeros((0, 3))
    sweep = np.array([
        _infill_angle_for_layer(layer, isovalue_rank) if params.infill_pattern is InfillPattern.CROSSED_ZIGZAG else 0.0
        for layer in layers], dtype=float)
    config = {
        "nozzle_width": params.nozzle_width, "wall_count": params.wall_count,
        "degenerate_area_ratio": params.degenerate_area_ratio, "infill_density": params.infill_density,
        "pattern": _NATIVE_PATTERN_CODES[params.infill_pattern],
        "bridge_contour_ratio_threshold": params.bridge_contour_ratio_threshold,
        "bridge_anchor_extra_walls": params.bridge_anchor_extra_walls,
        "power": params.out_of_plane_interpolation_power, "continuous_rings": params.continuous_rings,
        "split_rings_into_runs": params.split_rings_into_runs,
        "split_open_polylines": params.split_open_polylines,
        "run_split_min_amplitude": params.run_split_min_amplitude,
    }
    threads = _STEP6_THREADS if _STEP6_THREADS is not None else 0
    is_bridge = np.ascontiguousarray(overhang_map.is_bridge, dtype=np.uint8)
    span = np.ascontiguousarray(overhang_map.bridge_span_direction, dtype=float)
    # Les couches sont independantes : le noyau les remplit par TRANCHES, pour rendre la main (et la
    # progression, R-017) entre deux tranches. Assez de couches par tranche pour garder le parallelisme.
    n = len(layers)
    step = max(_NATIVE_MIN_CHUNK, -(-n // _NATIVE_CHUNKS))
    result: list = []
    for lo in range(0, n, step):
        hi = min(n, lo + step)
        p_lo, p_hi = int(loop_begin[lo]), int(loop_begin[hi])
        result += list(module.fill_layers(
            np.ascontiguousarray(frames[lo:hi]), np.ascontiguousarray(loop_begin[lo:hi + 1] - loop_begin[lo]),
            np.ascontiguousarray(point_begin[p_lo:p_hi + 1] - point_begin[p_lo]),
            np.ascontiguousarray(points[point_begin[p_lo]:point_begin[p_hi]]),
            np.ascontiguousarray(sweep[lo:hi]), config, kdtree, is_bridge, span, int(threads)))
        _progress.report(hi, n)
    return [replace(layer, toolpath_2d=list(polys)) for layer, polys in zip(layers, result)]


def _fill_layer(layer: Layer, mesh_data: MeshData, overhang_map: OverhangMap, params: Step6Params,
                isovalue_rank: dict) -> Layer:
    if True:
        contours_2d = [_project_to_local_2d(c, layer.local_frame) for c in layer.contour_3d]
        # offsets hors-plan REELS de chaque boucle, calcules une fois par
        # couche : ils servent a la fois au contour n=0 (restitution exacte,
        # limite #31) et de source d'interpolation pour les points generes.
        offsets_3d = [_out_of_plane_offset(c, layer.local_frame) for c in layer.contour_3d]
        open_loops, index_groups = _classify_contour_indices(contours_2d)

        toolpath: list[np.ndarray] = []

        for outer_idx, hole_idxs in index_groups:
            outer = open_loops[outer_idx]
            holes = [open_loops[h] for h in hole_idxs]
            members = {outer_idx, *hole_idxs}
            bridge_direction = None
            outer_out_of_plane = None
            source_2d_parts: list[np.ndarray] = []
            source_offset_parts: list[np.ndarray] = []

            for k, (contour_3d_loop, loop_2d, loop_offsets) in enumerate(zip(
                layer.contour_3d, contours_2d, offsets_3d
            )):
                # les TROUS du groupe portent eux aussi un offset reel connu :
                # les ignorer degraderait l'interpolation a leur voisinage.
                if k not in members:
                    continue
                matches_outer = k == outer_idx

                # contour_3d_loop porte toujours un point de fermeture
                # duplique (premier == dernier, garanti par step5) alors
                # que les boucles issues de _classify_contours sont toujours
                # ouvertes (open_loops) - retire ce doublon pour aligner les
                # longueurs avant d'associer un offset par point.
                open_2d = loop_2d[:-1] if np.allclose(loop_2d[0], loop_2d[-1]) else loop_2d
                open_offsets = loop_offsets[: len(open_2d)]
                source_2d_parts.append(open_2d)
                source_offset_parts.append(open_offsets)

                if matches_outer:
                    outer_out_of_plane = open_offsets
                    bridge_direction = _bridge_direction_for_loop(
                        contour_3d_loop, mesh_data, overhang_map, layer.local_frame,
                        params.bridge_contour_ratio_threshold,
                    )

            source_2d = np.concatenate(source_2d_parts, axis=0) if source_2d_parts else None
            source_offsets = (
                np.concatenate(source_offset_parts, axis=0) if source_offset_parts else None
            )

            if bridge_direction is not None:
                generated = _reorient_bridge_toolpath(outer, holes, bridge_direction, params)
            elif _native_fill_module(params) is not None and params.infill_pattern in _NATIVE_PATTERN_CODES:
                sweep_angle = (
                    _infill_angle_for_layer(layer, isovalue_rank)
                    if params.infill_pattern is InfillPattern.CROSSED_ZIGZAG
                    else 0.0
                )
                valid, infill, offset_walls = _native_fill_module(params).fill_region(
                    outer, holes, params.nozzle_width, params.wall_count, params.degenerate_area_ratio,
                    params.infill_density, _NATIVE_PATTERN_CODES[params.infill_pattern], sweep_angle,
                )
                if not valid:
                    continue
                generated = list(infill) + list(offset_walls) + [outer] + list(holes)
            else:
                walls, innermost = _offset_walls(
                    outer, holes, params.nozzle_width, params.wall_count, params.degenerate_area_ratio
                )
                if not walls:
                    continue
                sweep_angle = (
                    _infill_angle_for_layer(layer, isovalue_rank)
                    if params.infill_pattern is InfillPattern.CROSSED_ZIGZAG
                    else 0.0
                )
                if params.infill_pattern is InfillPattern.OSCILLATING_SPIRAL:
                    axis_local = (
                        _plate_axis_in_local_2d(layer.local_frame, params.spiral_axis_origin)
                        if params.spiral_axis_origin is not None else None
                    )
                    infill = _generate_spiral_infill(
                        innermost, params, _layer_rank(layer, isovalue_rank), axis_local
                    )
                else:
                    infill = _generate_infill(
                        innermost, params.infill_density, params.infill_pattern,
                        params.nozzle_width, sweep_angle,
                    )
                generated = infill + walls

            # contours d'origine n=0 (exterieur et trous), emis OUVERTS
            original_ids = {id(outer), *(id(h) for h in holes)}
            for poly in generated:
                arr = _to_toolpath_3col(
                    poly,
                    _out_of_plane_for_polyline(
                        poly, outer, outer_out_of_plane, source_2d, source_offsets,
                        params.out_of_plane_interpolation_power,
                    ),
                )
                if params.continuous_rings and id(poly) in original_ids and len(arr) >= 3                         and not np.array_equal(arr[0, :2], arr[-1, :2]):
                    arr = np.vstack([arr, arr[:1]])
                toolpath.append(arr)

        if params.split_rings_into_runs:
            toolpath = _split_rings_in_toolpath(
                toolpath, include_open_polylines=params.split_open_polylines,
                share_vertices=params.continuous_rings,
                min_amplitude=params.run_split_min_amplitude,
            )

        return replace(layer, toolpath_2d=toolpath)
