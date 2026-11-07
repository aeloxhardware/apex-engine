"""Etape 8 : assemblage des points bruts depuis l'ordonnancement, normales par point.

Découpé de `step8_machine_toolpath.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step8_machine_toolpath`.
"""

from __future__ import annotations

from dataclasses import field
from typing import Optional
import numpy as np
from core import _cpp
from core.mesh_types import Layer, PrintSchedule
from element_ardor.step6_toolpath2d import _lift_to_3d
from element_ardor.step8_params import (  # noqa: F401
    Step8Params,
    ToolpathResolutionError,
)
from element_ardor.step8_smoothing import (  # noqa: F401
    _RawPoints,
    _point_arrays,
)


def _layer_normal(local_frame: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    """
    Normale du repere local d'une couche = direction d'impression desiree en
    tout point de cette couche (meme lecture que step3/step4/step5 : le
    champ scalaire encode une direction d'impression par construction).
    """
    _origin, axis_u, axis_v = local_frame
    normal = np.cross(axis_u, axis_v)
    norm = np.linalg.norm(normal)
    if norm < 1e-12:
        return np.array([0.0, 0.0, 1.0])
    return normal / norm


def _iter_schedule_points_native(schedule, layers_by_id, params):
    """Assemblage natif des points bruts (cpp/src/core/step8_assembly.cpp), au bit
    identique. None = noyau natif inactif ou entree inhabituelle : le chemin Python
    ci-dessous s'applique (il leve aussi les erreurs de planning invalide)."""
    if not (_cpp.enabled() and hasattr(_cpp.module(), "build_raw_points")):
        return None
    table_index: dict[int, int] = {}
    layers_used: list[Layer] = []
    seg_layer: list[int] = []
    seg_index: list[int] = []
    seg_pass: list[int] = []
    pass_colors: list = []
    for pass_number, pass_ in enumerate(schedule.passes):
        pass_colors.append(pass_.color_id)
        for layer_id, segment_index in pass_.segments:
            slot = table_index.get(layer_id)
            if slot is None:
                layer = layers_by_id.get(layer_id)
                if layer is None or layer.toolpath_2d is None:
                    return None
                slot = table_index[layer_id] = len(layers_used)
                layers_used.append(layer)
            if segment_index >= len(layers_used[slot].toolpath_2d):
                return None
            seg_layer.append(slot)
            seg_index.append(segment_index)
            seg_pass.append(pass_number)
    if not seg_layer:
        return None
    columns = {np.asarray(p).shape[1] for layer in layers_used for p in layer.toolpath_2d
               if np.asarray(p).ndim == 2 and len(p)}
    if len(columns) != 1 or any(np.asarray(p).ndim != 2 for layer in layers_used for p in layer.toolpath_2d):
        return None
    ncol = columns.pop()
    if ncol not in (2, 3):
        return None
    frames = np.array([np.concatenate([np.asarray(f, dtype=float) for f in layer.local_frame])
                       for layer in layers_used], dtype=float)
    reversed_flags = np.array([bool(getattr(layer, "deposition_reversed", False)) for layer in layers_used],
                              dtype=np.uint8)
    ids = np.array([layer.id for layer in layers_used], dtype=np.int64)
    poly_begin = np.zeros(len(layers_used) + 1, dtype=np.int64)
    point_counts: list[int] = []
    chunks: list[np.ndarray] = []
    for slot, layer in enumerate(layers_used):
        poly_begin[slot + 1] = poly_begin[slot] + len(layer.toolpath_2d)
        for p in layer.toolpath_2d:
            arr = np.asarray(p, dtype=float).reshape(-1, ncol) if len(p) else np.zeros((0, ncol))
            point_counts.append(len(arr))
            chunks.append(arr)
    point_begin = np.zeros(len(point_counts) + 1, dtype=np.int64)
    np.cumsum(point_counts, out=point_begin[1:])
    points = np.ascontiguousarray(np.concatenate(chunks, axis=0) if chunks else np.zeros((0, ncol)))
    position, direction, layer_id, segment, pass_idx, travel = _cpp.module().build_raw_points(
        frames, reversed_flags, ids, poly_begin, point_begin, points,
        np.asarray(seg_layer, dtype=np.int64), np.asarray(seg_index, dtype=np.int64),
        np.asarray(seg_pass, dtype=np.int64))
    travel = travel.astype(bool)
    travel_feedrate = params.travel_feedrate if params.travel_feedrate > 0 else params.default_feedrate
    colors = [pass_colors[i] for i in pass_idx.tolist()]
    return _RawPoints(
        position=position, desired_direction=direction, layer_id=layer_id, segment_index=segment,
        color_id=colors, feedrate=np.where(travel, float(travel_feedrate), float(params.default_feedrate)),
        extrusion=np.where(travel, 0.0, 1.0), is_travel=travel,
    )


def _iter_schedule_points(
    schedule: PrintSchedule, layers_by_id: dict[int, Layer], params: Step8Params
) -> _RawPoints:
    """
    Parcourt `schedule.passes` DANS L'ORDRE (jamais re-derive : c'est l'ordre d'execution gele par step7),
    releve chaque segment 2D vers le 3D, et insere un point de trajet a vide (`extrusion=0.0`) a chaque
    frontiere de polyligne - y compris entre deux passes non contigues. Une frontiere est detectee des que la
    cle (layer_id, index) du segment courant differe de celle du segment precedent.

    Rend des TABLEAUX (`_RawPoints`), assembles par blocs (un par polyligne) : memes valeurs que la liste
    d'objets d'origine (2026-10-01).
    """
    native_result = _iter_schedule_points_native(schedule, layers_by_id, params)
    if native_result is not None:
        return native_result

    previous_key: Optional[tuple[int, int]] = None
    normal_cache: dict[int, np.ndarray] = {}
    blocks = _RawBlocks()
    for pass_ in schedule.passes:
        for layer_id, segment_index in pass_.segments:
            layer, polyline_2d = _schedule_segment_polyline(layers_by_id, layer_id, segment_index)
            if len(polyline_2d) == 0:
                continue

            points_3d = _lift_to_3d(polyline_2d, layer.local_frame)
            # une normale PAR COUCHE, pas par polyligne : 29 691 np.cross sur death @ 0,8 (0,84 s) pour 257
            # valeurs distinctes
            desired_direction = normal_cache.get(layer_id)
            if desired_direction is None:
                desired_direction = _layer_normal(layer.local_frame)
                normal_cache[layer_id] = desired_direction
            # Le sens de depot est donne par l'ORDRE de parcours des isosurfaces, pas par le signe de
            # grad(phi) : une couche atteinte en DESCENDANT le champ (Layer.deposition_reversed) se depose
            # sous la matiere existante, la tete pointe a l'oppose de la normale.
            if getattr(layer, "deposition_reversed", False):
                desired_direction = -desired_direction
            key = (layer_id, segment_index)

            if previous_key is not None and key != previous_key:
                blocks.add(points_3d[:1], desired_direction, layer_id, segment_index, pass_.color_id, True)
            blocks.add(points_3d, desired_direction, layer_id, segment_index, pass_.color_id, False)
            previous_key = key
    return blocks.to_raw_points(params)


def gradient_direction_field(tet_mesh, scalar_field, nodal: bool = False):
    """Champ de directions +grad(phi) normalise, lu au tetraedre dont le
    centroide est le plus proche (gradient P1 constant par tetraedre).
    Approximation assumee : le plus proche centroide n'est pas toujours le
    tetraedre contenant, mais la normale varie peu entre voisins, et le filtre
    median de B (`free_axis_median_window`) lisse le reste."""
    from core.spatial import KDTree
    centroids = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    grad = np.asarray(scalar_field.gradient, dtype=float)
    norm = np.linalg.norm(grad, axis=1)
    unit = np.full_like(grad, np.nan)
    ok = norm > 1e-12
    unit[ok] = grad[ok] / norm[ok, None]
    tree = KDTree(centroids)

    if not nodal:
        def field(points):
            _, idx = tree.query(np.asarray(points, dtype=float))
            return unit[np.asarray(idx, dtype=np.int64)]

        return field

    # Variante CONTINUE (`Step8Params.per_point_normal_nodal`) : gradient
    # moyenne aux noeuds (pondere par le volume des tetraedres incidents),
    # puis interpole en barycentrique dans le tetraedre le plus proche
    # (poids negatifs ramenes a 0 hors du tetraedre). Supprime la
    # discontinuite d'une maille a l'autre, cause des sauts de B (limite #63).
    tets = tet_mesh.tetrahedra
    verts = tet_mesh.tet_vertices
    corners = verts[tets]
    vol = np.abs(np.linalg.det(np.transpose(corners[:, 1:] - corners[:, :1], (0, 2, 1)))) / 6.0
    node_g = np.zeros((len(verts), 3))
    safe = np.where(ok[:, None], unit, 0.0) * vol[:, None]
    for k in range(4):
        np.add.at(node_g, tets[:, k], safe)
    nn = np.linalg.norm(node_g, axis=1)
    node_unit = np.full_like(node_g, np.nan)
    node_ok = nn > 1e-12
    node_unit[node_ok] = node_g[node_ok] / nn[node_ok, None]

    def field_nodal(points):
        p = np.asarray(points, dtype=float)
        _, idx = tree.query(p)
        idx = np.asarray(idx, dtype=np.int64)
        c = corners[idx]                                   # (M,4,3)
        m = np.transpose(c[:, 1:] - c[:, :1], (0, 2, 1))   # (M,3,3)
        rhs = (p - c[:, 0])[:, :, None]
        try:
            lam = np.linalg.solve(m, rhs)[:, :, 0]
        except np.linalg.LinAlgError:
            lam = np.full((len(p), 3), 0.25)
        w = np.column_stack([1.0 - lam.sum(axis=1), lam])
        w = np.clip(w, 0.0, None)
        w /= np.maximum(w.sum(axis=1), 1e-12)[:, None]
        g = np.einsum("mk,mkj->mj", w, np.nan_to_num(node_unit[tets[idx]]))
        n = np.linalg.norm(g, axis=1)
        out = np.full_like(g, np.nan)
        good = n > 1e-9
        out[good] = g[good] / n[good, None]
        fallback = ~good
        out[fallback] = unit[idx[fallback]]
        return out

    return field_nodal


def _apply_per_point_normals(raw_points, layers_by_id, direction_field) -> None:
    """Remplace la normale de couche par la normale locale, au signe pres de
    `deposition_reversed` ; garde la normale de couche la ou le champ est
    inexploitable ou oppose a elle (jamais un retournement de tete induit par
    un gradient local aberrant)."""
    positions, current = _point_arrays(raw_points)
    local = direction_field(positions)
    arrays = isinstance(raw_points, _RawPoints)
    if (arrays and _cpp.enabled() and hasattr(_cpp.module(), "apply_per_point_normals")
            and raw_points.desired_direction.dtype == np.float64
            and raw_points.desired_direction.flags["C_CONTIGUOUS"]
            and raw_points.desired_direction.flags["WRITEABLE"]):
        ids = np.asarray(raw_points.layer_id)
        reversed_by_id = {int(i): bool(getattr(layers_by_id[int(i)], "deposition_reversed", False))
                          for i in np.unique(ids)}
        flags = np.fromiter((reversed_by_id[int(i)] for i in ids), dtype=np.uint8, count=len(ids))
        _cpp.module().apply_per_point_normals(
            np.ascontiguousarray(local, dtype=np.float64), np.ascontiguousarray(current, dtype=np.float64),
            flags, raw_points.desired_direction)
        return
    for i, d in enumerate(local):
        if not np.all(np.isfinite(d)):
            continue
        layer = layers_by_id[int(raw_points.layer_id[i]) if arrays else raw_points[i].layer_id]
        if getattr(layer, "deposition_reversed", False):
            d = -d
        if float(np.dot(d, current[i])) <= 0.0:
            continue
        if arrays:
            raw_points.desired_direction[i] = d
        else:
            raw_points[i].desired_direction = d


def _schedule_segment_polyline(
    layers_by_id: dict[int, Layer], layer_id: int, segment_index: int
) -> tuple[Layer, np.ndarray]:
    layer = layers_by_id.get(layer_id)
    if layer is None:
        raise ToolpathResolutionError(
            f"PrintSchedule reference la couche {layer_id}, absente de `layers`."
        )
    toolpath = layer.toolpath_2d
    if toolpath is None:
        raise ToolpathResolutionError(
            f"La couche {layer_id} n'a pas de toolpath_2d - "
            f"element_ardor/step7_scheduling.py doit etre execute avant l'etape 8."
        )
    if segment_index >= len(toolpath):
        raise ToolpathResolutionError(
            f"PrintSchedule reference le segment {segment_index} de la couche "
            f"{layer_id}, hors limites ({len(toolpath)} segments disponibles)."
        )
    return layer, toolpath[segment_index]


class _RawBlocks:
    """Blocs de points bruts (un par polyligne ou point de trajet a vide), assembles en `_RawPoints`."""

    def __init__(self):
        self.positions: list[np.ndarray] = []
        self.direction: list[np.ndarray] = []
        self.count: list[int] = []
        self.layer: list[int] = []
        self.segment: list[int] = []
        self.color: list[Optional[int]] = []
        self.travel: list[bool] = []

    def add(self, points, direction, layer_id: int, segment_index: int, color, travel: bool) -> None:
        self.positions.append(points)
        self.direction.append(direction)
        self.count.append(len(points))
        self.layer.append(layer_id)
        self.segment.append(segment_index)
        self.color.append(color)
        self.travel.append(travel)

    def to_raw_points(self, params: Step8Params) -> _RawPoints:
        travel_feedrate = params.travel_feedrate if params.travel_feedrate > 0 else params.default_feedrate
        if not self.positions:
            empty3 = np.zeros((0, 3))
            return _RawPoints(empty3, empty3.copy(), np.zeros(0, np.int64), np.zeros(0, np.int64), [],
                              np.zeros(0), np.zeros(0), np.zeros(0, bool))
        counts = np.asarray(self.count, dtype=np.int64)
        travel = np.repeat(np.asarray(self.travel, dtype=bool), counts)
        color_ids: list[Optional[int]] = []
        for color, count in zip(self.color, self.count):
            color_ids.extend([color] * count)
        return _RawPoints(
            position=np.concatenate(self.positions, axis=0).astype(float, copy=False),
            desired_direction=np.repeat(np.asarray(self.direction, dtype=float), counts, axis=0),
            layer_id=np.repeat(np.asarray(self.layer, dtype=np.int64), counts),
            segment_index=np.repeat(np.asarray(self.segment, dtype=np.int64), counts),
            color_id=color_ids,
            feedrate=np.where(travel, float(travel_feedrate), float(params.default_feedrate)),
            extrusion=np.where(travel, 0.0, 1.0),
            is_travel=travel,
        )
