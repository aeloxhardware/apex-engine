"""
Trajectoire machine de la voie Bird (corps) : etape 8 d'Apex appliquee aux
couches tranchees dans l'espace deforme (porte de research/S4_Slicer/full_step89.py,
partie « etape 8 »).

Direction visee en chaque point = normale de la couche Bird DANS LA PIECE
REELLE, c'est-a-dire le gradient de phi (hauteur deformee) sur le tetraedre
d'ORIGINE qui contient le point, normalise. La resolution (IK couplee a la
position, filtre median de B, diagnostics) est exactement celle de la voie FEM
(`step8.resolve_raw_points`). La borne de B selon la hauteur et la
reorientation par paliers de la recherche existent deja dans Apex
(`step10.clamp_head_height`, `reorient_coupled_tiers`) : l'orchestrateur les
applique ensuite, comme pour la voie FEM.

Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer).
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from core.mesh_types import MachineProfile, Step8Diagnostic, ToolpathMachine
from pipeline.bird.deformed_slicing import barycentric, tet_gradients
from pipeline.step8_machine_toolpath import Step8Params, _RawPoints, resolve_raw_points


def locate(points: np.ndarray, original: np.ndarray, cells: np.ndarray, centroid_tree: cKDTree, k: int = 12):
    """Tetraedre d'origine contenant chaque point (sinon le moins exterieur des k candidats)."""
    _d, cand = centroid_tree.query(points, k=k)
    best = cand[:, 0].copy()
    score = np.full(len(points), -np.inf)
    for j in range(k):
        with np.errstate(all="ignore"):
            b = barycentric(points, original[cells[cand[:, j]]])
        s = np.nan_to_num(b.min(axis=1), nan=-np.inf)
        better = s > score
        best[better], score[better] = cand[better, j], s[better]
    return best


def layer_normals(cells: np.ndarray, original: np.ndarray, deformed: np.ndarray) -> np.ndarray:
    """Normale de couche Bird par tetraedre d'origine : gradient normalise de la hauteur deformee."""
    with np.errstate(all="ignore"):
        grad = tet_gradients(original, cells, deformed[:, 2])
    return grad / np.maximum(np.linalg.norm(grad, axis=1, keepdims=True), 1e-12)


def raw_points_from_layers(layers: list[dict], cells: np.ndarray, original: np.ndarray, deformed: np.ndarray,
                           params: Step8Params, origin_shift: np.ndarray = np.zeros(3)) -> _RawPoints:
    """Points bruts dans l'ORDRE D'IMPRESSION (ordre des couches puis des
    polylignes) : premier point de chaque polyligne = trajet a vide.
    `origin_shift` : decalage repere Bird -> repere de `toolpath_3d` (soustrait
    pour localiser les points dans le maillage d'origine)."""
    normals = layer_normals(cells, original, deformed)
    tree = cKDTree(original[cells].mean(axis=1))
    positions, directions, layer_ids, seg_ids, travel = [], [], [], [], []
    for layer_index, layer in enumerate(layers):
        for seg, poly in enumerate(layer["toolpath_3d"]):
            poly = np.asarray(poly, dtype=float)
            if len(poly) == 0:
                continue
            positions.append(poly)
            layer_ids.append(np.full(len(poly), layer_index))
            seg_ids.append(np.full(len(poly), seg))
            t = np.zeros(len(poly), dtype=bool)
            t[0] = True
            travel.append(t)
    if not positions:
        empty = np.zeros((0, 3))
        return _RawPoints(empty, empty, np.zeros(0, np.int64), np.zeros(0, np.int64), [], np.zeros(0), np.zeros(0),
                          np.zeros(0, bool))
    travel = np.concatenate(travel)
    # Localisation de TOUS les points en un appel (le calcul de chaque point ne depend
    # que de lui : requete k-d et resolution 3x3 independantes) -- identique au bit a
    # l'appel par polyligne, 11 s -> ~1 s sur le Benchy (3 685 polylignes).
    all_points = np.vstack(positions)
    directions = [normals[locate(all_points - origin_shift, original, cells, tree)]]
    travel_feed = params.travel_feedrate if params.travel_feedrate > 0 else params.default_feedrate
    return _RawPoints(
        position=np.vstack(positions), desired_direction=np.vstack(directions),
        layer_id=np.concatenate(layer_ids).astype(np.int64), segment_index=np.concatenate(seg_ids).astype(np.int64),
        color_id=[None] * len(travel),
        feedrate=np.where(travel, float(travel_feed), float(params.default_feedrate)),
        extrusion=np.where(travel, 0.0, 1.0), is_travel=travel)


def run_from_passes(passes: list[dict], cells: np.ndarray, original: np.ndarray, deformed: np.ndarray,
                    machine_profile: MachineProfile, params: Step8Params = Step8Params()
                    ) -> tuple[ToolpathMachine, list[Step8Diagnostic]]:
    """Trajectoire corps + lissage dans l'ordre de P5 (print_order.order_and_check,
    points dans le repere Bird). Corps : resolution d'Apex (normale de couche
    Bird). Lissage : angles de P5 tels quels (cible = normale de surface), comme
    full_step89 ; color_id 1 = lissage, 0 = corps (convention de la recherche,
    lue par l'export pour l'epaisseur de cordon)."""
    from core.mesh_types import ToolpathPoint
    body_layers = [{"toolpath_3d": p["points"]} for p in passes if p["kind"] == "body"]
    body_tp, body_diag = run(body_layers, cells, original, deformed, machine_profile, params)
    body_iter = iter(body_tp.points)
    points, index_map = [], {}
    body_counter = 0
    for p in passes:
        if p["kind"] == "body":
            for poly in p["points"]:
                for _ in range(len(poly)):
                    q = next(body_iter)
                    q.color_id = 0
                    index_map[body_counter] = len(points)
                    body_counter += 1
                    points.append(q)
        else:
            for poly, angles in zip(p["points"], p["angles"]):
                for i, (pos, a) in enumerate(zip(np.asarray(poly), np.asarray(angles))):
                    travel = i == 0
                    points.append(ToolpathPoint(
                        position=np.asarray(pos, dtype=float), axis_angles=[float(x) for x in a],
                        extrusion=0.0 if travel else 1.0,
                        feedrate=float(params.travel_feedrate if travel and params.travel_feedrate > 0
                                       else params.default_feedrate),
                        color_id=1))
    diagnostics = [Step8Diagnostic(point_index=index_map[d.point_index], layer_id=d.layer_id,
                                   segment_index=d.segment_index, residual_angle_error_deg=d.residual_angle_error_deg,
                                   position=d.position) for d in body_diag]
    return ToolpathMachine(points=points), diagnostics


def run(layers: list[dict], cells: np.ndarray, original: np.ndarray, deformed: np.ndarray,
        machine_profile: MachineProfile, params: Step8Params = Step8Params(),
        origin_shift: np.ndarray = np.zeros(3)) -> tuple[ToolpathMachine, list[Step8Diagnostic]]:
    raw = raw_points_from_layers(layers, cells, original, deformed, params, origin_shift)
    return resolve_raw_points(raw, machine_profile, params)
