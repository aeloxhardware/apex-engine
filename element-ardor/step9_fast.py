"""Etape 9 : contrôle rapide de collision (arbre k-d de matière passée, calcul par lots).

Découpé de `step9_collision_sim.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step9_collision_sim`.
"""

from __future__ import annotations

from typing import Optional
import numpy as np
from core import _cpp
from core import progress as _progress
from core.mesh_types import CollisionEvent, MachineProfile, ToolpathMachine, ToolpathPoint
from pipeline.step9_params import (  # noqa: F401
    CollisionSimParams,
    printed_margin_of,
)
from pipeline.step9_primitives import (  # noqa: F401
    PLATFORM_CENTER_OFFSET,
    PLATFORM_HALF_HEIGHT,
    PLATFORM_RADIUS,
    _check_platform_disc,
    _positioned_toolhead_vertices,
    _positioned_toolhead_vertices_batch,
)


def _densified_printed_material(points: list[ToolpathPoint], step: float) -> tuple[np.ndarray, np.ndarray]:
    """Points densifies le long de chaque cordon (meme regle que
    _grow_printed_index_to : segment i-1 -> i quand les deux points sont en
    extrusion) et indice `i` du point qui termine ce cordon."""
    pos = np.array([p.position for p in points], dtype=float)
    ext = np.array([p.extrusion for p in points], dtype=float)
    if len(pos) < 2:
        return np.zeros((0, 3)), np.zeros(0, dtype=int)
    ends = np.flatnonzero((ext[1:] > 0.0) & (ext[:-1] > 0.0)) + 1
    if len(ends) == 0:
        return np.zeros((0, 3)), np.zeros(0, dtype=int)
    a, b = pos[ends - 1], pos[ends]
    counts = np.maximum(1, np.ceil(np.linalg.norm(b - a, axis=1) / step).astype(int))
    seg = np.repeat(np.arange(len(ends)), counts + 1)
    # meme t que np.arange(c + 1) / c par segment, sans boucle Python
    offset = np.arange(len(seg)) - np.repeat(np.cumsum(counts + 1) - (counts + 1), counts + 1)
    t = offset / np.repeat(counts, counts + 1)
    dense = a[seg] + t[:, None] * (b[seg] - a[seg])
    return dense, ends[seg]


class FastCollisionChecker:
    """Controle rapide (voir CollisionSimParams.fast_kdtree), reutilisable
    point par point -- par `run` et par la reorientation de l'etape 10, qui
    re-teste d'autres angles au MEME point (la matiere ne depend que des
    positions, jamais des angles).

    Par point teste : les K plus proches points de matiere de chaque sommet de
    tete dans le rayon `reach`, filtres a la matiere deja deposee. Si les K
    sont tous dans le rayon (liste saturee) sans matiere passee, la requete
    exacte par boule est refaite pour ce sommet : la tete survole souvent de
    la matiere FUTURE, qui ne doit masquer aucune matiere passee.
    scipy directement (k > 1, que l'enveloppe core.spatial n'expose pas)."""

    K_NEAR = 16

    def __init__(self, toolpath: ToolpathMachine, machine_profile: MachineProfile, params: CollisionSimParams,
                 toolhead_base_samples: np.ndarray):
        from scipy.spatial import cKDTree

        self.machine, self.params, self.samples = machine_profile, params, toolhead_base_samples
        self.reach = printed_margin_of(params) + machine_profile.nozzle_diameter / 2.0
        self.dense, self.owner = _densified_printed_material(toolpath.points, params.fast_densify_step)
        self._tree = None       # arbre scipy construit a la demande : inutile en mode natif (2026-10-06)
        self.owner_pad = np.r_[self.owner, np.iinfo(np.int64).max]   # indice "absent" de cKDTree
        # Pas 19 du portage C++ : meme minimum exact, en natif, avec elagage
        # des sous-arbres entierement futurs (CLAUDE.md limite #56).
        self.native = None
        if len(self.dense) and _cpp.enabled() and hasattr(_cpp.module(), "PastMaterialIndex"):
            self.native = _cpp.module().PastMaterialIndex(
                np.ascontiguousarray(self.dense, dtype=float), np.ascontiguousarray(self.owner, dtype=np.int64))

    @property
    def tree(self):
        if self._tree is None and len(self.dense):
            from scipy.spatial import cKDTree
            self._tree = cKDTree(self.dense)
        return self._tree

    def penetrations(self, point_idx: int, position: np.ndarray, axis_angles: list[float], verts=None):
        """(penetration plateau ou None, penetration piece ou None). `verts` :
        tete deja positionnee (calcul par lots de `_run_fast_kdtree`)."""
        if verts is None:
            verts = _positioned_toolhead_vertices(self.samples, position, axis_angles, self.machine.rotation_axes)
        part = None
        if self.native is not None:
            best = self.native.min_past_distance(np.ascontiguousarray(verts, dtype=float), int(point_idx),
                                                 float(self.reach))
            if best < self.reach:
                part = self.reach - best
        elif self.tree is not None:
            k = min(self.K_NEAR, len(self.dense))
            dist, idx = self.tree.query(verts, k=k, distance_upper_bound=self.reach)
            dist, idx = dist.reshape(len(verts), k), idx.reshape(len(verts), k)
            past = self.owner_pad[idx] < point_idx
            d_past = np.where(past, dist, np.inf).min(axis=1)
            saturated = np.isfinite(dist[:, -1]) & ~np.isfinite(d_past)
            for vi in np.flatnonzero(saturated):
                near = np.asarray(self.tree.query_ball_point(verts[vi], self.reach), dtype=int)
                near = near[self.owner[near] < point_idx]
                if len(near):
                    d_past[vi] = float(np.linalg.norm(self.dense[near] - verts[vi], axis=1).min())
            best = float(d_past.min())
            if best < self.reach:
                part = self.reach - best
        plate = None
        if self.machine.platform_origin is not None and self.machine.rotation_axes:
            plate = _check_platform_disc(verts, self.machine.platform_origin, self.params.safety_margin,
                                         getattr(self.machine, "platform_radius", None))
        return plate, part


def _platform_disc_batch(verts_block: np.ndarray, machine_profile: MachineProfile,
                         params: CollisionSimParams) -> list:
    """`_check_platform_disc` pour un bloc de poses (P,S,3) : memes operations
    element par element, minimum par pose - memes valeurs que l'appel unitaire."""
    if machine_profile.platform_origin is None or not machine_profile.rotation_axes:
        return [None] * len(verts_block)
    margin = params.safety_margin
    origin = np.asarray(machine_profile.platform_origin, dtype=float)
    # Raccourci EXACT : distance >= ecart axial ; une pose dont tous les points
    # sont a plus de `margin` au-dessus (ou au-dessous) de l'epaisseur du
    # plateau ne peut pas le toucher. Seules les autres paient le calcul
    # complet (2 s sur 2 sur death sans ce tri : 86 M d'hypot).
    dz_all = np.maximum(np.abs(verts_block[..., 2] - origin[2] - PLATFORM_CENTER_OFFSET) - PLATFORM_HALF_HEIGHT, 0.0)
    near = np.flatnonzero(dz_all.min(axis=1) < margin)
    result: list = [None] * len(verts_block)
    if len(near) == 0:
        return result
    local = verts_block[near] - origin
    radial = np.hypot(local[..., 0], local[..., 1])
    radius_value = getattr(machine_profile, "platform_radius", None)
    radius = PLATFORM_RADIUS if radius_value is None else float(radius_value)
    dr = np.maximum(radial - radius, 0.0)
    dz = np.maximum(np.abs(local[..., 2] - PLATFORM_CENTER_OFFSET) - PLATFORM_HALF_HEIGHT, 0.0)
    minimum = np.hypot(dr, dz).min(axis=1)
    for k, m in zip(near, minimum):
        result[k] = None if float(m) >= margin else margin - float(m)
    return result


_FAST_CHECKER_CACHE: dict = {}


def fast_checker_for(toolpath: ToolpathMachine, machine_profile: MachineProfile, params: CollisionSimParams,
                     toolhead_base_samples: np.ndarray) -> "FastCollisionChecker":
    """FastCollisionChecker REUTILISE tant que la matiere est la meme (2026-10-06). La matiere ne
    depend que des positions et des extrusions, jamais des angles : step9 puis les iterations de
    step10 (qui ne changent que des angles) partagent le meme controleur au lieu de le
    reconstruire (5 constructions, 7,3 s sur death). Cache a une entree, compare par valeur
    (positions, extrusions, echantillons de tete, parametres) et par identite du profil machine."""
    pos = np.array([p.position for p in toolpath.points], dtype=float)
    ext = np.array([p.extrusion for p in toolpath.points], dtype=float)
    samples = np.asarray(toolhead_base_samples, dtype=float)
    hit = _FAST_CHECKER_CACHE.get("entry")
    if (hit is not None and hit[0] is machine_profile and hit[1] == params
            and hit[2].shape == pos.shape and np.array_equal(hit[2], pos) and np.array_equal(hit[3], ext)
            and hit[4].shape == samples.shape and np.array_equal(hit[4], samples)):
        return hit[5]
    checker = FastCollisionChecker(toolpath, machine_profile, params, samples)
    _FAST_CHECKER_CACHE["entry"] = (machine_profile, params, pos, ext, samples, checker)
    return checker


def fast_batch_events(checker: "FastCollisionChecker", toolpath: ToolpathMachine,
                      indices: list[int]) -> Optional[list[CollisionEvent]]:
    """Evenements du controle rapide pour une liste de points, en UN appel natif parallele
    (pose de la tete, plateau, matiere passee). None si le noyau natif ne s'applique pas
    (l'appelant retombe sur la boucle Python). Meme ordre que la boucle : par point,
    piece puis plateau."""
    if checker.native is None or not _cpp.enabled() or not hasattr(_cpp.module(), "fast_penetrations_batch"):
        return None
    axes = checker.machine.rotation_axes or []
    if not axes or not indices:
        return None
    pts = toolpath.points
    positions = np.array([pts[i].position for i in indices], dtype=float)
    angles = np.array([pts[i].axis_angles for i in indices], dtype=float).reshape(len(indices), len(axes))
    units = np.array([np.asarray(a.axis, dtype=float) / np.linalg.norm(a.axis) for a in axes])
    origin = checker.machine.platform_origin
    radius = getattr(checker.machine, "platform_radius", None)
    plate, part = _cpp.module().fast_penetrations_batch(
        checker.native, np.ascontiguousarray(checker.samples, dtype=float), np.ascontiguousarray(units),
        np.ascontiguousarray(positions), np.ascontiguousarray(angles), np.asarray(indices, dtype=np.int64),
        float(checker.reach), None if origin is None else np.asarray(origin, dtype=float),
        float(PLATFORM_RADIUS if radius is None else radius), float(PLATFORM_HALF_HEIGHT),
        float(PLATFORM_CENTER_OFFSET), float(checker.params.safety_margin))
    events: list[CollisionEvent] = []
    for k, point_idx in enumerate(indices):
        if not np.isnan(part[k]):
            events.append(CollisionEvent(point_index=point_idx, obstacle_type="printed_part",
                                         penetration_depth=float(part[k])))
        if not np.isnan(plate[k]):
            events.append(CollisionEvent(point_index=point_idx, obstacle_type="platform",
                                         penetration_depth=float(plate[k])))
    return events


_NATIVE_TRANCHES = 20       # tranches du controle rapide natif (une fraction de progression par tranche)


_NATIVE_MIN_TRANCHE = 5000  # jamais moins : un appel natif trop court ne parallelise plus


def _run_fast_kdtree(
    toolpath: ToolpathMachine,
    machine_profile: MachineProfile,
    params: CollisionSimParams,
    toolhead_base_samples: np.ndarray,
    sample_indices: list[int],
) -> list[CollisionEvent]:
    """Variante rapide de la boucle de `run` (voir FastCollisionChecker)."""
    checker = fast_checker_for(toolpath, machine_profile, params, toolhead_base_samples)
    # Points independants : le noyau les traite par tranches (rend la main, donc la progression, R-017).
    indices = list(sample_indices)
    tranche = max(_NATIVE_MIN_TRANCHE, -(-len(indices) // _NATIVE_TRANCHES))
    native_events = fast_batch_events(checker, toolpath, indices[:tranche])
    if native_events is not None:
        _progress.report(min(tranche, len(indices)), len(indices))
        for start in range(tranche, len(indices), tranche):
            native_events += fast_batch_events(checker, toolpath, indices[start:start + tranche])
            _progress.report(min(start + tranche, len(indices)), len(indices))
        return native_events
    events: list[CollisionEvent] = []
    axes = machine_profile.rotation_axes or []
    chunk = 1000
    for start in range(0, len(sample_indices), chunk):
        block = sample_indices[start:start + chunk]
        pts = [toolpath.points[i] for i in block]
        verts_block = _positioned_toolhead_vertices_batch(
            toolhead_base_samples, np.array([p.position for p in pts], dtype=float),
            np.array([p.axis_angles for p in pts], dtype=float), axes,
        ) if axes else None
        batch = None
        if verts_block is not None and checker.native is not None and hasattr(checker.native, "min_past_distance_batch"):
            # tout le bloc d'un coup, en parallele cote natif : meme calcul pose
            # par pose que `penetrations`, donc memes evenements au bit
            best = checker.native.min_past_distance_batch(
                np.ascontiguousarray(verts_block, dtype=float), np.asarray(block, dtype=np.int64), float(checker.reach))
            plates = _platform_disc_batch(verts_block, machine_profile, params)
            batch = (best, plates)
        _progress.report(start + len(block), len(sample_indices))
        for k, point_idx in enumerate(block):
            point = pts[k]
            verts = verts_block[k] if verts_block is not None else None
            if batch is not None:
                best_k = float(batch[0][k])
                part = checker.reach - best_k if best_k < checker.reach else None
                plate = batch[1][k]
            else:
                plate, part = checker.penetrations(point_idx, point.position, point.axis_angles, verts)
            if part is not None:
                events.append(CollisionEvent(point_index=point_idx, obstacle_type="printed_part",
                                             penetration_depth=part))
            if plate is not None:
                events.append(CollisionEvent(point_index=point_idx, obstacle_type="platform",
                                             penetration_depth=plate))
    return events
