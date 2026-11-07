"""Etape 10 : réorientation locale (paliers couplés, borne de hauteur, candidats, re-vérification).

Découpé de `step10_adjustment.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step10_adjustment`.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional
import numpy as np
import trimesh
from core import _cpp
from core.geometry_utils import _rodrigues_rotate_single, composed_rotation_angles_and_directions
from core.mesh_types import CollisionEvent, MachineProfile, RotationAxisSpec, ToolpathMachine, ToolpathPoint
from element_ardor.step9_collision_sim import CollisionSimParams, _evaluate_point, _grow_printed_index_to, printed_margin_of, _filter_nozzle_contact_zone, _PrintedGeometryIndex, _surface_sample_points, platform_disc
from element_ardor.step10_params import (  # noqa: F401
    Step10Params,
)


def clamp_head_height(
    toolpath: ToolpathMachine,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    axis_index: Optional[int] = None,
    grid_step_deg: float = 0.5,
) -> tuple[ToolpathMachine, int, list[float]]:
    """
    Borne l'angle d'orientation selon la hauteur de la pointe (porte de la
    recherche le 2026-09-30, `full_step89.py --height-clamp`) : pour chaque
    valeur de l'axe `axis_index`, point le plus bas de la VRAIE tete sous la
    pointe ; un point dont la tete passerait sous la marge du plateau est
    ramene a l'angle autorise le plus proche. Mesure en recherche : des
    contacts plateau jusqu'a 78 mm de penetration (tete retournee pres du
    plateau) etaient de vraies traversees. Redondant avec la reorientation
    de `run`, voulu par l'utilisateur (« la redondance est toujours bonne »).

    Table pre-calculee une fois : valide seulement si les AUTRES axes
    tournent autour de la verticale (sans effet sur z) -- sinon ValueError,
    jamais un resultat faux en silence. `axis_index` par defaut : le dernier
    axe non couple a la position (B sur la machine actuelle).

    Rend (trajectoire, nombre de points ramenes, ecarts angulaires).
    L'entree n'est jamais mutee.
    """
    axes = machine_profile.rotation_axes or []
    model = machine_profile.nozzle_model
    if not axes or not isinstance(model, trimesh.Trimesh):
        return ToolpathMachine(points=list(toolpath.points)), 0, []
    coupled = set(machine_profile.position_coupled_axes or [])
    if axis_index is None:
        free = [i for i in range(len(axes)) if i not in coupled]
        if not free:
            return ToolpathMachine(points=list(toolpath.points)), 0, []
        axis_index = free[-1]
    for i, spec in enumerate(axes):
        a = np.asarray(spec.axis, dtype=float)
        if i != axis_index and np.linalg.norm(np.cross(a / np.linalg.norm(a), [0.0, 0.0, 1.0])) > 1e-9:
            raise ValueError("clamp_head_height : les autres axes doivent tourner autour de la verticale")
    grid, zlow = _head_zlow_table(model, axes, axis_index, sim_params, grid_step_deg)
    plate_z = float(machine_profile.platform_origin[2]) if machine_profile.platform_origin is not None else 0.0

    if _cpp.enabled() and hasattr(_cpp.module(), "clamp_head_height_apply") and len(toolpath.points):
        return _clamp_apply_native(toolpath, axis_index, grid, zlow, plate_z, sim_params)

    points, moved, deviations = [], 0, []
    for p in toolpath.points:
        ok = p.position[2] - plate_z + zlow >= sim_params.safety_margin
        current = p.axis_angles[axis_index]
        if ok.any() and not ok[np.abs(grid - current).argmin()]:
            best = float(grid[ok][np.abs(grid[ok] - current).argmin()])
            angles = list(p.axis_angles)
            angles[axis_index] = best
            deviations.append(abs(best - current))
            p = replace(p, axis_angles=angles)
            moved += 1
        points.append(p)
    return ToolpathMachine(points=points), moved, deviations


def reorient_coupled_tiers(
    toolpath: ToolpathMachine,
    contact_indices: list[int],
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    tiers_deg: tuple[float, ...] = (45.0, 90.0, 181.0),
    grid_step_deg: float = 2.5,
    full_mm: float = 1.0,
    taper_mm: float = 3.0,
) -> tuple[ToolpathMachine, list[tuple[int, Optional[float], Optional[int]]]]:
    """
    Reorientation des points en contact pour une machine a axe COUPLE a la
    position (limite #59 : C impose, seul l'axe libre oriente) -- portee de la
    recherche le 2026-09-30 (`full_step89.py --reorient/--apply` : Benchy
    83 -> 9 contacts, dino 217 -> 0).

    Pour chaque point en contact : valeurs de l'axe libre de la plus proche a
    la plus eloignee, par PALIERS de tolerance sur l'ecart entre la direction
    de la tete obtenue et celle d'origine (un ecart de depot vaut mieux qu'une
    collision ; le dernier palier est sans limite). Controle : le controle
    rapide de step9 (`FastCollisionChecker`), plateau et piece. La correction
    trouvee est appliquee aux voisins de la meme polyligne : poids 1 a
    +-`full_mm`, decroissant lineairement jusqu'a 0 a +-`taper_mm` (sans cela,
    un seul point bascule et la tete saute).

    La reference de l'ecart est la direction de tete ACTUELLE du point (celle
    que step8 a resolue), faute de normale visee dans ToolpathPoint.
    Rend (trajectoire, journal [(indice, angle trouve ou None, palier)]).
    L'entree n'est jamais mutee.
    """
    from element_ardor.step9_collision_sim import fast_checker_for, _filter_nozzle_contact_zone

    axes = machine_profile.rotation_axes or []
    points = list(toolpath.points)
    if not contact_indices or not axes or not isinstance(machine_profile.nozzle_model, trimesh.Trimesh):
        return ToolpathMachine(points=points), []
    coupled = set(machine_profile.position_coupled_axes or [])
    free = [i for i in range(len(axes)) if i not in coupled]
    if len(free) != 1:
        raise ValueError("reorient_coupled_tiers : un seul axe d'orientation libre attendu")
    ax = free[0]
    samples = _filter_nozzle_contact_zone(
        _surface_sample_points(machine_profile.nozzle_model, sim_params.toolhead_sample_spacing),
        sim_params.nozzle_contact_clearance,
    )
    checker = fast_checker_for(toolpath, machine_profile, sim_params, samples)
    lo, hi = axes[ax].angle_range_deg
    grid = np.arange(lo, hi + 1e-9, grid_step_deg)
    setup = _TierSetup(axes, ax, grid, samples, checker, lo, hi, tiers_deg, full_mm, taper_mm, sim_params,
                       machine_profile)

    pos = np.array([p.position for p in points], dtype=float)
    ext = np.array([p.extrusion for p in points], dtype=float)
    new_angles = np.array([p.axis_angles for p in points], dtype=float)
    if checker.native is not None and hasattr(_cpp.module(), "reorient_coupled_tiers"):
        new_angles, log = _reorient_tiers_native(setup, pos, ext, new_angles, contact_indices)
    else:
        new_angles, log = _reorient_tiers_python(setup, pos, ext, new_angles, contact_indices)
    out = [replace(p, axis_angles=[float(v) for v in new_angles[i]])
           if not np.array_equal(new_angles[i], p.axis_angles) else p for i, p in enumerate(points)]
    return ToolpathMachine(points=out), log


def _direction_for_angles(
    axes: list[RotationAxisSpec], angles_deg: list[float], reference_direction: np.ndarray
) -> np.ndarray:
    """
    Direction de depot resultant de `angles_deg`, avec la MEME convention de
    composition que composed_rotation_angles_and_directions / step8 /
    step9 (rotations appliquees en ordre INVERSE de la liste) - jamais une
    convention propre a l'etape 10.
    """
    direction = np.asarray(reference_direction, dtype=float)
    for spec, angle_deg in zip(reversed(axes), reversed(list(angles_deg))):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        direction = _rodrigues_rotate_single(direction, axis, np.radians(angle_deg))
    return direction


def _collisions_at_points(
    toolpath: ToolpathMachine,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    point_indices: list[int],
) -> list[CollisionEvent]:
    """
    Re-teste la collision UNIQUEMENT aux indices demandes, avec les primitives de l'etape 9. Exact : l'index
    de matiere imprimee ne depend que de `position` (jamais de `axis_angles`), donc une reorientation locale
    laisse la geometrie imprimee identique et seul le point reoriente doit etre re-teste. L'index est
    reconstruit jusqu'au plus grand indice demande (croissance monotone, cordon en cours de depot exclu).
    """
    nozzle_model = machine_profile.nozzle_model
    if not isinstance(nozzle_model, trimesh.Trimesh) or len(nozzle_model.vertices) == 0:
        return []
    if not point_indices:
        return []

    # toolhead_sample_spacing (echelle de la TETE), jamais nozzle_diameter (limite #26 : 17734 points de test
    # par appel au lieu de ~1540 sur aetheris-4.glb).
    if sim_params.fast_kdtree:
        checker = _fast_checker(toolpath, machine_profile, sim_params)
        indices = [i for i in sorted(set(point_indices)) if i < len(toolpath.points)]
        from element_ardor.step9_collision_sim import fast_batch_events
        batch = fast_batch_events(checker, toolpath, indices)
        if batch is not None:
            return batch
        return [e for i in indices for e in _fast_events(checker, i, toolpath.points[i])]

    # Zone de contact de la buse retiree, comme step9 (sinon tout point re-teste est en collision, 2026-10-01).
    samples = _filter_nozzle_contact_zone(
        _surface_sample_points(nozzle_model, sim_params.toolhead_sample_spacing),
        sim_params.nozzle_contact_clearance,
    )
    printed_index = _PrintedGeometryIndex(
        machine_profile.nozzle_diameter, sim_params.grid_cell_size_multiplier
    )

    events: list[CollisionEvent] = []
    last_grown_idx = 0
    for point_idx in sorted(set(point_indices)):
        if point_idx >= len(toolpath.points):
            continue
        last_grown_idx = _grow_printed_index_to(
            printed_index, toolpath.points, point_idx, last_grown_idx
        )
        events.extend(
            _evaluate_point(
                point_idx,
                samples,
                toolpath.points[point_idx],
                machine_profile,
                printed_index,
                sim_params.safety_margin,
                printed_margin_of(sim_params),
            )
        )
    return events


def _violates_rate_limits(
    points: list[ToolpathPoint],
    point_idx: int,
    candidate_angles: list[float],
    machine_profile: MachineProfile,
    neighbor_margin: int,
) -> bool:
    """
    Une reorientation ponctuelle peut creer un saut d'angle que le lissage
    temporel de l'etape 8 avait precisement elimine. On rejette donc tout
    candidat qui violerait `max_angular_velocity` vis-a-vis de ses voisins
    immediats - sans cette verification, l'etape 10 "resoudrait" une
    collision en produisant une trajectoire mecaniquement inexecutable.

    Limite non traitee (assumee) : seule la vitesse est verifiee, pas
    l'acceleration (qui demanderait de raisonner sur trois points
    consecutifs et de re-lisser toute la serie).
    """
    max_velocity = machine_profile.max_angular_velocity
    if max_velocity <= 0:
        return False  # limite non specifiee : aucune contrainte a faire respecter

    max_velocity_deg_s = np.degrees(max_velocity)

    for offset in range(1, max(neighbor_margin, 1) + 1):
        for neighbor_idx in (point_idx - offset, point_idx + offset):
            if neighbor_idx < 0 or neighbor_idx >= len(points):
                continue
            distance = float(
                np.linalg.norm(points[neighbor_idx].position - points[point_idx].position)
            )
            feedrate = max(points[neighbor_idx].feedrate, 1e-6)
            dt = max(distance / feedrate, 1e-9)
            for current, neighbor in zip(candidate_angles, points[neighbor_idx].axis_angles):
                if abs(current - neighbor) / dt > max_velocity_deg_s:
                    return True
    return False


def _with_angles(toolpath: ToolpathMachine, point_idx: int, angles: list[float]) -> ToolpathMachine:
    """Copie de la trajectoire avec les angles d'un seul point remplaces (jamais de mutation en place)."""
    points = list(toolpath.points)
    points[point_idx] = replace(points[point_idx], axis_angles=list(angles))
    return ToolpathMachine(points=points)


@dataclass
class _ReorientationContext:
    """
    Etat partage entre TOUS les points traites par le niveau 1 d'une meme iteration externe de `run()`.
    Les echantillons de tete ne dependent que du modele ; l'index de matiere imprimee ne depend que de
    `position` et grandit de facon monotone (points traites en ordre CROISSANT, meme invariant que l'etape 9).
    Sans cet etat partage, l'index etait reconstruit a chaque point (2.5 appels/s -> 0.39 appels/s au 200e).
    """

    toolhead_samples: Optional[np.ndarray]
    printed_index: Optional[_PrintedGeometryIndex]
    last_grown_idx: int = 0
    # Mode rapide de step9 (`CollisionSimParams.fast_kdtree`) : la
    # re-verification utilise le MEME controleur que step9, pas l'index exact
    # par capsules (incoherence corrigee le 2026-10-01 : step10 restait sur
    # l'exact, des heures sur death.stl et ses 24 789 collisions).
    fast_checker: Optional[object] = None
    # Tableaux (positions, avances, angles COURANTS) du chemin natif de la
    # reorientation locale ; les angles suivent les reorientations reussies.
    positions: Optional[np.ndarray] = None
    feedrates: Optional[np.ndarray] = None
    angles: Optional[np.ndarray] = None

    def grow_to(self, points: list[ToolpathPoint], point_idx: int) -> None:
        if self.printed_index is None:
            return
        self.last_grown_idx = _grow_printed_index_to(
            self.printed_index, points, point_idx, self.last_grown_idx
        )


def _build_reorientation_context(
    toolpath: ToolpathMachine, machine_profile: MachineProfile, sim_params: CollisionSimParams
) -> _ReorientationContext:
    """
    Construit l'etat partage une seule fois par iteration externe de `run()`.
    `last_grown_idx=0` : l'appelant doit parcourir les points en ordre
    CROISSANT et appeler `grow_to` avant chaque test, jamais en desordre -
    meme invariant "jamais en avance sur le point teste" que l'etape 9.
    """
    nozzle_model = machine_profile.nozzle_model
    has_toolhead = isinstance(nozzle_model, trimesh.Trimesh) and len(nozzle_model.vertices) > 0
    if not has_toolhead:
        return _ReorientationContext(toolhead_samples=None, printed_index=None)

    # toolhead_sample_spacing (echelle de la TETE), jamais nozzle_diameter
    # (echelle du CORDON) - meme distinction que step9 (CLAUDE.md limite #26).
    # meme filtre de la zone de contact que step9 (voir _collisions_at_points)
    toolhead_samples = _filter_nozzle_contact_zone(
        _surface_sample_points(nozzle_model, sim_params.toolhead_sample_spacing),
        sim_params.nozzle_contact_clearance,
    )
    printed_index = _PrintedGeometryIndex(
        machine_profile.nozzle_diameter, sim_params.grid_cell_size_multiplier
    )
    fast_checker = _fast_checker(toolpath, machine_profile, sim_params) if sim_params.fast_kdtree else None
    context = _ReorientationContext(toolhead_samples=toolhead_samples, printed_index=printed_index,
                                    fast_checker=fast_checker)
    if (fast_checker is not None and getattr(fast_checker, "native", None) is not None
            and _cpp.enabled() and hasattr(_cpp.module(), "try_local_reorientation_coupled")):
        context.positions = np.array([p.position for p in toolpath.points], dtype=float)
        context.feedrates = np.array([p.feedrate for p in toolpath.points], dtype=float)
        context.angles = np.array([p.axis_angles for p in toolpath.points], dtype=float)
    return context


def _fast_checker(toolpath, machine_profile, sim_params):
    """Controleur rapide de step9, avec les MEMES echantillons de tete que
    step9 en mode rapide (zone de contact de la buse retiree)."""
    from element_ardor.step9_collision_sim import fast_checker_for, _filter_nozzle_contact_zone
    samples = _filter_nozzle_contact_zone(
        _surface_sample_points(machine_profile.nozzle_model, sim_params.toolhead_sample_spacing),
        sim_params.nozzle_contact_clearance,
    )
    return fast_checker_for(toolpath, machine_profile, sim_params, samples)


def _fast_events(checker, point_idx: int, point: ToolpathPoint) -> list[CollisionEvent]:
    plate, part = checker.penetrations(point_idx, point.position, point.axis_angles)
    events = []
    if part is not None:
        events.append(CollisionEvent(point_index=point_idx, obstacle_type="printed_part", penetration_depth=part))
    if plate is not None:
        events.append(CollisionEvent(point_index=point_idx, obstacle_type="platform", penetration_depth=plate))
    return events


def _coupled_candidates(axes, current_angles, reference_direction, samples):
    """Candidats de reorientation quand le premier axe est couple a la
    position : C fige a sa valeur courante, grille produit sur les axes
    INTERNES seulement. Rend (directions (K,3), angles (K, len(axes)))."""
    inner = axes[1:]
    grids = [np.linspace(spec.angle_range_deg[0], spec.angle_range_deg[1], samples) for spec in inner]
    mesh = np.meshgrid(*grids, indexing="ij")
    inner_angles = np.column_stack([m.ravel() for m in mesh])
    angles = np.column_stack([np.full(len(inner_angles), float(current_angles[0])), inner_angles])
    directions = np.array([_direction_for_angles(axes, list(a), reference_direction) for a in angles])
    return directions, angles


def _try_local_reorientation(
    toolpath: ToolpathMachine,
    point_idx: int,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    params: Step10Params,
    context: _ReorientationContext,
) -> Optional[list[float]]:
    """
    Cherche un jeu d'angles voisin qui elimine la collision au point
    `point_idx` sans depasser le seuil de surplomb.

    Candidats = grille de l'etape 8 (`composed_rotation_angles_and_directions`), essayes par ecart angulaire
    CROISSANT vis-a-vis de la direction d'origine ; arret des que l'ecart depasse
    `max_direction_deviation_deg` (le candidat reintroduirait un surplomb). Retourne les angles retenus, ou
    None. `context` est construit une seule fois par iteration externe de `run()` (voir
    `_ReorientationContext`), jamais ici.
    """
    axes = machine_profile.rotation_axes or []
    if not axes:
        return None  # aucune liberte de rotation : la reorientation n'a pas de sens

    point = toolpath.points[point_idx]
    # meme reference que step8 (sens d'impression a angles nuls = +Z) ; seuls des ECARTS entre directions
    # issues de cette meme reference sont compares ici.
    reference_direction = np.array([0.0, 0.0, 1.0])
    original_direction = _direction_for_angles(axes, point.axis_angles, reference_direction)

    coupled_case = list(getattr(machine_profile, "position_coupled_axes", None) or []) == [0] and len(axes) >= 2
    if coupled_case and context.angles is not None and context.fast_checker is not None:
        return _native_coupled_reorientation(point_idx, machine_profile, sim_params, params, context)

    if coupled_case:
        # Machine ACTUELLE (limite #59) : C est impose par la position du point, seul B (et les axes internes
        # suivants) peut varier ; grille 4x plus fine (un seul degre de liberte).
        directions, angles_grid = _coupled_candidates(
            axes, point.axis_angles, reference_direction, max(params.reorientation_samples_per_axis, 2) * 4
        )
    else:
        directions, angles_grid = composed_rotation_angles_and_directions(
            axes, reference_direction, samples_per_axis=max(params.reorientation_samples_per_axis, 2)
        )

    cosines = np.clip(directions @ original_direction, -1.0, 1.0)
    deviations_deg = np.degrees(np.arccos(cosines))

    has_toolhead = context.printed_index is not None
    if has_toolhead and context.fast_checker is None:
        context.grow_to(toolpath.points, point_idx)

    for candidate_idx in np.argsort(deviations_deg):
        if deviations_deg[candidate_idx] > params.max_direction_deviation_deg:
            break  # tries par ecart croissant : tout le reste depasse aussi le seuil

        candidate = [float(a) for a in angles_grid[candidate_idx]]
        if np.allclose(candidate, point.axis_angles, atol=1e-9):
            continue  # angles courants : c'est l'etat qui collisionne deja

        if _violates_rate_limits(
            toolpath.points, point_idx, candidate, machine_profile, params.rate_limit_neighbor_margin
        ):
            continue

        if not has_toolhead:
            return candidate  # aucun modele de tete : rien a re-verifier, meme garde que _collisions_at_points

        candidate_point = replace(point, axis_angles=list(candidate))
        if context.fast_checker is not None:
            events = _fast_events(context.fast_checker, point_idx, candidate_point)
        else:
            events = _evaluate_point(
                point_idx, context.toolhead_samples, candidate_point, machine_profile, context.printed_index,
                sim_params.safety_margin, printed_margin_of(sim_params)
            )
        if not events:
            return candidate

    return None


def _head_zlow_table(model, axes, axis_index: int, sim_params: CollisionSimParams, grid_step_deg: float):
    """Pour chaque valeur de l'axe `axis_index`, point le plus bas de la VRAIE tete sous la pointe.
    Retourne `(grid, zlow)`."""
    from element_ardor.step9_collision_sim import _filter_nozzle_contact_zone, _positioned_toolhead_vertices

    samples = _filter_nozzle_contact_zone(
        _surface_sample_points(model, sim_params.toolhead_sample_spacing), sim_params.nozzle_contact_clearance
    )
    lo, hi = axes[axis_index].angle_range_deg
    grid = np.arange(lo, hi + 1e-9, grid_step_deg)
    zlow = np.empty(len(grid))
    for g, value in enumerate(grid):
        angles = [0.0] * len(axes)
        angles[axis_index] = float(value)
        zlow[g] = _positioned_toolhead_vertices(samples, np.zeros(3), angles, axes)[:, 2].min()
    return grid, zlow


def _clamp_apply_native(toolpath, axis_index: int, grid, zlow, plate_z: float, sim_params: CollisionSimParams):
    zs = np.array([p.position[2] for p in toolpath.points], dtype=float)
    cur = np.array([p.axis_angles[axis_index] for p in toolpath.points], dtype=float)
    new_ax, mask = _cpp.module().clamp_head_height_apply(
        zs, cur, np.ascontiguousarray(grid, dtype=float), np.ascontiguousarray(zlow, dtype=float),
        plate_z, float(sim_params.safety_margin))
    points, deviations = [], []
    for i, p in enumerate(toolpath.points):
        if mask[i]:
            angles = list(p.axis_angles)
            angles[axis_index] = float(new_ax[i])
            deviations.append(abs(float(new_ax[i]) - float(cur[i])))
            p = replace(p, axis_angles=angles)
        points.append(p)
    return ToolpathMachine(points=points), int(mask.sum()), deviations


@dataclass
class _TierSetup:
    """Reglages communs de `reorient_coupled_tiers` (axes, grille, controleur, paliers, propagation)."""
    axes: list
    ax: int
    grid: np.ndarray
    samples: np.ndarray
    checker: object
    lo: float
    hi: float
    tiers_deg: tuple
    full_mm: float
    taper_mm: float
    sim_params: CollisionSimParams
    machine_profile: MachineProfile


def _reorient_tiers_native(setup: _TierSetup, pos, ext, new_angles, contact_indices):
    """Coeur de la boucle en natif (memes candidats, memes paliers, meme propagation)."""
    machine_profile, checker = setup.machine_profile, setup.checker
    axis_units = np.array([np.asarray(a.axis, dtype=float) / np.linalg.norm(a.axis) for a in setup.axes])
    disc = platform_disc(machine_profile)
    origin = machine_profile.platform_origin if disc is not None else None     # sans plateau decrit : test saute
    disc = disc or (0.0, 0.0, 0.0)
    new_angles, native_log = _cpp.module().reorient_coupled_tiers(
        np.ascontiguousarray(pos), np.ascontiguousarray(ext), np.ascontiguousarray(new_angles),
        np.ascontiguousarray(axis_units), int(setup.ax), np.ascontiguousarray(setup.grid, dtype=float),
        np.asarray(setup.tiers_deg, dtype=float), np.ascontiguousarray(setup.samples, dtype=float),
        checker.native, float(checker.reach),
        None if origin is None else np.asarray(origin, dtype=float),
        disc[0], disc[1], disc[2], float(setup.sim_params.safety_margin), float(setup.lo), float(setup.hi),
        float(setup.full_mm), float(setup.taper_mm),
        np.asarray(sorted(set(contact_indices)), dtype=np.int64),
    )
    log = [(int(i), float(a) if f else None, int(t) if f else None) for i, f, a, t in native_log]
    return new_angles, log


def _reorient_tiers_python(setup: _TierSetup, pos, ext, new_angles, contact_indices):
    """Boucle Python de reference : retourne `(new_angles, log)`."""
    from element_ardor.step9_collision_sim import _positioned_toolhead_vertices
    axes, ax, grid, checker = setup.axes, setup.ax, setup.grid, setup.checker
    z_unit = np.array([[0.0, 0.0, 1.0]])

    def head_dir(angles):
        return _positioned_toolhead_vertices(z_unit, np.zeros(3), angles, axes)[0]

    n_points = len(pos)
    original = new_angles.copy()   # a0 : angles d'ORIGINE du point, jamais ceux deja corriges par propagation
    log = []
    for idx in sorted(set(contact_indices)):
        a0 = list(original[idx])
        d0 = head_dir(a0)
        cands = sorted(grid, key=lambda b: abs(b - a0[ax]))
        devs = {}
        for b in cands:
            a = list(a0)
            a[ax] = float(b)
            devs[float(b)] = float(np.degrees(np.arccos(np.clip(head_dir(a) @ d0, -1.0, 1.0))))
        found = None
        for tier, tol in enumerate(setup.tiers_deg):
            for b in cands:
                dv = devs[float(b)]
                if dv > tol or (tier and dv <= setup.tiers_deg[tier - 1]):
                    continue
                a = list(a0)
                a[ax] = float(b)
                if checker.penetrations(idx, pos[idx], a) == (None, None):
                    found = (float(b), tier)
                    break
            if found:
                break
        log.append((idx, found[0] if found else None, found[1] if found else None))
        if not found:
            continue
        delta = found[0] - a0[ax]
        # voisins de la meme polyligne (points consecutifs en extrusion)
        start = idx
        while start > 0 and ext[start] > 0:
            start -= 1
        stop = idx
        while stop + 1 < n_points and ext[stop + 1] > 0:
            stop += 1
        seg = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pos[start:stop + 1], axis=0), axis=1))]
        s = np.abs(seg - seg[idx - start])
        w = np.clip((setup.taper_mm - s) / max(setup.taper_mm - setup.full_mm, 1e-9), 0.0, 1.0)
        near = np.flatnonzero(w > 0) + start
        new_angles[near, ax] = np.clip(new_angles[near, ax] + w[near - start] * delta, setup.lo, setup.hi)
        new_angles[idx, ax] = found[0]
    return new_angles, log


def _native_coupled_reorientation(
    point_idx: int, machine_profile: MachineProfile, sim_params: CollisionSimParams, params: Step10Params,
    context: _ReorientationContext,
) -> Optional[list[float]]:
    """Chemin natif (memes candidats par ecart croissant ; tri stable, la ou np.argsort ne garantit pas
    l'ordre des ex aequo)."""
    axes = machine_profile.rotation_axes or []
    axis_units = np.array([np.asarray(a.axis, dtype=float) / np.linalg.norm(a.axis) for a in axes])
    disc = platform_disc(machine_profile)
    origin = machine_profile.platform_origin if disc is not None else None
    disc = disc or (0.0, 0.0, 0.0)
    v_deg = float(np.degrees(machine_profile.max_angular_velocity)) if machine_profile.max_angular_velocity > 0 else 0.0
    found = _cpp.module().try_local_reorientation_coupled(
        context.positions, context.feedrates, context.angles, np.ascontiguousarray(axis_units),
        int(point_idx),
        np.array([a.angle_range_deg[0] for a in axes[1:]], dtype=float),
        np.array([a.angle_range_deg[1] for a in axes[1:]], dtype=float),
        int(max(params.reorientation_samples_per_axis, 2) * 4),
        float(params.max_direction_deviation_deg), v_deg, int(params.rate_limit_neighbor_margin),
        np.ascontiguousarray(context.toolhead_samples, dtype=float), context.fast_checker.native,
        float(context.fast_checker.reach),
        None if origin is None else np.asarray(origin, dtype=float),
        disc[0], disc[1], disc[2], float(sim_params.safety_margin),
    )
    return None if found is None else [float(v) for v in found]
