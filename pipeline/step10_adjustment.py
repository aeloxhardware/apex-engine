"""
pipeline/step10_adjustment.py

Etape 10 du pipeline : boucle d'ajustement iterative sur les collisions
rapportees par l'etape 9 et les orientations inatteignables rapportees par
l'etape 8.

Pour chaque probleme, tente dans l'ordre :
  1. REORIENTATION LOCALE : cherche un jeu d'angles voisin qui elimine la
     collision tout en restant sous le seuil de surplomb - via la MEME
     machinerie de recherche que step8 (composed_rotation_angles_and_directions
     + refine_axis_angles_local_search), jamais une nouvelle methode.
  2. RE-SEGMENTATION CIBLEE : emet une ResegmentationRequest vers l'appelant
     (hook injecte), avec la plage d'angles a exclure. Voir la docstring de
     ResegmentationRequest pour pourquoi l'etape 10 ne re-joue PAS
     elle-meme step3 -> step8.
  3. SUPPORT DE SECOURS : marque le point, avec log explicite.

Apres chaque ajustement, relance la verification de collision UNIQUEMENT
sur la portion modifiee (voir _collisions_at_points), pas la simulation
complete.

Sortie : (ToolpathMachine ajustee, AdjustmentReport).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable, Optional

import numpy as np
import trimesh

from core import _cpp
from core.geometry_utils import _rodrigues_rotate_single, composed_rotation_angles_and_directions
from core.mesh_types import (
    AdjustmentLogEntry,
    AdjustmentReport,
    CollisionEvent,
    CollisionReport,
    MachineProfile,
    ResegmentationRequest,
    RotationAxisSpec,
    Step8Diagnostic,
    ToolpathMachine,
    ToolpathPoint,
)
from pipeline.step9_collision_sim import (
    CollisionSimParams,
    _evaluate_point,
    _grow_printed_index_to,
    _filter_nozzle_contact_zone,
    _PrintedGeometryIndex,
    _surface_sample_points,
    run as collision_run,
)

# Un hook recoit une demande de re-segmentation et retourne la nouvelle
# trajectoire complete (step3 -> step8 re-joues par l'appelant), ou None
# s'il ne peut pas satisfaire la demande.
ResegmentationHook = Callable[[ResegmentationRequest], Optional[ToolpathMachine]]


@dataclass
class Step10Params:
    max_adjustment_iterations: int = 3
    # Grille de recherche de la reorientation locale - meme parametrage que
    # la recherche NN de step8 (grille produit sur les plages d'axes).
    reorientation_samples_per_axis: int = 12
    # Seuil de surplomb : ecart angulaire maximal tolere entre la direction
    # de depot d'origine et celle du candidat de reorientation. Un candidat
    # qui elimine la collision mais depasse ce seuil reintroduirait un
    # surplomb - il est rejete, c'est exactement la condition "tout en
    # respectant encore le seuil de surplomb" de la specification.
    max_direction_deviation_deg: float = 45.0
    # Demi-largeur de la plage d'angles signalee comme problematique dans
    # une ResegmentationRequest (autour de l'angle courant de chaque axe).
    excluded_range_half_width_deg: float = 15.0
    # Nombre de voisins de part et d'autre du point modifie sur lesquels
    # valider les limites de vitesse/acceleration angulaire.
    rate_limit_neighbor_margin: int = 1
    # Diagnostics d'ORIENTATION de l'etape 8 (ecart tete / direction visee,
    # sans collision) : seuls ceux dont l'ecart residuel depasse ce seuil
    # (degres) entrent dans la boucle (re-segmentation, puis support de
    # secours). Les autres sont journalises `orientation_deviation` et ne
    # comptent PAS comme supports. None (DEFAUT, decision de l'utilisateur
    # du 2026-10-07) = jamais : un support ne corrige pas une orientation.
    # 0.0 = comportement historique (tout diagnostic devient support :
    # 115 312 sur death contre 209 vraies collisions, limite #74).
    # Une collision reelle au meme point suit toujours la voie normale.
    orientation_support_min_residual_deg: Optional[float] = None


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
    from pipeline.step9_collision_sim import _filter_nozzle_contact_zone, _positioned_toolhead_vertices

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
    plate_z = float(machine_profile.platform_origin[2]) if machine_profile.platform_origin is not None else 0.0

    if _cpp.enabled() and hasattr(_cpp.module(), "clamp_head_height_apply") and len(toolpath.points):
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
    from pipeline.step9_collision_sim import (
        fast_checker_for,
        _filter_nozzle_contact_zone,
        _positioned_toolhead_vertices,
    )

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
    z_unit = np.array([[0.0, 0.0, 1.0]])

    def head_dir(angles):
        return _positioned_toolhead_vertices(z_unit, np.zeros(3), angles, axes)[0]

    pos = np.array([p.position for p in points], dtype=float)
    ext = np.array([p.extrusion for p in points], dtype=float)
    new_angles = np.array([p.axis_angles for p in points], dtype=float)
    log = []
    if checker.native is not None and hasattr(_cpp.module(), "reorient_coupled_tiers"):
        # Coeur de la boucle en natif (memes candidats, memes paliers, meme
        # propagation) ; seule la reconstruction des points reste ici.
        from pipeline.step9_collision_sim import (
            PLATFORM_CENTER_OFFSET, PLATFORM_HALF_HEIGHT, PLATFORM_RADIUS,
        )
        axis_units = np.array([np.asarray(a.axis, dtype=float) / np.linalg.norm(a.axis) for a in axes])
        radius = getattr(machine_profile, "platform_radius", None)
        origin = machine_profile.platform_origin
        new_angles, native_log = _cpp.module().reorient_coupled_tiers(
            np.ascontiguousarray(pos), np.ascontiguousarray(ext), np.ascontiguousarray(new_angles),
            np.ascontiguousarray(axis_units), int(ax), np.ascontiguousarray(grid, dtype=float),
            np.asarray(tiers_deg, dtype=float), np.ascontiguousarray(samples, dtype=float),
            checker.native, float(checker.reach),
            None if origin is None else np.asarray(origin, dtype=float),
            float(PLATFORM_RADIUS if radius is None else radius), float(PLATFORM_HALF_HEIGHT),
            float(PLATFORM_CENTER_OFFSET), float(sim_params.safety_margin), float(lo), float(hi),
            float(full_mm), float(taper_mm),
            np.asarray(sorted(set(contact_indices)), dtype=np.int64),
        )
        log = [(int(i), float(a) if f else None, int(t) if f else None) for i, f, a, t in native_log]
        out = [replace(p, axis_angles=[float(v) for v in new_angles[i]])
               if not np.array_equal(new_angles[i], p.axis_angles) else p for i, p in enumerate(points)]
        return ToolpathMachine(points=out), log
    for idx in sorted(set(contact_indices)):
        a0 = list(points[idx].axis_angles)
        d0 = head_dir(a0)
        cands = sorted(grid, key=lambda b: abs(b - a0[ax]))
        devs = {}
        for b in cands:
            a = list(a0)
            a[ax] = float(b)
            devs[float(b)] = float(np.degrees(np.arccos(np.clip(head_dir(a) @ d0, -1.0, 1.0))))
        found = None
        for tier, tol in enumerate(tiers_deg):
            for b in cands:
                dv = devs[float(b)]
                if dv > tol or (tier and dv <= tiers_deg[tier - 1]):
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
        while stop + 1 < len(points) and ext[stop + 1] > 0:
            stop += 1
        seg = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pos[start:stop + 1], axis=0), axis=1))]
        s = np.abs(seg - seg[idx - start])
        w = np.clip((taper_mm - s) / max(taper_mm - full_mm, 1e-9), 0.0, 1.0)
        near = np.flatnonzero(w > 0) + start
        new_angles[near, ax] = np.clip(new_angles[near, ax] + w[near - start] * delta, lo, hi)
        new_angles[idx, ax] = found[0]
    out = [replace(p, axis_angles=[float(v) for v in new_angles[i]]) if not np.array_equal(new_angles[i], p.axis_angles)
           else p for i, p in enumerate(points)]
    return ToolpathMachine(points=out), log


# ---------------------------------------------------------------------------
# Direction de depot pour un jeu d'angles donne
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Re-verification partielle (question 1 du plan)
# ---------------------------------------------------------------------------

def _collisions_at_points(
    toolpath: ToolpathMachine,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    point_indices: list[int],
) -> list[CollisionEvent]:
    """
    Re-teste la collision UNIQUEMENT aux indices demandes, en reutilisant
    les primitives de l'etape 9 (jamais une seconde implementation de la
    detection).

    Pourquoi c'est exact et pourquoi la "portion modifiee" se reduit aux
    seuls points reorientes : l'index de matiere imprimee de l'etape 9 est
    construit a partir de `points[i].position` UNIQUEMENT, jamais de
    `axis_angles` (voir _grow_printed_index_to). Une reorientation locale ne
    touche que les angles, donc la geometrie imprimee est strictement
    identique avant/apres, et aucun point en aval ne peut voir son etat de
    collision changer. Seul le point reoriente doit etre re-teste.

    L'index est reconstruit jusqu'au plus grand indice demande (croissance
    strictement monotone, meme invariant que l'etape 9 : jamais en avance
    sur le point teste, cordon en cours de depot exclu).
    """
    nozzle_model = machine_profile.nozzle_model
    if not isinstance(nozzle_model, trimesh.Trimesh) or len(nozzle_model.vertices) == 0:
        return []
    if not point_indices:
        return []

    # toolhead_sample_spacing (echelle de la TETE), jamais nozzle_diameter
    # (echelle du CORDON) - meme distinction que step9 (CLAUDE.md limite
    # #26). Utiliser nozzle_diameter ici reproduirait exactement le meme
    # blocage (17734 points de test par appel au lieu de ~1540, mesure
    # concretement sur aetheris-4.glb), simplement deplace de l'etape 9 a
    # l'etape 10 plutot qu'elimine.
    if sim_params.fast_kdtree:
        checker = _fast_checker(toolpath, machine_profile, sim_params)
        indices = [i for i in sorted(set(point_indices)) if i < len(toolpath.points)]
        from pipeline.step9_collision_sim import fast_batch_events
        batch = fast_batch_events(checker, toolpath, indices)
        if batch is not None:
            return batch
        return [e for i in indices for e in _fast_events(checker, i, toolpath.points[i])]

    # Zone de contact de la buse retiree, EXACTEMENT comme step9 : sans ce
    # filtre (oublie ici jusqu'au 2026-10-01), la pointe touche toujours le
    # cordon qu'elle depose et tout point re-teste etait declare en collision.
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
            )
        )
    return events


# ---------------------------------------------------------------------------
# Niveau 1 - Reorientation locale
# ---------------------------------------------------------------------------

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
    Etat partage entre TOUS les points traites par le niveau 1 (reorientation
    locale) d'une meme iteration externe de `run()` - jamais reconstruit par
    point ni par candidat.

    Les echantillons de tete ne dependent que du modele (fixe pour toute la
    boucle). L'index de matiere imprimee ne depend que de `position` (jamais
    de `axis_angles`, voir `_collisions_at_points`) et grandit de facon
    strictement monotone au fil de `problem_indices` traites en ordre
    CROISSANT par `run()` - exactement le meme invariant que la croissance
    de l'etape 9 (`_grow_printed_index_to` dans `step9_collision_sim.run`).

    Sans cet etat partage, `_try_local_reorientation` reconstruisait son
    propre index depuis 0 a CHAQUE point traite (jusqu'a 3588 fois sur la
    trajectoire de test), chacun re-inserant toutes les capsules jusqu'a son
    propre `point_idx` - mesure concretement : le taux d'appel chutait de
    ~2.5 appels/s au debut a 0.39 appels/s au 200e appel (elapsed 512s),
    une degradation croissante avec la position dans la trajectoire,
    incompatible avec le cout lineaire vise.
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
    from pipeline.step9_collision_sim import fast_checker_for, _filter_nozzle_contact_zone
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

    Reutilise la grille de recherche de l'etape 8
    (`composed_rotation_angles_and_directions`) plutot qu'une nouvelle
    exploration : les candidats sont exactement les orientations
    atteignables par la cinematique decrite, pour un nombre d'axes
    quelconque. Ils sont essayes par ecart angulaire CROISSANT vis-a-vis de
    la direction d'origine (le plus petit changement d'abord), et la
    recherche s'arrete des que l'ecart depasse `max_direction_deviation_deg`
    - au-dela, le candidat reintroduirait un surplomb.

    Retourne les angles retenus, ou None si aucun candidat ne convient.

    Performance : `context` (echantillons de tete + index de matiere
    imprimee) est fourni par l'appelant, construit UNE SEULE FOIS par
    iteration externe de `run()` et grossi de facon incrementale au fil des
    points traites (voir `_ReorientationContext`) - jamais reconstruit ici,
    ni par point ni par candidat.
    """
    axes = machine_profile.rotation_axes or []
    if not axes:
        return None  # aucune liberte de rotation : la reorientation n'a pas de sens

    point = toolpath.points[point_idx]
    # meme reference que step8 (sens d'impression a angles nuls = +Z). Seuls
    # des ECARTS entre directions issues de cette meme reference sont compares
    # ici : le comportement ne dependait pas de son signe, l'alignement est
    # une question de coherence.
    reference_direction = np.array([0.0, 0.0, 1.0])
    original_direction = _direction_for_angles(axes, point.axis_angles, reference_direction)

    coupled_case = list(getattr(machine_profile, "position_coupled_axes", None) or []) == [0] and len(axes) >= 2
    if coupled_case and context.angles is not None and context.fast_checker is not None:
        # Chemin natif (memes candidats par ecart croissant ; tri stable, la ou
        # np.argsort ne garantit pas l'ordre des ex aequo).
        from pipeline.step9_collision_sim import PLATFORM_CENTER_OFFSET, PLATFORM_HALF_HEIGHT, PLATFORM_RADIUS
        axis_units = np.array([np.asarray(a.axis, dtype=float) / np.linalg.norm(a.axis) for a in axes])
        radius = getattr(machine_profile, "platform_radius", None)
        origin = machine_profile.platform_origin
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
            float(PLATFORM_RADIUS if radius is None else radius), float(PLATFORM_HALF_HEIGHT),
            float(PLATFORM_CENTER_OFFSET), float(sim_params.safety_margin),
        )
        return None if found is None else [float(v) for v in found]

    if coupled_case:
        # Machine ACTUELLE (limite #59) : C est impose par la position du
        # point, seul B (et les axes internes suivants) peut varier. La grille
        # produit complete proposait des C inexecutables (corrige le
        # 2026-10-01). B seul : grille 4x plus fine (un seul degre de liberte).
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
    toolhead_samples = context.toolhead_samples
    printed_index = context.printed_index

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
                point_idx, toolhead_samples, candidate_point, machine_profile, printed_index,
                sim_params.safety_margin
            )
        if not events:
            return candidate

    return None


# ---------------------------------------------------------------------------
# Niveau 2 - Demande de re-segmentation ciblee
# ---------------------------------------------------------------------------

def _excluded_ranges_for_point(
    point: ToolpathPoint, axes: list[RotationAxisSpec], half_width_deg: float
) -> list[Optional[tuple[float, float]]]:
    """
    Plage d'angles a exclure par axe : une fenetre centree sur l'angle
    courant, celui qui a mene a la collision. Formule ainsi (et pas comme
    "un angle unique interdit") parce que la segmentation raisonne en
    directions candidates continues, pas en valeurs discretes.
    """
    return [
        (angle - half_width_deg, angle + half_width_deg)
        for angle, _spec in zip(point.axis_angles, axes)
    ]


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    toolpath: ToolpathMachine,
    collision_report: CollisionReport,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    params: Step10Params,
    step8_diagnostics: Optional[list[Step8Diagnostic]] = None,
    resegmentation_hook: Optional[ResegmentationHook] = None,
) -> tuple[ToolpathMachine, AdjustmentReport]:
    """
    Boucle d'ajustement sur les collisions (etape 9) et les orientations
    inatteignables (etape 8).

    Args:
        toolpath: sortie de l'etape 8 (jamais mutee - une copie ajustee est
            retournee).
        collision_report: sortie de l'etape 9 sur cette trajectoire.
        machine_profile: cinematique + `nozzle_model` deja charge.
        sim_params: parametres de l'etape 9, reutilises tels quels pour
            chaque re-verification (meme marge, meme echantillonnage).
        params: parametres propres a l'ajustement.
        step8_diagnostics: orientations inatteignables signalees par
            l'etape 8. Elles sautent le niveau 1 : la recherche IK locale
            a deja ete epuisee par l'etape 8 pour ces points, la relancer
            a l'identique ne peut rien donner de nouveau.
        resegmentation_hook: appele quand la reorientation locale echoue.
            Absent => le niveau 2 est saute et l'on passe directement au
            support de secours.

    Returns:
        (trajectoire ajustee, AdjustmentReport). Ne leve JAMAIS pour une
        collision non resolue : a l'epuisement de `max_adjustment_iterations`,
        les collisions restantes basculent en support de secours et
        `converged` passe a False - signal explicite, coherent avec
        `Step8Diagnostic` et `CollisionReport.chuck_checked`.
    """
    current = ToolpathMachine(points=list(toolpath.points))
    entries: list[AdjustmentLogEntry] = []
    fallback_marked: set[int] = set()
    axes = machine_profile.rotation_axes or []

    tracked_indices = sorted({event.point_index for event in collision_report.events})
    diagnostics = list(step8_diagnostics or [])
    threshold = params.orientation_support_min_residual_deg
    if threshold is None:
        threshold = float("inf")
    if diagnostics:
        # Ecart d'orientation sous le seuil : journalise, jamais un support
        # (voir Step10Params.orientation_support_min_residual_deg).
        logged: set[int] = set()
        for diagnostic in diagnostics:
            if diagnostic.residual_angle_error_deg <= threshold and diagnostic.point_index not in logged:
                logged.add(diagnostic.point_index)
                entries.append(
                    AdjustmentLogEntry(
                        point_index=diagnostic.point_index,
                        method="orientation_deviation",
                        detail=(
                            f"ecart d'orientation {diagnostic.residual_angle_error_deg:.1f} deg "
                            f"<= {threshold:g} deg : depot incline, sans support"
                        ),
                    )
                )
        diagnostics = [d for d in diagnostics if d.residual_angle_error_deg > threshold]
    diagnostic_indices = sorted({diagnostic.point_index for diagnostic in diagnostics})

    pending_events = list(collision_report.events)
    iterations_used = 0
    converged = False

    for iteration in range(1, max(params.max_adjustment_iterations, 0) + 1):
        iterations_used = iteration

        problem_indices = [
            index
            for index in sorted({event.point_index for event in pending_events} | set(diagnostic_indices))
            if index not in fallback_marked and index < len(current.points)
        ]
        if not problem_indices:
            converged = True
            break

        full_recheck_needed = False

        # Construit une seule fois par iteration externe, grossi ensuite de
        # facon monotone au fil de `problem_indices` (deja trie CROISSANT
        # ci-dessus) - voir _ReorientationContext.
        reorientation_context = _build_reorientation_context(current, machine_profile, sim_params)

        # Ensembles, pas listes : avec la normale par point, death compte ~115 000
        # diagnostics d'etape 8, et un `in` sur la liste triee rendait cette boucle
        # quadratique (73 s sur 81 s de step10, 2026-10-06). Meme resultat.
        diagnostic_set = set(diagnostic_indices)
        pending_event_indices = {event.point_index for event in pending_events}
        for point_index in problem_indices:
            is_diagnostic_only = point_index in diagnostic_set and point_index not in pending_event_indices

            # --- Niveau 1 : reorientation locale ------------------------------
            if not is_diagnostic_only:
                new_angles = _try_local_reorientation(
                    current, point_index, machine_profile, sim_params, params, reorientation_context
                )
                if new_angles is not None:
                    current = _with_angles(current, point_index, new_angles)
                    if reorientation_context.angles is not None:
                        reorientation_context.angles[point_index] = new_angles
                    entries.append(
                        AdjustmentLogEntry(
                            point_index=point_index,
                            method="local_reorientation",
                            detail=f"angles ajustes vers {[round(a, 3) for a in new_angles]}",
                        )
                    )
                    continue

            # --- Niveau 2 : re-segmentation ciblee ----------------------------
            if resegmentation_hook is not None:
                request = ResegmentationRequest(
                    point_index=point_index,
                    position=current.points[point_index].position,
                    excluded_axis_ranges_deg=_excluded_ranges_for_point(
                        current.points[point_index], axes, params.excluded_range_half_width_deg
                    ),
                    reason=(
                        "orientation inatteignable (etape 8)"
                        if is_diagnostic_only
                        else "collision non resolue par reorientation locale"
                    ),
                )
                new_toolpath = resegmentation_hook(request)
                if new_toolpath is not None:
                    current = new_toolpath
                    entries.append(
                        AdjustmentLogEntry(
                            point_index=point_index,
                            method="region_resegmentation",
                            detail=request.reason,
                        )
                    )
                    # Une re-segmentation invalide TOUTE la trajectoire (champ
                    # scalaire re-resolu globalement) : la re-verification
                    # partielle ne s'applique plus, il faut une passe complete.
                    full_recheck_needed = True
                    diagnostic_indices = []
                    break

            # --- Niveau 3 : support de secours --------------------------------
            fallback_marked.add(point_index)
            entries.append(
                AdjustmentLogEntry(
                    point_index=point_index,
                    method="fallback_support",
                    detail=(
                        "ni reorientation locale ni re-segmentation n'ont resolu "
                        "la collision - support classique requis"
                    ),
                )
            )

        if full_recheck_needed:
            pending_events = [
                event
                for event in collision_run(current, machine_profile, sim_params).events
                if event.point_index not in fallback_marked
            ]
            tracked_indices = sorted({event.point_index for event in pending_events})
        else:
            pending_events = [
                event
                for event in _collisions_at_points(
                    current, machine_profile, sim_params, tracked_indices
                )
                if event.point_index not in fallback_marked
            ]

        if not pending_events and not [
            index for index in diagnostic_indices if index not in fallback_marked
        ]:
            converged = True
            break

    # Epuisement des iterations avec des collisions encore ouvertes : bascule
    # forcee en support de secours plutot qu'une exception (voir docstring).
    unresolved = [event for event in pending_events if event.point_index not in fallback_marked]
    for event in unresolved:
        fallback_marked.add(event.point_index)
        entries.append(
            AdjustmentLogEntry(
                point_index=event.point_index,
                method="fallback_support",
                detail=(
                    f"non resolu apres {iterations_used} iteration(s) "
                    f"(max_adjustment_iterations atteint) - support classique requis"
                ),
            )
        )
    if unresolved:
        converged = False

    return current, AdjustmentReport(
        entries=entries,
        converged=converged,
        unresolved_collisions=unresolved,
        iterations_used=iterations_used,
    )
