"""
pipeline/step8_machine_toolpath.py

Etape 8 du pipeline : resolution de la trajectoire machine reelle a partir
de PrintSchedule (etape 7) et des couches recolorees qu'il reference.

Pour chaque point de trajectoire, ce fichier :
  1. Releve la position 2D locale vers le 3D (`step6._lift_to_3d`).
  2. Resout la direction d'impression desiree (normale du repere local de
     la couche) en angles concrets pour CHAQUE axe de
     `MachineProfile.rotation_axes` - generiquement, pour N axes, jamais 2
     codes en dur (voir CLAUDE.md limite de conception #18).
  3. Detecte et gere la singularite au centre de rotation (distance nulle a
     l'axe -> angle herite du point precedent plutot que recalcule).
  4. Signale (sans jamais lever d'exception) toute orientation inatteignable
     dans les plages d'axes disponibles - donnee pure vers l'etape 10.
  5. Lisse la serie temporelle d'angles pour respecter
     max_angular_velocity/max_angular_acceleration.
  6. Insere des trajets a vide (extrusion=0) aux frontieres de polyligne.

Ce fichier ne connait RIEN d'un firmware particulier : sa sortie est
position cartesienne + angles generiques + extrusion/feedrate. La
traduction vers une commande machine precise (ordre X/Y/Z/B, noms d'axes)
est repoussee vers un futur module d'export (io/gcode_export.py), hors
perimetre ici.

Depend de list[Layer] (retour de step7, toolpath_2d re-decoupe par couleur),
PrintSchedule (step7, ordre d'execution gele), MachineProfile.

Sortie : (ToolpathMachine, list[Step8Diagnostic]).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from core import _cpp
from core.geometry_utils import (
    _rodrigues_rotate,
    composed_rotation_angles_and_directions,
    refine_axis_angles_local_search,
    resolve_axis_angles_analytic,
)
from core.mesh_types import (
    Layer,
    MachineProfile,
    PrintPass,
    PrintSchedule,
    RotationAxisSpec,
    Step8Diagnostic,
    ToolpathMachine,
    ToolpathPoint,
)
from pipeline.step6_toolpath2d import _lift_to_3d


class ToolpathResolutionError(Exception):
    """
    Levee uniquement pour une entree structurellement invalide (toolpath_2d
    manquant, schedule pointant vers un segment/couche inexistant) - JAMAIS
    pour une orientation inatteignable, qui est un cas attendu couvert par
    Step8Diagnostic, pas une erreur bloquante.
    """


@dataclass
class Step8Params:
    # Nombre d'echantillons par axe pour la recherche du plus proche voisin
    # (composed_rotation_angles_and_directions) - cout total = ~samples**N.
    ik_samples_per_axis: int = 10
    ik_refine_iterations: int = 8
    # Multiplicateur du pas initial d'affinement, PAS une constante fixe en
    # degres : le pas de depart = ik_refine_step_multiplier * (360 /
    # ik_samples_per_axis), a l'echelle de l'ecart entre deux points voisins
    # de la grille NN (jusqu'a ~36deg pour samples_per_axis=10). Un pas fixe
    # trop petit laisse l'affinement piege pres du candidat NN au lieu de
    # converger vers le vrai optimum (mesure : residu ~15-20deg avec un pas
    # fixe de 2deg ; <0.1deg avec ce pas derive - voir tests).
    ik_refine_step_multiplier: float = 0.5
    # Distance en-deca de laquelle la position courante est consideree SUR
    # l'axe de rotation (singularite : angle indetermine, herite du point
    # precedent plutot que recalcule).
    singularity_radius_epsilon: float = 1e-3
    # Feedrate des trajets a vide inseres aux frontieres de polyligne. 0 =>
    # reutilise le feedrate du point suivant (pas de vitesse magique par
    # defaut non documentee).
    travel_feedrate: float = 0.0
    default_feedrate: float = 1000.0
    # Au-dela de ce seuil, l'erreur residuelle d'un point est jugee assez
    # grande pour meriter un Step8Diagnostic (signal vers l'etape 10).
    max_unreachable_angular_error_deg: float = 5.0
    # Troncature des angles par les limites de vitesse/acceleration du profil
    # (_clamp_rate). DESACTIVEE par defaut depuis le 2026-09-28 : le firmware
    # de la machine COORDONNE les axes (confirme par l'utilisateur) -- quand un
    # axe atteint sa limite, c'est tout le mouvement qui ralentit, l'angle
    # reste exact. Tronquer l'angle fait prendre a la tete un retard qu'elle
    # n'aura pas (mesure en recherche : jusqu'a 62 deg d'ecart a la normale sur
    # la coque du Benchy). Le bon modele est un ralentissement de l'avance, a
    # faire a l'export G-code (feuille de route). Les limites restent dans le
    # profil machine ; True retrouve l'ancien comportement.
    clamp_angle_rates: bool = False
    # Filtre MEDIAN (en nombre de points, 0 = aucun) des axes d'orientation
    # NON couples a la position, par polyligne (2026-09-30, mesure en
    # recherche). La normale de couche est constante par tetraedre : pres de
    # l'axe du plateau, C etant impose, B oscille 0 <-> ~14 deg d'un point a
    # l'autre (5,9 % des segments du corps de Squirtle), et le firmware, qui
    # coordonne les axes, ralentit tout le mouvement a chaque inversion
    # (x3,7 sur le temps du corps). 9 points : -41 % de temps. Opt-in : a
    # appliquer AVANT tout controle de collision (step9 valide les angles
    # filtres). L'ecart residuel a la direction visee est recalcule.
    free_axis_median_window: int = 0
    # Normale PAR POINT (phase 3.1, session du 2026-10-01). Par defaut, la cible
    # de l'IK est la normale du plan MOYEN de la couche, une seule pour toute la
    # couche. Sur la machine actuelle (C impose par la position, limite #59), une
    # normale inclinee fixe n'est dans le plan radial que d'une poignee de
    # points : mesure, 95 % des points de l_shape en `Step8Diagnostic`. Vrai =
    # la cible devient le gradient LOCAL du champ en chaque point, fourni par
    # l'orchestrateur (`run(..., direction_field=...)`) ; sans champ fourni,
    # comportement historique. Meme choix que la recherche (full_step89 :
    # « normale du corps = gradient de phi »).
    per_point_normal: bool = False
    # Variante continue de la normale par point (gradient moyenne aux noeuds,
    # interpole en barycentrique) ; sans effet si per_point_normal est faux.
    per_point_normal_nodal: bool = False


# ---------------------------------------------------------------------------
# Direction desiree et resolution cinematique inverse generique
# ---------------------------------------------------------------------------

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


def _resolve_axis_angles_via_grid_search(
    desired_direction: np.ndarray,
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    params: Step8Params,
) -> tuple[list[float], float]:
    """
    Recherche par grille + affinement local - methode HISTORIQUE, conservee
    comme repli de `_resolve_axis_angles` (voir sa docstring) pour le cas
    rare ou la resolution analytique ne trouve geometriquement aucune
    solution (paire externe incapable d'atteindre l'angle cible, verifie
    reel et non hypothetique pendant la conception - CLAUDE.md/plan).

      1. Recherche du plus proche voisin parmi les directions atteignables
         par composed_rotation_angles_and_directions (grille produit).
      2. Affinement local generique (refine_axis_angles_local_search),
         avec un pas initial derive de l'ecart de la grille NN - voir
         Step8Params.ik_refine_step_multiplier.

    Retourne (angles_deg, erreur_angulaire_residuelle_rad) - l'erreur est
    l'angle entre `desired_direction` et la direction reellement atteinte
    par les angles resolus (recalculee par cinematique directe), jamais
    supposee nulle.
    """
    if not axes:
        return [], 0.0

    directions, angles_grid = composed_rotation_angles_and_directions(
        axes, reference_direction, samples_per_axis=max(params.ik_samples_per_axis, 1)
    )

    desired = np.asarray(desired_direction, dtype=float)
    desired_norm = np.linalg.norm(desired)
    if desired_norm > 1e-12:
        desired = desired / desired_norm

    similarity = directions @ desired
    best_index = int(np.argmax(similarity))
    initial_angles = angles_grid[best_index]

    grid_spacing_deg = 360.0 / max(params.ik_samples_per_axis, 1)
    initial_step_deg = params.ik_refine_step_multiplier * grid_spacing_deg

    refined_angles, residual_error_rad = refine_axis_angles_local_search(
        axes,
        reference_direction,
        initial_angles,
        desired,
        iterations=params.ik_refine_iterations,
        initial_step_deg=initial_step_deg,
    )

    return [float(a) for a in refined_angles], float(residual_error_rad)


def _resolve_axis_angles(
    desired_direction: np.ndarray,
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    params: Step8Params,
    previous_angles: Optional[list[float]] = None,
) -> tuple[list[float], float]:
    """
    Resout `desired_direction` en un angle par axe de `axes`, generiquement
    pour N axes (0, 1, 2, ...), par cinematique inverse ANALYTIQUE
    (`core.geometry_utils.resolve_axis_angles_analytic` - bissection + atan2,
    voir sa docstring pour la justification geometrique complete et
    CLAUDE.md "Piste chiffree pour step8" pour le diagnostic qui a motive
    ce fix : la grille historique avait une erreur mediane mesuree de
    9.26 deg contre un seuil de rejet de 5 deg, echouant structurellement
    son propre critere).

    Repli vers `_resolve_axis_angles_via_grid_search` (methode historique,
    conservee intacte) UNIQUEMENT quand la resolution analytique ne trouve
    geometriquement aucune solution pour la paire d'axes externe - cas
    confirme reel (pas hypothetique) pendant la conception sur des axes
    non orthogonaux, jamais rencontre sur la machine a 2 axes orthogonaux
    actuelle.

    `previous_angles` (angles resolus au point de trajectoire precedent,
    None au premier point) sert de valeur de repliement pour les axes
    internes au-dela de la paire externe (axes[2:]) et de repli en cas de
    singularite locale - voir docstring de resolve_axis_angles_analytic.

    Retourne (angles_deg, erreur_angulaire_residuelle_rad) - l'erreur est
    l'angle entre `desired_direction` et la direction reellement atteinte
    par les angles resolus (recalculee par cinematique directe), jamais
    supposee nulle.
    """
    if not axes:
        return [], 0.0

    def _fallback() -> tuple[list, float]:
        return _resolve_axis_angles_via_grid_search(
            desired_direction, axes, reference_direction, params
        )

    angles, residual_error_rad = resolve_axis_angles_analytic(
        axes, reference_direction, desired_direction, previous_angles, _fallback
    )

    return [float(a) for a in angles], float(residual_error_rad)


# ---------------------------------------------------------------------------
# Singularite au centre de rotation
# ---------------------------------------------------------------------------

def _axis_chain_world_lines(
    axes: list[RotationAxisSpec], angles_deg: list[float], platform_origin: np.ndarray
) -> list[np.ndarray]:
    """
    Direction monde (unitaire) de l'axe de rotation de chaque DOF, apres
    application des rotations des DOF qui le precedent dans la chaine, aux
    angles courants `angles_deg`. Toutes les rotations partagent la meme
    origine `platform_origin` (RotationAxisSpec n'a pas d'origine propre par
    axe dans ce modele) : seule la DIRECTION de chaque droite de rotation
    varie le long de la chaine, pas son point d'ancrage.

    Necessaire generiquement (pas seulement pour l'axe 0) car un axe interne
    (ex: la tete B) peut lui-meme devenir singulier selon la geometrie de la
    machine decrite - voir _detect_singular_axes.
    """
    lines: list[np.ndarray] = []
    accumulated_rotation_angles: list[float] = []
    accumulated_axes: list[np.ndarray] = []

    for axis_index, spec in enumerate(axes):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        # applique les rotations des DOF PRECEDENTS (exterieurs) a cet axe
        world_axis = axis.reshape(1, 3)
        for prior_axis, prior_angle_deg in zip(
            reversed(accumulated_axes), reversed(accumulated_rotation_angles)
        ):
            world_axis = _rodrigues_rotate(
                world_axis, prior_axis, np.array([np.radians(prior_angle_deg)])
            )
        lines.append(world_axis[0])
        accumulated_axes.append(axis)
        accumulated_rotation_angles.append(angles_deg[axis_index])

    return lines


def _detect_singular_axes(
    position: np.ndarray,
    platform_origin: Optional[np.ndarray],
    axes: list[RotationAxisSpec],
    angles_deg: list[float],
    epsilon: float,
) -> list[int]:
    """
    Indices des axes dont la droite de rotation (ancree a `platform_origin`,
    orientee selon `_axis_chain_world_lines`) passe a moins de `epsilon` de
    `position` - generalise "R≈0 => C indetermine" a "distance a l'axe i≈0
    => angle i indetermine", pour n'importe quel axe de la chaine, pas
    seulement le premier.
    """
    if platform_origin is None or not axes:
        return []

    lines = _axis_chain_world_lines(axes, angles_deg, platform_origin)
    relative = np.asarray(position, dtype=float) - np.asarray(platform_origin, dtype=float)

    singular: list[int] = []
    for axis_index, line_direction in enumerate(lines):
        projection = np.dot(relative, line_direction) * line_direction
        perpendicular_distance = float(np.linalg.norm(relative - projection))
        if perpendicular_distance < epsilon:
            singular.append(axis_index)
    return singular


def _position_to_radius(
    position: np.ndarray, platform_origin: Optional[np.ndarray], axes: list[RotationAxisSpec]
) -> Optional[float]:
    """
    Distance perpendiculaire de `position` a la droite portee par l'axe de
    rotation LE PLUS EXTERNE (axes[0]), passant par `platform_origin` -
    concept de "rayon" bien defini des qu'il existe au moins un axe de
    rotation ancre a une origine. Pour une machine a 0 axe de rotation
    (orientation libre), aucun rayon n'a de sens : retourne None plutot
    qu'une valeur inventee.
    """
    if platform_origin is None or not axes:
        return None

    axis = np.asarray(axes[0].axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    relative = np.asarray(position, dtype=float) - np.asarray(platform_origin, dtype=float)
    projection = np.dot(relative, axis) * axis
    return float(np.linalg.norm(relative - projection))


def _inherit_angle_at_singularity(
    singular_axes: list[int], resolved_angles: list[float], previous_angles: Optional[list[float]]
) -> list[float]:
    """
    Remplace l'angle de chaque axe signale singulier par l'angle du POINT
    PRECEDENT plutot que par la valeur recalculee (potentiellement instable,
    type atan2(0,0)) - generalise l'heritage a autant d'axes simultanement
    singuliers que necessaire, pas seulement l'axe 0.
    """
    if not singular_axes or previous_angles is None:
        return resolved_angles

    angles = list(resolved_angles)
    for axis_index in singular_axes:
        if axis_index < len(previous_angles):
            angles[axis_index] = previous_angles[axis_index]
    return angles


# ---------------------------------------------------------------------------
# Lissage temporel (vitesse/acceleration angulaire)
# ---------------------------------------------------------------------------

def _is_full_turn_axis(spec: RotationAxisSpec) -> bool:
    lo, hi = spec.angle_range_deg
    return (hi - lo) >= 360.0


def _unwrap_angle_sequence(angles_deg: np.ndarray, is_full_turn: bool) -> np.ndarray:
    """
    Convertit une sequence brute en sequence continue (deplie) via le delta
    au plus court chemin - UNIQUEMENT pour un axe a tour complet (359deg ->
    1deg doit se lire comme un pas de 2deg, pas 358deg). Une plage bornee
    strictement (ex: B in [0,180]) ne peut par construction pas boucler :
    depli desactive (no-op), le depliage y introduirait une discontinuite
    qui n'existe pas physiquement.
    """
    if len(angles_deg) == 0:
        return angles_deg
    if not is_full_turn:
        return angles_deg.copy()

    unwrapped = np.empty_like(angles_deg)
    unwrapped[0] = angles_deg[0]
    for i in range(1, len(angles_deg)):
        delta = ((angles_deg[i] - angles_deg[i - 1] + 180.0) % 360.0) - 180.0
        unwrapped[i] = unwrapped[i - 1] + delta
    return unwrapped


def _rewrap_angle_sequence(
    unwrapped_deg: np.ndarray, angle_range_deg: tuple[float, float], is_full_turn: bool
) -> np.ndarray:
    """Inverse de _unwrap_angle_sequence : replie dans la plage physique valide apres lissage."""
    if len(unwrapped_deg) == 0:
        return unwrapped_deg
    if is_full_turn:
        lo, _hi = angle_range_deg
        return lo + np.mod(unwrapped_deg - lo, 360.0)
    lo, hi = angle_range_deg
    return np.clip(unwrapped_deg, lo, hi)


def _clamp_rate(
    unwrapped_deg: np.ndarray,
    dt_per_point: np.ndarray,
    max_velocity_rad_s: float,
    max_accel_rad_s2: float,
) -> np.ndarray:
    """
    Limiteur de taux causal (avant uniquement) : applique sequentiellement
    les contraintes |delta_angle/delta_t| <= max_velocity puis
    |delta_vitesse/delta_t| <= max_accel entre points consecutifs.

    Prefere a une spline : les limites machine sont des contraintes DURES
    que ce filtre respecte par construction (clamp direct), alors qu'une
    spline necessiterait un ajustement iteratif sans garantie de
    convergence sous contrainte de derivees. `dt_per_point[i]` est le temps
    entre le point i-1 et i (issu du feedrate/distance, pas suppose
    uniforme) - le filtre reste donc generique sur un espacement de points
    quelconque, coherent avec l'ordre gele et non re-echantillonne de
    PrintSchedule.

    Si max_velocity_rad_s ou max_accel_rad_s2 <= 0 (non specifie), aucune
    limite n'est appliquee sur la grandeur correspondante - une machine
    pour laquelle ces limites ne sont pas encore connues ne doit pas voir
    sa trajectoire tronquee a zero mouvement.
    """
    n = len(unwrapped_deg)
    if n <= 1:
        return unwrapped_deg.copy()

    clamped = unwrapped_deg.copy()
    previous_velocity_deg_s = 0.0

    for i in range(1, n):
        dt = max(float(dt_per_point[i]), 1e-9)
        desired_delta_deg = clamped[i] - clamped[i - 1]
        desired_velocity_deg_s = desired_delta_deg / dt

        if max_velocity_rad_s > 0:
            max_velocity_deg_s = np.degrees(max_velocity_rad_s)
            desired_velocity_deg_s = float(
                np.clip(desired_velocity_deg_s, -max_velocity_deg_s, max_velocity_deg_s)
            )

        if max_accel_rad_s2 > 0:
            max_accel_deg_s2 = np.degrees(max_accel_rad_s2)
            max_delta_velocity = max_accel_deg_s2 * dt
            desired_velocity_deg_s = float(
                np.clip(
                    desired_velocity_deg_s,
                    previous_velocity_deg_s - max_delta_velocity,
                    previous_velocity_deg_s + max_delta_velocity,
                )
            )

        clamped[i] = clamped[i - 1] + desired_velocity_deg_s * dt
        previous_velocity_deg_s = desired_velocity_deg_s

    return clamped


# ---------------------------------------------------------------------------
# Parcours du schedule et insertion des trajets a vide
# ---------------------------------------------------------------------------

@dataclass
class _RawPoint:
    position: np.ndarray
    desired_direction: np.ndarray
    layer_id: int
    segment_index: int
    color_id: Optional[int]
    feedrate: float
    extrusion: float
    is_travel: bool


def _iter_schedule_points(
    schedule: PrintSchedule, layers_by_id: dict[int, Layer], params: Step8Params
) -> list[_RawPoint]:
    """
    Parcourt `schedule.passes` DANS L'ORDRE (jamais re-derive : c'est
    l'ordre d'execution gele par step7), releve chaque segment 2D vers le
    3D, et insere un point de trajet a vide (`extrusion=0.0`) a chaque
    frontiere de polyligne - y compris entre deux passes non contigues
    (changement de couleur ou saut de couche). Une frontiere est detectee
    des que la cle (layer_id, index) du segment courant differe de celle du
    segment precedent, quelle que soit son origine (meme couche, couche
    differente, passe differente).
    """
    raw_points: list[_RawPoint] = []
    previous_key: Optional[tuple[int, int]] = None
    normal_cache: dict[int, np.ndarray] = {}

    def _segment_polyline(layer_id: int, segment_index: int) -> tuple[Layer, np.ndarray]:
        layer = layers_by_id.get(layer_id)
        if layer is None:
            raise ToolpathResolutionError(
                f"PrintSchedule reference la couche {layer_id}, absente de `layers`."
            )
        toolpath = layer.toolpath_2d
        if toolpath is None:
            raise ToolpathResolutionError(
                f"La couche {layer_id} n'a pas de toolpath_2d - "
                f"pipeline/step7_scheduling.py doit etre execute avant l'etape 8."
            )
        if segment_index >= len(toolpath):
            raise ToolpathResolutionError(
                f"PrintSchedule reference le segment {segment_index} de la couche "
                f"{layer_id}, hors limites ({len(toolpath)} segments disponibles)."
            )
        return layer, toolpath[segment_index]

    for pass_ in schedule.passes:
        for layer_id, segment_index in pass_.segments:
            layer, polyline_2d = _segment_polyline(layer_id, segment_index)
            if len(polyline_2d) == 0:
                continue

            points_3d = _lift_to_3d(polyline_2d, layer.local_frame)
            # une normale PAR COUCHE, pas par polyligne : 29 691 np.cross sur
            # death @ 0,8 (0,84 s) pour 257 valeurs distinctes
            desired_direction = normal_cache.get(layer_id)
            if desired_direction is None:
                desired_direction = _layer_normal(layer.local_frame)
                normal_cache[layer_id] = desired_direction
            # Le sens de depot est donne par l'ORDRE de parcours des
            # isosurfaces, pas par le signe de grad(phi). Une couche atteinte
            # en DESCENDANT le champ (voir Layer.deposition_reversed) se
            # depose sous la matiere existante : la tete doit pointer a
            # l'oppose de la normale. Tant que l'ordre etait "phi croissant",
            # ce cas n'existait pas et la supposition +grad(phi) restait vraie
            # par construction.
            if getattr(layer, "deposition_reversed", False):
                desired_direction = -desired_direction
            key = (layer_id, segment_index)

            if previous_key is not None and key != previous_key:
                raw_points.append(
                    _RawPoint(
                        position=points_3d[0],
                        desired_direction=desired_direction,
                        layer_id=layer_id,
                        segment_index=segment_index,
                        color_id=pass_.color_id,
                        feedrate=(
                            params.travel_feedrate
                            if params.travel_feedrate > 0
                            else params.default_feedrate
                        ),
                        extrusion=0.0,
                        is_travel=True,
                    )
                )

            for point_3d in points_3d:
                raw_points.append(
                    _RawPoint(
                        position=point_3d,
                        desired_direction=desired_direction,
                        layer_id=layer_id,
                        segment_index=segment_index,
                        color_id=pass_.color_id,
                        feedrate=params.default_feedrate,
                        extrusion=1.0,
                        is_travel=False,
                    )
                )

            previous_key = key

    return raw_points



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
    positions = np.array([p.position for p in raw_points], dtype=float)
    local = direction_field(positions)
    for p, d in zip(raw_points, local):
        if not np.all(np.isfinite(d)):
            continue
        layer = layers_by_id[p.layer_id]
        if getattr(layer, "deposition_reversed", False):
            d = -d
        if float(np.dot(d, p.desired_direction)) <= 0.0:
            continue
        p.desired_direction = d

# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    layers: list[Layer],
    schedule: PrintSchedule,
    machine_profile: MachineProfile,
    params: Step8Params,
    direction_field=None,
) -> tuple[ToolpathMachine, list[Step8Diagnostic]]:
    """
    Resout la trajectoire machine reelle a partir du planning d'impression.

    `direction_field` (optionnel, utilise si `params.per_point_normal`) :
    appelable `(M,3) positions -> (M,3)` directions unitaires +grad(phi),
    lignes NaN la ou le champ est inexploitable (repli sur la normale de
    couche). Construit par l'orchestrateur (`gradient_direction_field`).

    Args:
        layers: sortie de pipeline/step7_scheduling.py (`colored_layers`) -
            PAS les couches d'entree de step7 : `PrintPass.segments` pointe
            vers ces couches-la (toolpath_2d re-decoupe par couleur).
        schedule: sortie de pipeline/step7_scheduling.py. `schedule.passes`
            est consomme dans l'ordre exact, jamais re-derive.
        machine_profile: cinematique de la machine (rotation_axes generique,
            longueur arbitraire) et origine du plateau.
        params: parametres de resolution IK, de singularite et de lissage.

    Returns:
        (ToolpathMachine, diagnostics) - diagnostics est une liste plate de
        Step8Diagnostic, une entree par point dont l'orientation desiree
        n'a pas pu etre atteinte a moins de
        `params.max_unreachable_angular_error_deg` pres. Jamais
        d'exception pour ce cas : c'est une donnee pure vers l'etape 10.

    Raises:
        ToolpathResolutionError: si le schedule reference une couche ou un
            segment absent/non rempli (donnee structurellement invalide,
            distinct d'une orientation inatteignable).
    """
    layers_by_id = {layer.id: layer for layer in layers}
    raw_points = _iter_schedule_points(schedule, layers_by_id, params)
    if params.per_point_normal and direction_field is not None and raw_points:
        _apply_per_point_normals(raw_points, layers_by_id, direction_field)

    axes = machine_profile.rotation_axes or []
    # Direction a angles nuls = SENS D'IMPRESSION de la tete au repos : tete
    # vers le bas, donc impression vers +Z (axe du corps de la tete, +Z du
    # modele de tete, cf. CLAUDE.md). Meme convention que step3/step7
    # (composed_rotation_directions(..., reference_direction=up)), puisque la
    # cible passee a l'IK est la normale de couche, c'est-a-dire le sens
    # d'impression.
    # Correction (2026-09-27) : la reference valait -Z, si bien qu'une couche
    # plate (cible +Z) donnait B = 179,9 deg -- tete retournee, corps SOUS la
    # buse -- et toute orientation etait inversee (couche inclinee de 11,5 deg
    # -> B = 168,5 deg). Aucun G-code n'etant encore produit, seule l'etape 9
    # en voyait les effets.
    reference_direction = np.array([0.0, 0.0, 1.0])
    platform_origin = machine_profile.platform_origin

    if getattr(machine_profile, "position_coupled_axes", None) and axes:
        # limite #59 : premier axe impose par la position (machine R C Z B)
        resolved_angles, residual_errors = _resolve_all_axis_angles_coupled(
            raw_points, axes, reference_direction, machine_profile, params
        )
    else:
        resolved_angles, residual_errors = _resolve_all_axis_angles(
            raw_points, axes, reference_direction, platform_origin, params
        )

    if params.free_axis_median_window > 1 and axes:
        resolved_angles, residual_errors = _median_filter_free_axes(
            resolved_angles, raw_points, axes, machine_profile, reference_direction,
            params.free_axis_median_window,
        )

    if params.clamp_angle_rates:
        smoothed_angles = _smooth_all_axes(resolved_angles, raw_points, axes, machine_profile)
    else:
        # axes coordonnes par le firmware : angles exacts (voir clamp_angle_rates)
        smoothed_angles = [list(a) for a in resolved_angles]

    points: list[ToolpathPoint] = []
    diagnostics: list[Step8Diagnostic] = []

    points = [
        ToolpathPoint(
            position=raw_point.position,
            axis_angles=smoothed_angles[i],
            extrusion=raw_point.extrusion,
            feedrate=raw_point.feedrate,
            color_id=raw_point.color_id,
        )
        for i, raw_point in enumerate(raw_points)
    ]
    # conversion en degres VECTORISEE (meme operation element par element que
    # np.degrees scalaire) : 364 192 appels numpy evites sur death @ 0,8
    residual_deg = np.degrees(np.asarray(residual_errors, dtype=float))
    for i in np.nonzero(residual_deg > params.max_unreachable_angular_error_deg)[0]:
        raw_point = raw_points[i]
        diagnostics.append(
            Step8Diagnostic(
                point_index=int(i),
                layer_id=raw_point.layer_id,
                segment_index=raw_point.segment_index,
                residual_angle_error_deg=float(residual_deg[i]),
                position=raw_point.position,
            )
        )

    return ToolpathMachine(points=points), diagnostics


#: Issue du dernier appel a `_resolve_all_axis_angles` : `None` = chemin
#: Python (interrupteur eteint ou aucun point), sinon NOMBRE de points ayant
#: exige le repli grille cote Python (0 = natif de bout en bout).
#: Existe pour que les tests d'equivalence puissent VERIFIER que le natif a
#: bien produit le resultat - sans quoi un repli silencieux rendrait le test
#: vert en ne comparant que Python a lui-meme (trou reellement rencontre).
_LAST_NATIVE_FALLBACK_COUNT: Optional[int] = None


def _resolve_all_axis_angles_native(
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    platform_origin: Optional[np.ndarray],
    params: Step8Params,
) -> tuple[list[list[float]], list[float]]:
    """
    Chemin natif, avec REPRISE apres chaque repli.

    Le natif resout analytiquement et s'arrete au premier point que la
    resolution de la paire externe ne sait pas traiter (cible colineaire a
    l'axe externe : le residu du bracketing y est une tangence a zero, jamais
    un changement de signe - typiquement un depot vertical). Ce point-la est
    resolu par la recherche par grille Python, deja ecrite et testee, puis le
    natif reprend au suivant.

    **Reprendre au lieu d'abandonner est l'essentiel de ce mecanisme.** Une
    premiere version rejouait TOUTE la boucle en Python des le premier repli :
    le premier point d'une trajectoire a fond plat depose verticalement, donc
    UN point sur 297 994 suffisait a annuler l'integralite du gain (mesure :
    89,5 s contre 89,7 s, soit rien). La chaine `previous_angles` est passee
    explicitement d'un segment au suivant, donc l'equivalence est preservee.
    """
    global _LAST_NATIVE_FALLBACK_COUNT

    native = _cpp.module()
    directions = np.array([p.desired_direction for p in raw_points], dtype=float)
    positions = np.array([p.position for p in raw_points], dtype=float)
    axes_spec = [(np.asarray(spec.axis, dtype=float), *spec.angle_range_deg) for spec in axes]
    origin = None if platform_origin is None else np.asarray(platform_origin, dtype=float)

    resolved_angles: list[list[float]] = []
    residual_errors: list[float] = []
    previous_angles: Optional[list[float]] = None
    cursor = 0
    fallback_count = 0

    while cursor < len(raw_points):
        angles_matrix, residuals, fallback_index = native.resolve_axis_angles_batch(
            directions[cursor:],
            positions[cursor:],
            axes_spec,
            np.asarray(reference_direction, dtype=float),
            origin,
            params.singularity_radius_epsilon,
            None if previous_angles is None else np.asarray(previous_angles, dtype=float),
        )

        # Sur repli, les lignes AVANT `fallback_index` sont valides ; les
        # suivantes n'ont jamais ete calculees et ne doivent pas etre lues.
        resolved_count = len(angles_matrix) if fallback_index < 0 else int(fallback_index)
        for k in range(resolved_count):
            resolved_angles.append(list(angles_matrix[k]))
            residual_errors.append(float(residuals[k]))
        if resolved_count:
            previous_angles = resolved_angles[-1]

        if fallback_index < 0:
            break

        fallback_count += 1
        raw_point = raw_points[cursor + resolved_count]
        angles, residual_error_rad = _resolve_axis_angles(
            raw_point.desired_direction, axes, reference_direction, params, previous_angles
        )
        singular_axes = _detect_singular_axes(
            raw_point.position, platform_origin, axes, angles, params.singularity_radius_epsilon
        )
        angles = _inherit_angle_at_singularity(singular_axes, angles, previous_angles)

        resolved_angles.append(angles)
        residual_errors.append(residual_error_rad)
        previous_angles = angles
        cursor += resolved_count + 1

    _LAST_NATIVE_FALLBACK_COUNT = fallback_count
    return resolved_angles, residual_errors


def coupled_axis_angle(
    position: np.ndarray,
    platform_origin: np.ndarray,
    axis: np.ndarray,
    radial_direction: np.ndarray,
    singularity_radius: float,
) -> Optional[float]:
    """
    Angle (degres, [0, 360)) de l'axe COUPLE A LA POSITION (limite #59) : la
    rotation autour de `axis` qui oriente `radial_direction` (+R a angles
    nuls) vers le point. La tete ne se deplacant que dans ce demi-plan, c'est
    le seul angle qui permet d'atteindre le point. `None` sur l'axe lui-meme
    (distance < singularity_radius) : l'angle y est indetermine (limite #19).
    """
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    rel = np.asarray(position, dtype=float) - np.asarray(platform_origin, dtype=float)
    rel = rel - axis * float(rel @ axis)
    if np.linalg.norm(rel) < singularity_radius:
        return None
    e1 = np.asarray(radial_direction, dtype=float)
    e1 = e1 - axis * float(e1 @ axis)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    return float(np.degrees(np.arctan2(rel @ e2, rel @ e1)) % 360.0)


def _resolve_all_axis_angles_coupled(
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    machine_profile: MachineProfile,
    params: Step8Params,
) -> tuple[list[list[float]], list[float]]:
    """
    IK quand le premier axe est IMPOSE PAR LA POSITION (machine actuelle :
    R C Z B, pas d'axe Y ; limite #59). Pour chaque point : l'angle couple
    vient de la position (coupled_axis_angle), puis la direction cible est
    ramenee dans le repere de cet axe et les axes RESTANTS la resolvent (meme
    `_resolve_axis_angles` que le cas libre). L'erreur residuelle est celle de
    la chaine complete : sur la machine actuelle, B ne penche que dans le plan
    radial, donc une cible qui en sort garde une erreur egale a son angle au
    plan -- elle devient un Step8Diagnostic, jamais une exception.

    Sur l'axe (singularite), l'angle couple est herite du point precedent
    (0 au premier point), comme le cas libre (limite #19).

    Chemin Python uniquement, dans les deux positions de APEX_USE_CPP : le
    natif ne connait pas ce couplage (resultats identiques par construction).
    """
    coupled = list(machine_profile.position_coupled_axes)
    if coupled != [0]:
        raise ValueError(
            f"position_coupled_axes={coupled} : seul le premier axe peut etre couple a la position"
        )
    if machine_profile.radial_direction is None or machine_profile.platform_origin is None:
        raise ValueError("position_coupled_axes exige radial_direction et platform_origin")
    outer = axes[0]
    outer_axis = np.asarray(outer.axis, dtype=float)
    outer_axis = outer_axis / np.linalg.norm(outer_axis)
    inner_axes = axes[1:]

    if _cpp.enabled() and raw_points:
        return _resolve_all_axis_angles_coupled_native(
            raw_points, axes, outer_axis, reference_direction, machine_profile, params
        )

    resolved_angles: list[list[float]] = []
    residual_errors: list[float] = []
    previous: Optional[list[float]] = None
    for raw_point in raw_points:
        c = coupled_axis_angle(
            raw_point.position, machine_profile.platform_origin, outer_axis,
            machine_profile.radial_direction, params.singularity_radius_epsilon,
        )
        if c is None:
            c = previous[0] if previous is not None else 0.0
        target = np.asarray(raw_point.desired_direction, dtype=float)
        target = target / max(np.linalg.norm(target), 1e-12)
        local_target = _rodrigues_rotate(target[None, :], outer_axis, np.array([np.radians(-c)]))[0]
        inner, _res = _resolve_axis_angles(
            local_target, inner_axes, reference_direction, params,
            None if previous is None else previous[1:],
        )
        angles = [c] + list(inner)
        achieved = reference_direction
        for spec, angle_deg in zip(reversed(axes), reversed(angles)):
            ax = np.asarray(spec.axis, dtype=float)
            achieved = _rodrigues_rotate(achieved[None, :], ax / np.linalg.norm(ax),
                                         np.array([np.radians(angle_deg)]))[0]
        residual = float(np.arccos(np.clip(achieved @ target, -1.0, 1.0)))
        resolved_angles.append(angles)
        residual_errors.append(residual)
        previous = angles
    return resolved_angles, residual_errors



def _resolve_all_axis_angles_coupled_native(
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    outer_axis: np.ndarray,
    reference_direction: np.ndarray,
    machine_profile: MachineProfile,
    params: Step8Params,
) -> tuple[list[list[float]], list[float]]:
    """
    Pas 17 du portage C++ (CLAUDE.md limite #56) : meme calcul que la boucle
    Python de `_resolve_all_axis_angles_coupled`, reorganise en trois temps.

    1. Angle couple C de TOUS les points d'un coup (vectorise). Sa seule
       dependance sequentielle est l'heritage sur l'axe (singularite) : un
       remplissage vers l'avant, 0 au debut.
    2. Cibles ramenees dans le repere de C (Rodrigues par lots).
    3. Axes internes resolus par l'IK NATIVE du pas 2, SANS detection de
       singularite (la boucle Python n'en fait pas non plus pour les axes
       internes), avec la meme reprise au point fautif : la grille Python ne
       traite que ce point.

    Equivalence : TOLERANCE declaree, comme au pas 2 (angles a ~1e-9 deg).
    Les produits scalaires vectorises peuvent differer d'un ULP des produits
    scalaires point par point de la boucle Python.
    """
    global _LAST_NATIVE_FALLBACK_COUNT
    native = _cpp.module()
    positions = np.array([p.position for p in raw_points], dtype=float)
    targets = np.array([p.desired_direction for p in raw_points], dtype=float)
    targets = targets / np.maximum(np.linalg.norm(targets, axis=1), 1e-12)[:, None]

    rel = positions - np.asarray(machine_profile.platform_origin, dtype=float)
    rel = rel - np.outer(rel @ outer_axis, outer_axis)
    on_axis = np.linalg.norm(rel, axis=1) < params.singularity_radius_epsilon
    e1 = np.asarray(machine_profile.radial_direction, dtype=float)
    e1 = e1 - outer_axis * float(e1 @ outer_axis)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(outer_axis, e1)
    c = np.degrees(np.arctan2(rel @ e2, rel @ e1)) % 360.0
    if on_axis.any():
        c = c.copy()
        last = 0.0
        for i in range(len(c)):
            if on_axis[i]:
                c[i] = last
            last = c[i]

    local_targets = _rodrigues_rotate(targets, outer_axis, np.radians(-c))
    inner_axes = axes[1:]
    inner_spec = [(np.asarray(spec.axis, dtype=float), *spec.angle_range_deg) for spec in inner_axes]

    inner_rows: list[list[float]] = []
    previous: Optional[list[float]] = None
    cursor = 0
    fallback_count = 0
    n = len(raw_points)
    while cursor < n:
        angles_matrix, _res, fallback_index = native.resolve_axis_angles_batch(
            local_targets[cursor:], positions[cursor:], inner_spec,
            np.asarray(reference_direction, dtype=float), None,
            params.singularity_radius_epsilon,
            None if previous is None else np.asarray(previous, dtype=float),
        )
        count = len(angles_matrix) if fallback_index < 0 else int(fallback_index)
        inner_rows.extend(list(row) for row in angles_matrix[:count])
        if count:
            previous = inner_rows[-1]
        if fallback_index < 0:
            break
        fallback_count += 1
        k = cursor + count
        inner, _r = _resolve_axis_angles(local_targets[k], inner_axes, reference_direction, params, previous)
        inner_rows.append(list(inner))
        previous = inner_rows[-1]
        cursor = k + 1
    _LAST_NATIVE_FALLBACK_COUNT = fallback_count

    angles = np.column_stack([c, np.array(inner_rows, dtype=float).reshape(n, len(inner_axes))])
    achieved = np.tile(np.asarray(reference_direction, dtype=float), (n, 1))
    for k in range(len(axes) - 1, -1, -1):
        ax = np.asarray(axes[k].axis, dtype=float)
        achieved = _rodrigues_rotate(achieved, ax / np.linalg.norm(ax), np.radians(angles[:, k]))
    residuals = np.arccos(np.clip(np.einsum("ij,ij->i", achieved, targets), -1.0, 1.0))
    return [list(map(float, row)) for row in angles], [float(r) for r in residuals]

def _resolve_all_axis_angles(
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    platform_origin: Optional[np.ndarray],
    params: Step8Params,
) -> tuple[list[list[float]], list[float]]:
    """
    Resout l'angle de chaque axe pour TOUS les points de la trajectoire.

    La boucle porte une dependance SEQUENTIELLE reelle : `previous_angles`
    sert de valeur de repliement pour les axes internes (`axes[2:]`), de
    repli en cas de singularite locale, et d'angle herite quand un axe est
    detecte singulier. Elle n'est donc PAS parallelisable telle quelle -
    contrairement a ce que laissait croire son allure de boucle par point.

    Chemin natif (portage C++, voir cpp/src/core/kinematics.cpp) quand
    l'interrupteur est actif ; sinon la boucle Python historique, conservee
    intacte juste en dessous.

    Le natif ne porte PAS la recherche par grille de repli
    (`_resolve_axis_angles_via_grid_search`) : quand la resolution analytique
    de la paire externe echoue, il rend la main en signalant l'index fautif,
    et TOUTE la boucle est rejouee en Python. Rejouer tout, et non le seul
    point concerne, est necessaire a l'equivalence : reprendre au milieu avec
    une chaine `previous_angles` differente divergerait silencieusement.
    Ce cas n'est jamais atteint sur la machine a 2 axes orthogonaux actuelle.
    """
    global _LAST_NATIVE_FALLBACK_COUNT
    _LAST_NATIVE_FALLBACK_COUNT = None

    if _cpp.enabled() and raw_points:
        return _resolve_all_axis_angles_native(
            raw_points, axes, reference_direction, platform_origin, params
        )

    resolved_angles: list[list[float]] = []
    residual_errors: list[float] = []
    previous_angles: Optional[list[float]] = None

    for raw_point in raw_points:
        angles, residual_error_rad = _resolve_axis_angles(
            raw_point.desired_direction, axes, reference_direction, params, previous_angles
        )

        singular_axes = _detect_singular_axes(
            raw_point.position, platform_origin, axes, angles, params.singularity_radius_epsilon
        )
        angles = _inherit_angle_at_singularity(singular_axes, angles, previous_angles)

        resolved_angles.append(angles)
        residual_errors.append(residual_error_rad)
        previous_angles = angles

    return resolved_angles, residual_errors


def _median_filter_free_axes(resolved_angles, raw_points, axes, machine_profile, reference_direction, window):
    """Filtre median des axes non couples, par polyligne (voir
    Step8Params.free_axis_median_window). Une polyligne s'arrete a chaque
    trajet a vide et a chaque changement de (couche, segment). Rend (angles,
    ecarts residuels recalculés en radians)."""
    from scipy.ndimage import median_filter

    from core.geometry_utils import _rodrigues_rotate_single

    angles = np.array(resolved_angles, dtype=float)
    coupled = set(getattr(machine_profile, "position_coupled_axes", None) or [])
    free = [i for i in range(len(axes)) if i not in coupled]
    keys = [(p.layer_id, p.segment_index, p.is_travel) for p in raw_points]
    start = 0
    for i in range(1, len(raw_points) + 1):
        if i == len(raw_points) or keys[i] != keys[start] or raw_points[i].is_travel:
            if i - start > 2 and not raw_points[start].is_travel:
                for ax in free:
                    angles[start:i, ax] = median_filter(angles[start:i, ax], size=window, mode="nearest")
            start = i
    units = [np.asarray(s.axis, dtype=float) / np.linalg.norm(s.axis) for s in axes]
    residual = np.zeros(len(raw_points))
    for i, p in enumerate(raw_points):
        d = np.asarray(reference_direction, dtype=float)
        for u, a in zip(reversed(units), reversed(angles[i])):
            d = _rodrigues_rotate_single(d, u, np.radians(a))
        if p.desired_direction is None:
            continue
        t = np.asarray(p.desired_direction, dtype=float)
        t = t / max(np.linalg.norm(t), 1e-12)
        residual[i] = float(np.arccos(np.clip(d @ t, -1.0, 1.0)))
    return [list(a) for a in angles], list(residual)


def _smooth_all_axes(
    resolved_angles: list[list[float]],
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    machine_profile: MachineProfile,
) -> list[list[float]]:
    """
    Lisse la serie temporelle d'angles, axe par axe, sur l'ENSEMBLE des
    points (pas par passe : le lissage a besoin de la serie complete pour
    respecter les limites de vitesse/acceleration sur toute la trajectoire).
    Depli/repli (wraparound) pour les axes a tour complet, clamp de
    taux, puis reprojection - voir _unwrap_angle_sequence/_clamp_rate/
    _rewrap_angle_sequence.

    `dt_per_point[i]` derive du feedrate et de la distance parcourue entre
    les points i-1 et i, donnant un pas de temps physiquement significatif
    plutot qu'un taux d'echantillonnage suppose fixe.
    """
    n = len(raw_points)
    if n == 0 or not axes:
        return [list(a) for a in resolved_angles]

    dt_per_point = np.ones(n)
    for i in range(1, n):
        distance = float(np.linalg.norm(raw_points[i].position - raw_points[i - 1].position))
        feedrate = max(raw_points[i].feedrate, 1e-6)
        dt_per_point[i] = distance / feedrate

    angles_matrix = np.array(resolved_angles, dtype=float)  # (n, len(axes))
    smoothed_matrix = angles_matrix.copy()

    coupled = set(getattr(machine_profile, "position_coupled_axes", None) or [])
    for axis_index, spec in enumerate(axes):
        if axis_index in coupled:
            # angle impose par la position : le lisser deplacerait le point
            # atteint (limite #59) ; c'est la vitesse d'avance qui devra
            # respecter la vitesse du plateau (a traiter a l'export G-code)
            continue
        is_full_turn = _is_full_turn_axis(spec)
        column = angles_matrix[:, axis_index]

        unwrapped = _unwrap_angle_sequence(column, is_full_turn)
        clamped = _clamp_rate(
            unwrapped,
            dt_per_point,
            machine_profile.max_angular_velocity,
            machine_profile.max_angular_acceleration,
        )
        rewrapped = _rewrap_angle_sequence(clamped, spec.angle_range_deg, is_full_turn)
        smoothed_matrix[:, axis_index] = rewrapped

    return [list(row) for row in smoothed_matrix]
