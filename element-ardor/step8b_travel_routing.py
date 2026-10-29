"""
pipeline/step8b_travel_routing.py

Etape intermediaire entre step8 (trajectoire machine) et step9 (simulation
de collision) : route les trajets a vide LONGS en lift-and-traverse
(lever en Z, traverser a cette hauteur, redescendre) plutot que la ligne
droite brute inseree par `step8._iter_schedule_points` - voir CLAUDE.md
limite #46 pour le diagnostic complet qui a motive ce fix (trajets a vide
jusqu'a 44.5mm en ligne droite a travers le vide, mesures sur death.stl).

*** CE N'EST PAS UN VRAI PATH-PLANNING 3D. *** C'est un palliatif
deliberement simple (3 segments : lift vertical, traverse horizontale a
hauteur constante, descente verticale), choisi avec l'utilisateur comme
portee raisonnable pour cette session plutot qu'une recherche de chemin
complete (A*/RRT) dans l'espace libre - voir CLAUDE.md "limite du slicer a
corriger plus tard" pour la portee explicitement non couverte (une piece
tres tourmentee dont AUCUNE hauteur de lift ne degage la tete reste hors de
portee de ce mecanisme; TravelRoutingDiagnostic signale ce cas plutot que
de le masquer).

Ce module a BESOIN de la trajectoire ENTIERE (ToolpathMachine, sortie de
step8) et du modele de tete (MachineProfile.nozzle_model) pour valider un
lift par vraie detection de collision - contrairement a step8, qui produit
la trajectoire point par point sans connaitre l'etat de la matiere deja
imprimee en avance. D'ou sa position dans le pipeline : APRES step8, AVANT
step9 (qui simule la trajectoire deja routee).

Reutilise les PRIMITIVES de step9 (_PrintedGeometryIndex,
_grow_printed_index_to, _evaluate_point, _surface_sample_points) - jamais
une seconde implementation de la detection de collision, meme discipline
que step10 vis-a-vis de step9.

Sortie : (ToolpathMachine routee, list[TravelRoutingDiagnostic]).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np
import trimesh

from core.mesh_types import (
    MachineProfile,
    ToolpathMachine,
    ToolpathPoint,
    TravelRoutingDiagnostic,
)
from pipeline.step9_collision_sim import (
    CollisionSimParams,
    _aabb_overlaps,
    _grow_printed_index_to,
    _points_to_segment_distances,
    _PrintedGeometryIndex,
    _positioned_toolhead_vertices,
    _surface_sample_points,
    _toolhead_aabb,
)


@dataclass
class TravelRoutingParams:
    # En-deca de cette longueur (mm), un trajet a vide reste une ligne
    # droite (comportement historique inchange) - router un trajet court
    # n'a pas de benefice mesurable et ajoute un cout de calcul inutile.
    # Ordre de grandeur choisi au-dessus de l'espacement de point typique
    # en extrusion (~0.5mm mediane mesuree sur death.stl, CLAUDE.md limite
    # #37), pour ne cibler que les VRAIS trajets a vide de fin de polyligne,
    # pas un bruit d'echantillonnage local.
    length_threshold: float = 5.0
    # Hauteur de lift initiale testee (mm au-dessus du plus haut des deux
    # points start/end du trajet).
    lift_height: float = 3.0
    # Multiplicateur de plafond : si `lift_height` collisionne encore, la
    # hauteur est doublee successivement jusqu'a
    # `lift_height * max_lift_multiplier` avant d'abandonner et d'emettre
    # un TravelRoutingDiagnostic.
    max_lift_multiplier: float = 4.0
    # Nombre de points intermediaires generes le long de la traverse
    # horizontale (en plus des deux coudes lift-depart/lift-arrivee) - une
    # traverse longue merite plus d'un seul point de test de collision.
    # Defaut a 1 (un seul point median) plutot que plusieurs : chaque point
    # supplementaire multiplie directement le cout de validation (voir
    # routing_toolhead_sample_spacing ci-dessous) - un seul point median
    # suffit a detecter la grande majorite des collisions de traverse sur
    # les modeles mesures, le reste etant couvert par les deux coudes.
    traverse_samples: int = 1
    # Pas d'echantillonnage des ARETES du modele de tete, UNIQUEMENT pour la
    # validation de collision du routage (distinct de
    # CollisionSimParams.toolhead_sample_spacing, utilise par step9/step10
    # pour la simulation complete et precise). Un espacement plus grossier
    # ici est un compromis delibere et mesure : le routage teste bien plus
    # de candidats que step9 (jusqu'a plusieurs niveaux de lift x plusieurs
    # points de chemin x tous les trajets a vide longs), donc le cout de
    # `_positioned_toolhead_vertices`/`_check_printed_part` par appel doit
    # rester borne - voir CLAUDE.md limite #46 pour la mesure qui a motive
    # ce parametre (127 trajets longs sur l-triple-axis-overhanged-shape.stl,
    # ~1500 points d'echantillon a spacing=5.0 rendaient le routage complet
    # de l'ordre de plusieurs minutes, inutilisable en pratique).
    routing_toolhead_sample_spacing: float = 15.0
    # Multiplicateur de taille de cellule de la grille spatiale de
    # _PrintedGeometryIndex, UNIQUEMENT pour le routage (distinct de
    # CollisionSimParams.grid_cell_size_multiplier, utilise par step9/step10
    # pour des requetes ponctuelles autour d'UN SEUL point de trajectoire).
    # Une requete de routage porte sur l'AABB d'un CHEMIN ENTIER (plusieurs
    # dizaines de mm de large avec le lift), pas un point : avec le
    # multiplicateur par defaut de step9 (cellules de nozzle_diameter*4 =
    # 1.6mm), une telle requete dilate sur des centaines de cellules et
    # renvoie la quasi-totalite de l'index accumule comme candidats - mesure
    # concretement (script de diagnostic jetable, non conserve) sur
    # l-triple-axis-overhanged-shape.stl : le temps de routage passait de
    # 3s a 25s entre le 20e et le 80e trajet route, correle directement a la
    # taille de l'index (5000 -> 12000 capsules), pas au nombre de trajets
    # lui-meme - signature classique d'une grille trop fine pour la taille
    # de requete. Des cellules bien plus grandes (echelle de la TETE, pas du
    # CORDON - meme distinction que toolhead_sample_spacing vs
    # nozzle_diameter, limite #26) reduisent le nombre de cellules a
    # parcourir par requete sans changer le resultat (toujours un
    # sur-ensemble exact des vrais candidats, jamais un faux negatif).
    routing_grid_cell_size_multiplier: float = 40.0


def _travel_boundaries(points: list[ToolpathPoint]) -> list[int]:
    """
    Indices `i` tels que `points[i]` est le point de DESTINATION d'un trajet
    a vide inséré par step8 - reconnu par `extrusion == 0.0` sur `points[i]`
    ET sur `points[i-1]` n'etant PAS lui-meme ce meme point de trajet
    (`_iter_schedule_points` insere exactement UN point extrusion=0 par
    frontiere, jamais deux consecutifs) : concretement, `points[i]` est un
    point de trajet a vide si son extrusion est nulle - simple, coherent
    avec le contrat deja etabli par step8 (`is_travel`/`extrusion=0.0` non
    expose sur `ToolpathPoint`, seul `extrusion` l'est).
    """
    return [i for i in range(1, len(points)) if points[i].extrusion <= 0.0]


def _make_intermediate_point(
    position: np.ndarray, template: ToolpathPoint, feedrate: float
) -> ToolpathPoint:
    return replace(template, position=position, extrusion=0.0, feedrate=feedrate)


def _lifted_path_points(
    start: np.ndarray, end: np.ndarray, lift: float, traverse_samples: int
) -> list[np.ndarray]:
    """
    Positions intermediaires du chemin lift-and-traverse : monte
    verticalement de `lift` au-dessus de `start`, traverse horizontalement
    (echantillonnee) jusqu'a l'aplomb de `end`, redescend verticalement.
    Le point `start` lui-meme n'est PAS inclus (deja present dans la
    trajectoire) ; le point `end` final est ajoute par l'appelant (c'est le
    point de destination existant, jamais duplique ici).
    """
    up = np.array([0.0, 0.0, lift])
    lifted_start = start + up
    lifted_end = end + up

    traverse_points = []
    n = max(traverse_samples, 0)
    for k in range(1, n + 1):
        t = k / (n + 1)
        traverse_points.append(lifted_start + t * (lifted_end - lifted_start))

    return [lifted_start, *traverse_points, lifted_end]


def _path_collides(
    path_positions: list[np.ndarray],
    toolhead_base_samples: Optional[np.ndarray],
    rotation_axes,
    reference_angles: list[float],
    printed_index: Optional[_PrintedGeometryIndex],
    safety_margin: float,
) -> bool:
    """
    True si au moins un point du chemin candidat penetre la matiere deja
    imprimee.

    Positionne la tete a CHAQUE point du chemin candidat (comme step9),
    mais teste l'UNION de tous ces points positionnes contre l'index en UN
    SEUL passage de broad-phase + boucle de candidats, au lieu d'un appel
    complet a `_check_printed_part` par point. Necessaire en pratique : le
    routage teste bien plus de candidats que step9 (jusqu'a plusieurs
    niveaux de lift x plusieurs points de chemin x tous les trajets a vide
    longs d'une trajectoire), et repeter la boucle de candidats (avec sa
    propre requete `candidates_near`) par point individuel domine le temps
    total des lors que l'index de matiere imprimee est deja volumineux (meme
    pathologie que CLAUDE.md limites #27/#30, deja mesuree et corrigee pour
    step9/step10 - meme principe applique ici : un seul passage sur un plus
    gros tableau plutot que N passages sur N petits tableaux).

    Sans modele de tete/index, aucun test n'est possible : retourne False
    (rien a valider, meme garde que step10._try_local_reorientation quand
    `has_toolhead` est faux).
    """
    if toolhead_base_samples is None or printed_index is None or not path_positions:
        return False

    all_samples = np.concatenate(
        [
            _positioned_toolhead_vertices(toolhead_base_samples, position, reference_angles, rotation_axes)
            for position in path_positions
        ],
        axis=0,
    )

    aabb = _toolhead_aabb(all_samples)
    candidates = printed_index.candidates_near(aabb, safety_margin)
    if not candidates:
        return False

    for segment in candidates:
        if not _aabb_overlaps(aabb, segment.aabb, safety_margin):
            continue
        distances = _points_to_segment_distances(all_samples, segment.start, segment.end)
        distances -= segment.radius
        if float(distances.min()) < safety_margin:
            return True
    return False


def route_long_travels(
    toolpath: ToolpathMachine,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    params: TravelRoutingParams,
) -> tuple[ToolpathMachine, list[TravelRoutingDiagnostic]]:
    """
    Remplace chaque trajet a vide de longueur > `params.length_threshold`
    par un chemin lift-and-traverse valide par detection de collision
    contre la matiere deja imprimee (index incrementalement grossi, meme
    invariant "jamais en avance sur le point teste" que step9/step10).

    Ne mute jamais `toolpath` - retourne une copie. Les points en extrusion
    (`extrusion > 0`) sont copies tels quels et strictement dans le meme
    ordre : seuls les points de trajet a vide (`extrusion <= 0`) peuvent
    voir leur position dediee remplacee par une sequence de plusieurs
    points.

    Sans modele de tete charge (`machine_profile.nozzle_model` absent),
    aucune validation de collision n'est possible : le lift standard
    (`params.lift_height`) est applique sans essai d'augmentation ni
    diagnostic - un lift non valide reste strictement meilleur qu'une ligne
    droite non validee (memes soulevements deja utilises par les slicers
    3-axes classiques comme palliatif), mais n'est jamais presente comme
    "verifie sans collision" (aucun TravelRoutingDiagnostic n'est
    d'ailleurs emis dans ce cas, faute de pouvoir juger).
    """
    nozzle_model = machine_profile.nozzle_model
    has_toolhead = isinstance(nozzle_model, trimesh.Trimesh) and len(nozzle_model.vertices) > 0

    toolhead_base_samples = None
    if has_toolhead:
        toolhead_base_samples = _surface_sample_points(nozzle_model, params.routing_toolhead_sample_spacing)

    printed_index: Optional[_PrintedGeometryIndex] = None
    if has_toolhead:
        printed_index = _PrintedGeometryIndex(
            machine_profile.nozzle_diameter, params.routing_grid_cell_size_multiplier
        )

    points = list(toolpath.points)
    boundary_indices = set(_travel_boundaries(points))

    routed_points: list[ToolpathPoint] = []
    diagnostics: list[TravelRoutingDiagnostic] = []
    last_grown_idx = 0

    for i, point in enumerate(points):
        if has_toolhead:
            last_grown_idx = _grow_printed_index_to(printed_index, points, i, last_grown_idx)

        if i not in boundary_indices or i == 0:
            routed_points.append(point)
            continue

        start = points[i - 1].position
        end = point.position
        travel_length = float(np.linalg.norm(end - start))

        if travel_length <= params.length_threshold:
            routed_points.append(point)
            continue

        candidate_lift = params.lift_height
        max_lift = params.lift_height * max(params.max_lift_multiplier, 1.0)
        chosen_path: Optional[list[np.ndarray]] = None
        last_tested_lift = candidate_lift

        while True:
            path = _lifted_path_points(start, end, candidate_lift, params.traverse_samples)
            last_tested_lift = candidate_lift
            if not _path_collides(
                path, toolhead_base_samples, machine_profile.rotation_axes,
                point.axis_angles, printed_index, sim_params.safety_margin,
            ):
                chosen_path = path
                break
            if candidate_lift >= max_lift or not has_toolhead:
                break
            candidate_lift = min(candidate_lift * 2.0, max_lift)

        if chosen_path is None:
            # Aucune hauteur testee n'a evite la collision (ou pas de
            # modele de tete pour valider) : on route quand meme via le
            # DERNIER candidat essaye plutot que retomber sur la ligne
            # droite (strictement pire, jamais un progres) - mais on le
            # signale explicitement si un modele de tete etait disponible
            # et confirme la collision residuelle.
            chosen_path = _lifted_path_points(start, end, last_tested_lift, params.traverse_samples)
            if has_toolhead:
                diagnostics.append(
                    TravelRoutingDiagnostic(
                        point_index=i,
                        travel_length=travel_length,
                        max_lift_tested=last_tested_lift,
                        still_colliding=True,
                    )
                )

        travel_feedrate = points[i - 1].feedrate
        for intermediate_position in chosen_path[:-1]:
            routed_points.append(_make_intermediate_point(intermediate_position, point, travel_feedrate))
        routed_points.append(point)

    return ToolpathMachine(points=routed_points), diagnostics
