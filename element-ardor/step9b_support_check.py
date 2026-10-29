"""
pipeline/step9b_support_check.py

Etape 9b du pipeline : verifie, sur la trajectoire machine REELLE, que chaque
point de depot a de la matiere DEJA IMPRIMEE en dessous de lui au moment ou il
est extrude.

Positionne entre step8/step8b et step9. Ne depend ni de step9 ni de step10 : il
consomme la meme `ToolpathMachine` que step9, sans avoir besoin du
`CollisionReport`. Le placer tot donne le diagnostic le plus proche de la
sortie brute de step5/6/7/8 - et il reste valide apres step10, qui ne modifie
jamais `position`, seulement `axis_angles`.

## Ce qu'il mesure, et pourquoi ce n'est pas step2

`step2_overhang` detecte les surplombs a partir de la seule GEOMETRIE du
modele : il ne connait pas l'ordre d'impression. Or une face parfaitement
supportee dans la piece finie peut tres bien etre deposee dans le vide si la
matiere qui la porte est imprimee APRES elle. C'est exactement le defaut
mesure a l'origine de cette etape sur res/test_models/death.stl : 11.25% de la
matiere extrudee deposee sans rien en dessous, dont 99% avait pourtant un
support - imprime trop tard.

Cette etape separe donc explicitement les deux natures de defaut
(`UnsupportedDepositEvent.has_future_support`) :
  - support existant mais depose plus tard  -> defaut d'ORDRE, corrigible par
    un reordonnancement, sans toucher a la geometrie ;
  - aucun support a aucun moment            -> surplomb geometrique reel,
    qu'aucun ordre ne peut resoudre.

Sans ce distinguo, un taux global de "matiere en l'air" ne dit pas s'il faut
reordonner ou ajouter du support : les deux corrections sont incomparables.

## Validite du test vertical - et sa limite corrigee

La classification "au moment T" (defaut d'ordre vs non-supporte) cherche le
support VERS LE BAS (`-up_vector`) : c'est une question sur la geometrie
MONDE deja deposee a cet instant, independante de l'orientation de la tete a
ce moment precis - un test vertical y est valide.

Mais la classification FINALE "surplomb reel" (aucun support a AUCUN moment,
meme dans la piece finie) ne peut PAS se limiter a la verticale : mesure
directement sur res/test_models/death.stl, les 194 points que ce test
classait "surplomb reel" en ne regardant QUE la verticale ont TOUS (194/194)
un support des qu'on teste le long d'une direction atteignable par la
cinematique reelle de la machine (`machine_profile.rotation_axes`, meme
fonction `composed_rotation_directions` que step3) - la machine peut
deposer ce point en inclinant la tete/le plateau, ce n'est pas un vrai
surplomb geometrique. Voir `UnsupportedDepositEvent.supported_by_tilt` et
`_has_support_in_any_direction` : la classification finale teste maintenant
TOUTES les directions atteignables avant de conclure a un surplomb reel.

## Reutilisation, jamais duplication

L'indexation de la matiere imprimee (`_PrintedGeometryIndex`,
`_grow_printed_index_to`) et la distance point-capsule
(`_point_to_segment_distance`) sont importees telles quelles depuis step9 -
meme patron que `step10_adjustment`, qui importe deja ces primitives. Aucune
seconde implementation de l'etat incremental de la matiere imprimee n'existe
dans le projet.

Sortie : SupportCheckReport
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from core.geometry_utils import composed_rotation_directions, point_has_support_along_direction
from core.mesh_types import (
    AABB,
    MachineProfile,
    SupportCheckReport,
    ToolpathMachine,
    ToolpathPoint,
    UnsupportedDepositEvent,
)
from pipeline.step3_segmentation import _samples_per_axis
from pipeline.step9_collision_sim import (
    _PrintedGeometryIndex,
    _grow_printed_index_to,
    _point_to_segment_distance,
    _sample_point_indices,
)

#: Verticale monde, sous la forme `(1,3)` attendue par `index.has_support`.
#: Constante de module plutot qu'un `np.array` reconstruit a chaque appel :
#: `_has_support_below` est appelee 273 146 fois sur death.stl.
_WORLD_UP_SINGLETON = np.array([[0.0, 0.0, 1.0]])


@dataclass
class SupportCheckParams:
    """
    Les rayons/hauteurs de recherche sont exprimes en MULTIPLES de grandeurs
    physiques reelles (diametre de buse, hauteur de couche), jamais en
    constantes absolues : une constante universelle serait soit aberrante sur
    un modele de 5 mm, soit aveugle sur un modele de 500 mm. Meme discipline
    que `bridge_max_gap_distance` (step2), calibre relativement a
    `layer_height_base`.
    """

    # Rayon de recherche horizontal sous le point teste, en multiples de
    # `nozzle_diameter`. Le cordon a lui-meme cette largeur : un support
    # valide doit etre dans son emprise, pas a un rayon arbitraire.
    # 2.5 x 0.4 mm (defaut projet) = 1.0 mm.
    horizontal_search_radius_factor: float = 2.5

    # Hauteur de recherche verticale sous le point, en multiples de
    # `layer_height_base`. Doit couvrir plusieurs couches : sur une couche
    # courbe, la matiere porteuse n'est pas exactement une hauteur de couche
    # en dessous. 7.5 x 0.2 mm (defaut projet) = 1.5 mm.
    vertical_search_height_factor: float = 7.5

    # Tolerance "pose sur le plateau", en multiples de `nozzle_diameter` : un
    # point assez proche du plateau n'a structurellement besoin d'aucun
    # support en dessous. 2.5 x 0.4 mm = 1.0 mm.
    platform_contact_tolerance_factor: float = 2.5

    # Pas d'echantillonnage le long de la trajectoire (distance cumulee, meme
    # convention que `CollisionSimParams.path_sampling_step`).
    # 0.0 = teste TOUS les points en extrusion.
    #
    # Contrairement a step9 (ou sous-echantillonner est acceptable : la tete
    # a une emprise large, une collision manquee a un point est rattrapee au
    # suivant), ici c'est le POINT DE DEPOT lui-meme qui est teste. Sauter des
    # points sous-echantillonne DIRECTEMENT la mesure. Le defaut sur teste donc
    # tout ; un sous-echantillonnage reste possible, a la charge explicite de
    # l'appelant qui accepte alors une mesure approximative.
    path_sampling_step: float = 0.0

    # Localiser, pour chaque defaut d'ORDRE, le premier point de trajectoire
    # qui l'aurait supporte (`UnsupportedDepositEvent.supporting_point_index`).
    #
    # Desactive par defaut : cette localisation exige de rejouer la croissance
    # incrementale de la matiere, soit O(defauts x segments), la ou la
    # distinction ordre/surplomb - la donnee qui decide QUOI corriger - coute
    # O(defauts) et reste TOUJOURS calculee. Mesure sur death.stl (pire cas,
    # 22 393 evenements) : 9.7 s sans localisation, contre plusieurs minutes
    # avec. A activer pour une investigation ciblee, pas pour une mesure de
    # routine.
    locate_supporting_point: bool = False

    # Taille de cellule de la grille spatiale, en multiples de
    # `nozzle_diameter`. DEDIE, et volontairement 4x plus fin que le defaut de
    # step9 (4.0) : une requete de support porte sur une boite de ~2x2x1.5 mm
    # sous un point, la ou step9 interroge l'AABB d'une tete de 76 mm. Une
    # grille dimensionnee pour la tete retourne massivement trop de candidats
    # pour cet usage - mesure sur death.stl : 360 candidats en mediane a 4.0
    # contre 84 a 1.0, pour un cout de construction et de requete inchange.
    # Meme nature d'ajustement que `routing_grid_cell_size_multiplier` de
    # step8b, dedie pour la raison symetrique (requetes bien plus LARGES).
    grid_cell_size_multiplier: float = 1.0

    # Nombre de directions candidates testees pour la classification finale
    # "surplomb reel" (voir docstring de module, "Validite du test
    # vertical"). Meme defaut ET meme fonction que
    # `SegmentationParams.direction_candidate_count` (step3) - pas une
    # nouvelle densite inventee : mesure sur death.stl, cette densite
    # resout deja 194/194 des "surplombs reels" precedemment mal classes,
    # en 0.26s pour les 194 points concernes (la passe ne s'applique
    # qu'aux points deja identifies comme surplombs candidats, jamais a
    # toute la trajectoire).
    direction_candidate_count: int = 26


def _search_aabb_below(
    position: np.ndarray, horizontal_radius: float, vertical_height: float
) -> AABB:
    """
    Boite de recherche du support : etendue de `horizontal_radius` en X/Y, et
    de `vertical_height` VERS LE BAS uniquement (jamais vers le haut - une
    capsule au-dessus du point ne le supporte pas).

    Le sommet de la boite est le point lui-meme : une capsule strictement
    au-dessus ne peut jamais entrer dans les candidats.
    """
    return AABB(
        min_corner=np.array(
            [
                position[0] - horizontal_radius,
                position[1] - horizontal_radius,
                position[2] - vertical_height,
            ]
        ),
        max_corner=np.array(
            [position[0] + horizontal_radius, position[1] + horizontal_radius, position[2]]
        ),
    )


def _has_support_below(
    index: _PrintedGeometryIndex,
    position: np.ndarray,
    horizontal_radius: float,
    vertical_height: float,
) -> bool:
    """
    Existe-t-il, dans `index`, une capsule qui supporte `position` ?

    Un candidat supporte le point si son point le plus proche est
    STRICTEMENT en dessous (`dz > 0`), a moins de `vertical_height` en
    dessous, et a moins de `horizontal_radius` horizontalement.

    Le test porte sur le point le plus proche du SEGMENT, pas sur ses
    extremites : un cordon long passant sous le point le supporte, meme si
    aucune de ses deux extremites n'est proche. La projection scalaire clampee
    est la meme que celle de `step9._point_to_segment_distance` (verifie par
    test contre cette fonction), mais decomposee en composantes verticale et
    horizontale - ce qu'une distance scalaire unique ne permet pas.

    VECTORISE sur tous les candidats en une passe numpy, jamais une boucle
    Python par capsule. Mesure sur death.stl : la grille de step9 retourne
    432 candidats en MEDIANE par requete (jusqu'a 2338), pour ~199 000
    requetes - une boucle interpretee par candidat dominait entierement le
    cout (191 s pour la seule passe principale). Meme nature de correction que
    la vectorisation de `_points_to_segment_distances` deja faite dans step9,
    pour exactement la meme raison.

    Delegue a `core.geometry_utils.point_has_support_along_direction` avec
    `direction_up=[0,0,1]` (verticale monde) - la geometrie de projection
    elle-meme vit desormais dans la primitive partagee (reutilisee par la
    classification multi-directionnelle ci-dessous ET par le mecanisme de
    reordonnancement de step5), ce fichier ne fait plus que la recherche
    spatiale (candidates_near) et le cablage de l'AABB verticale.
    """
    return index.has_support(
        _search_aabb_below(position, horizontal_radius, vertical_height),
        position,
        _WORLD_UP_SINGLETON,
        horizontal_radius,
        vertical_height,
    )


def _has_support_in_any_direction(
    index: _PrintedGeometryIndex,
    position: np.ndarray,
    directions: np.ndarray,
    horizontal_radius: float,
    vertical_height: float,
) -> bool:
    """
    Existe-t-il, dans `index`, un support pour `position` le long d'AU MOINS
    UNE des `directions` (D,3) ? Voir docstring de module, "Validite du test
    vertical" : un point sans support VERTICAL peut quand meme etre depose
    sans surplomb si la machine incline la tete/le plateau dans une
    direction atteignable.

    Une SEULE requete spatiale (`candidates_near`, AABB centree sur
    `position` et etendue de `max(horizontal_radius, vertical_height)` dans
    les 3 axes - plus large que `_search_aabb_below`, qui suppose une seule
    direction verticale) : les memes candidats sont ensuite testes pour
    chaque direction via `point_has_support_along_direction`, sans reinterroger
    l'index a chaque direction. Court-circuite des la premiere direction
    trouvee. Mesure sur death.stl : 0.26s pour 194 points x 26 directions.
    """
    extent = max(horizontal_radius, vertical_height)
    aabb = AABB(min_corner=position - extent, max_corner=position + extent)
    return index.has_support(
        aabb, position, np.asarray(directions, dtype=float), horizontal_radius, vertical_height
    )


def _build_full_index(
    points: list[ToolpathPoint], nozzle_diameter: float, cell_size_multiplier: float
) -> _PrintedGeometryIndex:
    """
    Index contenant TOUTE la matiere extrudee de la trajectoire, sans egard a
    l'ordre - c'est-a-dire l'etat FINAL de la piece imprimee.

    Sert a distinguer un defaut d'ordre d'un surplomb reel : si un point non
    supporte au moment T trouve un support dans cet index complet, c'est
    qu'une matiere deposee plus tard aurait pu le porter.
    """
    index = _PrintedGeometryIndex(nozzle_diameter, cell_size_multiplier)

    # Selection vectorisee des paires en extrusion, puis UNE seule insertion en
    # lot : la boucle Python d'origine traversait la frontiere une fois par
    # capsule (744 324 fois sur death.stl) pour une operation triviale.
    extruding = np.array([p.extrusion > 0.0 for p in points], dtype=bool)
    if len(extruding) < 2:
        return index
    pair_mask = extruding[1:] & extruding[:-1]
    if not np.any(pair_mask):
        return index

    positions = np.array([p.position for p in points], dtype=float)
    index.append_segments(positions[:-1][pair_mask], positions[1:][pair_mask])
    return index


def run(
    toolpath: ToolpathMachine,
    machine_profile: MachineProfile,
    layer_height_base: float,
    platform_z: Optional[float],
    params: SupportCheckParams,
) -> SupportCheckReport:
    """
    Verifie le support de chaque point de depot sur la trajectoire machine.

    Args:
        toolpath: sortie de step8 (ou step8b si le routage est actif).
        machine_profile: pour `nozzle_diameter` (echelle du cordon).
        layer_height_base: `LayerParams.layer_height_base`, passe
            explicitement - step9b ne le redevine jamais depuis la
            trajectoire, deux sources de verite qui divergeraient en silence.
        platform_z: niveau du plateau, source de verite =
            `IngestReport.platform_z_used` (meme valeur que celle donnee a
            step3). `None` = plateau inconnu : aucun point n'est alors
            exempte au titre du contact plateau, ce qui SURESTIME le nombre
            d'evenements plutot que de masquer un vrai defaut.
        params: seuils de recherche.

    Returns:
        SupportCheckReport. `checked=False` si la trajectoire est vide : signal
        explicite, jamais une liste vide ambigue.

    Ne leve jamais : un depot non supporte est une DONNEE, pas une erreur -
    meme discipline que `Step8Diagnostic` et `CollisionReport`.
    """
    points = toolpath.points
    if not points:
        return SupportCheckReport(checked=False)

    nozzle_diameter = max(float(machine_profile.nozzle_diameter), 1e-6)
    horizontal_radius = params.horizontal_search_radius_factor * nozzle_diameter
    vertical_height = params.vertical_search_height_factor * max(float(layer_height_base), 1e-6)
    platform_tolerance = params.platform_contact_tolerance_factor * nozzle_diameter

    # Directions candidates pour la classification finale "surplomb reel"
    # (voir docstring de module) - calculees UNE SEULE FOIS, memes fonction
    # et defaut que step3 (`_generate_candidate_directions`), y compris le
    # prefixage manuel de `up` (composed_rotation_directions n'inclut la
    # direction de reference que si un axe echantillonne l'angle 0, jamais
    # garanti en general). `rotation_axes` vide/None -> `up` seul, repli
    # exact sur le comportement historique (verticale uniquement).
    up = np.array([0.0, 0.0, 1.0])
    rotation_axes = machine_profile.rotation_axes or []
    if rotation_axes:
        samples_per_axis = _samples_per_axis(params.direction_candidate_count, len(rotation_axes))
        candidate_directions = np.vstack(
            [up.reshape(1, 3), composed_rotation_directions(rotation_axes, up, samples_per_axis)]
        )
    else:
        candidate_directions = up.reshape(1, 3)

    sample_indices = _sample_point_indices(points, params.path_sampling_step)

    printed_index = _PrintedGeometryIndex(nozzle_diameter, params.grid_cell_size_multiplier)
    last_grown_idx = 0
    tested = 0
    unsupported_indices: list[int] = []

    for point_idx in sample_indices:
        point = points[point_idx]
        # L'index grossit AVANT le test, exactement comme step9 : il ne
        # contient alors que la matiere strictement anterieure, et jamais le
        # cordon qui se termine au point teste (celui-ci est confondu avec la
        # buse - le compter comme support serait un faux negatif structurel).
        last_grown_idx = _grow_printed_index_to(
            printed_index, points, point_idx, last_grown_idx
        )

        if point.extrusion <= 0.0:
            continue
        tested += 1

        position = np.asarray(point.position, dtype=float)
        if platform_z is not None and position[2] <= platform_z + platform_tolerance:
            continue  # pose sur le plateau : aucun support en dessous requis

        if not _has_support_below(
            printed_index, position, horizontal_radius, vertical_height
        ):
            unsupported_indices.append(point_idx)

    events = _classify_unsupported(
        points,
        unsupported_indices,
        nozzle_diameter,
        horizontal_radius,
        vertical_height,
        params.grid_cell_size_multiplier,
        candidate_directions,
        params.locate_supporting_point,
    )

    return SupportCheckReport(events=events, checked=True, tested_extruded_points=tested)


def _classify_unsupported(
    points: list[ToolpathPoint],
    unsupported_indices: list[int],
    nozzle_diameter: float,
    horizontal_radius: float,
    vertical_height: float,
    cell_size_multiplier: float,
    candidate_directions: np.ndarray,
    locate_supporting_point: bool = True,
) -> list[UnsupportedDepositEvent]:
    """
    Separe les depots non supportes en defauts d'ORDRE, surplombs RESOLUS PAR
    INCLINAISON, et surplombs REELS.

    Trois passes de couts tres differents, dans cet ordre precis :

    1. **Classification ordre/surplomb** (bon marche, cout borne) : un index
       contenant TOUTE la matiere extrudee - l'etat FINAL de la piece - est
       construit une fois, puis chaque evenement y est teste une fois
       (verticalement). Un evenement qui y trouve un support a bien de la
       matiere porteuse quelque part : elle est simplement deposee trop tard
       (defaut d'ORDRE). Un evenement qui n'en trouve pas n'a pas de support
       vertical, meme dans la piece finie - candidat surplomb, a affiner par
       la passe 3. Cout : O(evenements), pas O(evenements x segments).

    2. **Localisation** du premier point qui aurait supporte chaque defaut
       d'ordre (`supporting_point_index`) : celle-ci exige de rejouer la
       croissance incrementale, donc O(defauts d'ordre x segments). Elle n'est
       lancee QUE sur les defauts d'ordre deja identifies par la passe 1, et
       peut etre desactivee (`locate_supporting_point=False`) quand seuls les
       totaux importent - la distinction ordre/surplomb, elle, reste toujours
       calculee.

    3. **Classification surplomb reel vs resolu par inclinaison** (voir
       docstring de module, "Validite du test vertical") : les candidats
       surplomb de la passe 1 sont retestes contre le MEME `final_index`,
       mais le long de CHACUNE de `candidate_directions` - un point qui
       trouve un support dans une direction inclinee est un point que la
       machine peut deposer sans surplomb en orientant la tete, PAS un
       defaut d'ordre (aucun reordonnancement ne le concerne) ni un vrai
       surplomb geometrique. Cout : O(candidats surplomb x directions),
       negligeable car cette passe ne s'applique jamais a la masse des
       defauts d'ordre (deja ecartes par la passe 1).

    Cette separation est ce qui rend l'etape utilisable : la version initiale
    testait TOUS les evenements en attente a CHAQUE segment de la trajectoire,
    ce qui ne terminait pas en 10 minutes sur le pire cas mesure (22 393
    evenements x 199 000 segments sur death.stl avant correction du champ).
    """
    if not unsupported_indices:
        return []

    # --- Passe 1 : ordre vs surplomb, sur l'etat FINAL de la piece ---------
    final_index = _build_full_index(points, nozzle_diameter, cell_size_multiplier)
    ordering_defects = [
        event_idx
        for event_idx in unsupported_indices
        if _has_support_below(
            final_index,
            np.asarray(points[event_idx].position, dtype=float),
            horizontal_radius,
            vertical_height,
        )
    ]

    # --- Passe 3 : parmi les candidats surplomb, lesquels sont resolus par
    # une orientation de tete atteignable ? (voir docstring de module)
    overhang_candidates = [idx for idx in unsupported_indices if idx not in set(ordering_defects)]
    tilt_resolved = {
        event_idx
        for event_idx in overhang_candidates
        if _has_support_in_any_direction(
            final_index,
            np.asarray(points[event_idx].position, dtype=float),
            candidate_directions,
            horizontal_radius,
            vertical_height,
        )
    }

    # --- Passe 2 : quel point precis aurait supporte ce defaut d'ordre ? ---
    supporting: dict[int, int] = {}
    if locate_supporting_point and ordering_defects:
        pending = set(ordering_defects)
        index = _PrintedGeometryIndex(nozzle_diameter, cell_size_multiplier)
        for i in range(1, len(points)):
            if not pending:
                break
            if not (points[i].extrusion > 0.0 and points[i - 1].extrusion > 0.0):
                continue
            index.append_segment(points[i - 1].position, points[i].position)

            # Un evenement ne peut etre resolu que par de la matiere
            # POSTERIEURE a lui : la matiere anterieure a deja ete testee par
            # la passe principale, qui l'a precisement declare non supporte.
            resolved = [
                event_idx
                for event_idx in pending
                if i > event_idx
                and _has_support_below(
                    index,
                    np.asarray(points[event_idx].position, dtype=float),
                    horizontal_radius,
                    vertical_height,
                )
            ]
            for event_idx in resolved:
                supporting[event_idx] = i
                pending.discard(event_idx)

    ordering_set = set(ordering_defects)
    return [
        UnsupportedDepositEvent(
            point_index=event_idx,
            position=np.asarray(points[event_idx].position, dtype=float),
            has_future_support=event_idx in ordering_set,
            supporting_point_index=supporting.get(event_idx),
            supported_by_tilt=event_idx in tilt_resolved,
        )
        for event_idx in unsupported_indices
    ]
