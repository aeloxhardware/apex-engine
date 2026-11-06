"""Etape 6 : ponts (direction, réorientation) et découpage des anneaux en tronçons monotones.

Découpé de `step6_toolpath2d.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step6_toolpath2d`.
"""

from __future__ import annotations

import numpy as np
from core.mesh_types import MeshData, OverhangMap
from pipeline.step6_geometry import (  # noqa: F401
    _offset_walls,
)
from pipeline.step6_infill import (  # noqa: F401
    _generate_infill,
)
from pipeline.step6_params import (  # noqa: F401
    Step6Params,
)


def _bridge_direction_for_loop(
    contour_3d_loop: np.ndarray,
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    local_frame: tuple[np.ndarray, np.ndarray, np.ndarray],
    threshold: float,
) -> np.ndarray | None:
    """
    Determine si une boucle de contour 3D touche une zone is_bridge (step2),
    via le meme pont kdtree surface/volume deja utilise par step4 (ancrage
    plateau) et step5 (association region). Retourne la direction de portee
    du pont projetee dans le plan local de la couche (unitaire), ou None si
    la boucle ne touche aucune zone de pont.

    Seuil de PRESENCE, pas de majorite : un pont est une sous-region locale
    de la matiere d'une couche, jamais une majorite du perimetre total de sa
    boucle de contour (mesure directe sur la fixture pont : ratio entre 0.06
    et 0.21 sur toutes les couches concernees) - voir CORRECTION 1 du plan.
    """
    _, nearest_faces = mesh_data.kdtree.query(contour_3d_loop)
    nearest_faces = np.asarray(nearest_faces, dtype=np.int64)

    is_bridge_per_point = overhang_map.is_bridge[nearest_faces]
    ratio = float(is_bridge_per_point.mean()) if len(is_bridge_per_point) else 0.0
    if ratio < threshold:
        return None

    span_vectors = overhang_map.bridge_span_direction[nearest_faces[is_bridge_per_point]]
    valid = ~np.isnan(span_vectors).any(axis=1)
    if not np.any(valid):
        return None

    mean_span = span_vectors[valid].mean(axis=0)
    norm = np.linalg.norm(mean_span)
    if norm < 1e-12:
        return None
    mean_span = mean_span / norm

    _origin, axis_u, axis_v = local_frame
    span_2d = np.array([np.dot(mean_span, axis_u), np.dot(mean_span, axis_v)])
    span_norm = np.linalg.norm(span_2d)
    if span_norm < 1e-12:
        return None
    return span_2d / span_norm


def _reorient_bridge_toolpath(
    outer: np.ndarray,
    holes: list[np.ndarray],
    span_direction_2d: np.ndarray,
    params: Step6Params,
) -> list[np.ndarray]:
    """
    Trajectoire d'un contour pont : parois normales (wall_count), PUIS
    bridge_anchor_extra_walls parois supplementaires resserrees sur le
    perimetre (renfort d'ancrage sur les bords, cf. plan step6 Decision
    explicite : parois supplementaires plutot que densite d'infill
    augmentee), PUIS infill balaye perpendiculairement a span_direction_2d
    (la trajectoire doit TRAVERSER le pont, donc etre perpendiculaire a son
    axe de plus grande portee - balayer selon la portee reviendrait a
    longer le vide plutot qu'a le franchir).
    """
    total_wall_levels = params.wall_count + params.bridge_anchor_extra_walls
    walls, innermost = _offset_walls(
        outer, holes, params.nozzle_width, total_wall_levels, params.degenerate_area_ratio
    )
    if not walls:
        return [outer]

    sweep_angle = float(np.arctan2(span_direction_2d[1], span_direction_2d[0])) + np.pi / 2.0
    infill = _generate_infill(
        innermost, params.infill_density, params.infill_pattern, params.nozzle_width, sweep_angle
    )

    return infill + walls


def _monotonic_runs(offsets: np.ndarray, min_amplitude: float = 1e-12) -> list[tuple[int, int]]:
    """
    Decoupe un anneau en troncons [start, stop) D'INDICES DISJOINTS (jamais de
    sommet partage entre deux troncons - un decoupage "un sommet partage a
    chaque frontiere" semble naturel pour preserver la continuite du trace,
    mais duplique ce sommet des que les troncons deviennent des polylignes
    independantes : bug reel mesure, l'anneau ressortait avec PLUS de points
    que l'original).

    Troncons monotones en offset (croissant OU decroissant). Un changement de
    sens ne coupe que si l'offset a recule de plus de `min_amplitude` depuis
    l'extremum courant (hysteresis) : une oscillation plus petite est absorbee
    par le troncon en cours. La coupure est placee juste apres l'extremum,
    comme avant. Un anneau plat (tous offsets egaux, cas historique sans
    interpolation hors-plan) produit un seul troncon couvrant tout l'anneau -
    non-regression : `_split_ring_into_runs` renvoie alors l'anneau inchange,
    en une seule polyligne.
    """
    n = len(offsets)
    if n <= 2:
        return [(0, n)]

    runs: list[tuple[int, int]] = []
    start = 0
    direction = 0                      # 0 = indeterminee, 1 = montee, -1 = descente
    low = high = float(offsets[0])     # bornes tant que la direction est indeterminee
    ext_index = 0
    ext_value = float(offsets[0])      # extremum courant dans la direction suivie
    for i in range(1, n):
        value = float(offsets[i])
        if direction == 0:
            low, high = min(low, value), max(high, value)
            if value - low > min_amplitude:
                direction, ext_index, ext_value = 1, i, value
            elif high - value > min_amplitude:
                direction, ext_index, ext_value = -1, i, value
        elif direction == 1:
            if value >= ext_value:
                ext_index, ext_value = i, value
            elif ext_value - value > min_amplitude:
                runs.append((start, ext_index + 1))
                start = ext_index + 1
                direction, ext_index, ext_value = -1, i, value
        else:
            if value <= ext_value:
                ext_index, ext_value = i, value
            elif value - ext_value > min_amplitude:
                runs.append((start, ext_index + 1))
                start = ext_index + 1
                direction, ext_index, ext_value = 1, i, value
    runs.append((start, n))
    return runs


def _is_closed_ring(poly: np.ndarray) -> bool:
    """
    Une ligne d'infill ouverte (segment de balayage clippe, `_generate_infill`)
    n'est PAS un anneau. Un contour matiere l'est TOUJOURS par construction
    (`_offset_walls`/`_reorient_bridge_toolpath`), mais sa representation n'est
    pas uniforme : les parois n>=1 (generees par `shapely.buffer`) dupliquent
    leur premier point en derniere position (convention shapely), alors que le
    contour original n=0 (issu de `_classify_contours` -> `open_loops`,
    dedoublonne) ne le fait PAS - bug reel trouve en verifiant ce mecanisme :
    un test d'egalite stricte poly[0]==poly[-1] classait le contour n=0 comme
    une ligne OUVERTE et le laissait hors du decoupage, alors que c'est
    geometriquement le meme genre d'anneau que les autres parois.

    D'ou un test ROBUSTE aux deux representations : au moins 3 points, ET
    (premier == dernier point OU le dernier segment revient pres du premier
    point a une distance du meme ordre que le pas entre deux points
    consecutifs). Une ligne d'infill n'a structurellement aucune raison de
    revenir pres de son point de depart (droite clippee, extremites sur des
    bords opposes du polygone).
    """
    if len(poly) < 3:
        return False
    if np.allclose(poly[0, :2], poly[-1, :2]):
        return True
    closing_gap = float(np.hypot(*(poly[-1, :2] - poly[0, :2])))
    step_lengths = np.hypot(*(np.diff(poly[:, :2], axis=0)).T)
    typical_step = float(np.median(step_lengths)) if len(step_lengths) else 0.0
    return typical_step > 0.0 and closing_gap <= 2.0 * typical_step


def _merge_short_runs(
    runs: list[tuple[int, int]], offsets: np.ndarray | None = None
) -> list[tuple[int, int]]:
    """
    Fusionne tout troncon de moins de 2 points dans un voisin.

    Un troncon a UN SEUL point deviendrait une polyligne d'un point : step8 y
    insererait un trajet a vide puis un point isole, et
    `_grow_printed_index_to` (qui exige DEUX extremites en extrusion) n'en
    ferait aucun cordon - ce point cesserait purement et simplement d'etre
    imprime. Fusionner preserve la partition exacte (aucun point perdu ni
    duplique).

    LA DIRECTION DE LA FUSION N'EST PAS NEUTRE - c'est un bug reel mesure, pas
    une precaution theorique. Une version precedente fusionnait TOUJOURS en
    arriere (le point isole colle a la fin du troncon precedent) et
    documentait la consequence comme un detail acceptable ("un troncon qui
    n'est plus strictement monotone sur un seul de ses points"). Mesure sur
    res/test_models/death.stl : ce "detail" causait **425 des 1055 defauts
    d'ordre (40.3%)**. Le motif exact, verifie : le champ scalaire oscille
    faiblement (bruit du marching tetrahedra, dz median mesure -0.0141mm), un
    point isole forme un micro-creux, la fusion arriere le recolle a la FIN
    d'un troncon montant - il est donc depose APRES un point plus haut que
    lui, sans rien dessous. Sur 425 cas, 425 (100%) avaient leur predecesseur
    immediat plus haut qu'eux DANS LA MEME POLYLIGNE continue.

    D'ou la regle ici : fusionner en arriere UNIQUEMENT si le point isole est
    a une altitude >= au dernier point du troncon precedent (la montee reste
    vraie). Sinon, le reporter en AVANT - il devient le PREMIER point du
    troncon suivant, donc son point le plus bas, ce qui est exactement sa
    place physique. `offsets=None` conserve l'ancien comportement (fusion
    arriere systematique) pour les appelants qui n'ont pas l'altitude sous la
    main.
    """
    if offsets is None:
        merged: list[tuple[int, int]] = []
        for start, stop in runs:
            if merged and (stop - start) < 2:
                merged[-1] = (merged[-1][0], stop)
            else:
                merged.append((start, stop))
        if len(merged) >= 2 and (merged[0][1] - merged[0][0]) < 2:
            merged[1] = (merged[0][0], merged[1][1])
            merged.pop(0)
        return merged

    merged = []
    carried_start: int | None = None      # point isole reporte sur le troncon suivant
    for start, stop in runs:
        if carried_start is not None:
            start = carried_start
            carried_start = None
        if (stop - start) >= 2:
            merged.append((start, stop))
            continue
        # troncon trop court : choisir la direction qui preserve la montee
        if merged and offsets[start] >= offsets[merged[-1][1] - 1]:
            merged[-1] = (merged[-1][0], stop)
        else:
            carried_start = start
    if carried_start is not None:
        # plus aucun troncon apres : la fusion arriere est le seul recours
        # (jamais un point perdu, meme au prix de la monotonie sur ce point)
        if merged:
            merged[-1] = (merged[-1][0], runs[-1][1])
        else:
            merged.append((carried_start, runs[-1][1]))
    return merged


def _split_ring_into_runs(
    ring: np.ndarray, share_vertices: bool = False, min_amplitude: float = 1e-12
) -> list[np.ndarray]:
    """
    Decoupe un anneau FERME (K,3) (u, v, offset) en troncons monotones, chacun
    stocke dans le sens MONTANT - voir docstring de section pour la cause
    racine mesuree et pourquoi les troncons restent des polylignes SEPAREES.

    Ne cree ni ne supprime aucun point : la concatenation des troncons
    retournes contient exactement les memes points que `ring`, seul l'ordre
    INTERNE d'un troncon descendant est inverse. L'ordre DES troncons entre eux
    est celui de l'anneau d'origine - le choisir demande l'etat reel de la
    matiere deja deposee, que ce fichier n'a pas (cf. step7b).

    Un anneau plat ou trop court ressort en un seul troncon, identique a
    l'entree (non-regression stricte).
    """
    if len(ring) <= 2:
        return [ring]

    runs = _monotonic_runs(ring[:, 2], min_amplitude)
    if len(runs) <= 1:
        return [ring]

    segments: list[np.ndarray] = []
    for k, (start, stop) in enumerate(_merge_short_runs(runs, ring[:, 2])):
        # `share_vertices` : le troncon reprend le dernier sommet du precedent,
        # pour que le segment de jonction soit extrude (voir
        # Step6Params.continuous_rings) ; il suit le sens de ce troncon-ci.
        if share_vertices and k > 0 and start > 0:
            start -= 1
        segment = ring[start:stop]
        if len(segment) >= 2 and segment[-1, 2] < segment[0, 2]:
            segment = segment[::-1]
        segments.append(segment)
    return segments


def _split_rings_in_toolpath(
    toolpath: list[np.ndarray],
    include_open_polylines: bool = False,
    share_vertices: bool = False,
    min_amplitude: float = 1e-12,
) -> list[np.ndarray]:
    """
    Applique `_split_ring_into_runs` a chaque polyligne FERMEE de `toolpath`
    (parois, contour original).

    `include_open_polylines` etend le MEME traitement aux polylignes OUVERTES
    (lignes d'infill balayees, cf. `_generate_infill`). Une version precedente
    les excluait au motif qu'"une droite clippee n'a pas de retour au point de
    depart a gerer" - exact pour la FERMETURE, faux pour la MONOTONIE : sur une
    couche courbe, une droite clippee monte et descend le long de son parcours
    exactement comme un anneau.

    Mesure sur res/test_models/death.stl (5 couches les plus touchees,
    simulation directe sur donnees capturees) : les polylignes ouvertes
    portaient **211 des 243 paires descendantes** d'une couche (87%), contre 32
    pour les anneaux fermes - le decoupage historique ne couvrait donc que
    ~5% des paires de la trajectoire (174 sur 3797). Leur appliquer le meme
    decoupage en troncons monotones ramene 211 descentes a 9 (-95.7%).
    Retourner simplement la polyligne descendante ne suffit PAS (mesure : -4.3%
    seulement) : une ligne d'infill monte ET descend, l'inverser deplace le
    probleme sans le resoudre.

    Le nombre de polylignes AUGMENTE - c'est voulu : chaque frontiere de
    troncon devient une vraie frontiere de segment, donc un vrai trajet a vide
    insere par step8._iter_schedule_points, et non un cordon parasite extrude
    en travers de la piece (voir docstring de section). Le cout en trajets a
    vide supplementaires est donc a mesurer, jamais suppose negligeable.
    """
    result: list[np.ndarray] = []
    for poly in toolpath:
        if _is_closed_ring(poly) or include_open_polylines:
            result.extend(_split_ring_into_runs(poly, share_vertices, min_amplitude))
        else:
            result.append(poly)
    return result
