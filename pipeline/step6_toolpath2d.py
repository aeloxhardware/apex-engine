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
from core.geometry_utils import point_in_polygon_2d
from core.mesh_types import Layer, MeshData, OverhangMap


class ToolpathGenerationError(Exception):
    """Levee si une couche ne peut produire aucune trajectoire exploitable."""


class InfillPattern(Enum):
    # lignes paralleles independantes, meme orientation a toutes les couches
    # (comportement historique)
    LINES = "lines"
    # deux balayages croises a 90 deg dans CHAQUE couche, chacun a la moitie
    # de la densite
    GRID = "grid"
    # lignes croisees-alternees : un zigzag continu par zone, oriente selon
    # l'axe u du repere local une couche sur deux et selon l'axe v l'autre
    # (voir `_infill_angle_for_layer`)
    CROSSED_ZIGZAG = "crossed_zigzag"
    # spirale d'Archimede autour de l'AXE DU PLATEAU, oscillante : r(theta)
    # avec theta qui ne recule jamais (le plateau tourne toujours dans le meme
    # sens, les pics sont des mouvements radiaux). MEME ligne mediane a toutes
    # les couches (elles se superposent) ; seules les asperites different
    # d'une couche a l'autre (phase opposee, tirages differents).
    OSCILLATING_SPIRAL = "oscillating_spiral"


@dataclass
class Step6Params:
    nozzle_width: float = 0.4
    wall_count: int = 2
    infill_density: float = 0.2          # dans (0, 1]
    infill_pattern: InfillPattern = InfillPattern.CROSSED_ZIGZAG
    # seuil de PRESENCE (pas de majorite) : proportion minimale de points
    # d'une boucle de contour dont la face de surface la plus proche est
    # is_bridge, au-dela de laquelle la boucle entiere est traitee comme
    # trajectoire de pont. Un pont n'est qu'une petite portion locale du
    # perimetre total d'une couche (CLAUDE.md : zones is_bridge "qui
    # retombent dans cette couche", pas la couche entiere) - un seuil de
    # majorite ne peut structurellement jamais se declencher (mesure sur la
    # fixture pont : ratio reel entre 0.06 et 0.21 sur toutes les couches
    # concernees, jamais proche de 0.5). D'ou un defaut bas.
    bridge_contour_ratio_threshold: float = 0.05
    # parois supplementaires resserrees sur le perimetre d'un contour pont,
    # en plus de wall_count, avant de generer l'infill balaye - lecture
    # directe de CLAUDE.md "augmente la densite d'ancrage sur les bords"
    # (renfort de perimetre, pas de densite d'infill - cf. plan step6).
    bridge_anchor_extra_walls: int = 2
    # aire residuelle (relative a l'aire du contour d'origine) sous laquelle
    # un offset est considere comme ayant fait disparaitre la zone.
    degenerate_area_ratio: float = 1e-6
    # Exposant de la ponderation inverse-distance utilisee pour donner un
    # offset hors-plan aux points GENERES par shapely (parois n>=1, infill),
    # interpole depuis les offsets REELS connus sur les boucles de contour.
    # None = DESACTIVE (offset nul pour ces points, comportement historique) -
    # meme discipline opt-in que direction_blend_radius=0.0 (step4) ou
    # max_region_merge_angle_deg=None (step3).
    #
    # Sans ce mecanisme, ces points restent colles au plan local MOYEN de la
    # couche alors que la couche reelle s'en ecarte (mesure sur le modele en
    # T : jusqu'a 0.845 mm d'ecart), et peuvent donc tomber hors du solide -
    # residu que la limite #31 de CLAUDE.md documentait comme hors de son
    # perimetre. Le contour n=0 n'est JAMAIS concerne : il garde ses offsets
    # exacts, mesures, jamais interpoles.
    out_of_plane_interpolation_power: float | None = None

    # Decoupe chaque ANNEAU FERME (paroi, contour original) en troncons
    # monotones, chacun emis comme une POLYLIGNE SEPAREE et parcouru dans le
    # sens montant - voir Decision 5 plus bas. False = DESACTIVE (anneau
    # entier en une polyligne, comportement historique), meme discipline
    # opt-in que out_of_plane_interpolation_power ci-dessus. Necessite ce
    # dernier actif pour avoir un effet reel : sans lui l'offset hors-plan
    # des points generes est nul, chaque anneau est deja plat et ne produit
    # qu'un seul troncon.
    split_rings_into_runs: bool = False

    # Etend le decoupage en troncons monotones (ci-dessus) aux polylignes
    # OUVERTES - lignes d'infill balayees, cf. `_generate_infill` - et plus
    # seulement aux anneaux fermes. False = DESACTIVE (comportement
    # historique), meme discipline opt-in que les champs ci-dessus.
    #
    # Mesure sur death.stl (donnees capturees, 5 couches les plus touchees) :
    # les polylignes ouvertes portent 211 des 243 paires descendantes d'une
    # couche (87%) - le decoupage historique, limite aux anneaux, ne couvrait
    # que ~5% des paires de la trajectoire. Voir `_split_rings_in_toolpath`.
    split_open_polylines: bool = False

    # Continuite des parois (2026-10-01). Deux segments n'etaient JAMAIS
    # extrudes, mesure sur le pipeline par defaut :
    #   - la jonction entre deux troncons monotones (`_split_ring_into_runs`) :
    #     troncons a indices DISJOINTS, donc le segment qui les relie devient
    #     un trajet a vide - 7,7 % de la longueur deposee perdue sur death.stl,
    #     15,7 % sur l_shape_overhang.stl (des trous dans les parois) ;
    #   - la fermeture du contour d'origine n=0 et de ses trous, emis OUVERTS
    #     (0,1 % et 0,85 %).
    # True : chaque troncon reprend le dernier sommet du precedent (la
    # jonction reste monotone, elle prend le sens du troncon suivant), et les
    # contours n=0 recoivent leur point de fermeture. False = historique.
    continuous_rings: bool = False

    # Amplitude minimale (mm) d'un changement de sens pour COUPER un troncon
    # monotone. Le defaut 1e-12 reproduit l'ancien comportement (toute
    # oscillation, meme numerique, coupe). Mesure sur un cylindre plein
    # (couche plate) : l'anneau exterieur a 257 points et une etendue d'offset
    # de 3,8 um, soit 122 changements de sens -> 160 polylignes par couche au
    # lieu de 4, chacune avec son trajet a vide. Une oscillation sous cette
    # amplitude est du bruit du marching tetrahedra, bien en dessous de
    # l'epaisseur d'un cordon : elle ne cree aucun appui manquant.
    run_split_min_amplitude: float = 1e-12

    # --- OSCILLATING_SPIRAL (opt-in) --------------------------------------
    # Point de l'axe du plateau (x, y, z, repere machine) : centre de la
    # spirale. Source de verite unique = MachineProfile.platform_origin
    # (renseigne par main) ; None = centre de la zone (la monotonie de theta
    # n'est alors garantie qu'autour de ce centre, pas autour de l'axe).
    spiral_axis_origin: tuple[float, float, float] | None = None
    # "Seuil d'Alternance" (mm) : ecart radial des pics. None = 0,6 x la
    # limite (pas - buse) / 2. Toujours plafonne a (pas - buse) / 2 : deux
    # tours voisins ne se touchent jamais.
    spiral_alternation_threshold: float | None = None
    # Distance (mm, le long de la trajectoire) entre deux pics. None = 4 buses.
    spiral_peak_spacing: float | None = None
    # "Seuil de Randomisation" dans [0, 1] : 0 = oscillation reguliere,
    # 1 = amplitude, espacement, forme (pic net / courbe) et pics purement
    # radiaux tires au hasard. Graine fixe : resultat reproductible.
    spiral_randomization_threshold: float = 0.5
    spiral_seed: int = 0


# ---------------------------------------------------------------------------
# Projection / relevement repere local <-> 3D
# ---------------------------------------------------------------------------

def _project_to_local_2d(
    points_3d: np.ndarray, local_frame: tuple[np.ndarray, np.ndarray, np.ndarray]
) -> np.ndarray:
    """Projection PURE (u,v), utilisee par la geometrie interne 2D
    (classification, offset shapely, infill) - jamais de 3e colonne ici,
    ces algorithmes ne doivent voir que le plan."""
    origin, axis_u, axis_v = local_frame
    relative = points_3d - origin
    u = relative @ axis_u
    v = relative @ axis_v
    return np.stack([u, v], axis=1)


def _frame_normal(axis_u: np.ndarray, axis_v: np.ndarray) -> np.ndarray:
    """axis_u x axis_v, au bit pres comme np.cross (memes produits puis
    differences), sans son cout d'appel (~25 us) : appele une fois par
    polyligne relevee, soit ~29 000 fois sur death.stl (0,8 s mesurees)."""
    u0, u1, u2 = float(axis_u[0]), float(axis_u[1]), float(axis_u[2])
    v0, v1, v2 = float(axis_v[0]), float(axis_v[1]), float(axis_v[2])
    return np.array([u1 * v2 - u2 * v1, u2 * v0 - u0 * v2, u0 * v1 - u1 * v0])


def _out_of_plane_offset(
    points_3d: np.ndarray, local_frame: tuple[np.ndarray, np.ndarray, np.ndarray]
) -> np.ndarray:
    """
    Distance signee de chaque point au plan local (le long de la normale
    axis_u x axis_v), perdue par _project_to_local_2d. Un point 3D dont le
    contour reel ondule hors du plan approximatif de step5 (Decision 6, cf.
    CLAUDE.md) a une distance non nulle - preservee ici pour que
    _lift_to_3d puisse reconstruire EXACTEMENT le point d'origine plutot
    que sa projection sur le plan moyen (voir plan de correction du
    debordement, CLAUDE.md).
    """
    origin, axis_u, axis_v = local_frame
    normal = _frame_normal(axis_u, axis_v)
    return (points_3d - origin) @ normal


def _to_toolpath_3col(points_2d: np.ndarray, out_of_plane: np.ndarray | None = None) -> np.ndarray:
    """
    Assemble le format final (K,3) de Layer.toolpath_2d : (u, v, offset
    hors-plan). `out_of_plane` absent/None = colonne de zeros - le cas des
    parois offsettees et de l'infill quand
    `Step6Params.out_of_plane_interpolation_power` est None (ces points 2D
    nouveaux n'ont aucun equivalent 3D d'origine a restaurer) ; quand le
    mecanisme est actif, ils recoivent un offset INTERPOLE depuis les
    boucles de contour (`_interpolate_out_of_plane`). Nom du champ conserve
    ("toolpath_2d") malgre la 3e colonne : c'est toujours une trajectoire
    du repere LOCAL de la couche (pas un repere 3D machine), l'offset n'est
    qu'un complement pour un relevage fidele, pas une nouvelle dimension
    independante.
    """
    if out_of_plane is None:
        out_of_plane = np.zeros(len(points_2d))
    return np.concatenate([points_2d, out_of_plane[:, None]], axis=1)


def _lift_to_3d(
    points_2d_or_3col: np.ndarray, local_frame: tuple[np.ndarray, np.ndarray, np.ndarray]
) -> np.ndarray:
    """
    Releve une polyligne du repere local de la couche vers le 3D. Accepte
    (K,2) (comportement historique, offset hors-plan implicitement nul) ou
    (K,3) (u,v,offset - cf. _to_toolpath_3col) : la colonne d'offset,
    quand presente, reinjecte la composante hors-plan preservee par
    _out_of_plane_offset, pour que le point releve retombe EXACTEMENT sur
    le point 3D d'origine plutot que sur le plan moyen approximatif de
    step5 - voir CLAUDE.md, plan de correction du debordement.
    """
    origin, axis_u, axis_v = local_frame
    lifted = origin + points_2d_or_3col[:, 0:1] * axis_u + points_2d_or_3col[:, 1:2] * axis_v
    if points_2d_or_3col.shape[1] >= 3:
        normal = _frame_normal(axis_u, axis_v)
        lifted = lifted + points_2d_or_3col[:, 2:3] * normal
    return lifted


# Deux points 2D plus proches que cette distance sont traites comme
# confondus (on reprend alors la valeur source telle quelle, plutot que de
# diviser par une distance quasi nulle).
_COINCIDENT_EPS = 1e-9


def _interpolate_out_of_plane(
    query_2d: np.ndarray,
    source_2d: np.ndarray,
    source_offsets: np.ndarray,
    power: float,
) -> np.ndarray:
    """
    Interpole l'offset hors-plan en des points 2D GENERES (parois offsettees,
    infill - sans equivalent 3D d'origine), a partir des offsets REELS connus
    aux points des boucles de contour de la meme zone.

    Ponderation inverse-distance (`1/d**power`) : meme famille que
    `step4_field_solver._blend_tet_directions`, mais domaine different
    (scalaire dans le plan 2D d'une couche, vs direction 3D par tetraedre) -
    d'ou une fonction distincte plutot qu'une reutilisation forcee.

    `scipy.interpolate.LinearNDInterpolator` a ete ecarte : la triangulation
    de Delaunay d'un ANNEAU de points ne couvre que l'enveloppe CONVEXE, donc
    rendrait NaN a l'interieur de tout contour non convexe - cas courant ici.

    Args:
        query_2d: (K,2) points a evaluer.
        source_2d: (M,2) points de contour (u,v) aux offsets connus.
        source_offsets: (M,) offsets hors-plan reels, alignes sur source_2d.
        power: exposant de la ponderation (2.0 typiquement).

    Returns:
        (K,) offsets interpoles. Toujours borne par [min, max] des sources
        (propriete de la moyenne ponderee) : ce mecanisme ne peut pas creer
        un offset plus extreme que ceux reellement mesures sur le contour.
    """
    query_2d = np.asarray(query_2d, dtype=np.float64)
    source_2d = np.asarray(source_2d, dtype=np.float64)
    source_offsets = np.asarray(source_offsets, dtype=np.float64)

    if len(query_2d) == 0:
        return np.zeros(0)
    if len(source_2d) == 0:
        return np.zeros(len(query_2d))

    if _cpp.enabled():
        # Transcription bit-exacte de la suite - voir
        # cpp/include/apex/interpolate.hpp. Le gain vient de ne jamais
        # materialiser la matrice de distances (K,M) : une seule passe par
        # requete, plus proche voisin et les deux accumulations comprises.
        return _cpp.module().interpolate_out_of_plane(
            np.ascontiguousarray(query_2d),
            np.ascontiguousarray(source_2d),
            np.ascontiguousarray(source_offsets),
            float(power),
        )

    distances = np.linalg.norm(query_2d[:, None, :] - source_2d[None, :, :], axis=2)
    nearest = np.argmin(distances, axis=1)

    # 0**-power exploserait : on neutralise les distances quasi nulles ici et
    # on rend la valeur source exacte pour ces lignes juste apres.
    coincident = distances < _COINCIDENT_EPS
    safe_distances = np.where(coincident, 1.0, distances)
    weights = np.where(coincident, 0.0, safe_distances ** (-power))

    total = weights.sum(axis=1)
    weighted = (weights * source_offsets[None, :]).sum(axis=1)

    # `total <= 0` couvre le cas limite d'un sous-debordement numerique sur
    # des distances enormes : on retombe alors sur la source la plus proche,
    # jamais sur une division par zero.
    use_nearest = coincident.any(axis=1) | (total <= 0.0)
    return np.where(use_nearest, source_offsets[nearest], weighted / np.where(total > 0.0, total, 1.0))


def _out_of_plane_for_polyline(
    poly: np.ndarray,
    outer: np.ndarray,
    outer_out_of_plane: np.ndarray | None,
    source_2d: np.ndarray | None,
    source_offsets: np.ndarray | None,
    power: float | None,
) -> np.ndarray | None:
    """
    Offset hors-plan a attacher a une polyligne de trajectoire.

    Le contour n=0 (identique au contour d'origine) garde ses offsets EXACTS,
    mesures - jamais interpoles : c'est l'invariant du fix de la limite #31,
    et la seule polyligne dont la position 3D d'origine est connue.

    Tout le reste (parois n>=1, infill) est interpole si le mecanisme est
    actif, sinon `None` (colonne de zeros, comportement historique).
    """
    if _loop_matches(poly, outer):
        return outer_out_of_plane
    if power is None or source_2d is None or source_offsets is None or len(source_2d) == 0:
        return None
    return _interpolate_out_of_plane(poly[:, :2], source_2d, source_offsets, power)


# ---------------------------------------------------------------------------
# Decision 1 - Classification contour exterieur / trou par imbrication
# ---------------------------------------------------------------------------

def _signed_area_2d(polygon: np.ndarray) -> float:
    x = polygon[:, 0]
    y = polygon[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def _classify_contours(
    contours_2d: list[np.ndarray],
) -> list[tuple[np.ndarray, list[np.ndarray]]]:
    """Forme historique de `_classify_contour_indices` : les boucles OUVERTES
    plutot que leurs indices (conservee pour les appelants et les tests)."""
    open_loops, groups = _classify_contour_indices(contours_2d)
    return [(open_loops[o], [open_loops[h] for h in holes]) for o, holes in groups]


def _classify_contour_indices(
    contours_2d: list[np.ndarray],
) -> tuple[list[np.ndarray], list[tuple[int, list[int]]]]:
    """
    Rend (boucles ouvertes, [(indice matiere, [indices trous])]) - les indices
    sont ceux de `contours_2d`. Permet a `run` d'associer chaque boucle a son
    groupe SANS comparer toutes les boucles a tous les contours : la version
    par comparaison (`_loop_matches`) coutait 193 331 appels a np.allclose sur
    turkish.stl (O(groupes x boucles x trous)).

    Classe les boucles projetees en groupes (contour_matiere, [trous]) par
    imbrication geometrique (point_in_polygon_2d), PAS par signe d'aire
    (hypothese refutee empiriquement - cf. CLAUDE.md limite #11).

    Regle de parite : pour une boucle donnee, le nombre de boucles (parmi
    les plus grandes qu'elle) qui la contiennent determine son role -
    pair (0, 2, 4...) = matiere, impair (1, 3...) = trou. Generalise le cas
    a 2 boucles a un nombre arbitraire de boucles imbriquees (anneau dans
    anneau).

    Chaque trou est ensuite rattache au plus petit contour "matiere" qui le
    contient directement (son parent immediat), pas necessairement le plus
    grand contour de la couche - important des qu'il y a plusieurs ilots
    matiere disjoints.
    """
    if not contours_2d:
        return [], []

    # une seule boucle sans point de fermeture duplique pour les calculs
    open_loops = [c[:-1] if np.allclose(c[0], c[-1]) else c for c in contours_2d]
    areas = [abs(_signed_area_2d(loop)) for loop in open_loops]
    order = sorted(range(len(open_loops)), key=lambda i: areas[i], reverse=True)

    # containers[i] = liste des indices de boucles (parmi celles deja
    # traitees, donc plus grandes) qui contiennent la boucle i
    containers: dict[int, list[int]] = {i: [] for i in range(len(open_loops))}
    for rank, i in enumerate(order):
        probe_point = open_loops[i][0]
        for j in order[:rank]:
            if point_in_polygon_2d(probe_point, open_loops[j]):
                containers[i].append(j)

    is_hole = {i: (len(containers[i]) % 2 == 1) for i in range(len(open_loops))}

    # rattache chaque trou a son parent immediat DE MATIERE = le plus petit
    # englobant NON-TROU. Un trou peut etre imbrique a plus de 2 niveaux
    # (matiere > trou > matiere > trou...) : son plus petit conteneur au
    # sens de l'aire peut lui-meme etre classe trou (parite impaire), auquel
    # cas ce n'est pas un parent valide pour `groups` (qui n'a de cle que
    # pour les contours matiere) - remonte alors au conteneur non-trou de
    # plus petite aire. Bug reel mesure sur res/test_models/death.stl
    # (KeyError) : la version precedente prenait le plus petit conteneur
    # SANS filtrer par is_hole, plantant des qu'il etait lui-meme un trou.
    groups: dict[int, list[int]] = {i: [] for i in range(len(open_loops)) if not is_hole[i]}
    for i in range(len(open_loops)):
        if not is_hole[i]:
            continue
        matter_parents = [p for p in containers[i] if not is_hole[p]]
        if not matter_parents:
            # trou sans parent matiere (topologiquement incoherent pour une
            # isosurface bien formee) : ignore plutot que de planter -
            # meme esprit defensif que step5 (_map_tets_to_regions).
            continue
        parent = min(matter_parents, key=lambda p: areas[p])
        groups[parent].append(i)

    return open_loops, [(outer_idx, holes) for outer_idx, holes in groups.items()]


# ---------------------------------------------------------------------------
# Decision 3 - Offset polygonal progressif (parois), arret propre si degenere
# ---------------------------------------------------------------------------

def _to_shapely_polygon(outer: np.ndarray, holes: list[np.ndarray]) -> Polygon:
    return Polygon(outer, holes=[h for h in holes])


def _offset_walls(
    outer: np.ndarray, holes: list[np.ndarray], nozzle_width: float, wall_count: int, degenerate_area_ratio: float
) -> tuple[list[np.ndarray], Polygon | None]:
    """
    Genere les parois d'un contour matiere par offset progressif vers
    l'interieur (n = wall_count .. 1, puis le contour original n=0),
    ordonnees du centre vers l'exterieur dans la liste retournee.

    S'arrete proprement (sans erreur) des qu'un offset fait disparaitre la
    zone - une paroi plus fine que nozzle_width * wall_count ne peut pas
    recevoir wall_count parois pleines, cf. plan step6 Question 3. Retourne
    aussi la geometrie du dernier offset valide (pour servir de frontiere a
    l'infill), ou None si meme n=1 a disparu.

    Chaque niveau trace TOUS les bords de TOUS les morceaux de l'offset :
    exterieur ET trous. Correction (2026-09-26) de deux defauts mesures sur
    les couches Bird du Benchy :
      - seul le bord exterieur etait trace (`polygon.exterior`), et seul le
        contour exterieur d'origine etait ajoute a la fin : les trous ne
        recevaient AUCUNE paroi, ni leur contour, ni les parois autour ;
      - quand l'offset coupait la zone en plusieurs morceaux, seul le plus
        grand etait garde (`_largest_polygon`) : 99 couches du Benchy
        perdaient ainsi parois interieures et remplissage, 561 mm2 au total.
    Le seuil de degenerescence s'applique desormais morceau par morceau.
    """
    base_polygon = _to_shapely_polygon(outer, holes)
    if base_polygon.is_empty or base_polygon.area <= 0:
        return [], None

    reference_area = base_polygon.area
    walls: list[tuple[int, np.ndarray]] = []
    innermost = None

    for n in range(wall_count, 0, -1):
        offset = base_polygon.buffer(-nozzle_width * n, join_style="mitre")
        pieces = [p for p in _polygon_pieces(offset) if p.area > degenerate_area_ratio * reference_area]
        if not pieces:
            continue
        for piece in pieces:
            walls.append((n, np.asarray(piece.exterior.coords)))
            walls.extend((n, np.asarray(ring.coords)) for ring in piece.interiors)
        if innermost is None:
            innermost = pieces[0] if len(pieces) == 1 else MultiPolygon(pieces)

    # parois du centre vers l'exterieur : n decroissant a deja ce sens
    ordered = [pts for _n, pts in walls]
    ordered.append(outer)  # n=0, contour original, en dernier (le plus exterieur)
    ordered.extend(holes)  # n=0 des trous : le bord visible de chaque cavite

    return ordered, innermost


def _polygon_pieces(geometry) -> list[Polygon]:
    """Tous les morceaux polygonaux d'un resultat de `buffer` (Polygon,
    MultiPolygon ou vide)."""
    if geometry.is_empty:
        return []
    if geometry.geom_type == "Polygon":
        return [geometry]
    if geometry.geom_type == "MultiPolygon":
        return list(geometry.geoms)
    return []



# ---------------------------------------------------------------------------
# Decision 4 - Infill
# ---------------------------------------------------------------------------

_MAX_INFILL_LINES = 5000
# longueur maximale d'un raccord de zigzag, en espacements de balayage
_ZIGZAG_MAX_CONNECTOR_SPACINGS = 2.5

def _generate_infill(
    boundary: Polygon,
    density: float,
    pattern: InfillPattern,
    nozzle_width: float,
    sweep_angle_rad: float = 0.0,
) -> list[np.ndarray]:
    """
    Genere l'infill a l'interieur de `boundary` par balayage de lignes
    paralleles clippees par le polygone (shapely.intersection), selon
    `pattern` :
      - LINES : lignes independantes (une polyligne de 2 points chacune),
        orientees selon `sweep_angle_rad` ;
      - GRID : deux balayages croises a 90 deg, chacun a densite / 2 ;
      - CROSSED_ZIGZAG : les lignes d'un meme balayage sont reliees bout a
        bout en zigzag (`_connect_zigzag`) ; l'alternance entre couches est
        portee par l'appelant via `sweep_angle_rad`
        (`_infill_angle_for_layer`).

    L'espacement entre lignes d'un balayage vaut nozzle_width / density :
    1.0 = lignes jointives (une largeur de buse), 0.2 = une ligne toutes les
    cinq largeurs de buse.

    Correction (2026-09-26) : l'espacement etait derive de la TAILLE de la
    couche (`diagonal / density * 0.05`), jamais de la buse - si bien qu'une
    couche de 35 mm a 20 % recevait une ligne tous les 8,75 mm au lieu de
    2 mm, et que plus une couche etait grande, plus elle etait creuse.
    Mesure sur les couches Bird du Benchy (43 mm de remplissage au lieu
    d'environ 160 sur une couche de 424 mm2). Et `pattern` etait recu sans
    jamais etre lu : GRID produisait des lignes simples.

    sweep_angle_rad force l'orientation du balayage (ponts : traverser
    perpendiculairement a la portee, cf. Decision 2 du plan ; couches
    normales en CROSSED_ZIGZAG : 0 ou pi/2 selon la couche).
    """
    if boundary is None or boundary.is_empty or density <= 0:
        return []

    if pattern is InfillPattern.GRID:
        half = min(density, 1.0) / 2.0
        return (
            _generate_infill(boundary, half, InfillPattern.LINES, nozzle_width, sweep_angle_rad)
            + _generate_infill(
                boundary, half, InfillPattern.LINES, nozzle_width, sweep_angle_rad + np.pi / 2.0
            )
        )

    if pattern is InfillPattern.CROSSED_ZIGZAG:
        # `boundary` est le polygone de la paroi la plus interieure : son bord
        # EST l'axe de cette paroi. Les raccords du zigzag suivant ce bord, ils
        # repasseraient sur la paroi (double depot) : on les en ecarte d'une
        # demi-largeur de buse, cordons jointifs. Une zone plus etroite
        # qu'une buse ne recoit pas de remplissage (les parois s'y touchent).
        zone = boundary.buffer(-nozzle_width / 2.0, join_style="mitre")
        if zone.is_empty:
            return []
        scans, spacing = _clipped_scanlines(zone, density, nozzle_width, sweep_angle_rad)
        return _connect_zigzag(scans, zone, spacing)

    scans, spacing = _clipped_scanlines(boundary, density, nozzle_width, sweep_angle_rad)
    return [segment for scan in scans for segment in scan]


def _clipped_scanlines(
    boundary: Polygon, density: float, nozzle_width: float, sweep_angle_rad: float
) -> tuple[list[list[np.ndarray]], float]:
    """
    Lignes de balayage clippees par `boundary`, regroupees PAR ligne de
    balayage (dans l'ordre de balayage) et, au sein d'une ligne, triees et
    orientees dans le sens du balayage. Rend aussi l'espacement retenu.
    """
    minx, miny, maxx, maxy = boundary.bounds
    diagonal = float(np.hypot(maxx - minx, maxy - miny))
    if diagonal <= 0:
        return [], 0.0

    spacing = nozzle_width / min(density, 1.0)
    # borne basse pour eviter un nombre de lignes explosif (buse ou densite
    # aberrante) : le balayage couvre 2 x diagonal, donc au plus
    # _MAX_INFILL_LINES lignes par couche.
    spacing = max(spacing, 2.0 * diagonal / _MAX_INFILL_LINES, 1e-6)

    center = np.array([(minx + maxx) / 2.0, (miny + maxy) / 2.0])
    half_extent = diagonal  # sur-couvre largement pour rester valide apres rotation

    offsets = np.arange(-half_extent, half_extent, spacing)
    cos_a, sin_a = np.cos(sweep_angle_rad), np.sin(sweep_angle_rad)
    direction = np.array([cos_a, sin_a])
    perpendicular = np.array([-sin_a, cos_a])

    scans: list[list[np.ndarray]] = []
    for offset in offsets:
        base = center + perpendicular * offset
        p0 = base - direction * half_extent
        p1 = base + direction * half_extent
        clipped = LineString([p0, p1]).intersection(boundary)
        if clipped.is_empty:
            continue
        parts = clipped.geoms if hasattr(clipped, "geoms") else [clipped]
        segments = []
        for part in parts:
            if part.geom_type != "LineString":
                continue
            coords = np.asarray(part.coords)
            if len(coords) < 2:
                continue
            if np.dot(coords[-1] - coords[0], direction) < 0:
                coords = coords[::-1]
            segments.append(coords)
        if segments:
            segments.sort(key=lambda c: float(np.dot(c[0], direction)))
            scans.append(segments)
    return scans, spacing


def _connect_zigzag(
    scans: list[list[np.ndarray]], boundary: Polygon, spacing: float
) -> list[np.ndarray]:
    """
    Relie les segments de lignes de balayage SUCCESSIVES en zigzag : depuis
    la fin d'un segment, on enchaine sur l'extremite la plus proche d'un
    segment libre de la ligne suivante, que l'on parcourt alors en sens
    inverse. Le raccord SUIT LE BORD de la zone entre les deux extremites
    (`_boundary_connector`), comme dans un slicer classique ; a defaut (les
    deux extremites ne sont pas sur un meme bord), un raccord droit est
    accepte s'il reste dans `boundary`. Dans les deux cas il mesure au plus
    `_ZIGZAG_MAX_CONNECTOR_SPACINGS` espacements - au-dela il longerait un
    bord presque parallele au balayage et deposerait une ligne en trop.
    Sinon le zigzag s'arrete et un nouveau commence (une zone concave ou
    trouee donne donc plusieurs zigzags, relies par des trajets a vide).

    Un raccord droit seul ne suffit pas : mesure sur les couches Bird du
    Benchy, dont les bords sont des polygones irreguliers issus du maillage,
    le raccord droit sortait legerement de la zone a presque chaque ligne et
    la mediane tombait a UNE ligne par zigzag.
    """
    from shapely.prepared import prep

    inside = prep(boundary.buffer(max(1e-6, spacing * 1e-4)))
    max_connector = _ZIGZAG_MAX_CONNECTOR_SPACINGS * spacing
    rings = [
        LineString(ring.coords)
        for piece in _polygon_pieces(boundary)
        for ring in (piece.exterior, *piece.interiors)
    ]
    used = [[False] * len(scan) for scan in scans]
    chains: list[np.ndarray] = []

    for k, scan in enumerate(scans):
        for j, segment in enumerate(scan):
            if used[k][j]:
                continue
            used[k][j] = True
            parts = [segment]
            end = segment[-1]
            current = k
            while current + 1 < len(scans):
                best = None
                for jj, candidate in enumerate(scans[current + 1]):
                    if used[current + 1][jj]:
                        continue
                    for reverse, start in ((False, candidate[0]), (True, candidate[-1])):
                        distance = float(np.linalg.norm(start - end))
                        if best is None or distance < best[0]:
                            best = (distance, jj, reverse)
                if best is None or best[0] > max_connector:
                    break
                distance, jj, reverse = best
                candidate = scans[current + 1][jj]
                start = candidate[-1] if reverse else candidate[0]
                connector = _boundary_connector(end, start, rings, max_connector)
                if connector is None:
                    if distance > 0 and not inside.covers(LineString([end, start])):
                        break
                    connector = np.empty((0, 2))
                used[current + 1][jj] = True
                oriented = candidate[::-1] if reverse else candidate
                parts.append(connector)
                parts.append(oriented)
                end = oriented[-1]
                current += 1
            chains.append(np.vstack(parts))
    return chains


def _boundary_connector(
    a: np.ndarray, b: np.ndarray, rings: list[LineString], max_length: float
) -> np.ndarray | None:
    """
    Points intermediaires du plus court chemin LE LONG D'UN BORD (anneau
    exterieur ou trou) entre deux points situes sur ce bord, extremites
    exclues. None si aucun anneau ne porte les deux points, ou si le chemin
    depasse `max_length`.
    """
    from shapely.geometry import Point
    from shapely.ops import substring

    pa, pb = Point(a), Point(b)
    tolerance = 1e-6 * max(1.0, max_length)
    for ring in rings:
        if ring.distance(pa) > tolerance or ring.distance(pb) > tolerance:
            continue
        length = ring.length
        sa, sb = ring.project(pa), ring.project(pb)
        forward = (sb - sa) % length
        backward = (sa - sb) % length
        if min(forward, backward) > max_length:
            return None

        def along(s0, s1):
            if s1 >= s0:
                return np.asarray(substring(ring, s0, s1).coords)
            return np.vstack([np.asarray(substring(ring, s0, length).coords),
                              np.asarray(substring(ring, 0.0, s1).coords)])

        path = along(sa, sb) if forward <= backward else along(sb, sa)[::-1]
        return path[1:-1] if len(path) > 2 else np.empty((0, 2))
    return None


# ---------------------------------------------------------------------------
# OSCILLATING_SPIRAL
# ---------------------------------------------------------------------------

_SPIRAL_SAMPLE_STEP_FACTOR = 0.5       # pas d'echantillonnage = 0,5 buse
_SPIRAL_MIN_RADIUS_PITCHES = 0.5       # rayon de depart : pres de l'axe, C tourne tres vite
_SPIRAL_DENSITY_TOLERANCE = 0.05
_SPIRAL_DENSITY_ITERATIONS = 6
_SPIRAL_STOP_PROBABILITY = 0.3         # part (x randomisation) de pics purement radiaux


def _layer_rank(layer: Layer, isovalue_rank: dict) -> int:
    if layer.isovalue is not None and layer.isovalue in isovalue_rank:
        return int(isovalue_rank[layer.isovalue])
    return int(layer.id)


def _plate_axis_in_local_2d(
    local_frame: tuple[np.ndarray, np.ndarray, np.ndarray], axis_origin: tuple[float, float, float]
) -> np.ndarray:
    """
    Intersection de l'axe du plateau (vertical, passant par `axis_origin`) avec
    le plan de la couche, en coordonnees locales (u, v). Couche quasi
    verticale (pas d'intersection stable) : projection du point lui-meme.
    """
    origin, axis_u, axis_v = (np.asarray(a, dtype=float) for a in local_frame)
    normal = np.cross(axis_u, axis_v)
    normal /= max(float(np.linalg.norm(normal)), 1e-12)
    a0 = np.asarray(axis_origin, dtype=float)
    denom = float(normal[2])
    point = a0 if abs(denom) < 0.1 else a0 + np.array([0.0, 0.0, float(normal @ (origin - a0)) / denom])
    d = point - origin
    return np.array([float(d @ axis_u), float(d @ axis_v)])


def _oscillating_spiral_polyline(
    center: np.ndarray, r_min_zone: float, r_max: float, pitch: float, amplitude: float,
    wavelength: float, randomization: float, nozzle_width: float, rng: np.random.Generator,
    first_sign: float, angle_offset: float,
) -> np.ndarray:
    """
    Spirale r = pitch * theta / 2pi, oscillante, en coordonnees polaires
    (theta, r) autour de `center` : une suite de points de controle dont theta
    ne decroit JAMAIS (deux points au meme theta = pic qui n'avance pas,
    mouvement purement radial). Entre deux points, la deviation radiale suit
    une droite (pic net) ou un cosinus (courbe), toujours monotone en theta.
    `angle_offset` pivote toute la spirale (pi = decalage d'un demi-pas).
    """
    rho = float(np.clip(randomization, 0.0, 1.0))
    cap = max(0.0, 0.5 * (pitch - nozzle_width))   # deux tours voisins restent a >= 1 buse
    amp = min(amplitude, cap)
    step = nozzle_width * _SPIRAL_SAMPLE_STEP_FACTOR
    r_start = max(_SPIRAL_MIN_RADIUS_PITCHES * pitch, r_min_zone - 2.0 * amp)
    theta = 2.0 * np.pi * r_start / pitch
    theta_end = 2.0 * np.pi * (r_max + 2.0 * amp) / pitch
    base = lambda th: pitch * th / (2.0 * np.pi)

    ctrl_theta = [theta]
    ctrl_dev = [0.0]
    ctrl_curve = [False]
    sign = first_sign
    min_half = nozzle_width
    k = 0
    prev_stop = False
    while theta < theta_end and amp > 0.0:
        r_eff = max(base(theta), r_start)
        half = max(min_half, 0.5 * wavelength * (1.0 + rho * 0.6 * rng.uniform(-1.0, 1.0)))
        # jamais deux pics radiaux de suite : la tete repasserait sur la meme ligne
        stop = k > 0 and not prev_stop and rng.random() < _SPIRAL_STOP_PROBABILITY * rho
        prev_stop = stop
        a_k = float(np.clip(amp * (1.0 + rho * rng.uniform(-0.6, 0.0)), 0.25 * amp, amp))
        curve = bool(rho > 0.0 and rng.random() < 0.5)
        theta = theta if stop else theta + half / r_eff
        ctrl_theta.append(theta)
        ctrl_dev.append(sign * a_k)
        ctrl_curve.append(curve)
        sign = -sign
        k += 1
    if amp <= 0.0:
        ctrl_theta.append(theta_end)
        ctrl_dev.append(0.0)
        ctrl_curve.append(False)

    thetas: list[np.ndarray] = [np.array([ctrl_theta[0]])]
    devs: list[np.ndarray] = [np.array([ctrl_dev[0]])]
    for i in range(1, len(ctrl_theta)):
        t0, t1 = ctrl_theta[i - 1], ctrl_theta[i]
        d0, d1 = ctrl_dev[i - 1], ctrl_dev[i]
        r_eff = max(base(t1), r_start)
        length = max((t1 - t0) * r_eff, abs(d1 - d0), 1e-9)
        n = max(1, int(np.ceil(length / step)))
        s = np.linspace(0.0, 1.0, n + 1)[1:]
        thetas.append(t0 + (t1 - t0) * s)
        ease = (1.0 - np.cos(np.pi * s)) / 2.0 if ctrl_curve[i] else s
        devs.append(d0 + (d1 - d0) * ease)
    th = np.concatenate(thetas)
    r = base(th) + np.concatenate(devs)
    ang = th + angle_offset
    return np.column_stack([center[0] + r * np.cos(ang), center[1] + r * np.sin(ang)])


def _angle_nondecreasing(points: np.ndarray, center: np.ndarray, tolerance: float = 1e-9) -> bool:
    """Angle polaire autour de `center` (deroule) jamais decroissant le long de `points`."""
    d = points - center
    if np.any(np.hypot(d[:, 0], d[:, 1]) < 1e-9):
        return True            # passage par l'axe : angle indetermine, rien a verifier
    angles = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
    return bool(np.all(np.diff(angles) >= -tolerance))


def _generate_spiral_infill(
    boundary: Polygon, params: "Step6Params", rank: int, axis_local: np.ndarray | None
) -> list[np.ndarray]:
    """
    Remplissage OSCILLATING_SPIRAL d'une zone. Toutes les couches suivent la
    MEME spirale mediane (elles se superposent) ; seules les asperites
    different : oscillation en opposition de phase une couche sur deux et
    tirages differents (graine + rang). Le pas est corrige
    d'apres la longueur REELLEMENT deposee pour tenir `infill_density`.
    Les raccords entre morceaux suivent le bord (comme le zigzag) et ne sont
    pas garantis monotones en theta.
    """
    from shapely.prepared import prep

    if boundary is None or boundary.is_empty or params.infill_density <= 0:
        return []
    width = params.nozzle_width
    zone = boundary.buffer(-width / 2.0, join_style="mitre")
    if zone.is_empty or zone.area <= 0:
        return []
    density = min(params.infill_density, 1.0)
    coords_all = np.vstack([np.asarray(piece.exterior.coords) for piece in _polygon_pieces(zone)])
    center = np.asarray(axis_local, dtype=float) if axis_local is not None else np.array(zone.centroid.coords[0])
    radii = np.linalg.norm(coords_all - center, axis=1)
    r_max = float(radii.max())
    r_min_zone = 0.0 if zone.contains(shapely.Point(center)) else float(zone.distance(shapely.Point(center)))
    minx, miny, maxx, maxy = zone.bounds
    diagonal = float(np.hypot(maxx - minx, maxy - miny))
    min_pitch = max(2.0 * diagonal / _MAX_INFILL_LINES, 1e-6)

    pitch = max(width / density, min_pitch)
    odd = rank % 2 == 1
    best: LineString | None = None
    # amplitude et longueur d'onde fixees UNE fois (a partir du pas nominal) :
    # sinon le pas corrige les ferait varier et la correction de densite ne
    # convergerait pas.
    amplitude = (
        params.spiral_alternation_threshold
        if params.spiral_alternation_threshold is not None
        else 0.6 * max(0.0, 0.5 * (pitch - width))
    )
    wavelength = params.spiral_peak_spacing if params.spiral_peak_spacing is not None else 4.0 * width
    for _ in range(_SPIRAL_DENSITY_ITERATIONS):
        rng = np.random.default_rng([int(params.spiral_seed) & 0xFFFFFFFF, int(rank)])
        line = _oscillating_spiral_polyline(
            center, r_min_zone, r_max, pitch, amplitude, wavelength,
            params.spiral_randomization_threshold, width, rng,
            first_sign=-1.0 if odd else 1.0, angle_offset=0.0,
        )
        clipped = LineString(line).intersection(zone)
        length = float(clipped.length)
        best = clipped
        achieved = length * width / zone.area
        ratio = achieved / density if density > 0 else 1.0
        if abs(ratio - 1.0) <= _SPIRAL_DENSITY_TOLERANCE or length <= 0:
            break
        pitch = max(pitch * float(np.clip(ratio, 0.5, 2.0)), min_pitch)
    if best is None or best.is_empty:
        return []

    pieces = [np.asarray(g.coords) for g in (best.geoms if hasattr(best, "geoms") else [best])
              if g.geom_type == "LineString" and len(g.coords) >= 2]
    if not pieces:
        return []
    # ordre le long de la spirale : position du premier point sur la ligne
    line_geom = LineString(line)
    order = shapely.line_locate_point(line_geom, shapely.points(np.array([p[0] for p in pieces])))
    pieces = [pieces[i] for i in np.argsort(order, kind="stable")]

    inside = prep(zone.buffer(max(1e-6, pitch * 1e-4)))
    max_connector = _ZIGZAG_MAX_CONNECTOR_SPACINGS * pitch
    rings = [LineString(ring.coords) for piece in _polygon_pieces(zone) for ring in (piece.exterior, *piece.interiors)]
    chains: list[np.ndarray] = []
    current = pieces[0]
    for nxt in pieces[1:]:
        end, start = current[-1], nxt[0]
        gap = float(np.linalg.norm(start - end))
        connector = _boundary_connector(end, start, rings, max_connector) if gap <= max_connector else None
        if connector is None and gap > 0 and (gap > max_connector or not inside.covers(LineString([end, start]))):
            chains.append(current)
            current = nxt
            continue
        joint = np.vstack([end[None, :], connector if connector is not None else np.empty((0, 2)), start[None, :]])
        if not _angle_nondecreasing(joint, center):
            # un raccord qui ferait reculer le plateau : trajet a vide a la place
            chains.append(current)
            current = nxt
            continue
        parts = [current] + ([connector] if connector is not None and len(connector) else []) + [nxt]
        current = np.vstack(parts)
    chains.append(current)
    return chains



def _infill_angle_for_layer(layer: Layer, isovalue_rank: dict) -> float:
    """
    Orientation du balayage d'une couche normale (hors pont) en
    CROSSED_ZIGZAG : axe u du repere local (0) une couche sur deux, axe v
    (pi/2) l'autre. La parite est prise sur le RANG de l'isovaleur, pas sur
    la position dans la liste : l'ordre d'impression (step5, en arbre) fait
    se succeder des couches de branches differentes, alors que deux couches
    EMPILEES ont des isovaleurs consecutives - ce sont elles qui doivent se
    croiser. Les morceaux d'une meme isovaleur (clusters) partagent donc la
    meme orientation. Sans isovaleur, repli sur l'identifiant de couche.
    """
    if layer.isovalue is not None and layer.isovalue in isovalue_rank:
        rank = isovalue_rank[layer.isovalue]
    else:
        rank = layer.id
    return 0.0 if rank % 2 == 0 else np.pi / 2.0


# ---------------------------------------------------------------------------
# Decision 2 - Pont surface (is_bridge) <-> boucle 2D
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

# Nombre de fils de step6 (None = min(8, coeurs)) ; 1 = serie historique.
_STEP6_THREADS = None

# Pas 15 (limite #56) : parois et remplissage d'une zone par le noyau natif
# (Clipper2 au lieu de shapely/GEOS). PAS au bit pres - meme offset en
# onglet, autre arithmetique et autre point de depart d'anneau -, d'ou
# l'option. ACTIVE PAR DEFAUT quand le noyau natif est actif (2026-10-26) :
# APEX_NATIVE_FILL=0 rend shapely. Mesure sur les 4 modeles (limite #82).
# Les ponts restent en Python.
_NATIVE_FILL = os.environ.get("APEX_NATIVE_FILL", "1") == "1"
_NATIVE_PATTERN_CODES = {InfillPattern.LINES: 0, InfillPattern.GRID: 1, InfillPattern.CROSSED_ZIGZAG: 2}


def _native_fill_module():
    if not _NATIVE_FILL:
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

    def one_layer(layer: Layer) -> Layer:
        return _fill_layer(layer, mesh_data, overhang_map, params, isovalue_rank)

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
# APEX_NATIVE_LAYER=0 rend le chemin couche par couche de Python.
_NATIVE_LAYER = os.environ.get("APEX_NATIVE_LAYER", "1") == "1"


def _fill_layers_native(layers, mesh_data, overhang_map, params, isovalue_rank):
    """Liste de Layer avec `toolpath_2d` rempli, ou None si le chemin natif ne s'applique pas
    (noyau inactif, motif non natif, kd-tree sans version native)."""
    module = _native_fill_module()
    if (not _NATIVE_LAYER or module is None or not hasattr(module, "fill_layers")
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
    result = module.fill_layers(
        frames, loop_begin, point_begin, points, sweep, config, kdtree,
        np.ascontiguousarray(overhang_map.is_bridge, dtype=np.uint8),
        np.ascontiguousarray(overhang_map.bridge_span_direction, dtype=float), int(threads))
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
            elif _native_fill_module() is not None and params.infill_pattern in _NATIVE_PATTERN_CODES:
                sweep_angle = (
                    _infill_angle_for_layer(layer, isovalue_rank)
                    if params.infill_pattern is InfillPattern.CROSSED_ZIGZAG
                    else 0.0
                )
                valid, infill, offset_walls = _native_fill_module().fill_region(
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


def _loop_matches(loop_2d: np.ndarray, outer_open: np.ndarray) -> bool:
    """
    Rapproche une boucle projetee (avec point de fermeture duplique) du
    contour "matiere" ouvert issu de _classify_contours (meme boucle,
    representation differente) - comparaison sur le premier point, unique
    par construction (les boucles d'une meme couche ne se croisent pas).
    """
    open_loop = loop_2d[:-1] if np.allclose(loop_2d[0], loop_2d[-1]) else loop_2d
    return open_loop.shape == outer_open.shape and np.allclose(open_loop, outer_open)
