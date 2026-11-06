"""Etape 6 : repère local 2D, décalage hors plan, classification des contours et parois.

Découpé de `step6_toolpath2d.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step6_toolpath2d`.
"""

from __future__ import annotations

import numpy as np
from shapely.geometry import MultiPolygon, Polygon
from core import _cpp
from core.geometry_utils import point_in_polygon_2d


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
    les couches du Benchy :
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


def _loop_matches(loop_2d: np.ndarray, outer_open: np.ndarray) -> bool:
    """
    Rapproche une boucle projetee (avec point de fermeture duplique) du
    contour "matiere" ouvert issu de _classify_contours (meme boucle,
    representation differente) - comparaison sur le premier point, unique
    par construction (les boucles d'une meme couche ne se croisent pas).
    """
    open_loop = loop_2d[:-1] if np.allclose(loop_2d[0], loop_2d[-1]) else loop_2d
    return open_loop.shape == outer_open.shape and np.allclose(open_loop, outer_open)
