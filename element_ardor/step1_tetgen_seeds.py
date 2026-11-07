"""Etape 1 : points intérieurs, cavités et découpe du raffinement près du plateau (TetGen).

Découpé de `step1_ingest.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step1_ingest`.
"""

from __future__ import annotations

import numpy as np
import trimesh
import trimesh.transformations as trimesh_tf
from trimesh import geometry as trimesh_geometry
from trimesh import grouping as trimesh_grouping
from trimesh.creation import triangulate_polygon
from trimesh.path import polygons as trimesh_polygons
from core.mesh_types import JunctionZone
from element_ardor.step1_types import (  # noqa: F401
    MeshIngestError,
)


_INTERIOR_POINT_SEED = 0  # deterministe - reproductibilite des seeds de region


_INTERIOR_POINT_INITIAL_SAMPLES = 64


_INTERIOR_POINT_MAX_SAMPLES = 4096  # plafond - au-dela, echec explicite


def _find_guaranteed_interior_point(
    shell: trimesh.Trimesh, fallback_point: np.ndarray, context_label: str
) -> np.ndarray:
    """
    Retourne un point GARANTI a l'interieur de `shell` (verifie via
    `shell.contains()`), pour toute coque watertight arbitraire (convexe ou
    non). Essaie d'abord `fallback_point` (typiquement un centroide deja
    calcule par l'appelant, chemin rapide - reste correct sur la tres grande
    majorite des coques, notamment tout ce qui est convexe ou proche de
    convexe). Si `fallback_point` echoue le test d'inclusion (mesure
    concrete : forme non convexe ou coque a plusieurs composantes
    disjointes, cf. CLAUDE.md limite #42 - un centroide de coque "rest" hors
    du solide reel sur un modele en T, qui casse SILENCIEUSEMENT
    `TetGen.add_region` : la region dont le seed est hors-coque n'est plus
    localisable dans la triangulation, et son volume cible est ignore au
    profit de l'autre region, appliquee alors a tout le maillage), retombe
    sur un tirage par rejet (`trimesh.sample.volume_mesh`, seed FIXE pour la
    reproductibilite - meme discipline que le determinisme deja exige du tri
    DFS spatial de step5), avec une escalade deterministe du nombre
    d'essais si le tirage initial ne produit aucun point accepte (mesure :
    `volume_mesh` peut renvoyer 0 point sur une coque fine par rapport a sa
    bounding box).

    Retourne TOUJOURS un seul point INDIVIDUELLEMENT valide, jamais une
    moyenne de plusieurs points acceptes - piege mesure concretement : la
    moyenne de plusieurs points valides peut elle-meme retomber hors d'une
    coque non convexe, meme si chacun des points moyennes est
    individuellement a l'interieur.

    Args:
        shell: coque watertight (suppose deja etanche cote appelant, meme
            convention que le reste de ce fichier a ce stade du pipeline).
        fallback_point: point candidat a essayer en premier (ex: centroide).
        context_label: texte court identifiant l'appelant, uniquement pour
            un message d'erreur lisible (ex: "coque 'rest'").

    Returns:
        Point (3,) verifie a l'interieur de `shell`.

    Raises:
        MeshIngestError: coque vide/degeneree (aucune face, bounding box
            indefinie - sans ce garde-fou, `trimesh.sample.volume_mesh`
            plante avec un `TypeError` cryptique, mesure concretement) avant
            meme d'essayer l'echantillonnage, ou aucun point trouve apres
            escalade jusqu'a `_INTERIOR_POINT_MAX_SAMPLES` essais (coque
            quasi degeneree, ratio volume/bounding-box extremement faible).
    """
    if shell.faces.shape[0] == 0 or shell.bounds is None:
        raise MeshIngestError(
            f"{context_label} : coque vide ou degeneree, impossible d'y "
            "chercher un point interieur."
        )

    if shell.contains([fallback_point])[0]:
        return np.asarray(fallback_point, dtype=np.float64)

    sample_count = _INTERIOR_POINT_INITIAL_SAMPLES
    while sample_count <= _INTERIOR_POINT_MAX_SAMPLES:
        points = trimesh.sample.volume_mesh(shell, count=sample_count, seed=_INTERIOR_POINT_SEED)
        if len(points) > 0:
            return np.asarray(points[0], dtype=np.float64)
        sample_count *= 8

    raise MeshIngestError(
        f"{context_label} : impossible de trouver un point garanti a "
        f"l'interieur apres {_INTERIOR_POINT_MAX_SAMPLES} essais "
        "d'echantillonnage volumique (coque quasi degeneree - ratio "
        "volume/bounding-box extremement faible)."
    )


def _find_internal_cavities(mesh: trimesh.Trimesh) -> list[np.ndarray]:
    """
    Detecte les cavites internes du mesh (ex: boite creuse) et retourne un
    point representatif a l'interieur de chacune, a passer a
    tetgen.TetGen.add_hole(). Sans ça, tetgen ne sait pas qu'un volume
    interne doit rester vide et peut generer des tetraedres qui traversent
    la cavite (mesure : jusqu'a ~6% des tetraedres sur une boite creuse
    simple, centroides verifies hors du solide reel).

    Methode : separe le mesh en coques connexes etanches
    (`mesh.split(only_watertight=True)`). Le signe du volume d'une coque
    (propriete trimesh, deja orientee) distingue frontiere solide (positif)
    de cavite (negatif - ses normales, heritees d'une operation booleenne
    de difference, pointent vers l'interieur du vide). Le centroide de
    chaque coque-cavite sert de point representatif ; si ce centroide tombe
    hors de sa propre coque (cavite non convexe, ex: forme en tore), repli
    via `_find_guaranteed_interior_point` (rejection sampling deterministe)
    plutot qu'un echec immediat.
    """
    shells = mesh.split(only_watertight=True)
    if len(shells) <= 1:
        return []

    hole_points = []
    for shell in shells:
        if shell.volume >= 0:
            continue  # coque solide (frontiere exterieure ou paroi normale), pas une cavite

        centroid = shell.vertices.mean(axis=0)
        hole_points.append(
            _find_guaranteed_interior_point(shell, centroid, "cavite interne (volume signe negatif)")
        )

    return hole_points


def _shared_cut_cap(
    open_mesh: trimesh.Trimesh, plane_origin: np.ndarray, plane_normal: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Triangule UNE SEULE FOIS le contour de bord cree par une coupe
    `slice_plane(..., cap=False)`, pour etre reutilisee identique des deux
    cotes de la coupe - voir `_split_for_platform_refinement` pour pourquoi
    deux triangulations independantes (une par `slice_plane(cap=True)`)
    divergent sur un contour non convexe.

    Reproduit la logique interne du parametre `cap=True` de
    `trimesh.intersections.slice_mesh_plane` (dedup des sommets, extraction
    des aretes de bord sur le plan, `triangulate_polygon`), mais appliquee
    UNE FOIS a la place de deux fois independamment.

    Retourne (cap_vertices (K,3), cap_faces (T,3)) dans le repere du monde,
    triangulation orientee de facon a ce que la normale des triangles suive
    `plane_normal` (meme convention que `slice_plane(cap=True)` : le cap
    ferme le mesh du cote NEGATIF de `plane_normal`, cote qui a ete tranche).
    """
    vertices = open_mesh.vertices
    faces = open_mesh.faces

    # slice_plane(cap=False) laisse des sommets numeriquement dupliques a
    # chaque intersection arete/plan, sans les fusionner - meme etape de
    # dedup que le code interne de slice_mesh_plane avant de calculer un cap.
    unique, inverse = trimesh_grouping.unique_rows(vertices)
    vertices = vertices[unique]
    faces = inverse[faces]
    faces = faces[(faces[:, :1] != faces[:, 1:]).all(axis=1)]

    dots = np.dot(vertices - plane_origin, plane_normal)
    on_plane = np.abs(dots) < 1e-8

    edges = trimesh_geometry.faces_to_edges(faces)
    edges.sort(axis=1)
    edges = edges[on_plane[edges].all(axis=1)]
    edges = edges[edges[:, 0] != edges[:, 1]]

    # arete de bord du contour de coupe = apparait exactement une fois
    # (une arete interieure a une face qui serait par coincidence entierement
    # sur le plan apparaitrait deux fois, une par face adjacente).
    unique_edge_idx = trimesh_grouping.group_rows(edges, require_count=1)
    edges_boundary = edges[unique_edge_idx]

    to_2D = trimesh_geometry.plane_transform(origin=plane_origin, normal=-plane_normal)
    to_3D = np.linalg.inv(to_2D)
    vertices_2D = trimesh_tf.transform_points(vertices, to_2D)

    all_v2 = []
    all_f = []
    offset = 0
    for poly in trimesh_polygons.edges_to_polygons(edges_boundary, vertices_2D[:, :2]):
        v2, f2 = triangulate_polygon(poly, force_vertices=True)
        all_v2.append(v2)
        all_f.append(f2 + offset)
        offset += len(v2)

    if not all_v2:
        return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64)

    v2_all = np.vstack(all_v2)
    f_all = np.vstack(all_f)
    v3_all = trimesh_tf.transform_points(trimesh.util.stack_3D(v2_all), to_3D)

    return np.asarray(v3_all, dtype=np.float64), np.asarray(f_all, dtype=np.int64)


def _split_for_platform_refinement(
    mesh: trimesh.Trimesh,
    near_platform_margin: float,
    platform_z: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Decoupe `mesh` en un PLC combine, deux regions de part et d'autre de `z = platform_z +
    near_platform_margin` ("slab" pres du plateau, "rest" au-dela), pour alimenter `TetGen.add_region` avec un
    seed par region (le mecanisme bgmesh/metric est ecarte, voir limite #35).

    Le contour de coupe est triangule UNE SEULE FOIS (`_shared_cut_cap`) et n'apparait qu'UNE fois dans le PLC,
    comme facette INTERNE partagee. Bug mesure sur death.stl : deux `slice_plane(cap=True)` independants
    triangulent le meme contour non convexe de facon differente (12 triangles divergents sur 47), ce que TetGen
    refuse ("self-intersections") ; le repli de `_tetrahedralize` (limite #43) masquait l'echec en desactivant
    le raffinement plateau (couches en « escalier »). Apres correction : death.stl 112698 tetraedres avec
    raffinement reel (107037 en mode degrade), l_shape_overhang.stl 1677 (inchange).

    Retourne (combined_vertices (N,3), combined_faces (M,3), slab_seed (3,), rest_seed (3,)). Chaque seed part
    du centroide de sa coque ouverte, repousse vers l'interieur de sa region, VERIFIE via `shell.contains()`
    (`_find_guaranteed_interior_point`) sur une coque watertight temporaire (jamais le PLC final) ; tirage par
    rejet deterministe sur forme non convexe (limite #42 : un seed hors-coque casse silencieusement
    `TetGen.add_region`).
    """
    cut_z = platform_z + near_platform_margin
    plane_origin = np.array([0.0, 0.0, cut_z])
    plane_normal_up = np.array([0.0, 0.0, 1.0])

    slab_open = mesh.slice_plane(plane_origin=plane_origin, plane_normal=-plane_normal_up, cap=False)
    rest_open = mesh.slice_plane(plane_origin=plane_origin, plane_normal=plane_normal_up, cap=False)

    cap_vertices, cap_faces = _shared_cut_cap(rest_open, plane_origin, plane_normal_up)

    n_slab = len(slab_open.vertices)
    n_rest = len(rest_open.vertices)

    combined_vertices_dup = np.vstack([slab_open.vertices, rest_open.vertices, cap_vertices])
    combined_faces_dup = np.vstack([
        slab_open.faces,
        rest_open.faces + n_slab,
        cap_faces + n_slab + n_rest,  # une seule copie - facette interne partagee
    ])

    # slab_open/rest_open portent chacun leur copie des sommets du contour de coupe : sans cette fusion
    # globale, le cap resterait deconnecte des bords ouverts malgre des coordonnees identiques.
    unique, inverse = trimesh_grouping.unique_rows(combined_vertices_dup)
    combined_vertices = combined_vertices_dup[unique]
    combined_faces = inverse[combined_faces_dup]
    combined_faces = combined_faces[(combined_faces[:, :1] != combined_faces[:, 1:]).all(axis=1)]

    # Coques watertight temporaires (coque ouverte + cap), uniquement pour le test d'inclusion des seeds.
    slab_for_seed = trimesh.Trimesh(
        vertices=np.vstack([slab_open.vertices, cap_vertices]),
        faces=np.vstack([slab_open.faces, cap_faces[:, ::-1] + n_slab]),
        process=True,
    )
    rest_for_seed = trimesh.Trimesh(
        vertices=np.vstack([rest_open.vertices, cap_vertices]),
        faces=np.vstack([rest_open.faces, cap_faces + n_rest]),
        process=True,
    )

    slab_seed_candidate = slab_open.vertices.mean(axis=0)
    slab_seed_candidate[2] = min(slab_seed_candidate[2], platform_z + near_platform_margin * 0.5)
    rest_seed_candidate = rest_open.vertices.mean(axis=0)
    rest_seed_candidate[2] = max(rest_seed_candidate[2], platform_z + near_platform_margin * 1.5)

    slab_seed = _find_guaranteed_interior_point(
        slab_for_seed, slab_seed_candidate, "coque 'slab' (pres du plateau)"
    )
    rest_seed = _find_guaranteed_interior_point(
        rest_for_seed, rest_seed_candidate, "coque 'rest' (reste du modele)"
    )

    return (
        np.asarray(combined_vertices, dtype=np.float64),
        np.asarray(combined_faces, dtype=np.int64),
        np.asarray(slab_seed, dtype=np.float64),
        np.asarray(rest_seed, dtype=np.float64),
    )


def _merged_junction_seed(mesh: trimesh.Trimesh, junction_zones: list[JunctionZone]) -> np.ndarray:
    """
    Fusionne toutes les `JunctionZone` detectees en UN SEUL seed interne
    (centroide pondere par rayon), pour alimenter une UNIQUE region
    `TetGen.add_region` supplementaire plutot que N regions separees.

    Necessaire car `TetGen.add_region` (binding 0.8.4) ne se generalise PAS
    correctement a N regions actives simultanement - mesure directe
    (scripts/measure_add_region_n_zones.py, session de la limite #48) : sur
    4 regions (1 "rest" + 3 zones fines), AUCUNE des 3 zones fines
    n'atteint le facteur de densite attendu face a "rest" (l'une d'elles,
    junction_y, se retrouve meme MOINS dense que la region censee etre
    grossiere) - une generalisation du bug deja documente a 2 regions
    (limite #42), pas juste une variante mineure. Le fallback a 2 regions
    (1 "rest" + 1 zone fine fusionnee) fonctionne correctement, mesure sur
    le meme modele (scripts/measure_add_region_merged_zone.py) : facteur
    ~48x de densite locale pour la zone fusionnee contre la zone "rest" et
    un point de controle hors de toute zone fine, aucune fuite du volume
    fin hors de son rayon d'influence.

    Consequence acceptee : deux jonctions multi-region tres eloignees
    l'une de l'autre partagent le meme volume cible fin, meme si l'une
    seule en aurait besoin - imprecision mineure comparee au fait de payer
    le cout d'un maillage fin partout (`resolution` globalement raffine),
    et strictement meilleure que le comportement actuel (aucun raffinement
    du tout aux jonctions, limite #47).

    Le point retourne n'est PAS garanti a l'interieur du mesh (une moyenne
    de centroides individuellement valides peut retomber hors d'un solide
    non convexe, piege deja documente pour `_find_guaranteed_interior_point`)
    - l'appelant doit le passer par cette fonction avant de l'utiliser comme
    seed `add_region`.
    """
    centroids = np.array([zone.centroid for zone in junction_zones])
    radii = np.array([zone.radius for zone in junction_zones])
    weights = radii / radii.sum()
    return np.average(centroids, axis=0, weights=weights)


def _find_rest_seed_outside_junction_zones(
    mesh: trimesh.Trimesh, junction_zones: list[JunctionZone]
) -> np.ndarray:
    """
    Seed "rest" garanti a la fois A L'INTERIEUR du mesh ET LE PLUS ELOIGNE
    POSSIBLE de toute `JunctionZone` - bug reel mesure en ecrivant les tests
    de cette fonction (pas suppose) : le centroide naif du mesh (`mesh.
    vertices.mean(axis=0)`) peut tomber A L'INTERIEUR du rayon d'une zone de
    jonction (mesure sur le tripode synthetique : centroide a 3.42mm du
    centre d'une zone de rayon 6.0mm), auquel cas TetGen applique le volume
    FIN de la region jonction a la region "rest" tout entiere - explosion
    de densite mesuree sur tout le maillage (905 -> 44985 tetraedres sur le
    tripode, loin de toute jonction : 34 -> 2302 tetraedres locaux), pas
    seulement pres de la jonction ciblee.

    Essaie d'abord le centroide naif (chemin rapide, correct des que le
    mesh est assez grand par rapport aux zones de jonction). Sinon, cherche
    parmi TOUS les sommets du mesh (deterministe, pas un tirage aleatoire)
    celui qui MAXIMISE la distance minimale a toutes les zones - garantit de
    trouver le meilleur candidat possible plutot que d'esperer un tirage
    aleatoire chanceux, cas mesure directement sur res/test_models/death.stl
    (motif festonne tres etendu, une zone fusionnee de rayon 47mm couvre une
    portion importante du volume) : un rejection-sampling par tirage
    aleatoire y echouait apres 4096 essais malgre l'existence de sommets
    valides (les coins/extremites du mesh, hors de la zone dense mais rares
    relativement au volume total echantillonne par `trimesh.sample.
    volume_mesh`).
    """
    candidate = mesh.vertices.mean(axis=0)

    def _min_dist_to_zones(point: np.ndarray) -> float:
        return min(
            (np.linalg.norm(point - zone.centroid) - zone.radius for zone in junction_zones),
            default=np.inf,
        )

    if mesh.contains([candidate])[0] and _min_dist_to_zones(candidate) > 0.0:
        return np.asarray(candidate, dtype=np.float64)

    # Les sommets du mesh sont sur la SURFACE, jamais utilisables tels
    # quels comme seed volumique TetGen - on cherche le sommet le plus
    # eloigne de toutes les zones, puis on le POUSSE vers le centroide
    # global du mesh (garanti interieur pour un solide simplement connexe,
    # meme discipline que le "candidate" naif ci-dessus) jusqu'a ce que le
    # point pousse soit lui-meme a l'interieur du mesh.
    vertices = mesh.vertices
    min_dists = np.array([
        min((np.linalg.norm(v - zone.centroid) - zone.radius for zone in junction_zones), default=np.inf)
        for v in vertices
    ])
    best_idx = int(np.argmax(min_dists))

    if min_dists[best_idx] <= 0.0:
        raise MeshIngestError(
            "Impossible de trouver un seed 'rest' hors de toutes les zones de "
            "jonction : meme le sommet du mesh le plus eloigne "
            f"(distance={min_dists[best_idx]:.2f}mm au bord de la zone la plus "
            "proche) reste a l'interieur d'une zone - les zones de jonction "
            "couvrent la quasi-totalite du volume sur ce modele."
        )

    surface_point = vertices[best_idx]
    for fraction in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5):
        pushed = surface_point + fraction * (candidate - surface_point)
        if mesh.contains([pushed])[0] and _min_dist_to_zones(pushed) > 0.0:
            return np.asarray(pushed, dtype=np.float64)

    raise MeshIngestError(
        "Impossible de pousser le seed 'rest' (sommet le plus eloigne des "
        "zones de jonction) a l'interieur du mesh sans re-tomber dans une "
        "zone de jonction - geometrie non convexe non couverte par cette "
        "heuristique."
    )
