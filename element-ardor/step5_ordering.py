"""Etape 5 : regroupement spatial des contours et ordre de parcours par continuité (DFS).

Découpé de `step5_layers.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step5_layers`.
"""

from __future__ import annotations

from typing import Optional
import networkx as nx
import numpy as np
from core.geometry_utils import point_in_polygon_2d


def _min_pairwise_distance(loop_a: np.ndarray, loop_b: np.ndarray) -> float:
    """Distance minimale point-a-point entre deux boucles 3D (enveloppes brutes, pas de projection)."""
    diffs = loop_a[:, None, :] - loop_b[None, :, :]
    return float(np.linalg.norm(diffs, axis=-1).min())


def _loops_connected(
    proj_a: np.ndarray,
    proj_b: np.ndarray,
    loop_a_3d: np.ndarray,
    loop_b_3d: np.ndarray,
    gap_threshold: float,
) -> bool:
    """
    Deux boucles sont "connectees" si elles sont soit IMBRIQUEES (une
    contient l'autre - ex: paroi exterieure + cavite interne, ne doit JAMAIS
    etre scinde, cf. CLAUDE.md), soit a une distance de bord-a-bord
    inferieure a `gap_threshold` (evite une discontinuite artificielle juste
    avant une fusion geometrique reelle a une isovaleur voisine).

    `proj_a`/`proj_b` : boucles projetees en 2D dans un repere de reference
    commun (utilise uniquement pour le test d'imbrication - la topologie
    relative de deux boucles franchement disjointes ne peut pas s'imbriquer
    accidentellement par un mauvais choix de plan, cf. CLAUDE.md).
    `loop_a_3d`/`loop_b_3d` : boucles brutes 3D (test de distance, pas de
    projection necessaire).

    Factorisee depuis `_cluster_contours` (clustering intra-isovaleur) pour
    reutilisation par `_build_continuity_graph` (continuite inter-isovaleur,
    Decision 8) - meme critere, une seule source de verite.
    """
    nested = (
        point_in_polygon_2d(proj_b[0], proj_a)
        or point_in_polygon_2d(proj_a[0], proj_b)
    )
    if nested:
        return True
    return _min_pairwise_distance(loop_a_3d, loop_b_3d) <= gap_threshold


def _cluster_contours(
    contours: list[np.ndarray],
    provisional_frame: tuple[np.ndarray, np.ndarray, np.ndarray],
    gap_threshold: float,
) -> list[list[int]]:
    """
    Partitionne `contours` (boucles fermees d'une meme couche) en clusters
    spatiaux, via `_loops_connected` (imbrication ou distance de bord-a-bord).

    `provisional_frame` sert UNIQUEMENT de reference de projection pour le
    test d'imbrication (calcule sur TOUTES les boucles de la couche, avant
    clustering - le meme plan "approximatif" qui motive ce fix en aval,
    mais les erreurs d'angle de ce plan moyen n'affectent pas la topologie
    relative de deux boucles tant que la projection ne degenere pas).

    Retourne une liste de clusters, chaque cluster etant une liste d'indices
    dans `contours` (triee, ordre stable).
    """
    n = len(contours)
    graph = nx.Graph()
    graph.add_nodes_from(range(n))

    origin, axis_u, axis_v = provisional_frame
    projected = []
    for c in contours:
        relative = c[:-1] - origin  # sans le point de fermeture duplique
        projected.append(np.stack([relative @ axis_u, relative @ axis_v], axis=1))

    for i in range(n):
        for j in range(i + 1, n):
            if _loops_connected(projected[i], projected[j], contours[i], contours[j], gap_threshold):
                graph.add_edge(i, j)

    return [sorted(component) for component in nx.connected_components(graph)]


def _tets_for_cluster(
    unique_tets: np.ndarray,
    tet_centers: np.ndarray,
    cluster_points: np.ndarray,
    tolerance: float,
) -> np.ndarray:
    """
    Sous-ensemble de `unique_tets` dont le centroide est proche (a
    `tolerance` pres) d'au moins un point du cluster - approximation de la
    correspondance triangle->tetraedre-source, PERDUE par
    `trimesh.outline()` (verifie : les entites de `outline.entities` ne
    portent aucune reference aux indices de `faces` d'entree). Un tetraedre
    ambigu (proche de plusieurs clusters) peut legitimement contribuer a
    plusieurs d'entre eux sans consequence : region_ids ne sert qu'a un
    rang topologique MAX pour le tri final (_sort_key), pas a une
    partition exacte.
    """
    if len(unique_tets) == 0:
        return unique_tets
    dists = np.linalg.norm(
        tet_centers[:, None, :] - cluster_points[None, :, :], axis=-1
    ).min(axis=1)
    return unique_tets[dists <= tolerance]


def _build_continuity_graph(
    raw_layers: list[tuple[float, list[np.ndarray], list[int], tuple]],
    gap_threshold: float,
    progress: Optional[np.ndarray] = None,
) -> nx.DiGraph:
    """
    Graphe oriente : un noeud par entree de `raw_layers` (indice dans cette
    liste), une arete i->j si le cluster i (isovaleur v_i) et le cluster j
    (isovaleur v_j) sont connectes (`_loops_connected`) ET v_i/v_j sont des
    isovaleurs distinctes ADJACENTES dans la sequence reellement presente
    dans `raw_layers` (pas la grille theorique de `_generate_isovalues` -
    une isovaleur sans contour exploitable est simplement absente et ne doit
    pas casser la continuite de ses voisines survivantes).

    Le seuil de connexion est `gap_threshold`, le MEME que celui utilise par
    `_cluster_contours` pour le clustering intra-isovaleur : deux isovaleurs
    successives sont physiquement proches en 3D (c'est la definition meme de
    "couches successives"), un cluster et son prolongement direct a
    l'isovaleur suivante doit etre au moins aussi proche que deux boucles de
    la meme couche qui ne fusionnent pas encore - pas une nouvelle constante
    arbitraire.

    La projection necessaire au test d'imbrication de `_loops_connected`
    reutilise le `local_frame` DEJA calcule par cluster dans `raw_layers`
    (contrairement a `_cluster_contours`, qui doit encore utiliser un plan
    provisoire commun avant que les clusters n'existent) : chaque cluster a
    deja son propre plan local coherent depuis le clustering. Projection
    faite dans le repere du cluster i (choix arbitraire mais coherent pour
    le test de distance ; le test d'imbrication teste deja les deux sens).

    `progress` (optionnel, une valeur par entree de `raw_layers`) remplace
    l'isovaleur pour ORIENTER les aretes, sans jamais toucher au critere
    d'ADJACENCE, qui reste l'isovaleur : "quelles couches se touchent" est une
    question de geometrie (deux isosurfaces successives sont voisines), tandis
    que "laquelle vient avant l'autre" est une question de progression
    d'impression. Confondre les deux est precisement ce qui fait echouer un
    tri global par distance geodesique seule (mesure : x4 de degradation,
    l'ordre entrelacant des couches spatialement eloignees qui n'ont en commun
    que leur distance au plateau).

    Avec `progress`, une arete peut donc etre INVERSEE par rapport a
    l'isovaleur : c'est exactement le cas d'une coque suspendue dont le champ
    croit de l'extremite vers l'attache, ou l'ordre correct est l'inverse de
    l'isovaleur. Une paire de progression strictement egale ne produit aucune
    arete (aucun des deux ne precede l'autre).

    Retourne un graphe acyclique par construction : les aretes vont toujours
    d'une valeur strictement inferieure vers une strictement superieure, que
    cette valeur soit l'isovaleur ou `progress`.
    """
    graph = nx.DiGraph()
    graph.add_nodes_from(range(len(raw_layers)))

    groups: dict[float, list[int]] = {}
    for idx, (isovalue, *_rest) in enumerate(raw_layers):
        groups.setdefault(round(isovalue, 6), []).append(idx)

    sorted_isovalues = sorted(groups.keys())

    for k in range(len(sorted_isovalues) - 1):
        indices_a = groups[sorted_isovalues[k]]
        indices_b = groups[sorted_isovalues[k + 1]]

        for i in indices_a:
            _, contours_i, _region_ids_i, frame_i = raw_layers[i]
            origin_i, axis_u_i, axis_v_i = frame_i
            points_i = np.concatenate(contours_i, axis=0)
            proj_i = np.stack(
                [
                    (points_i - origin_i) @ axis_u_i,
                    (points_i - origin_i) @ axis_v_i,
                ],
                axis=1,
            )

            for j in indices_b:
                _, contours_j, _region_ids_j, _frame_j = raw_layers[j]
                points_j = np.concatenate(contours_j, axis=0)
                proj_j = np.stack(
                    [
                        (points_j - origin_i) @ axis_u_i,
                        (points_j - origin_i) @ axis_v_i,
                    ],
                    axis=1,
                )
                if not _loops_connected(
                    proj_i, proj_j, points_i, points_j, gap_threshold
                ):
                    continue

                weight = _min_pairwise_distance(points_i, points_j)
                if progress is None:
                    graph.add_edge(i, j, weight=weight)
                elif progress[i] < progress[j]:
                    graph.add_edge(i, j, weight=weight)
                elif progress[j] < progress[i]:
                    # La progression d'impression s'oppose a l'isovaleur ici :
                    # l'arete est inversee, pas ignoree.
                    graph.add_edge(j, i, weight=weight)

    return graph


def _subtree_sizes(graph: nx.DiGraph, nodes: set[int]) -> dict[int, int]:
    """
    Taille du sous-arbre atteignable depuis chaque noeud (nombre de
    descendants, soi-meme inclus), restreint a `nodes`. Utilise
    `nx.descendants` (pas une somme naive des tailles d'enfants) car le
    graphe de continuite peut avoir des recombinaisons (in-degree >= 2) : un
    sous-arbre partage entre deux parents ne doit pas etre compte en double.
    Cout negligeable a l'echelle reelle de ces graphes (quelques dizaines de
    noeuds par couche).
    """
    subgraph = graph.subgraph(nodes)
    return {node: len(nx.descendants(subgraph, node)) + 1 for node in nodes}


def _dfs_order_within_rank(
    graph: nx.DiGraph,
    nodes: set[int],
    raw_layers: list[tuple[float, list[np.ndarray], list[int], tuple]],
    progress: Optional[np.ndarray] = None,
) -> list[int]:
    """
    Ordre de visite DFS des `nodes` (sous-ensemble de raw_layers appartenant
    a un meme rang topologique), respectant la continuite spatiale : une
    branche est entierement terminee avant de sauter a une autre.

    A un noeud de branchement, visite la branche de plus PETITE taille de
    sous-arbre en premier (`_subtree_sizes`) : termine les branches courtes
    d'abord, jamais de retour en arriere sur une grande branche restee en
    suspens - robuste a un nombre arbitraire d'enfants, contrairement a un
    critere de proximite au point de sortie du parent qui degenere en
    quasi-egalite au moment meme du split (cf. CLAUDE.md).

    Entre composantes faiblement connexes independantes d'un meme rang,
    visite par isovaleur minimale croissante (generalisation directe du tri
    historique par isovaleur, applique a l'unite "composante" plutot qu'a un
    cluster individuel).

    `progress` (optionnel) remplace l'isovaleur comme mesure de PROGRESSION
    d'impression - ordre des composantes et choix des racines. La continuite
    spatiale (structure du graphe, DFS, priorite aux petites branches) est
    strictement inchangee : ce parametre ne dit pas "ou aller", seulement "par
    quel bout commencer".

    Tri explicite partout (tie-break par indice de noeud) : deux appels
    identiques produisent toujours exactement le meme ordre.
    """
    def key_of(node: int) -> float:
        return float(progress[node]) if progress is not None else raw_layers[node][0]

    subgraph = graph.subgraph(nodes)
    components = sorted(
        nx.weakly_connected_components(subgraph),
        key=lambda comp: (min(key_of(n) for n in comp), min(comp)),
    )

    sizes = _subtree_sizes(graph, nodes)
    visited: set[int] = set()
    order: list[int] = []

    def visit(node: int) -> None:
        if node in visited:
            return
        visited.add(node)
        order.append(node)
        children = sorted(
            (c for c in subgraph.successors(node) if c not in visited),
            key=lambda c: (sizes[c], c),
        )
        for child in children:
            visit(child)

    for component in components:
        roots = sorted(
            (n for n in component if subgraph.in_degree(n) == 0),
            key=lambda n: (key_of(n), n),
        )
        if not roots:
            # garde-fou defensif (un DAG fini a toujours une source, ne
            # devrait pas arriver) : repli sur le noeud de plus petite
            # isovaleur (ou progression) de la composante.
            roots = [min(component, key=lambda n: (key_of(n), n))]
        for root in roots:
            visit(root)
        # filet de securite : un noeud non atteint depuis les racines
        # detectees est quand meme ajoute, en fin de composante, par
        # isovaleur croissante - ne perd jamais un Layer silencieusement.
        for node in sorted(component, key=lambda n: (raw_layers[n][0], n)):
            if node not in visited:
                visit(node)

    return order


def _spatial_dfs_order(
    raw_layers: list[tuple[float, list[np.ndarray], list[int], tuple]],
    region_ranks: dict[int, int],
    gap_threshold: float,
    progress: Optional[np.ndarray] = None,
) -> list[int]:
    """
    Ordre final des indices de `raw_layers` : rang topologique de region
    croissant (contrainte physique dure, JAMAIS traversee - partition
    stricte des noeuds par rang avant tout parcours, pas une ponderation),
    puis continuite spatiale DFS a l'interieur de chaque rang
    (`_dfs_order_within_rank`). Retourne une permutation complete des
    indices de `raw_layers`.
    """
    graph = _build_continuity_graph(raw_layers, gap_threshold, progress)

    by_rank: dict[int, set[int]] = {}
    for idx, (_isovalue, _contours, region_ids, _frame) in enumerate(raw_layers):
        rank = max((region_ranks.get(r, 0) for r in region_ids), default=0)
        by_rank.setdefault(rank, set()).add(idx)

    order: list[int] = []
    for rank in sorted(by_rank):
        order.extend(
            _dfs_order_within_rank(graph, by_rank[rank], raw_layers, progress)
        )
    return order
