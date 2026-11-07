"""Suite découpée de `test_step5_layers.py` (R-030 point 1, limite de 800 lignes par fichier) : tests déplacés tels quels."""

from __future__ import annotations

import numpy as np
import trimesh
from core.mesh_types import RegionDependencyGraph
from element_ardor.step4_field_solver import FieldSolverParams, run as field_run
from element_ardor.step5_layers import LayerParams, _march_tetrahedra_batch, _march_tetrahedron, _spatial_dfs_order, run as layers_run
from test_step5_layers import (  # noqa: F401
    _pipeline_data_from_trimesh_default,
    hollow_box_data,
    portal_data,
    t_shape_single_region_data,
    tripod_data,
)


def test_march_tetrahedra_batch_matches_scalar_loop():
    """
    Compare, triangle par triangle et EN ORDRE (pas de tri qui masquerait
    un bug de winding), la sortie de _march_tetrahedra_batch a celle de
    _march_tetrahedron appelee ligne par ligne - sur un melange de
    tetraedres a valeurs continues (cas generique) et a valeurs discretes
    (force des egalites/limites n_above, y compris n_above in (0,4) et
    champ constant degenere hi==lo) pour couvrir tous les chemins.
    """
    rng = np.random.default_rng(2024)

    for trial in range(30):
        m = int(rng.integers(50, 400))
        tet_points = rng.uniform(-5.0, 5.0, size=(m, 4, 3))

        # ~30% valeurs discretes (force des egalites/cas limites), ~5%
        # constantes par ligne (degenere hi==lo), le reste continu.
        tet_values = rng.uniform(0.0, 4.0, size=(m, 4))
        discrete_mask = rng.random(m) < 0.30
        tet_values[discrete_mask] = rng.choice(
            [0.0, 1.0, 2.0, 3.0, 1.0], size=(discrete_mask.sum(), 4)
        )
        constant_mask = rng.random(m) < 0.05
        tet_values[constant_mask] = rng.uniform(0.0, 4.0, size=(constant_mask.sum(), 1))

        isovalue = float(rng.uniform(0.0, 4.0))

        batch_vertices, batch_local_idx = _march_tetrahedra_batch(tet_points, tet_values, isovalue)
        batch_triangles = batch_vertices.reshape(-1, 3, 3)

        for row in range(m):
            expected = _march_tetrahedron(tet_points[row], tet_values[row], isovalue)
            actual = batch_triangles[batch_local_idx == row]

            assert len(actual) == len(expected), (
                f"trial={trial} row={row}: {len(actual)} triangles obtenus, "
                f"{len(expected)} attendus (isovalue={isovalue}, "
                f"tet_values={tet_values[row]})"
            )
            for tri_actual, tri_expected in zip(actual, expected):
                assert np.allclose(tri_actual, tri_expected, atol=1e-9), (
                    f"trial={trial} row={row}: triangle different (ordre/winding "
                    f"inclus) - obtenu {tri_actual}, attendu {tri_expected}"
                )


def test_march_tetrahedra_batch_all_degenerate_returns_empty():
    """Toutes les lignes non-croisantes (n_above in (0,4)) ou a champ
    constant : aucun triangle, aucune exception - le repli 'vertices_parts
    vide' doit etre exerce."""
    tet_points = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ]
    )
    tet_values = np.array(
        [
            [10.0, 10.0, 10.0, 10.0],  # n_above == 4 (tout au-dessus)
            [-10.0, -10.0, -10.0, -10.0],  # n_above == 0 (tout en-dessous)
            [2.0, 2.0, 2.0, 2.0],  # champ constant, degenere (hi==lo)
        ]
    )

    vertices, local_idx = _march_tetrahedra_batch(tet_points, tet_values, isovalue=2.0)

    assert len(vertices) == 0
    assert len(local_idx) == 0


def test_march_tetrahedra_batch_empty_input_returns_empty():
    vertices, local_idx = _march_tetrahedra_batch(
        np.zeros((0, 4, 3)), np.zeros((0, 4)), isovalue=1.0
    )
    assert len(vertices) == 0
    assert len(local_idx) == 0


def test_split_disjoint_clusters_default_is_disabled(hollow_box_data):
    """
    Non-regression stricte : sans split_disjoint_clusters (defaut None), le
    comportement est identique a aujourd'hui - la couche mi-hauteur de la
    boite creuse produit toujours un seul Layer avec 2 contour_3d.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = hollow_box_data
    params = LayerParams(layer_height_base=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    mid_layers = [
        layer for layer in layers
        if abs(np.mean([c[:, 2].mean() for c in layer.contour_3d])) < 2.0
    ]
    multi_contour_layers = [layer for layer in mid_layers if len(layer.contour_3d) >= 2]
    assert multi_contour_layers, (
        "Comportement historique attendu : la couche mi-hauteur garde ses "
        "2 contours dans un seul Layer quand le clustering est desactive"
    )


def test_hollow_box_nested_contours_not_split_even_when_enabled(hollow_box_data):
    """
    Cas le plus critique : meme avec split_disjoint_clusters ACTIVE, deux
    boucles IMBRIQUEES (paroi exterieure + cavite interne) ne doivent
    JAMAIS etre scindees en 2 Layer distincts - ce sont une seule section
    transversale avec un trou, pas deux clusters separes (CLAUDE.md).
    """
    mesh_data, tet_mesh, overhang_map, regions, field = hollow_box_data
    params = LayerParams(layer_height_base=1.0, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    mid_layers = [
        layer for layer in layers
        if abs(np.mean([c[:, 2].mean() for c in layer.contour_3d])) < 2.0
    ]
    multi_contour_layers = [layer for layer in mid_layers if len(layer.contour_3d) >= 2]
    assert multi_contour_layers, (
        "Le clustering actif ne doit PAS scinder les boucles imbriquees de "
        "la boite creuse - elles doivent rester dans un seul Layer"
    )


def test_tripod_split_disjoint_clusters_produces_distinct_local_frames(tripod_data):
    """
    Cas cible du fix : sur le modele tripode (jonction reelle a 4 regions),
    avec split_disjoint_clusters active, au moins une isovaleur partagee
    par plusieurs Layer distincts doit exister, chacun avec un local_frame
    dont la normale differe sensiblement des autres - preuve directe que le
    plan moyen incoherent (cause du faisceau convergent observe
    visuellement, CLAUDE.md) est resolu par construction.
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph,
        FieldSolverParams(),
    )
    params = LayerParams(layer_height_base=0.2, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, region_graph, params)

    by_isovalue: dict[float, list] = {}
    for layer in layers:
        key = round(layer.isovalue, 6)
        by_isovalue.setdefault(key, []).append(layer)

    split_groups = [group for group in by_isovalue.values() if len(group) >= 2]
    assert split_groups, (
        "Au moins une isovaleur devrait avoir ete scindee en plusieurs Layer "
        "sur le modele tripode avec split_disjoint_clusters actif"
    )

    found_distinct_normals = False
    for group in split_groups:
        normals = [np.cross(l.local_frame[1], l.local_frame[2]) for l in group]
        for i in range(len(normals)):
            for j in range(i + 1, len(normals)):
                if abs(np.dot(normals[i], normals[j])) < 0.9:
                    found_distinct_normals = True
    assert found_distinct_normals, (
        "Au moins deux Layer issus d'une meme isovaleur scindee devraient "
        "avoir des plans locaux nettement distincts (produit scalaire des "
        "normales < 0.9), preuve que chaque cluster a son propre plan "
        "coherent plutot qu'un plan moyen incoherent"
    )


def test_tripod_split_disjoint_clusters_region_ids_are_distinct_subsets(tripod_data):
    """
    Les region_ids des Layer resultants d'une isovaleur scindee doivent
    etre des sous-ensembles distincts (pas tous identiques) - preuve que le
    recalcul par cluster (approximation par proximite spatiale, CLAUDE.md)
    associe bien chaque cluster a SA region, pas a l'union de toutes.
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph,
        FieldSolverParams(),
    )
    params = LayerParams(layer_height_base=0.2, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, region_graph, params)

    by_isovalue: dict[float, list] = {}
    for layer in layers:
        key = round(layer.isovalue, 6)
        by_isovalue.setdefault(key, []).append(layer)

    split_groups = [group for group in by_isovalue.values() if len(group) >= 2]
    assert split_groups, "Au moins une isovaleur scindee est necessaire pour ce test"

    found_distinct_region_sets = False
    for group in split_groups:
        region_sets = [frozenset(l.region_ids) for l in group]
        if len(set(region_sets)) > 1:
            found_distinct_region_sets = True
    assert found_distinct_region_sets, (
        "Au moins un groupe de Layer scindes devrait avoir des region_ids "
        "distincts entre clusters (pas tous identiques)"
    )


def test_close_disjoint_contours_not_split_within_tolerance():
    """
    Deux boites reliees par un pont fin, avec un ecart controle sous le
    seuil de tolerance (split_disjoint_clusters * layer_height_base) : les
    2 boites doivent rester dans un SEUL Layer a l'endroit du gap - preuve
    que le mecanisme n'introduit pas de discontinuite artificielle avant
    une fusion geometrique reelle a une isovaleur voisine.
    """
    box_a = trimesh.creation.box(extents=[4, 4, 4])
    box_a.apply_translation([-2.5, 0, 10])
    box_b = trimesh.creation.box(extents=[4, 4, 4])
    box_b.apply_translation([2.5, 0, 10])
    bridge = trimesh.creation.box(extents=[5, 1, 1])
    bridge.apply_translation([0, 0, 10])

    mesh = box_a.union(bridge, engine="manifold").union(box_b, engine="manifold")
    mesh_data, tet_mesh, overhang_map, regions, field = _pipeline_data_from_trimesh_default(
        mesh, resolution=1.0, platform_z=8.0
    )

    layer_height_base = 0.5
    # tolerance large (10x layer_height_base) : l'objectif de ce test est de
    # prouver l'absence de scission artificielle quand le gap reel est petit
    # devant le seuil, pas de calibrer un seuil de production - la
    # calibration fine se fait sur le modele reel (CLAUDE.md).
    params = LayerParams(layer_height_base=layer_height_base, split_disjoint_clusters=10.0)
    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    # les couches proches de z=10 (hauteur des 2 boites + pont) doivent
    # rester des Layer uniques (pas de scission artificielle malgre la
    # geometrie "en haltere" - le pont, bien que fin, garde les boucles
    # suffisamment proches en distance de bord-a-bord).
    near_bridge = [layer for layer in layers if any(abs(c[:, 2].mean() - 10.0) < 1.0 for c in layer.contour_3d)]
    assert near_bridge, "Au moins une couche proche du pont est necessaire pour ce test"

    by_isovalue: dict[float, list] = {}
    for layer in near_bridge:
        key = round(layer.isovalue, 6)
        by_isovalue.setdefault(key, []).append(layer)

    assert all(len(group) == 1 for group in by_isovalue.values()), (
        "Aucune isovaleur proche du pont ne devrait etre scindee en "
        "plusieurs Layer : l'ecart entre les 2 boites doit rester sous la "
        "tolerance de distance"
    )


def _rank_sequence(layers, region_ranks_by_id: dict[int, int]) -> list[int]:
    ranks = []
    for layer in layers:
        rank = max((region_ranks_by_id.get(r, 0) for r in layer.region_ids), default=0)
        ranks.append(rank)
    return ranks


def _branch_switch_count(layers, jump_threshold: float) -> int:
    """
    Nombre de fois ou l'origine du local_frame "saute" de plus de
    `jump_threshold` entre deux couches consecutives - proxy direct du
    nombre de trajets a vide inseres par step7/step8 entre deux clusters
    spatialement distincts (CLAUDE.md, faisceau convergent mesure sur le
    modele en T reel).

    Se base sur la POSITION du cluster (local_frame.origin), pas sur
    region_ids : deux clusters spatialement disjoints peuvent legitimement
    partager les memes region_ids (une seule region step3 peut couvrir
    plusieurs composantes spatiales du solide, cf. CLAUDE.md - c'est
    exactement le cas cible de ce fix), donc region_ids seul ne peut pas
    detecter la bascule.
    """
    switches = 0
    for i in range(1, len(layers)):
        prev_origin = layers[i - 1].local_frame[0]
        cur_origin = layers[i].local_frame[0]
        if np.linalg.norm(cur_origin - prev_origin) > jump_threshold:
            switches += 1
    return switches


def test_t_shape_spatial_order_groups_branch_before_switching(t_shape_single_region_data):
    """
    Cas cible exact du fix (reproduit le bug mesure sur
    res/test_models/l-triple-axis-overhanged-shape.stl) : pied et bras
    partagent la meme region step3 (meme nappe de surface connexe), donc le
    tri (rank, isovalue) seul les fait alterner a quasi chaque isovaleur.
    Le tri spatial DFS doit produire strictement moins de bascules de
    branche que ce tri naif reconstruit independamment sur les memes
    couches.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = t_shape_single_region_data
    params = LayerParams(layer_height_base=0.5, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    split_isovalues = {
        round(layer.isovalue, 6)
        for layer in layers
        if sum(1 for other in layers if round(other.isovalue, 6) == round(layer.isovalue, 6)) >= 2
    }
    assert split_isovalues, (
        "Au moins une isovaleur scindee est necessaire pour que ce test "
        "prouve quoi que ce soit sur l'ordre - verifier la geometrie de la "
        "fixture si ceci echoue"
    )

    # seuil de saut : la moitie de la distance separant les 2 branches (16
    # unites entre pied et bras dans cette fixture) - largement au-dessus
    # du deplacement normal d'une couche a l'autre au sein d'une meme
    # branche (~layer_height_base), largement en-dessous du saut reel
    # pied<->bras.
    jump_threshold = 5.0
    naive_order = sorted(layers, key=lambda l: (l.isovalue,))
    naive_switches = _branch_switch_count(naive_order, jump_threshold)
    spatial_switches = _branch_switch_count(layers, jump_threshold)

    assert spatial_switches < naive_switches, (
        f"Le tri spatial DFS ({spatial_switches} bascules) devrait produire "
        f"nettement moins de bascules de branche que l'ordre par isovaleur "
        f"seule ({naive_switches} bascules)"
    )


def test_tripod_spatial_order_never_violates_region_rank(tripod_data):
    """
    Non-regression stricte de la contrainte dure : le tri spatial DFS ne
    doit JAMAIS reordonner un cluster de rang superieur avant un cluster de
    rang inferieur, meme actif. Rangs non-decroissants dans la liste
    retournee.
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph,
        FieldSolverParams(),
    )
    params = LayerParams(layer_height_base=0.2, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, region_graph, params)

    region_ranks_by_id = {0: 0, 1: 1, 2: 1, 3: 1}
    ranks = _rank_sequence(layers, region_ranks_by_id)

    assert ranks == sorted(ranks), (
        "Le tri spatial DFS ne doit jamais violer le rang topologique de "
        "region - sequence obtenue non triee"
    )


def test_tripod_spatial_order_deterministic(tripod_data):
    """Deux appels identiques doivent produire exactement le meme ordre - garantit
    l'absence de dependance a un ordre d'iteration non deterministe de networkx."""
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph,
        FieldSolverParams(),
    )
    params = LayerParams(layer_height_base=0.2, split_disjoint_clusters=3.0)

    layers_a = layers_run(tet_mesh, mesh_data, field, region_graph, params)
    layers_b = layers_run(tet_mesh, mesh_data, field, region_graph, params)

    assert [round(l.isovalue, 9) for l in layers_a] == [round(l.isovalue, 9) for l in layers_b]
    assert [frozenset(l.region_ids) for l in layers_a] == [frozenset(l.region_ids) for l in layers_b]


def test_hollow_box_spatial_order_keeps_nested_contours_single_layer(hollow_box_data):
    """
    Non-regression : le tri spatial DFS (actif automatiquement des que
    split_disjoint_clusters l'est) ne doit introduire aucune scission
    indirecte des boucles imbriquees de la boite creuse - meme nombre de
    Layer qu'avec le seul clustering (sans consequence de l'ordre).
    """
    mesh_data, tet_mesh, overhang_map, regions, field = hollow_box_data
    params = LayerParams(layer_height_base=1.0, split_disjoint_clusters=3.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    mid_layers = [
        layer for layer in layers
        if abs(np.mean([c[:, 2].mean() for c in layer.contour_3d])) < 2.0
    ]
    assert all(len(layer.contour_3d) == 2 for layer in mid_layers), (
        "Les boucles imbriquees de la boite creuse doivent rester dans un "
        "seul Layer (2 contour_3d), le graphe de continuite ne doit pas "
        "introduire de scission indirecte"
    )


def test_portal_fronts_stay_split_until_they_meet(portal_data):
    """
    Le champ est une distance au plateau a travers la matiere : les deux montants
    sont deux fronts qui montent en parallele, traversent la traverse chacun de
    leur cote et ne se rejoignent qu'a la fin. Le clustering doit donc scinder
    chaque isovaleur en 2 Layer (un par front) presque partout, sans erreur ni
    duplication.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = portal_data
    params = LayerParams(layer_height_base=1.0, split_disjoint_clusters=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    by_isovalue: dict[float, list] = {}
    for layer in layers:
        by_isovalue.setdefault(round(layer.isovalue, 6), []).append(layer)

    split = [g for g in by_isovalue.values() if len(g) >= 2]
    assert len(split) >= 0.9 * len(by_isovalue), (
        "Les deux fronts devraient rester scindes en 2 Layer a presque toutes les isovaleurs"
    )


def test_portal_spatial_order_visits_one_leg_fully_before_switching(portal_data):
    """
    L'ordre spatial DFS doit regrouper TOUTES les couches d'un meme montant
    avant de passer a l'autre - au plus une poignee de bascules (une par
    racine independante visitee - les 2 montants sont 2 composantes non
    connectees sous la traverse, donc 2 racines - pas une bascule a CHAQUE
    isovaleur comme le ferait un tri (rank, isovalue) naif). Compare
    directement au tri naif pour le prouver, plutot que de supposer un
    chiffre absolu (le nombre exact de racines/artefacts de maillage n'est
    pas le point teste ici).
    """
    mesh_data, tet_mesh, overhang_map, regions, field = portal_data
    params = LayerParams(layer_height_base=1.0, split_disjoint_clusters=1.0)

    layers = layers_run(tet_mesh, mesh_data, field, regions, params)

    by_isovalue: dict[float, list] = {}
    for layer in layers:
        key = round(layer.isovalue, 6)
        by_isovalue.setdefault(key, []).append(layer)
    split_isovalues = {iv for iv, g in by_isovalue.items() if len(g) >= 2}
    assert split_isovalues, "Au moins une isovaleur scindee est necessaire pour ce test"

    # seuil de saut : la moitie de la distance separant les 2 montants (16
    # unites, cf. _portal_mesh) - au-dessus du deplacement normal d'une
    # couche a l'autre au sein d'un meme montant.
    jump_threshold = 5.0
    naive_order = sorted(layers, key=lambda l: (l.isovalue,))
    naive_switches = _branch_switch_count(naive_order, jump_threshold)
    spatial_switches = _branch_switch_count(layers, jump_threshold)

    assert spatial_switches < naive_switches, (
        f"L'ordre spatial DFS ({spatial_switches} bascules) devrait produire "
        f"nettement moins de bascules que l'ordre par isovaleur seule "
        f"({naive_switches} bascules) sur le portique"
    )


def test_portal_merge_node_not_duplicated(portal_data):
    """
    Garde-fou direct contre une regression future de _dfs_order_within_rank :
    aucune duplication du noeud de recombinaison malgre ses 2 aretes
    entrantes. Se base sur (isovalue, origine du local_frame arrondie) -
    PAS (isovalue, region_ids) seul, qui coincide legitimement pour les 2
    montants du portique (une seule region pour tout le modele, cf.
    portal_data) sans que ce soit une duplication.
    """
    mesh_data, tet_mesh, overhang_map, regions, field = portal_data

    params_with_sort = LayerParams(layer_height_base=1.0, split_disjoint_clusters=1.0)
    layers_with_sort = layers_run(tet_mesh, mesh_data, field, regions, params_with_sort)

    keys_seen = [
        (round(l.isovalue, 6), tuple(np.round(l.local_frame[0], 3)))
        for l in layers_with_sort
    ]
    assert len(keys_seen) == len(set(keys_seen)), (
        "Un meme (isovalue, origine de cluster) apparait plusieurs fois - "
        "signe d'une duplication du noeud de recombinaison par le DFS"
    )


def test_spatial_order_is_a_permutation():
    """
    Test unitaire direct sur _spatial_dfs_order, isole du pipeline complet -
    construit des raw_layers/region_ranks synthetiques (style tripode ET
    style portique) et verifie que l'ordre retourne est une vraie
    permutation (aucun index perdu, aucun duplique).
    """
    identity_frame = (np.zeros(3), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))

    def _square_loop(cx: float, cy: float, cz: float, half: float = 1.0) -> np.ndarray:
        pts = np.array([
            [cx - half, cy - half, cz],
            [cx + half, cy - half, cz],
            [cx + half, cy + half, cz],
            [cx - half, cy + half, cz],
            [cx - half, cy - half, cz],
        ])
        return pts

    # Style tripode : 1 pied (rank 0) qui se scinde en 2 branches (rank 1),
    # chacune une chaine simple sur 2 isovaleurs supplementaires.
    raw_layers_tripod = [
        (0.0, [_square_loop(0, 0, 0)], [0], identity_frame),
        (1.0, [_square_loop(-3, 0, 1)], [1], identity_frame),
        (1.0, [_square_loop(3, 0, 1)], [1], identity_frame),
        (2.0, [_square_loop(-3, 0, 2)], [1], identity_frame),
        (2.0, [_square_loop(3, 0, 2)], [1], identity_frame),
    ]
    region_ranks_tripod = {0: 0, 1: 1}
    order_tripod = _spatial_dfs_order(raw_layers_tripod, region_ranks_tripod, gap_threshold=1.5)
    assert sorted(order_tripod) == list(range(len(raw_layers_tripod)))

    # Style portique : 2 montants separes (rank 0) qui fusionnent en 1 seul
    # cluster (rank 0 aussi - meme region pour tout le portique).
    raw_layers_portal = [
        (0.0, [_square_loop(-5, 0, 0)], [0], identity_frame),
        (0.0, [_square_loop(5, 0, 0)], [0], identity_frame),
        (1.0, [_square_loop(-5, 0, 1)], [0], identity_frame),
        (1.0, [_square_loop(5, 0, 1)], [0], identity_frame),
        (2.0, [_square_loop(0, 0, 2, half=6.0)], [0], identity_frame),
    ]
    region_ranks_portal = {0: 0}
    order_portal = _spatial_dfs_order(raw_layers_portal, region_ranks_portal, gap_threshold=1.5)
    assert sorted(order_portal) == list(range(len(raw_layers_portal)))


def test_isovalue_order_ignores_region_rank(tripod_data):
    """
    Champ qui EST la progression (distance au plateau) : le rang de region ne
    doit plus primer. Dependances inversees (bras -> pied) pour que le tri
    historique (rang, isovaleur) imprime les bras, d'isovaleur haute, avant le
    pied : sans ce desaccord le test ne prouverait rien.
    """
    mesh_data, tet_mesh, overhang_map, region_graph = tripod_data
    field = field_run(
        tet_mesh, mesh_data, overhang_map, region_graph,
        FieldSolverParams(),
    )
    inverted = RegionDependencyGraph(
        regions=region_graph.regions,
        dependencies=[(1, 0), (2, 0), (3, 0)],
        platform_contact_faces=region_graph.platform_contact_faces,
    )

    historic = layers_run(tet_mesh, mesh_data, field, inverted,
                          LayerParams(layer_height_base=0.5, isovalue_order=False))
    iso_historic = [layer.isovalue for layer in historic]
    assert iso_historic != sorted(iso_historic), "le tri par rang doit ici contredire l'isovaleur"

    ordered = layers_run(
        tet_mesh, mesh_data, field, inverted, LayerParams(layer_height_base=0.5, isovalue_order=True)
    )
    iso_ordered = [layer.isovalue for layer in ordered]
    assert iso_ordered == sorted(iso_ordered)
    assert sorted(iso_ordered) == sorted(iso_historic), "meme ensemble de couches, seul l'ordre change"


def test_isovalue_order_is_the_default():
    """Le champ est une distance au plateau : le tri par isovaleur est le defaut."""
    assert LayerParams().isovalue_order is True
