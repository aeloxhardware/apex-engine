"""
pipeline/step3_segmentation.py

Etape 3 du pipeline : segmentation du mesh en sous-volumes (Region) partageant
une direction d'impression viable, et construction du graphe de dependances
physiques entre ces regions.

Depend de MeshData (step1) et OverhangMap (step2).

Cette etape resout aussi une limite laissee ouverte a l'etape 2 : la
distinction entre "face reposant sur le plateau d'impression" et "face
reellement en surplomb dans le vide". step2_overhang.py n'a aucune
connaissance du plateau ; c'est ici, ou l'on raisonne pour la premiere fois
en termes de region/impression reelle, que cette information est injectee.

Sortie : RegionDependencyGraph
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import networkx as nx
import numpy as np
import scipy.sparse as sp

import scipy.sparse.csgraph as spgraph
import scipy.spatial
import scipy.spatial.distance

from core.geometry_utils import (
    build_face_adjacency,
    composed_rotation_directions,
    fibonacci_sphere_directions,
)
from core.mesh_types import (
    MeshData,
    OverhangMap,
    Region,
    RegionDependencyGraph,
    RotationAxisSpec,
)
from core.spatial import KDTree


class SegmentationError(Exception):
    """Levee si la segmentation produit un resultat physiquement invalide (cycle de dependances)."""


@dataclass
class SegmentationParams:
    direction_candidate_count: int = 26
    max_angle_from_direction_deg: float = 45.0
    # Seuil angulaire (degres) au-dela duquel une fusion de deux regions
    # adjacentes est refusee, si la direction optimale de la region
    # combinee s'ecarte trop de la direction optimale d'au moins une des
    # deux regions prises separement (voir _region_growing). None (defaut)
    # = DESACTIVE, comportement historique inchange : une fusion est
    # acceptee des qu'au moins une direction candidate reste compatible
    # avec les deux regions, peu importe l'angle. Opt-in explicite, meme
    # discipline que direction_blend_radius=0.0 (step4) et
    # near_platform_resolution=None (step1).
    #
    # PRECEDEMMENT tente : un critere de tolerance sur le RISQUE RESIDUEL
    # scalaire (region_count_penalty, ratio de degradation) plutot qu'un
    # angle - ECARTE apres mesure directe sur un modele en T reel (2 bras
    # perpendiculaires) : avec direction_candidate_count~26, il existe
    # presque toujours une direction candidate de residu quasi nul meme
    # pour une region combinant deux bras incompatibles (la sphere de
    # candidats est assez dense pour qu'un compromis numerique existe sans
    # etre geometriquement coherent) - mesure : own_best=shared_best=0.0
    # pour TOUTES les fusions testees, le critere ne se declenchait donc
    # JAMAIS meme a tolerance 100%. Ce defaut S'AGGRAVE avec la densite de
    # candidats, donc particulierement dangereux pour une cinematique
    # future dense/continue (voir CLAUDE.md limite #18). Le critere
    # angulaire ci-dessus ne depend PAS du nombre de candidats (un angle
    # entre deux directions reste un angle), donc robuste par construction
    # a cette classe de probleme - voir CLAUDE.md pour le detail complet.
    max_region_merge_angle_deg: Optional[float] = None
    # None (defaut) = DETECTION AUTOMATIQUE depuis le mesh (min(z) des
    # sommets) - meme discipline opt-in que max_region_merge_angle_deg
    # ci-dessus. Une valeur explicite reste un override total. Corrige un
    # bug reel : un ancien defaut fige a 0.0 laissait
    # `_compute_platform_contact` ne detecter AUCUNE face de contact sur
    # tout modele dont le point le plus bas n'est pas exactement a z=0 (ex:
    # res/test_models/death.stl, min(z)=1.245) - cassant en cascade
    # l'ancrage plateau du solveur de champ (step4). Voir CLAUDE.md.
    platform_z: Optional[float] = None
    platform_contact_tolerance: float = 0.5
    # Si None ou liste vide : utilise fibonacci_sphere_directions (orientation
    # libre, 5 axes). Si fourni : decrit une cinematique de rotation
    # arbitraire (N axes composes, chacun avec sa propre plage angulaire) via
    # composed_rotation_directions - PAS figee sur une seule machine. Un seul
    # axe a plage complete (0,360) degenere exactement vers l'ancien
    # comportement a 1 DOF (cercle) ; plusieurs axes decrivent une cinematique
    # composee (ex: plateau + tete d'une machine a 2 DOF).
    rotation_axes: Optional[list[RotationAxisSpec]] = None
    up_vector: np.ndarray = None

    def __post_init__(self):
        if self.up_vector is None:
            self.up_vector = np.array([0.0, 0.0, 1.0])


# ---------------------------------------------------------------------------
# Etape 1 - Filtrage des faces au contact du plateau
# ---------------------------------------------------------------------------

def _compute_platform_contact(
    mesh_data: MeshData, platform_z: float, tolerance: float
) -> np.ndarray:
    """
    Marque les faces dont tous les sommets sont a moins de `tolerance` du
    plan du plateau (z = platform_z). Ces faces reposent physiquement sur
    le plateau et ne doivent pas etre traitees comme des surplombs par la
    suite de la segmentation, meme si leur normale pointe vers le bas.
    """
    vertex_z = mesh_data.vertices[:, 2]
    face_vertex_z = vertex_z[mesh_data.faces]  # shape (M, 3)
    max_z_per_face = face_vertex_z.max(axis=1)
    return max_z_per_face < (platform_z + tolerance)


def _effective_risk_score(
    overhang_map: OverhangMap, is_platform_contact: np.ndarray
) -> np.ndarray:
    """
    Score de risque a utiliser pour la segmentation : identique a
    OverhangMap.risk_score, sauf pour les faces au contact du plateau, dont
    le risque est force a 0 (elles n'ont pas besoin d'une direction
    d'impression alternative, elles reposent deja sur une surface stable).
    """
    effective = overhang_map.risk_score.copy()
    effective[is_platform_contact] = 0.0
    return effective


# ---------------------------------------------------------------------------
# Etape 3 - Directions candidates et compatibilite par face
# ---------------------------------------------------------------------------

def _generate_candidate_directions(params: SegmentationParams) -> np.ndarray:
    """
    Genere l'ensemble des directions d'impression candidates a tester,
    toujours avec la verticale (up_vector) en premiere position (direction
    par defaut, prioritaire en cas d'egalite de cout).

    Dispatch generique sur la cinematique decrite par `rotation_axes` : ni
    None/vide (orientation libre) ni une liste d'axes quelconque n'est un
    cas special code en dur ici - composed_rotation_directions absorbe
    uniformement 0, 1 ou N axes (voir CLAUDE.md limite #18 : ne pas figer
    step3 sur une seule cinematique de machine).
    """
    up = params.up_vector / np.linalg.norm(params.up_vector)

    if not params.rotation_axes:
        sphere = fibonacci_sphere_directions(params.direction_candidate_count)
        directions = np.vstack([up.reshape(1, 3), sphere])
    else:
        composed = composed_rotation_directions(
            params.rotation_axes,
            reference_direction=up,
            samples_per_axis=_samples_per_axis(
                params.direction_candidate_count, len(params.rotation_axes)
            ),
        )
        directions = np.vstack([up.reshape(1, 3), composed])

    return directions


def _samples_per_axis(direction_candidate_count: int, axis_count: int) -> int:
    """
    Convertit un nombre de directions cible en nombre d'echantillons PAR
    AXE, puisque composed_rotation_directions echantillonne une grille
    produit (total = samples_per_axis ** axis_count, modulo les DOF figes).
    Racine axis_count-ieme, borne a 2 au minimum pour qu'un axe a plage non
    degeneree produise toujours au moins 2 points distincts.
    """
    if axis_count <= 0:
        return direction_candidate_count
    return max(int(round(direction_candidate_count ** (1.0 / axis_count))), 2)


def _face_direction_compatibility(
    face_normals: np.ndarray, directions: np.ndarray, max_angle_deg: float
) -> np.ndarray:
    """
    Pour chaque face et chaque direction candidate, determine si imprimer
    dans cette direction rendrait la face imprimable (i.e. elle n'est PAS
    en surplomb selon cette direction).

    Convention : "imprimer dans la direction d" signifie que d joue le role
    de la verticale locale (la buse depose selon -d, donc "vers le bas
    local" = -d). Une face est en surplomb si sa normale pointe a moins de
    max_angle_deg de "vers le bas local" (-d) - exactement le meme critere
    que step2_overhang._compute_is_overhang, generalise a une direction
    quelconque plutot que la seule verticale machine.

    Important : une face dont la normale est perpendiculaire a `d` (angle
    proche de 90 deg, typiquement un mur vertical relativement a la
    verticale machine) n'est en surplomb selon AUCUNE direction candidate -
    elle est donc compatible partout. C'est le comportement voulu : un mur
    ne contraint jamais le choix de direction d'une region, contrairement a
    une face en surplomb franc qui n'est compatible qu'avec un sous-
    ensemble restreint de directions.

    Retourne un masque booleen (M, D) : compatible[f, d] = True si la face f
    n'est PAS en surplomb dans la direction candidate d (donc imprimable).
    """
    norms = np.linalg.norm(face_normals, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    normalized_normals = face_normals / norms

    # angle entre chaque normale et "vers le bas local" (-d) pour chaque
    # direction candidate d
    cos_angles = normalized_normals @ (-directions).T  # shape (M, D)
    cos_angles = np.clip(cos_angles, -1.0, 1.0)
    angle_from_down_deg = np.degrees(np.arccos(cos_angles))

    is_overhang_for_direction = angle_from_down_deg < max_angle_deg
    return ~is_overhang_for_direction


def _angle_between_deg(a: np.ndarray, b: np.ndarray) -> float:
    """Angle en degres entre deux vecteurs unitaires (produit scalaire, clip pour stabilite numerique)."""
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0))))


def _region_growing(
    face_adjacency: sp.csr_matrix,
    compatibility: np.ndarray,
    risk_score: np.ndarray,
    directions: np.ndarray,
    max_region_merge_angle_deg: Optional[float],
) -> np.ndarray:
    """
    Regroupe les faces en regions par croissance iterative : deux faces
    adjacentes fusionnent si elles partagent au moins une direction
    candidate compatible commune. Les faces a risque eleve sont traitees en
    priorite (elles contraignent le plus fortement le choix de direction).

    `max_region_merge_angle_deg` (None = desactive, comportement
    historique) : si fourni, une fusion n'est acceptee que si la direction
    optimale de la region combinee (argmin du risque residuel, restreint
    aux directions deja compatibles avec les deux) reste a moins de
    `max_region_merge_angle_deg` de la direction optimale de CHACUNE des
    deux regions prises separement - empeche deux sous-volumes
    geometriquement incompatibles (ex: deux bras perpendiculaires) de
    fusionner sur la seule base qu'une direction de compromis reste
    techniquement compatible avec les deux, meme si elle est loin de
    l'optimum de chacune (CLAUDE.md, cas mesure sur un modele en T).

    Un critere de TOLERANCE SUR LE RISQUE RESIDUEL (plutot qu'un angle) a
    ete tente en premier et ECARTE apres mesure directe : avec un nombre
    de candidats meme modere (~26), il existe presque toujours une
    direction de residu quasi nul pour n'importe quelle region, y compris
    une fusion geometriquement incoherente - le residu ne capture pas la
    divergence directionnelle. Un angle entre deux directions ne depend
    PAS de la densite de candidats (contrairement a un minimum sur un grand
    ensemble), donc robuste par construction a une cinematique dense ou
    continue (voir CLAUDE.md limite #18).

    Retourne un tableau region_id (M,) - un identifiant de region par face.
    """
    face_count = compatibility.shape[0]
    region_id = np.arange(face_count)  # chaque face demarre comme sa propre region
    # ensemble de directions compatibles restantes par region (intersection
    # cumulative des faces qui la composent)
    region_compat = {i: compatibility[i].copy() for i in range(face_count)}

    quality_gate_enabled = max_region_merge_angle_deg is not None
    region_residual: dict[int, np.ndarray] = {}
    if quality_gate_enabled:
        # risque residuel par direction candidate, cumule par region : une
        # face contribue son risque a chaque direction pour laquelle elle
        # N'EST PAS compatible (meme quantite que _select_region_directions,
        # mais maintenue incrementalement plutot que recalculee a chaque
        # paire testee - cout O(D) par fusion, pas O(taille de region * D)).
        for i in range(face_count):
            region_residual[i] = (~compatibility[i]).astype(np.float64) * risk_score[i]

    adjacency_coo = face_adjacency.tocoo()
    # traite les paires de faces les plus a risque en premier : une fusion
    # impliquant une face a fort risque contraint plus tot le choix de
    # direction, ce qui limite les retours en arriere.
    edge_risk = risk_score[adjacency_coo.row] + risk_score[adjacency_coo.col]
    # Tri STABLE (2026-10-01) : chaque arete figure deux fois avec le meme
    # risque, les ex-aequo sont donc certains, et l'ordre d'un tri NON stable
    # (le `quicksort` de numpy, en realite introsort) depend de son
    # implementation. Mesure : sur turkish, le seul passage au tri stable fait
    # passer la segmentation de 11 a 24 regions. Ex-aequo departages dans
    # l'ordre (ligne, colonne) de la matrice d'adjacence, reproductible (et
    # transcrit tel quel par le natif). step9b sur les 4 modeles : death et
    # les 2 modeles fins inchanges, turkish 0,198 -> 0,172 % de matiere sans
    # appui, collisions 0,21 -> 0,13 %.
    order = np.argsort(-edge_risk, kind="stable")

    def find(f: int) -> int:
        while region_id[f] != f:
            region_id[f] = region_id[region_id[f]]  # path compression
            f = region_id[f]
        return f

    for idx in order:
        f1, f2 = int(adjacency_coo.row[idx]), int(adjacency_coo.col[idx])
        r1, r2 = find(f1), find(f2)
        if r1 == r2:
            continue

        shared_compat = region_compat[r1] & region_compat[r2]
        if not np.any(shared_compat):
            # aucune direction commune ne rend les deux regions imprimables
            # ensemble : fusion refusee.
            continue

        if quality_gate_enabled:
            combined_residual = region_residual[r1] + region_residual[r2]
            idx_combined = int(np.argmin(np.where(shared_compat, combined_residual, np.inf)))
            dir_combined = directions[idx_combined]
            idx_r1 = int(np.argmin(region_residual[r1]))
            idx_r2 = int(np.argmin(region_residual[r2]))
            angle_r1 = _angle_between_deg(directions[idx_r1], dir_combined)
            angle_r2 = _angle_between_deg(directions[idx_r2], dir_combined)
            if max(angle_r1, angle_r2) > max_region_merge_angle_deg:
                # la direction combinee s'ecarte trop de l'optimum d'au
                # moins une des deux regions : fusion refusee, la scission
                # est preferable.
                continue

        # fusionne r2 dans r1
        region_id[r2] = r1
        region_compat[r1] = shared_compat
        if quality_gate_enabled:
            region_residual[r1] = region_residual[r1] + region_residual[r2]
            del region_residual[r2]
        del region_compat[r2]

    # normalise les ids finaux (path compression complete + relabeling contigu)
    final_ids = np.array([find(f) for f in range(face_count)])
    unique_ids, relabeled = np.unique(final_ids, return_inverse=True)
    return relabeled


# ---------------------------------------------------------------------------
# Etape 5 - Selection de la direction finale par region
# ---------------------------------------------------------------------------

def _platform_root_region_ids(
    region_labels: np.ndarray, is_platform_contact: np.ndarray
) -> set[int]:
    """Regions ayant au moins une face au contact du plateau (racines de croissance)."""
    return {
        int(region)
        for region in np.unique(region_labels)
        if bool(np.any(is_platform_contact[region_labels == region]))
    }


def _region_adjacency_graph(region_labels: np.ndarray, face_adjacency) -> nx.Graph:
    """
    Adjacence NON ORIENTEE entre regions : la connectivite reelle de la
    matiere, independante de toute notion de haut/bas. L'orientation vient
    ensuite, et uniquement, de la profondeur depuis le plateau
    (`_compute_print_depths`).
    """
    adjacency = nx.Graph()
    adjacency.add_nodes_from(int(region) for region in np.unique(region_labels))

    coo = face_adjacency.tocoo()
    labels_row = region_labels[coo.row]
    labels_col = region_labels[coo.col]
    cross = labels_row != labels_col
    for r1, r2 in zip(labels_row[cross], labels_col[cross]):
        adjacency.add_edge(int(r1), int(r2))

    return adjacency


def _select_region_directions(
    region_labels: np.ndarray,
    compatibility: np.ndarray,
    risk_score: np.ndarray,
    directions: np.ndarray,
) -> dict[int, np.ndarray]:
    """
    Pour chaque region, choisit parmi les directions candidates compatibles
    avec l'ensemble de ses faces celle qui minimise le risque residuel
    total. Si aucune direction n'est compatible avec toutes les faces de la
    region (ne devrait pas arriver si _region_growing a fonctionne
    correctement, mais gere en filet de securite), choisit la direction qui
    minimise le risque residuel meme partiellement.
    """
    region_direction: dict[int, np.ndarray] = {}

    for region in np.unique(region_labels):
        face_mask = region_labels == region
        region_compat = compatibility[face_mask]  # (n_faces_in_region, D)
        region_risk = risk_score[face_mask]

        # risque residuel par direction candidate = somme des risques des
        # faces NON compatibles avec cette direction
        incompatible = ~region_compat  # (n_faces_in_region, D)
        residual_risk = (incompatible * region_risk[:, None]).sum(axis=0)  # (D,)

        best_direction_idx = int(np.argmin(residual_risk))

        region_direction[int(region)] = directions[best_direction_idx]

    return region_direction


# ---------------------------------------------------------------------------
# Etape 6 - Graphe de dependances entre regions
# ---------------------------------------------------------------------------

def _lowest_region(
    region_ids, region_labels: np.ndarray, face_centers: np.ndarray, up: np.ndarray
) -> int:
    """
    Region dont le point le plus bas le long de `up` est le plus bas de
    l'ensemble `region_ids` - racine de repli quand aucune region du groupe
    n'est au contact du plateau.

    Utilise le MINIMUM par region (pas la moyenne) : ce qui compte pour
    demarrer une croissance, c'est le point le plus proche du plateau, pas
    le centre de masse - une region large et haute avec une seule pointe
    basse est un meilleur point de depart qu'une petite region entierement
    situee plus haut.
    """
    return min(
        region_ids,
        key=lambda rid: float((face_centers[region_labels == rid] @ up).min()),
    )


def _compute_print_depths(
    adjacency: nx.Graph,
    root_ids: set[int],
    region_labels: np.ndarray,
    face_centers: np.ndarray,
    up: np.ndarray,
) -> dict[int, int]:
    """
    Profondeur de chaque region dans l'arborescence de croissance depuis le
    plateau : 0 pour une racine (region au contact du plateau), n pour une
    region atteinte en n sauts de region a travers la matiere.

    C'est la question physiquement pertinente pour l'ordre d'impression -
    "par ou la matiere rejoint-elle le plateau" - et NON "qui est le plus
    bas". Une piece pendante (excroissance suspendue sous une plateforme)
    est plus basse que ce a quoi elle est accrochee, mais ne peut exister
    qu'apres : elle n'atteint le plateau qu'A TRAVERS son point
    d'accrochage, donc a une profondeur strictement superieure.

    `nx.single_source_shortest_path_length` sur un noeud source VIRTUEL
    relie a toutes les racines, plutot qu'un parcours ecrit a la main :
    donne directement un parcours en largeur multi-sources, et networkx
    est deja la dependance de graphe du projet (detection de cycle ici
    meme, tri topologique en step5). Le noeud virtuel n'existe que le temps
    du calcul, il n'apparait jamais dans le graphe de dependances rendu.

    Deux cas limites, tous deux REELS (pas theoriques), traites
    explicitement :
      - aucune racine plateau : le critere de contact exige que les 3
        sommets d'une face tiennent sous le plateau + tolerance, donc un
        objet en contact TANGENT (sphere, arete, pointe) en produit
        legitimement zero. Repli : la region la plus basse sert de racine.
      - regions inatteignables depuis les racines : une cavite interne
        forme une nappe de surface qui ne partage aucune arete avec la
        coque exterieure (cas nominal, cf. _find_internal_cavities en
        step1) ; un mesh multi-corps aussi. Chaque composante connexe sans
        racine est traitee independamment, avec sa propre racine de repli.
        Les profondeurs ne sont comparables qu'A L'INTERIEUR d'une
        composante - seule l'orientation des aretes en depend.
    """
    depths: dict[int, int] = {}

    for component in nx.connected_components(adjacency):
        component_roots = component & root_ids
        if not component_roots:
            component_roots = {
                _lowest_region(component, region_labels, face_centers, up)
            }

        virtual_source = "__platform__"
        search_graph = adjacency.subgraph(component).copy()
        search_graph.add_node(virtual_source)
        for root_id in component_roots:
            search_graph.add_edge(virtual_source, root_id)

        hops = nx.single_source_shortest_path_length(search_graph, virtual_source)
        for region_id, hop_count in hops.items():
            if region_id == virtual_source:
                continue
            depths[int(region_id)] = int(hop_count) - 1  # racines a 0, pas 1

    return depths


def _build_dependency_graph(
    mesh_data: MeshData,
    region_labels: np.ndarray,
    region_directions: dict[int, np.ndarray],
    is_platform_contact: np.ndarray,
    up_vector: np.ndarray,
) -> RegionDependencyGraph:
    """
    Construit le graphe de dependances entre regions sous forme
    d'ARBORESCENCE DE CROISSANCE DEPUIS LE PLATEAU : une region A doit etre
    imprimee avant une region B si B ne rejoint le plateau qu'a travers A.

    L'ordre suit la facon dont la matiere se depose reellement - de proche
    en proche depuis le plateau - et NON une comparaison de hauteur. Le
    critere de hauteur (historique : centroide global le plus bas imprime
    en premier) est faux pour toute piece PENDANTE : une excroissance
    suspendue sous une plateforme est plus basse que ce a quoi elle est
    accrochee, mais ne peut exister qu'apres elle. Deux regions voisines se
    touchent d'ailleurs par une ARETE, cote a cote et jamais empilees :
    aucune comparaison de hauteur entre elles, globale ou locale, ne peut
    exprimer "qui supporte qui".

    Mecanisme : parcours en largeur multi-sources depuis les regions au
    contact du plateau (`_compute_print_depths`), puis orientation de
    chaque arete d'adjacence de la profondeur FAIBLE vers la profondeur
    FORTE. Deux regions de meme profondeur ne portent aucune arete : elles
    poussent en parallele depuis le plateau, aucune contrainte physique ne
    les ordonne.

    ACYCLICITE GARANTIE PAR CONSTRUCTION : la profondeur croit strictement
    le long de chaque arete conservee. C'est une garantie plus forte que
    celle de la reference globale `up_vector` historique, qui existait
    precisement pour eviter les cycles artificiels sur un anneau de regions
    (sphere segmentee en tranches) - ce cas est desormais couvert sans
    dependre d'une reference de direction. Le controle de cycle ci-dessous
    est conserve (contrat documente de `run()`), mais il est devenu
    structurellement inatteignable.

    Corollaire : une region au contact du plateau est a profondeur 0, le
    minimum, donc aucune arete ne peut pointer vers elle - la suppression
    forcee des aretes entrantes des regions plateau, necessaire avec
    l'ancien critere de hauteur, n'a plus lieu d'etre.
    """
    face_adjacency = build_face_adjacency(mesh_data.faces, mesh_data.face_count)

    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    up = up_vector / np.linalg.norm(up_vector)

    graph = nx.DiGraph()
    regions: dict[int, Region] = {}
    for region_id in np.unique(region_labels):
        face_ids = np.where(region_labels == region_id)[0]
        regions[int(region_id)] = Region(
            id=int(region_id),
            face_ids=face_ids,
            target_direction=region_directions[int(region_id)],
        )
        graph.add_node(int(region_id))

    root_ids = _platform_root_region_ids(region_labels, is_platform_contact)
    adjacency = _region_adjacency_graph(region_labels, face_adjacency)
    depths = _compute_print_depths(adjacency, root_ids, region_labels, face_centers, up)

    for r1, r2 in adjacency.edges():
        depth_1, depth_2 = depths[r1], depths[r2]
        if depth_1 < depth_2:
            graph.add_edge(r1, r2)
        elif depth_2 < depth_1:
            graph.add_edge(r2, r1)
        # profondeurs egales : regions soeurs, aucune contrainte d'ordre

    if not nx.is_directed_acyclic_graph(graph):
        cycle = nx.find_cycle(graph)
        raise SegmentationError(
            f"Le graphe de dependances entre regions contient un cycle "
            f"({cycle}) - segmentation physiquement inimprimable dans un "
            f"ordre valide. Ajustez les parametres de segmentation ou "
            f"revoyez la geometrie source."
        )

    dependencies = list(graph.edges())
    return RegionDependencyGraph(
        regions=regions,
        dependencies=dependencies,
        # expose au reste du pipeline : l'etape 4 en a besoin pour ancrer le
        # champ scalaire au plateau. Le deduire du graphe (regions sans arete
        # entrante) serait faux des que la segmentation ne produit qu'une
        # seule region - cette region est alors racine ET couvre tout le
        # modele, plateau compris.
        platform_contact_faces=is_platform_contact,
        region_print_depth=depths,
    )


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def _native_segmentation(mesh_data, effective_risk, is_platform_contact, directions, params):
    """Pas 21b du portage (cpp/src/core/segmentation.cpp) : chemin par defaut
    en natif, IDENTIQUE au Python (ordre des dependances compris). None si le
    noyau natif est inactif."""
    from core import _cpp
    if not _cpp.enabled():
        return None
    module = _cpp.module()
    if module is None or not hasattr(module, "segment_default"):
        return None
    labels, direction_index, depth, deps = module.segment_default(
        np.ascontiguousarray(mesh_data.vertices, dtype=float),
        np.ascontiguousarray(mesh_data.faces, dtype=np.int64),
        np.ascontiguousarray(mesh_data.face_normals, dtype=float),
        np.ascontiguousarray(effective_risk, dtype=float),
        np.ascontiguousarray(is_platform_contact, dtype=np.uint8),
        np.ascontiguousarray(directions, dtype=float),
        float(params.max_angle_from_direction_deg),
        -1.0 if params.max_region_merge_angle_deg is None else float(params.max_region_merge_angle_deg),
        np.ascontiguousarray(params.up_vector, dtype=float),
    )
    regions = {
        int(r): Region(id=int(r), face_ids=np.where(labels == r)[0],
                       target_direction=directions[int(direction_index[r])])
        for r in range(len(direction_index))
    }
    pairs = np.asarray(deps, dtype=np.int64).reshape(-1, 2)
    return RegionDependencyGraph(
        regions=regions,
        dependencies=[(int(a), int(b)) for a, b in pairs],
        platform_contact_faces=is_platform_contact,
        region_print_depth={int(r): int(d) for r, d in enumerate(depth)},
    )


def run(
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    params: SegmentationParams,
) -> RegionDependencyGraph:
    """
    Segmente le mesh en sous-volumes partageant une direction d'impression
    viable, et construit le graphe de dependances physiques entre eux.

    Args:
        mesh_data: sortie de pipeline/step1_ingest.py.
        overhang_map: sortie de pipeline/step2_overhang.py.
        params: seuils et parametres de segmentation.

    Returns:
        RegionDependencyGraph avec les regions, leurs dependances (garanties
        acycliques) et `region_print_depth`, la profondeur de chaque region
        dans l'arborescence de croissance depuis le plateau qui a servi a
        orienter ces dependances.

    Raises:
        SegmentationError: si le graphe de dependances resultant contient
            un cycle (segmentation physiquement inimprimable). Devenu
            structurellement inatteignable depuis que l'ordre derive d'une
            profondeur d'arborescence strictement croissante le long de
            chaque arete - conserve comme garde-fou de contrat.
    """
    resolved_platform_z = (
        float(mesh_data.vertices[:, 2].min()) if params.platform_z is None else params.platform_z
    )
    is_platform_contact = _compute_platform_contact(
        mesh_data, resolved_platform_z, params.platform_contact_tolerance
    )
    effective_risk = _effective_risk_score(overhang_map, is_platform_contact)

    directions = _generate_candidate_directions(params)

    native = _native_segmentation(
        mesh_data, effective_risk, is_platform_contact, directions, params
    )
    if native is not None:
        return native
    compatibility = _face_direction_compatibility(
        mesh_data.face_normals, directions, params.max_angle_from_direction_deg
    )

    face_adjacency = build_face_adjacency(mesh_data.faces, mesh_data.face_count)
    region_labels = _region_growing(
        face_adjacency, compatibility, effective_risk, directions,
        params.max_region_merge_angle_deg,
    )
    region_directions = _select_region_directions(
        region_labels, compatibility, effective_risk, directions
    )

    return _build_dependency_graph(
        mesh_data, region_labels, region_directions, is_platform_contact, params.up_vector
    )
