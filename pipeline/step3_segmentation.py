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

import networkx as nx
import numpy as np
import scipy.sparse as sp

from core.geometry_utils import (
    build_face_adjacency,
    fibonacci_sphere_directions,
    rotation_axis_circle_directions,
)
from core.mesh_types import MeshData, OverhangMap, Region, RegionDependencyGraph


class SegmentationError(Exception):
    """Levee si la segmentation produit un resultat physiquement invalide (cycle de dependances)."""


@dataclass
class SegmentationParams:
    direction_candidate_count: int = 26
    max_angle_from_direction_deg: float = 45.0
    region_count_penalty: float = 0.1
    platform_z: float = 0.0
    platform_contact_tolerance: float = 0.5
    # Si None : utilise fibonacci_sphere_directions (orientation libre, 5 axes).
    # Si fourni : utilise rotation_axis_circle_directions, qui reflete la
    # cinematique reelle d'une machine 4 axes (3 lineaires + 1 rotatif) - un
    # seul DOF de rotation donne un CERCLE de directions atteignables, pas
    # une sphere complete.
    rotation_axis_direction: np.ndarray = None
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
    """
    up = params.up_vector / np.linalg.norm(params.up_vector)

    if params.rotation_axis_direction is not None:
        circle = rotation_axis_circle_directions(
            params.rotation_axis_direction, params.direction_candidate_count
        )
        directions = np.vstack([up.reshape(1, 3), circle])
    else:
        sphere = fibonacci_sphere_directions(params.direction_candidate_count)
        directions = np.vstack([up.reshape(1, 3), sphere])

    return directions


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


# ---------------------------------------------------------------------------
# Etape 4 - Region growing pondere
# ---------------------------------------------------------------------------

def _region_growing(
    face_adjacency: sp.csr_matrix,
    compatibility: np.ndarray,
    risk_score: np.ndarray,
    region_count_penalty: float,
) -> np.ndarray:
    """
    Regroupe les faces en regions par croissance iterative : deux faces
    adjacentes fusionnent si elles partagent au moins une direction
    candidate compatible commune. Les faces a risque eleve sont traitees en
    priorite (elles contraignent le plus fortement le choix de direction).

    Retourne un tableau region_id (M,) - un identifiant de region par face.
    """
    face_count = compatibility.shape[0]
    region_id = np.arange(face_count)  # chaque face demarre comme sa propre region
    # ensemble de directions compatibles restantes par region (intersection
    # cumulative des faces qui la composent)
    region_compat = {i: compatibility[i].copy() for i in range(face_count)}

    adjacency_coo = face_adjacency.tocoo()
    # traite les paires de faces les plus a risque en premier : une fusion
    # impliquant une face a fort risque contraint plus tot le choix de
    # direction, ce qui limite les retours en arriere.
    edge_risk = risk_score[adjacency_coo.row] + risk_score[adjacency_coo.col]
    order = np.argsort(-edge_risk)

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
            # ensemble : fusion refusee, la regularisation anti-fragmentation
            # (region_count_penalty) module indirectement ce refus en amont
            # via l'ordre de traitement (regions a fort risque fusionnees
            # en priorite), pas ici par un terme de cout explicite.
            continue

        # fusionne r2 dans r1
        region_id[r2] = r1
        region_compat[r1] = shared_compat
        del region_compat[r2]

    # normalise les ids finaux (path compression complete + relabeling contigu)
    final_ids = np.array([find(f) for f in range(face_count)])
    unique_ids, relabeled = np.unique(final_ids, return_inverse=True)
    return relabeled


# ---------------------------------------------------------------------------
# Etape 5 - Selection de la direction finale par region
# ---------------------------------------------------------------------------

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

def _build_dependency_graph(
    mesh_data: MeshData,
    region_labels: np.ndarray,
    region_directions: dict[int, np.ndarray],
    is_platform_contact: np.ndarray,
    up_vector: np.ndarray,
) -> RegionDependencyGraph:
    """
    Construit le graphe de dependances entre regions adjacentes. Une region
    A doit etre imprimee avant une region B si B se trouve "au-dessus" de A
    au niveau de leur frontiere commune.

    Important : la projection utilisee pour determiner "au-dessus" se fait
    sur `up_vector` (la verticale machine globale), PAS sur une moyenne
    locale des directions d'impression de chaque paire de regions. Utiliser
    une reference locale par paire n'est pas transitif d'une paire a
    l'autre et peut produire des cycles artificiels des qu'un ensemble de
    regions forme un anneau (typiquement une sphere segmentee en tranches
    autour d'un axe) - chaque paire semble localement coherente, mais
    l'ensemble ne l'est pas. Une reference globale unique garantit un ordre
    partiel valide.
    """
    face_adjacency = build_face_adjacency(mesh_data.faces, mesh_data.face_count)
    adjacency_coo = face_adjacency.tocoo()

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

    region_platform_contact: dict[int, bool] = {}
    for region_id, region in regions.items():
        region_platform_contact[region_id] = bool(
            np.any(is_platform_contact[region.face_ids])
        )

    seen_pairs = set()
    for idx in range(len(adjacency_coo.row)):
        f1, f2 = int(adjacency_coo.row[idx]), int(adjacency_coo.col[idx])
        r1, r2 = int(region_labels[f1]), int(region_labels[f2])
        if r1 == r2:
            continue
        pair = (min(r1, r2), max(r1, r2))
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)

        # heuristique de dependance : la region dont le centre de faces
        # partagees est "en dessous" le long de la verticale machine
        # globale (up_vector) doit etre imprimee en premier. Reference
        # globale et non locale, voir docstring pour la justification.
        center_r1 = face_centers[region_labels == r1].mean(axis=0)
        center_r2 = face_centers[region_labels == r2].mean(axis=0)

        projection_r1 = np.dot(center_r1, up)
        projection_r2 = np.dot(center_r2, up)

        if projection_r1 <= projection_r2:
            graph.add_edge(r1, r2)
        else:
            graph.add_edge(r2, r1)

    # les regions au contact du plateau ne devraient avoir aucune dependance
    # amont (arete entrante) - si le graphe en genere une malgre tout (cas
    # geometrique ambigu), elle est retiree explicitement : le plateau est
    # par definition la reference la plus basse.
    for region_id, is_contact in region_platform_contact.items():
        if is_contact:
            incoming = list(graph.in_edges(region_id))
            graph.remove_edges_from(incoming)

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
    )


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

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
        RegionDependencyGraph avec les regions et leurs dependances,
        garanti acyclique.

    Raises:
        SegmentationError: si le graphe de dependances resultant contient
            un cycle (segmentation physiquement inimprimable).
    """
    is_platform_contact = _compute_platform_contact(
        mesh_data, params.platform_z, params.platform_contact_tolerance
    )
    effective_risk = _effective_risk_score(overhang_map, is_platform_contact)

    directions = _generate_candidate_directions(params)
    compatibility = _face_direction_compatibility(
        mesh_data.face_normals, directions, params.max_angle_from_direction_deg
    )

    face_adjacency = build_face_adjacency(mesh_data.faces, mesh_data.face_count)
    region_labels = _region_growing(
        face_adjacency, compatibility, effective_risk, params.region_count_penalty
    )

    region_directions = _select_region_directions(
        region_labels, compatibility, effective_risk, directions
    )

    return _build_dependency_graph(
        mesh_data, region_labels, region_directions, is_platform_contact, params.up_vector
    )