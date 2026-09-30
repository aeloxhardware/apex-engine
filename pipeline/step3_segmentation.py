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

from dataclasses import dataclass, field, replace
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
    TetMesh,
)
from core.spatial import KDTree
# Reutilise les primitives VOLUMIQUES de step4, jamais reimplementees ici :
# aucun cycle d'import (step4 n'importe rien de step3). Voir
# `_branch_labels_from_tet_mesh` ci-dessous pour pourquoi step3 en a besoin.
from pipeline.step4_field_solver import (
    FieldSolverParams,
    _branch_indicator,
    _build_tet_tet_adjacency,
    _compute_shape_gradients,
    _identify_platform_nodes,
)


class SegmentationError(Exception):
    """Levee si la segmentation produit un resultat physiquement invalide (cycle de dependances)."""


@dataclass
class SegmentationParams:
    # AVEC `branch_cut_threshold` ET `anchor_growth_weight` actifs ensemble,
    # 26 candidats (grille par defaut, fibonacci_sphere_directions) se sont
    # mesures TROP CLAIRSEMES pour representer un axe de croissance
    # presque horizontal (`_anchor_growth_directions`, PCA volumique) : le
    # meilleur candidat compatible avait systematiquement une composante Z
    # parasite. Mesure sur real_overhang.stl (bras horizontal) :
    #   26  candidats -> [-0.04, -0.87, -0.50]  (Z parasite, 50%)
    #   60  candidats -> [+0.17, -0.98, -0.08]  (quasi plat)
    #   120 candidats -> [-0.09, -0.97, -0.24]  (NON MONOTONE, ne pas
    #                                            interpoler entre valeurs)
    # Reste 26 par defaut ici (comportement historique inchange pour tout
    # usage sans segmentation par branche) - `run()` releve ce plancher
    # UNIQUEMENT quand les deux mecanismes sont actifs ensemble (voir
    # `_effective_direction_candidate_count`), jamais silencieusement pour
    # les autres usages de ce parametre.
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
    # l'ancrage plateau du solveur FEM (step4). Voir CLAUDE.md.
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

    # Penalise les directions de croissance qui s'opposent a l'ANCRAGE reel
    # de la region - c'est-a-dire au sens dans lequel sa matiere s'eloigne
    # de la ou elle rejoint le reste de la piece (plateau, ou frontiere avec
    # une region voisine). None (defaut) = DESACTIVE, comportement
    # historique inchange - meme discipline que max_region_merge_angle_deg.
    #
    # Corrige un defaut STRUCTUREL de `_select_region_directions`, mesure sur
    # res/test_models/death.stl : le choix ne minimisait QUE le risque de
    # surplomb residuel, sans jamais considerer par ou la region rejoint la
    # matiere deja imprimee. Pour un plafond de cavite (le capot du Benchy
    # inverse, region 3 : 717 faces, la plus grande du modele), la direction
    # de moindre surplomb est VERS LE BAS ([-0.26, 0.08, -0.96] mesure) -
    # optimale au sens du critere, et pourtant physiquement absurde : la
    # suivre voudrait dire commencer a imprimer en l'air au sommet du capot
    # (z~100) pour descendre, alors que le capot ne tient au pilier que par
    # sa BASE (z~59). Une telle cible fait empiler les couches du haut vers
    # le bas, chacune atterrissant dans le vide.
    #
    # Pourquoi ce critere plutot que `FieldSolverParams.enforce_upward_growth`
    # (step4), qui traite le meme symptome : ce dernier retourne de force
    # TOUTE direction descendante, y compris pour une piece REELLEMENT
    # suspendue par le haut, ou descendre est le bon sens de croissance -
    # exactement le cas d'usage 4 axes que le projet vise. Le critere
    # d'ancrage, lui, accepte une direction descendante quand l'ancrage est
    # en haut, et ne la refuse que lorsqu'elle s'oppose a la croissance
    # reelle. Il agit en amont (sur la CIBLE, pas sur le champ resolu), la ou
    # le defaut est reellement introduit.
    #
    # Valeur = poids de la penalite, en multiples du risque total de la
    # region : une direction parfaitement opposee a la croissance d'ancrage
    # se voit ajouter `anchor_growth_weight * risque_total_region` a son
    # residu, une direction orthogonale ou alignee n'est jamais penalisee.
    anchor_growth_weight: Optional[float] = None

    # --- Segmentation par BRANCHE (structure volumique), pas par face -----
    # `_region_growing` ci-dessous decide par ORIENTATION DE FACE : deux
    # faces fusionnent si une direction candidate reste compatible avec les
    # deux. Mesure sur `res/test_models/real_overhang.stl` (pilier + bras
    # horizontal + pilier suspendu, trois branches distinctes) : ce critere
    # ne les separe JAMAIS, a aucun seuil de `max_region_merge_angle_deg`
    # teste (30, 20, 10, 5 degres) - les parois laterales des trois branches
    # sont toutes a peu pres verticales localement, donc compatibles entre
    # elles au sens de ce critere. La plus grosse region (22/28 faces) reste
    # perpetuellement les trois branches reunies, avec UNE direction unique.
    #
    # Ce parametre ajoute une contrainte de fusion supplementaire, orthogonale
    # a l'angle : deux faces adjacentes ne fusionnent que si elles
    # appartiennent a la MEME composante connexe du volume UNE FOIS COUPE aux
    # zones de ramification (voir `_branch_labels_from_tet_mesh`). La
    # ramification est detectee par la divergence du flux geodesique -
    # |div(X)|, X = gradient normalise de la distance lisse au plateau
    # (methode de la chaleur) - qui vaut ~0 dans un corps de tube (le front
    # de distance n'accelere pas) et s'eleve nettement la ou l'aire du front
    # change d'un coup (un branchement). Mesure sur real_overhang.stl : 3
    # composantes propres (pilier/bras/suspendu) au seuil de coupe adopte.
    #
    # None (defaut) = DESACTIVE, comportement historique inchange - la
    # segmentation reste purement par orientation de face. Nombre plutot
    # qu'un booleen : c'est le SEUIL de coupe de l'indicateur
    # (`branch_cut_threshold`), pas un simple interrupteur - voir cette
    # fonction pour la sensibilite mesuree (0.05 separe les 3 branches,
    # 0.08 les refusionne deja).
    branch_cut_threshold: Optional[float] = None

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


def _merge_tiny_branch_components(
    tet_mesh: TetMesh,
    tet_labels: np.ndarray,
    tet_centers: np.ndarray,
    min_branch_tets: int = 20,
) -> np.ndarray:
    """
    Fusionne toute composante microscopique dans le voisin ADJACENT (au sens
    du graphe tetraedre-tetraedre) dont un tetraedre est le plus proche d'un
    tetraedre de la composante - meme critere que documente historiquement
    (distance point-le-plus-proche-a-point-le-plus-proche, jamais un
    centroide global ni le voisin le plus frequent), mais VECTORISE.

    Bug de performance reel trouve en executant sur `death.stl` (mesh dense,
    2418 faces de surface reelles) : le seuil `branch_cut_threshold` par
    defaut fragmente ce modele en 318 239 micro-composantes (contre
    seulement 3-4 significatives sur `real_overhang.stl`, le modele
    synthetique sur lequel ce mecanisme avait ete mis au point). L'ancienne
    version bouclait en Python sur CHAQUE micro-composante et reconstruisait
    un `cKDTree` par candidat voisin a l'interieur de cette boucle - avec
    318k composantes, ce n'etait tout simplement jamais fini (tue apres
    plusieurs minutes sans terminer, meme a resolution grossiere).

    Fix : un seul passage vectorise par ITERATION plutot qu'un cKDTree par
    composante. A chaque iteration, pour CHAQUE arete de `full_adjacency`
    reliant un tetraedre "petit" (composante < seuil) a un tetraedre "grand"
    (composante >= seuil), calcule la distance euclidienne des deux centres
    en un seul appel numpy vectorise sur TOUTES les aretes candidates a la
    fois ; chaque petite composante adopte alors le label du grand voisin
    connecte par l'arete de distance minimale parmi les siennes. Repete
    jusqu'a stabilisation (les grandes composantes elles-memes ne bougent
    jamais, donc le nombre de petites composantes ne peut que decroitre a
    chaque iteration - boucle bornee par le nombre de petites composantes
    initial, converge en pratique en tres peu de passes car l'immense
    majorite des petites composantes touche directement une grande
    composante des la premiere iteration).
    """
    full_adjacency = _build_tet_tet_adjacency(tet_mesh).tocoo()
    row, col = full_adjacency.row, full_adjacency.col
    edge_dist = np.linalg.norm(tet_centers[row] - tet_centers[col], axis=1)

    label_of_tet = tet_labels.copy()
    for _ in range(64):  # borne dure de securite, jamais atteinte en pratique
        sizes = np.bincount(label_of_tet, minlength=int(label_of_tet.max()) + 1)
        is_tiny_tet = sizes[label_of_tet] < min_branch_tets
        if not np.any(is_tiny_tet):
            break

        row_tiny = is_tiny_tet[row]
        col_tiny = is_tiny_tet[col]
        # arete reliant un tetraedre "petit" a un tetraedre "grand" (jamais
        # petit-petit : une petite composante doit rejoindre une VRAIE
        # branche, pas fusionner avec une autre miette) - direction row->col
        # ET col->row couvertes separement, la matrice d'adjacence etant
        # symetrique mais stockee une seule fois par paire non ordonnee.
        candidate_mask = row_tiny & ~col_tiny
        candidate_mask |= col_tiny & ~row_tiny
        if not np.any(candidate_mask):
            # des petites composantes subsistent mais aucune ne touche
            # directement une grande - ne peut arriver qu'a un ilot
            # entierement entoure d'autres miettes, cas degenere non mesure ;
            # les tetraedres restants gardent leur label, pas de boucle infinie.
            break

        cand_small = np.where(row_tiny[candidate_mask], row[candidate_mask], col[candidate_mask])
        cand_large = np.where(row_tiny[candidate_mask], col[candidate_mask], row[candidate_mask])
        cand_dist = edge_dist[candidate_mask]

        # pour chaque tetraedre "petit" source, garde l'arete de distance
        # minimale vers un tetraedre "grand" (argsort stable + prise du
        # premier par groupe, vectorise plutot qu'une boucle par tetraedre).
        order = np.argsort(cand_dist, kind="stable")
        cand_small_sorted = cand_small[order]
        cand_large_sorted = cand_large[order]
        _, first_occurrence = np.unique(cand_small_sorted, return_index=True)
        best_small = cand_small_sorted[first_occurrence]
        best_large_target = cand_large_sorted[first_occurrence]

        # chaque tetraedre petit adopte le label ACTUEL de son grand voisin
        # le plus proche - propage a tous les tetraedres de sa composante
        # (pas seulement lui-meme) pour que la composante entiere bascule
        # d'un coup plutot que tetraedre par tetraedre.
        new_label_for_small_tet = label_of_tet[best_large_target]
        remap = dict(zip(
            label_of_tet[best_small].tolist(), new_label_for_small_tet.tolist()
        ))
        # une composante ne recoit qu'UNE cible (la premiere rencontree parmi
        # ses membres dans `best_small` suffit : tous ses membres partagent
        # le meme label actuel, donc le meme groupe dans `remap`).
        mask_to_remap = np.isin(label_of_tet, list(remap.keys()))
        if np.any(mask_to_remap):
            old_keys = np.array(list(remap.keys()))
            new_values = np.array(list(remap.values()))
            order_keys = np.argsort(old_keys)
            positions = np.searchsorted(old_keys[order_keys], label_of_tet[mask_to_remap])
            label_of_tet[mask_to_remap] = new_values[order_keys][positions]

    return label_of_tet


def _branch_labels_from_tet_mesh(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    is_platform_contact: np.ndarray,
    cut_threshold: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Etiquette de BRANCHE par face de surface, derivee de la structure
    VOLUMIQUE de la piece - independante de l'orientation des faces (voir
    `SegmentationParams.branch_cut_threshold`).

    Retourne `(face_labels, tet_labels, tet_centers)` : `face_labels` (M,)
    est ce que `_region_growing` consomme (contrainte de fusion par face) ;
    `tet_labels`/`tet_centers` (T,) sont les MEMES etiquettes mais au niveau
    du volume, dans le meme espace d'indices que `face_labels` - destines a
    `_anchor_growth_directions`, qui calcule l'axe principal d'une branche
    sur ses tetraedres plutot que sur ses quelques faces de surface (voir sa
    docstring : sur un maillage grossier, un axe calcule sur les faces peut
    lui-meme etre corrompu par une face mal classee).

    Methode : calcule |div(X)| par tetraedre (`_branch_indicator`, deja
    ecrit pour la modulation de lissage - reutilise ici tel quel, pas
    reimplemente), coupe le graphe d'adjacence tetraedre-tetraedre a chaque
    paire dont l'indicateur depasse `cut_threshold`, prend les composantes
    connexes resultantes. Chaque face de surface herite de l'etiquette du
    tetraedre le plus proche de son centre (meme technique de plus-proche-
    voisin que `_map_boundary_faces_to_surface`, step4, mais dans le sens
    inverse : de tetraedre vers face plutot que de face vers tetraedre).

    Les nœuds plateau sont retrouves directement depuis `is_platform_contact`
    (deja calcule par l'appelant a ce stade de `run()`) plutot que via
    `_identify_platform_nodes` : celle-ci exige
    `RegionDependencyGraph.platform_contact_faces`, qui n'existe qu'APRES la
    segmentation - la region par cette fonction meme. Utiliser cette source
    aurait ete circulaire ; le seuil de z est reconstruit ici a l'identique
    (meme tolerance numerique que `_identify_platform_nodes`).

    Coupure DELIBEREMENT bruitee pres des FRONTIERES du domaine (source
    plateau, extremites libres) - meme faux positif deja documente pour
    `_branch_indicator` (limite de conception, pas un bug) : ces zones
    produisent des milliers de micro-composantes (mesure : 9539 composantes
    au total sur real_overhang.stl, dont seulement 3-4 de taille
    significative). Filtrees ici par une taille minimale absolue plutot
    qu'un ratio, pour rester stable independamment du nombre total de
    tetraedres du modele.
    """
    platform_faces_z = mesh_data.vertices[mesh_data.faces[is_platform_contact]]
    if platform_faces_z.size > 0:
        platform_level = float(platform_faces_z[:, :, 2].max())
    else:
        platform_level = float(tet_mesh.tet_vertices[:, 2].min())
    scale = float(np.ptp(tet_mesh.tet_vertices[:, 2])) or 1.0
    epsilon = 1e-9 * scale
    platform_nodes = np.where(tet_mesh.tet_vertices[:, 2] <= platform_level + epsilon)[0]

    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)
    indicator, _boundary = _branch_indicator(
        tet_mesh, shape_gradients, volumes, platform_nodes, FieldSolverParams()
    )

    adjacency = _build_tet_tet_adjacency(tet_mesh).tocoo()
    edge_indicator = np.maximum(indicator[adjacency.row], indicator[adjacency.col])
    keep = edge_indicator < cut_threshold

    cut_graph = sp.coo_matrix(
        (np.ones(int(keep.sum())), (adjacency.row[keep], adjacency.col[keep])),
        shape=(tet_mesh.tet_count, tet_mesh.tet_count),
    ).tocsr()
    _n_components, tet_labels = spgraph.connected_components(cut_graph, directed=False)
    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)

    # Fusionne toute composante microscopique dans le voisin dont le POINT LE
    # PLUS PROCHE (pas le centroide GLOBAL) est le plus proche de la
    # composante, jamais le voisin le plus FREQUENT en nombre de tetraedres.
    #
    # Bug reel corrige ici, trouve par inspection visuelle du GLB des
    # parties : une micro-composante isolee peut etre une vraie JONCTION
    # (pas seulement un artefact de bord plateau/extremite) - mesure sur
    # real_overhang.stl, 23 tetraedres a y~-2/z~19.5, exactement au coin
    # pilier/bras. La fusionner par frequence de voisinage la rattachait au
    # PLUS GROS voisin (le pilier), donnant au pilier des tetraedres proches
    # du sommet de la piece a z=22 - ce qui faisait ensuite classer a tort,
    # par plus-proche-voisin, une grande face de surface du bras dans la
    # region du pilier.
    #
    # PREMIER correctif tente (centroide GLOBAL a centroide GLOBAL) :
    # MESURE COMME INSUFFISANT ET ECARTE. Le pilier s'etend sur 33mm de
    # hauteur (z de -15 a 18.7) ; son centroide global est tire vers le bas
    # et se retrouve, en distance euclidienne brute, PLUS LOIN de la
    # jonction que le bras alors meme que le HAUT du pilier y touche
    # directement - mesure : dist(jonction, centroide pilier) = 18.7 contre
    # dist(jonction, centroide bras) = 8.2, un rapport qui ne reflete que la
    # taille du pilier, pas sa proximite reelle. Deux nouvelles instabilites
    # sont apparues avec ce critere (6 regions au lieu de 3, une region a une
    # seule face isolee) : le retenir aurait ete pire que l'origine.
    #
    # Le bon critere est la distance MINIMALE (point le plus proche a point
    # le plus proche), independante de la taille ou de la forme globale
    # d'une composante - une jonction rejoint la branche qu'elle TOUCHE
    # reellement, pas celle dont le centre est statistiquement le plus proche.
    tet_labels = _merge_tiny_branch_components(tet_mesh, tet_labels, tet_centers)
    tree = KDTree(tet_centers)

    face_labels = _majority_branch_label_per_face(
        mesh_data.vertices, mesh_data.faces, tree, tet_labels
    )
    unique_labels, relabeled = np.unique(face_labels, return_inverse=True)

    # tet_labels/tet_centers, RENUMEROTES dans le meme espace d'indices que
    # `relabeled` (via le meme `unique_labels`), sont retournes en plus des
    # labels par face : `_anchor_growth_directions` en a besoin pour calculer
    # l'axe principal d'une branche sur le VOLUME plutot que sur les
    # quelques faces de surface qui la representent (voir sa docstring -
    # sur un maillage grossier, l'axe des faces peut lui-meme etre corrompu
    # par des faces mal classees).
    #
    # BUG REEL CORRIGE ICI, trouve par mesure directe (le fix
    # anchor_growth_weight semblait sans effet - la cause etait ici, pas
    # dans le calcul de l'axe principal lui-meme) : `tet_labels` (avant
    # relabeling) peut contenir une etiquette de branche qui n'apparait SUR
    # AUCUNE face de surface - un ilot de tetraedres purement internes,
    # jamais le plus proche d'aucun centre de face. `searchsorted` sur une
    # valeur absente de `unique_labels` renvoie un indice arbitraire (le
    # point d'insertion), corrompant silencieusement le mapping pour TOUS
    # les tetraedres de cet ilot ET, en aval, pour toute etiquette qui suit
    # dans l'ordre trie. Mesure sur real_overhang.stl : tet_labels avait 4
    # valeurs uniques (0,1,2,3) contre 3 pour face_labels, decalant le
    # mapping et faisant echouer silencieusement la PCA volumique.
    #
    # Fix : les tetraedres SANS etiquette correspondante parmi les faces (un
    # ilot purement interne n'a, par definition, aucune utilite pour choisir
    # UNE direction d'impression de SURFACE) recoivent -1, exclu
    # explicitement par `_anchor_growth_directions` (qui replie alors sur les
    # faces) plutot que de leur inventer un rattachement arbitraire.
    label_lookup = {int(old): int(new) for new, old in enumerate(unique_labels)}
    tet_labels_remapped = np.array(
        [label_lookup.get(int(lbl), -1) for lbl in tet_labels], dtype=np.int64
    )
    return relabeled.astype(np.int64), tet_labels_remapped, tet_centers


def _majority_branch_label_per_face(
    vertices: np.ndarray,
    faces: np.ndarray,
    tet_center_tree: KDTree,
    tet_labels: np.ndarray,
) -> np.ndarray:
    """
    Etiquette de branche par face de surface, par VOTE MAJORITAIRE sur des
    points echantillonnes A L'INTERIEUR de chaque triangle - jamais le seul
    centroide de la face.

    Bug reel corrige ici, trouve par inspection visuelle du GLB des parties
    (export_parts_glb) : sur un maillage grossier (real_overhang.stl, 28
    faces au total, fourni par un utilisateur - PAS genere par ce projet,
    donc sa finesse ne peut jamais etre supposee), une face triangulaire
    peut a elle seule couvrir plusieurs branches. Exemple mesure : le
    triangle de sommets (2,0,22)/(-2,0,22)/(2,-20,22) traverse en diagonale
    TOUTE la face superieure de la piece, de l'extremite du bras (y=0)
    jusqu'a celle du pilier suspendu (y=-20). Son CENTROIDE (l'ancien point
    de test unique) tombe dans le bras, alors qu'un tiers du triangle
    appartient reellement au pilier suspendu - la face entiere heritait a
    tort de la branche du bras, produisant des faces de region a
    l'etendue absurde (mesure : region du bras avec des sommets a y=-20,
    alors que le bras lui-meme s'arrete a y=-16) et donc des directions
    cibles incoherentes.

    Echantillonnage barycentrique DENSITE PROPORTIONNELLE A L'AIRE (pas un
    nombre fixe de points) : une face 100x plus grande qu'une autre doit
    recevoir plus de points, sous peine de sous-echantillonner exactement
    les grandes faces qui posent le probleme. Densite calibree sur la
    taille caracteristique des tetraedres (aire d'un triangle equilateral
    de cote = arete mediane du TetMesh), pas un nombre absolu qui n'aurait
    aucun sens d'un modele a l'autre.
    """
    triangles = vertices[faces]  # (M,3,3)
    edges_a = triangles[:, 1] - triangles[:, 0]
    edges_b = triangles[:, 2] - triangles[:, 0]
    areas = 0.5 * np.linalg.norm(np.cross(edges_a, edges_b), axis=1)

    # Taille de reference : aire d'un triangle equilateral dont le cote est
    # l'arete mediane des tetraedres du domaine - meme echelle que le
    # maillage volumique auquel on compare, jamais une constante absolue.
    tet_points = tet_center_tree.data
    if len(tet_points) > 1:
        sample_edges = tet_points[: min(200, len(tet_points))]
        median_spacing = float(np.median(
            scipy.spatial.distance.pdist(sample_edges)
        )) if len(sample_edges) > 1 else 1.0
    else:
        median_spacing = 1.0
    reference_area = max((np.sqrt(3.0) / 4.0) * median_spacing ** 2, 1e-9)

    labels = np.zeros(len(faces), dtype=np.int64)
    rng = np.random.default_rng(0)  # deterministe

    for i in range(len(faces)):
        v0, v1, v2 = triangles[i]
        sample_count = max(4, int(np.ceil(areas[i] / reference_area)) + 1)
        # Coordonnees barycentriques uniformes (methode standard : deux
        # tirages uniformes, replies sur le triangle si leur somme depasse 1).
        u = rng.uniform(0.0, 1.0, sample_count)
        v = rng.uniform(0.0, 1.0, sample_count)
        over = u + v > 1.0
        u[over] = 1.0 - u[over]
        v[over] = 1.0 - v[over]
        points = v0 + u[:, None] * (v1 - v0) + v[:, None] * (v2 - v0)
        # Sommets et centroide inclus explicitement : garantit une bonne
        # couverture meme pour sample_count minimal.
        points = np.vstack([points, triangles[i], triangles[i].mean(axis=0)])

        _, nearest_tets = tet_center_tree.query(points)
        sample_labels = tet_labels[nearest_tets]
        labels[i] = np.bincount(sample_labels).argmax()

    return labels


def _region_growing(
    face_adjacency: sp.csr_matrix,
    compatibility: np.ndarray,
    risk_score: np.ndarray,
    directions: np.ndarray,
    max_region_merge_angle_deg: Optional[float],
    branch_labels: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Regroupe les faces en regions par croissance iterative : deux faces
    adjacentes fusionnent si elles partagent au moins une direction
    candidate compatible commune. Les faces a risque eleve sont traitees en
    priorite (elles contraignent le plus fortement le choix de direction).

    `branch_labels` (optionnel, un entier par face) : contrainte de fusion
    SUPPLEMENTAIRE et orthogonale a l'angle - deux faces de branches
    differentes (voir `SegmentationParams.branch_cut_threshold`) ne
    fusionnent jamais, quelle que soit leur compatibilite de direction.
    Corrige un defaut mesure du critere par angle seul : sur une piece a
    parois lateralement toutes verticales (pilier/bras/pilier suspendu), la
    compatibilite d'orientation ne separe jamais les branches, a aucun
    seuil - `branch_labels` regarde la structure du VOLUME plutot que
    l'orientation des faces.

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

        if branch_labels is not None and branch_labels[f1] != branch_labels[f2]:
            # Faces de branches VOLUMIQUES differentes : fusion refusee sans
            # meme regarder la compatibilite de direction. C'est la seule
            # facon de separer real_overhang.stl (voir docstring).
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


def _anchor_growth_directions(
    region_labels: np.ndarray,
    face_centers: np.ndarray,
    is_platform_contact: np.ndarray,
    face_adjacency,
    depths: dict[int, int],
    up: np.ndarray,
    branch_tet_labels: Optional[np.ndarray] = None,
    branch_tet_centers: Optional[np.ndarray] = None,
) -> dict[int, np.ndarray]:
    """
    Sens de croissance NATUREL de chaque region : de son ANCRAGE vers le
    reste de sa matiere.

    L'ancrage d'une region, c'est la ou elle rejoint ce qui la PORTE :
      - ses faces au contact du plateau, si elle en a (c'est une racine) ;
      - sinon, ses faces frontieres avec une region STRICTEMENT MOINS
        PROFONDE - son ou ses parents dans l'arborescence de croissance, par
        lesquels seuls la matiere peut lui parvenir ;
      - a defaut (region isolee, ou sans parent identifiable), `up` en repli.

    La direction rendue est `centroide(region) - centroide(ancrage)`,
    normalisee : elle pointe du point d'attache vers le volume a construire.
    C'est la question physiquement pertinente - "dans quel sens cette
    matiere s'eloigne-t-elle de ce qui la porte" - et NON une comparaison de
    hauteur, qui se tromperait sur toute piece suspendue (plus basse que son
    support, mais devant croitre vers le bas depuis lui).

    Le filtrage sur la PROFONDEUR est essentiel, et non un raffinement :
    prendre toutes les frontieres sans distinguer parents et enfants a ete
    mesure comme AGGRAVANT le probleme sur death.stl (4 -> 6 regions a
    direction descendante) - une region bordee surtout par ses enfants
    voyait son "ancrage" place au-dessus d'elle, et sa croissance donc
    pointee vers le bas. `_compute_print_depths` ne depend d'aucune
    direction (uniquement adjacence, racines plateau et centres de faces),
    il n'y a donc aucune circularite a s'en servir ici.
    """
    growth: dict[int, np.ndarray] = {}
    up_unit = up / (np.linalg.norm(up) or 1.0)

    coo = face_adjacency.tocoo()
    labels_row = region_labels[coo.row]
    labels_col = region_labels[coo.col]
    cross = labels_row != labels_col

    # Faces de r touchant une region STRICTEMENT moins profonde que r.
    parent_boundary: dict[int, set[int]] = {}
    for face_id, own_region, other_region in zip(
        coo.row[cross], labels_row[cross], labels_col[cross]
    ):
        own_id, other_id = int(own_region), int(other_region)
        if depths.get(other_id, 0) < depths.get(own_id, 0):
            parent_boundary.setdefault(own_id, set()).add(int(face_id))

    for region in np.unique(region_labels):
        region_id = int(region)
        region_face_ids = np.where(region_labels == region)[0]

        platform_faces = region_face_ids[is_platform_contact[region_face_ids]]
        if len(platform_faces) > 0:
            anchor_faces = platform_faces
        else:
            anchor_faces = np.array(
                sorted(parent_boundary.get(region_id, set())), dtype=np.int64
            )

        if len(anchor_faces) == 0:
            growth[region_id] = up_unit
            continue

        offset = face_centers[region_face_ids].mean(axis=0) - face_centers[anchor_faces].mean(axis=0)
        norm = float(np.linalg.norm(offset))
        # Ancrage confondu avec le centre de la region (region reduite a sa
        # propre frontiere) : aucun sens de croissance lisible, repli sur up
        # plutot qu'une direction issue d'un vecteur nul normalise.
        if norm <= 1e-9:
            growth[region_id] = up_unit
            continue

        # PROJECTION sur l'axe d'inertie PRINCIPAL de la region, plutot que
        # le vecteur brut centroide-a-centroide tel quel.
        #
        # Bug reel corrige ici, trouve par inspection visuelle
        # (export_parts_glb) : sur un maillage grossier (real_overhang.stl,
        # 28 faces), les quelques faces d'ANCRAGE ne sont pas forcement
        # representatives du reste de la region - mesure sur le bras de ce
        # modele, ancre par 4 faces de z moyen 20.33 (le HAUT de la
        # jonction), alors que le centroide de TOUTE la region (qui inclut
        # aussi le dessous du bras) tombe a z=18.70, legerement plus bas.
        # Le vecteur brut portait donc une composante Z parasite (-0.35)
        # sans rapport avec l'orientation reelle du bras (horizontal, z
        # quasi constant sur toute sa longueur) : artefact d'echantillonnage
        # de quelques faces de bord, pas la vraie direction de la branche.
        #
        # L'axe principal (premiere composante de la SVD des centres de
        # faces recentres) capture la direction ou la region s'ETEND
        # REELLEMENT, independamment de l'endroit precis ou son ancrage a
        # ete echantillonne. Le SENS reste celui du vecteur brut (l'axe
        # principal n'a pas de sens propre, seulement une orientation) : on
        # projette `offset` dessus puis on prend le signe qui va dans le
        # meme sens que `offset`, jamais l'inverse.
        #
        # PCA sur les TETRAEDRES de la branche quand ils sont disponibles
        # (segmentation par branche active), PAS sur les faces de surface.
        # Bug reel corrige ici : sur real_overhang.stl (28 faces), l'axe
        # calcule sur les faces du bras restait lui-meme corrompu par des
        # faces mal classees appartenant en realite au pilier suspendu
        # (mesure : axe [0.09, 0.85, 0.52], une composante Z tout aussi
        # parasite que celle qu'on cherchait a corriger). Les tetraedres
        # d'une branche volumique, bien plus nombreux (des milliers contre
        # une poignee de faces), donnent un axe fiable independamment de la
        # grossierete du maillage de SURFACE.
        pca_points = None
        if branch_tet_labels is not None and branch_tet_centers is not None:
            # La region_id EST l'etiquette de branche quand la segmentation
            # par branche est active (une region = une composante connexe du
            # graphe tetraedre-tetraedre coupe) : meme espace d'indices,
            # verifie par construction dans `_branch_labels_from_tet_mesh`
            # (relabeling contigu partage entre face_labels et tet_labels).
            tet_mask = branch_tet_labels == region_id
            if np.any(tet_mask):
                pca_points = branch_tet_centers[tet_mask]

        if pca_points is None:
            pca_points = face_centers[region_face_ids]

        centred = pca_points - pca_points.mean(axis=0)
        if centred.shape[0] >= 2:
            _u, _s, vt = np.linalg.svd(centred, full_matrices=False)
            principal_axis = vt[0]
            projection = float(np.dot(offset, principal_axis))
            if abs(projection) > 1e-9:
                aligned_axis = principal_axis if projection > 0 else -principal_axis
                growth[region_id] = aligned_axis
                continue

        growth[region_id] = offset / norm

    return growth


def _select_region_directions(
    region_labels: np.ndarray,
    compatibility: np.ndarray,
    risk_score: np.ndarray,
    directions: np.ndarray,
    anchor_growth: Optional[dict[int, np.ndarray]] = None,
    anchor_growth_weight: Optional[float] = None,
) -> dict[int, np.ndarray]:
    """
    Pour chaque region, choisit parmi les directions candidates compatibles
    avec l'ensemble de ses faces celle qui minimise le risque residuel
    total. Si aucune direction n'est compatible avec toutes les faces de la
    region (ne devrait pas arriver si _region_growing a fonctionne
    correctement, mais gere en filet de securite), choisit la direction qui
    minimise le risque residuel meme partiellement.

    `anchor_growth`/`anchor_growth_weight` (opt-in, voir
    `SegmentationParams.anchor_growth_weight`) ajoutent au residu une
    penalite proportionnelle a l'OPPOSITION entre la direction candidate et
    le sens de croissance depuis l'ancrage de la region. Sans elle, le seul
    critere est le surplomb residuel - qui, pour un plafond de cavite, elit
    une direction descendante physiquement irrealisable (il faudrait
    commencer a imprimer en l'air). La penalite ne borne rien et n'interdit
    rien : une direction opposee reste elue si elle reste, meme penalisee,
    la moins risquee - le mecanisme corrige un critere incomplet, il ne le
    remplace pas.
    """
    region_direction: dict[int, np.ndarray] = {}
    use_anchor = anchor_growth is not None and anchor_growth_weight is not None

    for region in np.unique(region_labels):
        face_mask = region_labels == region
        region_compat = compatibility[face_mask]  # (n_faces_in_region, D)
        region_risk = risk_score[face_mask]

        # risque residuel par direction candidate = somme des risques des
        # faces NON compatibles avec cette direction
        incompatible = ~region_compat  # (n_faces_in_region, D)
        residual_risk = (incompatible * region_risk[:, None]).sum(axis=0)  # (D,)

        if use_anchor:
            growth_direction = anchor_growth.get(int(region))
            if growth_direction is not None:
                # Penalite CONTINUE, croissante avec le desalignement :
                # 0 pour une direction alignee sur la croissance, 0.5 pour
                # une orthogonale, 1 pour une exactement opposee. Un simple
                # `max(0, -alignement)` a ete mesure insuffisant : il rend
                # TOUTES les directions non-opposees stricement equivalentes
                # (penalite nulle), si bien que sur une region sans risque -
                # ou le residu est deja nul partout - l'argmin restait
                # arbitraire et pouvait elire une direction quasi-opposee au
                # sens de croissance (mesure sur death.stl : 3 regions du
                # sommet basculaient ainsi de montantes a descendantes).
                misalignment = (1.0 - (directions @ growth_direction)) * 0.5
                # Echelle = risque total de la region, pour que la penalite
                # reste commensurable au residu qu'elle corrige, quelle que
                # soit la taille de la region ou l'amplitude des risques.
                scale = float(region_risk.sum())
                # Region sans aucun risque (typiquement posee a plat sur le
                # plateau) : le residu est nul pour toutes les directions et
                # l'argmin est arbitraire. Une echelle unitaire laisse alors
                # l'ancrage trancher, plutot que le premier indice venu.
                if scale <= 0.0:
                    scale = 1.0
                residual_risk = residual_risk + anchor_growth_weight * scale * misalignment

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

def run(
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    params: SegmentationParams,
    tet_mesh: Optional[TetMesh] = None,
) -> RegionDependencyGraph:
    """
    Segmente le mesh en sous-volumes partageant une direction d'impression
    viable, et construit le graphe de dependances physiques entre eux.

    Args:
        mesh_data: sortie de pipeline/step1_ingest.py.
        overhang_map: sortie de pipeline/step2_overhang.py.
        params: seuils et parametres de segmentation.
        tet_mesh: sortie de pipeline/step1_ingest.py (le meme TetMesh que
            step4 consommera). REQUIS uniquement si
            `params.branch_cut_threshold is not None` - la segmentation par
            branche a besoin de la structure VOLUMIQUE, indisponible sur le
            seul maillage de surface. `None` (defaut) : `branch_cut_threshold`
            doit rester `None` sous peine de `SegmentationError` explicite,
            plutot que d'ignorer silencieusement un parametre demande.

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

    # Voir SegmentationParams.direction_candidate_count : la grille par
    # defaut (26 candidats) s'est mesuree trop clairsemee pour representer
    # un axe de croissance presque horizontal quand `branch_cut_threshold`
    # ET `anchor_growth_weight` sont actifs ensemble. Releve UNIQUEMENT dans
    # ce cas, jamais pour un usage qui n'exerce pas les deux mecanismes -
    # `params` lui-meme n'est jamais mute, seule cette copie locale l'est.
    effective_params = params
    if (
        params.branch_cut_threshold is not None
        and params.anchor_growth_weight is not None
        and params.direction_candidate_count < 60
    ):
        effective_params = replace(params, direction_candidate_count=60)

    directions = _generate_candidate_directions(effective_params)
    compatibility = _face_direction_compatibility(
        mesh_data.face_normals, directions, params.max_angle_from_direction_deg
    )

    branch_labels = None
    branch_tet_labels = None
    branch_tet_centers = None
    if params.branch_cut_threshold is not None:
        if tet_mesh is None:
            raise SegmentationError(
                "SegmentationParams.branch_cut_threshold est renseigne mais "
                "tet_mesh est absent : la segmentation par branche a besoin "
                "du TetMesh (sortie de step1_ingest.run()). Passez-le a "
                "segmentation_run(), ou remettez branch_cut_threshold a None."
            )
        branch_labels, branch_tet_labels, branch_tet_centers = _branch_labels_from_tet_mesh(
            tet_mesh, mesh_data, is_platform_contact, params.branch_cut_threshold
        )

    face_adjacency = build_face_adjacency(mesh_data.faces, mesh_data.face_count)
    region_labels = _region_growing(
        face_adjacency, compatibility, effective_risk, directions,
        params.max_region_merge_angle_deg, branch_labels,
    )

    anchor_growth = None
    if params.anchor_growth_weight is not None:
        # Profondeurs de croissance depuis le plateau - calculees ICI, avant
        # le choix des directions, car elles n'en dependent pas (uniquement
        # adjacence, racines plateau et centres de faces). Les MEMES
        # primitives servent ensuite a `_build_dependency_graph` : une seule
        # source de verite sur "qui porte qui", jamais deux calculs paralleles
        # qui pourraient diverger en silence.
        face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
        up_unit = params.up_vector / np.linalg.norm(params.up_vector)
        anchor_growth = _anchor_growth_directions(
            region_labels,
            face_centers,
            is_platform_contact,
            face_adjacency,
            _compute_print_depths(
                _region_adjacency_graph(region_labels, face_adjacency),
                _platform_root_region_ids(region_labels, is_platform_contact),
                region_labels,
                face_centers,
                up_unit,
            ),
            params.up_vector,
            branch_tet_labels,
            branch_tet_centers,
        )

    region_directions = _select_region_directions(
        region_labels,
        compatibility,
        effective_risk,
        directions,
        anchor_growth,
        params.anchor_growth_weight,
    )

    return _build_dependency_graph(
        mesh_data, region_labels, region_directions, is_platform_contact, params.up_vector
    )