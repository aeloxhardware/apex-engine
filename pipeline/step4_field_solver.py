"""
pipeline/step4_field_solver.py

Etape 4 du pipeline : assemblage et resolution du champ scalaire continu sur
le TetMesh, dont les isosurfaces (etape 5) definiront les couches courbes.

Depend de TetMesh + MeshData (step1), OverhangMap (step2),
RegionDependencyGraph (step3), et des Constraint[] envoyees dynamiquement
par l'UI (points d'attache, contraintes esthetiques/structurelles).

Le systeme assemble est A.phi = b, symetrique defini positif, somme de :
  - un terme de lissage (rigidite FEM P1) ;
  - un terme de surplomb (le gradient doit suivre la direction d'impression
    de chaque region, pres des faces a risque) ;
  - un ancrage automatique du plateau (leve la singularite du systeme) ;
  - les contraintes utilisateur, assemblees par core/constraints.py.

Structure en trois fonctions (cf. plan, decision 6) : `build_static_system()`
assemble tout ce qui NE depend PAS des contraintes,
`solve_with_constraints()` ajoute les contraintes et resout, et `run()` est
l'enveloppe des deux pour un slicing complet en une passe. Ce decoupage
permet a l'orchestrateur de re-resoudre rapidement apres l'edition d'un
point d'attache, sans reassembler la partie structurelle ni relancer
step1/2/3.

Sortie : ScalarField
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.spatial

from core.constraints import ConstraintAssemblyContext, assemble_all_constraints
from core.mesh_types import (
    Constraint,
    MeshData,
    OverhangMap,
    RegionDependencyGraph,
    ScalarField,
    TetMesh,
)


class FieldSolverError(Exception):
    """Levee si le systeme assemble n'est pas exploitable ou si le solveur ne converge pas."""


@dataclass
class FieldSolverParams:
    laplacian_smoothness_weight: float = 1.0
    gradient_overhang_weight: float = 5.0
    resolution: float = 1.5
    max_solver_iterations: int = 2000
    convergence_tolerance: float = 1e-8

    # Poids de base du terme "le gradient suit la direction d'impression de
    # la region", applique a TOUS les tetraedres. Sans lui, le champ ne
    # serait pilote que pres des surplombs et resterait indetermine (donc
    # nul, l'ancrage plateau valant 0) dans le reste du volume : une piece
    # sans surplomb n'aurait aucune couche. C'est ce terme qui produit le
    # slicing plan classique par defaut ; gradient_overhang_weight ne fait
    # que le renforcer la ou le risque est eleve.
    gradient_baseline_weight: float = 1.0

    # Nombre de sauts d'adjacence tetraedre-tetraedre sur lesquels propager
    # l'influence d'une face a risque. 1 = seulement les tetraedres au
    # contact direct de la face frontiere concernee.
    overhang_influence_depth: int = 1

    # Ancrage automatique des regions racines (contact plateau) a
    # platform_anchor_value. Necessaire pour lever la singularite du
    # systeme (un Laplacien seul a un noyau constant : phi = c est toujours
    # solution) quand aucune Constraint[] utilisateur n'est fournie.
    enable_platform_anchor: bool = True
    platform_anchor_weight: float = 1e6
    platform_anchor_value: float = 0.0

    # Seuil (nombre de noeuds du TetMesh) en-dessous duquel un solveur
    # direct (spsolve) est utilise plutot qu'un solveur iteratif (cg).
    # Le remplissage d'une factorisation directe 3D croit vite, d'ou un
    # seuil bien inferieur a ce qu'on tolererait en 2D.
    direct_solve_node_threshold: int = 30000
    use_ilu_preconditioner: bool = False


@dataclass
class FieldSolverStaticSystem:
    """
    Partie du systeme lineaire qui NE depend PAS des Constraint[] : elle
    reste valide tant que tet_mesh, mesh_data, overhang_map,
    region_dependency_graph et les poids/seuils de FieldSolverParams sont
    inchanges.

    L'orchestrateur la construit une fois apres step3, la met en cache, et
    la reutilise a chaque edition de point d'attache via
    `solve_with_constraints()`. Si un parametre non lie aux contraintes
    change (poids de lissage, de surplomb, d'ancrage, seuils solveur), le
    cache doit etre invalide et `build_static_system()` rappele.
    """
    A_structural: sp.csr_matrix
    b_structural: np.ndarray
    ctx: ConstraintAssemblyContext


# ---------------------------------------------------------------------------
# Decision 1 - Matrice de rigidite FEM P1 (Laplacien discret geometrique)
# ---------------------------------------------------------------------------

def _compute_shape_gradients(tet_mesh: TetMesh) -> tuple[np.ndarray, np.ndarray]:
    """
    Pour chaque tetraedre, calcule G (3,4) - le gradient constant de chacune
    des 4 fonctions de forme barycentriques - et son volume.

    Methode standard P1 : les colonnes de l'inverse transpose de la matrice
    des aretes (x1-x0, x2-x0, x3-x0) donnent les gradients des fonctions de
    forme 1..3 ; le gradient de la fonction de forme 0 est l'oppose de leur
    somme (partition de l'unite : la somme des 4 fonctions vaut 1 partout,
    donc la somme de leurs gradients est nulle).

    Retourne (shape_gradients (T,3,4), volumes (T,)). Les tetraedres
    degeneres (volume ~ 0) ont un gradient nul et un volume nul : leur
    contribution disparait naturellement des assemblages.
    """
    tet_count = tet_mesh.tet_count
    shape_gradients = np.zeros((tet_count, 3, 4))
    volumes = np.zeros(tet_count)

    vertices_per_tet = tet_mesh.tet_vertices[tet_mesh.tetrahedra]  # (T,4,3)
    edge_matrices = np.transpose(
        vertices_per_tet[:, 1:, :] - vertices_per_tet[:, 0:1, :], (0, 2, 1)
    )  # (T,3,3), colonnes = x1-x0, x2-x0, x3-x0

    determinants = np.linalg.det(edge_matrices)
    signed_volumes = determinants / 6.0
    valid = np.abs(signed_volumes) > 1e-14

    if np.any(valid):
        inverses = np.linalg.inv(edge_matrices[valid])          # (V,3,3)
        grad_1_3 = np.transpose(inverses, (0, 2, 1))            # (V,3,3)
        grad_0 = -grad_1_3.sum(axis=2, keepdims=True)           # (V,3,1)
        shape_gradients[valid] = np.concatenate([grad_0, grad_1_3], axis=2)
        volumes[valid] = np.abs(signed_volumes[valid])

    return shape_gradients, volumes


def _assemble_stiffness_matrix(
    tet_mesh: TetMesh, shape_gradients: np.ndarray, volumes: np.ndarray
) -> sp.csr_matrix:
    """
    Assemble la matrice de rigidite FEM P1 globale (Laplacien discret
    geometriquement correct) : K_e = V_e * G_e^T G_e par element, sommee
    aux indices de sommets globaux. Symetrique semi-definie positive, de
    noyau les champs constants.
    """
    local_matrices = volumes[:, None, None] * np.einsum(
        "tki,tkj->tij", shape_gradients, shape_gradients
    )  # (T,4,4)
    return _scatter_local_matrices(tet_mesh, np.arange(tet_mesh.tet_count), local_matrices)


def _scatter_local_matrices(
    tet_mesh: TetMesh, tet_indices: np.ndarray, local_matrices: np.ndarray
) -> sp.csr_matrix:
    """
    Disperse des matrices locales (n,4,4) vers la matrice creuse globale,
    aux indices de sommets des tetraedres `tet_indices`. Les doublons
    d'indices sont sommes par coo_matrix, ce qui realise l'assemblage FEM.
    """
    node_count = tet_mesh.node_count
    if len(tet_indices) == 0:
        return sp.csr_matrix((node_count, node_count))

    tets = tet_mesh.tetrahedra[tet_indices]                     # (n,4)
    rows = np.repeat(tets, 4, axis=1).ravel()                   # (n*16,)
    cols = np.tile(tets, (1, 4)).ravel()                        # (n*16,)
    data = local_matrices.reshape(-1)

    return sp.coo_matrix(
        (data, (rows, cols)), shape=(node_count, node_count)
    ).tocsr()


def _compute_per_tet_gradient(
    tet_mesh: TetMesh, shape_gradients: np.ndarray, values: np.ndarray
) -> np.ndarray:
    """
    Gradient constant par tetraedre du champ P1 defini par `values` (une
    valeur par noeud) : grad(phi)|_e = G_e @ phi_e. Representation exacte
    d'un champ lineaire par element, pas une approximation.
    """
    values_per_tet = values[tet_mesh.tetrahedra]                # (T,4)
    return np.einsum("tij,tj->ti", shape_gradients, values_per_tet)


# ---------------------------------------------------------------------------
# Decision 2 - Pont face de surface -> tetraedre volumique
# ---------------------------------------------------------------------------

def _tet_face_table(tet_mesh: TetMesh) -> dict[tuple[int, int, int], list[int]]:
    """
    Table face (triplet trie de noeuds) -> liste des tetraedres proprietaires.
    Une face interne a 2 proprietaires, une face frontiere un seul.
    """
    face_to_tets: dict[tuple[int, int, int], list[int]] = {}
    local_faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))

    for t_idx, tet in enumerate(tet_mesh.tetrahedra):
        for local_face in local_faces:
            key = tuple(sorted(int(tet[i]) for i in local_face))
            face_to_tets.setdefault(key, []).append(t_idx)

    return face_to_tets


def _map_boundary_faces_to_surface(
    tet_mesh: TetMesh, mesh_data: MeshData
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Identifie les faces frontiere du TetMesh (celles n'appartenant qu'a un
    seul tetraedre) et, pour chacune, retrouve la face du mesh de surface
    d'origine la plus proche via MeshData.kdtree (construit a l'etape 1 sur
    les centres de faces).

    Un mapping par index direct serait fragile : tetgen est appele sans le
    switch `Y`, donc il peut inserer des points de Steiner sur la frontiere
    et celle-ci n'est pas identique noeud-a-noeud au mesh source.

    Retourne (boundary_face_nodes (K,3), boundary_face_tets (K,),
    nearest_surface_face_id (K,)).
    """
    face_to_tets = _tet_face_table(tet_mesh)

    boundary_nodes = []
    boundary_tets = []
    for key, owners in face_to_tets.items():
        if len(owners) == 1:
            boundary_nodes.append(key)
            boundary_tets.append(owners[0])

    if not boundary_nodes:
        empty_int = np.zeros((0, 3), dtype=np.int64)
        return empty_int, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)

    boundary_face_nodes = np.array(boundary_nodes, dtype=np.int64)
    boundary_face_tets = np.array(boundary_tets, dtype=np.int64)

    centroids = tet_mesh.tet_vertices[boundary_face_nodes].mean(axis=1)
    _, nearest_surface_face_id = mesh_data.kdtree.query(centroids)

    return boundary_face_nodes, boundary_face_tets, np.asarray(nearest_surface_face_id, dtype=np.int64)


def _build_tet_tet_adjacency(tet_mesh: TetMesh) -> sp.csr_matrix:
    """
    Matrice d'adjacence tetraedre-tetraedre (deux tetraedres sont adjacents
    s'ils partagent une face). Distincte de TetMesh.tet_adjacency, qui est
    une adjacence noeud-noeud : ici on propage un masque defini PAR
    TETRAEDRE, pas par noeud.
    """
    tet_count = tet_mesh.tet_count
    face_to_tets = _tet_face_table(tet_mesh)

    rows, cols = [], []
    for owners in face_to_tets.values():
        if len(owners) == 2:
            t1, t2 = owners
            rows.extend([t1, t2])
            cols.extend([t2, t1])

    if not rows:
        return sp.csr_matrix((tet_count, tet_count))

    data = np.ones(len(rows))
    return sp.csr_matrix((data, (rows, cols)), shape=(tet_count, tet_count))


# ---------------------------------------------------------------------------
# Decision 3 - Terme de direction de gradient (base + sur-ponderation surplomb)
# ---------------------------------------------------------------------------

def _spread_risk_influence(
    tet_tet_adjacency: sp.csr_matrix, seed_tets: np.ndarray, hops: int
) -> np.ndarray:
    """
    Propage un masque booleen de tetraedres "a risque" sur `hops` sauts
    d'adjacence tetraedre-tetraedre, en partant de `seed_tets`.
    """
    tet_count = tet_tet_adjacency.shape[0]
    mask = np.zeros(tet_count, dtype=bool)
    if len(seed_tets) == 0:
        return mask
    mask[seed_tets] = True

    frontier = mask.copy()
    for _ in range(max(0, hops)):
        reached = (tet_tet_adjacency @ frontier.astype(np.float64)) > 0
        new_frontier = reached & (~mask)
        if not np.any(new_frontier):
            break
        mask |= new_frontier
        frontier = new_frontier

    return mask


def _assemble_gradient_direction_term(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    active_tets: np.ndarray,
    tet_directions: np.ndarray,
    tet_weights: np.ndarray,
) -> tuple[sp.csr_matrix, np.ndarray]:
    """
    Aligne le gradient du champ sur la direction d'impression de la region,
    via l'energie de correspondance de gradient
    E_e = w_e * V_e * ||grad(phi) - d||^2. Sa minimisation donne :
        A_local = w_e * V_e * G_e^T G_e
        b_local = w_e * V_e * G_e^T d

    Le terme en b est la force motrice du champ : une penalite purement
    homogene (sur la seule composante perpendiculaire, G_e^T (I - d d^T) G_e)
    n'injecterait rien dans le second membre et laisserait le systeme
    s'effondrer sur phi = 0 partout en l'absence de contrainte de valeur.

    Applique a tous les tetraedres de `active_tets`, avec un poids par
    tetraedre : poids de base partout (le champ doit definir des couches
    dans tout le volume, pas seulement pres des surplombs), majore autour
    des zones a risque.
    """
    node_count = tet_mesh.node_count
    active = np.where(active_tets & (tet_weights > 0.0) & (volumes > 0.0))[0]
    if len(active) == 0:
        return sp.csr_matrix((node_count, node_count)), np.zeros(node_count)

    gradients = shape_gradients[active]                          # (n,3,4)
    directions = tet_directions[active]                          # (n,3)
    norms = np.linalg.norm(directions, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    directions = directions / norms

    scale = tet_weights[active] * volumes[active]                # (n,)

    local_matrices = scale[:, None, None] * np.einsum(
        "tki,tkj->tij", gradients, gradients
    )                                                            # (n,4,4)
    local_rhs = scale[:, None] * np.einsum("tki,tk->ti", gradients, directions)  # (n,4)

    A = _scatter_local_matrices(tet_mesh, active, local_matrices)

    b = np.zeros(node_count)
    np.add.at(b, tet_mesh.tetrahedra[active].ravel(), local_rhs.ravel())

    return A, b


# ---------------------------------------------------------------------------
# Decision 4 - Ancrage automatique du plateau
# ---------------------------------------------------------------------------

def _identify_platform_nodes(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    region_dependency_graph: RegionDependencyGraph,
) -> np.ndarray:
    """
    Retrouve les noeuds du TetMesh reposant sur le plateau.

    `RegionDependencyGraph.platform_contact_faces` (calcule par step3, seule
    etape qui connait la position du plateau) designe QUELLES faces touchent
    le plateau ; le niveau geometrique du plateau en est deduit directement
    (z maximal de leurs sommets), puis les noeuds du TetMesh sont testes
    contre ce plan. Aucun parametre `platform_z` n'est redefini ici : la
    source de verite reste step3.

    Deux pieges evites :
      - NE PAS deduire le plateau de la topologie du graphe ("regions sans
        arete entrante") : des que la segmentation ne produit qu'une seule
        region - cas de toute forme simple, cube compris - cette region est
        a la fois racine et couvre tout le modele, ce qui ancrerait toute la
        surface a 0 et empecherait le champ de croitre.
      - NE PAS passer par la face de surface la plus proche de chaque face
        frontiere du TetMesh : ce mapping (correct pour heriter d'une
        direction de region, ou une petite erreur est lissee) est trop
        grossier pour un predicat binaire. Mesure sur un cube : il ancrait
        des noeuds de paroi laterale jusqu'a z=-1.82 tout en manquant de
        vrais noeuds de base a z=-5.
    """
    platform_faces = region_dependency_graph.platform_contact_faces
    if platform_faces is None:
        raise FieldSolverError(
            "RegionDependencyGraph.platform_contact_faces est absent : "
            "impossible d'ancrer le champ au plateau. Utilisez la sortie de "
            "pipeline/step3_segmentation.run(), ou desactivez explicitement "
            "l'ancrage via FieldSolverParams(enable_platform_anchor=False)."
        )

    platform_faces = np.asarray(platform_faces, dtype=bool)
    if not np.any(platform_faces):
        return np.zeros(0, dtype=np.int64)

    # niveau du plateau = sommet le plus haut des faces en contact (la bande
    # de tolerance de step3 est donc deja incluse dans cette valeur)
    platform_vertex_z = mesh_data.vertices[mesh_data.faces[platform_faces]][:, :, 2]
    platform_level = float(platform_vertex_z.max())

    scale = float(np.ptp(tet_mesh.tet_vertices[:, 2])) or 1.0
    epsilon = 1e-9 * scale

    return np.where(tet_mesh.tet_vertices[:, 2] <= platform_level + epsilon)[0]


def _assemble_platform_anchor(
    node_count: int, platform_nodes: np.ndarray, weight: float, anchor_value: float
) -> tuple[sp.csr_matrix, np.ndarray]:
    """
    Penalite molle forcant field(node) ~= anchor_value sur les noeuds du
    plateau. Meme mecanique que les contraintes utilisateur (un seul
    mecanisme d'assemblage dans tout le fichier), avec un poids eleve
    dedie plutot qu'une elimination Dirichlet dure.
    """
    A = sp.csr_matrix((node_count, node_count))
    b = np.zeros(node_count)
    if len(platform_nodes) == 0:
        return A, b

    A = sp.csr_matrix(
        (np.full(len(platform_nodes), weight), (platform_nodes, platform_nodes)),
        shape=(node_count, node_count),
    )
    b[platform_nodes] = weight * anchor_value
    return A, b


# ---------------------------------------------------------------------------
# Decision 5 - Resolution du systeme lineaire
# ---------------------------------------------------------------------------

def _solve_linear_system(
    A: sp.csr_matrix,
    b: np.ndarray,
    params: FieldSolverParams,
    initial_guess: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Resout A.phi = b. Le systeme est symetrique defini positif par
    construction (rigidite FEM + penalites quadratiques), donc CG est le
    solveur iteratif approprie - pas GMRES, reserve aux matrices non
    symetriques.
    """
    node_count = A.shape[0]

    if node_count < params.direct_solve_node_threshold:
        try:
            return spla.spsolve(A.tocsc(), b)
        except Exception as exc:
            raise FieldSolverError(f"Le solveur direct (spsolve) a echoue: {exc}") from exc

    if params.use_ilu_preconditioner:
        try:
            ilu = spla.spilu(A.tocsc())
        except Exception as exc:
            raise FieldSolverError(f"La factorisation ILU a echoue: {exc}") from exc
        preconditioner = spla.LinearOperator(A.shape, ilu.solve)
    else:
        diagonal = A.diagonal().copy()
        diagonal[diagonal == 0.0] = 1.0
        preconditioner = spla.LinearOperator(A.shape, lambda x: x / diagonal)

    values, info = spla.cg(
        A,
        b,
        x0=initial_guess,
        rtol=params.convergence_tolerance,
        maxiter=params.max_solver_iterations,
        M=preconditioner,
    )

    if info > 0:
        raise FieldSolverError(
            f"Le solveur iteratif (cg) n'a pas converge apres "
            f"{params.max_solver_iterations} iterations (info={info}). Augmentez "
            f"max_solver_iterations ou convergence_tolerance, ou activez "
            f"use_ilu_preconditioner."
        )
    if info < 0:
        raise FieldSolverError(f"Le solveur iteratif (cg) a recu une entree illegale (info={info}).")

    return values


# ---------------------------------------------------------------------------
# Decision 6 - Systeme statique (invariant aux contraintes) / resolution
# ---------------------------------------------------------------------------

def build_static_system(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    region_dependency_graph: RegionDependencyGraph,
    params: FieldSolverParams,
) -> FieldSolverStaticSystem:
    """
    Assemble la partie du systeme qui NE depend PAS des Constraint[] :
    lissage (rigidite FEM), terme de surplomb, ancrage plateau, plus les
    structures partagees (KDTree des noeuds, gradients de fonctions de
    forme, volumes).

    A appeler UNE FOIS apres step3. Le resultat reste valide tant que la
    geometrie et les poids/seuils ne changent pas ; seule l'edition des
    Constraint[] utilisateur passe ensuite par `solve_with_constraints()`.
    """
    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)

    A = params.laplacian_smoothness_weight * _assemble_stiffness_matrix(
        tet_mesh, shape_gradients, volumes
    )
    b = np.zeros(tet_mesh.node_count)

    boundary_face_nodes, boundary_face_tets, nearest_surface_face_id = (
        _map_boundary_faces_to_surface(tet_mesh, mesh_data)
    )

    # Direction d'impression et risque par TETRAEDRE (pas seulement les
    # tetraedres frontiere) : chaque tetraedre herite de la face de surface
    # la plus proche de son centroide, via le kdtree deja construit a
    # l'etape 1 sur les centres de faces.
    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    _, nearest_face_per_tet = mesh_data.kdtree.query(tet_centers)
    nearest_face_per_tet = np.asarray(nearest_face_per_tet, dtype=np.int64)

    face_direction = np.zeros((mesh_data.face_count, 3))
    face_direction[:, 2] = 1.0  # defaut : verticale, si une face n'a pas de region
    for region in region_dependency_graph.regions.values():
        face_direction[region.face_ids] = region.target_direction

    tet_directions = face_direction[nearest_face_per_tet]

    # Risque effectif : une face posee sur le plateau n'est pas un surplomb,
    # meme si sa normale pointe vers le bas (limite de conception #1 de
    # CLAUDE.md). Meme raisonnement que step3._effective_risk_score().
    effective_risk = np.asarray(overhang_map.risk_score, dtype=np.float64).copy()
    if region_dependency_graph.platform_contact_faces is not None:
        effective_risk[np.asarray(region_dependency_graph.platform_contact_faces, dtype=bool)] = 0.0

    tet_risk = effective_risk[nearest_face_per_tet]

    # Sur-ponderation autour des surplombs : les tetraedres frontiere a
    # risque servent de germes, propages sur overhang_influence_depth sauts.
    risky_tets = np.zeros(tet_mesh.tet_count, dtype=bool)
    if len(boundary_face_tets) > 0:
        risky_boundary = effective_risk[nearest_surface_face_id] > 0.0
        seed_tets = np.unique(boundary_face_tets[risky_boundary])
        if len(seed_tets) > 0:
            risky_tets = _spread_risk_influence(
                _build_tet_tet_adjacency(tet_mesh),
                seed_tets,
                params.overhang_influence_depth,
            )

    # Tous les tetraedres recoivent le poids de base ; ceux marques a risque
    # recoivent en plus la sur-ponderation proportionnelle au risque.
    tet_weights = np.full(tet_mesh.tet_count, params.gradient_baseline_weight, dtype=np.float64)
    tet_weights[risky_tets] += params.gradient_overhang_weight * tet_risk[risky_tets]

    gradient_A, gradient_b = _assemble_gradient_direction_term(
        tet_mesh,
        shape_gradients,
        volumes,
        np.ones(tet_mesh.tet_count, dtype=bool),
        tet_directions,
        tet_weights,
    )
    A = A + gradient_A
    b = b + gradient_b

    if params.enable_platform_anchor:
        platform_nodes = _identify_platform_nodes(
            tet_mesh, mesh_data, region_dependency_graph
        )
        anchor_A, anchor_b = _assemble_platform_anchor(
            tet_mesh.node_count,
            platform_nodes,
            params.platform_anchor_weight,
            params.platform_anchor_value,
        )
        A = A + anchor_A
        b = b + anchor_b

    ctx = ConstraintAssemblyContext(
        tet_mesh=tet_mesh,
        mesh_data=mesh_data,
        overhang_map=overhang_map,
        region_dependency_graph=region_dependency_graph,
        node_kdtree=scipy.spatial.cKDTree(tet_mesh.tet_vertices),
        shape_gradients=shape_gradients,
        tet_volumes=volumes,
        params=params,
    )

    return FieldSolverStaticSystem(A_structural=A.tocsr(), b_structural=b, ctx=ctx)


def solve_with_constraints(
    static_system: FieldSolverStaticSystem,
    constraints: Optional[list[Constraint]],
    params: FieldSolverParams,
    initial_guess: Optional[np.ndarray] = None,
) -> ScalarField:
    """
    Ajoute les contraintes utilisateur au systeme statique et resout.
    C'est le chemin rapide apres l'edition d'un point d'attache : seules
    les contraintes sont reassemblees, la partie structurelle est reutilisee
    telle quelle.

    `initial_guess` (typiquement le champ precedent) accelere la convergence
    du solveur iteratif quand l'edition est une petite perturbation ; sans
    effet sur le chemin direct.
    """
    constraint_A, constraint_b = assemble_all_constraints(constraints, static_system.ctx)

    A = (static_system.A_structural + constraint_A).tocsr()
    b = static_system.b_structural + constraint_b

    values = _solve_linear_system(A, b, params, initial_guess=initial_guess)
    gradient = _compute_per_tet_gradient(
        static_system.ctx.tet_mesh, static_system.ctx.shape_gradients, values
    )

    return ScalarField(values=values, gradient=gradient)


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    region_dependency_graph: RegionDependencyGraph,
    constraints: Optional[list[Constraint]],
    params: FieldSolverParams,
) -> ScalarField:
    """
    Assemble et resout le champ scalaire continu sur le TetMesh, en une
    passe. Enveloppe de `build_static_system()` + `solve_with_constraints()`
    pour un slicing complet ; pour une re-resolution rapide apres edition
    d'un point d'attache, appeler ces deux fonctions separement et mettre
    le systeme statique en cache.

    Args:
        tet_mesh: sortie de pipeline/step1_ingest.py.
        mesh_data: sortie de pipeline/step1_ingest.py (pour le mapping
            surface -> volume via son kdtree).
        overhang_map: sortie de pipeline/step2_overhang.py.
        region_dependency_graph: sortie de pipeline/step3_segmentation.py.
        constraints: contraintes utilisateur, potentiellement vide/None.
        params: poids et parametres de resolution.

    Returns:
        ScalarField : une valeur par noeud du TetMesh, un gradient constant
        par tetraedre.

    Raises:
        FieldSolverError: si le solveur ne converge pas ou si le systeme
            assemble n'est pas exploitable.
    """
    static_system = build_static_system(
        tet_mesh, mesh_data, overhang_map, region_dependency_graph, params
    )
    return solve_with_constraints(static_system, constraints, params)
