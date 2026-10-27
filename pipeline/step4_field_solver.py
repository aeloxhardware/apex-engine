"""
pipeline/step4_field_solver.py

Etape 4 du pipeline : le champ scalaire continu sur le TetMesh, dont les
isosurfaces (etape 5) definissent les couches courbes.

Le champ EST la distance au plateau a travers la matiere, calculee par la
methode de la chaleur (Crane, Weischedel, Wardetzky) : chaque couche est une
coquille equidistante du plateau, donc geodesiquement homogene, et l'ordre
d'impression par isovaleur croissante est la progression de la matiere. Une
piece suspendue vient apres ce qui la porte parce qu'elle en est plus loin.

Trois etapes, sur les matrices FEM P1 (rigidite et masse concentree) :
diffusion de la chaleur depuis le plateau, direction normalisee de son flux
(eventuellement ramenee au demi-plan atteignable par la tete), puis
integration de Poisson `L*phi = div(X)`.

Depend de TetMesh + MeshData (step1) et de RegionDependencyGraph (step3, qui
seul connait le contact avec le plateau). Sortie : ScalarField.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from core import _cpp
from core import progress as _progress
from core.mesh_types import (
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
    # Taille cible des tetraedres (mm), lue par step1.
    resolution: float = 1.5
    # Solveur lineaire : direct (`spsolve`) sous `direct_solve_node_threshold`
    # noeuds, sinon gradient conjugue preconditionne (Jacobi, ou ILU).
    max_solver_iterations: int = 2000
    convergence_tolerance: float = 1e-8
    direct_solve_node_threshold: int = 30000
    use_ilu_preconditioner: bool = False
    # Penalite molle de poids eleve qui ancre le champ a 0 sur les noeuds
    # source (plateau) : leve la singularite du Laplacien pur.
    platform_anchor_weight: float = 1e6
    up_vector: Optional[np.ndarray] = None
    # Temps de diffusion = heat_method_time_factor * (arete mediane)^2.
    # 100 retenu : les valeurs sont convergees de 10 a 1000 (55.30 contre
    # 55.36, soit 0.1 % d'ecart) ; a 1.0 le champ sature loin du plateau et la
    # direction de son gradient n'est plus que du bruit numerique.
    heat_method_time_factor: float = 100.0
    # Origine du plateau : renseignee, la direction du flux de la chaleur est
    # projetee sur le demi-plan (+radial, vertical) atteignable par la tete
    # (C est impose par la position, B ne penche que vers +R) AVANT
    # l'integration. `main._effective_field_solver_params` la renseigne depuis
    # le profil machine quand `position_coupled_axes == [0]`.
    radial_tilt_origin: Optional[np.ndarray] = None
    heat_project_reachable: bool = True

    def __post_init__(self):
        if self.up_vector is None:
            self.up_vector = np.array([0.0, 0.0, 1.0])


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

    if _cpp.enabled():
        # Formule des cofacteurs en C++ plutot que `np.linalg.inv` (LAPACK) :
        # deux algorithmes differents pour le meme inverse 3x3, donc des
        # arrondis differents. Ecart RELATIF mesure sur death.stl : 7,2e-14,
        # soit la precision machine - voir la note d'equivalence de step4,
        # ou la tolerance remplace l'egalite bit a bit (limite #56, pas 7).
        gradients, volumes = _cpp.module().compute_shape_gradients(
            np.ascontiguousarray(tet_mesh.tet_vertices, dtype=np.float64),
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64),
        )
        return gradients, volumes

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


class _FemPattern:
    """
    Motif de sparsite CSR de l'assemblage FEM, partageable entre matrices.

    Les deux matrices de step4 (rigidite et terme de direction) ont
    EXACTEMENT le meme motif : celui de l'adjacence noeud-noeud du maillage,
    qui ne depend que de la topologie, jamais des poids. Le chemin scipy le
    reconstruisait pourtant deux fois, en triant a chaque fois les 6,2 M
    entrees COO (16 par tetraedre) puis en sommant les doublons.

    Construit explicitement par l'appelant et passe aux deux assemblages -
    jamais memorise sur le TetMesh : un cache invisible survivrait a une
    modification de `tetrahedra` et produirait une matrice silencieusement
    fausse.
    """

    __slots__ = ("_native", "indptr", "indices", "node_count")

    def __init__(self, tet_mesh: TetMesh):
        self.node_count = tet_mesh.node_count
        self._native = _cpp.module().FemSparsity(
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64), tet_mesh.node_count
        )
        self.indptr = self._native.indptr()
        self.indices = self._native.indices()

    def gram(
        self,
        tet_mesh: TetMesh,
        shape_gradients: np.ndarray,
        scale: np.ndarray,
        active: Optional[np.ndarray],
    ) -> sp.csr_matrix:
        """Matrice `sum_e scale_e * G_e^T G_e`, assemblee DIRECTEMENT en CSR."""
        data = self._native.assemble_gram(
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64),
            np.ascontiguousarray(shape_gradients, dtype=np.float64),
            np.ascontiguousarray(scale, dtype=np.float64),
            None if active is None else np.ascontiguousarray(active, dtype=np.int32),
        )
        return sp.csr_matrix(
            (data, self.indices, self.indptr), shape=(self.node_count, self.node_count)
        )


def _assemble_stiffness_matrix(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    tet_weights: Optional[np.ndarray] = None,
    pattern: Optional["_FemPattern"] = None,
) -> sp.csr_matrix:
    """
    Assemble la matrice de rigidite FEM P1 globale (Laplacien discret
    geometriquement correct) : K_e = w_e * V_e * G_e^T G_e par element, sommee
    aux indices de sommets globaux. Symetrique semi-definie positive, de
    noyau les champs constants.

    `tet_weights` (T,) module le poids de lissage PAR ELEMENT (voir
    `FieldSolverParams.smoothness_modulation_strength`). `None` = poids
    uniforme, strictement identique au comportement historique - et pas
    seulement a une tolerance pres : c'est la MEME expression, sans
    multiplication supplementaire.

    Le MOTIF DE SPARSITE ne depend pas des poids, seulement de la topologie
    du maillage. Une factorisation directe (`spsolve`) coute donc exactement
    la meme chose avec ou sans modulation ; seul un solveur iteratif peut
    voir son conditionnement change.
    """
    scale = volumes if tet_weights is None else volumes * tet_weights

    if _cpp.enabled():
        # Assemblage DIRECT en CSR : ni matrices locales materialisees, ni
        # COO intermediaire de 6,2 M entrees a trier puis dedupliquer.
        pattern = pattern if pattern is not None else _FemPattern(tet_mesh)
        return pattern.gram(tet_mesh, shape_gradients, scale, None)

    local_matrices = scale[:, None, None] * np.einsum(
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


def _assemble_lumped_mass_matrix(
    tet_mesh: TetMesh, volumes: np.ndarray
) -> sp.csr_matrix:
    """
    Matrice de MASSE P1 condensee (lumped) : chaque tetraedre distribue son
    volume a parts egales sur ses 4 sommets, et la matrice est diagonale.

    Condensee et non consistante : la methode de la chaleur ne demande qu'un
    terme de masse bien conditionne, et une diagonale garantit que
    `M + t*L` reste definie positive sans hypothese sur la qualite du
    maillage - point qui compte ici, `tetgen` ne garantissant aucune
    regularite (meme raison que le `abs(volume)` deja utilise ailleurs dans
    ce fichier).
    """
    node_count = tet_mesh.node_count
    contributions = np.repeat(np.abs(volumes) / 4.0, 4)
    indices = tet_mesh.tetrahedra.ravel()
    diagonal = np.zeros(node_count)
    np.add.at(diagonal, indices, contributions)
    return sp.diags(diagonal, format="csr")


def _assemble_divergence_rhs(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    tet_field: np.ndarray,
) -> np.ndarray:
    """
    Second membre de l'equation de Poisson `L*phi = div(X)` pour un champ
    vectoriel `X` constant par tetraedre.

    En P1, `integral(grad(w_i) . X)` sur un element vaut `V_e * (G_e^T X_e)_i` :
    c'est exactement la meme contraction que
    `_assemble_gradient_direction_term` utilise pour son second membre, avec
    le volume comme poids au lieu du poids de surplomb. Aucune nouvelle
    formule d'integration n'est introduite.
    """
    local = volumes[:, None] * np.einsum("tki,tk->ti", shape_gradients, tet_field)
    rhs = np.zeros(tet_mesh.node_count)
    np.add.at(rhs, tet_mesh.tetrahedra.ravel(), local.ravel())
    return rhs


def _heat_method_distance(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    source_nodes: np.ndarray,
    params: FieldSolverParams,
) -> np.ndarray:
    """
    Distance geodesique LISSE au plateau, par la methode de la chaleur
    (Crane, Weischedel, Wardetzky).

    Pourquoi elle plutot que Dijkstra sur les aretes : Dijkstra suit les
    aretes du maillage en zigzag, ce qui donne une distance dont la norme du
    gradient vaut 4,81 en mediane sur death.stl au lieu de 1,0 attendu pour
    une distance. Cette irregularite la rend inutilisable comme champ
    generateur de couches. La methode de la chaleur produit une distance
    CONTINUE dont `||grad|| ~ 1` par construction, sur le meme domaine et
    avec les memes matrices FEM deja assemblees ici.

    Trois etapes :
      1. diffusion `(M + t*L) u = source` sur un temps court `t` ;
      2. direction du flux `X = -grad(u)/||grad(u)||` par tetraedre - seule
         la DIRECTION est retenue, c'est ce qui elimine l'irregularite
         d'amplitude ;
      3. integration `L*phi = div(X)`, puis recalage a 0 sur la source.

    `t` est derive de la taille caracteristique du maillage (carre de la
    longueur d'arete mediane), jamais une constante absolue : la formulation
    de la chaleur est dimensionnee en longueur au carre, une valeur fixe
    n'aurait pas de sens d'un modele a l'autre.

    Retourne un tableau par NOEUD. Les deux systemes sont singuliers du seul
    champ constant (Laplacien pur), ce que la source pour le premier et le
    recalage final pour le second suffisent a lever.
    """
    node_count = tet_mesh.node_count
    if len(source_nodes) == 0:
        return np.full(node_count, np.inf)

    if _cpp.enabled() and not params.use_ilu_preconditioner and hasattr(_cpp.module(), "heat_method_distance"):
        # Chemin natif (tolerance, pas le bit : LDLT/CG d'Eigen contre
        # SuperLU/CG de scipy). None = non convergence : on retombe sur Python,
        # qui leve l'erreur habituelle.
        project = params.heat_project_reachable and params.radial_tilt_origin is not None
        up = np.asarray(params.up_vector if params.up_vector is not None else [0.0, 0.0, 1.0], dtype=float)
        native = _cpp.module().heat_method_distance(
            np.ascontiguousarray(tet_mesh.tet_vertices, dtype=np.float64),
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64),
            np.ascontiguousarray(shape_gradients, dtype=np.float64),
            np.ascontiguousarray(volumes, dtype=np.float64),
            np.ascontiguousarray(source_nodes, dtype=np.int64),
            float(params.heat_method_time_factor),
            float(params.platform_anchor_weight),
            float(params.convergence_tolerance),
            int(params.max_solver_iterations),
            int(params.direct_solve_node_threshold),
            np.asarray(params.radial_tilt_origin, dtype=float) if project else None,
            up,
        )
        if native is not None:
            return native

    L = _assemble_stiffness_matrix(tet_mesh, shape_gradients, volumes)
    M = _assemble_lumped_mass_matrix(tet_mesh, volumes)

    # Temps de diffusion, homogene a une longueur au carre.
    edges = tet_mesh.tet_vertices[tet_mesh.tetrahedra[:, [0, 1]]]
    h = float(np.median(np.linalg.norm(edges[:, 0] - edges[:, 1], axis=1))) or 1.0
    t = params.heat_method_time_factor * h * h

    source = np.zeros(node_count)
    source[source_nodes] = 1.0
    def solve(A, b):
        # Le systeme de diffusion (temps t = 100 h^2) est mal conditionne sur un
        # grand modele (fig.stl, 165 mm, 309 000 noeuds : 2 000 iterations de CG
        # ne suffisent pas). Repli UNIQUE avec 10x plus d'iterations, jamais en
        # premiere intention : le chemin habituel reste inchange.
        try:
            return _solve_linear_system(A, b, params)
        except FieldSolverError:
            more = dataclasses.replace(params, max_solver_iterations=max(20000, 10 * params.max_solver_iterations))
            return _solve_linear_system(A, b, more)

    u = solve((M + t * L).tocsr(), source)

    gradients = _compute_per_tet_gradient(tet_mesh, shape_gradients, u)
    norms = np.linalg.norm(gradients, axis=1)
    flow = np.zeros_like(gradients)
    usable = norms > 1e-12
    flow[usable] = -gradients[usable] / norms[usable][:, None]

    if params.heat_project_reachable and params.radial_tilt_origin is not None:
        projected, _changed = _project_to_reachable_radial_tilt(
            tet_mesh, flow, params.radial_tilt_origin,
            params.up_vector if params.up_vector is not None else np.array([0.0, 0.0, 1.0]),
        )
        flow[usable] = projected[usable]

    rhs = _assemble_divergence_rhs(tet_mesh, shape_gradients, volumes, flow)

    # Le Laplacien pur est singulier (noyau = constantes). On ancre les
    # noeuds source a 0 par penalite forte plutot que d'eliminer des lignes :
    # meme mecanisme que l'ancrage plateau deja en place, et le systeme reste
    # symetrique defini positif.
    anchor = sp.diags(
        np.where(np.isin(np.arange(node_count), source_nodes),
                 params.platform_anchor_weight, 0.0),
        format="csr",
    )
    phi = solve((L + anchor).tocsr(), rhs)

    # Une distance est positive et nulle a la source ; le signe de la
    # solution de Poisson n'est pas garanti par le systeme lui-meme.
    if phi[source_nodes].mean() > phi.mean():
        phi = -phi
    return phi - phi[source_nodes].min()


def _compute_platform_level(
    mesh_data: MeshData, region_dependency_graph: RegionDependencyGraph
) -> Optional[float]:
    """
    Niveau geometrique (z) du plateau, deduit de
    `RegionDependencyGraph.platform_contact_faces` (step3, seule etape qui
    connait la position reelle du plateau) - sommet le plus haut des faces
    en contact, la bande de tolerance de step3 est donc deja incluse.
    Retourne None si aucune face de contact n'est disponible (segmentation
    sans plateau identifie) - source de verite UNIQUE, reutilisee a la fois
    par `_identify_platform_nodes` (ancrage de phi) et
    `_blend_tet_directions` (source d'influence "plateau" du melange de
    direction) pour ne jamais faire diverger silencieusement deux calculs
    du meme niveau.
    """
    platform_faces = region_dependency_graph.platform_contact_faces
    if platform_faces is None:
        return None
    platform_faces = np.asarray(platform_faces, dtype=bool)
    if not np.any(platform_faces):
        return None
    platform_vertex_z = mesh_data.vertices[mesh_data.faces[platform_faces]][:, :, 2]
    return float(platform_vertex_z.max())


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

    platform_level = _compute_platform_level(mesh_data, region_dependency_graph)
    if platform_level is None:
        return np.zeros(0, dtype=np.int64)

    scale = float(np.ptp(tet_mesh.tet_vertices[:, 2])) or 1.0
    epsilon = 1e-9 * scale

    return np.where(tet_mesh.tet_vertices[:, 2] <= platform_level + epsilon)[0]


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


def _project_to_reachable_radial_tilt(
    tet_mesh: TetMesh, directions: np.ndarray, origin: np.ndarray, up_vector: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Projette chaque direction sur le demi-plan (+radial, vertical) de son
    tetraedre : seul degre d'orientation de la machine actuelle quand C est
    impose par la position (limite #59). Rend (directions, masque des
    tetraedres modifies de plus de 1e-6). Pres de l'axe (rayon < 1e-6), la
    radiale est indefinie et C y est libre : direction laissee intacte."""
    up = np.asarray(up_vector, dtype=float)
    up = up / np.linalg.norm(up)
    centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1) - np.asarray(origin, dtype=float)
    radial = centers - np.outer(centers @ up, up)
    rn = np.linalg.norm(radial, axis=1)
    ok = rn > 1e-6
    rhat = np.zeros_like(radial)
    rhat[ok] = radial[ok] / rn[ok, None]
    d_r = np.einsum("ij,ij->i", directions, rhat)
    d_z = directions @ up
    d_r = np.maximum(d_r, 0.0)
    projected = d_r[:, None] * rhat + d_z[:, None] * up
    norm = np.linalg.norm(projected, axis=1)
    degenerate = norm < 1e-9          # direction purement tangentielle : verticale
    projected[degenerate] = up
    norm[degenerate] = 1.0
    projected /= norm[:, None]
    result = np.where(ok[:, None], projected, directions)
    changed = np.linalg.norm(result - directions, axis=1) > 1e-6
    return result, changed


def run(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    region_dependency_graph: RegionDependencyGraph,
    params: FieldSolverParams,
) -> ScalarField:
    """
    Calcule le champ scalaire : la distance au plateau a travers la matiere,
    par la methode de la chaleur. `overhang_map` n'intervient pas (l'orientation
    des couches vient de la geometrie du champ) ; il reste dans la signature
    pour garder le contrat entre etapes.

    Raises:
        FieldSolverError: si le contact avec le plateau est inconnu ou si le
            solveur ne converge pas.
    """
    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)
    platform_nodes = _identify_platform_nodes(tet_mesh, mesh_data, region_dependency_graph)
    _progress.milestone(0.1)                 # R-017 : gradients de forme et nœuds du plateau
    values = _heat_method_distance(tet_mesh, shape_gradients, volumes, platform_nodes, params)
    _progress.milestone(0.9)                 # diffusion, flux et integration resolus
    return ScalarField(
        values=values,
        gradient=_compute_per_tet_gradient(tet_mesh, shape_gradients, values),
    )
