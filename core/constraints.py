"""
core/constraints.py

Logique de fusion des contraintes utilisateur dans le systeme lineaire de
l'etape 4 (resolution du champ scalaire).

Ce fichier est le POINT D'EXTENSION du solveur de champ : ajouter un nouveau
type de contrainte (renfort structurel, nouvelle regle esthetique, etc.) se
fait ici, en ajoutant une fonction d'assemblage + une entree dans
CONSTRAINT_ASSEMBLERS. `pipeline/step4_field_solver.py` n'a jamais besoin
d'etre modifie : il appelle uniquement `assemble_all_constraints()` et ne
connait aucun ConstraintType par son nom.

Sens de dependance : pipeline/step4_field_solver.py importe ce fichier,
jamais l'inverse (meme convention que core/mesh_types.py et
core/geometry_utils.py).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Optional

import numpy as np
import scipy.sparse as sp

from core.mesh_types import (
    Constraint,
    ConstraintFalloff,
    ConstraintType,
    MeshData,
    OverhangMap,
    RegionDependencyGraph,
    TetMesh,
)


class ConstraintAssemblyError(Exception):
    """Levee si une contrainte ne peut pas etre assemblee (type non enregistre, falloff inconnu...)."""


# ---------------------------------------------------------------------------
# Contexte partage passe a chaque assembleur
# ---------------------------------------------------------------------------

@dataclass
class ConstraintAssemblyContext:
    """
    Donnees partagees dont un assembleur de contrainte peut avoir besoin.
    Aucune logique de calcul ici - uniquement un regroupement, pour eviter
    que chaque assembleur recalcule les memes structures (KDTree sur les
    noeuds, gradients de fonctions de forme...).

    `shape_gradients` et `tet_volumes` sont precalcules une seule fois par
    `build_static_system()` (etape 4) : ils ne dependent que de la geometrie
    du TetMesh, jamais des contraintes elles-memes.
    """
    tet_mesh: TetMesh
    mesh_data: MeshData
    overhang_map: OverhangMap
    region_dependency_graph: RegionDependencyGraph
    node_kdtree: Any                     # core.spatial.KDTree sur tet_vertices (scipy ou natif)
    shape_gradients: np.ndarray          # shape (T, 3, 4) - gradients des fonctions de forme P1
    tet_volumes: np.ndarray              # shape (T,)
    params: Any = None                   # FieldSolverParams, si un assembleur en a besoin


# ---------------------------------------------------------------------------
# Noyaux d'attenuation
# ---------------------------------------------------------------------------

def falloff_weight(
    distance: np.ndarray, influence_radius: float, falloff: ConstraintFalloff
) -> np.ndarray:
    """
    Noyau d'attenuation en fonction de la distance au centre de la
    contrainte. Vaut 1 au centre et decroit jusqu'a 0 au bord de
    `influence_radius`.
    """
    normalized = np.clip(distance / max(influence_radius, 1e-12), 0.0, 1.0)

    if falloff == ConstraintFalloff.LINEAR:
        return 1.0 - normalized
    if falloff == ConstraintFalloff.GAUSSIAN:
        # exp(-4x^2) : ~0.018 au bord, decroissance douce et centree
        return np.exp(-4.0 * normalized ** 2)
    if falloff == ConstraintFalloff.CUBIC:
        return (1.0 - normalized ** 3) ** 3

    raise ConstraintAssemblyError(f"Falloff non reconnu: {falloff}")


# ---------------------------------------------------------------------------
# Assembleurs concrets
# ---------------------------------------------------------------------------

def _assemble_field_value_penalty(
    constraint: Constraint, ctx: ConstraintAssemblyContext
) -> tuple[sp.csr_matrix, np.ndarray]:
    """
    Penalite molle sur la VALEUR du champ : force field(node) ~=
    target_field_value pour les noeuds situes dans influence_radius de la
    position de la contrainte, ponderes par le noyau de falloff.

    Contribue a la fois a A (terme diagonal) et a b (second membre) -
    c'est ce qui permet a un point d'attache de reellement tirer le champ
    vers une valeur donnee.
    """
    node_count = ctx.tet_mesh.node_count
    A = sp.csr_matrix((node_count, node_count))
    b = np.zeros(node_count)

    if constraint.target_field_value is None:
        return A, b

    nearby = ctx.node_kdtree.query_ball_point(constraint.position, constraint.influence_radius)
    if not nearby:
        return A, b
    nearby = np.asarray(nearby, dtype=np.int64)

    distances = np.linalg.norm(ctx.tet_mesh.tet_vertices[nearby] - constraint.position, axis=1)
    weights = constraint.weight * falloff_weight(
        distances, constraint.influence_radius, constraint.falloff
    )

    positive = weights > 0.0
    if not np.any(positive):
        return A, b
    nearby = nearby[positive]
    weights = weights[positive]

    A = sp.csr_matrix((weights, (nearby, nearby)), shape=(node_count, node_count))
    np.add.at(b, nearby, weights * constraint.target_field_value)

    return A, b


def _assemble_gradient_direction_penalty(
    constraint: Constraint, ctx: ConstraintAssemblyContext
) -> tuple[sp.csr_matrix, np.ndarray]:
    """
    Penalite molle sur la DIRECTION du gradient : energie de correspondance
    E_e = w_e * V_e * ||grad(phi) - d||^2 sur les tetraedres proches de la
    contrainte, ou d = target_normal normalisee.

    Sa minimisation donne, par element :
        A_local = w_e * V_e * G_e^T G_e
        b_local = w_e * V_e * G_e^T d

    Le terme en b est essentiel : une penalite purement homogene (par
    exemple sur la seule composante perpendiculaire, G_e^T (I - d d^T) G_e)
    n'injecterait aucune force motrice et laisserait le systeme s'effondrer
    sur le champ nul en l'absence d'autre contrainte de valeur.
    """
    node_count = ctx.tet_mesh.node_count
    A = sp.csr_matrix((node_count, node_count))
    b = np.zeros(node_count)

    if constraint.target_normal is None:
        return A, b

    direction_norm = np.linalg.norm(constraint.target_normal)
    if direction_norm < 1e-12:
        raise ConstraintAssemblyError(
            f"target_normal de norme nulle sur une contrainte {constraint.type} "
            f"en position {constraint.position} - direction indefinie."
        )
    d = np.asarray(constraint.target_normal, dtype=np.float64) / direction_norm

    # tetraedres concernes : ceux dont le centroide est dans le rayon
    # d'influence (test au centre plutot que par sommet - un tetraedre est
    # l'unite sur laquelle le gradient est defini, pas ses sommets pris
    # individuellement).
    tet_centers = ctx.tet_mesh.tet_vertices[ctx.tet_mesh.tetrahedra].mean(axis=1)
    distances = np.linalg.norm(tet_centers - constraint.position, axis=1)
    affected = np.where(distances <= constraint.influence_radius)[0]
    if len(affected) == 0:
        return A, b

    weights = constraint.weight * falloff_weight(
        distances[affected], constraint.influence_radius, constraint.falloff
    )

    rows, cols, data = [], [], []
    for local_idx, t_idx in enumerate(affected):
        weight = weights[local_idx]
        volume = ctx.tet_volumes[t_idx]
        if weight <= 0.0 or volume <= 0.0:
            continue

        gradients = ctx.shape_gradients[t_idx]          # (3,4)
        local_matrix = weight * volume * (gradients.T @ gradients)   # (4,4)
        local_rhs = weight * volume * (gradients.T @ d)              # (4,)

        tet = ctx.tet_mesh.tetrahedra[t_idx]
        for i in range(4):
            b[int(tet[i])] += local_rhs[i]
            for j in range(4):
                rows.append(int(tet[i]))
                cols.append(int(tet[j]))
                data.append(local_matrix[i, j])

    if rows:
        A = sp.csr_matrix((data, (rows, cols)), shape=(node_count, node_count))

    return A, b


def _assemble_user_constraint(
    constraint: Constraint, ctx: ConstraintAssemblyContext
) -> tuple[sp.csr_matrix, np.ndarray]:
    """
    Assembleur generique des contraintes utilisateur : applique la penalite
    de valeur si `target_field_value` est fourni, ET/OU la penalite de
    direction si `target_normal` est fourni. Les deux sont cumulables sur
    une meme contrainte (l'UI peut envoyer les deux champs).
    """
    value_A, value_b = _assemble_field_value_penalty(constraint, ctx)
    normal_A, normal_b = _assemble_gradient_direction_penalty(constraint, ctx)
    return value_A + normal_A, value_b + normal_b


# ---------------------------------------------------------------------------
# Registre : le point d'extension
# ---------------------------------------------------------------------------

ConstraintAssembler = Callable[
    [Constraint, ConstraintAssemblyContext], "tuple[sp.csr_matrix, np.ndarray]"
]

# Pour ajouter un nouveau type de contrainte : ecrire sa fonction
# d'assemblage ci-dessus, puis l'enregistrer ici. Aucun autre fichier du
# pipeline n'a besoin d'etre modifie.
CONSTRAINT_ASSEMBLERS: dict[ConstraintType, ConstraintAssembler] = {
    ConstraintType.ATTACH_POINT: _assemble_user_constraint,
    ConstraintType.AESTHETIC: _assemble_user_constraint,
    ConstraintType.STRUCTURAL: _assemble_user_constraint,
}


def assemble_all_constraints(
    constraints: Optional[list[Constraint]], ctx: ConstraintAssemblyContext
) -> tuple[sp.csr_matrix, np.ndarray]:
    """
    Somme les contributions (A, b) de toutes les contraintes utilisateur,
    en dispatchant chacune vers son assembleur enregistre.

    Args:
        constraints: contraintes envoyees par l'UI, potentiellement vide/None.
        ctx: contexte partage precalcule par l'etape 4.

    Returns:
        Tuple (A_contraintes, b_contraintes), a additionner au systeme
        structurel de l'etape 4.

    Raises:
        ConstraintAssemblyError: si un type de contrainte n'a pas
            d'assembleur enregistre dans CONSTRAINT_ASSEMBLERS.
    """
    node_count = ctx.tet_mesh.node_count
    A = sp.csr_matrix((node_count, node_count))
    b = np.zeros(node_count)

    if not constraints:
        return A, b

    for constraint in constraints:
        assembler = CONSTRAINT_ASSEMBLERS.get(constraint.type)
        if assembler is None:
            raise ConstraintAssemblyError(
                f"Aucun assembleur enregistre pour le type de contrainte "
                f"{constraint.type}. Ajoutez-le a CONSTRAINT_ASSEMBLERS dans "
                f"core/constraints.py."
            )
        constraint_A, constraint_b = assembler(constraint, ctx)
        A = A + constraint_A
        b = b + constraint_b

    return A, b
