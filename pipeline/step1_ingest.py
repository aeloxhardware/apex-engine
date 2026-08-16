"""
pipeline/step1_ingest.py

Etape 1 du pipeline : ingestion du mesh source (STL/3MF) et preparation de
toutes les structures dont les etapes suivantes ont besoin.

Cette etape est la seule du pipeline qui touche a l'I/O disque cote geometrie.
Elle ne depend d'aucune autre etape - c'est le point d'entree du pipeline et
le fichier le plus simple a tester en isolation (aucune donnee amont requise,
juste un chemin de fichier).

Sortie : (MeshData, TetMesh, IngestReport)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.spatial
import trimesh

from core.mesh_types import MeshData, TetMesh


@dataclass
class IngestReport:
    """
    Resume de l'ingestion, remonte a l'UI TypeScript pour donner un retour
    sur la qualite du mesh d'entree avant de lancer le reste du pipeline.
    """
    was_watertight_initially: bool
    is_watertight_after_repair: bool
    faces_removed_degenerate: int
    vertices_merged: int
    tetrahedra_generated: int
    source_face_count: int
    source_vertex_count: int


class MeshIngestError(Exception):
    """Levee quand le mesh source ne peut pas etre charge ou exploite."""


# ---------------------------------------------------------------------------
# Etape 1 - Chargement du fichier brut
# ---------------------------------------------------------------------------

def _load_raw_mesh(mesh_path: str, allow_multi_mesh: bool) -> trimesh.Trimesh:
    loaded = trimesh.load(mesh_path, force="mesh" if not allow_multi_mesh else None)

    if isinstance(loaded, trimesh.Scene):
        if not allow_multi_mesh:
            raise MeshIngestError(
                f"'{mesh_path}' contient plusieurs objets (Scene) mais "
                "allow_multi_mesh=False. Fusionnez le fichier ou activez ce parametre."
            )
        geometries = list(loaded.geometry.values())
        if len(geometries) == 0:
            raise MeshIngestError(f"'{mesh_path}' ne contient aucune geometrie.")
        loaded = trimesh.util.concatenate(geometries)

    if not isinstance(loaded, trimesh.Trimesh):
        raise MeshIngestError(
            f"'{mesh_path}' n'a pas produit un Trimesh exploitable "
            f"(type obtenu: {type(loaded)})."
        )

    if loaded.faces.shape[0] == 0:
        raise MeshIngestError(f"'{mesh_path}' est charge mais ne contient aucune face.")

    return loaded


# ---------------------------------------------------------------------------
# Etape 2 - Nettoyage et reparation
# ---------------------------------------------------------------------------

def _clean_and_repair(mesh: trimesh.Trimesh) -> tuple[trimesh.Trimesh, IngestReport]:
    was_watertight_initially = bool(mesh.is_watertight)

    vertex_count_before = mesh.vertices.shape[0]
    face_count_before = mesh.faces.shape[0]

    mesh.merge_vertices()
    # API trimesh >= 5.0 : ces methodes retournent un masque booleen a
    # appliquer via update_faces, plutot que de modifier le mesh en place.
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.update_faces(mesh.unique_faces())

    if not mesh.is_watertight:
        trimesh.repair.fill_holes(mesh)
        trimesh.repair.fix_normals(mesh)
        trimesh.repair.fix_winding(mesh)

    vertices_merged = max(0, vertex_count_before - mesh.vertices.shape[0])
    faces_removed = max(0, face_count_before - mesh.faces.shape[0])

    report = IngestReport(
        was_watertight_initially=was_watertight_initially,
        is_watertight_after_repair=bool(mesh.is_watertight),
        faces_removed_degenerate=faces_removed,
        vertices_merged=vertices_merged,
        tetrahedra_generated=0,  # rempli plus tard, apres tetraedrisation
        source_face_count=mesh.faces.shape[0],
        source_vertex_count=mesh.vertices.shape[0],
    )
    return mesh, report


# ---------------------------------------------------------------------------
# Etape 4 - Structures d'acceleration (BVH + KDTree)
# ---------------------------------------------------------------------------

def _build_acceleration_structures(mesh: trimesh.Trimesh):
    # RayMeshIntersector choisit automatiquement pyembree si disponible dans
    # l'environnement, sinon retombe sur l'implementation pure numpy - pas
    # besoin de brancher explicitement, trimesh gere le fallback lui-meme.
    bvh = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    kdtree = scipy.spatial.cKDTree(mesh.triangles_center)
    return bvh, kdtree


# ---------------------------------------------------------------------------
# Etape 5 - Assemblage du MeshData
# ---------------------------------------------------------------------------

def _build_mesh_data(mesh: trimesh.Trimesh) -> MeshData:
    bvh, kdtree = _build_acceleration_structures(mesh)
    return MeshData(
        vertices=np.asarray(mesh.vertices, dtype=np.float64),
        faces=np.asarray(mesh.faces, dtype=np.int64),
        face_normals=np.asarray(mesh.face_normals, dtype=np.float64),
        vertex_normals=np.asarray(mesh.vertex_normals, dtype=np.float64),
        is_watertight=bool(mesh.is_watertight),
        bvh=bvh,
        kdtree=kdtree,
    )


# ---------------------------------------------------------------------------
# Etape 6 - Tetraedrisation volumique (tetgen)
# ---------------------------------------------------------------------------

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
    chaque coque-cavite sert de point representatif, verifie via
    `shell.contains()` avant utilisation.
    """
    shells = mesh.split(only_watertight=True)
    if len(shells) <= 1:
        return []

    hole_points = []
    for shell in shells:
        if shell.volume >= 0:
            continue  # coque solide (frontiere exterieure ou paroi normale), pas une cavite

        centroid = shell.vertices.mean(axis=0)
        if not shell.contains([centroid])[0]:
            raise MeshIngestError(
                "Cavite interne detectee (volume signe negatif) mais son "
                "centroide n'est pas a l'interieur de sa propre coque "
                "(cavite non-convexe) - impossible de deriver un point de "
                "trou fiable pour tetgen sans heuristique plus poussee."
            )
        hole_points.append(centroid)

    return hole_points


def _tetrahedralize(mesh: trimesh.Trimesh, resolution: float) -> TetMesh:
    """
    Genere un maillage tetraedrique volumique respectant la frontiere du mesh
    de surface, via TetGen. `resolution` est converti en volume maximal
    autorise par tetraedre (parametre -a de tetgen).
    """
    import tetgen

    if not mesh.is_watertight:
        raise MeshIngestError(
            "Impossible de tetraedriser un mesh non etanche apres reparation. "
            "Verifiez la geometrie source."
        )

    max_tet_volume = resolution ** 3

    tet = tetgen.TetGen(mesh.vertices, mesh.faces)
    for hole_point in _find_internal_cavities(mesh):
        tet.add_hole(hole_point)
    tet.tetrahedralize(switches=f"pq1.2a{max_tet_volume:.10f}")
    # Lire depuis les attributs de l'objet plutot que la valeur de retour :
    # le nombre de valeurs retournees par tetrahedralize() varie selon la
    # version de tetgen installee, tandis que .node/.elem sont stables.
    tet_vertices = tet.node
    tetrahedra = tet.elem

    tet_adjacency = _build_tet_adjacency(tetrahedra, len(tet_vertices))

    return TetMesh(
        tet_vertices=np.asarray(tet_vertices, dtype=np.float64),
        tetrahedra=np.asarray(tetrahedra, dtype=np.int64),
        tet_adjacency=tet_adjacency,
    )


def _build_tet_adjacency(tetrahedra: np.ndarray, node_count: int) -> sp.csr_matrix:
    """
    Construit la matrice d'adjacence creuse noeud-noeud du maillage
    tetraedrique (1 si deux noeuds partagent une arete d'un tetraedre).
    Reutilisee par l'assemblage du Laplacien discret a l'etape 4 du pipeline.
    """
    edges = set()
    for tet in tetrahedra:
        for i in range(4):
            for j in range(i + 1, 4):
                a, b = int(tet[i]), int(tet[j])
                edges.add((min(a, b), max(a, b)))

    if len(edges) == 0:
        return sp.csr_matrix((node_count, node_count))

    rows, cols = zip(*edges)
    rows = np.array(rows)
    cols = np.array(cols)

    data = np.ones(len(rows) * 2)
    all_rows = np.concatenate([rows, cols])
    all_cols = np.concatenate([cols, rows])

    adjacency = sp.csr_matrix(
        (data, (all_rows, all_cols)), shape=(node_count, node_count)
    )
    return adjacency


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    mesh_path: str,
    resolution: float,
    allow_multi_mesh: bool = False,
) -> tuple[MeshData, TetMesh, IngestReport]:
    """
    Charge, nettoie et prepare le mesh source pour le reste du pipeline.

    Args:
        mesh_path: chemin vers le fichier STL/3MF source.
        resolution: taille cible de tetraedre (FieldSolverParams.resolution),
            utilisee pour controler la densite de la tetraedrisation.
        allow_multi_mesh: si True, fusionne les geometries multiples d'une
            Scene plutot que de lever une erreur.

    Returns:
        Tuple (MeshData, TetMesh, IngestReport).

    Raises:
        MeshIngestError: si le fichier est illisible, vide, ou ne peut pas
            etre reduit a un volume tetraedrisable meme apres reparation.
    """
    raw_mesh = _load_raw_mesh(mesh_path, allow_multi_mesh)
    cleaned_mesh, report = _clean_and_repair(raw_mesh)

    mesh_data = _build_mesh_data(cleaned_mesh)
    tet_mesh = _tetrahedralize(cleaned_mesh, resolution)

    report.tetrahedra_generated = tet_mesh.tet_count

    return mesh_data, tet_mesh, report