"""
pipeline/step2_overhang.py

Etape 2 du pipeline : detection des surplombs et des ponts internes.

Ne depend que de MeshData (sortie de step1_ingest.py) - testable en isolation
avec des mesh synthetiques simples.

Deux detections distinctes et complementaires :
  - surplomb (overhang) : une face dont la normale pointe trop vers le bas
    par rapport a la verticale de reference.
  - pont interne (bridge) : une face quasi-horizontale (donc non detectee
    comme surplomb) mais qui n'a rien en dessous d'elle dans le mesh -
    typiquement le "plafond" d'une cavite interne.

Sortie : OverhangMap
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from core.mesh_types import MeshData, OverhangMap


@dataclass
class OverhangParams:
    overhang_angle_threshold_deg: float = 45.0
    bridge_detection_ray_density: int = 4          # subdivisions barycentriques par face candidate
    bridge_support_ratio_threshold: float = 0.7    # proportion de points non supportes -> is_bridge
    bridge_max_gap_distance: float = 2.0            # distance max de ray-cast avant "non supporte"
    # NOTE: ce defaut (2.0, generalement en mm) est arbitraire et depend
    # fortement de l'echelle du modele et de la hauteur de couche cible.
    # A ajuster explicitement par l'appelant plutot que de compter sur ce
    # defaut - il doit typiquement etre de l'ordre de quelques hauteurs de
    # couche (layer_height_base), pas une valeur absolue fixe.
    bridge_horizontal_angle_tolerance_deg: float = 10.0  # tolerance autour de l'horizontale pour etre "candidat pont"
    up_vector: np.ndarray = None                    # defaut [0,0,1], voir __post_init__

    def __post_init__(self):
        if self.up_vector is None:
            self.up_vector = np.array([0.0, 0.0, 1.0])


# ---------------------------------------------------------------------------
# Etape 1 - Angle de surplomb par face
# ---------------------------------------------------------------------------

def _compute_overhang_angles(mesh_data: MeshData, up_vector: np.ndarray) -> np.ndarray:
    """
    Retourne, pour chaque face, l'angle (en degres) entre sa normale et la
    verticale descendante (-up_vector). Plus cet angle est petit, plus la
    face pointe vers le bas (surplomb severe). 0 deg = pointe droit vers le
    bas, 180 deg = pointe droit vers le haut.
    """
    down_vector = -up_vector / np.linalg.norm(up_vector)
    normals = mesh_data.face_normals
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    norms[norms == 0] = 1.0  # evite division par zero sur des normales degenerees
    normalized = normals / norms

    dot = np.clip(normalized @ down_vector, -1.0, 1.0)
    angle_from_down_deg = np.degrees(np.arccos(dot))
    return angle_from_down_deg


def _compute_is_overhang(angle_from_down_deg: np.ndarray, threshold_deg: float) -> np.ndarray:
    """
    Une face est en surplomb si elle pointe a moins de `threshold_deg` de la
    verticale descendante (angle_from_down petit = pointe vers le bas).
    """
    return angle_from_down_deg < threshold_deg


# ---------------------------------------------------------------------------
# Etape 2 - Score de risque continu
# ---------------------------------------------------------------------------

def _compute_risk_score(angle_from_down_deg: np.ndarray, threshold_deg: float) -> np.ndarray:
    """
    Score continu dans [0, 1] via une sigmoide centree sur le seuil, plutot
    qu'une coupure dure - evite les regions instables au niveau du seuil
    lors de la segmentation (etape 3 du pipeline).

    risk = 1  -> surplomb severe (pointe droit vers le bas)
    risk = 0  -> aucun risque (pointe droit vers le haut)
    risk = 0.5 -> exactement au seuil
    """
    # Largeur de transition de la sigmoide : ~1/4 du seuil, purement empirique,
    # suffisant pour eviter une coupure brutale sans lisser excessivement.
    transition_width = max(threshold_deg * 0.25, 1e-3)
    x = (threshold_deg - angle_from_down_deg) / transition_width
    risk = 1.0 / (1.0 + np.exp(-x))
    return risk


# ---------------------------------------------------------------------------
# Etape 3 - Detection des ponts internes par ray-casting
# ---------------------------------------------------------------------------

def _barycentric_sample_points(vertices: np.ndarray, faces: np.ndarray, subdivisions: int) -> list[np.ndarray]:
    """
    Pour chaque face (triangle), genere une grille de points par coordonnees
    barycentriques a `subdivisions` niveaux. Retourne une liste (une entree
    par face) de tableaux (K, 3) de points echantillons.
    """
    tri_verts = vertices[faces]  # shape (M, 3, 3)
    n = max(subdivisions, 1)

    bary_coords = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            bary_coords.append((i / n, j / n, k / n))
    bary_coords = np.array(bary_coords)  # shape (K, 3)

    # points[f, p, :] = bary_coords[p] @ tri_verts[f]
    points_per_face = np.einsum("pk,fkc->fpc", bary_coords, tri_verts)
    return [points_per_face[f] for f in range(points_per_face.shape[0])]


def _is_face_horizontal_candidate(
    angle_from_down_deg: np.ndarray, tolerance_deg: float
) -> np.ndarray:
    """
    Une face est candidate 'pont' si elle est quasi-horizontale ET que sa
    normale pointe vers le BAS (angle proche de 0 deg par rapport a la
    verticale descendante) - c'est le cas d'un plafond de cavite, qui a
    potentiellement besoin d'un support par en dessous.

    Une face quasi-horizontale dont la normale pointe vers le HAUT (angle
    proche de 180 deg, comme le dessus d'un solide plein) n'est PAS
    candidate : elle regarde vers le ciel, elle n'a par definition pas
    besoin d'etre supportee par en dessous - ce qu'il y a en dessous d'elle
    est necessairement du solide (ou le plateau), jamais du vide, sinon la
    face elle-meme n'existerait pas a cet endroit du mesh.
    """
    return angle_from_down_deg < tolerance_deg


def _detect_bridges(
    mesh_data: MeshData,
    angle_from_down_deg: np.ndarray,
    params: OverhangParams,
) -> np.ndarray:
    """
    Pour chaque face candidate (quasi-horizontale, normale vers le bas),
    echantillonne des points et lance un ray-cast vertical vers le bas
    depuis chacun. Si une proportion suffisante des points n'a rien en
    dessous, la face est marquee is_bridge.
    """
    # NOTE: on ne filtre PAS sur (~is_overhang) ici. Un plafond de cavite
    # quasi-horizontal (angle proche de 0) est aussi, par un seuil d'angle
    # simple, un "surplomb" classique - mais c'est precisement le cas que le
    # mecanisme de detection de pont (ray-cast, vide en dessous) doit
    # traiter differemment d'un surplomb en pente franche. Les deux
    # detections ne sont donc pas mutuellement exclusives par construction :
    # une face peut etre a la fois is_overhang (angle < seuil) et testee
    # comme candidate pont.
    candidates = _is_face_horizontal_candidate(
        angle_from_down_deg, params.bridge_horizontal_angle_tolerance_deg
    )

    is_bridge = np.zeros(mesh_data.face_count, dtype=bool)
    candidate_indices = np.where(candidates)[0]

    if len(candidate_indices) == 0:
        return is_bridge

    sample_points_per_face = _barycentric_sample_points(
        mesh_data.vertices, mesh_data.faces[candidate_indices], params.bridge_detection_ray_density
    )

    down_direction = np.array([0.0, 0.0, -1.0])
    # petit decalage pour eviter que le rayon ne s'auto-intersecte avec sa
    # propre face d'origine au point de depart
    epsilon = 1e-6

    for local_idx, face_idx in enumerate(candidate_indices):
        points = sample_points_per_face[local_idx]
        origins = points + down_direction * epsilon
        directions = np.tile(down_direction, (len(points), 1))

        locations, index_ray, _ = mesh_data.bvh.intersects_location(origins, directions)

        unsupported_count = 0
        for p_idx in range(len(points)):
            hits_for_point = locations[index_ray == p_idx]
            if len(hits_for_point) == 0:
                unsupported_count += 1
                continue
            distances = np.linalg.norm(hits_for_point - origins[p_idx], axis=1)
            if distances.min() > params.bridge_max_gap_distance:
                unsupported_count += 1

        unsupported_ratio = unsupported_count / len(points)
        if unsupported_ratio >= params.bridge_support_ratio_threshold:
            is_bridge[face_idx] = True

    return is_bridge


# ---------------------------------------------------------------------------
# Etape 4 - Direction de portee du pont (PCA locale par composante connexe)
# ---------------------------------------------------------------------------

def _build_face_adjacency(faces: np.ndarray, face_count: int) -> sp.csr_matrix:
    """Graphe d'adjacence face-face : deux faces sont adjacentes si elles partagent une arete."""
    edge_to_faces: dict[tuple[int, int], list[int]] = {}
    for f_idx, face in enumerate(faces):
        for i in range(3):
            a, b = int(face[i]), int(face[(i + 1) % 3])
            edge = (min(a, b), max(a, b))
            edge_to_faces.setdefault(edge, []).append(f_idx)

    rows, cols = [], []
    for shared_faces in edge_to_faces.values():
        if len(shared_faces) == 2:
            f1, f2 = shared_faces
            rows.extend([f1, f2])
            cols.extend([f2, f1])

    if len(rows) == 0:
        return sp.csr_matrix((face_count, face_count))

    data = np.ones(len(rows))
    return sp.csr_matrix((data, (rows, cols)), shape=(face_count, face_count))


def _connected_components(mask: np.ndarray, adjacency: sp.csr_matrix) -> list[np.ndarray]:
    """
    Regroupe les indices de faces marquees `mask` en composantes connexes,
    en ne traversant que les aretes reliant deux faces elles-memes marquees.
    """
    n_components, labels = sp.csgraph.connected_components(
        adjacency, directed=False
    )
    components = []
    for label in range(n_components):
        face_indices = np.where((labels == label) & mask)[0]
        if len(face_indices) > 0:
            components.append(face_indices)
    return components


def _compute_bridge_span_directions(
    mesh_data: MeshData, is_bridge: np.ndarray, face_adjacency: sp.csr_matrix
) -> np.ndarray:
    """
    Pour chaque composante connexe de faces-pont, calcule par PCA l'axe de
    plus grande extension dans le plan approximatif de la zone. Cette
    direction sera utilisee (etape 6 du pipeline) pour orienter la
    trajectoire perpendiculairement a la portee du pont.
    """
    span_direction = np.full((mesh_data.face_count, 3), np.nan)

    if not np.any(is_bridge):
        return span_direction

    components = _connected_components(is_bridge, face_adjacency)

    for face_indices in components:
        vert_indices = np.unique(mesh_data.faces[face_indices].ravel())
        points = mesh_data.vertices[vert_indices]

        centroid = points.mean(axis=0)
        centered = points - centroid

        # PCA via SVD : le premier vecteur singulier droit donne l'axe de
        # plus grande variance, i.e. la direction de plus grande extension.
        try:
            _, _, vt = np.linalg.svd(centered, full_matrices=False)
            main_axis = vt[0]
        except np.linalg.LinAlgError:
            main_axis = np.array([1.0, 0.0, 0.0])

        span_direction[face_indices] = main_axis

    return span_direction


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(mesh_data: MeshData, params: OverhangParams) -> OverhangMap:
    """
    Detecte les surplombs et les ponts internes du mesh, et calcule un score
    de risque continu par face utilise comme poids de contrainte a l'etape 3
    (segmentation) et l'etape 4 (champ scalaire).

    Args:
        mesh_data: sortie de pipeline/step1_ingest.py.
        params: seuils et parametres de detection.

    Returns:
        OverhangMap avec risk_score, is_overhang, is_bridge et
        bridge_span_direction alignes sur l'indexation des faces.
    """
    angle_from_down_deg = _compute_overhang_angles(mesh_data, params.up_vector)

    is_overhang = _compute_is_overhang(angle_from_down_deg, params.overhang_angle_threshold_deg)
    risk_score = _compute_risk_score(angle_from_down_deg, params.overhang_angle_threshold_deg)

    is_bridge = _detect_bridges(mesh_data, angle_from_down_deg, params)

    face_adjacency = _build_face_adjacency(mesh_data.faces, mesh_data.face_count)
    bridge_span_direction = _compute_bridge_span_directions(mesh_data, is_bridge, face_adjacency)

    return OverhangMap(
        risk_score=risk_score,
        is_overhang=is_overhang,
        is_bridge=is_bridge,
        bridge_span_direction=bridge_span_direction,
    )
