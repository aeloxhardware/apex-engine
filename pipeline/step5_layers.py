"""
pipeline/step5_layers.py

Etape 5 du pipeline : extraction des couches courbes (isosurfaces) a partir
du ScalarField resolu par step4, sur le TetMesh de step1.

Depend de TetMesh + MeshData (step1), ScalarField (step4),
RegionDependencyGraph (step3).

Ce que ce fichier fait, dans l'ordre :
  1. Convertit `layer_height_base` (une distance physique) en un espacement
     de valeurs de champ, via la mediane ponderee par volume de ||grad(phi)||
     sur tout le TetMesh (une seule valeur representative - un champ non
     eikonal ne permet pas un espacement exactement constant partout, cf.
     Decision 2 du plan).
  2. Pour chaque isovaleur, extrait le patch triangule par marching
     tetrahedra (implementation maison - aucune bibliotheque du projet ne
     fournit cette operation sur un maillage tetraedrique arbitraire).
  3. En extrait les boucles de contour fermees via le trace de bord du patch
     (`trimesh.Trimesh.outline()`) : comme le TetMesh est conforme a la
     frontiere (tetgen appele avec le switch `p`), le bord du patch EST deja
     geometriquement l'intersection avec le mesh original - pas besoin d'une
     deuxieme operation geometrique explicite (cf. Decision 3 du plan).
  4. Associe chaque couche aux regions de step3 qu'elle traverse, et
     l'ordonne par ISOVALEUR SEULE - jamais par rang topologique de region
     (voir _compute_region_ranks et la docstring de run() plus bas pour la
     raison : le rang, discret et assigne par tetraedre, n'est PAS toujours
     coherent spatialement avec l'isovaleur, continue, pres d'une frontiere
     de region).

Sortie : liste ordonnee de Layer.
"""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np
import trimesh

from core.mesh_types import Layer, MeshData, RegionDependencyGraph, ScalarField, TetMesh


class LayerExtractionError(Exception):
    """Levee si le ScalarField ne peut pas etre converti en couches (gradient manquant, champ degenere)."""


@dataclass
class LayerParams:
    layer_height_base: float = 0.2
    # tolerance de fusion des sommets du patch triangule, avant extraction du
    # contour, relative a l'etendue du TetMesh - evite des micro-trous de
    # bord dus au bruit flottant plutot qu'a une vraie discontinuite
    # geometrique (cf. Decision 5 du plan).
    weld_tolerance_factor: float = 1e-6


# ---------------------------------------------------------------------------
# Decision 1 - Marching tetrahedra (3 cas geometriques)
# ---------------------------------------------------------------------------

def _interpolate_edge(
    p_from: np.ndarray, v_from: float, p_to: np.ndarray, v_to: float, isovalue: float
) -> np.ndarray:
    t = (isovalue - v_from) / (v_to - v_from)
    return p_from + t * (p_to - p_from)


def _orient_triangle(
    p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, reference_direction: np.ndarray
) -> np.ndarray:
    """
    Oriente le triangle (p0,p1,p2) pour que sa normale ait un produit
    scalaire positif avec `reference_direction`. Utilise pour toutes les
    variantes de marching tetrahedra (cf. Decision 3bis du plan) : garantit
    une orientation coherente sur tout le patch extrait, sans dependre de
    l'ordre local (arbitraire) des sommets du tetraedre source.
    """
    normal = np.cross(p1 - p0, p2 - p0)
    if np.dot(normal, reference_direction) < 0.0:
        return np.array([p0, p2, p1])
    return np.array([p0, p1, p2])


def _march_tetrahedron(
    tet_points: np.ndarray, tet_values: np.ndarray, isovalue: float
) -> list[np.ndarray]:
    """
    Extrait les triangles de l'isosurface `isovalue` pour UN tetraedre.

    tet_points: (4,3), tet_values: (4,). Retourne une liste de triangles
    (chacun (3,3)), toujours orientes normale pointant du cote "en dessous"
    vers le cote "au dessus" de l'isovaleur.

    La direction de reference pour l'orientation (`p[hi] - p[lo]`, ou hi/lo
    sont les sommets de valeur max/min du tetraedre) est GARANTIE d'avoir un
    produit scalaire positif avec le vrai gradient local : le champ est P1
    (lineaire par tetraedre), donc phi(hi) - phi(lo) = grad(phi) . (p[hi] -
    p[lo]) exactement, et phi(hi) > phi(lo) par construction. C'est ce qui
    rend l'orientation independante de l'ordre local (arbitraire, cf.
    step4_field_solver._compute_shape_gradients qui prend deja abs(volume)
    pour la meme raison) des sommets du tetraedre.

    Convention de degenerescence : un sommet dont la valeur egale exactement
    l'isovaleur est classe "au-dessus" (>=), pour eliminer toute ambiguite
    de classement sans biaiser les valeurs elles-memes.
    """
    above = tet_values >= isovalue
    n_above = int(above.sum())

    if n_above == 0 or n_above == 4:
        return []

    hi_idx = int(np.argmax(tet_values))
    lo_idx = int(np.argmin(tet_values))
    if tet_values[hi_idx] == tet_values[lo_idx]:
        # tetraedre a champ constant (degenere pour cette isovaleur) : pas
        # de direction de reference exploitable, on l'ignore proprement.
        return []
    reference_direction = tet_points[hi_idx] - tet_points[lo_idx]

    if n_above in (1, 3):
        lone_is_above = n_above == 1
        lone_idx = int(np.where(above == lone_is_above)[0][0])
        rest_idx = [i for i in range(4) if i != lone_idx]

        pts = []
        for r in rest_idx:
            if lone_is_above:
                p = _interpolate_edge(
                    tet_points[r], tet_values[r], tet_points[lone_idx], tet_values[lone_idx], isovalue
                )
            else:
                p = _interpolate_edge(
                    tet_points[lone_idx], tet_values[lone_idx], tet_points[r], tet_values[r], isovalue
                )
            pts.append(p)

        return [_orient_triangle(pts[0], pts[1], pts[2], reference_direction)]

    # n_above == 2 : quadrilatere -> 2 triangles. Les 4 aretes qui
    # traversent l'isovaleur relient chacun des 2 sommets "au-dessus" a
    # chacun des 2 sommets "en-dessous" (les aretes above-above et
    # below-below ne traversent pas).
    a1, a2 = [int(i) for i in np.where(above)[0]]
    b1, b2 = [int(i) for i in np.where(~above)[0]]

    p_a1b1 = _interpolate_edge(tet_points[b1], tet_values[b1], tet_points[a1], tet_values[a1], isovalue)
    p_a1b2 = _interpolate_edge(tet_points[b2], tet_values[b2], tet_points[a1], tet_values[a1], isovalue)
    p_a2b1 = _interpolate_edge(tet_points[b1], tet_values[b1], tet_points[a2], tet_values[a2], isovalue)
    p_a2b2 = _interpolate_edge(tet_points[b2], tet_values[b2], tet_points[a2], tet_values[a2], isovalue)

    return [
        _orient_triangle(p_a1b1, p_a1b2, p_a2b2, reference_direction),
        _orient_triangle(p_a1b1, p_a2b2, p_a2b1, reference_direction),
    ]


def _tet_volumes(tet_mesh: TetMesh) -> np.ndarray:
    """Volume de chaque tetraedre (meme formule que step4_field_solver, non partagee - trop petit pour justifier un utilitaire commun)."""
    verts = tet_mesh.tet_vertices[tet_mesh.tetrahedra]  # (T,4,3)
    edge_matrices = np.transpose(verts[:, 1:, :] - verts[:, 0:1, :], (0, 2, 1))  # (T,3,3)
    return np.abs(np.linalg.det(edge_matrices)) / 6.0


def _extract_isosurface_triangles(
    tet_mesh: TetMesh,
    values: np.ndarray,
    isovalue: float,
    tet_min: np.ndarray,
    tet_max: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extrait le patch triangule complet de l'isosurface `isovalue` sur tout
    le TetMesh. Filtre d'abord vectoriellement les tetraedres candidats
    (`tet_min <= isovalue <= tet_max`), pour ne classifier triangle par
    triangle que les tetraedres realistiquement traverses.

    Retourne (vertices (3N,3), faces (N,3), source_tet_ids (N,)) - un
    triplet de sommets et l'indice du tetraedre source par triangle extrait
    (necessaire pour l'association couche->regions de la decision 4).
    """
    candidates = np.where((tet_min <= isovalue) & (isovalue <= tet_max))[0]

    triangles: list[np.ndarray] = []
    source_tets: list[int] = []

    for t_idx in candidates:
        tet = tet_mesh.tetrahedra[t_idx]
        tris = _march_tetrahedron(tet_mesh.tet_vertices[tet], values[tet], isovalue)
        for tri in tris:
            triangles.append(tri)
            source_tets.append(int(t_idx))

    if not triangles:
        return np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64), np.zeros(0, dtype=np.int64)

    vertices = np.concatenate(triangles, axis=0)  # (3N,3)
    faces = np.arange(3 * len(triangles), dtype=np.int64).reshape(-1, 3)
    return vertices, faces, np.array(source_tets, dtype=np.int64)


def _weld_vertices(
    vertices: np.ndarray, faces: np.ndarray, tolerance: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fusionne les sommets proches (a `tolerance` pres) avant extraction du
    contour. Deux tetraedres adjacents interpolent independamment le meme
    point d'arete partagee ; ce devrait etre bit-identique en theorie (memes
    entrees, meme formule), mais fusionner par tolerance absorbe tout ecart
    residuel sans hypothese fragile sur l'arithmetique flottante.
    """
    if tolerance <= 0 or len(vertices) == 0:
        return vertices, faces

    rounded = np.round(vertices / tolerance).astype(np.int64)
    unique_rounded, unique_indices, inverse = np.unique(
        rounded, axis=0, return_index=True, return_inverse=True
    )
    welded_vertices = vertices[unique_indices]
    remapped_faces = inverse[faces].reshape(faces.shape)
    return welded_vertices, remapped_faces


# ---------------------------------------------------------------------------
# Decision 3 - Contours fermes via trace de bord (trimesh.outline())
# ---------------------------------------------------------------------------

def _extract_closed_contours(vertices: np.ndarray, faces: np.ndarray) -> list[np.ndarray]:
    """
    Extrait les boucles de contour fermees du patch triangule (vertices,
    faces), via le trace de bord de trimesh. `process=False` : on ne laisse
    pas trimesh reordonner/fusionner les sommets a notre place (deja fait
    explicitement par `_weld_vertices`), pour garder le controle complet de
    la correspondance triangle -> tetraedre source en amont.

    Gere nativement les patchs a plusieurs boucles disjointes (ex: paroi
    exterieure + paroi interieure d'une cavite a la meme hauteur de champ) -
    chaque boucle devient un element de la liste retournee.
    """
    if len(faces) == 0:
        return []

    try:
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
        outline = mesh.outline()
    except Exception:
        # patch degenere (triangles nuls, topologie non-variete locale) :
        # aucun contour exploitable plutot qu'un crash - cf. Decision 5.
        return []

    contours = []
    for entity in outline.entities:
        if not getattr(entity, "closed", True):
            continue
        points = np.asarray(entity.discrete(outline.vertices))
        if len(points) < 4:  # minimum : 3 sommets distincts + point de fermeture
            continue
        contours.append(points)

    return contours


# ---------------------------------------------------------------------------
# Decision 2 - Espacement des isovaleurs (mediane ponderee par volume)
# ---------------------------------------------------------------------------

def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    values_sorted = values[order]
    cum_weights = np.cumsum(weights[order])
    total = cum_weights[-1]
    if total <= 0.0:
        return float(np.median(values))
    idx = int(np.searchsorted(cum_weights, 0.5 * total))
    idx = min(idx, len(values_sorted) - 1)
    return float(values_sorted[idx])


def _representative_gradient_norm(scalar_field: ScalarField, tet_mesh: TetMesh) -> float:
    """
    Mediane ponderee par volume de ||grad(phi)|| sur tout le TetMesh - le
    scalaire "representatif" utilise pour convertir layer_height_base en
    espacement de champ (cf. Decision 2 du plan).

    Mediane plutot que moyenne : l'ancrage plateau de step4
    (`platform_anchor_weight`, poids eleve sur un sous-ensemble restreint de
    noeuds frontiere) tend a produire un gradient localement plus fort sur
    la fine couche de tetraedres adjacents a ces noeuds ancres. Une moyenne
    brute se laisse fausser par cette poignee de valeurs extremes ; la
    mediane ponderee reste representative du comportement typique du champ.
    """
    if scalar_field.gradient is None:
        raise LayerExtractionError(
            "ScalarField.gradient est requis (fourni par pipeline/step4_field_solver.py) "
            "pour convertir layer_height_base en espacement de champ."
        )

    norms = np.linalg.norm(scalar_field.gradient, axis=1)
    volumes = _tet_volumes(tet_mesh)
    return _weighted_median(norms, volumes)


def _generate_isovalues(scalar_field: ScalarField, tet_mesh: TetMesh, params: LayerParams) -> np.ndarray:
    """
    Genere la sequence d'isovaleurs espacees de `delta = layer_height_base *
    representatif(||grad(phi)||)`, entre les bornes REELLES du champ resolu
    (min/max mesures, pas supposees etre exactement [0, max] - l'ancrage
    plateau de step4 est une penalite molle, pas une elimination dure).
    """
    gradient_norm = _representative_gradient_norm(scalar_field, tet_mesh)
    if gradient_norm <= 1e-12:
        raise LayerExtractionError(
            "La norme representative du gradient du champ est quasi nulle - le "
            "champ scalaire semble degenere (constant), impossible d'en deriver "
            "un espacement d'isovaleurs."
        )

    delta = params.layer_height_base * gradient_norm
    v_min = float(scalar_field.values.min())
    v_max = float(scalar_field.values.max())

    if v_max - v_min < delta:
        # modele plus fin qu'une seule couche : une couche unique au centre
        return np.array([0.5 * (v_min + v_max)])

    first = v_min + 0.5 * delta
    count = int(np.floor((v_max - first) / delta)) + 1
    return first + delta * np.arange(count)


# ---------------------------------------------------------------------------
# Decision 4 - Association couche -> regions, rang topologique
# ---------------------------------------------------------------------------

def _map_tets_to_regions(
    tet_mesh: TetMesh, mesh_data: MeshData, region_dependency_graph: RegionDependencyGraph
) -> np.ndarray:
    """
    Pour chaque tetraedre, retrouve la region de sa face de surface la plus
    proche (via mesh_data.kdtree sur le centroide du tetraedre) - meme
    technique que step4_field_solver.build_static_system pour le meme
    pont surface/volume, pas reinventee ici.

    Retourne (T,) region_id, ou -1 si aucune region ne couvre la face
    trouvee (ne devrait pas arriver : step3 garantit que chaque face a
    exactement une region - garde-fou defensif, pas un cas attendu).
    """
    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    _, nearest_face = mesh_data.kdtree.query(tet_centers)
    nearest_face = np.asarray(nearest_face, dtype=np.int64)

    face_to_region = np.full(mesh_data.face_count, -1, dtype=np.int64)
    for region in region_dependency_graph.regions.values():
        face_to_region[region.face_ids] = region.id

    return face_to_region[nearest_face]


def _compute_region_ranks(region_dependency_graph: RegionDependencyGraph) -> dict[int, int]:
    """
    Rang topologique de chaque region : plus long chemin depuis une region
    racine (aucune arete entrante) dans le DAG de dependances. Une region
    racine a rang 0. Accumulation en ordre topologique (meme outil -
    networkx - deja utilise par step3_segmentation pour la detection de
    cycle, nouvel usage ici).

    NE SERT PLUS A TRIER LES COUCHES dans run() (voir sa docstring) : un
    tetraedre proche d'une frontiere de region peut avoir une isovaleur
    proche de celle de la region voisine sans que son rang topologique le
    reflete, ce qui produisait un ordre d'impression aberrant (mesure
    concretement sur un L a 2 regions - une couche de rang 0 isolee tres
    haute en isovaleur se retrouvait triee avant tout le corps de la region
    de rang 1, qui la precede pourtant physiquement). Conservee : fonction
    pure, testee independamment, utile a un futur besoin qui raisonnerait
    reellement sur l'ordre des regions (pas les couches elles-memes).
    """
    graph = nx.DiGraph()
    graph.add_nodes_from(region_dependency_graph.regions.keys())
    graph.add_edges_from(region_dependency_graph.dependencies)

    ranks: dict[int, int] = {}
    for region_id in nx.topological_sort(graph):
        predecessors = list(graph.predecessors(region_id))
        ranks[region_id] = 0 if not predecessors else 1 + max(ranks[p] for p in predecessors)
    return ranks


# ---------------------------------------------------------------------------
# Decision 6 - Repere local (origine, axe_u, axe_v)
# ---------------------------------------------------------------------------

def _compute_local_frame(
    contours: list[np.ndarray], tet_gradients: np.ndarray, tet_volumes: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Repere local approximatif de la couche : origine = centroide des points
    de contour, normale = gradient moyen (pondere par volume) des
    tetraedres contribuant a la couche, axe_u/axe_v = base orthonormee du
    plan perpendiculaire (meme construction que
    core.geometry_utils.rotation_axis_circle_directions, reutilisee ici).
    """
    all_points = np.concatenate(contours, axis=0)
    origin = all_points.mean(axis=0)

    weighted_gradient = (tet_gradients * tet_volumes[:, None]).sum(axis=0)
    norm = np.linalg.norm(weighted_gradient)
    normal = np.array([0.0, 0.0, 1.0]) if norm < 1e-12 else weighted_gradient / norm

    reference = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(reference, normal)) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])

    axis_u = reference - np.dot(reference, normal) * normal
    axis_u = axis_u / np.linalg.norm(axis_u)
    axis_v = np.cross(normal, axis_u)

    return origin, axis_u, axis_v


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    scalar_field: ScalarField,
    region_dependency_graph: RegionDependencyGraph,
    params: LayerParams,
) -> list[Layer]:
    """
    Extrait les couches courbes du ScalarField et les ordonne.

    Args:
        tet_mesh: sortie de pipeline/step1_ingest.py.
        mesh_data: sortie de pipeline/step1_ingest.py (kdtree pour le pont
            surface/volume).
        scalar_field: sortie de pipeline/step4_field_solver.py (gradient
            requis, pas seulement les valeurs).
        region_dependency_graph: sortie de pipeline/step3_segmentation.py.
        params: hauteur de couche cible et tolerances.

    Returns:
        Liste ordonnee de Layer (id reassigne selon l'ordre d'impression
        final : ISOVALEUR SEULE, jamais le rang topologique des regions
        traversees. Le champ scalaire est resolu globalement (FEM,
        lissage laplacien) et son isovaleur encode deja par construction
        la progression spatiale coherente de l'ancrage plateau vers
        l'exterieur - un tetraedre proche d'une frontiere de region peut
        legitimement avoir une isovaleur proche de celle de la region
        voisine (continuite du lissage), sans que cela remette en cause
        l'ordre physique reel. Trier par rang d'abord produisait un ordre
        d'impression aberrant : une couche de rang 0 isolee tres haute en
        isovaleur (donc physiquement situee pres du sommet du bras) se
        retrouvait triee AVANT tout le corps de la region de rang 1 (le
        bras lui-meme), qui la precede pourtant physiquement - mesure
        concretement sur res/test_models/l_shape_overhang.stl une fois
        step3 corrige pour scinder base/bras en 2 regions reelles (voir
        CLAUDE.md). `pipeline/step7_scheduling.py` consomme cette liste
        via `enumerate(layers)` comme ordre d'impression reel - ce n'est
        pas qu'un souci d'affichage, un mauvais ordre ici casse
        l'ordonnancement en aval.

    Raises:
        LayerExtractionError: si le champ n'a pas de gradient ou est
            degenere (norme de gradient quasi nulle partout).
    """
    isovalues = _generate_isovalues(scalar_field, tet_mesh, params)

    tet_values_per_tet = scalar_field.values[tet_mesh.tetrahedra]  # (T,4)
    tet_min = tet_values_per_tet.min(axis=1)
    tet_max = tet_values_per_tet.max(axis=1)

    tet_regions = _map_tets_to_regions(tet_mesh, mesh_data, region_dependency_graph)
    volumes = _tet_volumes(tet_mesh)

    scale = float(np.ptp(tet_mesh.tet_vertices, axis=0).max()) or 1.0
    weld_tolerance = params.weld_tolerance_factor * scale

    raw_layers: list[tuple[float, list[np.ndarray], list[int], tuple[np.ndarray, np.ndarray, np.ndarray]]] = []

    for isovalue in isovalues:
        vertices, faces, source_tets = _extract_isosurface_triangles(
            tet_mesh, scalar_field.values, float(isovalue), tet_min, tet_max
        )
        if len(faces) == 0:
            continue

        welded_vertices, welded_faces = _weld_vertices(vertices, faces, weld_tolerance)
        contours = _extract_closed_contours(welded_vertices, welded_faces)
        if not contours:
            continue

        unique_tets = np.unique(source_tets)
        region_ids = sorted({int(r) for r in tet_regions[unique_tets] if r >= 0})

        local_frame = _compute_local_frame(
            contours, scalar_field.gradient[unique_tets], volumes[unique_tets]
        )

        raw_layers.append((float(isovalue), contours, region_ids, local_frame))

    # Isovaleur SEULE - voir la docstring de run() pour pourquoi le rang
    # topologique de region ne doit plus intervenir ici.
    raw_layers.sort(key=lambda entry: entry[0])

    return [
        Layer(
            id=i,
            contour_3d=contours,
            region_ids=region_ids,
            local_frame=local_frame,
            isovalue=isovalue,
        )
        for i, (isovalue, contours, region_ids, local_frame) in enumerate(raw_layers)
    ]
