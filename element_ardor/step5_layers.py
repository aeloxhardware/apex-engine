"""
element_ardor/step5_layers.py

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
     l'ordonne par rang topologique maximal (la contrainte la plus stricte
     parmi les regions traversees, cf. Decision 4) puis par isovaleur - ou,
     si `split_disjoint_clusters` est actif, par un ordre DFS de continuite
     spatiale a l'interieur de chaque rang (cf. Decision 8) : necessaire des
     qu'une meme region step3 (une nappe de surface connexe) couvre
     plusieurs composantes spatialement disjointes du solide (ex: pied et
     pointe de bras d'un T partageant la meme region "verticale" - sinon le
     tri par isovaleur seule les fait alterner a chaque couche, CLAUDE.md).

Sortie : liste ordonnee de Layer.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import networkx as nx
import numpy as np
import trimesh

from core.geometry_utils import (
    point_in_polygon_2d,
)
from core.mesh_types import (
    AABB,
    Layer,
    MeshData,
    RegionDependencyGraph,
    ScalarField,
    TetMesh,
)
from core import _cpp, runtime
from core import progress as _progress
from core.spatial import KDTree

from element_ardor.step5_params import (  # noqa: F401  (réexport : imports existants)
    LayerExtractionError,
    LayerParams,
)
from element_ardor.step5_marching import (  # noqa: F401  (réexport : imports existants)
    _CAP_EXCLUDE_TABLE,
    _QUAD_BIT_WEIGHTS,
    _QUAD_TABLE_A,
    _QUAD_TABLE_B,
    _build_quad_tables,
    _extract_isosurface_triangles,
    _interpolate_edge,
    _interpolate_edge_batch,
    _march_tetrahedra_batch,
    _march_tetrahedron,
    _orient_triangle,
    _orient_triangle_batch,
    _tet_volumes,
    _weld_vertices,
)
from element_ardor.step5_ordering import (  # noqa: F401  (réexport : imports existants)
    _build_continuity_graph,
    _cluster_contours,
    _dfs_order_within_rank,
    _loops_connected,
    _min_pairwise_distance,
    _spatial_dfs_order,
    _subtree_sizes,
    _tets_for_cluster,
)


# ---------------------------------------------------------------------------
# Decision 1bis - Version VECTORISEE de _march_tetrahedron, appliquee a TOUS
# les tetraedres candidats d'une isovaleur a la fois (zero boucle Python sur
# les tetraedres/triangles). Remplace le goulot mesure : _march_tetrahedron
# etait appelee ~321000 fois sur death.stl (89% du temps de step5), dont
# 421545 appels a np.cross sur des vecteurs 3x3 individuels via
# _orient_triangle - meme anti-pattern de dispatch numpy generique deja
# corrige pour l'IK (CLAUDE.md limite #37). Logique mathematique
# STRICTEMENT identique a _march_tetrahedron/_interpolate_edge/
# _orient_triangle (conservees, jamais supprimees - _march_tetrahedron sert
# d'oracle de test), seule l'application passe du scalaire au batch.
# ---------------------------------------------------------------------------


# Nombre de couches ou le trace natif a rendu la main a trimesh (bord non-variete) :
# `runtime.state().outline_fallbacks`, verifie par les tests d'equivalence.


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
    # Triangles DEGENERES (deux sommets confondus par `_weld_vertices`) : quand
    # un noeud est a moins de la tolerance de soudure du plan, les points de
    # coupe de toutes ses aretes fusionnent. Sans aire, ils ne portent aucun
    # bord, mais leurs aretes repetees (multiplicite 6 a 22 mesuree) cassent le
    # trace de bord : boucles ouvertes, donc jetees. Mesure sur un Benchy
    # (noeuds a 2,6e-5 et 9,2e-5 mm du plan) : couche entiere perdue (192 mm2)
    # et 23 % d'une autre.
    if len(faces):
        faces = faces[(faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2])
                      & (faces[:, 0] != faces[:, 2])]
    if len(faces) == 0:
        return []

    # Pas 14 du portage C++ : trace de bord natif, identique a outline() dans le
    # cas des boucles simples ; None = bord non-variete -> trimesh ci-dessous.
    if _cpp.enabled() and hasattr(_cpp.module(), "boundary_loops"):
        loops = _cpp.module().boundary_loops(np.ascontiguousarray(faces, dtype=np.int64), len(vertices))
        if loops is not None:
            return [vertices[loop] for loop in loops if len(loop) >= 4]
        runtime.state().outline_fallbacks += 1

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
            "ScalarField.gradient est requis (fourni par element_ardor/step4_field_solver.py) "
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
# Decision 7 - Clustering spatial des boucles disjointes non imbriquees
# (CLAUDE.md : jonction en T/etoile, une isovaleur peut traverser
# simultanement deux parties du solide tres eloignees l'une de l'autre)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Decision 8 - Continuite spatiale inter-isovaleur (ordre DFS des clusters)
#
# Le clustering (Decision 7) resout le cas ou une couche a une isovaleur
# DONNEE a plusieurs boucles disjointes. Mais deux clusters issus d'une
# meme region topologique step3 (une region est une nappe de surface
# connexe, PAS necessairement une seule composante spatiale du solide - ex:
# le pied et la pointe d'un bras d'un modele en T partagent souvent la
# region "verticale" unique, cf. CLAUDE.md) ont le meme rang ET les memes
# region_ids : le tri _sort_key=(rank, isovalue) les fait alterner a chaque
# isovaleur (pied->pointe->pied->pointe...), produisant un trajet a vide a
# chaque bascule. Ce mecanisme construit un graphe de continuite spatiale
# entre clusters d'isovaleurs successives, puis un parcours DFS qui termine
# une branche avant de sauter a une autre - generalise a un nombre
# arbitraire de branches, aux recombinaisons (deux branches qui
# fusionnent), sans jamais violer le rang topologique de region (partition
# stricte avant le DFS, pas une ponderation).
# ---------------------------------------------------------------------------


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
        tet_mesh: sortie de element_ardor/step1_ingest.py.
        mesh_data: sortie de element_ardor/step1_ingest.py (kdtree pour le pont
            surface/volume).
        scalar_field: sortie de element_ardor/step4_field_solver.py (gradient
            requis, pas seulement les valeurs).
        region_dependency_graph: sortie de element_ardor/step3_segmentation.py.
        params: hauteur de couche cible et tolerances.

    Returns:
        Liste ordonnee de Layer (id reassigne selon l'ordre d'impression
        final : rang topologique des regions traversees, puis isovaleur).

    Raises:
        LayerExtractionError: si le champ n'a pas de gradient ou est
            degenere (norme de gradient quasi nulle partout).
    """
    isovalues = _generate_isovalues(scalar_field, tet_mesh, params)

    tet_values_per_tet = scalar_field.values[tet_mesh.tetrahedra]  # (T,4)
    tet_regions = _map_tets_to_regions(tet_mesh, mesh_data, region_dependency_graph)
    region_ranks = _compute_region_ranks(region_dependency_graph)
    volumes = _tet_volumes(tet_mesh)
    scale = float(np.ptp(tet_mesh.tet_vertices, axis=0).max()) or 1.0
    # Seuil de connexion spatiale entre boucles (clustering et tri par continuite) : quelques hauteurs de
    # couche separent le prolongement direct d'un cluster de ce qui ne l'est pas (mesure sur le modele en T).
    gap_threshold = (params.split_disjoint_clusters if params.split_disjoint_clusters is not None else 3.0) \
        * params.layer_height_base
    ctx = _ExtractContext(
        tet_mesh=tet_mesh, scalar_field=scalar_field, params=params, tet_regions=tet_regions, volumes=volumes,
        tet_min=tet_values_per_tet.min(axis=1), tet_max=tet_values_per_tet.max(axis=1),
        weld_tolerance=params.weld_tolerance_factor * scale, gap_threshold=gap_threshold,
        cluster_tet_tolerance=_cluster_tet_tolerance(tet_mesh, params, scale),
        gradient_c=np.ascontiguousarray(scalar_field.gradient, dtype=float),
        volumes_c=np.ascontiguousarray(volumes, dtype=float),
        tet_regions_c=np.ascontiguousarray(tet_regions, dtype=np.int64),
    )
    # Pas 13 du portage C++ (limite #56) : marching tetrahedra + soudure en natif, construit UNE fois par champ.
    if _cpp.enabled() and hasattr(_cpp.module(), "IsosurfaceExtractor"):
        ctx.native_extractor = _cpp.module().IsosurfaceExtractor(
            np.ascontiguousarray(tet_mesh.tet_vertices, dtype=float),
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64),
            np.ascontiguousarray(scalar_field.values, dtype=float),
        )

    raw_layers: list[tuple[float, list[np.ndarray], list[int], tuple[np.ndarray, np.ndarray, np.ndarray]]] = []
    for done, isovalue in enumerate(isovalues):
        _progress.report(done, len(isovalues))   # R-017 : couches deja extraites
        raw_layers.extend(_extract_layer_entries(float(isovalue), ctx))
    raw_layers = _order_raw_layers(raw_layers, region_ranks, gap_threshold, params)

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


@dataclass
class _ExtractContext:
    """Donnees communes a l'extraction de toutes les couches d'un champ."""
    tet_mesh: TetMesh
    scalar_field: ScalarField
    params: LayerParams
    tet_regions: np.ndarray
    volumes: np.ndarray
    tet_min: np.ndarray
    tet_max: np.ndarray
    weld_tolerance: float
    gap_threshold: float
    cluster_tet_tolerance: float
    native_extractor: object = None
    gradient_c: object = None
    volumes_c: object = None
    tet_regions_c: object = None


def _cluster_tet_tolerance(tet_mesh: TetMesh, params: LayerParams, scale: float) -> float:
    """Taille caracteristique d'un tetraedre (mediane des longueurs d'arete x2), tolerance de proximite de
    `_tets_for_cluster` (Decision 7). PAS une fraction de `scale` : mesure sur le modele tripode, elle tombait
    sous la distance reelle centroide-tetraedre -> point de cluster le plus proche (0.328), filtrait TOUS les
    tetraedres et forcait un repli sur unique_tets entier (region_ids identiques pour tous les clusters)."""
    if params.split_disjoint_clusters is None:
        return 0.0
    tet_verts_per_tet = tet_mesh.tet_vertices[tet_mesh.tetrahedra]
    edge_lengths = np.linalg.norm(
        tet_verts_per_tet[:, 1:, :] - tet_verts_per_tet[:, 0:1, :], axis=-1
    )
    return float(np.median(edge_lengths)) * 2.0 if edge_lengths.size else scale


def _extract_layer_entries(isovalue: float, ctx: _ExtractContext) -> list:
    """Couche(s) brute(s) d'une isovaleur : liste de `(isovalue, contours, region_ids, local_frame)`."""
    native_extractor, params = ctx.native_extractor, ctx.params
    if (native_extractor is not None and (params.split_disjoint_clusters is None)
            and hasattr(native_extractor, "extract_layer")):
        # Chemin par defaut en un appel natif : isosurface, contours, regions et repere local.
        # Statut 2 (bord non simple) : on retombe sur trimesh ci-dessous.
        layer = native_extractor.extract_layer(float(isovalue), ctx.weld_tolerance, ctx.gradient_c,
                                               ctx.volumes_c, ctx.tet_regions_c)
        if isinstance(layer, tuple):
            contours, region_ids, local_frame = layer
            return [(float(isovalue), contours, list(region_ids), local_frame)]
        if layer == 0:
            return []
    if native_extractor is not None:
        welded_vertices, welded_faces, source_tets = native_extractor.extract(float(isovalue), ctx.weld_tolerance)
        if len(welded_faces) == 0:
            return []
    else:
        vertices, faces, source_tets = _extract_isosurface_triangles(
            ctx.tet_mesh, ctx.scalar_field.values, float(isovalue), ctx.tet_min, ctx.tet_max
        )
        if len(faces) == 0:
            return []
        welded_vertices, welded_faces = _weld_vertices(vertices, faces, ctx.weld_tolerance)
    contours = _extract_closed_contours(welded_vertices, welded_faces)
    if not contours:
        return []

    unique_tets = np.unique(source_tets)
    if params.split_disjoint_clusters is None or len(contours) < 2:
        # chemin historique inchange (mot pour mot) - non-regression par construction.
        region_ids = sorted({int(r) for r in ctx.tet_regions[unique_tets] if r >= 0})
        local_frame = _compute_local_frame(
            contours, ctx.scalar_field.gradient[unique_tets], ctx.volumes[unique_tets]
        )
        return [(float(isovalue), contours, region_ids, local_frame)]
    return _split_layer_clusters(isovalue, contours, unique_tets, ctx)


def _order_raw_layers(raw_layers: list, region_ranks: dict, gap_threshold: float, params: LayerParams) -> list:
    if params.split_disjoint_clusters is not None:
        # tri spatial DFS (Decision 8) : degenere exactement vers le tri historique (rank, isovalue) quand le
        # clustering ne produit jamais plus d'un cluster par isovaleur - pas de parametre separe.
        spatial_order = _spatial_dfs_order(raw_layers, region_ranks, gap_threshold)
        return [raw_layers[i] for i in spatial_order]

    def _sort_key(entry: tuple) -> tuple[int, float]:
        isovalue, _contours, region_ids, _frame = entry
        if params.isovalue_order:
            return (0, isovalue)
        rank = max((region_ranks.get(r, 0) for r in region_ids), default=0)
        return (rank, isovalue)

    return sorted(raw_layers, key=_sort_key)


def _split_layer_clusters(isovalue: float, contours: list, unique_tets: np.ndarray, ctx: _ExtractContext) -> list:
    """Une couche scindee en clusters spatialement disjoints (`split_disjoint_clusters`)."""
    tet_mesh, scalar_field = ctx.tet_mesh, ctx.scalar_field
    provisional_frame = _compute_local_frame(
        contours, scalar_field.gradient[unique_tets], ctx.volumes[unique_tets]
    )
    clusters = _cluster_contours(contours, provisional_frame, ctx.gap_threshold)
    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra[unique_tets]].mean(axis=1)
    entries = []
    for cluster_indices in clusters:
        cluster_contours = [contours[i] for i in cluster_indices]
        cluster_points = np.concatenate(cluster_contours, axis=0)
        cluster_tets = _tets_for_cluster(
            unique_tets, tet_centers, cluster_points, tolerance=ctx.cluster_tet_tolerance
        )
        if len(cluster_tets) == 0:
            # aucun tetraedre proche (tolerance trop stricte sur un maillage tres grossier) : repli sur
            # unique_tets entier plutot qu'un local_frame/region_ids vides.
            cluster_tets = unique_tets
        region_ids = sorted({int(r) for r in ctx.tet_regions[cluster_tets] if r >= 0})
        local_frame = _compute_local_frame(
            cluster_contours, scalar_field.gradient[cluster_tets], ctx.volumes[cluster_tets]
        )
        entries.append((float(isovalue), cluster_contours, region_ids, local_frame))
    return entries

from element_ardor.step5_marching import (  # noqa: E402,F401  (reexport : imports et patchs existants)
    _march_cap_batch,
    _march_quad_batch,
)
