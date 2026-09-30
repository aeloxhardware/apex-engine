"""
pipeline/step1b_junction_detection.py

Detection des zones de jonction multi-region a forte variance directionnelle
(CLAUDE.md limite #48, resout la limite #47 : non-monotonie du champ
scalaire dans une jonction ou plusieurs regions a directions divergentes
sont adjacentes sur une zone geometrique dense).

Ne fait AUCUN calcul de tetraedrisation ni de champ scalaire : prend en
entree le MeshData et le RegionDependencyGraph d'une premiere passe LEGERE
(step1 mesh+surface seul -> step2 -> step3, PAS de tetraedrisation fine),
et retourne des zones candidates a raffiner localement lors de la VRAIE
tetraedrisation (pipeline/step1_ingest.py, add_region generalise).

Fonction pure, testable en isolation - pas de dependance a TetMesh/
ScalarField, contrairement au reste du diagnostic (step5).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from core.geometry_utils import build_face_adjacency
from core.mesh_types import JunctionZone, MeshData, RegionDependencyGraph


def _angle_between_deg(a: np.ndarray, b: np.ndarray) -> float:
    """Angle en degres entre deux vecteurs unitaires - meme formule que
    pipeline/step3_segmentation.py::_angle_between_deg, dupliquee ici plutot
    qu'importee (fonction privee non exportee, 3 lignes, pas de dependance
    circulaire step1b<->step3 introduite pour si peu)."""
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0))))


def _merge_nearby_zones_once(zones: list[JunctionZone], merge_distance: float) -> list[JunctionZone]:
    """Une seule passe de fusion (voir `_merge_nearby_zones` pour le contexte
    complet et la discipline d'iteration jusqu'a convergence)."""
    n = len(zones)
    if n <= 1:
        return zones

    centroids = np.array([zone.centroid for zone in zones])
    radii = np.array([zone.radius for zone in zones])

    # Distance centre-a-centre entre toutes les paires (n petit compare au
    # nombre de faces du mesh - une matrice dense O(n^2) reste negligeable
    # meme pour plusieurs centaines de zones, contrairement au graphe de
    # faces qui, lui, doit rester creux).
    diff = centroids[:, None, :] - centroids[None, :, :]
    dist = np.linalg.norm(diff, axis=2)
    radius_sum = radii[:, None] + radii[None, :]
    connected = dist <= (radius_sum + merge_distance)
    np.fill_diagonal(connected, False)

    adjacency = sp.csr_matrix(connected)
    n_components, labels = sp.csgraph.connected_components(adjacency, directed=False)

    merged: list[JunctionZone] = []
    for label in range(n_components):
        group_indices = np.where(labels == label)[0]
        if len(group_indices) == 1:
            merged.append(zones[group_indices[0]])
            continue

        group_centroids = centroids[group_indices]
        group_radii = radii[group_indices]
        weights = group_radii / group_radii.sum()
        centroid = np.average(group_centroids, axis=0, weights=weights)
        # Rayon du groupe fusionne : englobe chaque zone membre a sa
        # distance reelle du nouveau centroide (pas une simple somme de
        # rayons, qui surestimerait largement un groupe etale).
        radius = float(
            np.max(np.linalg.norm(group_centroids - centroid, axis=1) + group_radii)
        )
        merged.append(JunctionZone(centroid=centroid, radius=radius))

    return merged


def _merge_nearby_zones(zones: list[JunctionZone], merge_distance: float) -> list[JunctionZone]:
    """
    Fusionne les `JunctionZone` dont les spheres d'influence se recouvrent
    ou sont proches (distance centre-a-centre < somme des rayons +
    `merge_distance`) en une seule zone par groupe connexe - ITERE jusqu'a
    convergence (nombre de zones stable), pas une seule passe.

    Necessaire sur une geometrie a jonctions denses (motif festonne, mesure
    sur death.stl : 301 zones detectees, 24 regions, distance mediane au
    plus proche voisin 1.18mm) - le graphe d'adjacence de FACES utilise par
    `detect_junction_zones` se fragmente en de nombreuses composantes des
    que des faces non marquees s'intercalent entre deux clusters de
    jonction physiquement tres proches, alors qu'un seul volume de
    raffinement les couvrirait tous. Sans cette fusion,
    `_merged_junction_seed` (pipeline/step1_ingest.py) moyenne 301
    centroides disperses sur toute une region du modele en un seul point
    qui ne represente plus aucune zone reelle - mesure directement : 0
    effet sur le nombre d'anomalies de non-monotonie (CLAUDE.md limite #48).

    L'iteration est necessaire (pas optionnelle) : une seule passe de
    `_merge_nearby_zones_once` peut laisser deux GROUPES fusionnes se
    chevaucher malgre eux, parce que le rayon d'un groupe fusionne
    (englobant chaque membre a sa distance reelle) peut rapprocher deux
    groupes qui ne se touchaient pas via leurs zones D'ORIGINE - mesure
    directement sur death.stl : une premiere passe (301 -> 2) laissait 2
    zones dont la distance centre-a-centre (28.8mm) restait BIEN EN DECA de
    la somme de leurs rayons (57.3mm, chevauchement net), simplement parce
    qu'elles provenaient de deux composantes de faces jamais directement
    adjacentes AVANT fusion. Sans re-verifier apres fusion, ce chevauchement
    residuel produisait un seed 'rest' impossible a placer hors des deux
    zones (`_find_rest_seed_outside_junction_zones` echouait avec
    `MeshIngestError` sur ce meme modele).
    """
    current = zones
    while True:
        merged = _merge_nearby_zones_once(current, merge_distance)
        if len(merged) == len(current):
            return merged
        current = merged


def detect_junction_zones(
    mesh_data: MeshData,
    region_dependency_graph: RegionDependencyGraph,
    angle_threshold_deg: float,
    margin: float,
    merge_distance: float = 5.0,
) -> list[JunctionZone]:
    """
    Identifie les clusters de faces frontiere entre deux regions dont les
    `target_direction` divergent au-dela de `angle_threshold_deg`, et
    retourne une zone (centroide + rayon) par cluster spatial - candidate a
    un raffinement local de maillage (`junction_refinement_resolution`,
    pipeline/step1_ingest.py) pour corriger l'erreur de discretisation FEM
    documentee limite #34/#36/#47 (un tetraedre P1 assez grand pour
    s'etendre sur une zone ou la direction cible change beaucoup ne peut
    representer qu'un seul gradient constant).

    Args:
        mesh_data: sortie d'une passe LEGERE de step1 (mesh de surface
            seul, tetraedrisation fine non necessaire ici).
        region_dependency_graph: sortie de step3 sur cette meme passe legere.
        angle_threshold_deg: seuil de divergence directionnelle entre deux
            regions adjacentes pour marquer leur frontiere comme jonction
            (meme role que SegmentationParams.max_region_merge_angle_deg,
            valeur separee car les deux poursuivent des objectifs
            differents : l'un empeche une fusion de region, l'autre
            declenche un raffinement de maillage).
        margin: marge (mm) ajoutee au rayon de chaque zone detectee, pour
            couvrir les tetraedres immediatement voisins de la frontiere
            elle-meme, pas seulement les faces exactement dessus.
        merge_distance: distance (mm) sous laquelle deux zones proches sont
            fusionnees en une seule (`_merge_nearby_zones`) - necessaire sur
            une geometrie a jonctions denses ou le graphe d'adjacence de
            FACES se fragmente en de nombreuses composantes separees par
            des faces non marquees, alors qu'un seul volume de raffinement
            les couvrirait tous physiquement (mesure sur death.stl : 301
            zones brutes, distance mediane au plus proche voisin 1.18mm -
            sans fusion, `_merged_junction_seed` moyenne des centroides
            disperses sur toute une region du modele en un point qui ne
            represente plus aucune zone reelle, mesure : 0 effet sur le
            nombre d'anomalies de non-monotonie).

    Returns:
        Liste de `JunctionZone` (centroide + rayon), APRES fusion des zones
        proches - une par cluster spatial de faces frontiere a forte
        divergence, regroupe si necessaire. Liste vide si aucune frontiere
        ne depasse le seuil (cas normal sur un modele a une seule region, ou
        a regions dont les directions restent proches).
    """
    face_ids_to_region: dict[int, int] = {}
    for region in region_dependency_graph.regions.values():
        for face_id in region.face_ids:
            face_ids_to_region[int(face_id)] = region.id

    face_count = mesh_data.face_count
    face_adjacency = build_face_adjacency(mesh_data.faces, face_count)
    face_adjacency_coo = face_adjacency.tocoo()

    # Marque les faces situees sur une frontiere a forte divergence
    # directionnelle (au moins une face adjacente d'une AUTRE region dont
    # l'angle entre target_direction depasse le seuil).
    is_junction_face = np.zeros(face_count, dtype=bool)
    junction_rows: list[int] = []
    junction_cols: list[int] = []

    for f1, f2 in zip(face_adjacency_coo.row, face_adjacency_coo.col):
        if f1 >= f2:
            continue  # matrice symetrique, ne traiter chaque paire qu'une fois
        r1 = face_ids_to_region.get(int(f1))
        r2 = face_ids_to_region.get(int(f2))
        if r1 is None or r2 is None or r1 == r2:
            continue
        dir1 = region_dependency_graph.regions[r1].target_direction
        dir2 = region_dependency_graph.regions[r2].target_direction
        if _angle_between_deg(dir1, dir2) > angle_threshold_deg:
            is_junction_face[f1] = True
            is_junction_face[f2] = True
            junction_rows.extend([f1, f2])
            junction_cols.extend([f2, f1])

    if not np.any(is_junction_face):
        return []

    # Regroupe les faces marquees en composantes connexes SUR LE SOUS-GRAPHE
    # des seules aretes de jonction detectees (pas l'adjacence complete) :
    # deux faces marquees mais reliees seulement par des faces non marquees
    # ne doivent pas fusionner en une seule zone.
    junction_adjacency = sp.csr_matrix(
        (np.ones(len(junction_rows)), (junction_rows, junction_cols)),
        shape=(face_count, face_count),
    )
    n_components, labels = sp.csgraph.connected_components(junction_adjacency, directed=False)

    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)

    zones: list[JunctionZone] = []
    for label in range(n_components):
        cluster_faces = np.where((labels == label) & is_junction_face)[0]
        if len(cluster_faces) == 0:
            continue
        cluster_centers = face_centers[cluster_faces]
        centroid = cluster_centers.mean(axis=0)
        radius = float(np.linalg.norm(cluster_centers - centroid, axis=1).max()) + margin
        zones.append(JunctionZone(centroid=centroid, radius=radius))

    return _merge_nearby_zones(zones, merge_distance)
