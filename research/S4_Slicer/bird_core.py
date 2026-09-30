"""
Coeur algorithmique de la methode S4-Slicer (Joshua Bird, GPL-3.0), extrait
des cellules 2, 3, 4, 7 du notebook main.ipynb en fonctions reutilisables,
sans etat module-level -- pour etre pilote depuis un serveur interactif
(server.py) qui garde le maillage courant en memoire entre les iterations.

SEUL changement par rapport a l'algorithme original :
- open3d (o3d.io.read_triangle_mesh) remplace par trimesh.load(...), car
  open3d n'a pas de wheel pour Python 3.14. Tableaux numeriquement
  equivalents.
- np.cross complete explicitement en 3D (numpy>=2.0 a retire le broadcast
  implicite 2D->3D).

Rien d'autre n'est modifie dans la logique de deformation elle-meme.
"""
import base64
import pickle

import networkx as nx
import numpy as np
import pyvista as pv
import tetgen
import trimesh
import scipy.sparse as sp
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation as R
from scipy.sparse import lil_matrix

UP_VECTOR = np.array([0, 0, 1])


def encode_object(obj):
    return base64.b64encode(pickle.dumps(obj)).decode("utf-8")


def decode_object(encoded_str):
    return pickle.loads(base64.b64decode(encoded_str))


# ============================================================================
# Etat du modele courant : regroupe tout ce que les fonctions de la cellule 3
# (calculate_path_length_to_base_gradient) et de la cellule 2 (voisinages)
# doivent partager -- equivalent des variables globales du notebook, mais
# explicite et passe en argument plutot qu'implicite.
# ============================================================================


class ModelState:
    """Tout ce qui ne change PAS d'une iteration a l'autre : topologie du
    maillage tetraedrique (voisinages, graphe de distances, cellules au
    contact du plateau). Construit une seule fois a l'ingestion."""

    def __init__(self, cell_neighbour_dict, cell_neighbour_graph, bottom_cells, bottom_cells_mask, cell_to_face):
        self.cell_neighbour_dict = cell_neighbour_dict
        self.cell_neighbour_graph = cell_neighbour_graph
        self.bottom_cells = bottom_cells
        self.bottom_cells_mask = bottom_cells_mask
        self.cell_to_face = cell_to_face


# ============================================================================
# CELLULE 2 (adaptee : trimesh au lieu d'open3d) -- ingestion + topologie
# ============================================================================


def load_and_tetrahedralize(stl_path, max_tet_volume=None, surface_subdivisions=0, flatten_bottom=None):
    """Charge un STL et le tetraedrise. Retourne le UnstructuredGrid brut
    (avant tout calcul d'attribut), recentre en XY et pose au sol en Z.

    `max_tet_volume` : plafond de volume par tetraedre, passe a TetGen via le
    switch `-a` (meme mecanisme que step1_ingest._tetrahedralize cote Apex).
    `None` = defaut TetGen, comportement historique inchange. Il raffine le
    CHAMP DE ROTATION (defini par tetraedre), donc la granularite avec
    laquelle la tete peut s'incliner le long de la piece.

    `surface_subdivisions` : nombre de passes de subdivision des triangles
    d'ENTREE avant tetraedrisation. Levier DISTINCT du precedent : TetGen
    ajoute des noeuds a l'interieur du volume mais ne subdivise jamais la
    peau qu'on lui donne, donc raffiner le volume ne rend pas la MESURE
    par face plus fine. Sur un STL a 28 faces (real_overhang), une rotation
    unique s'applique a des faces de 24 mm2 -- c'est la mesure qui est
    grossiere, pas seulement le champ.

    `flatten_bottom` : epaisseur (mm) sous laquelle les sommets du dessous
    sont ramenes sur le plateau AVANT tetraedrisation (mesh_flatten.py, qui
    ecarte les sommets dont l'ecrasement serait faux) -- en pratique la
    hauteur de la premiere couche, qui est de toute facon ecrasee a
    l'impression. Un dessous presque plan mais etage (tilted benchy
    undertriangulated : sommets a 0 ; 0,004 ; 0,01 ; ... 0,13 mm) faisait
    generer a TetGen 21 004 tetraedres microscopiques (2,6e-6 mm3 en mediane,
    6 600 fois plus petits qu'ailleurs, 43 % du maillage) pour remplir ce
    coin extremement mince. None = geometrie inchangee."""
    mesh = trimesh.load(stl_path)
    if flatten_bottom is not None:
        from mesh_flatten import flatten_bottom as _flatten_bottom
        mesh, _displacement = _flatten_bottom(mesh, flatten_bottom)

    for _ in range(surface_subdivisions):
        mesh = mesh.subdivide()

    input_tet = tetgen.TetGen(np.asarray(mesh.vertices), np.asarray(mesh.faces))
    if max_tet_volume is not None:
        input_tet.tetrahedralize(switches=f"pq1.414a{max_tet_volume}")
    else:
        input_tet.tetrahedralize()
    input_tet = input_tet.grid

    x_min, x_max, y_min, y_max, z_min, z_max = input_tet.bounds
    input_tet.points -= np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min])

    return input_tet


def build_topology(tet):
    """Construit les voisinages (point/edge/face) et le graphe de distances
    euclidiennes entre centres de cellules -- une seule fois par maillage,
    jamais recalcule entre deux iterations de deformation."""
    cell_neighbour_dict = {
        neighbour_type: {face: [] for face in range(tet.number_of_cells)}
        for neighbour_type in ["point", "edge", "face"]
    }
    for neighbour_type in ["point", "edge", "face"]:
        cell_neighbours = []
        for cell_index in range(tet.number_of_cells):
            neighbours = tet.cell_neighbors(cell_index, f"{neighbour_type}s")
            for neighbour in neighbours:
                if neighbour > cell_index:
                    cell_neighbours.append((cell_index, neighbour))
        for face_1, face_2 in np.array(cell_neighbours):
            cell_neighbour_dict[neighbour_type][face_1].append(face_2)
            cell_neighbour_dict[neighbour_type][face_2].append(face_1)

        tet.field_data[f"cell_{neighbour_type}_neighbours"] = np.array(cell_neighbours)

    cell_neighbour_graph = nx.Graph()
    cell_centers = tet.cell_centers().points
    for edge in tet.field_data["cell_point_neighbours"]:
        distance = np.linalg.norm(cell_centers[edge[0]] - cell_centers[edge[1]])
        cell_neighbour_graph.add_weighted_edges_from([(edge[0], edge[1], distance)])

    return cell_neighbour_dict, cell_neighbour_graph


def update_tet_attributes(tet, state):
    """CELLULE 2 (update_tet_attributes) -- inchangee. Recalcule normales,
    centres, overhang_angle, in_air a partir de la geometrie COURANTE du
    tet (donc a rappeler apres chaque deformation)."""
    surface_mesh = tet.extract_surface()
    cell_to_face = state.cell_to_face

    cells = tet.cells.reshape(-1, 5)[:, 1:]
    tet.add_field_data(cells, "cells")
    cell_vertices = tet.points
    tet.add_field_data(cell_vertices, "cell_vertices")
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]
    tet.add_field_data(faces, "faces")
    face_vertices = surface_mesh.points
    tet.add_field_data(face_vertices, "face_vertices")

    tet.cell_data["face_normal"] = np.full((tet.number_of_cells, 3), np.nan)
    surface_mesh_face_normals = surface_mesh.face_normals
    for cell_index, face_indices in cell_to_face.items():
        face_normals = surface_mesh_face_normals[face_indices]
        most_down_normal_index = np.argmin(face_normals[:, 2])
        tet.cell_data["face_normal"][cell_index] = face_normals[most_down_normal_index]
    tet.cell_data["face_normal"] = (
        tet.cell_data["face_normal"] / np.linalg.norm(tet.cell_data["face_normal"], axis=1)[:, None]
    )

    tet.cell_data["face_center"] = np.empty((tet.number_of_cells, 3))
    tet.cell_data["face_center"][:, :] = np.nan
    surface_mesh_cell_centers = surface_mesh.cell_centers().points
    for cell_index, face_indices in cell_to_face.items():
        face_centers = surface_mesh_cell_centers[face_indices]
        most_down_center_index = np.argmin(face_centers[:, 2])
        tet.cell_data["face_center"][cell_index] = face_centers[most_down_center_index]

    tet.cell_data["cell_center"] = tet.cell_centers().points

    bottom_cell_threshold = np.nanmin(tet.cell_data["face_center"][:, 2]) + 0.3
    bottom_cells_mask_local = tet.cell_data["face_center"][:, 2] < bottom_cell_threshold
    tet.cell_data["is_bottom"] = bottom_cells_mask_local
    bottom_cells_local = np.where(bottom_cells_mask_local)[0]

    face_normals = tet.cell_data["face_normal"].copy()
    face_normals[bottom_cells_mask_local] = np.nan
    overhang_angle = np.arccos(np.dot(face_normals, UP_VECTOR))
    tet.cell_data["overhang_angle"] = overhang_angle

    overhang_direction = face_normals[:, :2].copy()
    overhang_direction /= np.linalg.norm(overhang_direction, axis=1)[:, None]
    tet.cell_data["overhang_direction"] = overhang_direction

    IN_AIR_THRESHOLD = 1
    tet.cell_data["in_air"] = np.full(tet.number_of_cells, False)

    _, paths_to_bottom = nx.multi_source_dijkstra(state.cell_neighbour_graph, set(bottom_cells_local))

    tet.cell_data["path_to_bottom"] = np.full(
        (tet.number_of_cells, np.max([len(x) for x in paths_to_bottom.values()])), -1
    )
    for cell_index, path_to_bottom in paths_to_bottom.items():
        tet.cell_data["path_to_bottom"][cell_index, : len(path_to_bottom)] = path_to_bottom

    for cell_index in range(tet.number_of_cells):
        path_to_bottom = paths_to_bottom[cell_index]
        if len(path_to_bottom) > 1:
            cell_heights = tet.cell_data["cell_center"][path_to_bottom, 2]
            if np.any(cell_heights > tet.cell_data["cell_center"][cell_index, 2] + IN_AIR_THRESHOLD):
                tet.cell_data["in_air"][cell_index] = True

    return tet


def calculate_tet_attributes(tet):
    """CELLULE 2 (calculate_tet_attributes) -- inchangee. Calcule la table
    cell_to_face/face_to_cell (UNE SEULE FOIS, topologie invariante), puis
    delegue a update_tet_attributes pour le reste. Retourne aussi un
    ModelState pret a l'emploi pour les iterations suivantes."""
    surface_mesh = tet.extract_surface()

    cells = tet.cells.reshape(-1, 5)[:, 1:]
    tet.add_field_data(cells, "cells")
    cell_vertices = tet.points
    tet.add_field_data(cell_vertices, "cell_vertices")
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]
    tet.add_field_data(faces, "faces")
    face_vertices = surface_mesh.points
    tet.add_field_data(face_vertices, "face_vertices")

    cell_to_face = {}
    face_to_cell = {face_index: [] for face_index in range(len(faces))}
    cell_to_face_vertices = {}
    face_to_cell_vertices = {}
    for cell_vertex_index, cell_vertex in enumerate(tet.field_data["cell_vertices"].reshape(-1, 3)):
        face_vertex_index = np.where((face_vertices == cell_vertex).all(axis=1))[0]
        if len(face_vertex_index) == 1:
            cell_to_face_vertices[cell_vertex_index] = face_vertex_index[0]
            face_to_cell_vertices[face_vertex_index[0]] = cell_vertex_index

    for cell_index, cell in enumerate(tet.field_data["cells"]):
        face_vertex_indices = [
            cell_to_face_vertices[cell_vertex_index]
            for cell_vertex_index in cell
            if cell_vertex_index in cell_to_face_vertices
        ]
        if len(face_vertex_indices) >= 3:
            extracted = surface_mesh.extract_points(face_vertex_indices, adjacent_cells=False)
            if extracted.number_of_cells >= 1:
                cell_to_face[cell_index] = list(extracted.cell_data["vtkOriginalCellIds"])
                for face_index in extracted.cell_data["vtkOriginalCellIds"]:
                    face_to_cell[face_index].append(cell_index)

    tet.add_field_data(encode_object(cell_to_face), "cell_to_face")
    tet.add_field_data(encode_object(face_to_cell), "face_to_cell")

    tet.cell_data["has_face"] = np.zeros(tet.number_of_cells)
    for cell_index, face_indices in cell_to_face.items():
        tet.cell_data["has_face"][cell_index] = 1

    cell_neighbour_dict, cell_neighbour_graph = build_topology(tet)

    # bottom_cells calcule a blanc pour amorcer update_tet_attributes (qui en
    # a besoin via state.cell_neighbour_graph, deja pret) -- meme sequence
    # que le notebook (calculate_tet_attributes appelle update_tet_attributes
    # en interne).
    dummy_state = ModelState(cell_neighbour_dict, cell_neighbour_graph, None, None, cell_to_face)
    tet = update_tet_attributes(tet, dummy_state)

    bottom_cells_mask = tet.cell_data["is_bottom"]
    bottom_cells = np.where(bottom_cells_mask)[0]

    tet.cell_data["overhang_angle"][bottom_cells] = np.nan

    state = ModelState(cell_neighbour_dict, cell_neighbour_graph, bottom_cells, bottom_cells_mask, cell_to_face)
    return tet, state


# ============================================================================
# CELLULE 3 -- inchangee
# ============================================================================


def planeFit(points):
    from numpy.linalg import svd

    points = np.reshape(points, (np.shape(points)[0], -1))
    ctr = points.mean(axis=1)
    x = points - ctr[:, np.newaxis]
    M = np.dot(x, x.T)
    return ctr, svd(M)[0][:, -1]


def calculate_path_length_to_base_gradient(tet, state, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO):
    cell_neighbour_dict = state.cell_neighbour_dict
    cell_neighbour_graph = state.cell_neighbour_graph
    bottom_cells = state.bottom_cells

    path_length_to_base_gradient = np.zeros((tet.number_of_cells))

    cell_distance_to_bottom = np.empty((tet.number_of_cells))
    cell_distance_to_bottom[:] = np.nan
    distances_to_bottom, paths_to_bottom = nx.multi_source_dijkstra(cell_neighbour_graph, set(bottom_cells))
    closest_bottom_cell_indices = np.zeros((tet.number_of_cells), dtype=int)
    for cell_index in range(tet.number_of_cells):
        face_normal = tet.cell_data["face_normal"][cell_index]

        cell_is_overhang = np.arccos(np.dot(face_normal, [0, 0, 1])) > np.deg2rad(90 + MAX_OVERHANG)
        if cell_is_overhang and cell_index not in bottom_cells:
            closest_bottom_cell_indices[cell_index] = paths_to_bottom[cell_index][0]
            cell_distance_to_bottom[cell_index] = distances_to_bottom[cell_index]

    tet.cell_data["cell_distance_to_bottom"] = cell_distance_to_bottom

    for cell_index in range(tet.number_of_cells):
        if not np.isnan(cell_distance_to_bottom[cell_index]):
            local_cells = cell_neighbour_dict["edge"][cell_index]
            local_cells = np.hstack((local_cells, cell_index))

            local_cell_path_lengths = [cell_distance_to_bottom[local_cell] for local_cell in local_cells]
            local_cell_path_lengths = np.array(local_cell_path_lengths)

            local_cells = np.array(local_cells)[~np.isnan(local_cell_path_lengths)]
            local_cell_path_lengths = local_cell_path_lengths[~np.isnan(local_cell_path_lengths)]

            if len(local_cell_path_lengths) < 3:
                location_to_roll_to = tet.cell_data["cell_center"][closest_bottom_cell_indices[cell_index], :2]

                direction_to_bottom = location_to_roll_to - tet.cell_data["cell_center"][cell_index, :2]
                direction_to_bottom /= np.linalg.norm(direction_to_bottom)

                cell_center = tet.cell_data["cell_center"][cell_index, :2].copy()
                cell_center /= np.linalg.norm(cell_center)

                optimal_rotation_direction = np.dot(cell_center, direction_to_bottom) / np.abs(
                    np.dot(cell_center, direction_to_bottom)
                )
                if np.isnan(optimal_rotation_direction):
                    optimal_rotation_direction = 0

                path_length_to_base_gradient[cell_index] = optimal_rotation_direction

            else:
                points = np.hstack((tet.cell_data["cell_center"][local_cells, :2], local_cell_path_lengths[:, None]))
                _, plane_normal = planeFit(points.T)

                cell_center_direction_normalized = tet.cell_data["cell_center"][cell_index, :2] / np.linalg.norm(
                    tet.cell_data["cell_center"][cell_index, :2]
                )
                gradient_in_radial_direction = np.dot(cell_center_direction_normalized, plane_normal[:2])

                if np.isnan(gradient_in_radial_direction):
                    gradient_in_radial_direction = np.mean(
                        path_length_to_base_gradient[local_cells][~np.isnan(path_length_to_base_gradient[local_cells])]
                    )
                    if np.isnan(gradient_in_radial_direction):
                        gradient_in_radial_direction = 0

                path_length_to_base_gradient[cell_index] = gradient_in_radial_direction

    if INITIAL_ROTATION_FIELD_SMOOTHING != 0:
        for _ in range(INITIAL_ROTATION_FIELD_SMOOTHING):
            smoothed_path_length_to_base_gradient = np.zeros((tet.number_of_cells))
            for cell_index in range(tet.number_of_cells):
                if path_length_to_base_gradient[cell_index] != 0:
                    neighbours = cell_neighbour_dict["point"][cell_index]
                    local_cells = neighbours.copy()
                    for neighbour in neighbours:
                        local_cells.extend(cell_neighbour_dict["point"][neighbour])
                    local_cells = np.array(list(set(local_cells)))
                    local_cells = local_cells[path_length_to_base_gradient[local_cells] != 0]
                    smoothed_path_length_to_base_gradient[cell_index] = np.mean(
                        path_length_to_base_gradient[local_cells]
                    )

        path_length_to_base_gradient = smoothed_path_length_to_base_gradient

    if not SET_INITIAL_ROTATION_TO_ZERO:
        path_length_to_base_gradient[path_length_to_base_gradient == 0] = np.nan
    tet.cell_data["path_length_to_base_gradient"] = path_length_to_base_gradient

    return path_length_to_base_gradient


# ============================================================================
# CELLULE 4 -- inchangee (fix numpy>=2.0 sur np.cross)
# ============================================================================


def calculate_initial_rotation_field(
    tet, state, MAX_OVERHANG, ROTATION_MULTIPLIER, STEEP_OVERHANG_COMPENSATION,
    INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION,
):
    initial_rotation_field = np.abs(np.deg2rad(90 + MAX_OVERHANG) - tet.cell_data["overhang_angle"])

    path_length_to_base_gradient = calculate_path_length_to_base_gradient(
        tet, state, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO
    )

    if STEEP_OVERHANG_COMPENSATION:
        initial_rotation_field[tet.cell_data["in_air"]] += 2 * (
            np.deg2rad(180) - tet.cell_data["overhang_angle"][tet.cell_data["in_air"]]
        )

    initial_rotation_field *= path_length_to_base_gradient

    initial_rotation_field = np.clip(initial_rotation_field * ROTATION_MULTIPLIER, -np.deg2rad(360), np.deg2rad(360))

    initial_rotation_field = np.clip(initial_rotation_field, MAX_NEG_ROTATION, MAX_POS_ROTATION)

    tet.cell_data["initial_rotation_field"] = initial_rotation_field

    return initial_rotation_field


def calculate_rotation_matrices(tet, rotation_field):
    # FIX numpy>=2.0 : np.cross exige des vecteurs 3D des deux cotes (le
    # broadcast implicite 2D->3D, z=0, a ete retire).
    xy = tet.cell_data["cell_center"][:, :2]
    xy_3d = np.hstack([xy, np.zeros((xy.shape[0], 1))])
    tangential_vectors = np.cross(np.array([0, 0, 1]), xy_3d)
    tangential_vectors /= np.linalg.norm(tangential_vectors, axis=1)[:, None]
    tangential_vectors[np.isnan(tangential_vectors).any(axis=1)] = [1, 0, 0]

    rotation_matrices = R.from_rotvec(rotation_field[:, None] * tangential_vectors).as_matrix()

    return rotation_matrices


def calculate_unique_vertices_rotated(tet, rotation_field):
    rotation_matrices = calculate_rotation_matrices(tet, rotation_field)

    unique_vertices = np.zeros((tet.number_of_cells, 4, 3))
    for cell_index, cell in enumerate(tet.field_data["cells"]):
        unique_vertices[cell_index] = tet.field_data["cell_vertices"][cell]

    cell_centers = tet.cell_data["cell_center"]

    unique_vertices_rotated = cell_centers.reshape(-1, 1, 3, 1) + rotation_matrices.reshape(-1, 1, 3, 3) @ (
        unique_vertices.reshape(-1, 4, 3, 1) - cell_centers.reshape(-1, 1, 3, 1)
    )

    return unique_vertices_rotated


def apply_rotation_field_unique_vertices(tet, rotation_field):
    unique_vertices_rotated = calculate_unique_vertices_rotated(tet, rotation_field)

    unique_cells = np.zeros((tet.number_of_cells, 5), dtype=int)
    unique_cells[:, 0] = 4
    unique_cells[:, 1:] = np.arange(tet.number_of_cells * 4).reshape(-1, 4)

    new_tet = pv.UnstructuredGrid(
        unique_cells.flatten(), np.full(tet.number_of_cells, pv.CellType.TETRA), unique_vertices_rotated.reshape(-1, 3)
    )

    return new_tet


def planar_base_mask(tet, overhang_deg, layer_height, max_height=None):
    """
    Base a imprimer A PLAT (non deformee) : tout le volume de la piece
    d'ORIGINE situe sous son surplomb le plus bas, avec au minimum les
    cellules posees sur le plateau.

    Surplomb = meme critere que Bird : face exterieure dont l'angle a la
    verticale montante depasse 90 deg + overhang_deg, faces posees sur le
    plateau exclues. La hauteur du surplomb est celle du SOMMET le plus bas
    des faces en surplomb (la ou elles commencent), pas de leur centre ; une
    cellule appartient a la base si TOUS ses noeuds sont a cette hauteur ou
    en dessous. Aucune constante propre a un modele : la hauteur de la base
    decoule de la geometrie de chaque piece.

    Mesure qui a fixe ce choix (Benchy incline) : le rond fait 0,5 mm et le
    pilier incline qui le surmonte est un vrai surplomb des 0,5 mm. Prendre
    le centre des faces, moins une marge d'une couche, coupait le rond en
    deux (base de 0,3 mm) ; le sommet le plus bas, sans marge, le garde entier.

    La base fait au moins une hauteur de couche (`layer_height`), meme si un
    surplomb commence plus bas.

    `max_height` (optionnel) plafonne la base a cette hauteur au-dessus du
    plateau : variante "contact + marge", ou seules les premieres couches sont
    figees et Bird reste libre au-dessus, meme sous le premier surplomb.

    A appeler UNE FOIS, sur le tet d'origine (calculate_tet_attributes).
    Rend (masque des cellules, masque des noeuds, z du plus bas surplomb).
    """
    cells = np.asarray(tet.field_data["cells"])
    points = np.asarray(tet.points)
    tri = np.asarray(tet.field_data["face_vertices"])[np.asarray(tet.field_data["faces"])]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.maximum(np.linalg.norm(normals, axis=1), 1e-12)[:, None]
    z0 = points[:, 2].min()
    # Face posee sur le plateau = entierement dans la PREMIERE COUCHE et
    # tournee vers le bas : elle s'imprime de toute facon comme premiere
    # couche. Une tolerance de 1e-3 mm ne suffisait pas : le dessous d'un STL
    # n'est pas toujours plan (tilted benchy undertriangulated : sommets
    # etages de 0 a 0,12 mm) ; aucune face n'y etait reconnue et le dessous
    # passait pour un surplomb a 0,00 mm.
    on_plate = (tri[:, :, 2].max(axis=1) < z0 + layer_height) & (normals[:, 2] < 0.0)
    angle = np.degrees(np.arccos(np.clip(normals[:, 2], -1.0, 1.0)))
    overhanging = (angle > 90.0 + overhang_deg) & ~on_plate
    if overhanging.any():
        z_first = float(tri[overhanging][:, :, 2].min())
    else:
        z_first = float(points[:, 2].max())
    top = max(z_first, z0 + layer_height)
    if max_height is not None:
        top = min(top, z0 + max(max_height, layer_height))
    # Ce sont les NOEUDS sous la hauteur de base qui sont figes ; une cellule
    # recoit une rotation nulle si tous ses noeuds le sont. Ajouter en bloc
    # les cellules "du bas" de Bird (tolerance 0,3 mm) entrainait des
    # tetraedres hauts de plusieurs mm (45 % de la piece figee sur ce modele).
    node_mask = points[:, 2] <= top + 1e-6
    cell_mask = node_mask[cells].all(axis=1)
    return cell_mask, node_mask, z_first


def optimize_rotations(
    tet, state, NEIGHBOUR_LOSS_WEIGHT, MAX_OVERHANG, ROTATION_MULTIPLIER, ITERATIONS,
    STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO,
    MAX_POS_ROTATION, MAX_NEG_ROTATION, progress_callback=None, fixed_zero_cells=None,
):
    """fixed_zero_cells : masque de cellules dont la rotation est FIXEE a 0
    (base plane, cf. planar_base_mask). Contrainte dure : ces cellules sont
    retirees des inconnues, elles n'entrent plus que comme constantes dans le
    terme de lissage entre voisins. Sans elle, la base -- qui n'a aucune
    cible propre (surplomb NaN) -- heritait par lissage de la rotation du
    pilier et se creusait en cuvette (mesure sur le Benchy : 6,75 mm)."""
    initial_rotation_field = calculate_initial_rotation_field(
        tet, state, MAX_OVERHANG, ROTATION_MULTIPLIER, STEEP_OVERHANG_COMPENSATION,
        INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION,
    )
    if fixed_zero_cells is not None and fixed_zero_cells.any():
        return _optimize_rotations_with_fixed(
            tet, initial_rotation_field, NEIGHBOUR_LOSS_WEIGHT, ITERATIONS,
            np.asarray(fixed_zero_cells, dtype=bool), progress_callback,
        )
    num_cells_with_initial_rotation = np.sum(~np.isnan(initial_rotation_field))

    # FIX PERF (meme defaut que calculate_deformation ci-dessous, cf. session
    # de profilage) : jac_sparsity() bouclait meme en PUR PYTHON
    # ("for i, (cell_1, cell_2) in enumerate(...)") sur potentiellement des
    # dizaines de milliers d'aretes -- pire que le lil_matrix fantaisie de
    # objective_jacobian. La sparsite (structure des lignes/colonnes non
    # nulles) est INVARIANTE d'une evaluation a l'autre ; seules les valeurs
    # changent -- indices COO precalcules UNE SEULE FOIS, hors des closures.
    cell_face_neighbours = tet.field_data["cell_face_neighbours"]
    _n_neighbours = len(cell_face_neighbours)
    _cell_1 = cell_face_neighbours[:, 0]
    _cell_2 = cell_face_neighbours[:, 1]
    _valid_cell_indices = np.where(~np.isnan(initial_rotation_field))[0]

    # ligne i (0..n_neighbours-1) : deux entrees non nulles, cell_1 et cell_2.
    # lignes n_neighbours..fin : une seule entree (contrainte de cible directe).
    _rot_rows = np.concatenate([
        np.arange(_n_neighbours), np.arange(_n_neighbours),
        _n_neighbours + np.arange(len(_valid_cell_indices)),
    ])
    _rot_cols = np.concatenate([_cell_1, _cell_2, _valid_cell_indices])

    def objective_function(rotation_field):
        neighbour_differences = rotation_field[cell_face_neighbours[:, 0]] - rotation_field[cell_face_neighbours[:, 1]]
        neighbour_losses = NEIGHBOUR_LOSS_WEIGHT * neighbour_differences**2

        initial_rotation_losses = (
            rotation_field[_valid_cell_indices] - initial_rotation_field[_valid_cell_indices]
        ) ** 2

        return np.concatenate((neighbour_losses, initial_rotation_losses))

    def objective_jacobian(rotation_field):
        differences = rotation_field[_cell_1] - rotation_field[_cell_2]
        target_diff = 2 * (rotation_field[_valid_cell_indices] - initial_rotation_field[_valid_cell_indices])

        jac_data = np.concatenate([
            2 * NEIGHBOUR_LOSS_WEIGHT * differences,
            -2 * NEIGHBOUR_LOSS_WEIGHT * differences,
            target_diff,
        ]).astype(np.float32)

        jac = sp.coo_matrix(
            (jac_data, (_rot_rows, _rot_cols)),
            shape=(_n_neighbours + num_cells_with_initial_rotation, tet.number_of_cells),
        )
        return jac.tocsr()

    def jac_sparsity():
        sparsity_data = np.ones(len(_rot_rows), dtype=np.int8)
        sparsity = sp.coo_matrix(
            (sparsity_data, (_rot_rows, _rot_cols)),
            shape=(_n_neighbours + num_cells_with_initial_rotation, tet.number_of_cells),
        )
        return sparsity.tocsr()

    smoothed_rotation_field = np.zeros((tet.number_of_cells))

    def _callback(intermediate_result):
        # Signature scipy>=1.16 : intermediate_result est un OptimizeResult
        # portant nfev/cost a CETTE iteration externe de trf -- appele une
        # fois par iteration, jamais a l'interieur d'une evaluation de
        # fonction/jacobienne, donc cout negligeable meme sur un maillage dense.
        if progress_callback is not None:
            progress_callback("rotation", intermediate_result.nfev, ITERATIONS, float(intermediate_result.cost))

    result = least_squares(
        objective_function,
        smoothed_rotation_field,
        jac=objective_jacobian,
        max_nfev=ITERATIONS,
        jac_sparsity=jac_sparsity(),
        verbose=0,
        method="trf",
        ftol=1e-6,
        callback=_callback if progress_callback is not None else None,
    )

    return result.x, result.cost


def _optimize_rotations_with_fixed(tet, initial_rotation_field, NEIGHBOUR_LOSS_WEIGHT, ITERATIONS,
                                   fixed, progress_callback=None):
    """Meme probleme que optimize_rotations (memes residus, meme solveur),
    avec les cellules `fixed` retirees des inconnues et figees a 0."""
    n_cells = tet.number_of_cells
    free = np.flatnonzero(~fixed)
    column = np.full(n_cells, -1)
    column[free] = np.arange(len(free))

    pairs = tet.field_data["cell_face_neighbours"]
    pairs = pairs[~(fixed[pairs[:, 0]] & fixed[pairs[:, 1]])]      # paires figees : constantes
    c1, c2 = pairs[:, 0], pairs[:, 1]
    targets = np.flatnonzero(~np.isnan(initial_rotation_field) & ~fixed)
    n_pairs = len(pairs)

    rows, cols, which = [], [], []
    for side, cells in ((0, c1), (1, c2)):
        keep = column[cells] >= 0
        rows.append(np.flatnonzero(keep))
        cols.append(column[cells[keep]])
        which.append(np.full(keep.sum(), side))
    rows.append(n_pairs + np.arange(len(targets)))
    cols.append(column[targets])
    rows, cols, which = np.concatenate(rows), np.concatenate(cols), np.concatenate(which)
    n_rows = n_pairs + len(targets)
    shape = (n_rows, len(free))

    def full(x):
        field = np.zeros(n_cells)
        field[free] = x
        return field

    def objective_function(x):
        field = full(x)
        d = field[c1] - field[c2]
        return np.concatenate((NEIGHBOUR_LOSS_WEIGHT * d**2,
                               (field[targets] - initial_rotation_field[targets]) ** 2))

    def objective_jacobian(x):
        field = full(x)
        d = field[c1] - field[c2]
        pair_rows = rows[: len(rows) - len(targets)]
        pair_which = which[: len(rows) - len(targets)]
        pair_values = np.where(pair_which == 0, 2 * NEIGHBOUR_LOSS_WEIGHT, -2 * NEIGHBOUR_LOSS_WEIGHT) * d[pair_rows]
        target_values = 2 * (field[targets] - initial_rotation_field[targets])
        data = np.concatenate([pair_values, target_values]).astype(np.float32)
        return sp.coo_matrix((data, (rows, cols)), shape=shape).tocsr()

    sparsity = sp.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=shape).tocsr()

    def _callback(intermediate_result):
        if progress_callback is not None:
            progress_callback("rotation", intermediate_result.nfev, ITERATIONS, float(intermediate_result.cost))

    result = least_squares(
        objective_function, np.zeros(len(free)), jac=objective_jacobian, max_nfev=ITERATIONS,
        jac_sparsity=sparsity, verbose=0, method="trf", ftol=1e-6,
        callback=_callback if progress_callback is not None else None,
    )
    return full(result.x), result.cost


# ============================================================================
# CELLULE 7 -- inchangee
# ============================================================================

_N = np.eye(4) - 1 / 4 * np.ones((4, 4))


def calculate_deformation(tet, rotation_field, ITERATIONS, progress_callback=None, pinned=None):
    """pinned : (masque des noeuds, positions imposees) -- noeuds FIXES a ces
    positions (base plane : positions d'ORIGINE, a chaque passe). Contrainte
    dure : retires des inconnues, ils n'entrent plus que comme constantes."""
    if pinned is not None and np.asarray(pinned[0]).any():
        return _calculate_deformation_pinned(tet, rotation_field, ITERATIONS, pinned, progress_callback)
    new_vertices = tet.points.copy()
    params = new_vertices.flatten()

    rotation_matrices = calculate_rotation_matrices(tet, rotation_field)

    old_vertices = tet.field_data["cell_vertices"][tet.field_data["cells"]]
    old_vertices_transformed = np.einsum("ijk,ikl->ijl", rotation_matrices, (_N @ old_vertices).transpose(0, 2, 1))

    # FIX PERF (mesure : ~36% du temps de chaque evaluation, cf. session de
    # profilage) : objective_jacobian/jac_sparsity construisaient un
    # lil_matrix rempli par assignation fantaisie ("J[rows, cols] = valeurs"),
    # le point faible connu de ce format -- l'indexation fantaisie de lil_matrix
    # n'est PAS vectorisee en interne, elle boucle element par element cote
    # scipy. La sparsite (quels (ligne, colonne) sont non nuls) est INVARIANTE
    # d'une evaluation a l'autre -- seules les VALEURS changent -- donc les
    # indices COO sont precalcules UNE SEULE FOIS ici, en dehors des deux
    # fonctions rappelees par least_squares a chaque iteration.
    num_cells_per_vertex = len(tet.field_data["cells"][0])  # = 4 (tetraedre)
    _cell_indices = np.repeat(np.arange(tet.number_of_cells), num_cells_per_vertex)
    _vertex_indices = np.ravel(tet.field_data["cells"])
    # Un triplet (ligne=cellule, colonne=vertex*3+dim) par dimension x,y,z :
    # meme structure que la boucle "for dim in range(3)" d'origine, mais les
    # 3 blocs de colonnes sont concatenes en un seul assemblage COO vectorise.
    _jac_rows = np.tile(_cell_indices, 3)
    _jac_cols = np.concatenate([_vertex_indices * 3 + dim for dim in range(3)])

    def objective_function(params):
        new_vertices = params[: tet.number_of_points * 3].reshape(-1, 3)
        new_vertices_transformed = (_N @ new_vertices[tet.field_data["cells"]]).transpose(0, 2, 1)
        position_losses = np.linalg.norm(new_vertices_transformed - old_vertices_transformed, axis=(1, 2)) ** 2
        return position_losses

    def objective_jacobian(params):
        new_vertices = params[: tet.number_of_points * 3].reshape(-1, 3)
        new_vertices_transformed = (_N @ new_vertices[tet.field_data["cells"]]).transpose(0, 2, 1)
        diff = new_vertices_transformed - old_vertices_transformed
        diff = diff.transpose(0, 2, 1)
        # meme valeurs que l'original (2 * diff[:, :, dim].ravel() pour
        # dim=0,1,2), assemblees dans le meme ordre que _jac_cols ci-dessus.
        jac_data = np.concatenate([2 * diff[:, :, dim].ravel() for dim in range(3)]).astype(np.float32)
        J = sp.coo_matrix(
            (jac_data, (_jac_rows, _jac_cols)), shape=(tet.number_of_cells, len(params))
        )
        return J.tocsr()

    def jac_sparsity():
        sparsity_data = np.ones(len(_jac_rows), dtype=np.int8)
        sparsity = sp.coo_matrix(
            (sparsity_data, (_jac_rows, _jac_cols)), shape=(tet.number_of_cells, len(params))
        )
        return sparsity.tocsr()

    def _callback(intermediate_result):
        if progress_callback is not None:
            progress_callback("deformation", intermediate_result.nfev, ITERATIONS, float(intermediate_result.cost))

    result = least_squares(
        objective_function,
        params,
        max_nfev=ITERATIONS,
        verbose=0,
        jac=objective_jacobian,
        jac_sparsity=jac_sparsity(),
        method="trf",
        x_scale="jac",
        callback=_callback if progress_callback is not None else None,
    )

    return result.x[: tet.number_of_points * 3].reshape(-1, 3), result.cost


def normal_constraint_rows(original, cells):
    """Pour chaque tetraedre, les coefficients (4,) qui donnent, a partir des
    4 hauteurs DEFORMEES de ses noeuds, la composante radiale (a_r) et
    tangentielle (a_t) du gradient de la hauteur deformee par rapport a la
    geometrie d'ORIGINE -- c'est la normale de couche finale dans la piece.
    g = M^-1 D z4 (M : aretes d'origine en lignes, D : differences)."""
    p = original[cells]
    m = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=1)
    d = np.array([[-1.0, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]])
    g_op = np.linalg.solve(m, np.broadcast_to(d, (len(cells), 3, 4)))          # (T,3,4)
    c = p.mean(axis=1)[:, :2]
    r = np.column_stack([c, np.zeros(len(c))])
    r /= np.maximum(np.linalg.norm(r, axis=1, keepdims=True), 1e-9)
    t = np.cross([0.0, 0.0, 1.0], r)
    return np.einsum("tk,tkj->tj", r, g_op), np.einsum("tk,tkj->tj", t, g_op)


def _calculate_deformation_pinned(tet, rotation_field, ITERATIONS, pinned, progress_callback=None,
                                  normal_weights=None):
    """Meme probleme que calculate_deformation (memes residus, meme
    solveur), les noeuds epingles retires des inconnues.

    normal_weights = (w_out, w_tan) : contrainte de NORMALE FINALE
    atteignable par la machine actuelle (B d'un seul cote, C impose par la
    position ; parametre selectionnable, 2026-09-29). Deux residus de plus
    par tetraedre, lineaires en hauteurs deformees, mesures sur la geometrie
    d'ORIGINE (pinned[1]) : w_out * min(0, g.r) (la couche ne doit pas
    pencher vers le centre) et w_tan * (g.t) (pas de pente de cote).
    Borner le signe de la rotation par cellule ne suffisait pas : la
    deformation ajustee cree ses propres pentes (mesure sur dino et hors
    centre)."""
    node_mask, pinned_positions = np.asarray(pinned[0], dtype=bool), np.asarray(pinned[1], dtype=float)
    n_points = tet.number_of_points
    cells = tet.field_data["cells"]
    free = np.flatnonzero(~node_mask)
    column = np.full(n_points, -1)
    column[free] = np.arange(len(free))

    base_vertices = np.asarray(tet.points, dtype=float).copy()
    base_vertices[node_mask] = pinned_positions[node_mask]

    rotation_matrices = calculate_rotation_matrices(tet, rotation_field)
    old_vertices = tet.field_data["cell_vertices"][cells]
    old_vertices_transformed = np.einsum("ijk,ikl->ijl", rotation_matrices, (_N @ old_vertices).transpose(0, 2, 1))

    cell_rows = np.repeat(np.arange(tet.number_of_cells), cells.shape[1])
    vertex_ids = np.ravel(cells)
    keep = column[vertex_ids] >= 0
    jac_rows = np.tile(cell_rows[keep], 3)
    jac_cols = np.concatenate([column[vertex_ids[keep]] * 3 + dim for dim in range(3)])
    shape = (tet.number_of_cells, len(free) * 3)

    def vertices(x):
        v = base_vertices.copy()
        v[free] = x.reshape(-1, 3)
        return v

    def objective_function(x):
        transformed = (_N @ vertices(x)[cells]).transpose(0, 2, 1)
        return np.linalg.norm(transformed - old_vertices_transformed, axis=(1, 2)) ** 2

    def objective_jacobian(x):
        diff = ((_N @ vertices(x)[cells]).transpose(0, 2, 1) - old_vertices_transformed).transpose(0, 2, 1)
        data = np.concatenate([2 * diff[:, :, dim].ravel()[keep] for dim in range(3)]).astype(np.float32)
        return sp.coo_matrix((data, (jac_rows, jac_cols)), shape=shape).tocsr()

    sparsity = sp.coo_matrix((np.ones(len(jac_rows), dtype=np.int8), (jac_rows, jac_cols)), shape=shape).tocsr()

    if normal_weights is not None:
        w_out, w_tan = normal_weights
        a_r, a_t = normal_constraint_rows(pinned_positions, cells)
        n_cells = tet.number_of_cells
        z_cols = column[cells] * 3 + 2                                        # (T,4), -1*3+2 si epingle
        free4 = column[cells] >= 0
        nr_rows = np.repeat(np.arange(n_cells), 4).reshape(-1, 4)
        base_objective, base_jacobian = objective_function, objective_jacobian

        def objective_function(x):                                           # noqa: F811
            z4 = vertices(x)[cells][:, :, 2]
            sr = (a_r * z4).sum(axis=1)
            st = (a_t * z4).sum(axis=1)
            return np.concatenate([base_objective(x), w_out * np.minimum(0.0, sr), w_tan * st])

        def objective_jacobian(x):                                           # noqa: F811
            z4 = vertices(x)[cells][:, :, 2]
            active = ((a_r * z4).sum(axis=1) < 0)[:, None]
            rows = np.concatenate([nr_rows[free4] + n_cells, nr_rows[free4] + 2 * n_cells])
            cols = np.concatenate([z_cols[free4], z_cols[free4]])
            data = np.concatenate([(w_out * a_r * active)[free4], (w_tan * a_t)[free4]])
            extra = sp.coo_matrix((data, (rows - n_cells, cols)), shape=(2 * n_cells, shape[1]))
            return sp.vstack([base_jacobian(x), extra]).tocsr()

        rows = np.concatenate([nr_rows[free4], nr_rows[free4] + n_cells])
        cols = np.concatenate([z_cols[free4], z_cols[free4]])
        extra_s = sp.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(2 * n_cells, shape[1]))
        sparsity = sp.vstack([sparsity, extra_s]).tocsr()

    def _callback(intermediate_result):
        if progress_callback is not None:
            progress_callback("deformation", intermediate_result.nfev, ITERATIONS, float(intermediate_result.cost))

    result = least_squares(
        objective_function, base_vertices[free].ravel(), max_nfev=ITERATIONS, verbose=0,
        jac=objective_jacobian, jac_sparsity=sparsity, method="trf", x_scale="jac",
        callback=_callback if progress_callback is not None else None,
    )
    return vertices(result.x), result.cost


# ============================================================================
# Une iteration complete (optimize_rotations + calculate_deformation +
# update_tet_attributes), pour l'usage interactif du serveur.
# ============================================================================


def run_one_iteration(tet, state, params, progress_callback=None, planar_base=None):
    """params: dict avec les cles NEIGHBOUR_LOSS_WEIGHT, MAX_OVERHANG,
    ROTATION_MULTIPLIER, SET_INITIAL_ROTATION_TO_ZERO,
    INITIAL_ROTATION_FIELD_SMOOTHING, MAX_POS_ROTATION_DEG, MAX_NEG_ROTATION_DEG,
    ROTATE_ITERATIONS, DEFORM_ITERATIONS, STEEP_OVERHANG_COMPENSATION.

    progress_callback(phase, nfev, max_nfev, cost), phase in
    {"rotation", "deformation"} -- appele une fois par iteration EXTERNE de
    least_squares (pas par evaluation de fonction/jacobienne), donc cout
    negligeable meme sur un maillage dense.

    Retourne (nouveau_tet, info) ou info contient les couts et la fraction
    de surplomb residuelle -- jamais d'etat cache, tout est dans la valeur
    de retour.

    planar_base : None (comportement d'origine) ou (masque des cellules,
    masque des noeuds, positions d'ORIGINE des noeuds) issu de
    planar_base_mask -- la base est alors figee : rotation nulle, noeuds a
    leur position d'origine, a chaque passe."""
    cell_mask = planar_base[0] if planar_base is not None else None
    pinned = (planar_base[1], planar_base[2]) if planar_base is not None else None
    rotation_field, rotate_cost = optimize_rotations(
        tet, state,
        params["NEIGHBOUR_LOSS_WEIGHT"], params["MAX_OVERHANG"], params["ROTATION_MULTIPLIER"],
        params["ROTATE_ITERATIONS"], params["STEEP_OVERHANG_COMPENSATION"],
        params["INITIAL_ROTATION_FIELD_SMOOTHING"], params["SET_INITIAL_ROTATION_TO_ZERO"],
        np.deg2rad(params["MAX_POS_ROTATION_DEG"]), np.deg2rad(params["MAX_NEG_ROTATION_DEG"]),
        progress_callback=progress_callback, fixed_zero_cells=cell_mask,
    )

    if params.get("REACHABLE_NORMAL_WEIGHT"):
        if pinned is None:
            raise ValueError("REACHABLE_NORMAL_WEIGHT exige la base figee (positions d'origine)")
        w = float(params["REACHABLE_NORMAL_WEIGHT"])
        new_vertices, deform_cost = _calculate_deformation_pinned(
            tet, rotation_field, params["DEFORM_ITERATIONS"], pinned, progress_callback,
            normal_weights=(w, 0.5 * w))
    else:
        new_vertices, deform_cost = calculate_deformation(
            tet, rotation_field, params["DEFORM_ITERATIONS"], progress_callback=progress_callback,
            pinned=pinned,
        )

    deformed_tet = pv.UnstructuredGrid(tet.cells, np.full(tet.number_of_cells, pv.CellType.TETRA), new_vertices)
    for key in tet.field_data.keys():
        deformed_tet.field_data[key] = tet.field_data[key]
    for key in tet.cell_data.keys():
        deformed_tet.cell_data[key] = tet.cell_data[key]
    deformed_tet = update_tet_attributes(deformed_tet, state)

    overhang_fraction = float(
        np.nanmean(deformed_tet.cell_data["overhang_angle"] > np.deg2rad(90 + params["MAX_OVERHANG"]))
    )

    info = {
        "rotate_cost": float(rotate_cost),
        "deform_cost": float(deform_cost),
        "overhang_fraction": overhang_fraction,
        "number_of_cells": int(tet.number_of_cells),
    }
    return deformed_tet, info


def export_surface_stl(tet, path):
    """Recentre en bas de bounding box (comme la cellule 11 du notebook) et
    exporte la surface. Modifie tet.points en place (comportement identique
    au notebook original)."""
    x_min, x_max, y_min, y_max, z_min, z_max = tet.bounds
    offsets = np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min])
    tet.points -= offsets
    tet.extract_surface().save(path)
