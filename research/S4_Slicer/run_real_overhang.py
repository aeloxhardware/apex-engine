"""
Extraction fidele des cellules 2, 3, 4, 7, 11 du notebook main.ipynb (S4-Slicer,
Joshua Bird, GPL-3.0) -- executee sur res/test_models/real_overhang.stl du
projet Apex, avec les parametres PAR DEFAUT de la cellule 4 (une seule passe,
pas d'iteration manuelle).

SEUL changement par rapport au notebook original : open3d (o3d.io.read_triangle_mesh)
remplace par trimesh.load(...), car open3d n'a pas de wheel pour Python 3.14.
Les tableaux .vertices/.faces produits sont numeriquement equivalents a
mesh.vertices/mesh.triangles d'open3d -- aucune autre ligne de l'algorithme
de deformation n'est modifiee.

SAVE_GIF=False partout : on ne veut que le maillage deforme en sortie, pas
d'animation.
"""
import sys
import time
import base64
import pickle

import networkx as nx
import numpy as np
import pyvista as pv
import tetgen
import trimesh
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation as R
from scipy.sparse import lil_matrix

up_vector = np.array([0, 0, 1])


def encode_object(obj):
    return base64.b64encode(pickle.dumps(obj)).decode("utf-8")


def decode_object(encoded_str):
    return pickle.loads(base64.b64decode(encoded_str))


# ============================================================================
# CELLULE 2 (adaptee : trimesh au lieu d'open3d, model_name = real_overhang)
# ============================================================================

model_name = sys.argv[1] if len(sys.argv) > 1 else "real_overhang"
mesh = trimesh.load(f"input_models/{model_name}.stl")

input_tet = tetgen.TetGen(np.asarray(mesh.vertices), np.asarray(mesh.faces))
input_tet.tetrahedralize()
input_tet = input_tet.grid

PART_OFFSET = np.array([0.0, 0.0, 0.0])
x_min, x_max, y_min, y_max, z_min, z_max = input_tet.bounds
input_tet.points -= np.array(
    [(x_min + x_max) / 2, (y_min + y_max) / 2, z_min]
) + PART_OFFSET

print(f"[real_overhang] tetraedres: {input_tet.number_of_cells}, "
      f"points: {input_tet.number_of_points}")

# find neighbours
cell_neighbour_dict = {
    neighbour_type: {face: [] for face in range(input_tet.number_of_cells)}
    for neighbour_type in ["point", "edge", "face"]
}
for neighbour_type in ["point", "edge", "face"]:
    cell_neighbours = []
    for cell_index in range(input_tet.number_of_cells):
        neighbours = input_tet.cell_neighbors(cell_index, f"{neighbour_type}s")
        for neighbour in neighbours:
            if neighbour > cell_index:
                cell_neighbours.append((cell_index, neighbour))
    for face_1, face_2 in np.array(cell_neighbours):
        cell_neighbour_dict[neighbour_type][face_1].append(face_2)
        cell_neighbour_dict[neighbour_type][face_2].append(face_1)

    input_tet.field_data[f"cell_{neighbour_type}_neighbours"] = np.array(cell_neighbours)

cell_neighbour_graph = nx.Graph()
cell_centers = input_tet.cell_centers().points
for edge in input_tet.field_data["cell_point_neighbours"]:
    distance = np.linalg.norm(cell_centers[edge[0]] - cell_centers[edge[1]])
    cell_neighbour_graph.add_weighted_edges_from([(edge[0], edge[1], distance)])


def update_tet_attributes(tet):
    surface_mesh = tet.extract_surface()
    cell_to_face = decode_object(tet.field_data["cell_to_face"])

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
    bottom_cells_mask = tet.cell_data["face_center"][:, 2] < bottom_cell_threshold
    tet.cell_data["is_bottom"] = bottom_cells_mask
    bottom_cells = np.where(bottom_cells_mask)[0]

    face_normals = tet.cell_data["face_normal"].copy()
    face_normals[bottom_cells_mask] = np.nan
    overhang_angle = np.arccos(np.dot(face_normals, up_vector))
    tet.cell_data["overhang_angle"] = overhang_angle

    overhang_direction = face_normals[:, :2].copy()
    overhang_direction /= np.linalg.norm(overhang_direction, axis=1)[:, None]
    tet.cell_data["overhang_direction"] = overhang_direction

    IN_AIR_THRESHOLD = 1
    tet.cell_data["in_air"] = np.full(tet.number_of_cells, False)

    _, paths_to_bottom = nx.multi_source_dijkstra(cell_neighbour_graph, set(bottom_cells))

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
    surface_mesh = tet.extract_surface()

    cells = tet.cells.reshape(-1, 5)[:, 1:]
    tet.add_field_data(cells, "cells")
    cell_vertices = tet.points
    tet.add_field_data(cell_vertices, "cell_vertices")
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]
    tet.add_field_data(faces, "faces")
    face_vertices = surface_mesh.points
    tet.add_field_data(face_vertices, "face_vertices")

    shared_vertices = []
    for cell_1, cell_2 in tet.field_data["cell_point_neighbours"]:
        shared_vertices_these_faces = np.intersect1d(cells[cell_1], cells[cell_2])
        for vertex in shared_vertices_these_faces:
            shared_vertices.append(
                {
                    "cell_1_index": cell_1,
                    "cell_2_index": cell_2,
                    "cell_1_vertex_index": np.where(cells[cell_1] == vertex)[0][0],
                    "cell_2_vertex_index": np.where(cells[cell_2] == vertex)[0][0],
                }
            )

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

    tet = update_tet_attributes(tet)

    bottom_cells_mask = tet.cell_data["is_bottom"]
    bottom_cells = np.where(bottom_cells_mask)[0]

    tet.cell_data["overhang_angle"][bottom_cells] = np.nan

    return tet, bottom_cells_mask, bottom_cells


bottom_cells_mask = None
bottom_cells = None
input_tet, bottom_cells_mask, bottom_cells = calculate_tet_attributes(input_tet)

bottom_cell_graph = nx.Graph()
for cell_index in bottom_cells:
    bottom_cell_graph.add_node(cell_index)
cell_point_neighbour_dict = cell_neighbour_dict["point"]
for cell_index in bottom_cells:
    for neighbour in cell_point_neighbour_dict[cell_index]:
        if neighbour in bottom_cells:
            bottom_cell_graph.add_edge(cell_index, neighbour)

bottom_cell_groups = [list(x) for x in list(nx.connected_components(bottom_cell_graph))]

undeformed_tet = input_tet.copy()

print(f"[real_overhang] cellules au contact plateau: {int(bottom_cells_mask.sum())}")

# ============================================================================
# CELLULE 3 (inchangee)
# ============================================================================


def planeFit(points):
    import numpy as np
    from numpy.linalg import svd

    points = np.reshape(points, (np.shape(points)[0], -1))
    assert points.shape[0] <= points.shape[1], (
        f"There are only {points.shape[1]} points in {points.shape[0]} dimensions."
    )
    ctr = points.mean(axis=1)
    x = points - ctr[:, np.newaxis]
    M = np.dot(x, x.T)
    return ctr, svd(M)[0][:, -1]


def calculate_path_length_to_base_gradient(tet, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO):
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
        for i in range(INITIAL_ROTATION_FIELD_SMOOTHING):
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
# CELLULE 4 (inchangee, SAVE_GIF=False, MAX_POS/NEG_ROTATION corriges a 360
# comme le dit le commentaire "normally set to 360 unless you get collisions")
# ============================================================================


def calculate_initial_rotation_field(
    tet, MAX_OVERHANG, ROTATION_MULTIPLIER, STEEP_OVERHANG_COMPENSATION,
    INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION,
):
    initial_rotation_field = np.full((tet.number_of_cells), np.nan)
    initial_rotation_field = np.abs(np.deg2rad(90 + MAX_OVERHANG) - tet.cell_data["overhang_angle"])

    path_length_to_base_gradient = calculate_path_length_to_base_gradient(
        tet, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO
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
    # FIX numpy>=2.0 : np.cross exige desormais des vecteurs 3D des deux
    # cotes (le broadcast implicite d'un vecteur 2D -> 3D, z=0, a ete
    # retire). Complete explicitement par une 3e composante nulle,
    # mathematiquement identique a ce que l'ancien comportement produisait.
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


def apply_rotation_field(tet, rotation_field):
    new_vertices = np.zeros((tet.number_of_points, 3))
    vertices_count = np.zeros((tet.number_of_points))
    for cell in tet.field_data["cells"]:
        vertices_count[cell] += 1

    unique_vertices_rotated = calculate_unique_vertices_rotated(tet, rotation_field)

    for cell_index, vertices in enumerate(unique_vertices_rotated):
        for i, vertex in enumerate(vertices):
            new_vertices[tet.field_data["cells"][cell_index, i]] += vertex.T[0] / vertices_count[
                tet.field_data["cells"][cell_index][i]
            ]

    new_tet = pv.UnstructuredGrid(tet.cells, np.full(tet.number_of_cells, pv.CellType.TETRA), new_vertices)

    return new_tet


def optimize_rotations(
    tet, NEIGHBOUR_LOSS_WEIGHT, MAX_OVERHANG, ROTATION_MULTIPLIER, ITERATIONS, SAVE_GIF,
    STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO,
    MAX_POS_ROTATION, MAX_NEG_ROTATION,
):
    initial_rotation_field = calculate_initial_rotation_field(
        tet, MAX_OVERHANG, ROTATION_MULTIPLIER, STEEP_OVERHANG_COMPENSATION,
        INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION,
    )
    num_cells_with_initial_rotation = np.sum(~np.isnan(initial_rotation_field))

    def objective_function(rotation_field):
        cell_face_neighbours = tet.field_data["cell_face_neighbours"]
        neighbour_differences = rotation_field[cell_face_neighbours[:, 0]] - rotation_field[cell_face_neighbours[:, 1]]
        neighbour_losses = NEIGHBOUR_LOSS_WEIGHT * neighbour_differences**2

        overhanging_mask = tet.cell_data["overhang_angle"] > np.deg2rad(90 + MAX_OVERHANG)
        valid_cell_indices = np.where(~np.isnan(initial_rotation_field))[0]
        initial_rotation_losses = (rotation_field[valid_cell_indices] - initial_rotation_field[valid_cell_indices]) ** 2

        return np.concatenate((neighbour_losses, initial_rotation_losses))

    def objective_jacobian(rotation_field):
        cell_face_neighbours = tet.field_data["cell_face_neighbours"]
        jac = lil_matrix((len(cell_face_neighbours) + num_cells_with_initial_rotation, tet.number_of_cells), dtype=np.float32)

        cell_1 = cell_face_neighbours[:, 0]
        cell_2 = cell_face_neighbours[:, 1]

        differences = rotation_field[cell_1] - rotation_field[cell_2]

        jac[range(len(cell_face_neighbours)), cell_1] = 2 * NEIGHBOUR_LOSS_WEIGHT * differences
        jac[range(len(cell_face_neighbours)), cell_2] = -2 * NEIGHBOUR_LOSS_WEIGHT * differences

        valid_cell_indices = np.where(~np.isnan(initial_rotation_field))[0]

        jac[len(cell_face_neighbours) + np.arange(len(valid_cell_indices)), valid_cell_indices] = (
            2 * (rotation_field[valid_cell_indices] - initial_rotation_field[valid_cell_indices])
        )

        return jac.tocsr()

    def jac_sparsity():
        cell_face_neighbours = tet.field_data["cell_face_neighbours"]
        sparsity = lil_matrix((len(cell_face_neighbours) + num_cells_with_initial_rotation, tet.number_of_cells), dtype=np.int8)

        for i, (cell_1, cell_2) in enumerate(cell_face_neighbours):
            sparsity[i, cell_1] = 1
            sparsity[i, cell_2] = 1

        valid_cell_indices = np.where(~np.isnan(initial_rotation_field))[0]
        i = 0
        for cell_index, initial_rotation in enumerate(initial_rotation_field):
            if cell_index in valid_cell_indices:
                sparsity[len(cell_face_neighbours) + i, cell_index] = 1
                i += 1

        return sparsity.tocsr()

    smoothed_rotation_field = np.zeros((tet.number_of_cells))

    result = least_squares(
        objective_function,
        smoothed_rotation_field,
        jac=objective_jacobian,
        max_nfev=ITERATIONS,
        jac_sparsity=jac_sparsity(),
        verbose=2,
        method="trf",
        ftol=1e-6,
    )

    return result.x


NEIGHBOUR_LOSS_WEIGHT = 20
MAX_OVERHANG = 30
ROTATION_MULTIPLIER = 1
SET_INITIAL_ROTATION_TO_ZERO = False
INITIAL_ROTATION_FIELD_SMOOTHING = 30
MAX_POS_ROTATION = np.deg2rad(360)   # corrige : le commentaire du notebook dit
MAX_NEG_ROTATION = np.deg2rad(-360)  # "normally set to 360 unless you get collisions"
ITERATIONS = 100
SAVE_GIF = False
STEEP_OVERHANG_COMPENSATION = True

t0 = time.time()
rotation_field = optimize_rotations(
    undeformed_tet, NEIGHBOUR_LOSS_WEIGHT, MAX_OVERHANG, ROTATION_MULTIPLIER, ITERATIONS, SAVE_GIF,
    STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO,
    MAX_POS_ROTATION, MAX_NEG_ROTATION,
)
print(f"[real_overhang] optimize_rotations: {time.time() - t0:.2f}s")

undeformed_tet_with_rotated_tetrahedrons = apply_rotation_field_unique_vertices(undeformed_tet, rotation_field)
undeformed_tet_with_rotated_tetrahedrons.cell_data["rotation_field"] = rotation_field

# ============================================================================
# CELLULE 7 (inchangee, SAVE_GIF=False)
# ============================================================================

N = np.eye(4) - 1 / 4 * np.ones((4, 4))


def calculate_deformation(tet, rotation_field, ITERATIONS, SAVE_GIF):
    new_vertices = tet.points.copy()

    params = new_vertices.flatten()

    rotation_matrices = calculate_rotation_matrices(tet, rotation_field)

    old_vertices = tet.field_data["cell_vertices"][tet.field_data["cells"]]
    old_vertices_transformed = np.einsum("ijk,ikl->ijl", rotation_matrices, (N @ old_vertices).transpose(0, 2, 1))

    def objective_function(params):
        new_vertices = params[: tet.number_of_points * 3].reshape(-1, 3)

        new_vertices_transformed = (N @ new_vertices[tet.field_data["cells"]]).transpose(0, 2, 1)

        position_losses = np.linalg.norm(new_vertices_transformed - old_vertices_transformed, axis=(1, 2)) ** 2

        return position_losses

    def objective_jacobian(params):
        J = lil_matrix((tet.number_of_cells, len(params)), dtype=np.float32)

        new_vertices = params[: tet.number_of_points * 3].reshape(-1, 3)

        old_vertices = tet.field_data["cell_vertices"][tet.field_data["cells"]]

        new_vertices_transformed = (N @ new_vertices[tet.field_data["cells"]]).transpose(0, 2, 1)

        diff = new_vertices_transformed - old_vertices_transformed

        diff = diff.transpose(0, 2, 1)

        cell_indices = np.repeat(np.arange(tet.number_of_cells), len(tet.field_data["cells"][0]))
        vertex_indices = np.ravel(tet.field_data["cells"])

        for dim in range(3):
            J[cell_indices, vertex_indices * 3 + dim] = 2 * diff[:, :, dim].ravel()

        return J.tocsr()

    def jac_sparsity():
        sparsity = lil_matrix((tet.number_of_cells, len(params)), dtype=np.int8)

        cell_indices = np.repeat(np.arange(tet.number_of_cells), len(tet.field_data["cells"][0]))
        vertex_indices = np.ravel(tet.field_data["cells"])

        for dim in range(3):
            sparsity[cell_indices, vertex_indices * 3 + dim] = 1

        return sparsity.tocsr()

    result = least_squares(
        objective_function,
        params,
        max_nfev=ITERATIONS,
        verbose=2,
        jac=objective_jacobian,
        jac_sparsity=jac_sparsity(),
        method="trf",
        x_scale="jac",
    )

    return result.x[: tet.number_of_points * 3].reshape(-1, 3)


ITERATIONS = 1000
SAVE_GIF = False
t0 = time.time()
new_vertices = calculate_deformation(undeformed_tet, rotation_field, ITERATIONS, SAVE_GIF)
print(f"[real_overhang] calculate_deformation: {time.time() - t0:.2f}s")
deformed_tet = pv.UnstructuredGrid(undeformed_tet.cells, np.full(undeformed_tet.number_of_cells, pv.CellType.TETRA), new_vertices)

for key in undeformed_tet.field_data.keys():
    deformed_tet.field_data[key] = undeformed_tet.field_data[key]
for key in undeformed_tet.cell_data.keys():
    deformed_tet.cell_data[key] = undeformed_tet.cell_data[key]
deformed_tet = update_tet_attributes(deformed_tet)

# ============================================================================
# CELLULE 11 (inchangee) : export STL du mesh deforme
# ============================================================================

x_min, x_max, y_min, y_max, z_min, z_max = deformed_tet.bounds
offsets_applied = np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min])
deformed_tet.points -= offsets_applied

import os
os.makedirs("output_models", exist_ok=True)
deformed_tet.extract_surface().save(f"output_models/{model_name}_deformed_tet.stl")

print(f"[real_overhang] mesh deforme exporte -> output_models/{model_name}_deformed_tet.stl")
print(f"[real_overhang] fraction de cellules encore en surplomb (>{90+MAX_OVERHANG} deg): "
      f"{np.nanmean(deformed_tet.cell_data['overhang_angle'] > np.deg2rad(90 + MAX_OVERHANG)):.4f}")

# Sauvegarde pour reprise manuelle (equivalent cellule 12 du notebook) :
# permet a iterate.py de faire "undeformed_tet = deformed_tet.copy()" sans
# tout relancer depuis l'ingestion STL.
os.makedirs("pickle_files", exist_ok=True)
with open(f"pickle_files/deformed_{model_name}.pkl", "wb") as f:
    pickle.dump(deformed_tet, f)
print(f"[real_overhang] etat sauvegarde -> pickle_files/deformed_{model_name}.pkl (iteration 1)")
