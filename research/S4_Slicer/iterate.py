"""
Une iteration manuelle supplementaire de la methode Bird -- equivalent des
cellules markdown 8/9 du notebook original ("Run below to do another
iteration" -> undeformed_tet = deformed_tet.copy()).

Charge l'etat sauvegarde par run_real_overhang.py (ou une iteration
precedente de ce script), relance optimize_rotations + calculate_deformation
sur ce nouveau point de depart, exporte le STL et re-sauvegarde l'etat.

Usage : python iterate.py [NEIGHBOUR_LOSS_WEIGHT] [MAX_OVERHANG] [ROTATION_MULTIPLIER] [DEFORM_ITERATIONS]
Sans argument : reprend les memes parametres que la passe precedente
(NEIGHBOUR_LOSS_WEIGHT=20, MAX_OVERHANG=30, ROTATION_MULTIPLIER=2, DEFORM_ITERATIONS=1000).
"""
import sys
import os
import time
import base64
import pickle

import networkx as nx
import numpy as np
import pyvista as pv
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation as R
from scipy.sparse import lil_matrix

model_name = "real_overhang"
up_vector = np.array([0, 0, 1])


def encode_object(obj):
    return base64.b64encode(pickle.dumps(obj)).decode("utf-8")


def decode_object(encoded_str):
    return pickle.loads(base64.b64decode(encoded_str))


# ---------------------------------------------------------------------------
# Charger l'etat (deformed_tet de la passe precedente devient undeformed_tet)
# ---------------------------------------------------------------------------

with open(f"pickle_files/deformed_{model_name}.pkl", "rb") as f:
    undeformed_tet = pickle.load(f)

# Reconstruire ce dont optimize_rotations/calculate_deformation ont besoin
# et qui n'est pas persiste avec le tet lui-meme (graphes, dicts de voisinage).
cell_neighbour_dict = {
    neighbour_type: {face: [] for face in range(undeformed_tet.number_of_cells)}
    for neighbour_type in ["point", "edge", "face"]
}
for neighbour_type in ["point", "edge", "face"]:
    for cell_1, cell_2 in undeformed_tet.field_data[f"cell_{neighbour_type}_neighbours"]:
        cell_neighbour_dict[neighbour_type][cell_1].append(cell_2)
        cell_neighbour_dict[neighbour_type][cell_2].append(cell_1)

cell_neighbour_graph = nx.Graph()
cell_centers = undeformed_tet.cell_centers().points
for edge in undeformed_tet.field_data["cell_point_neighbours"]:
    distance = np.linalg.norm(cell_centers[edge[0]] - cell_centers[edge[1]])
    cell_neighbour_graph.add_weighted_edges_from([(edge[0], edge[1], distance)])

bottom_cells_mask = undeformed_tet.cell_data["is_bottom"]
bottom_cells = np.where(bottom_cells_mask)[0]

# ---------------------------------------------------------------------------
# CELLULE 3 (planeFit, calculate_path_length_to_base_gradient) -- inchangee
# ---------------------------------------------------------------------------


def planeFit(points):
    from numpy.linalg import svd

    points = np.reshape(points, (np.shape(points)[0], -1))
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


# ---------------------------------------------------------------------------
# CELLULE 4 -- inchangee (fix numpy>=2.0 sur np.cross deja applique)
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# CELLULE 7 -- inchangee
# ---------------------------------------------------------------------------

N = np.eye(4) - 1 / 4 * np.ones((4, 4))


def calculate_deformation(tet, rotation_field, ITERATIONS):
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
    bottom_cells_mask_local = tet.cell_data["face_center"][:, 2] < bottom_cell_threshold
    tet.cell_data["is_bottom"] = bottom_cells_mask_local
    bottom_cells_local = np.where(bottom_cells_mask_local)[0]

    face_normals = tet.cell_data["face_normal"].copy()
    face_normals[bottom_cells_mask_local] = np.nan
    overhang_angle = np.arccos(np.dot(face_normals, up_vector))
    tet.cell_data["overhang_angle"] = overhang_angle

    overhang_direction = face_normals[:, :2].copy()
    overhang_direction /= np.linalg.norm(overhang_direction, axis=1)[:, None]
    tet.cell_data["overhang_direction"] = overhang_direction

    IN_AIR_THRESHOLD = 1
    tet.cell_data["in_air"] = np.full(tet.number_of_cells, False)

    _, paths_to_bottom = nx.multi_source_dijkstra(cell_neighbour_graph, set(bottom_cells_local))

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


# ---------------------------------------------------------------------------
# Parametres de cette iteration (identiques par defaut a l'iteration precedente)
# ---------------------------------------------------------------------------

NEIGHBOUR_LOSS_WEIGHT = float(sys.argv[1]) if len(sys.argv) > 1 else 20
MAX_OVERHANG = float(sys.argv[2]) if len(sys.argv) > 2 else 30
ROTATION_MULTIPLIER = float(sys.argv[3]) if len(sys.argv) > 3 else 2
DEFORM_ITERATIONS = int(sys.argv[4]) if len(sys.argv) > 4 else 1000

SET_INITIAL_ROTATION_TO_ZERO = False
INITIAL_ROTATION_FIELD_SMOOTHING = 30
MAX_POS_ROTATION = np.deg2rad(360)
MAX_NEG_ROTATION = np.deg2rad(-360)
ROTATE_ITERATIONS = 100
STEEP_OVERHANG_COMPENSATION = True

print(f"[iterate] parametres: NEIGHBOUR_LOSS_WEIGHT={NEIGHBOUR_LOSS_WEIGHT}, "
      f"MAX_OVERHANG={MAX_OVERHANG}, ROTATION_MULTIPLIER={ROTATION_MULTIPLIER}, "
      f"DEFORM_ITERATIONS={DEFORM_ITERATIONS}")

t0 = time.time()
rotation_field = optimize_rotations(
    undeformed_tet, NEIGHBOUR_LOSS_WEIGHT, MAX_OVERHANG, ROTATION_MULTIPLIER, ROTATE_ITERATIONS, False,
    STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO,
    MAX_POS_ROTATION, MAX_NEG_ROTATION,
)
print(f"[iterate] optimize_rotations: {time.time() - t0:.2f}s")

t0 = time.time()
new_vertices = calculate_deformation(undeformed_tet, rotation_field, DEFORM_ITERATIONS)
print(f"[iterate] calculate_deformation: {time.time() - t0:.2f}s")

deformed_tet = pv.UnstructuredGrid(undeformed_tet.cells, np.full(undeformed_tet.number_of_cells, pv.CellType.TETRA), new_vertices)

for key in undeformed_tet.field_data.keys():
    deformed_tet.field_data[key] = undeformed_tet.field_data[key]
for key in undeformed_tet.cell_data.keys():
    deformed_tet.cell_data[key] = undeformed_tet.cell_data[key]
deformed_tet = update_tet_attributes(deformed_tet)

# ---------------------------------------------------------------------------
# Export + sauvegarde d'etat pour la prochaine iteration
# ---------------------------------------------------------------------------

x_min, x_max, y_min, y_max, z_min, z_max = deformed_tet.bounds
offsets_applied = np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min])
deformed_tet.points -= offsets_applied

os.makedirs("output_models", exist_ok=True)
deformed_tet.extract_surface().save(f"output_models/{model_name}_deformed_tet.stl")

with open(f"pickle_files/deformed_{model_name}.pkl", "wb") as f:
    pickle.dump(deformed_tet, f)

overhang_fraction = np.nanmean(deformed_tet.cell_data['overhang_angle'] > np.deg2rad(90 + MAX_OVERHANG))
print(f"[iterate] mesh exporte -> output_models/{model_name}_deformed_tet.stl")
print(f"[iterate] fraction en surplomb (>{90+MAX_OVERHANG} deg): {overhang_fraction:.4f}")
