import json
import sys

# Liste des modules à importer
imports_to_process = [
    {"task": "IMPORT_NETWORKX", "code": "import networkx as nx"},
    {"task": "IMPORT_NUMPY", "code": "import numpy as np"},
    {"task": "IMPORT_PYVISTA", "code": "import pyvista as pv"},
    {"task": "IMPORT_TETGEN", "code": "import tetgen"},
    {"task": "IMPORT_SCIPY_OPTIMIZE", "code": "from scipy.optimize import minimize, least_squares"},
    {"task": "IMPORT_SCIPY_ROTATION", "code": "from scipy.spatial.transform import Rotation as R"},
    {"task": "IMPORT_OPEN3D", "code": "import open3d as o3d"},
    {"task": "IMPORT_TIME", "code": "import time"},
    {"task": "IMPORT_OS", "code": "import os"},
    {"task": "IMPORT_SYS", "code": "import sys"},
    {"task": "IMPORT_JSON", "code": "import json"},
    {"task": "IMPORT_PICKLE", "code": "import pickle"},
    {"task": "IMPORT_BASE64", "code": "import base64"},
    {"task": "IMPORT_CUPY", "code": "import cupy as cp"},
    {"task": "IMPORT_CUPYX", "code": "import cupyx"},
    {"task": "IMPORT_CUPYX_SCIPY_SPARSE", "code": "import cupyx.scipy.sparse as cusparse"},
    {"task": "IMPORT_CUPYX_SCIPY_SPARSE_LINALG", "code": "import cupyx.scipy.sparse.linalg as culinalg"},
    {"task": "IMPORT_SCIPY_SPARSE_LIL", "code": "from scipy.sparse import lil_matrix"},
    {"task": "IMPORT_PYGCODE", "code": "from pygcode import Line"},
    {"task": "IMPORT_CUPYX_SCATTER_ADD", "code": "from cupyx import scatter_add"},
    {"task": "IMPORT_CONFIG", "code": "import config; from config import *"}
]

total_imports = len(imports_to_process)

# Exécution avec une progression calculée de 0 à 100
for i, item in enumerate(imports_to_process):
    # Calcul du pourcentage (de 0 à 100)
    percentage = int((i / total_imports) * 100)
    
    print(json.dumps({
        "type": "progress", 
        "task": item["task"], 
        "current": percentage, 
        "total": 100
    }))
    sys.stdout.flush()
    
    # Exécution de l'import
    exec(item["code"])

# Envoi du statut final à 100%
print(json.dumps({
    "type": "progress", 
    "task": "IMPORTS_COMPLETED", 
    "current": 100, 
    "total": 100
}))
sys.stdout.flush()

#Variables globales
input_tet = None
up_vector = np.array([0, 0, 1])
cell_neighbour_graph = None
cell_neighbour_dict = None
undeformed_tet = None
deformed_tet = None
bottom_cells = None
bottom_cells_mask = None
rotation_field = None
initial_rotation_field = None
N = None
gcode_points = None
vertex_transformations = None
vertex_rotations = None
cells_array = None
precalc_new_positions = None
precalc_rotations = None
z_squish_scales = None
new_gcode_points = None
temp = None
new_vertices = None
tet = None
up_vector_gpu = None

def encode_object(obj):
    return base64.b64encode(pickle.dumps(obj)).decode('utf-8')

def decode_object(encoded_str):
    return pickle.loads(base64.b64decode(encoded_str))
    
def SETUP():
    global input_tet, cell_neighbour_graph, cell_neighbour_dict, up_vector_gpu

    up_vector_np = np.array([0, 0, 1])
    up_vector_gpu = cp.array([0, 0, 1], dtype=cp.float64)
    
    print("[GLOBAL] Succès du démarrage de l'exécution. Si ce message venait à ne pas s'afficher après 10 secondes, merci de recommencer.")

    pv.set_jupyter_backend('trame')
    pv.global_theme.notebook = False


    # -----------------------------------------------------------------------
    # Utilitaires (inchangés — pickle/base64 sont CPU par nature, encoder/
    # décoder un objet Python arbitraire sur GPU n'a pas de sens : ce sont des
    # dictionnaires Python, pas des tenseurs numériques)
    # -----------------------------------------------------------------------

    def encode_object(obj):
        return base64.b64encode(pickle.dumps(obj)).decode('utf-8')


    def decode_object(encoded_str):
        return pickle.loads(base64.b64decode(encoded_str))


    up_vector_np = np.array([0, 0, 1])
    up_vector_gpu = cp.array([0, 0, 1], dtype=cp.float64)

    print(f"[GPU] CuPy détecté — device actif : {cp.cuda.Device()}")
    print(f"[GPU] Mémoire GPU libre / totale : {[x / 1e9 for x in cp.cuda.Device().mem_info]} Go")


    # -----------------------------------------------------------------------
    # Chargement du maillage et tétraédrisation (reste CPU — voir note en tête
    # de fichier : TetGen n'a pas d'équivalent GPU)
    # -----------------------------------------------------------------------

    model_name = "model"
    mesh = o3d.io.read_triangle_mesh(f'input_model/{model_name}.stl')
    print(f"[CPU] Maillage chargé : {len(mesh.vertices)} sommets, {len(mesh.triangles)} triangles")

    input_tet = tetgen.TetGen(np.asarray(mesh.vertices), np.asarray(mesh.triangles))
    input_tet.tetrahedralize()
    input_tet = input_tet.grid
    print(f"[CPU] Tétraédrisation terminée : {input_tet.number_of_cells} cellules")

    PART_OFFSET = np.array([0., 0., 0.])
    x_min, x_max, y_min, y_max, z_min, z_max = input_tet.bounds
    input_tet.points -= np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min]) + PART_OFFSET


    # -----------------------------------------------------------------------
    # Recherche de voisins (point/edge/face) et construction du graphe de
    # proximité pondéré.
    #
    # Portage GPU : la boucle Python `for cell_index in range(...)` appelant
    # `input_tet.cell_neighbors(...)` reste sur CPU car c'est une requête de
    # connectivité topologique VTK (pas un calcul numérique, pas vectorisable
    # par nature — VTK n'expose pas cette API en batch). En revanche, le
    # calcul des DISTANCES entre centres de cellules voisines (utilisé pour
    # pondérer le graphe de proximité) est un calcul numérique pur sur des
    # tableaux de taille O(n_edges) : c'est vectorisé et déplacé sur GPU.
    # -----------------------------------------------------------------------

    cell_neighbour_dict = {neighbour_type: {face: [] for face in range(input_tet.number_of_cells)} for neighbour_type in ["point", "edge", "face"]}
    
    neighbour_types = ["point", "edge", "face"]
    total_types = len(neighbour_types)

    for i, neighbour_type in enumerate(neighbour_types):
        # Émission du début de tâche pour cette catégorie
        progress_data = {
            "type": "progress",
            "task": f"search_neighbours",
            "current": i,
            "total": total_types
        }
        sys.stdout.write(json.dumps(progress_data) + "\n")
        sys.stdout.flush()

        cell_neighbours = []
        for cell_index in range(input_tet.number_of_cells):
            neighbours = input_tet.cell_neighbors(cell_index, f"{neighbour_type}s")
            for neighbour in neighbours:
                if neighbour > cell_index:
                    cell_neighbours.append((cell_index, neighbour))
            
            # Optionnel : Émettre une mise à jour tous les X cellules pour éviter la saturation
            if cell_index % 1 == 0:
                progress_update = {
                    "type": "progress",
                    "task": f"search_neighbours_{neighbour_type}",
                    "current": cell_index,
                    "total": input_tet.number_of_cells
                }
                sys.stdout.write(json.dumps(progress_update) + "\n")
                sys.stdout.flush()

        for face_1, face_2 in np.array(cell_neighbours):
            cell_neighbour_dict[neighbour_type][face_1].append(face_2)
            cell_neighbour_dict[neighbour_type][face_2].append(face_1)

        input_tet.field_data[f"cell_{neighbour_type}_neighbours"] = np.array(cell_neighbours)

    print(f"[GPU] Transfert des centres de cellules et des paires de voisins vers le GPU")
    cell_centers_cpu = input_tet.cell_centers().points
    cell_centers_gpu = cp.asarray(cell_centers_cpu)  # <-- unique transfert CPU->GPU pour cette section

    point_edges_gpu = cp.asarray(input_tet.field_data["cell_point_neighbours"])
    # indices des deux extrémités de chaque arête, vectorisé
    edge_src = point_edges_gpu[:, 0]
    edge_dst = point_edges_gpu[:, 1]

    # calcul vectorisé de la distance euclidienne pour TOUTES les arêtes en une
    # seule opération GPU (remplace la boucle Python + nx.add_weighted_edges_from
    # implicite ligne par ligne de l'original)
    diff_gpu = cell_centers_gpu[edge_src] - cell_centers_gpu[edge_dst]
    distances_gpu = cp.linalg.norm(diff_gpu, axis=1)
    print(f"[GPU] {distances_gpu.shape[0]} distances de voisinage calculées sur GPU")

    # Le graphe NetworkX lui-même doit rester un objet Python CPU : c'est une
    # structure de données CPU-native utilisée ensuite pour Dijkstra "de
    # référence" si besoin de debug, mais la vraie traversée se fera en GPU
    # (voir gpu_multi_source_dijkstra plus bas). On construit quand même le
    # graphe NetworkX pour compatibilité totale de l'API / du comportement en
    # cas d'inspection manuelle, mais AUCUN calcul numérique lourd n'est fait
    # via NetworkX dans la suite de ce fichier.
    distances_cpu = cp.asnumpy(distances_gpu)  # transfert de sortie unique pour cette section
    edge_src_cpu = cp.asnumpy(edge_src)
    edge_dst_cpu = cp.asnumpy(edge_dst)

    cell_neighbour_graph = nx.Graph()
    weighted_edges = list(zip(edge_src_cpu.tolist(), edge_dst_cpu.tolist(), distances_cpu.tolist()))
    cell_neighbour_graph.add_weighted_edges_from(weighted_edges)
    print(f"[CPU] Graphe de proximité construit : {cell_neighbour_graph.number_of_nodes()} noeuds, {cell_neighbour_graph.number_of_edges()} arêtes")


    # -----------------------------------------------------------------------
    # update_tet_attributes : normales de face, centres de face, centres de
    # cellule, bottom cells, overhang_angle/direction, in_air.
    #
    # Portage GPU : `tet.extract_surface(...)` reste CPU (VTK, topologie pure,
    # cf. note en tête de fichier). Tout le calcul numérique qui en découle
    # (indexation des normales/centres par cell_to_face, produit scalaire avec
    # up_vector, normalisation, comparaison de hauteurs pour in_air) est
    # vectorisé et déplacé sur GPU. Les boucles Python
    # `for cell_index, face_indices in cell_to_face.items()` de l'original
    # étaient nécessaires seulement à cause de la taille variable de
    # face_indices par cellule (ragged array) — remplacées ici par un mapping
    # GPU vectorisé (un indice de face "most-down" par cellule), qui produit
    # exactement le même résultat sans boucle Python ni triple aller-retour
    # CPU/GPU.
    # -----------------------------------------------------------------------
        
def update_tet_attributes(tet, up_vector_gpu, cell_neighbour_graph):

    print("DEBUG: Entrée dans update_tet_attributes...")
    
    surface_mesh = tet.extract_surface(algorithm=None)
    
    # 1. EXTRACTION ET DÉCODAGE UNIQUE ET ROBUSTE
    cell_to_face = tet.field_data["cell_to_face"]
    
    if hasattr(cell_to_face, 'ndim') and cell_to_face.ndim > 0:
        cell_to_face = cell_to_face[0]
        
    if hasattr(cell_to_face, 'item'): 
        cell_to_face = cell_to_face.item()
        
    if isinstance(cell_to_face, bytes):
        cell_to_face = cell_to_face.decode('utf-8')
        
    if isinstance(cell_to_face, str):
        cell_to_face = decode_object(cell_to_face)
    
    print("DEBUG: Assignation des données du maillage...")
    # 2. ASSIGNATION DES DONNÉES DU MAILLAGE
    cells = tet.cells.reshape(-1, 5)[:, 1:]
    tet.add_field_data(cells, "cells")
    cell_vertices = tet.points
    tet.add_field_data(cell_vertices, "cell_vertices")
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]
    tet.add_field_data(faces, "faces")
    face_vertices = surface_mesh.points
    tet.add_field_data(face_vertices, "face_vertices")

    print("DEBUG: Vectorisation des normales...")
    
    # --- VECTORISATION GPU : Normales de face "most-down" ---
    cell_indices = np.array(list(cell_to_face.keys()), dtype=np.int64) # <-- Fini l'erreur ici !
    face_lists = list(cell_to_face.values())

    print("DEBUG: Padding des indices...")

    # Padding des indices pour avoir une matrice 2D homogène (remplie de -1)
    max_faces = max(len(f) for f in face_lists) if face_lists else 0
    padded_faces = np.full((len(face_lists), max_faces), -1, dtype=np.int64)
    for i, f in enumerate(face_lists):
        padded_faces[i, :len(f)] = f
        
    cell_indices_gpu = cp.asarray(cell_indices)
    padded_faces_gpu = cp.asarray(padded_faces)

    # Astuce du -1 : On ajoute une ligne avec Z = +inf à la fin.
    # Le padding (-1) pointera vers cette ligne, et cp.argmin l'ignorera.
    normals_np = np.vstack([surface_mesh.face_normals, [0.0, 0.0, np.inf]])
    normals_gpu = cp.asarray(normals_np)

    print("Extraction vectorisée...")
    
    # Extraction vectorisée et argmin sur l'axe Z
    z_normals_gpu = normals_gpu[padded_faces_gpu, 2]
    best_face_local_idx = cp.argmin(z_normals_gpu, axis=1)
    
    row_indices_gpu = cp.arange(len(face_lists))
    most_down_face_idx_gpu = padded_faces_gpu[row_indices_gpu, best_face_local_idx]

    face_normal_gpu = cp.full((tet.number_of_cells, 3), cp.nan, dtype=cp.float64)
    face_normal_gpu[cell_indices_gpu] = normals_gpu[most_down_face_idx_gpu]
    
    print("DEBUG: Normalisation...")

    # Normalisation sécurisée
    norms = cp.linalg.norm(face_normal_gpu, axis=1)[:, None]
    safe_norms = cp.where(norms == 0, 1.0, norms)
    face_normal_gpu = face_normal_gpu / safe_norms
    tet.cell_data['face_normal'] = face_normal_gpu.get()

    print("DEBUG: Vectorisation des centres de faces...")

    # --- VECTORISATION GPU : Centres de face "most-down" ---
    centers_np = np.vstack([surface_mesh.cell_centers().points, [0.0, 0.0, np.inf]])
    centers_gpu = cp.asarray(centers_np)

    z_centers_gpu = centers_gpu[padded_faces_gpu, 2]
    best_center_local_idx = cp.argmin(z_centers_gpu, axis=1)
    most_down_center_idx_gpu = padded_faces_gpu[row_indices_gpu, best_center_local_idx]

    face_center_gpu = cp.full((tet.number_of_cells, 3), cp.nan, dtype=cp.float64)
    face_center_gpu[cell_indices_gpu] = centers_gpu[most_down_center_idx_gpu]
    tet.cell_data['face_center'] = face_center_gpu.get()

    print("DEBUG: Centre des cellules...")
    
    # --- Centres de cellule ---
    cell_center_cpu = tet.cell_centers().points
    tet.cell_data["cell_center"] = cell_center_cpu
    cell_center_gpu = cp.asarray(cell_center_cpu)

    print("DEBUG: Bottom..")

    # --- Bottom cells ---
    bottom_cell_threshold_gpu = cp.nanmin(face_center_gpu[:, 2]) + 0.3
    bottom_cells_mask_gpu = face_center_gpu[:, 2] < bottom_cell_threshold_gpu
    bottom_cells_mask = bottom_cells_mask_gpu.get()
    tet.cell_data['is_bottom'] = bottom_cells_mask
    bottom_cells = np.where(bottom_cells_mask)[0]

    # --- overhang_angle / overhang_direction ---
    face_normal_masked_gpu = face_normal_gpu.copy()
    face_normal_masked_gpu[bottom_cells_mask_gpu] = cp.nan

    dot_vals = cp.clip(cp.dot(face_normal_masked_gpu, up_vector_gpu), -1.0, 1.0)  # évite les NaN dus à l'arrondi flottant
    overhang_angle_gpu = cp.arccos(dot_vals)
    tet.cell_data['overhang_angle'] = overhang_angle_gpu.get()

    overhang_direction_gpu = face_normal_masked_gpu[:, :2].copy()
    norms_gpu = cp.linalg.norm(overhang_direction_gpu, axis=1)[:, None]

    divided_gpu = overhang_direction_gpu / cp.where(norms_gpu > 0, norms_gpu, 1.0)
    zeros_gpu = cp.zeros_like(overhang_direction_gpu)
    positive_norm_mask_gpu = cp.broadcast_to(norms_gpu > 0, overhang_direction_gpu.shape)
    overhang_direction_gpu = cp.where(positive_norm_mask_gpu, divided_gpu, zeros_gpu)

    tet.cell_data['overhang_direction'] = overhang_direction_gpu.get()

    print("DEBUG: in_air...")
    
    # --- in_air : Dijkstra reste sur CPU pour garantir l'exactitude stricte des chemins ---
    IN_AIR_THRESHOLD = 1
    
    _, paths_to_bottom = nx.multi_source_dijkstra(cell_neighbour_graph, set(bottom_cells))

    # Optimisation de la construction de la matrice des chemins
    paths_list = list(paths_to_bottom.values())
    max_path_len = max(len(p) for p in paths_list) if paths_list else 0
    padded_paths = np.full((tet.number_of_cells, max_path_len), -1, dtype=np.int32)
    path_lengths = np.zeros(tet.number_of_cells, dtype=np.int32)

    for cell_index, path in paths_to_bottom.items():
        path_len = len(path)
        padded_paths[cell_index, :path_len] = path
        path_lengths[cell_index] = path_len
        
    tet.cell_data['path_to_bottom'] = padded_paths

    padded_paths_gpu = cp.asarray(padded_paths)
    path_lengths_gpu = cp.asarray(path_lengths)

    valid_mask_gpu = padded_paths_gpu >= 0
    safe_indices_gpu = cp.where(valid_mask_gpu, padded_paths_gpu, 0)
    
    path_heights_gpu = cell_center_gpu[safe_indices_gpu, 2]
    path_heights_gpu = cp.where(valid_mask_gpu, path_heights_gpu, -cp.inf)

    own_height_gpu = cell_center_gpu[:, 2][:, None] + IN_AIR_THRESHOLD
    higher_than_self_gpu = cp.any(path_heights_gpu > own_height_gpu, axis=1)
    has_path_gpu = path_lengths_gpu > 1

    tet.cell_data['in_air'] = (higher_than_self_gpu & has_path_gpu).get()

    return tet
    
def calculate_tet_attributes(tet):
    '''
    Calculate shared vertices between cells, cell to face & face to cell relations, and bottom cells of the tetrahedral mesh.
    '''
    surface_mesh = tet.extract_surface(algorithm=None)  # [CPU] VTK

    cells = tet.cells.reshape(-1, 5)[:, 1:]
    tet.add_field_data(cells, "cells")
    cell_vertices = tet.points
    tet.add_field_data(cell_vertices, "cell_vertices")
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]
    tet.add_field_data(faces, "faces")
    face_vertices = surface_mesh.points
    tet.add_field_data(face_vertices, "face_vertices")

    shared_vertices = []
    # Conversion du tqdm en boucle manuelle pour "Extraction sommets partagés"
    neighbours = tet.field_data["cell_point_neighbours"]
    total_neigh = len(neighbours)
    for i, (cell_1, cell_2) in enumerate(neighbours):
        if i % 1 == 0:
            print(json.dumps({"type": "progress", "task": "shared_vertices_extraction", "current": i, "total": total_neigh}))
            sys.stdout.flush()
        
        shared_vertices_these_faces = np.intersect1d(cells[cell_1], cells[cell_2])
        for vertex in shared_vertices_these_faces:
            shared_vertices.append({
                    "cell_1_index": cell_1,
                    "cell_2_index": cell_2,
                    "cell_1_vertex_index": np.where(cells[cell_1] == vertex)[0][0],
                    "cell_2_vertex_index": np.where(cells[cell_2] == vertex)[0][0],
                })

    cell_to_face = {}
    face_to_cell = {face_index: [] for face_index in range(len(faces))}
    cell_to_face_vertices = {}
    face_to_cell_vertices = {}

    # Mappage sommets/faces
    cv_data = tet.field_data["cell_vertices"].reshape(-1, 3)
    total_cv = len(cv_data)
    for i, cell_vertex in enumerate(cv_data):
        if i % 1 == 0:
            print(json.dumps({"type": "progress", "task": "vertex_face_mapping", "current": i, "total": total_cv}))
            sys.stdout.flush()
            
        face_vertex_index = np.where((face_vertices == cell_vertex).all(axis=1))[0]
        if len(face_vertex_index) == 1:
            cell_to_face_vertices[i] = face_vertex_index[0]
            face_to_cell_vertices[face_vertex_index[0]] = i

    # Mappage cellules/faces
    cells_data = tet.field_data["cells"]
    total_cells = len(cells_data)
    for i, cell in enumerate(cells_data):
        if i % 1 == 0:
            print(json.dumps({"type": "progress", "task": "cell_face_mapping", "current": i, "total": total_cells}))
            sys.stdout.flush()
            
        face_vertex_indices = [cell_to_face_vertices[cv_idx] for cv_idx in cell if cv_idx in cell_to_face_vertices]
        if len(face_vertex_indices) >= 3:
            extracted = surface_mesh.extract_points(face_vertex_indices, adjacent_cells=False)
            if extracted.number_of_cells >= 1:
                cell_to_face[i] = list(extracted.cell_data['vtkOriginalCellIds'])
                for face_index in extracted.cell_data['vtkOriginalCellIds']:
                    face_to_cell[face_index].append(i)

    tet.add_field_data(encode_object(cell_to_face), "cell_to_face")
    tet.add_field_data(encode_object(face_to_cell), "face_to_cell")

    tet.cell_data['has_face'] = np.zeros(tet.number_of_cells)
    for cell_index, face_indices in cell_to_face.items():
        tet.cell_data['has_face'][cell_index] = 1

    tet = update_tet_attributes(tet, up_vector_gpu, cell_neighbour_graph)
    
    bottom_cells_mask = tet.cell_data['is_bottom']
    bottom_cells = np.where(bottom_cells_mask)[0]
    tet.cell_data['overhang_angle'][bottom_cells] = np.nan

    return tet, bottom_cells_mask, bottom_cells
    
def SETUP_END():
    global bottom_cells, bottom_cells_mask, input_tet, undeformed_tet
    bottom_cells_mask = None
    bottom_cells = None
    input_tet, bottom_cells_mask, bottom_cells = calculate_tet_attributes(input_tet)
    print(f"[GPU/CPU] Attributs de cellule calculés : {int(bottom_cells_mask.sum())} bottom cells sur {input_tet.number_of_cells}")

    # -----------------------------------------------------------------------
    # Groupes de bottom cells connectés — reste CPU (nx.connected_components,
    # structure topologique, pas de calcul numérique vectorisable ici).
    # -----------------------------------------------------------------------

    bottom_cell_graph = nx.Graph()
    for cell_index in bottom_cells:
        bottom_cell_graph.add_node(cell_index)
    cell_point_neighbour_dict = cell_neighbour_dict["point"]
    for cell_index in bottom_cells:
        for neighbour in cell_point_neighbour_dict[cell_index]:
            if neighbour in bottom_cells:
                bottom_cell_graph.add_edge(cell_index, neighbour)

    bottom_cell_groups = [list(x) for x in list(nx.connected_components(bottom_cell_graph))]
    input_tet.field_data["bottom_cell_groups"] = encode_object(bottom_cell_groups)

    undeformed_tet = input_tet.copy()
    print(f"[CPU] undeformed_tet prêt : {undeformed_tet.number_of_cells} cellules")
    print(f"[GLOBAL] Import des données en mémoire flash à accès aléatoire et fermeture du noyau...")
    print(f"[GLOBAL] Terminé. Opération réussie.")
    
def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    import numpy as np
    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1))
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T)
    return ctr, svd(M)[0][:,-1]

def calculate_rotation_matrices(tet, rotation_field):
    centers_xy = tet.cell_data["cell_center"].copy()
    centers_xy[:, 2] = 0  # On force Z à 0 pour rester en 3D pure
    tangential_vectors = np.cross(np.array([0, 0, 1]), centers_xy)
    
    norms = np.linalg.norm(tangential_vectors, axis=1)
    # Eviter la division par zéro proprement
    mask = norms > 0
    tangential_vectors[mask] /= norms[mask, None]
    tangential_vectors[~mask] = [1, 0, 0]

    rotation_matrices = R.from_rotvec(rotation_field[:, None] * tangential_vectors).as_matrix()
    return rotation_matrices
    
def calculate_path_length_to_base_gradient(tet, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO):
    '''
    Calculate the path length to base gradient for each cell in the tetrahedral mesh with respect to the radial direction. This is used to determine the optimal rotation direction for each cell.

    returns: path_length_to_base_gradient. A scalar for each cell in the tetrahedral mesh. This is the gradient in the radial direction of the path length to the closest bottom cell.
    '''
    global bottom_cells
    
    # calculate initial rotation direction for each face
    path_length_to_base_gradient = np.zeros((tet.number_of_cells)) # this is a scalar with respect to the radial direction. ie the vector pointing to the cell center

    # find the path length for every overhang cell to a bottom cell
    cell_distance_to_bottom = np.empty((tet.number_of_cells))
    cell_distance_to_bottom[:] = np.nan
    distances_to_bottom, paths_to_bottom = nx.multi_source_dijkstra(cell_neighbour_graph, set(bottom_cells))# set([x[0] for x in tet.field_data["bottom_cell_groups"]]))
    closest_bottom_cell_indices = np.zeros((tet.number_of_cells), dtype=int)
    for cell_index in range(tet.number_of_cells):
        face_normal = tet.cell_data["face_normal"][cell_index]

        cell_is_overhang = np.arccos(np.dot(face_normal, [0,0,1])) > np.deg2rad(90+MAX_OVERHANG)
        if cell_is_overhang and cell_index not in bottom_cells:
            closest_bottom_cell_indices[cell_index] = paths_to_bottom[cell_index][0]
            cell_distance_to_bottom[cell_index] = distances_to_bottom[cell_index]

    tet.cell_data["cell_distance_to_bottom"] = cell_distance_to_bottom

    # calculate the gradient of path length to base for each cell
    for cell_index in range(tet.number_of_cells):
        if not np.isnan(cell_distance_to_bottom[cell_index]):
            local_cells = cell_neighbour_dict["edge"][cell_index]
            local_cells = np.hstack((local_cells, cell_index))
            # add neighbours neighbours
            # local_cells = neighbours.copy()
            # for neighbour in neighbours:
            #     local_cells.extend(cell_neighbour_dict["point"][neighbour])
            # local_cells = np.array(list(set(local_cells)))

            local_cell_path_lengths = [cell_distance_to_bottom[local_cell] for local_cell in local_cells]
            local_cell_path_lengths = np.array(local_cell_path_lengths)

            # remove neighbours with path length of nan
            local_cells = np.array(local_cells)[~np.isnan(local_cell_path_lengths)]
            local_cell_path_lengths = local_cell_path_lengths[~np.isnan(local_cell_path_lengths)]

            # if there are less than 3 neighbours with path length, roll to the closest bottom cell
            if len(local_cell_path_lengths) < 3:
                location_to_roll_to = tet.cell_data["cell_center"][closest_bottom_cell_indices[cell_index], :2]

                direction_to_bottom = location_to_roll_to - tet.cell_data["cell_center"][cell_index, :2]
                direction_to_bottom /= np.linalg.norm(direction_to_bottom)

                cell_center = tet.cell_data["cell_center"][cell_index, :2].copy()
                cell_center /= np.linalg.norm(cell_center)

                optimal_rotation_direction = np.dot(cell_center, direction_to_bottom) / np.abs(np.dot(cell_center, direction_to_bottom))
                if np.isnan(optimal_rotation_direction):
                    optimal_rotation_direction = 0

                path_length_to_base_gradient[cell_index] = optimal_rotation_direction

            # if there are 3 or more neighbours with path length, calculate the gradient in the radial direction
            # and use that as the optimal rotation direction
            else:
                points = np.hstack((tet.cell_data["cell_center"][local_cells, :2], local_cell_path_lengths[:, None]))
                _, plane_normal = planeFit(points.T)

                cell_center_direction_normalized = tet.cell_data["cell_center"][cell_index, :2] / np.linalg.norm(tet.cell_data["cell_center"][cell_index, :2])
                gradient_in_radial_direction = np.dot(cell_center_direction_normalized, plane_normal[:2])

                # if the gradient is nan, use the average of the neighbours
                if np.isnan(gradient_in_radial_direction):
                    gradient_in_radial_direction = np.mean(path_length_to_base_gradient[local_cells][~np.isnan(path_length_to_base_gradient[local_cells])])
                    if np.isnan(gradient_in_radial_direction):
                        gradient_in_radial_direction = 0

                path_length_to_base_gradient[cell_index] = gradient_in_radial_direction

    # smooth path_length_to_base_gradient with neighbours
    # not needed because we do neighbour difference minimization in the optimization step?
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
                    local_cells = local_cells[path_length_to_base_gradient[local_cells]!=0]
                    smoothed_path_length_to_base_gradient[cell_index] = np.mean(path_length_to_base_gradient[local_cells])

        path_length_to_base_gradient = smoothed_path_length_to_base_gradient

    # replace 0 with nan
    if not SET_INITIAL_ROTATION_TO_ZERO:
        path_length_to_base_gradient[path_length_to_base_gradient == 0] = np.nan
    tet.cell_data["path_length_to_base_gradient"] = path_length_to_base_gradient # very sexy

    return path_length_to_base_gradient
    
def calculate_initial_rotation_field_gpu(tet, MAX_OVERHANG, ROTATION_MULTIPLIER, STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION):
    '''
    Calcule le champ initial entièrement sur GPU.
    '''
    overhang_angle_gpu = cp.asarray(tet.cell_data['overhang_angle'])
    in_air_gpu = cp.asarray(tet.cell_data["in_air"])
    
    initial_rotation_field_gpu = cp.abs(cp.deg2rad(90+MAX_OVERHANG) - overhang_angle_gpu)

    # calculate_path_length_to_base_gradient est CPU (utilise networkx) et écrit
    # elle-même tet.cell_data["cell_distance_to_bottom"] et tet.cell_data["path_length_to_base_gradient"]
    # en effet de bord sur l'objet tet. On récupère juste sa valeur de retour ici.
    path_length_to_base_gradient = calculate_path_length_to_base_gradient(tet, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO)
    path_length_to_base_gradient_gpu = cp.asarray(path_length_to_base_gradient)

    if STEEP_OVERHANG_COMPENSATION:
        initial_rotation_field_gpu[in_air_gpu] += 2.0 * (cp.deg2rad(180) - overhang_angle_gpu[in_air_gpu])

    initial_rotation_field_gpu *= path_length_to_base_gradient_gpu
    initial_rotation_field_gpu = cp.clip(initial_rotation_field_gpu * ROTATION_MULTIPLIER, -cp.deg2rad(360), cp.deg2rad(360))
    initial_rotation_field_gpu = cp.clip(initial_rotation_field_gpu, MAX_NEG_ROTATION, MAX_POS_ROTATION)

    # Ligne réintégrée : la version CPU d'origine la faisait (tet.cell_data["initial_rotation_field"] = ...),
    # et Noyau 2 en a besoin pour certains plots. On la restaure ici pour garder le même comportement observable.
    tet.cell_data["initial_rotation_field"] = initial_rotation_field_gpu.get()

    return initial_rotation_field_gpu

def calculate_unique_vertices_rotated_gpu(tet, rotation_field):
    '''
    Rotation des sommets via la formule de Rodrigues 100% CuPy.
    '''
    rotation_field_gpu = cp.asarray(rotation_field)
    cell_centers_gpu = cp.asarray(tet.cell_data["cell_center"])
    
    centers_xy = cell_centers_gpu.copy()
    centers_xy[:, 2] = 0
    Z_axis = cp.array([0.0, 0.0, 1.0])
    
    tangential_vectors = cp.cross(Z_axis, centers_xy)
    norms = cp.linalg.norm(tangential_vectors, axis=1)
    
    mask = norms > 0
    tangential_vectors[mask] /= norms[mask, None]
    tangential_vectors[~mask] = cp.array([1.0, 0.0, 0.0])

    cells = tet.cells.reshape(-1, 5)[:, 1:]
    vertices_gpu = cp.asarray(tet.points[cells])
    
    shifted_gpu = vertices_gpu - cell_centers_gpu[:, None, :]
    
    theta = rotation_field_gpu[:, None, None]
    k = tangential_vectors[:, None, :]
    
    cos_theta = cp.cos(theta)
    sin_theta = cp.sin(theta)
    
    k_broadcast = cp.broadcast_to(k, shifted_gpu.shape)
    
    k_cross_v = cp.cross(k_broadcast, shifted_gpu, axis=2)
    k_dot_v = cp.sum(k_broadcast * shifted_gpu, axis=2, keepdims=True)
    
    v_rot = shifted_gpu * cos_theta + k_cross_v * sin_theta + k_broadcast * k_dot_v * (1.0 - cos_theta)
    
    unique_vertices_rotated_gpu = cell_centers_gpu[:, None, :] + v_rot
    
    return unique_vertices_rotated_gpu.get()

def apply_rotation_field_unique_vertices(tet, rotation_field):
    unique_vertices_rotated = calculate_unique_vertices_rotated_gpu(tet, rotation_field)
    
    num_cells = tet.number_of_cells
    unique_cells = np.empty((num_cells, 5), dtype=int)
    unique_cells[:, 0] = 4
    unique_cells[:, 1:] = np.arange(num_cells * 4).reshape(-1, 4)

    new_tet = pv.UnstructuredGrid(unique_cells.flatten(), np.full(num_cells, pv.CellType.TETRA), unique_vertices_rotated.reshape(-1, 3))
    
    for key in tet.cell_data.keys():
        new_tet.cell_data[key] = tet.cell_data[key].copy()
        
    return new_tet

def apply_rotation_field_unique_vertices(tet, rotation_field):
    unique_vertices_rotated = calculate_unique_vertices_rotated_gpu(tet, rotation_field)
    
    num_cells = tet.number_of_cells
    unique_cells = np.empty((num_cells, 5), dtype=int)
    unique_cells[:, 0] = 4
    unique_cells[:, 1:] = np.arange(num_cells * 4).reshape(-1, 4)

    new_tet = pv.UnstructuredGrid(unique_cells.flatten(), np.full(num_cells, pv.CellType.TETRA), unique_vertices_rotated.reshape(-1, 3))
    
    for key in tet.cell_data.keys():
        new_tet.cell_data[key] = tet.cell_data[key].copy()
        
    return new_tet

# ==========================================
# FONCTION D'OPTIMISATION PRINCIPALE
# ==========================================

def optimize_rotations(tet, NEIGHBOUR_LOSS_WEIGHT, MAX_OVERHANG, ROTATION_MULTIPLIER, ITERATIONS, SAVE_GIF, STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION):
    
    initial_rotation_field_gpu = calculate_initial_rotation_field_gpu(tet, MAX_OVERHANG, ROTATION_MULTIPLIER, STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION)
    
    cell_face_neighbours_gpu = cp.asarray(tet.field_data["cell_face_neighbours"], dtype=cp.int32)
    num_cells = tet.number_of_cells
    W = cp.float32(NEIGHBOUR_LOSS_WEIGHT)
    
    cell_1 = cell_face_neighbours_gpu[:, 0]
    cell_2 = cell_face_neighbours_gpu[:, 1]
    
    rows_A = cp.concatenate((cell_1, cell_2))
    cols_A = cp.concatenate((cell_2, cell_1))
    data_A = cp.full(rows_A.shape, -W, dtype=cp.float32)
    
    degrees = cp.zeros(num_cells, dtype=cp.float32)
    scatter_add(degrees, cell_1, W)
    scatter_add(degrees, cell_2, W)
    
    valid_mask = ~cp.isnan(initial_rotation_field_gpu)
    degrees[valid_mask] += 1.0
    
    rows_D = cp.arange(num_cells, dtype=cp.int32)
    cols_D = rows_D
    data_D = degrees
    
    all_rows = cp.concatenate((rows_A, rows_D))
    all_cols = cp.concatenate((cols_A, cols_D))
    all_data = cp.concatenate((data_A, data_D))
    
    M_gpu = cusparse.coo_matrix((all_data, (all_rows, all_cols)), shape=(num_cells, num_cells)).tocsr()
    
    b_gpu = cp.zeros(num_cells, dtype=cp.float32)
    b_gpu[valid_mask] = initial_rotation_field_gpu[valid_mask]

    x0 = cp.zeros(num_cells, dtype=cp.float32)
    x_gpu, info = culinalg.cg(M_gpu, b_gpu, x0=x0, maxiter=ITERATIONS, atol=1e-6)

    if info > 0:
        print(f"⚠️ Tolérance non atteinte, arrêt à {ITERATIONS} itérations (le résultat reste exploitable).")
    elif info < 0:
        print("❌ Erreur critique dans le solveur CuPy.")
    
    smoothed_rotation_field = x_gpu.get()
    
    return smoothed_rotation_field
    
def START_TET_MESH_DEFORM():
    global undeformed_tet
    global initial_rotation_field
    global rotation_field          # <-- ajoute cette ligne
    
    rotation_field = optimize_rotations(
        undeformed_tet,
        NEIGHBOUR_LOSS_WEIGHT,
        MAX_OVERHANG,
        ROTATION_MULTIPLIER,
        ITERATIONS,
        SAVE_GIF,
        STEEP_OVERHANG_COMPENSATION,
        INITIAL_ROTATION_FIELD_SMOOTHING,
        SET_INITIAL_ROTATION_TO_ZERO,
        MAX_POS_ROTATION,
        MAX_NEG_ROTATION
    )
    undeformed_tet_with_rotated_tetrahedrons = apply_rotation_field_unique_vertices(undeformed_tet, rotation_field)
    undeformed_tet_with_rotated_tetrahedrons.cell_data["rotation_field"] = rotation_field
    # undeformed_tet_with_rotated_tetrahedrons.plot(scalars="rotation_field")
    
def PATH_PROCESSING_SETUP():
    print("Entrée dans PATH_PROCESSING_SETUP...")
    # view the initial rotation field we are trying to optimize towards
    # tet.extract_cells(tet.cell_data['overhang_angle'] > np.deg2rad(90 + MAX_OVERHANG)).plot(scalars="initial_rotation_field")
    # tet.cell_data['overhang_angle'] > np.deg2rad(90 + MAX_OVERHANG)
    # undeformed_tet.plot(scalars="initial_rotation_field")
    # undeformed_tet.plot(scalars="in_air")
    # undeformed_tet.plot(scalars="overhang_angle")
    # undeformed_tet.plot(scalars="path_length_to_base_gradient")
    # undeformed_tet.plot(scalars=new_tet1.cell_data['rotation_field'])
    # undeformed_tet.plot(scalars="cell_distance_to_bottom", cpos=[-0.5, -1, -1]) 
    
def show_path_to_base(tet, cell_index, plotter=pv.Plotter()):
    path_to_bottom = tet.cell_data['path_to_bottom'][cell_index]
    first_negative_index = np.where(path_to_bottom == -1)[0][0]
    path_to_bottom = path_to_bottom[:first_negative_index]
    print(path_to_bottom)

    for i in range(len(path_to_bottom)-1):
        path_to_base = pv.Line(tet.cell_data["cell_center"][path_to_bottom[i]], tet.cell_data["cell_center"][path_to_bottom[i+1]])
        plotter.add_mesh(path_to_base, color="red")

    plotter.add_mesh(tet, opacity=0.2)
    plotter.show()
    
def show_path_to_base_gradient_calculation(tet, cell_indices):
    plotter = pv.Plotter()
    scalar = np.full(tet.number_of_cells, 0.0)

    for cell_index in cell_indices:
        local_cells = cell_neighbour_dict["edge"][cell_index]
        local_cells = np.hstack((local_cells, cell_index))
        local_cells_with_path_lengths = [x for x in local_cells if not np.isnan(tet.cell_data['cell_distance_to_bottom'][x])]
        path_lengths = tet.cell_data['cell_distance_to_bottom'][local_cells_with_path_lengths]

        cell_centers = tet.cell_data["cell_center"][local_cells_with_path_lengths]
        cell_centers_z_is_path_length = cell_centers.copy()
        cell_centers_z_is_path_length[:, 2] = path_lengths - np.min(path_lengths) + np.min(cell_centers[:, 2])
        points = pv.PolyData(cell_centers_z_is_path_length)
        glyph = points.glyph(geom=pv.Sphere(theta_resolution=8, phi_resolution=8, radius=0.1))

        if len(cell_centers_z_is_path_length) < 3:
            continue

        p, n = planeFit(cell_centers_z_is_path_length.T)
        plane = pv.Plane(center=p, direction=n, i_size=4, j_size=4)
        normal_arrow = pv.Arrow(start=p, direction=n, scale=1.5)

        # reflect arrow across xy plane
        reflected_n = -n - 2 * (np.dot(-n, up_vector)) * up_vector

        # extract radial component
        cell_center_direction_normalized = tet.cell_data["cell_center"][cell_index, :2] / np.linalg.norm(tet.cell_data["cell_center"][cell_index, :2])
        gradient_in_radial_direction = np.dot(cell_center_direction_normalized, reflected_n[:2]) * cell_center_direction_normalized
        nozzle_arrow = pv.Arrow(start=tet.cell_data["cell_center"][cell_index], direction=np.hstack((gradient_in_radial_direction, reflected_n[2])), scale=1)


        # plotter.add_mesh(glyph, color="blue", opacity=0.5, )
        # plotter.add_mesh(plane, color="green", opacity=0.2)
        # plotter.add_mesh(normal_arrow, color="orange")
        plotter.add_mesh(nozzle_arrow, color="red")
    scalar[local_cells_with_path_lengths] = 0.5
    scalar[cell_index] = 1
    plotter.add_mesh(tet, opacity=0.2, scalars=scalar, cmap="binary")
    plotter.show()
    
def show_dijkstras(tet, cell_index):
    plotter = pv.Plotter()#window_size=[3840, 2160])
    lines = []
    for neighbour in tet.field_data["cell_face_neighbours"]:
        lines += [tet.cell_data["cell_center"][neighbour[0]], tet.cell_data["cell_center"][neighbour[1]]]
    mesh = pv.line_segments_from_points(lines)
    plotter.add_mesh(mesh, color="grey", opacity=0.4)

    points = pv.PolyData(tet.cell_data["cell_center"])
    glyph = points.glyph(geom=pv.Sphere(theta_resolution=8, phi_resolution=8, radius=0.1))
    plotter.add_mesh(glyph, color="red", opacity=0.2)

    if cell_index is None:
        plotter.show()
        return

    plotter.camera_position = "xz"

    # run dijkstra's algorithm and visualize
    plotter.open_gif(f'gifs/dijkstra.gif')
    distances, paths = nx.single_source_dijkstra(cell_neighbour_graph, cell_index)
    dijkstra_actors = []
    for i in np.arange(0, tet.cell_data['cell_distance_to_bottom'][cell_index], 0.2):
        nodes_in_range = [node for node, distance in distances.items() if distance < i]
        if len(nodes_in_range) == 0:
            continue
        if set(nodes_in_range) & set(bottom_cells):
            break
        points = pv.PolyData(tet.cell_data["cell_center"][nodes_in_range])
        glyph = points.glyph(geom=pv.Sphere(theta_resolution=8, phi_resolution=8, radius=0.2))
        actor = plotter.add_mesh(glyph, color="blue", opacity=0.4)
        dijkstra_actors.append(actor)
        plotter.write_frame()


    path_to_bottom = tet.cell_data['path_to_bottom'][cell_index]
    first_negative_index = np.where(path_to_bottom == -1)[0][0]
    path_to_bottom = path_to_bottom[:first_negative_index]
    print(path_to_bottom)

    for i in range(len(path_to_bottom)-1):
        path_to_base = pv.Line(tet.cell_data["cell_center"][path_to_bottom[i]], tet.cell_data["cell_center"][path_to_bottom[i+1]])
        plotter.add_mesh(path_to_base, color="blue", line_width=5)

    for actor in dijkstra_actors:
        plotter.remove_actor(actor)

    plotter.write_frame()

    plotter.show()
    plotter.close()
    
def N_MATRIX_SETUP():
    global undeformed_tet

    # show_dijkstras(undeformed_tet, np.where(undeformed_tet.cell_data['has_face'].astype(bool) & (undeformed_tet.cell_data["overhang_angle"] > 3))[0][1])
    # show_path_to_base_gradient_calculation(undeformed_tet, [np.where(undeformed_tet.cell_data['has_face'].astype(bool) & (undeformed_tet.cell_data["overhang_angle"] > 3))[0][1]])
    # show_path_to_base_gradient_calculation(undeformed_tet, np.where(undeformed_tet.cell_data['has_face'].astype(bool) & (undeformed_tet.cell_data["overhang_angle"] > np.deg2rad(90 + MAX_OVERHANG)))[0])
    global N
    
    N = cp.eye(4) - 1/4 * cp.ones((4, 4))
    save_gif_i = 0
    
def calculate_deformation(tet, rotation_field, ITERATIONS, SAVE_GIF):
    global N
    
    '''
    Find the optimal deformation of the tetrahedral mesh to make cells have the same rotation as
    the given rotation field.
    
    100% GPU Custom Adam Optimizer.
    '''
    # 1. On charge tout en VRAM immédiatement
    vertices_gpu = cp.asarray(tet.points.copy())
    cells_gpu = cp.asarray(tet.field_data["cells"])
    old_vertices_gpu = cp.asarray(tet.field_data["cell_vertices"][tet.field_data["cells"]])
    
    rot_mat_cpu = calculate_rotation_matrices(tet, rotation_field)
    rotation_matrices_gpu = cp.asarray(rot_mat_cpu)

    # 2. Précalcul de la cible sur le GPU
    # (N @ old_vertices).transpose(0, 2, 1)
    transformed_old = cp.einsum('ij,kjl->kil', N, old_vertices_gpu)
    transformed_old_T = cp.transpose(transformed_old, (0, 2, 1))
    old_vertices_transformed_gpu = cp.einsum('ijk,ikl->ijl', rotation_matrices_gpu, transformed_old_T)

    plotter = pv.Plotter(off_screen=True) if SAVE_GIF else None
    if SAVE_GIF:
        plotter.open_gif(f'gifs/calculate_deformation.gif')

    # 3. Paramètres de l'optimiseur Adam
    lr = 0.05 
    beta1 = 0.9
    beta2 = 0.999
    m = cp.zeros_like(vertices_gpu)
    v = cp.zeros_like(vertices_gpu)
    
    flat_cells = cells_gpu.flatten()

    for i in range(DEFORM_ITERATIONS):
        # --- FORWARD PASS ---
        v_cells = vertices_gpu[cells_gpu] # (C, 4, 3)
        N_v_cells = cp.einsum('ij,kjl->kil', N, v_cells) # (C, 4, 3)
        pred = cp.transpose(N_v_cells, (0, 2, 1)) # (C, 3, 4)
        
        diff = pred - old_vertices_transformed_gpu # (C, 3, 4)
        
        # --- BACKWARD PASS ---
        # dL/dV = 2 * N^T @ diff^T
        diff_T = cp.transpose(diff, (0, 2, 1)) # (C, 4, 3)
        grad_v_cells = 2 * cp.einsum('ij,kjl->kil', N.T, diff_T) # (C, 4, 3)
        
        grad = cp.zeros_like(vertices_gpu)
        flat_grad = grad_v_cells.reshape(-1, 3)
        
        # Accumulation super rapide sur GPU
        for dim in range(3):
            cupyx.scatter_add(grad[:, dim], flat_cells, flat_grad[:, dim])
        
        # --- ADAM STEP ---
        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * (grad ** 2)
        m_hat = m / (1 - beta1 ** (i + 1))
        v_hat = v / (1 - beta2 ** (i + 1))
        
        vertices_gpu -= lr * m_hat / (cp.sqrt(v_hat) + 1e-8)

        # --- GESTION DU GIF (Seul transfert CPU) ---
        if SAVE_GIF:
            global save_gif_i
            save_gif_i += 1
            if save_gif_i % 10 == 0:
                new_tet = pv.UnstructuredGrid(tet.cells, np.full(tet.number_of_cells, pv.CellType.TETRA), vertices_gpu.get())
                mesh_actor = plotter.add_mesh(new_tet)
                plotter.write_frame()
                plotter.remove_actor(mesh_actor)

    if SAVE_GIF:
        plotter.close()

    return vertices_gpu.get()
 
def EXEC_DEFORMATIONS():
    global undeformed_tet
    global deformed_tet
    global rotation_field
    global new_vertices
    
    new_vertices = calculate_deformation(undeformed_tet, rotation_field, ITERATIONS, SAVE_GIF)
    deformed_tet = pv.UnstructuredGrid(undeformed_tet.cells, np.full(undeformed_tet.number_of_cells, pv.CellType.TETRA), new_vertices)
    # deformed_tet.plot()

    for key in undeformed_tet.field_data.keys():
        deformed_tet.field_data[key] = undeformed_tet.field_data[key]
    for key in undeformed_tet.cell_data.keys():
        deformed_tet.cell_data[key] = undeformed_tet.cell_data[key]
    deformed_tet = update_tet_attributes(deformed_tet, up_vector_gpu, cell_neighbour_graph)
    
def UT_COPY():
   global undeformed_tet
   global deformed_tet
   undeformed_tet = deformed_tet.copy()
   
def SAVE_STL():
    x_min, x_max, y_min, y_max, z_min, z_max = deformed_tet.bounds
    offsets_applied = np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min])
    deformed_tet.points -= offsets_applied

    deformed_tet.extract_surface().save(f'{OUTPUT_MODEL_DIR}/model_deformed_tet.stl')
    
def SAVE_DEFORM():
    global deformed_tet
    
    #DO NOT RUN before finished all iterations
    # save to pickle
    with open(f'{PKL_DIR}/deformed_model.pkl', 'wb') as f:
        pickle.dump(deformed_tet, f)
        
#--------------------------------------------------------------------------------------
#    ATTENTION, DÉCOUPER LE MODÈLE AVEC UNE TRANCHEUSE CLASSIQUE AVANT DE CONTINUER
#--------------------------------------------------------------------------------------

def OPEN_DEFORM():
    global deformed_tet
    
    print(json.dumps({"type": "info", "message": "Chargement du modèle déformé..."}))
    sys.stdout.flush()
    
    deformed_tet = pickle.load(open(f'{PKL_DIR}/deformed_model.pkl', 'rb'))
    
    print(json.dumps({"type": "info", "message": "Modèle chargé avec succès."}))
    sys.stdout.flush()
    
def EXEC_CALC_TET_ATTR():
    global deformed_tet
    
    print(json.dumps({"type": "info", "message": "Début du calcul des attributs du tétraèdre..."}))
    sys.stdout.flush()
    
    deformed_tet, _, _ = calculate_tet_attributes(deformed_tet)
    
    print(json.dumps({"type": "info", "message": "Calcul des attributs terminé."}))
    sys.stdout.flush()
    
def calc_barycentric_coordinates(tet_a, tet_b, tet_c, tet_d, point):
    '''
    Calculate the barycentric coordinates of a point in a tetrahedron. This is used to interpolate
    parameters from the vertices of the tetrahedron to a point within the tetrhedron.
    '''

    total_volume = tetrahedron_volume(tet_a, tet_b, tet_c, tet_d)

    if total_volume == 0:
        raise ValueError("The points do not form a valid tetrahedron (zero volume).")

    # Calculate the sub-volumes for each face
    vol_a = tetrahedron_volume(point, tet_b, tet_c, tet_d)
    vol_b = tetrahedron_volume(point, tet_a, tet_c, tet_d)
    vol_c = tetrahedron_volume(point, tet_a, tet_b, tet_d)
    vol_d = tetrahedron_volume(point, tet_a, tet_b, tet_c)

    # Calculate barycentric coordinates as the ratio of sub-volumes to total volume
    lambda_a = vol_a / total_volume
    lambda_b = vol_b / total_volume
    lambda_c = vol_c / total_volume
    lambda_d = vol_d / total_volume

    # The barycentric coordinates should sum to 1
    return np.array([lambda_a, lambda_b, lambda_c, lambda_d])
    
def project_point_onto_plane(plane_x_axis, plane_y_axis, point):
    projected_x = np.sum(plane_x_axis * point, axis=1)
    projected_y = np.sum(plane_y_axis * point, axis=1)

    return np.array([projected_x, projected_y]).T
    
def EXEC_CALC_TET_ATTR():
    global deformed_tet

    deformed_tet, _, _ = calculate_tet_attributes(deformed_tet)
    
def PREPARE_ROTATE():
    global deformed_tet, input_tet
    global vertex_transformations
    global vertex_rotations
    
    # ---------------------------------------------------------
    # 1. PRÉPARATION ET ROTATIONS (Identique à ton code)
    # ---------------------------------------------------------
    vertex_transformations = deformed_tet.points - input_tet.points

    # calculate tangential vectors (axis of rotation) for each cell
    centers_xy = np.zeros_like(input_tet.cell_data["cell_center"])
    centers_xy[:, :2] = input_tet.cell_data["cell_center"][:, :2]  # On garde Z à 0 mais en restant en 3D

    tangential_vectors = np.cross(np.array([0, 0, 1]), centers_xy)
    # normalize
    tangential_vectors /= np.linalg.norm(tangential_vectors, axis=1)[:, None]
    # replace nan with [1,0,0]
    tangential_vectors[np.isnan(tangential_vectors).any(axis=1)] = [1, 0, 0]

    num_cells_per_vertex = np.zeros((input_tet.number_of_points))
    for cell_index, cell in enumerate(input_tet.field_data["cells"]):
        num_cells_per_vertex[cell] += 1
        
    vertex_rotations = np.zeros((deformed_tet.number_of_points))
    cell_rotations = np.zeros((deformed_tet.number_of_cells))

    total_cells = len(deformed_tet.field_data["cells"])

    for cell_index, cell in enumerate(deformed_tet.field_data["cells"]):
        # 1. Progression
        if cell_index % 50 == 0:
            msg = {"type": "progress", "task": "overhang_check", "current": cell_index, "total": total_cells}
            sys.stdout.write(json.dumps(msg) + "\n")
            sys.stdout.flush()
            
        # 2. RÉCUPÉRATION LOCALE (Fini les globales !)
        # On extrait les sommets de la cellule actuelle uniquement
        new_vertices = deformed_tet.field_data["cell_vertices"][cell] 
        old_vertices = input_tet.field_data["cell_vertices"][cell]
        
        new_cell_center = deformed_tet.cell_data["cell_center"][cell_index]
        old_cell_center = input_tet.cell_data["cell_center"][cell_index]

        # 3. Calcul local (ne modifie pas l'original, travaille sur une copie)
        new_verts_centered = new_vertices - new_cell_center
        old_verts_centered = old_vertices - old_cell_center

        plane_x_vector = old_cell_center[:2] / (np.linalg.norm(old_cell_center[:2]) + 1e-9)
        plane_x_vector = np.array([plane_x_vector[0], plane_x_vector[1], 0])
        plane_y_vector = np.array([0, 0, 1])

        # 4. Projection
        new_vertices_projected = project_point_onto_plane(plane_x_vector, plane_y_vector, new_verts_centered)
        old_vertices_projected = project_point_onto_plane(plane_x_vector, plane_y_vector, old_verts_centered)
        
        # 5. Calcul de rotation
        covariance_matrix = np.dot(new_vertices_projected.T, old_vertices_projected)
        U, _, Vt = np.linalg.svd(covariance_matrix)
        rotation_matrix = np.dot(U, Vt)

        rotation = -np.arccos(min(max(rotation_matrix[0, 0], -1), 1))
        if rotation_matrix[1, 0] < 0:
            rotation = -rotation

        rotation = max(min(rotation, np.deg2rad(MAX_ROTATION)), np.deg2rad(MIN_ROTATION))
        cell_rotations[cell_index] = rotation

        for vertex_index in cell:
            vertex_rotations[vertex_index] += rotation / num_cells_per_vertex[vertex_index]
            
def CALCULATE_VECTORIZED_ZSQUISH():
    global deformed_tet, input_tet
    global z_squish_scales
    
    # ---------------------------------------------------------
    # 2. CALCUL VECTORISÉ DU Z-SQUISH (Oubli de la boucle lente !)
    # ---------------------------------------------------------
    print("Calcul vectorisé du z-squish...")
    def batch_tet_volume(v):
        # v shape: (N, 4, 3)
        v03 = v[:, 0] - v[:, 3]
        v13 = v[:, 1] - v[:, 3]
        v23 = v[:, 2] - v[:, 3]
        # |(V0-V3) . ((V1-V3) x (V2-V3))| / 6
        return np.abs(np.sum(v03 * np.cross(v13, v23), axis=1)) / 6.0

    warped_v = deformed_tet.field_data["cell_vertices"][deformed_tet.field_data["cells"]]
    unwarped_v = input_tet.field_data["cell_vertices"][input_tet.field_data["cells"]]
    z_squish_scales = batch_tet_volume(unwarped_v) / batch_tet_volume(warped_v)

def OPEN_READ_GCODE():
    global deformed_tet
    global gcode_points
    # ---------------------------------------------------------
    # 3. LECTURE DU G-CODE (Identique)
    # ---------------------------------------------------------
    pos = np.array([0., 0., 20.])
    feed = 5000
    gcode_points = []
    with open(f'{INPUT_GCODE_DIR}/model_deformed_tet.gcode', 'r') as fh:
        lines = fh.readlines() 
        
        total_lines = len(lines)

        for i, line_text in enumerate(lines):
            # Envoi de la progression toutes les 500 lignes (pour ne pas saturer le canal)
            if i % 1 == 0:
                progress_data = {
                    "type": "progress", 
                    "task": "gcode_parsing", 
                    "current": i, 
                    "total": total_lines
                }
                sys.stdout.write(json.dumps(progress_data) + "\n")
                sys.stdout.flush()
        
            line = Line(line_text)

            if not line.block.gcodes:
                continue

            for gcode in sorted(line.block.gcodes):
                if gcode.word == "G01" or gcode.word == "G00":
                    prev_pos = pos.copy()

                    if gcode.X is not None: pos[0] = gcode.X
                    if gcode.Y is not None: pos[1] = gcode.Y
                    if gcode.Z is not None: pos[2] = gcode.Z

                    inv_time_feed = None
                    for word in line.block.words:
                        if word.letter == "F":
                            feed = word.value

                    extrusion = None
                    for param in line.block.modal_params:
                        if param.letter == "E":
                            extrusion = param.value

                    delta_pos = pos - prev_pos
                    distance = np.linalg.norm(delta_pos)
                    if distance > 0:
                        num_segments = -(-distance // SEG_SIZE) 
                        seg_distance = distance/num_segments

                        time_to_complete_move = (1/feed) * seg_distance
                        inv_time_feed = 1/time_to_complete_move if time_to_complete_move != 0 else None

                        for i in range(int(num_segments)):
                            gcode_points.append({
                                "position": (prev_pos + delta_pos * (i+1) / num_segments),
                                "command": gcode.word,
                                "extrusion": extrusion/num_segments if extrusion is not None else None,
                                "inv_time_feed": inv_time_feed,
                                "move_length": seg_distance,
                                "start_position": prev_pos,
                                "end_position": pos,
                                "unsegmented_move_length": distance,
                                "after_retract": False,
                                "feed": feed
                            })
                    else:
                        time_to_complete_move = (1/feed) * distance
                        inv_time_feed = 1/time_to_complete_move if time_to_complete_move != 0 else None

                        gcode_points.append({
                            "position": pos.copy(),
                            "command": gcode.word,
                            "extrusion": extrusion,
                            "inv_time_feed": inv_time_feed,
                            "move_length": distance,
                            "unsegmented_move_length": distance,
                            "after_retract": False,
                            "feed": feed
                        })
                        
def BRC_INTERPOL_PREPROCESS():
    global deformed_tet
    global gcode_points
    global vertex_transformations
    global vertex_rotations
    global cells_array
    global precalc_new_positions
    global precalc_rotations
    # ---------------------------------------------------------
    # 4. PRÉ-CALCUL VECTORIEL DES INTERPOLATIONS
    # ---------------------------------------------------------
    print("Pré-calcul vectoriel des cellules et interpolations barycentriques... Ceci peut prendre un certain temps.")
    positions_array = np.array([pt["position"] for pt in gcode_points])
    commands_array = np.array([pt["command"] for pt in gcode_points])

    cells_array = np.array(deformed_tet.find_containing_cell(positions_array))
    closest_cells_array = np.array(deformed_tet.find_closest_cell(positions_array))

    # Remplacer les -1 par closest_cell pour les commandes G01 uniquement
    is_g01 = commands_array == "G01"
    needs_closest = (cells_array == -1) & is_g01
    cells_array[needs_closest] = closest_cells_array[needs_closest]

    N_points = len(gcode_points)
    precalc_new_positions = np.full((N_points, 3), np.nan)
    precalc_rotations = np.full(N_points, np.nan)

    valid_mask = cells_array != -1
    valid_indices = np.where(valid_mask)[0]
    valid_cells = cells_array[valid_mask]
    valid_positions = positions_array[valid_mask]

    # Extraction des sommets des cellules valides
    v_indices = deformed_tet.field_data["cells"][valid_cells] # (M, 4)
    v = deformed_tet.field_data["cell_vertices"][v_indices] # (M, 4, 3)

    # Résolution par inversion matricielle de masse (linalg.solve)
    r = valid_positions - v[:, 3, :]
    T = np.empty((len(valid_cells), 3, 3))
    T[:, :, 0] = v[:, 0, :] - v[:, 3, :]
    T[:, :, 1] = v[:, 1, :] - v[:, 3, :]
    T[:, :, 2] = v[:, 2, :] - v[:, 3, :]

    try:
        # On ajoute une dimension à r pour obtenir (M, 3, 1), puis on l'enlève du résultat avec .squeeze(-1)
        L012 = np.linalg.solve(T, r[..., np.newaxis]).squeeze(-1)
        L3 = 1.0 - L012.sum(axis=1)
        L = np.column_stack((L012, L3))
        
        # Validation du critère 1.01
        good_bary = np.sum(L, axis=1) <= 1.01
        
        final_indices = valid_indices[good_bary]
        final_v_indices = v_indices[good_bary]
        final_L = L[good_bary]
        final_positions = valid_positions[good_bary]
        
        # Transformation (somme pondérée vectorisée)
        v_transf_subset = vertex_transformations[final_v_indices]
        transformations = np.sum(v_transf_subset * final_L[:, :, None], axis=1)
        precalc_new_positions[final_indices] = final_positions - transformations
        
        # Rotation (somme pondérée vectorisée)
        v_rot_subset = vertex_rotations[final_v_indices]
        precalc_rotations[final_indices] = np.sum(v_rot_subset * final_L, axis=1)

    except np.linalg.LinAlgError:
        # Sécurité au cas où des tétraèdres dégénérés feraient crasher le solveur
        print("Matrices singulières détectées. Utilisation du fallback classique")
        total_steps = len(valid_indices)

        for i, idx in enumerate(valid_indices):
            # 1. Envoi de la télémétrie à ton interface TypeScript
            # On envoie une mise à jour tous les 100 pas pour ne pas saturer le canal
            if i % 1 == 0 or i == total_steps - 1:
                progress = {
                    "type": "progress", 
                    "task": "overhang_check", 
                    "current": i, 
                    "total": total_steps
                }
                sys.stdout.write(json.dumps(progress) + "\n")
                sys.stdout.flush()
            
            # 2. Ton calcul original
            c_idx = valid_cells[i]
            c_v = v[i]
            pos_i = valid_positions[i]
            b_coords = calc_barycentric_coordinates(c_v[0], c_v[1], c_v[2], c_v[3], pos_i)
            
            if np.sum(b_coords) <= 1.01:
                transf = np.sum(vertex_transformations[v_indices[i]] * b_coords[:, None], axis=0)
                precalc_new_positions[idx] = pos_i - transf
                precalc_rotations[idx] = np.sum(vertex_rotations[v_indices[i]] * b_coords)

def ASSEMBLY_LOOP():
    global gcode_points
    global cells_array
    global precalc_new_positions
    global precalc_rotations
    global z_squish_scales
    global new_gcode_points
    # ---------------------------------------------------------
    # 5. BOUCLE D'ASSEMBLAGE (Super rapide maintenant !)
    # ---------------------------------------------------------
    new_gcode_points = []
    prev_new_position = None
    travelling_over_air = False
    travelling = False
    prev_rotation = 0
    prev_travelling = False
    prev_command = "G00"
    ROTATION_AVERAGING_ALPHA = 0.2
    RETRACTION_LENGTH = 1.0
    ROTATION_MAX_DELTA = np.deg2rad(1)
    MAX_EXTRUSION_MULTIPLIER = 10
    lost_vertices = []
    highest_printed_point = 0

    total_gcode = len(gcode_points)

    # Remplacement de la barre de prog par une boucle enumerate pure
    for i, gcode_point in enumerate(gcode_points):
        
        # Envoi du message JSON toutes les 500 itérations (ou autre fréquence de ton choix)
        # pour ne pas noyer ton interface sous les messages.
        if i % 1 == 0:
            progress_msg = {"type": "progress", "task": "dynamic_transform", "current": i, "total": total_gcode}
            sys.stdout.write(json.dumps(progress_msg) + "\n")
            sys.stdout.flush()

        # ... tout ton code de transformation reste INCHANGÉ en dessous ...
        position = gcode_point["position"]
        command = gcode_point["command"]
        inv_time_feed = gcode_point["inv_time_feed"]
        extrusion = gcode_point["extrusion"]
        containing_cell_index = cells_array[i]

        new_position = precalc_new_positions[i]
        rotation = precalc_rotations[i]
        
        if np.isnan(new_position).any() or np.isnan(rotation):
            new_position = None
            rotation = None

        dont_smooth_rotation = False
        
        if new_position is None:
            if command == "G01":
                lost_vertices.append(position)
                continue
            elif command == "G00" and not travelling_over_air and prev_new_position is not None:
                new_position = np.array([prev_new_position[0], prev_new_position[1], highest_printed_point])
                rotation = max(min(prev_rotation, np.deg2rad(45)), np.deg2rad(-45))
                dont_smooth_rotation = True 
                travelling_over_air = True
            elif travelling_over_air:
                continue
            else:
                continue
        else:
            if travelling_over_air:
                new_position[2] = highest_printed_point
                rotation = max(min(rotation, np.deg2rad(45)), np.deg2rad(-45))
                dont_smooth_rotation = True 
            travelling_over_air = False

        extrusion_multiplier = 1
        if extrusion is not None and extrusion != RETRACTION_LENGTH and extrusion != -RETRACTION_LENGTH:
            extrusion_multiplier = extrusion_multiplier * z_squish_scales[containing_cell_index]
            extrusion = extrusion * min(extrusion_multiplier, MAX_EXTRUSION_MULTIPLIER)
        elif extrusion == -RETRACTION_LENGTH:
            travelling = True
        elif extrusion == RETRACTION_LENGTH:
            travelling = False
            
        if prev_rotation is not None and not dont_smooth_rotation:
            rotation = ROTATION_AVERAGING_ALPHA * rotation + (1 - ROTATION_AVERAGING_ALPHA) * prev_rotation

        if prev_rotation is not None and prev_new_position is not None and np.abs(rotation - prev_rotation) > ROTATION_MAX_DELTA:
            delta_rotation = rotation - prev_rotation
            num_interpolations = int(np.abs(delta_rotation) / ROTATION_MAX_DELTA) + 1
            delta_pos = new_position - prev_new_position
            for j in range(num_interpolations):
                new_gcode_points.append({
                    "position": prev_new_position + (delta_pos * ((j+1) / num_interpolations)),
                    "original_position": position,
                    "rotation": prev_rotation + (delta_rotation * ((j+1) / num_interpolations)),
                    "command": prev_command,
                    "extrusion": extrusion/num_interpolations if extrusion is not None else None,
                    "inv_time_feed": inv_time_feed * num_interpolations if inv_time_feed is not None else None,
                    "extrusion_multiplier": extrusion_multiplier,
                    "feed": gcode_point["feed"],
                    "travelling": prev_travelling
                })
        else:
            new_gcode_points.append({
                "position": new_position,
                "original_position": position,
                "rotation": rotation,
                "command": command,
                "extrusion": extrusion,
                "inv_time_feed": inv_time_feed,
                "extrusion_multiplier": extrusion_multiplier,
                "feed": gcode_point["feed"],
                "travelling": travelling
            })

        prev_rotation = rotation
        prev_new_position = new_position.copy()
        prev_travelling = travelling
        prev_command = command

        if command == "G01" and extrusion is not None and extrusion > 0 and (highest_printed_point != 0 or new_position[2] < 1):
            highest_printed_point = max(highest_printed_point, new_position[2])

    print(f"Lost {len(lost_vertices)} vertices")
    
def SAVE_GCODE():
    global gcode_points
    global new_gcode_points
    
    prev_r = 0
    prev_theta = 0
    prev_z = 20
    theta_accum = 0

    # save transformed gcode
    with open(f'{OUTPUT_GCODE_DIR}/.gcode', 'w') as fh:
        # write header
        fh.write("G94 ; mm/min feed  \n")
        fh.write("G28 ; home \n")
        fh.write("M83 ; relative extrusion \n")
        fh.write("G1 E10 ; prime extruder \n")
        fh.write("G94 ; mm/min feed \n")
        fh.write("G90 ; absolute positioning \n")
        fh.write(f"G0 C{prev_theta} X{prev_r} Z{prev_z} B0 ; go to start \n")
        fh.write("G93 ; inverse time feed \n")

        for i, point in enumerate(new_gcode_points):
            position = point["position"]
            rotation = point["rotation"]

            if np.all(np.isnan(position)):
                continue

            if position[2] < 0:
                continue

            z_hop = 0
            if point["travelling"]:
                z_hop = 1

            # convert to polar coordinates
            r = np.linalg.norm(position[:2])
            theta = np.arctan2(position[1], position[0])
            z = position[2]

            # compensate for nozzle offset
            r += -np.sin(rotation) * (NOZZLE_OFFSET + z_hop)
            z += (np.cos(rotation) - 1) * (NOZZLE_OFFSET + z_hop) + z_hop

            delta_theta = theta - prev_theta
            if delta_theta > np.pi:
                delta_theta -= 2*np.pi
            if delta_theta < -np.pi:
                delta_theta += 2*np.pi

            theta_accum += delta_theta

            string = f"{point['command']} Y{np.rad2deg(theta_accum):.5f} X{r:.5f} Z{z:.5f} B{np.rad2deg(rotation):.5f}" # polar printer
            # string = f"{point['command']} X{position[0]:.5f} Y{position[1]:.5f} Z{position[2]} B{np.rad2deg(rotation):.5f}" # cartesian printer (3 axis)

            if point["extrusion"] is not None:
                string += f" E{point['extrusion']:.4f}"

            no_feed_value = False
            if point["inv_time_feed"] is not None:
                string += f" F{(point['inv_time_feed']):.4f}"
            else:
                string += f" F20000"
                fh.write(f"G94\n")
                no_feed_value = True

            fh.write(string + "\n")

            if no_feed_value:
                fh.write(f"G93\n") # back to inv feed

            # update previous values
            prev_r = r
            prev_theta = theta
            prev_z = z
            
def PREV_TRAJECTORY_DOTS():
    global gcode_points
    global new_gcode_points
    global temp

    # 1. Initialisation des listes pour tout collecter en un seul passage
    positions = []
    rotations = []
    travellings = []
    commands = []
    feeds = []
    original_zs = []
    extrusion_multipliers = []

    total_points = len(new_gcode_points)

    for i, point in enumerate(new_gcode_points):
        # Envoi de la télémétrie toutes les X itérations (ex: toutes les 500) pour ne pas saturer le flux
        if i % 3 == 0 or i == total_points - 1:
            progress = {"type": "progress", "task": "gcode_analysis", "current": i, "total": total_points}
            sys.stdout.write(json.dumps(progress) + "\n")
            sys.stdout.flush()

        positions.append(point["position"])
        rotations.append(np.rad2deg(point["rotation"]))
        travellings.append(point["travelling"])
        commands.append(point["command"])
        feeds.append(min(point["feed"], 10000))
        original_zs.append(point["original_position"][2])
        
        em = point["extrusion_multiplier"]
        extrusion_multipliers.append(min(em, 15) if em is not None else np.nan)

    # 3. Conversion de notre structure en PolyData PyVista
    temp = pv.PolyData(np.array(positions))

    temp.cell_data["rotation"] = np.array(rotations)
    temp.cell_data["travelling"] = np.array(travellings)
    temp.cell_data["command"] = np.array(commands)
    temp.cell_data["feed"] = np.array(feeds)
    temp.cell_data["original_z"] = np.array(original_zs)
    temp.cell_data["original_z_bands"] = temp.cell_data["original_z"] % 1
    temp.cell_data["extrusion_multiplier"] = np.array(extrusion_multipliers)

    # 4. Calcul vectoriel ultra-rapide pour le delta_rotation (hors de la boucle)
    if len(new_gcode_points) > 1:
        # np.diff calcule directement la différence entre éléments consécutifs
        delta = np.diff(rotations)
        # On insère la valeur initiale '1' au début et on applique le clip
        temp.cell_data["delta_rotation"] = np.clip(np.insert(delta, 0, 1), -10, 10)
    else:
        temp.cell_data["delta_rotation"] = np.array([1])

    # 5. Extraction et affichage du résultat
    mask = temp.cell_data["command"] == "G01"
    temp.extract_cells(mask).plot(scalars="original_z_bands", cpos=[-0.5, -1, 0.5], point_size=10)
    
def PREV_TRAJECTORY_PATH():
    import plotly.graph_objects as go
    global temp

    fig = go.Figure()

    # On ne garde que G01, adieu les G00 !
    commands_config = [
        {"cmd": "G01", "color": "blue", "name": "G01 (Travail)"}
    ]

    total_configs = len(commands_config) # Le total pour ton interface
        
    for i, config in enumerate(commands_config):
        
        # 1. Envoi de l'état de progression au format JSON
        progress = {
            "type": "progress", 
            "task": "Extraction et injection des trajectoires G-code", 
            "current": i, 
            "total": total_configs
        }
        sys.stdout.write(json.dumps(progress) + "\n")
        sys.stdout.flush()

        # 2. Ton calcul
        indices = np.where(temp.cell_data["command"] == config["cmd"])[0]
        
        if len(indices) > 0:
            fig.add_trace(go.Scatter3d(
                x=temp.points[indices, 0],
                y=temp.points[indices, 1],
                z=temp.points[indices, 2],
                mode='lines+markers',
                marker=dict(size=2, color=config["color"], opacity=0.5),
                line=dict(color=config["color"], width=1),
                name=config["name"]
            ))
                
    # Ajustement des proportions (aspect ratio)
    max_range = np.array([np.ptp(temp.points[:,0]), np.ptp(temp.points[:,1]), np.ptp(temp.points[:,2])])
    aspect_ratio = max_range / np.max(max_range)

    fig.update_layout(
        scene=dict(
            aspectmode='manual',
            aspectratio=dict(x=aspect_ratio[0], y=aspect_ratio[1], z=aspect_ratio[2]),
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        ),
        margin=dict(r=0, l=0, b=0, t=0),
        width=800,
        height=800
    )

    fig.show(renderer="browser")