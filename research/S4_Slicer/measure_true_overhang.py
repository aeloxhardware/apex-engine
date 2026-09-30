"""
Mesure du surplomb REEL d'un modele traite par la methode Bird, dans
l'espace NON-DEFORME -- c'est-a-dire sur la piece telle qu'elle sera
reellement imprimee.

POURQUOI CE SCRIPT EXISTE (erreur de mesure commise puis corrigee) :
measure_real_overhang.py mesurait l'angle des faces du mesh DEFORME. C'est
faux : le mesh deforme n'est jamais imprime, c'est un espace de travail
intermediaire. Une zone plate de la piece reelle (ex: la base circulaire au
contact plateau) apparait courbee -- donc "en surplomb" -- dans l'espace
deforme, uniquement parce que la deformation d'une zone voisine l'a
entrainee par continuite. Apres reprojection (deform_back), elle redevient
plate. Mesurer sur le mesh deforme revient donc a mesurer l'AMPLEUR DE LA
DEFORMATION, pas le surplomb imprimable.

CE QUE MESURE CE SCRIPT :
pour chaque face de la surface d'ORIGINE, l'angle entre sa normale et la
verticale APRES compensation par l'inclinaison de tete que la machine
appliquera a cet endroit. Une face dont la normale pointe vers le bas dans
le repere piece n'est PAS en surplomb si la tete est inclinee de facon a
deposer perpendiculairement a elle.

La rotation par tetraedre est MESUREE sur la deformation (algorithme de
Kabsch entre sommets d'origine et sommets deformes, projetes dans le plan
radial), pas reprise du rotation_field de l'optimisation : c'est la
transcription exacte de ce que fait la cellule 17 du notebook original au
moment de reprojeter le G-code. Elle est ensuite bornee a la plage
physique de la tete [MIN_ROTATION, MAX_ROTATION], comme dans le notebook.

Usage:
  python measure_true_overhang.py <modele> [seuil_apex_deg]
  (<modele> = nom dans input_models/, sans .stl)
"""
import sys

import numpy as np

import bird_core

# Plage physique de la tete, valeurs du notebook original (cellule 17).
# Bird : MAX_ROTATION = 30, MIN_ROTATION = -130 degres.
MAX_ROTATION_DEG = 30.0
MIN_ROTATION_DEG = -130.0

model_name = sys.argv[1] if len(sys.argv) > 1 else "real_overhang"
threshold_apex_deg = float(sys.argv[2]) if len(sys.argv) > 2 else 45.0
# Deux leviers DISTINCTS, a ne pas confondre (c'est tout l'objet de la mesure) :
#   arg 3 -- plafond de volume par tetraedre (switch `-a` de TetGen). Raffine
#            le CHAMP DE ROTATION, defini par tetraedre.
#   arg 4 -- passes de subdivision des triangles d'ENTREE. Raffine la MESURE,
#            faite par face de surface. TetGen ne subdivise jamais la peau.
max_tet_volume = float(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != "-" else None
surface_subdivisions = int(sys.argv[4]) if len(sys.argv) > 4 else 0

# Le barème de continuation valide en session : ROTATION_MULTIPLIER=1 fixe,
# MAX_OVERHANG decroissant, chaque passe repartant du mesh deja deforme.
BAREME = [30, 20, 10, 5]


def project_point_onto_plane(plane_x_axis, plane_y_axis, points):
    """Transcription de la fonction du notebook (cellule 15)."""
    projected_x = np.sum(plane_x_axis * points, axis=1)
    projected_y = np.sum(plane_y_axis * points, axis=1)
    return np.array([projected_x, projected_y]).T


def measure_cell_rotations(input_tet, deformed_tet):
    """
    Rotation effective par tetraedre, MESUREE sur la deformation.

    Transcription de la cellule 17 du notebook : Kabsch entre les sommets
    d'origine et les sommets deformes, tous deux centres puis projetes dans
    le plan radial (le plan contenant l'axe Z et la direction radiale du
    tetraedre) -- c'est le plan dans lequel la tete de la machine peut
    physiquement s'incliner.
    """
    n_cells = deformed_tet.number_of_cells
    cell_rotations = np.zeros(n_cells)

    cells = deformed_tet.field_data["cells"]
    new_verts_all = deformed_tet.field_data["cell_vertices"]
    old_verts_all = input_tet.field_data["cell_vertices"]
    new_centers = deformed_tet.cell_data["cell_center"]
    old_centers = input_tet.cell_data["cell_center"]

    plane_y_vector = np.array([0.0, 0.0, 1.0])

    for cell_index, cell in enumerate(cells):
        new_vertices = new_verts_all[cell] - new_centers[cell_index]
        old_vertices = old_verts_all[cell] - old_centers[cell_index]

        old_center_xy = old_centers[cell_index][:2]
        norm_xy = np.linalg.norm(old_center_xy)
        if norm_xy < 1e-12:
            # tetraedre sur l'axe : direction radiale indeterminee, meme
            # singularite que R=0 documentee en limite #19 cote Apex.
            plane_x_vector = np.array([1.0, 0.0, 0.0])
        else:
            pxy = old_center_xy / norm_xy
            plane_x_vector = np.array([pxy[0], pxy[1], 0.0])

        new_proj = project_point_onto_plane(plane_x_vector, plane_y_vector, new_vertices)
        old_proj = project_point_onto_plane(plane_x_vector, plane_y_vector, old_vertices)

        covariance_matrix = np.dot(new_proj.T, old_proj)
        U, _, Vt = np.linalg.svd(covariance_matrix)
        rotation_matrix = np.dot(U, Vt)

        rotation = -np.arccos(min(max(rotation_matrix[0, 0], -1), 1))
        if rotation_matrix[1, 0] < 0:
            rotation = -rotation

        # bornage a la plage physique de la tete, comme le notebook
        rotation = max(min(rotation, np.deg2rad(MAX_ROTATION_DEG)), np.deg2rad(MIN_ROTATION_DEG))
        cell_rotations[cell_index] = rotation

    return cell_rotations


def face_rotations_from_cells(tet, cell_rotations):
    """
    Rotation de tete applicable a chaque face de la surface d'origine.

    cell_to_face donne, pour chaque tetraedre de bord, les faces de surface
    qu'il porte -- on remonte donc de la cellule vers ses faces (meme table
    que bird_core utilise pour calculer face_normal / overhang_angle).
    """
    cell_to_face = bird_core.decode_object(tet.field_data["cell_to_face"])
    n_faces = len(tet.field_data["faces"])

    face_rot_sum = np.zeros(n_faces)
    face_rot_count = np.zeros(n_faces)
    for cell_index, face_indices in cell_to_face.items():
        for face_index in face_indices:
            face_rot_sum[face_index] += cell_rotations[cell_index]
            face_rot_count[face_index] += 1

    with np.errstate(invalid="ignore", divide="ignore"):
        face_rot = np.where(face_rot_count > 0, face_rot_sum / np.maximum(face_rot_count, 1), np.nan)
    return face_rot, face_rot_count


def rotate_about_tangential(normals, centers, angles):
    """
    Applique a chaque normale la rotation d'angle `angle` autour de l'axe
    tangentiel local (perpendiculaire au rayon et a Z) -- la MEME convention
    d'axe que calculate_rotation_matrices dans bird_core, pour que la
    compensation soit coherente avec la deformation qui l'a produite.
    """
    xy = centers[:, :2]
    xy_3d = np.hstack([xy, np.zeros((xy.shape[0], 1))])
    tangential = np.cross(np.array([0.0, 0.0, 1.0]), xy_3d)
    norms = np.linalg.norm(tangential, axis=1)
    safe = norms > 1e-12
    tangential[safe] /= norms[safe][:, None]
    tangential[~safe] = [1.0, 0.0, 0.0]

    # Rodrigues, vectorise
    k = tangential
    cos_a = np.cos(angles)[:, None]
    sin_a = np.sin(angles)[:, None]
    k_cross_n = np.cross(k, normals)
    k_dot_n = np.sum(k * normals, axis=1)[:, None]
    return normals * cos_a + k_cross_n * sin_a + k * k_dot_n * (1 - cos_a)


def report(label, normals, areas, threshold_from_vertical, total_area, exclude=None):
    up = np.array([0.0, 0.0, 1.0])
    cos_angle = np.clip(normals @ up, -1, 1)
    angle_deg = np.degrees(np.arccos(cos_angle))
    is_overhang = angle_deg > threshold_from_vertical
    if exclude is not None:
        is_overhang = is_overhang & ~exclude
    area_overhang = areas[is_overhang].sum()
    print(f"  {label:<46} {is_overhang.sum():>6} faces  {area_overhang / total_area * 100:>6.2f} % de la surface")
    return is_overhang


def platform_contact_mask(centers, tolerance=1e-3):
    """
    Faces au contact du plateau : elles ont une normale vers le bas (donc
    "en surplomb" au sens de l'angle) mais reposent sur le plateau, jamais
    dans le vide. Meme exclusion que step3_segmentation.py cote Apex, qui
    force leur risque a zero (limite #1).
    """
    z_min = centers[:, 2].min()
    return centers[:, 2] < z_min + tolerance


def main():
    threshold_from_vertical = 90 + threshold_apex_deg

    print(f"Modele          : {model_name}")
    print(f"Seuil Apex      : {threshold_apex_deg} deg  (angle(normale, +Z) > {threshold_from_vertical} deg)")
    print(f"Plage tete      : [{MIN_ROTATION_DEG}, {MAX_ROTATION_DEG}] deg (valeurs du notebook Bird)")
    print(f"Bareme          : MAX_OVERHANG {BAREME}, ROTATION_MULTIPLIER=1")
    print(f"Densite tetra   : max_tet_volume={max_tet_volume}  (None = defaut TetGen)")
    print(f"Subdiv surface  : {surface_subdivisions} passe(s)  (0 = peau d'origine)")
    print()

    input_tet = bird_core.load_and_tetrahedralize(
        f"input_models/{model_name}.stl",
        max_tet_volume=max_tet_volume,
        surface_subdivisions=surface_subdivisions,
    )
    input_tet, state = bird_core.calculate_tet_attributes(input_tet)
    original_tet = input_tet.copy(deep=True)
    print(f"tetraedres : {input_tet.number_of_cells}, faces de surface : {len(input_tet.field_data['faces'])}")
    print()

    tet = input_tet
    for i, max_overhang in enumerate(BAREME):
        params = {
            "NEIGHBOUR_LOSS_WEIGHT": 20,
            "MAX_OVERHANG": max_overhang,
            "ROTATION_MULTIPLIER": 1,
            "SET_INITIAL_ROTATION_TO_ZERO": False,
            "INITIAL_ROTATION_FIELD_SMOOTHING": 30,
            "MAX_POS_ROTATION_DEG": 360,
            "MAX_NEG_ROTATION_DEG": -360,
            "ROTATE_ITERATIONS": 100,
            "DEFORM_ITERATIONS": 1000,
            "STEEP_OVERHANG_COMPENSATION": True,
        }
        tet, info = bird_core.run_one_iteration(tet, state, params)
        print(f"  passe {i + 1}/{len(BAREME)} (MAX_OVERHANG={max_overhang:>2}) : "
              f"deform_cost={info['deform_cost']:.3f}  overhang_fraction_interne={info['overhang_fraction'] * 100:.2f}%")
    deformed_tet = tet
    print()

    # --- rotations effectives, mesurees sur la deformation ------------------
    cell_rotations = measure_cell_rotations(original_tet, deformed_tet)
    face_rot, face_rot_count = face_rotations_from_cells(original_tet, cell_rotations)

    covered = face_rot_count > 0
    print(f"faces de surface portees par un tetraedre de bord : {covered.sum()} / {len(covered)}")
    rot_deg = np.degrees(cell_rotations)
    print(f"rotation de tete mesuree (par tetraedre) : "
          f"min={rot_deg.min():.1f} deg  max={rot_deg.max():.1f} deg  mediane={np.median(rot_deg):.1f} deg")
    print()

    # --- comparaison des trois mesures --------------------------------------
    surface = original_tet.extract_surface()
    normals_orig = surface.face_normals
    areas = surface.compute_cell_sizes(length=False, area=True, volume=False)["Area"]
    centers_orig = surface.cell_centers().points
    total_area = areas.sum()

    surface_def = deformed_tet.extract_surface()
    normals_def = surface_def.face_normals
    areas_def = surface_def.compute_cell_sizes(length=False, area=True, volume=False)["Area"]

    on_platform = platform_contact_mask(centers_orig)
    print(f"faces exclues comme contact plateau : {on_platform.sum()}")
    print()

    angles = np.where(covered, face_rot, 0.0)
    normals_compensated = rotate_about_tangential(normals_orig, centers_orig, angles)
    norms = np.linalg.norm(normals_compensated, axis=1)
    normals_compensated /= np.maximum(norms, 1e-12)[:, None]

    print("MESURES (contact plateau exclu)")
    before = report("1. piece d'origine, sans compensation", normals_orig, areas,
                    threshold_from_vertical, total_area, exclude=on_platform)
    report("2. mesh DEFORME (ce que je mesurais a tort)", normals_def, areas_def,
           threshold_from_vertical, areas_def.sum())
    after = report("3. piece reelle AVEC inclinaison de tete", normals_compensated, areas,
                   threshold_from_vertical, total_area, exclude=on_platform)

    # --- detail par face, pour voir QUELLES faces restent en surplomb -------
    up = np.array([0.0, 0.0, 1.0])
    ang_before = np.degrees(np.arccos(np.clip(normals_orig @ up, -1, 1)))
    ang_after = np.degrees(np.arccos(np.clip(normals_compensated @ up, -1, 1)))
    concerned = np.where(before | after)[0]
    if len(concerned):
        print()
        print("DETAIL des faces en surplomb (avant et/ou apres compensation)")
        print(f"  {'face':>5} {'z':>7} {'aire':>7} {'angle avant':>12} {'rot tete':>9} {'angle apres':>12}  statut")
        for i in concerned:
            statut = "RESOLU" if (before[i] and not after[i]) else ("reste" if after[i] else "-")
            print(f"  {i:>5} {centers_orig[i][2]:>7.1f} {areas[i]:>7.2f} "
                  f"{ang_before[i]:>11.1f} deg {np.degrees(angles[i]):>8.1f} deg {ang_after[i]:>11.1f} deg  {statut}")


if __name__ == "__main__":
    main()
