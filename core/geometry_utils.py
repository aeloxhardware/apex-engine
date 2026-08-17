"""
core/geometry_utils.py

Utilitaires geometriques partages entre plusieurs etapes du pipeline.

- build_face_adjacency : graphe dual du mesh (deja duplique localement dans
  step2_overhang.py ; step3_segmentation.py consomme la version partagee ici).
- fibonacci_sphere_directions : repartition uniforme de directions sur la
  sphere complete, pour une machine a orientation libre (5 axes).
- rotation_axis_circle_directions : repartition de directions sur le cercle
  atteignable par une machine 4 axes reelle (3 lineaires + 1 seul DOF de
  rotation) - PAS une sphere complete, voir CLAUDE.md limite de conception #6.
- point_in_polygon_2d : test d'imbrication geometrique standard
  (ray-casting), utilise pour distinguer un contour exterieur d'un trou
  (step5 pour la verification de coherence, step6 pour la classification
  des contours avant remplissage) - PAS le signe de l'aire projetee, dont
  l'hypothese s'est revelee fausse empiriquement (voir CLAUDE.md limite #11).
- composed_rotation_directions : repartition de directions atteignables par
  la composition ordonnee de N axes de rotation (core.mesh_types.RotationAxisSpec),
  chacun avec sa propre plage angulaire - generalisation qui englobe
  rotation_axis_circle_directions (cas a 1 axe, plage complete) SANS le
  remplacer, pour decrire une cinematique de machine arbitraire (pas figee
  sur une machine en particulier - voir CLAUDE.md limite de conception #18).
- composed_rotation_angles_and_directions : meme grille que ci-dessus, mais
  retourne aussi les angles par axe associes a chaque direction - requis par
  step8 (resolution cinematique inverse), qui doit repartir d'angles
  concrets, pas seulement de la direction resultante.
- refine_axis_angles_local_search : affinement generique par coordinate-
  descent d'un jeu d'angles vers une direction cible, utilise par step8
  apres une recherche de plus proche voisin sur composed_rotation_angles_and_directions.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def build_face_adjacency(faces: np.ndarray, face_count: int) -> sp.csr_matrix:
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


def fibonacci_sphere_directions(count: int) -> np.ndarray:
    """
    Repartit `count` directions unitaires quasi-uniformement sur la sphere
    complete, via une spirale de Fibonacci. Utilise quand
    SegmentationParams.rotation_axis_direction est None (orientation libre,
    type 5 axes).
    """
    if count <= 0:
        return np.zeros((0, 3))

    indices = np.arange(count)
    golden_ratio = (1.0 + np.sqrt(5.0)) / 2.0

    # z reparti uniformement sur [-1, 1], theta = angle azimutal
    z = 1.0 - 2.0 * (indices + 0.5) / count
    radius = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    theta = 2.0 * np.pi * indices / golden_ratio

    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    return np.stack([x, y, z], axis=1)


def rotation_axis_circle_directions(
    rotation_axis_direction: np.ndarray, count: int
) -> np.ndarray:
    """
    Repartit `count` directions unitaires uniformement sur le cercle des
    directions atteignables par une machine 4 axes reelle : 3 axes lineaires
    + un seul DOF de rotation autour de `rotation_axis_direction`. Un seul
    DOF de rotation ne balaie qu'un cercle perpendiculaire a l'axe de
    rotation, jamais la sphere complete (reservee a rotation_axis_direction
    is None -> fibonacci_sphere_directions).
    """
    if count <= 0:
        return np.zeros((0, 3))

    axis = np.asarray(rotation_axis_direction, dtype=float)
    axis = axis / np.linalg.norm(axis)

    # construit une base orthonormee (u, v) du plan perpendiculaire a axis
    reference = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(reference, axis)) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])

    u = reference - np.dot(reference, axis) * axis
    u = u / np.linalg.norm(u)
    v = np.cross(axis, u)

    angles = 2.0 * np.pi * np.arange(count) / count
    directions = (
        np.cos(angles)[:, None] * u[None, :] + np.sin(angles)[:, None] * v[None, :]
    )

    return directions


def _rodrigues_rotate(vectors: np.ndarray, axis: np.ndarray, angles_rad: np.ndarray) -> np.ndarray:
    """
    Tourne chaque ligne de `vectors` (N,3) autour de `axis` (3,, unitaire)
    d'un angle correspondant dans `angles_rad` (N,), via la formule de
    rotation de Rodrigues. Vectorise sur N pour eviter une boucle Python
    lors de l'echantillonnage d'une grille produit d'angles.
    """
    cos_a = np.cos(angles_rad)[:, None]
    sin_a = np.sin(angles_rad)[:, None]
    dot = (vectors @ axis)[:, None]
    cross = np.cross(np.broadcast_to(axis, vectors.shape), vectors)
    return vectors * cos_a + cross * sin_a + axis[None, :] * dot * (1.0 - cos_a)


def composed_rotation_angles_and_directions(
    axes: "list", reference_direction: np.ndarray, samples_per_axis: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Meme grille produit et meme composition que composed_rotation_directions,
    mais retourne AUSSI les angles par axe (degres, shape (P, len(axes)))
    associes a chaque direction (shape (P, 3)) de l'echantillonnage.

    composed_rotation_directions ne retournait que les directions,
    suffisant pour step3/step7 (qui ne testent qu'une compatibilite
    d'angle), mais pas pour step8 : une resolution cinematique inverse doit
    repartir d'angles CONCRETS une fois le plus proche voisin trouve, pas
    seulement de la direction resultante. Factorise ici pour ne pas
    dupliquer la grille produit et la composition Rodrigues.

    `axes[0]` est la rotation la plus "externe" (ex: le plateau C), `axes[-1]`
    la plus "interne" (ex: la tete B) - meme convention que l'ordre
    d'application matricielle.

    Generalisation directe et verifiable de rotation_axis_circle_directions,
    pas un remplacement :
      - axes vide -> `reference_direction` seule (aucune rotation), angles (P,0).
      - un seul axe a plage complete (0,360) -> memes directions que
        rotation_axis_circle_directions (verifie par test).
      - N axes -> grille produit des echantillons de chaque plage, chaque
        point de la grille compose par rotations de Rodrigues successives.

    Ne fait AUCUNE hypothese sur le nombre d'axes ni leurs plages : une
    machine a 2 DOF dont les deux plages sont restreintes ne couvre en
    general PAS la sphere complete, contrairement a une machine dont un axe
    a plage complete et l'autre une demi-plage (cas mesure separement, pas
    suppose ici).
    """
    reference_direction = np.asarray(reference_direction, dtype=float)
    reference_direction = reference_direction / np.linalg.norm(reference_direction)

    if not axes:
        return reference_direction.reshape(1, 3), np.zeros((1, 0))

    # grille produit des angles echantillonnes par axe, un DOF fige
    # (min==max) ne produit qu'un seul echantillon plutot que
    # samples_per_axis redondants. Un tour complet (plage >= 360deg) exclut
    # l'extremite superieure (0deg et 360deg sont le meme point, cf.
    # rotation_axis_circle_directions) ; une plage bornee strictement (ex:
    # B in [0,180], deux valeurs physiquement distinctes) inclut les deux
    # bornes - confondre les deux cas produirait un espacement d'angle faux
    # (mesure : linspace(0,360,8,endpoint=True) espace de 360/7 au lieu de
    # 360/8, dupliquant le point 0/360).
    per_axis_angles_deg = []
    per_axis_angles_rad = []
    for spec in axes:
        lo, hi = spec.angle_range_deg
        is_full_turn = (hi - lo) >= 360.0
        n = 1 if hi <= lo else max(samples_per_axis, 1)
        angles_deg = np.linspace(lo, hi, n, endpoint=not is_full_turn)
        per_axis_angles_deg.append(angles_deg)
        per_axis_angles_rad.append(np.radians(angles_deg))

    grids_deg = np.meshgrid(*per_axis_angles_deg, indexing="ij")
    grids_rad = np.meshgrid(*per_axis_angles_rad, indexing="ij")
    flat_angles_deg = [g.reshape(-1) for g in grids_deg]  # une colonne d'angles (deg) par axe
    flat_angles_rad = [g.reshape(-1) for g in grids_rad]
    total_points = flat_angles_rad[0].shape[0]

    angles_per_point = np.stack(flat_angles_deg, axis=1)  # (P, len(axes))

    directions = np.tile(reference_direction, (total_points, 1))
    # applique les rotations dans l'ordre INVERSE de la liste (la plus
    # interne d'abord, cf. axes[-1] = tete, la plus proche de la buse) pour
    # que la composition matricielle corresponde a R_1 . R_2 ... . R_n . d_ref
    for spec, angles in zip(reversed(axes), reversed(flat_angles_rad)):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        directions = _rodrigues_rotate(directions, axis, angles)

    return directions, angles_per_point


def composed_rotation_directions(
    axes: "list", reference_direction: np.ndarray, samples_per_axis: int
) -> np.ndarray:
    """
    Directions atteignables par la composition ORDONNEE de N rotations
    (`axes`, une liste de `core.mesh_types.RotationAxisSpec`), chacune
    echantillonnee sur sa propre plage angulaire.

    Wrapper mince autour de composed_rotation_angles_and_directions (voir
    sa docstring pour le detail de la generalisation) qui ne retourne que
    les directions - conserve pour les appelants (step3, step7) qui n'ont
    besoin que de tester une compatibilite d'angle, jamais de repartir des
    angles concrets.
    """
    directions, _angles = composed_rotation_angles_and_directions(
        axes, reference_direction, samples_per_axis
    )
    return directions


def refine_axis_angles_local_search(
    axes: "list",
    reference_direction: np.ndarray,
    initial_angles_deg: np.ndarray,
    target_direction: np.ndarray,
    iterations: int,
    initial_step_deg: float,
) -> tuple[np.ndarray, float]:
    """
    Affine `initial_angles_deg` (typiquement le meilleur candidat d'une
    grille composed_rotation_angles_and_directions) pour approcher au plus
    pres `target_direction`, par coordinate-descent generique : a chaque
    iteration, chaque axe est tour a tour nudge de +-`step` (clampe a son
    `angle_range_deg`), le mouvement n'est garde que s'il reduit l'erreur
    angulaire ; `step` est divise par 2 a chaque iteration.

    Seule primitive par-axe utilisee : _rodrigues_rotate (deja existante) -
    generalise a N axes SANS jacobien analytique par-N (qui devrait etre
    re-derive a chaque nombre/configuration d'axes, contraire a l'objectif
    de cinematique arbitraire). Un axe verrouille (min==max) ne produit
    aucun nudge valide et reste immobile - degenerescence geree par la
    meme boucle, pas un cas special.

    IMPORTANT : `initial_step_deg` doit demarrer a l'echelle de l'ecart
    entre deux points voisins de la grille NN d'origine (ex: 360/samples_per_axis),
    PAS une petite constante fixe - sinon l'affinement reste piege pres du
    candidat NN au lieu de converger vers le vrai optimum local (mesure :
    un pas fixe de 2deg laisse un residu de 15-20deg quand l'ecart de grille
    est ~36deg ; un pas derive de l'ecart de grille converge sous 0.1deg).

    Retourne (angles_affines_deg, erreur_angulaire_residuelle_rad).
    """
    reference_direction = np.asarray(reference_direction, dtype=float)
    reference_direction = reference_direction / np.linalg.norm(reference_direction)
    target_direction = np.asarray(target_direction, dtype=float)
    target_norm = np.linalg.norm(target_direction)
    if target_norm > 1e-12:
        target_direction = target_direction / target_norm

    angles = np.array(initial_angles_deg, dtype=float).copy()

    def _direction_for(angles_deg: np.ndarray) -> np.ndarray:
        direction = reference_direction.reshape(1, 3)
        for spec, angle_deg in zip(reversed(axes), reversed(list(angles_deg))):
            axis = np.asarray(spec.axis, dtype=float)
            axis = axis / np.linalg.norm(axis)
            direction = _rodrigues_rotate(direction, axis, np.array([np.radians(angle_deg)]))
        return direction[0]

    def _angular_error(direction: np.ndarray) -> float:
        cos_angle = float(np.clip(np.dot(direction, target_direction), -1.0, 1.0))
        return float(np.arccos(cos_angle))

    if not axes:
        return angles, _angular_error(reference_direction)

    best_error = _angular_error(_direction_for(angles))
    step = float(initial_step_deg)

    for _ in range(max(iterations, 0)):
        for axis_index, spec in enumerate(axes):
            lo, hi = spec.angle_range_deg
            if hi <= lo:
                continue  # axe verrouille : aucun nudge possible
            for delta in (step, -step):
                candidate = angles.copy()
                candidate[axis_index] = float(np.clip(candidate[axis_index] + delta, lo, hi))
                error = _angular_error(_direction_for(candidate))
                if error < best_error:
                    angles = candidate
                    best_error = error
        step *= 0.5

    return angles, best_error


def point_in_polygon_2d(point: np.ndarray, polygon: np.ndarray) -> bool:
    """
    Ray-casting standard : `point` (2,) est-il a l'interieur de `polygon`
    (N,2) ? Technique robuste et independante de toute convention de
    winding, contrairement au signe de l'aire projetee (verifie empiriquement
    faux sur un anneau plat synthetique - CLAUDE.md limite #11).
    """
    x, y = point
    inside = False
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            x_intersect = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < x_intersect:
                inside = not inside
    return inside
