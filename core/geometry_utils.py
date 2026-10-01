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
- point_has_support_along_direction : test "y a-t-il un segment support sous
  ce point, le long d'une direction quelconque" - generalise le test
  vertical historique de step9b_support_check (implicite [0,0,1]) a une
  direction arbitraire, pour verifier qu'un point sans support VERTICAL
  n'est pas simplement un point que la machine peut deposer en inclinant la
  tete. Reutilisee par step5_layers (verification de support d'un cluster
  contre les clusters deja places, direction = up_vector uniquement).
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
  Conservee comme repli (voir resolve_axis_angles_analytic) - jamais
  supprimee, toujours utilisee independamment par step10 (enumeration de
  candidats proches pour la reorientation post-collision).
- resolve_axis_angles_analytic : resolution CINEMATIQUE INVERSE generique
  N-axes par bissection/atan2, remplace la recherche par grille de
  _resolve_axis_angles (step8) - voir CLAUDE.md "Piste chiffree pour
  step8" pour le diagnostic (grille : erreur mediane mesuree 9.26 deg sur
  un modele reel dense, ~15000x plus lent que la fermeture analytique
  mesuree sur le cas 2 axes). Repose sur un invariant geometrique verifie
  numeriquement (angle(direction, axes[0]) independant de l'angle
  d'axes[0] lui-meme, pour tout N) qui ne suffit PAS a determiner une
  solution unique pour N>=3 (comptage de degres de liberte : cible = 2
  contraintes sur S^2, N inconnues => famille de solutions de dimension
  N-2 des que N>=3) - la fonction resout donc la PAIRE EXTERNE (axes[0],
  axes[1]) exactement, apres avoir replie les axes internes (axes[2:]) a
  une valeur fixe (continuite avec le point precedent). Retombe
  identiquement (repliement vide) sur la methode 2-axes deja validee
  quand N=2.
"""

from __future__ import annotations

import math
from typing import Callable, Optional

import numpy as np
import scipy.sparse as sp

from core.mesh_types import RotationAxisSpec


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


def _rodrigues_rotate_single(vector: np.ndarray, axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """
    Meme formule que _rodrigues_rotate, pour UN SEUL vecteur (3,) plutot
    qu'un lot (N,3). Ecrite en scalaire pur (indexation directe, pas de
    np.cross/broadcast_to) : np.cross sur un tableau a une seule ligne
    passe par tout le chemin generique numpy N-dimensionnel (moveaxis,
    normalize_axis_tuple, broadcast_to), un surcout enorme mesure pour un
    produit vectoriel 3x3 trivial (~10x plus lent que ce chemin scalaire,
    profile sur step8 - voir CLAUDE.md). Reservee aux appelants qui
    tournent un vecteur a la fois (refine_axis_angles_local_search,
    step10_adjustment._direction_for_angles) - _rodrigues_rotate reste
    inchangee pour les appelants batch legitimes (grille produit,
    sommets d'un mesh de tete en une fois).
    """
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)
    dot = vector[0] * axis[0] + vector[1] * axis[1] + vector[2] * axis[2]
    cross_x = axis[1] * vector[2] - axis[2] * vector[1]
    cross_y = axis[2] * vector[0] - axis[0] * vector[2]
    cross_z = axis[0] * vector[1] - axis[1] * vector[0]
    factor = dot * (1.0 - cos_a)
    return np.array([
        vector[0] * cos_a + cross_x * sin_a + axis[0] * factor,
        vector[1] * cos_a + cross_y * sin_a + axis[1] * factor,
        vector[2] * cos_a + cross_z * sin_a + axis[2] * factor,
    ])


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

    # Axes unitaires precalcules UNE FOIS (pas a chaque nudge - spec.axis ne
    # change jamais pendant la recherche) : evite de refaire N normalisations
    # par appel a _direction_for, redondance pure, pas seulement un surcout
    # de dispatch numpy (voir CLAUDE.md, perf step8).
    reversed_axes_unit = [
        np.asarray(spec.axis, dtype=float) / np.linalg.norm(spec.axis) for spec in reversed(axes)
    ]
    tx, ty, tz = float(target_direction[0]), float(target_direction[1]), float(target_direction[2])

    def _direction_for(angles_deg: np.ndarray) -> np.ndarray:
        direction = reference_direction
        for axis, angle_deg in zip(reversed_axes_unit, reversed(list(angles_deg))):
            direction = _rodrigues_rotate_single(direction, axis, np.radians(angle_deg))
        return direction

    def _angular_error(direction: np.ndarray) -> float:
        # dot/clip/arccos scalaires purs : np.dot/np.clip sur des scalaires
        # passent par le dispatch numpy generique, un surcout mesure comme
        # dominant sur ~700000 appels (CLAUDE.md, perf step8) pour une
        # operation triviale.
        cos_angle = direction[0] * tx + direction[1] * ty + direction[2] * tz
        if cos_angle > 1.0:
            cos_angle = 1.0
        elif cos_angle < -1.0:
            cos_angle = -1.0
        return math.acos(cos_angle)

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
                nudged = candidate[axis_index] + delta
                if nudged > hi:
                    nudged = hi
                elif nudged < lo:
                    nudged = lo
                candidate[axis_index] = nudged
                error = _angular_error(_direction_for(candidate))
                if error < best_error:
                    angles = candidate
                    best_error = error
        step *= 0.5

    return angles, best_error


# ---------------------------------------------------------------------------
# Cinematique inverse ANALYTIQUE generique N-axes (remplace la recherche par
# grille de step8 pour le cas courant) - voir docstring de module et
# CLAUDE.md "Piste chiffree pour step8" pour le diagnostic complet.
# ---------------------------------------------------------------------------

# En-deca, un vecteur est considere (anti)parallele a un axe : sa composante
# perpendiculaire est trop petite pour definir un angle stable par atan2
# (singularite deja documentee, CLAUDE.md limite #19 - "R=0 => C indetermine").
_PARALLEL_EPS = 1e-9
# En-deca, deux axes sont consideres paralleles entre eux (1 seul DOF
# effectif pour la paire) plutot que quasi-paralleles avec un atan2 instable.
_AXES_PARALLEL_EPS = 1e-9
# 1e-5 deg reste tres largement sous l'erreur de la grille remplacee
# (mediane mesuree 9.26 deg) tout en coupant nettement le nombre
# d'iterations de bissection necessaires (chaque iteration divise
# l'intervalle par 2 : 1e-7deg sur une plage de 360deg demande ~32
# iterations, 1e-5deg ~25 - la difference est marginale par iteration mais
# multipliee par le nombre de points de trajectoire).
_BISECTION_TOLERANCE_DEG = 1e-5
_BISECTION_MAX_ITERATIONS = 60
# Nombre de points de balayage grossier pour BRACKETER une racine de g(t1)
# avant de bissecter - rien ne garantit la monotonie de g sur un tour
# complet (axes non orthogonaux), un simple test aux deux bornes suffit
# pour le cas orthogonal mais pas en general (voir plan/CLAUDE.md).
_BRACKETING_SAMPLES = 64


def _cross3(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """np.cross scalaire pur (3,)x(3,) - meme raison que _rodrigues_rotate_single
    (limite #37, CLAUDE.md) : np.cross generique a un surcout mesure enorme
    sur un produit vectoriel 3x3 trivial appele des dizaines de milliers de
    fois dans la bissection de _resolve_axis_pair_analytic."""
    return np.array([
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ])


def _dot3(a: np.ndarray, b: np.ndarray) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm3(a: np.ndarray) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def _perpendicular_basis(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Base orthonormee (u, v) du plan perpendiculaire a `axis` (unitaire) -
    meme construction que rotation_axis_circle_directions (reference fixe
    [1,0,0], repli [0,1,0] si quasi-colineaire), pour rester coherente avec
    le reste du fichier plutot que d'inventer une troisieme convention.
    """
    reference = np.array([1.0, 0.0, 0.0])
    if abs(_dot3(reference, axis)) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    u = reference - _dot3(reference, axis) * axis
    u = u / _norm3(u)
    v = _cross3(axis, u)
    return u, v


def _signed_angle_in_plane(
    vector_from: np.ndarray, vector_to: np.ndarray, axis: np.ndarray
) -> Optional[float]:
    """
    Angle signe (degres) de `vector_from` vers `vector_to`, mesure dans le
    plan perpendiculaire a `axis` - None si l'une des projections est quasi
    nulle (vecteur (anti)parallele a `axis`, angle indetermine par
    construction, meme singularite que #19).
    """
    u, v = _perpendicular_basis(axis)
    from_perp = vector_from - _dot3(vector_from, axis) * axis
    to_perp = vector_to - _dot3(vector_to, axis) * axis
    from_norm = _norm3(from_perp)
    to_norm = _norm3(to_perp)
    if from_norm < _PARALLEL_EPS or to_norm < _PARALLEL_EPS:
        return None
    fu, fv = _dot3(from_perp, u), _dot3(from_perp, v)
    tu, tv = _dot3(to_perp, u), _dot3(to_perp, v)
    return float(np.degrees(np.arctan2(tv * fu - tu * fv, tu * fu + tv * fv)))


def _clamp_into_range(angle_deg: float, lo: float, hi: float) -> float:
    """
    Ramene `angle_deg` dans `[lo, hi]`. Un tour complet (plage >= 360deg)
    est TOUJOURS exactement atteignable : replie par modulo. Une plage
    bornee strictement clampe a la borne la plus proche au sens CIRCULAIRE
    (ex: 350deg est plus proche de lo=0 que de hi=20 par enroulement, pas
    au sens lineaire naif) - meme piege que le depli de step8
    (_unwrap_angle_sequence), evite ici a la source.
    """
    if hi - lo >= 360.0:
        return lo + ((angle_deg - lo) % 360.0)
    if lo <= angle_deg <= hi:
        return angle_deg
    dist_lo = min(abs(angle_deg - lo), 360.0 - abs(angle_deg - lo))
    dist_hi = min(abs(angle_deg - hi), 360.0 - abs(angle_deg - hi))
    return lo if dist_lo <= dist_hi else hi


def _forward_direction(
    axes: list, angles_deg: list, reference_direction: np.ndarray
) -> np.ndarray:
    """Cinematique DIRECTE : meme ordre de composition (le plus interne
    d'abord) que composed_rotation_angles_and_directions/step8, via la
    primitive scalaire (un seul point a la fois, pas une grille)."""
    direction = reference_direction
    for spec, angle_deg in zip(reversed(axes), reversed(angles_deg)):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        direction = _rodrigues_rotate_single(direction, axis, np.radians(angle_deg))
    return direction


def _angle_between(a: np.ndarray, b: np.ndarray) -> float:
    """
    Angle (degres) entre deux vecteurs UNITAIRES. Scalaire Python pur
    (math.acos, pas np.dot/np.clip) - meme correction que
    refine_axis_angles_local_search._angular_error (CLAUDE.md, perf step8) :
    np.clip/np.dot sur des scalaires triviaux passe par tout le dispatch
    numpy generique, un surcout mesure comme dominant sur cette fonction
    appelee ~90 fois par point de trajectoire (bracketing + bissection).
    """
    cos_angle = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    if cos_angle > 1.0:
        cos_angle = 1.0
    elif cos_angle < -1.0:
        cos_angle = -1.0
    return math.degrees(math.acos(cos_angle))


def _resolve_single_axis_analytic(
    axis_spec: RotationAxisSpec,
    reference: np.ndarray,
    target: np.ndarray,
    fallback_angle_deg: float,
) -> float:
    """
    Resout l'UNIQUE angle d'un axe pour amener `reference` le plus pres
    possible de `target` (cas de base N=1, et sous-probleme interne de la
    paire externe pour N>=2). Verrouille (lo==hi) : retourne la valeur
    verrouillee sans chercher. Singularite (`reference` ou `target`
    (anti)parallele a l'axe) : retourne `fallback_angle_deg` (l'angle
    precedent, coherent avec l'heritage de singularite deja fait en aval
    par step8 - ne PAS dupliquer cette logique ici, juste ne jamais
    produire de NaN).
    """
    lo, hi = axis_spec.angle_range_deg
    if hi <= lo:
        return lo

    axis = np.asarray(axis_spec.axis, dtype=float)
    axis = axis / np.linalg.norm(axis)

    delta = _signed_angle_in_plane(reference, target, axis)
    if delta is None:
        return _clamp_into_range(fallback_angle_deg, lo, hi)

    return _clamp_into_range(delta, lo, hi)


def _resolve_axis_pair_analytic(
    outer_spec: RotationAxisSpec,
    inner_spec: RotationAxisSpec,
    reference_effective: np.ndarray,
    target: np.ndarray,
    fallback_angles_deg: tuple[float, float],
) -> Optional[tuple[float, float]]:
    """
    Resout la PAIRE (outer_spec, inner_spec) pour amener
    `reference_effective` sur `target`, par bissection sur l'angle interne
    puis atan2 pour l'angle externe (voir CLAUDE.md/plan pour la
    justification geometrique complete de l'invariant utilise).

    Retourne None si aucune racine n'a pu etre bracketee sur la plage de
    l'axe interne - CA ARRIVE REELLEMENT (pas juste en theorie, verifie
    numeriquement pendant la conception) quand la paire externe est
    geometriquement incapable d'atteindre l'angle cible par rapport a
    outer_spec, meme en balayant tout l'axe interne (l'ensemble atteignable
    de angle(.,outer_axis) est un ANNEAU borne par les angles entre les
    deux axes et la reference effective, pas [0,180] en general) - dans ce
    cas l'appelant doit retomber sur le repli grille+raffinement.
    """
    outer_axis = np.asarray(outer_spec.axis, dtype=float)
    outer_axis = outer_axis / np.linalg.norm(outer_axis)
    inner_axis = np.asarray(inner_spec.axis, dtype=float)
    inner_axis = inner_axis / np.linalg.norm(inner_axis)

    outer_lo, outer_hi = outer_spec.angle_range_deg
    inner_lo, inner_hi = inner_spec.angle_range_deg

    outer_locked = outer_hi <= outer_lo
    inner_locked = inner_hi <= inner_lo

    if outer_locked and inner_locked:
        return outer_lo, inner_lo

    if inner_locked:
        # un seul DOF libre (outer) : reduit au cas N=1 contre une
        # reference deja tournee par l'angle verrouille de l'axe interne.
        rotated_reference = _rodrigues_rotate_single(
            reference_effective, inner_axis, np.radians(inner_lo)
        )
        angle = _resolve_single_axis_analytic(
            outer_spec, rotated_reference, target, fallback_angles_deg[0]
        )
        return angle, inner_lo

    if outer_locked:
        # un seul DOF libre (inner) : la cible doit d'abord etre "devissee"
        # de la rotation externe verrouillee avant de resoudre l'interne.
        unwound_target = _rodrigues_rotate_single(target, outer_axis, np.radians(-outer_lo))
        angle = _resolve_single_axis_analytic(
            inner_spec, reference_effective, unwound_target, fallback_angles_deg[1]
        )
        return outer_lo, angle

    if abs(np.dot(outer_axis, inner_axis)) > 1.0 - _AXES_PARALLEL_EPS:
        # axes paralleles : 1 seul DOF effectif pour la paire (deux
        # rotations autour de la meme droite s'additionnent). Reduit a
        # l'axe interne seul, avec la plage COMBINEE (seule la somme des
        # deux angles compte pour la direction resultante) ; l'angle
        # externe est fixe a sa borne basse, la difference est absorbee
        # par l'interne.
        combined_spec = RotationAxisSpec(
            axis=inner_axis, angle_range_deg=(inner_lo + outer_lo, inner_hi + outer_hi)
        )
        combined = _resolve_single_axis_analytic(
            combined_spec, reference_effective, target, fallback_angles_deg[0] + fallback_angles_deg[1]
        )
        inner_angle = _clamp_into_range(combined - outer_lo, inner_lo, inner_hi)
        return outer_lo, inner_angle

    target_angle_to_outer = _angle_between(target, outer_axis)

    def _achieved_angle_to_outer(inner_angle_deg: float) -> float:
        rotated = _rodrigues_rotate_single(
            reference_effective, inner_axis, np.radians(inner_angle_deg)
        )
        return _angle_between(rotated, outer_axis)

    # Bracketing VECTORISE (un seul appel _rodrigues_rotate batch sur tous
    # les echantillons, pas une boucle Python de _BRACKETING_SAMPLES appels
    # scalaires) - meme discipline de performance que le reste du fichier
    # (limite #37, CLAUDE.md) : cette fonction est appelee une fois par
    # point de trajectoire, le cout du bracketing doit rester marginal face
    # a la bissection elle-meme.
    samples = np.linspace(inner_lo, inner_hi, max(_BRACKETING_SAMPLES, 2))
    rotated_batch = _rodrigues_rotate(
        np.tile(reference_effective, (len(samples), 1)), inner_axis, np.radians(samples)
    )
    cos_angles = np.clip(rotated_batch @ outer_axis, -1.0, 1.0)
    achieved_batch = np.degrees(np.arccos(cos_angles))
    residuals = achieved_batch - target_angle_to_outer

    bracket = None
    sign_changes = np.where(np.diff(np.sign(residuals)) != 0)[0]
    if len(sign_changes) > 0:
        i = int(sign_changes[0])
        bracket = (samples[i], samples[i + 1], residuals[i], residuals[i + 1])
    else:
        exact = np.where(residuals == 0.0)[0]
        if len(exact) > 0:
            i = int(exact[0])
            bracket = (samples[i], samples[i], residuals[i], residuals[i])

    if bracket is None:
        return None

    lo_t, hi_t, lo_r, _hi_r = bracket
    inner_angle = lo_t
    for _ in range(_BISECTION_MAX_ITERATIONS):
        mid_t = 0.5 * (lo_t + hi_t)
        mid_r = _achieved_angle_to_outer(mid_t) - target_angle_to_outer
        inner_angle = mid_t
        if abs(mid_r) < _BISECTION_TOLERANCE_DEG:
            break
        if (lo_r > 0) == (mid_r > 0):
            lo_t, lo_r = mid_t, mid_r
        else:
            hi_t = mid_t

    intermediate = _rodrigues_rotate_single(
        reference_effective, inner_axis, np.radians(inner_angle)
    )
    delta = _signed_angle_in_plane(intermediate, target, outer_axis)
    if delta is None:
        outer_angle = _clamp_into_range(fallback_angles_deg[0], outer_lo, outer_hi)
    else:
        outer_angle = _clamp_into_range(delta, outer_lo, outer_hi)

    return outer_angle, inner_angle


def resolve_axis_angles_analytic(
    axes: "list",
    reference_direction: np.ndarray,
    target_direction: np.ndarray,
    previous_angles: Optional[list] = None,
    fallback: Optional[Callable[[], tuple[list, float]]] = None,
) -> tuple[list, float]:
    """
    Resolution cinematique inverse ANALYTIQUE, generique pour N axes (voir
    docstring de module). Remplace la recherche par grille de step8 pour le
    chemin principal ; ne leve jamais d'exception et ne renvoie jamais de
    NaN, quelle que soit la configuration d'axes.

    Args:
        axes: list[RotationAxisSpec], meme convention d'ordre que partout
            ailleurs dans ce fichier (axes[0] = le plus externe/applique en
            dernier, axes[-1] = le plus interne/applique en premier).
        reference_direction: direction de depart (angles tous nuls).
        target_direction: direction desiree (n'a pas besoin d'etre unitaire,
            normalisee ici).
        previous_angles: angles resolus au point de trajectoire precedent,
            utilises comme (a) valeur de repliement des axes internes
            (axes[2:], regle de continuite explicite - voir plan) et (b)
            repli en cas de singularite locale (vecteur (anti)parallele a
            un axe). None au premier point : repli sur 0/milieu de plage.
        fallback: fonction `() -> (list[float], float)` appelee UNIQUEMENT
            si la paire externe est geometriquement incapable d'atteindre
            la cible (verifie que ca arrive reellement, pas juste en
            theorie - voir _resolve_axis_pair_analytic). Si None, une
            erreur de config serait un bug appelant : ne devrait jamais
            arriver en pratique cote step8, qui fournit toujours le repli
            grille+raffinement existant.

    Returns:
        (angles_deg, residual_error_rad) - meme contrat que
        refine_axis_angles_local_search : la residuelle est TOUJOURS
        recalculee par cinematique directe, jamais supposee nulle.
    """
    reference_direction = np.asarray(reference_direction, dtype=float)
    reference_direction = reference_direction / np.linalg.norm(reference_direction)
    target = np.asarray(target_direction, dtype=float)
    target_norm = np.linalg.norm(target)
    if target_norm > 1e-12:
        target = target / target_norm

    n = len(axes)

    if n == 0:
        return [], _angle_between(reference_direction, target) * np.pi / 180.0

    def _previous_for(index: int) -> float:
        if previous_angles is not None and index < len(previous_angles):
            return previous_angles[index]
        lo, hi = axes[index].angle_range_deg
        if hi - lo >= 360.0:
            return 0.0
        return lo if lo <= 0.0 <= hi else 0.5 * (lo + hi)

    if n == 1:
        angle = _resolve_single_axis_analytic(
            axes[0], reference_direction, target, _previous_for(0)
        )
        achieved = _forward_direction(axes, [angle], reference_direction)
        residual_rad = np.radians(_angle_between(achieved, target))
        return [angle], float(residual_rad)

    # Repliement des axes internes (axes[2:]) a une valeur fixe - regle de
    # continuite explicite (herite du point precedent), PAS une necessite
    # mathematique : voir plan/CLAUDE.md pour la discussion de ce choix de
    # conception et de son risque (la paire externe peut alors se
    # retrouver incapable d'atteindre la cible - d'ou le repli).
    inner_angles = [_previous_for(i) for i in range(2, n)]
    reference_effective = reference_direction
    for spec, angle_deg in zip(reversed(axes[2:]), reversed(inner_angles)):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        reference_effective = _rodrigues_rotate_single(
            reference_effective, axis, np.radians(angle_deg)
        )

    fallback_angles = (_previous_for(0), _previous_for(1))
    solved_pair = _resolve_axis_pair_analytic(
        axes[0], axes[1], reference_effective, target, fallback_angles
    )

    if solved_pair is None:
        if fallback is None:
            # Ne devrait jamais arriver (step8 fournit toujours un repli) -
            # geler la paire a son etat precedent plutot que planter, meme
            # discipline "jamais d'exception pour une orientation
            # inatteignable" que le reste du fichier.
            outer_angle, inner_angle = fallback_angles
        else:
            return fallback()
    else:
        outer_angle, inner_angle = solved_pair

    angles = [outer_angle, inner_angle] + inner_angles
    achieved = _forward_direction(axes, angles, reference_direction)
    residual_rad = np.radians(_angle_between(achieved, target))
    return angles, float(residual_rad)


def point_has_support_along_direction(
    position: np.ndarray,
    direction_up: np.ndarray,
    candidate_starts: np.ndarray,
    candidate_ends: np.ndarray,
    horizontal_radius: float,
    vertical_height: float,
) -> bool:
    """
    `position` a-t-il un support (un point d'un des segments
    `candidate_starts[i]`-`candidate_ends[i]`) STRICTEMENT en dessous de lui
    le long de `direction_up` - c'est-a-dire un point du segment a moins de
    `horizontal_radius` transversalement et entre 0 et `vertical_height`
    "sous" lui le long de `-direction_up` ?

    Generalise le test vertical historique (implicitement `direction_up =
    [0,0,1]`, teste via position[2]/position[0:2] directement) a une
    direction QUELCONQUE : construit une base orthonormee locale (u, v,
    direction_up) et reprojette position + segments dedans avant d'appliquer
    exactement la meme geometrie (projection scalaire clampee sur le
    segment, plus proche point, test horizontal/vertical separe). Utilisee a
    la fois par pipeline/step9b_support_check.py (support d'un point de
    trajectoire contre des capsules imprimees, testee sur PLUSIEURS
    directions atteignables par la machine) et par pipeline/step5_layers.py
    (support d'un point de contour contre des clusters deja places, testee
    UNIQUEMENT selon up_vector - step5 n'a pas encore resolu axis_angles).

    `candidate_starts`/`candidate_ends` doivent deja etre le resultat d'un
    filtrage grossier (AABB) fait par l'appelant - cette fonction ne fait
    aucune recherche spatiale, seulement la geometrie precise sur les
    candidats deja retenus. Vectorise sur tous les candidats en une passe
    numpy (pas de boucle Python par candidat), meme discipline que
    step9._points_to_segment_distances.

    Retourne False si `candidate_starts` est vide (pas de candidat a
    tester).
    """
    if len(candidate_starts) == 0:
        return False

    d = np.asarray(direction_up, dtype=float)
    norm = np.linalg.norm(d)
    d = d / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0])

    # base orthonormee (u, v, d) - meme construction que _perpendicular_basis
    # (reference fixe [1,0,0], repli [0,1,0] si quasi-colineaire), pour
    # rester coherente avec le reste du fichier plutot qu'une 3e convention.
    reference = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(reference, d)) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    u = reference - np.dot(reference, d) * d
    u = u / np.linalg.norm(u)
    v = np.cross(d, u)
    basis = np.stack([u, v, d], axis=1)  # colonnes = base locale

    position_local = basis.T @ np.asarray(position, dtype=float)
    starts_local = np.asarray(candidate_starts, dtype=float) @ basis
    ends_local = np.asarray(candidate_ends, dtype=float) @ basis

    segments = ends_local - starts_local
    length_sq = np.einsum("ij,ij->i", segments, segments)
    degenerate = length_sq < 1e-18
    safe_length = np.where(degenerate, 1.0, length_sq)
    t = np.clip(
        np.einsum("ij,ij->i", position_local[None, :] - starts_local, segments) / safe_length,
        0.0,
        1.0,
    )
    t = np.where(degenerate, 0.0, t)
    closest = starts_local + t[:, None] * segments

    delta_up = position_local[2] - closest[:, 2]
    horizontal = np.hypot(position_local[0] - closest[:, 0], position_local[1] - closest[:, 1])

    return bool(
        np.any(
            (delta_up > 1e-9) & (delta_up <= vertical_height) & (horizontal <= horizontal_radius)
        )
    )


def point_in_polygon_2d(point: np.ndarray, polygon: np.ndarray) -> bool:
    """
    Ray-casting standard : `point` (2,) est-il a l'interieur de `polygon`
    (N,2) ? Technique robuste et independante de toute convention de
    winding, contrairement au signe de l'aire projetee (verifie empiriquement
    faux sur un anneau plat synthetique - CLAUDE.md limite #11).
    """
    from core import _cpp  # import local : core._cpp importe deja ce module indirectement
    if _cpp.enabled() and hasattr(_cpp.module(), "point_in_polygon_2d"):
        return bool(_cpp.module().point_in_polygon_2d(float(point[0]), float(point[1]),
                                                      np.ascontiguousarray(polygon, dtype=float)))
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
