"""Cinématique inverse analytique N axes (rotations de Rodrigues, résolution axe par axe et par paire).

Découpé de `geometry_utils.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `geometry_utils`.
"""

from __future__ import annotations

import math
from typing import Callable, Optional
import numpy as np
from core.mesh_types import RotationAxisSpec


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
