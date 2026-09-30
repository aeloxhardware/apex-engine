"""
tests/test_geometry_utils.py

Tests pour core/geometry_utils.py, en particulier composed_rotation_directions
(generalisation a une cinematique de rotation arbitraire a N axes composes,
voir CLAUDE.md limite de conception #18).

Execution : pytest tests/test_geometry_utils.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from core.geometry_utils import (
    _rodrigues_rotate,
    _rodrigues_rotate_single,
    composed_rotation_directions,
    resolve_axis_angles_analytic,
    rotation_axis_circle_directions,
)
from core.mesh_types import RotationAxisSpec


def _forward(axes: list, angles_deg: list, reference: np.ndarray) -> np.ndarray:
    """Cinematique directe de reference pour les tests - meme ordre de
    composition (le plus interne d'abord) que le reste du fichier."""
    direction = reference.copy()
    for spec, angle_deg in zip(reversed(axes), reversed(angles_deg)):
        axis = spec.axis / np.linalg.norm(spec.axis)
        direction = _rodrigues_rotate_single(direction, axis, np.radians(angle_deg))
    return direction


# ---------------------------------------------------------------------------
# Degenerescences : le modele general doit reproduire les cas simples existants
# ---------------------------------------------------------------------------

def test_empty_axes_returns_reference_direction_only():
    """Aucune rotation -> la direction de reference seule, inchangee."""
    reference = np.array([0.3, 0.4, 0.5])
    directions = composed_rotation_directions([], reference, samples_per_axis=8)

    assert directions.shape == (1, 3)
    normalized = reference / np.linalg.norm(reference)
    assert np.allclose(directions[0], normalized)


def test_single_axis_full_turn_matches_rotation_axis_circle_directions():
    """
    Un seul axe a plage complete (0,360) doit produire EXACTEMENT le meme
    ensemble de directions que rotation_axis_circle_directions - preuve que
    le modele general n'est pas une approximation du cas a 1 DOF, mais une
    vraie generalisation qui le contient.
    """
    axis = np.array([0.0, 0.0, 1.0])
    reference = np.array([1.0, 0.0, 0.0])
    count = 8

    composed = composed_rotation_directions(
        [RotationAxisSpec(axis=axis, angle_range_deg=(0, 360))], reference, samples_per_axis=count
    )
    circle = rotation_axis_circle_directions(axis, count)

    assert composed.shape == circle.shape
    assert np.allclose(sorted(composed.tolist()), sorted(circle.tolist()), atol=1e-9)


def test_single_axis_locked_produces_one_sample():
    """Un DOF fige (min == max) ne doit produire qu'un seul echantillon, pas samples_per_axis."""
    axis = np.array([0.0, 0.0, 1.0])
    reference = np.array([1.0, 0.0, 0.0])

    directions = composed_rotation_directions(
        [RotationAxisSpec(axis=axis, angle_range_deg=(45.0, 45.0))], reference, samples_per_axis=16
    )

    assert directions.shape == (1, 3)
    # rotation de 45deg de (1,0,0) autour de Z
    expected = np.array([np.cos(np.radians(45.0)), np.sin(np.radians(45.0)), 0.0])
    assert np.allclose(directions[0], expected, atol=1e-9)


# ---------------------------------------------------------------------------
# Cas machine actuelle : couverture spherique complete MESUREE, pas supposee
# ---------------------------------------------------------------------------

def test_two_axis_current_machine_kinematics_covers_full_sphere():
    """
    Machine actuelle : axe plateau (tour complet) + axe tete perpendiculaire
    (demi-plage 0-180deg). CLAUDE.md limite #18 affirme que cette
    configuration couvre la sphere complete - verifie ici par le calcul
    (dispersion des z, x, y sur toute leur plage possible), pas suppose vrai
    dans le code de production.
    """
    platform_axis = np.array([0.0, 0.0, 1.0])
    head_axis = np.array([1.0, 0.0, 0.0])
    reference = np.array([0.0, 0.0, -1.0])  # buse vers le bas a B=0

    axes = [
        RotationAxisSpec(axis=platform_axis, angle_range_deg=(0, 360)),
        RotationAxisSpec(axis=head_axis, angle_range_deg=(0, 180)),
    ]
    directions = composed_rotation_directions(axes, reference, samples_per_axis=24)

    norms = np.linalg.norm(directions, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-9), "Toutes les directions doivent rester unitaires"

    # couverture spherique complete : les 3 composantes doivent chacune
    # approcher leurs bornes theoriques [-1, 1]
    for axis_idx in range(3):
        assert directions[:, axis_idx].min() < -0.95, f"Composante {axis_idx} ne couvre pas le bas de sa plage"
        assert directions[:, axis_idx].max() > 0.95, f"Composante {axis_idx} ne couvre pas le haut de sa plage"


def test_two_axis_restricted_kinematics_does_not_cover_full_sphere():
    """
    Contre-preuve explicitement requise : une cinematique a 2 axes tous deux
    a plage restreinte (0-90deg chacun, PAS la configuration de la machine
    actuelle) ne
    doit PAS couvrir la sphere complete - prouve que le systeme ne fait
    aucune hypothese cachee du type "2 axes => sphere complete", ce qui
    serait retomber sur une reference figee deguisee.
    """
    axis_a = np.array([0.0, 0.0, 1.0])
    axis_b = np.array([0.0, 1.0, 0.0])
    reference = np.array([0.0, 0.0, -1.0])

    axes = [
        RotationAxisSpec(axis=axis_a, angle_range_deg=(0, 90)),
        RotationAxisSpec(axis=axis_b, angle_range_deg=(0, 90)),
    ]
    directions = composed_rotation_directions(axes, reference, samples_per_axis=16)

    # une cinematique restreinte sur les deux DOF ne peut pas atteindre tout
    # l'hemisphere oppose a la reference - la composante correspondante ne
    # doit pas approcher +1 comme dans le cas de la machine actuelle
    assert directions[:, 2].max() < 0.95, (
        "Une cinematique a 2 axes tous deux restreints ne devrait pas "
        "couvrir la sphere complete, contrairement au cas de la machine actuelle"
    )


# ---------------------------------------------------------------------------
# Coherence generale
# ---------------------------------------------------------------------------

def test_three_axes_all_unit_length():
    """Verifie que la composition reste stable (directions unitaires) au-dela de 2 axes."""
    axes = [
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0, 360)),
        RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0, 180)),
        RotationAxisSpec(axis=np.array([0.0, 1.0, 0.0]), angle_range_deg=(-30, 30)),
    ]
    directions = composed_rotation_directions(axes, np.array([0.0, 0.0, -1.0]), samples_per_axis=6)

    norms = np.linalg.norm(directions, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-9)
    assert directions.shape[0] == 6 ** 3


# ---------------------------------------------------------------------------
# _rodrigues_rotate_single : chemin scalaire dedie (CLAUDE.md, perf step8)
# ---------------------------------------------------------------------------

def test_rodrigues_rotate_single_matches_batched_path():
    """
    _rodrigues_rotate_single (un seul vecteur, sans np.cross generique) doit
    produire un resultat identique a _rodrigues_rotate (chemin batch
    existant, inchange) pour la meme rotation - seule la mecanique interne
    change, jamais le resultat mathematique.
    """
    cases = [
        (np.array([0.1, 0.2, 0.97]), np.array([0.0, 0.0, 1.0]), 0.3),
        (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.pi / 2),
        (np.array([0.5, -0.5, 0.7071]), np.array([1.0, 1.0, 0.0]) / np.sqrt(2), 1.234),
        (np.array([0.0, 0.0, -1.0]), np.array([1.0, 0.0, 0.0]), 0.0),  # angle nul
    ]
    for vector, axis, angle_rad in cases:
        expected = _rodrigues_rotate(vector.reshape(1, 3), axis, np.array([angle_rad]))[0]
        actual = _rodrigues_rotate_single(vector, axis, angle_rad)
        assert np.allclose(actual, expected, atol=1e-12), (vector, axis, angle_rad)


# ---------------------------------------------------------------------------
# resolve_axis_angles_analytic : IK generique N-axes par bissection/atan2
# (remplace la recherche par grille de step8 - voir CLAUDE.md "IK analytique
# implementee" pour le diagnostic complet et le risque de N>=3).
# ---------------------------------------------------------------------------

_PLATFORM_AXIS = np.array([0.0, 0.0, 1.0])
_HEAD_AXIS = np.array([1.0, 0.0, 0.0])
_MACHINE_AXES = [
    RotationAxisSpec(axis=_PLATFORM_AXIS, angle_range_deg=(0.0, 360.0)),
    RotationAxisSpec(axis=_HEAD_AXIS, angle_range_deg=(0.0, 180.0)),
]
_REFERENCE = np.array([0.0, 0.0, -1.0])  # meme reference_direction que step8.run()


def test_resolve_two_axis_real_machine_round_trip_is_near_exact():
    """
    Round-trip sur la machine reelle (2 axes orthogonaux) : le repliement
    interne est vide pour N=2, la methode retombe IDENTIQUEMENT sur la
    fermeture 2-axes deja validee - erreur attendue de l'ordre de la
    tolerance de bissection (1e-5 deg), tres largement sous le seuil de
    rejet de step8 (5 deg) et sous l'erreur mediane mesuree de la grille
    historique (9.26 deg).
    """
    rng = np.random.default_rng(1)
    max_error = 0.0
    for _ in range(200):
        platform_deg = rng.uniform(0.0, 360.0)
        head_deg = rng.uniform(0.0, 180.0)
        target = _forward(_MACHINE_AXES, [platform_deg, head_deg], _REFERENCE)

        angles, residual_rad = resolve_axis_angles_analytic(_MACHINE_AXES, _REFERENCE, target)
        achieved = _forward(_MACHINE_AXES, angles, _REFERENCE)
        error_deg = np.degrees(np.arccos(np.clip(np.dot(achieved, target), -1.0, 1.0)))

        max_error = max(max_error, error_deg)
        assert abs(np.degrees(residual_rad) - error_deg) < 1e-4, (
            "la residuelle retournee doit correspondre a l'ecart REEL recalcule"
        )

    assert max_error < 1e-3, f"erreur max {max_error} deg - attendu proche de 1e-5 deg"


def test_resolve_zero_axes_returns_empty_list():
    angles, residual_rad = resolve_axis_angles_analytic([], _REFERENCE, np.array([1.0, 0.0, 0.0]))
    assert angles == []
    assert residual_rad >= 0.0


def test_resolve_single_axis_matches_direct_atan2():
    """N=1 : un seul DOF, doit converger vers l'angle exact quand la cible
    est atteignable (meme angle par rapport a l'axe que la reference)."""
    axis_spec = RotationAxisSpec(axis=_PLATFORM_AXIS, angle_range_deg=(0.0, 360.0))
    target_angle_deg = 137.0
    target = _forward([axis_spec], [target_angle_deg], _REFERENCE)

    angles, residual_rad = resolve_axis_angles_analytic([axis_spec], _REFERENCE, target)

    assert len(angles) == 1
    assert np.degrees(residual_rad) < 1e-3
    achieved = _forward([axis_spec], angles, _REFERENCE)
    assert np.allclose(achieved, target, atol=1e-6)


def test_resolve_locked_single_axis_returns_locked_value_without_search():
    """Axe verrouille (lo==hi) : retourne la valeur verrouillee sans
    recherche, meme si la cible est inatteignable - meme contrat que
    refine_axis_angles_local_search (residuelle non nulle, pas d'exception)."""
    axis_spec = RotationAxisSpec(axis=_PLATFORM_AXIS, angle_range_deg=(42.0, 42.0))
    target = np.array([1.0, 0.0, 0.0])  # tres different de ce que 42deg peut produire

    angles, residual_rad = resolve_axis_angles_analytic([axis_spec], _REFERENCE, target)

    assert angles == [42.0]
    assert residual_rad > 0.0


def test_resolve_three_axis_nonorthogonal_never_crashes_or_produces_nan():
    """
    N=3, axes non orthogonaux : le repliement interne (axes[2:]) est une
    regle de continuite explicite, pas une necessite mathematique (CLAUDE.md)
    - la residuelle n'est donc PAS forcement proche de 0 pour une cible
    choisie au hasard. Ce test verifie l'invariant de robustesse (jamais de
    NaN/exception, residuelle toujours coherente avec un recalcul direct),
    pas une precision particuliere.
    """
    axes = [
        RotationAxisSpec(axis=np.array([0.5, 0.3, 0.8]) / np.linalg.norm([0.5, 0.3, 0.8]), angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=np.array([0.9, -0.2, 0.1]) / np.linalg.norm([0.9, -0.2, 0.1]), angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=np.array([0.1, 0.9, 0.3]) / np.linalg.norm([0.1, 0.9, 0.3]), angle_range_deg=(0.0, 360.0)),
    ]
    rng = np.random.default_rng(9)
    for _ in range(100):
        target = rng.normal(size=3)
        target = target / np.linalg.norm(target)
        previous = [rng.uniform(0.0, 360.0) for _ in axes]

        angles, residual_rad = resolve_axis_angles_analytic(axes, _REFERENCE, target, previous)

        assert len(angles) == 3
        assert all(np.isfinite(a) for a in angles)
        assert np.isfinite(residual_rad) and residual_rad >= 0.0

        achieved = _forward(axes, angles, _REFERENCE)
        recomputed_deg = np.degrees(
            np.arccos(np.clip(np.dot(achieved, target), -1.0, 1.0))
        )
        assert abs(np.degrees(residual_rad) - recomputed_deg) < 1e-3, (
            "la residuelle doit toujours correspondre a l'ecart REEL, jamais suppose nul"
        )


def test_resolve_three_axis_reduces_to_two_axis_case_when_inner_folded_to_target():
    """
    Cas ou le repliement interne EST la valeur qui a produit la cible
    (previous_angles == les vrais angles internes utilises) : la paire
    externe doit alors pouvoir reconstruire la cible avec une erreur
    proche de zero, exactement comme le cas N=2 - preuve que le
    repliement n'est pas juste "ne plante pas" mais fonctionne
    correctement quand il correspond a la geometrie reelle.
    """
    inner_axis = np.array([0.0, 1.0, 0.0])
    axes = [
        RotationAxisSpec(axis=_PLATFORM_AXIS, angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=_HEAD_AXIS, angle_range_deg=(0.0, 180.0)),
        RotationAxisSpec(axis=inner_axis, angle_range_deg=(-30.0, 30.0)),
    ]
    rng = np.random.default_rng(3)
    max_error = 0.0
    for _ in range(50):
        true_angles = [rng.uniform(0.0, 360.0), rng.uniform(0.0, 180.0), rng.uniform(-30.0, 30.0)]
        target = _forward(axes, true_angles, _REFERENCE)

        angles, residual_rad = resolve_axis_angles_analytic(
            axes, _REFERENCE, target, previous_angles=true_angles
        )
        achieved = _forward(axes, angles, _REFERENCE)
        error_deg = np.degrees(np.arccos(np.clip(np.dot(achieved, target), -1.0, 1.0)))
        max_error = max(max_error, error_deg)

    assert max_error < 1e-2, f"erreur max {max_error} deg avec repliement exact"


def test_resolve_parallel_axis_pair_reduces_to_single_dof():
    """Paire externe avec axes PARALLELES : un seul DOF effectif, doit
    quand meme resoudre sans exception ni NaN (plage combinee)."""
    axes = [
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 180.0)),
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 180.0)),
    ]
    target = _forward(axes, [40.0, 30.0], _REFERENCE)  # somme = 70deg autour de Z

    angles, residual_rad = resolve_axis_angles_analytic(axes, _REFERENCE, target)

    assert len(angles) == 2
    assert all(np.isfinite(a) for a in angles)
    assert np.isfinite(residual_rad)


def test_resolve_falls_back_when_outer_pair_cannot_reach_target():
    """
    Cas confirme REEL pendant la conception (pas hypothetique) : sur des
    axes non orthogonaux, la paire externe seule peut etre geometriquement
    incapable d'atteindre la cible pour un repliement interne donne -
    l'angle atteignable par rapport a axes[0] est borne a un ANNEAU, pas
    [0,180] en general. Le fallback doit alors etre appele et son resultat
    utilise tel quel, jamais une exception.
    """
    axes = [
        RotationAxisSpec(axis=np.array([0.5, 0.3, 0.8]) / np.linalg.norm([0.5, 0.3, 0.8]), angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=np.array([0.9, -0.2, 0.1]) / np.linalg.norm([0.9, -0.2, 0.1]), angle_range_deg=(0.0, 360.0)),
    ]
    # cible choisie hors de l'anneau atteignable pour ref=_REFERENCE (verifie
    # empiriquement pendant la conception - voir plan) : angle tres proche
    # de l'axe externe, difficile a atteindre pour cette paire specifique.
    target = axes[0].axis

    fallback_called = []

    def fallback():
        fallback_called.append(True)
        return [11.0, 22.0], 0.5

    angles, residual_rad = resolve_axis_angles_analytic(
        axes, _REFERENCE, target, fallback=fallback
    )

    # Cette cible/ces axes DECLENCHENT reellement le repli (verifie
    # explicitement, pas suppose) : sinon ce test ne prouverait rien sur
    # le mecanisme de secours.
    assert fallback_called, "cette configuration devrait declencher le repli grille"
    assert angles == [11.0, 22.0]
    assert residual_rad == 0.5
