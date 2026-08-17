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

from core.geometry_utils import composed_rotation_directions, rotation_axis_circle_directions
from core.mesh_types import RotationAxisSpec


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
