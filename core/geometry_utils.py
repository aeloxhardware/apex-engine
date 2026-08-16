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
