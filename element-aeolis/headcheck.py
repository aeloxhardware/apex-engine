"""Contrôle de collision de la vraie tête (buse verticale) contre la matière déjà déposée.

La tête (GLB, pointe à l'origine, corps vers +z) est réduite à son ENVELOPPE INFÉRIEURE sur une
grille polaire (rayon × azimut) : pour chaque cellule, le point le plus bas. C'est exact pour
le test « un point de la tête passe-t-il sous la surface imprimée ? » puisque seul le plus bas
d'une colonne verticale peut toucher une surface en hauteur.

Pour un point de dépôt p, la matière déjà posée a pour hauteur, en (x, y), la surface de la
couche déposée juste avant ; la couche en cours n'est pas comptée (la buse est confondue avec
le cordon qu'elle dépose, comme à l'étape 9 du pipeline principal). Hors de l'empreinte, il n'y
a que le plateau.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import trimesh

DEFAULT_HEAD = Path(__file__).resolve().parent.parent / "res" / "3d" / "toolheads" / "aetheris-4" / "aetheris-4.glb"


@dataclass
class HeadEnvelope:
    dx: np.ndarray                # (M,) décalage horizontal par rapport à la pointe
    dy: np.ndarray
    h: np.ndarray                 # (M,) hauteur au-dessus de la pointe
    clearance_slope_deg: float    # pente d'un plan montant tolérée dans le pire azimut


@dataclass
class HeadReport:
    checked: bool = False
    points: int = 0
    colliding_points: int = 0
    max_penetration: float = 0.0
    clearance_slope_deg: float = 0.0


def load_envelope(path: str | Path | None = None, contact_clearance: float = 3.0,
                  r_step: float = 0.5, n_az: int = 72, samples: int = 60000) -> HeadEnvelope:
    from file_io.mesh_loader import load_toolhead_mesh
    mesh = load_toolhead_mesh(str(path or DEFAULT_HEAD))
    pts = np.vstack([mesh.vertices,
                     trimesh.sample.sample_surface(mesh, samples, seed=0)[0]])
    pts = pts[np.linalg.norm(pts, axis=1) >= contact_clearance]   # zone de contact buse/cordon
    r = np.hypot(pts[:, 0], pts[:, 1])
    az = np.arctan2(pts[:, 1], pts[:, 0])
    ir = np.minimum((r / r_step).astype(int), int(r.max() / r_step))
    ia = ((az + np.pi) / (2 * np.pi) * n_az).astype(int) % n_az
    key = ir * n_az + ia
    order = np.lexsort((pts[:, 2], key))                           # par cellule, z croissant
    first = order[np.r_[True, np.diff(key[order]) != 0]]
    low = pts[first]
    rr = np.hypot(low[:, 0], low[:, 1])
    slope = np.degrees(np.arctan2(low[:, 2], np.maximum(rr, 1e-9)))
    return HeadEnvelope(low[:, 0], low[:, 1], low[:, 2], float(slope.min()))


def check_points(env: HeadEnvelope, field, z_fn, t_fn, xy: np.ndarray, z_tip: np.ndarray,
                 margin: float = 0.0, printed_fn=None) -> tuple[int, float]:
    """Nombre de points de dépôt en collision et pénétration maximale (mm).

    `z_fn(x, y)` : surface de la couche en cours ; `t_fn(x, y)` : son épaisseur locale. La matière
    déjà posée vaut z_fn − t_fn. `printed_fn(x, y)` (facultatif) rend directement cette différence : il
    évite d'interpoler deux fois la grille aux mêmes points.
    """
    if len(xy) == 0:
        return 0, 0.0
    # Un point de la tête plus haut que toute la pièce ne peut toucher aucune surface imprimée.
    reach = float(np.nanmax(field.top) - np.nanmin(field.bottom)) + margin
    keep = env.h <= reach
    env = HeadEnvelope(env.dx[keep], env.dy[keep], env.h[keep], env.clearance_slope_deg)
    colliding, worst = 0, 0.0
    x0, x1, y0, y1 = field.xs[0], field.xs[-1], field.ys[0], field.ys[-1]
    for lo in range(0, len(xy), 400):
        p = xy[lo:lo + 400]
        zt = z_tip[lo:lo + 400]
        hx = p[:, 0:1] + env.dx[None, :]
        hy = p[:, 1:2] + env.dy[None, :]
        inside = (hx >= x0) & (hx <= x1) & (hy >= y0) & (hy <= y1)
        fx, fy = np.clip(hx, x0, x1), np.clip(hy, y0, y1)
        if printed_fn is not None:
            printed = printed_fn(fx.ravel(), fy.ravel()).reshape(hx.shape)
        else:
            printed = z_fn(fx.ravel(), fy.ravel()).reshape(hx.shape) - t_fn(fx.ravel(), fy.ravel()).reshape(hx.shape)
        gap = (zt[:, None] + env.h[None, :]) - printed - margin       # < 0 : la tête est dans la matière
        gap = np.where(inside & ~np.isnan(gap), gap, np.inf)
        pen = -gap.min(axis=1)
        hit = pen > 1e-6
        colliding += int(hit.sum())
        if hit.any():
            worst = max(worst, float(pen[hit].max()))
    return colliding, worst
