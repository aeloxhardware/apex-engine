"""Cinematique d'Aetheris-4 Gen 1 : portique XYZ, deux rotations A et B PORTEES PAR LA TETE.

Le firmware ne fait AUCUNE cinematique (decision de l'utilisateur, 2026-10-31) : le slicer commande
la position des axes physiques, pivot compense. Ce module est la geometrie directe et inverse de cette
compensation, pour une chaine de N rotations de tete (N = 2 ici : A autour de X, B autour de Y, convention
ISO / RS-274) et des decalages fixes entre leurs centres.

Modele (meme convention que le reste du moteur, `core/geometry_utils` et `step9`) :
  - les angles s'appliquent dans l'ordre de `MachineProfile.rotation_axes` : axes[0] est le PARENT
    (le plus proche du chariot), axes[-1] le plus proche de la buse ; R = R_0(a_0) · R_1(a_1) · ... ;
  - a angles nuls la tete pend verticalement, direction de depot (axe du corps, pointe -> pivot) = +Z ;
  - pointe = c + o_0 + R_1 (o_1 + R_2 (... + R_n o_n)) ou `c` est le point commande du chariot
    (X, Y, Z firmware) et o_0 .. o_n les decalages fixes. Cas par defaut (UN centre commun a toutes les
    rotations) : o_0 = o_1 = 0, o_n = (0, 0, -L), L = `pivot_to_tip`.
  - convention de sortie « z_reference = tip » : a angles nuls le point commande EST la pointe.
    Le point commande est donc  fw = pointe - (chaine(angles) - chaine(0)).

VALEURS A CONFIRMER SUR LA MACHINE (l'utilisateur ne les connait pas encore, 2026-10-31) : l'ordre de
montage de A et B, L, d'eventuels decalages distincts entre les deux centres, les plages, le signe de
chaque axe. Rien ici n'est mesure ; tout est parametre.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

AXIS_X = (1.0, 0.0, 0.0)
AXIS_Y = (0.0, 1.0, 0.0)


@dataclass
class HeadChain:
    """Chaine de rotations de tete. `axes` et `letters` suivent l'ordre parent -> buse."""
    axes: list = field(default_factory=lambda: [AXIS_X, AXIS_Y])
    letters: list = field(default_factory=lambda: ["A", "B"])
    pivot_to_tip: float = 46.0                     # L, mm (valeur de la machine polaire, PAS mesuree pour Gen 1)
    # Decalages o_0 .. o_{n-1} (mm, avant le dernier). Le dernier, o_n = (0, 0, -L), vient de pivot_to_tip.
    offsets: Optional[Sequence] = None

    def n(self) -> int:
        return len(self.axes)

    def offset_list(self) -> np.ndarray:
        n = self.n()
        o = np.zeros((n + 1, 3))
        if self.offsets is not None:
            given = np.asarray(self.offsets, dtype=float).reshape(-1, 3)
            if len(given) != n:
                raise ValueError(f"`offsets` attend {n} decalages (un par axe, avant la buse), recu {len(given)}")
            o[:n] = given
        o[n] = (0.0, 0.0, -float(self.pivot_to_tip))
        return o


def _rot(axis, theta):
    """(N,3,3) : rotations main droite autour de `axis` (unitaire) d'angles `theta` (rad, (N,))."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    c = np.cos(theta)[:, None, None]
    s = np.sin(theta)[:, None, None]
    cross = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
    return c * np.eye(3) + s * cross + (1.0 - c) * np.outer(axis, axis)


def rotations(chain: HeadChain, angles_deg) -> list:
    """Liste des n rotations individuelles R_k, chacune (N,3,3)."""
    ang = np.atleast_2d(np.asarray(angles_deg, dtype=float))
    if ang.shape[1] != chain.n():
        raise ValueError(f"{chain.n()} angles attendus par point, recu {ang.shape[1]}")
    return [_rot(a, np.radians(ang[:, k])) for k, a in enumerate(chain.axes)]


def body_axis(chain: HeadChain, angles_deg) -> np.ndarray:
    """(N,3) direction du corps de la tete (pointe -> pivot) = direction de depot, R · (0,0,1)."""
    rs = rotations(chain, angles_deg)
    total = rs[0]
    for r in rs[1:]:
        total = total @ r
    return total[:, :, 2]


def _chain_offset(chain: HeadChain, angles_deg) -> np.ndarray:
    """(N,3) : o_0 + R_1(o_1 + R_2(... + R_n o_n))."""
    rs = rotations(chain, angles_deg)
    o = chain.offset_list()
    n = chain.n()
    p = np.broadcast_to(o[n], (len(rs[0]), 3)).copy()
    for k in range(n, 0, -1):
        p = o[k - 1] + np.einsum("nij,nj->ni", rs[k - 1], p)
    return p


def tip_to_firmware(chain: HeadChain, tip, angles_deg) -> np.ndarray:
    """Pointe (N,3) + angles d'Apex (N,n) -> point commande du chariot (X, Y, Z firmware)."""
    tip = np.atleast_2d(np.asarray(tip, dtype=float))
    return tip - (_chain_offset(chain, angles_deg) - chain.offset_list().sum(axis=0))


def firmware_to_tip(chain: HeadChain, fw, angles_deg) -> np.ndarray:
    """Inverse de `tip_to_firmware` (apercu, controle aller-retour, ecart de pointe a mi-mouvement)."""
    fw = np.atleast_2d(np.asarray(fw, dtype=float))
    return fw + (_chain_offset(chain, angles_deg) - chain.offset_list().sum(axis=0))
