"""Etape 7 : paramètres et erreur d'ordonnancement.

Découpé de `step7_scheduling.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step7_scheduling`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
from core.mesh_types import RotationAxisSpec


class SchedulingError(Exception):
    """Levee si le planning ne peut pas etre construit (donnees incoherentes)."""


@dataclass
class HoleCostModel:
    """Cout d'un trou temporaire (contrat UI : SchedulingParams.hole_cost_model)."""
    reheat_penalty: float = 1.0
    adhesion_risk_weight: float = 1.0


@dataclass
class PurgeCostModel:
    """Cout d'une purge classique (contrat UI : SchedulingParams.purge_cost_model)."""
    volume_per_purge: float = 1.0
    time_per_purge: float = 1.0


@dataclass
class Step7Params:
    # Distance max entre un point de trajectoire et la surface du mesh
    # au-dela de laquelle ce point est considere INTERNE (invisible), donc
    # sans couleur imposee. Typiquement de l'ordre d'une largeur de buse.
    surface_color_distance: float = 0.5
    hole_cost_model: HoleCostModel = field(default_factory=HoleCostModel)
    purge_cost_model: PurgeCostModel = field(default_factory=PurgeCostModel)
    max_open_holes: int = 4
    # Directions d'echappement testees pour la deadline d'accessibilite.
    escape_direction_count: int = 16
    # Cinematique de rotation de la machine : si None/vide, les directions
    # d'echappement sont prises sur la sphere complete (orientation libre).
    # Si fourni, decrit une cinematique arbitraire (N axes composes, chacun
    # avec sa propre plage) via composed_rotation_directions - PAS figee sur
    # une seule machine (voir CLAUDE.md limite de conception #18). Un seul
    # axe a plage complete degenere vers l'ancien comportement a 1 DOF.
    rotation_axes: Optional[list[RotationAxisSpec]] = None
    # Nombre d'echantillons le long d'un rayon d'echappement.
    deadline_ray_samples: int = 24
