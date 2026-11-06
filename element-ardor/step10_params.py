"""Etape 10 : paramètres et type du hook de re-segmentation.

Découpé de `step10_adjustment.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step10_adjustment`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional
from core.mesh_types import ResegmentationRequest, ToolpathMachine


# Un hook recoit une demande de re-segmentation et retourne la nouvelle
# trajectoire complete (step3 -> step8 re-joues par l'appelant), ou None
# s'il ne peut pas satisfaire la demande.
ResegmentationHook = Callable[[ResegmentationRequest], Optional[ToolpathMachine]]


@dataclass
class Step10Params:
    max_adjustment_iterations: int = 3
    # Grille de recherche de la reorientation locale - meme parametrage que
    # la recherche NN de step8 (grille produit sur les plages d'axes).
    reorientation_samples_per_axis: int = 12
    # Seuil de surplomb : ecart angulaire maximal tolere entre la direction
    # de depot d'origine et celle du candidat de reorientation. Un candidat
    # qui elimine la collision mais depasse ce seuil reintroduirait un
    # surplomb - il est rejete, c'est exactement la condition "tout en
    # respectant encore le seuil de surplomb" de la specification.
    max_direction_deviation_deg: float = 45.0
    # Demi-largeur de la plage d'angles signalee comme problematique dans
    # une ResegmentationRequest (autour de l'angle courant de chaque axe).
    excluded_range_half_width_deg: float = 15.0
    # Nombre de voisins de part et d'autre du point modifie sur lesquels
    # valider les limites de vitesse/acceleration angulaire.
    rate_limit_neighbor_margin: int = 1
    # Diagnostics d'ORIENTATION de l'etape 8 (ecart tete / direction visee,
    # sans collision) : seuls ceux dont l'ecart residuel depasse ce seuil
    # (degres) entrent dans la boucle (re-segmentation, puis support de
    # secours). Les autres sont journalises `orientation_deviation` et ne
    # comptent PAS comme supports. None (DEFAUT, decision de l'utilisateur
    # du 2026-10-07) = jamais : un support ne corrige pas une orientation.
    # 0.0 = comportement historique (tout diagnostic devient support :
    # 115 312 sur death contre 209 vraies collisions, limite #74).
    # Une collision reelle au meme point suit toujours la voie normale.
    orientation_support_min_residual_deg: Optional[float] = None
