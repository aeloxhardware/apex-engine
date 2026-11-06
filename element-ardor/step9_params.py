"""Etape 9 : paramètres de la simulation de collision.

Découpé de `step9_collision_sim.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step9_collision_sim`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


def printed_margin_of(params: "CollisionSimParams") -> float:
    """Marge appliquee a la matiere deja deposee."""
    return params.safety_margin if params.printed_part_margin is None else params.printed_part_margin


@dataclass
class CollisionSimParams:
    path_sampling_step: float
    safety_margin: float
    # Marge propre a la MATIERE DEJA DEPOSEE, mesuree a la surface du cordon
    # (rayon de buse). None = `safety_margin`. NEGATIVE = recouvrement toleré :
    # le modele de tete est exact, donc la pointe peut effleurer ou mordre un peu
    # la couche sans que ce soit une collision (decision de l'utilisateur,
    # 2026-10-27). Mesure : sur deux G-code (un cartesien plat, un Benchy polaire)
    # les contacts a marge 0 recouvrent la matiere de 0,18 mm au plus. Le plateau
    # garde `safety_margin` : la pointe ne descend jamais sous lui.
    printed_part_margin: Optional[float] = None
    # Non utilise par step9 lui-meme - porte pour le contrat CollisionSimParams
    # (CLAUDE.md), consomme par l'etape 10 qui boucle sur step9.
    max_adjustment_iterations: int = 3
    # Taille de cellule de la grille spatiale des capsules imprimees, en
    # multiple de nozzle_diameter (voir _PrintedGeometryIndex).
    grid_cell_size_multiplier: float = 4.0
    # Pas d'echantillonnage des ARETES du modele de tete (_surface_sample_points),
    # en mm - DELIBEREMENT distinct de nozzle_diameter (largeur de buse,
    # rayon des capsules imprimees : un concept different, la taille du
    # CORDON, pas la taille de la TETE). Les deux etaient confondus a
    # l'origine (limite #21) ; mesure concretement en executant main.py sur
    # le modele de tete reel (aetheris-4.glb, aretes de plusieurs dizaines
    # de mm) avec nozzle_diameter=0.4 par defaut : 17734 points de test par
    # point de trajectoire echantillonne, simulation totale de l'ordre de
    # la dizaine de minutes sur un modele trivial - inutilisable. Defaut
    # choisi a l'echelle typique d'une tete d'impression (quelques mm),
    # PAS a l'echelle d'un cordon d'extrusion.
    toolhead_sample_spacing: float = 5.0
    # Rayon (mm, repere local du modele de tete) autour de la pointe de buse
    # (origine locale) a l'interieur duquel les points d'echantillonnage sont
    # EXCLUS des tests de collision. La pointe de buse est par construction
    # toujours au contact du cordon qu'elle vient de deposer et proche du
    # plateau sur les premieres couches - sans cette zone de garde, ce seul
    # point declenche une collision a CHAQUE point de trajectoire teste,
    # sans exception (mesure sur l_shape_overhang.stl : exactement 2
    # evenements par point echantillonne, printed_part + platform - pas un
    # motif de collision reel, un artefact structurel du test). A l'echelle
    # de la buse/marge de securite (quelques mm), PAS de nozzle_diameter
    # (echelle du cordon, meme piege que toolhead_sample_spacing ci-dessus,
    # limites #21/#26). Defaut non nul : ce champ corrige un bug de faux
    # positif systematique, pas une fonctionnalite optionnelle a activer.
    nozzle_contact_clearance: float = 3.0
    # Controle RAPIDE (opt-in, porte de la recherche le 2026-09-30) : la matiere
    # imprimee est densifiee le long des cordons (pas `fast_densify_step`) dans
    # UN arbre k-d construit une fois, chaque point densifie portant l'indice
    # du point de trajectoire qui le depose ; un point teste ne voit que la
    # matiere deposee AVANT lui (meme invariant que _grow_printed_index_to).
    # Remplace la boucle de capsules, dont le debit s'effondre avec la matiere
    # accumulee (limite #27) : 30 min -> ~1 min sur le Benchy en recherche.
    # Densification : la distance a l'axe du cordon est SURESTIMEE d'au plus
    # un demi-pas (0,125 mm) -- a couvrir par safety_margin (1 mm).
    fast_kdtree: bool = False
    fast_densify_step: float = 0.25
