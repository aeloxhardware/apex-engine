"""Etape 5 : paramètres et erreur d'extraction des couches.

Découpé de `step5_layers.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step5_layers`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


class LayerExtractionError(Exception):
    """Levee si le ScalarField ne peut pas etre converti en couches (gradient manquant, champ degenere)."""


@dataclass
class LayerParams:
    layer_height_base: float = 0.2
    # tolerance de fusion des sommets du patch triangule, avant extraction du
    # contour, relative a l'etendue du TetMesh - evite des micro-trous de
    # bord dus au bruit flottant plutot qu'a une vraie discontinuite
    # geometrique (cf. Decision 5 du plan).
    weld_tolerance_factor: float = 1e-6
    # Scinde une couche en plusieurs Layer distincts (meme isovalue,
    # region_ids et local_frame recalcules par sous-ensemble) quand ses
    # boucles de contour forment plusieurs composantes spatialement
    # disjointes ET non imbriquees (ex: jonction en T ou en etoile, ou une
    # isovaleur traverse simultanement le pied et l'extremite d'un bras
    # eloigne - CLAUDE.md). None (defaut) = DESACTIVE, comportement
    # historique inchange (un seul Layer par isovaleur, local_frame moyenne
    # sur toutes les boucles) - meme discipline que
    # direction_blend_radius=0.0 (step4), near_platform_resolution=None
    # (step1), max_region_merge_angle_deg=None (step3). Valeur = facteur
    # multiplicatif de layer_height_base utilise comme tolerance de
    # distance (cf. _cluster_contours) : deux boucles non imbriquees
    # distantes de moins de split_disjoint_clusters * layer_height_base
    # restent dans le meme cluster (evite une discontinuite artificielle
    # juste avant qu'elles ne fusionnent geometriquement a une isovaleur
    # voisine).
    #
    # Active AUSSI, automatiquement, le tri spatial par continuite (DFS sur
    # un graphe de proximite entre clusters d'isovaleurs successives,
    # respectant STRICTEMENT le rang topologique de region comme tri
    # primaire - jamais traverse, cf. Decision 8) : sans ce tri, deux
    # clusters partageant le meme rang ET les memes region_ids (typiquement
    # pied et pointe de bras d'un T, meme region topologique au sens step3
    # malgre une topologie SOLIDE separee) restent ordonnes uniquement par
    # isovaleur croissante, ce qui les fait alterner a chaque isovaleur
    # (pied->pointe->pied->pointe...) et produit un faisceau de trajets a
    # vide entre les deux a chaque couche - mesure directement sur
    # res/test_models/l-triple-axis-overhanged-shape.stl (CLAUDE.md). Pas de
    # parametre separe : ce tri n'a de sens/d'effet que si le clustering
    # produit deja plusieurs Layer par isovaleur, et degenere exactement
    # vers le tri historique (rank, isovalue) en leur absence - un flag
    # independant ajouterait une combinatoire sans benefice reel.
    split_disjoint_clusters: Optional[float] = None

    # --- TRI PAR ISOVALEUR SEULE (2026-10-22) ------------------------------
    # Le tri historique `(rang de region, isovaleur)` suppose que l'isovaleur
    # ne dit pas a elle seule l'ordre de croissance, et s'en remet au rang de
    # step3 pour le dire. C'est l'inverse quand le champ EST une distance au
    # plateau a travers la matiere (`use_heat_geodesic_field`) : phi est alors
    # la progression elle-meme, et le rang de region la casse. Mesure sur
    # death.stl (champ de la chaleur) : les couches de rang 2 (phi 23,6 a 39,6)
    # passaient avant celles de rang 3 a phi 8,3, donc des pieces imprimees
    # avant la matiere qui les relie au reste (correlation ordre/isovaleur 0,78).
    #
    # True (defaut) : le champ est la distance a travers la matiere. False =
    # tri (rang de region, isovaleur).
    isovalue_order: bool = True
