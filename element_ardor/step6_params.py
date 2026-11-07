"""Etape 6 : paramètres, motifs de remplissage et erreur de génération.

Découpé de `step6_toolpath2d.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step6_toolpath2d`.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ToolpathGenerationError(Exception):
    """Levee si une couche ne peut produire aucune trajectoire exploitable."""


class InfillPattern(Enum):
    # lignes paralleles independantes, meme orientation a toutes les couches
    # (comportement historique)
    LINES = "lines"
    # deux balayages croises a 90 deg dans CHAQUE couche, chacun a la moitie
    # de la densite
    GRID = "grid"
    # lignes croisees-alternees : un zigzag continu par zone, oriente selon
    # l'axe u du repere local une couche sur deux et selon l'axe v l'autre
    # (voir `_infill_angle_for_layer`)
    CROSSED_ZIGZAG = "crossed_zigzag"
    # spirale d'Archimede autour de l'AXE DU PLATEAU, oscillante : r(theta)
    # avec theta qui ne recule jamais (le plateau tourne toujours dans le meme
    # sens, les pics sont des mouvements radiaux). MEME ligne mediane a toutes
    # les couches (elles se superposent) ; seules les asperites different
    # d'une couche a l'autre (phase opposee, tirages differents).
    OSCILLATING_SPIRAL = "oscillating_spiral"


@dataclass
class Step6Params:
    # R-030 pt 3 : False = shapely / couche par couche ; sans effet si le noyau natif est inactif (#82).
    native_fill: bool = True
    native_layer: bool = True
    nozzle_width: float = 0.4
    wall_count: int = 2
    infill_density: float = 0.2          # dans (0, 1]
    infill_pattern: InfillPattern = InfillPattern.CROSSED_ZIGZAG
    # seuil de PRESENCE (pas de majorite) : proportion minimale de points
    # d'une boucle de contour dont la face de surface la plus proche est
    # is_bridge, au-dela de laquelle la boucle entiere est traitee comme
    # trajectoire de pont. Un pont n'est qu'une petite portion locale du
    # perimetre total d'une couche (CLAUDE.md : zones is_bridge "qui
    # retombent dans cette couche", pas la couche entiere) - un seuil de
    # majorite ne peut structurellement jamais se declencher (mesure sur la
    # fixture pont : ratio reel entre 0.06 et 0.21 sur toutes les couches
    # concernees, jamais proche de 0.5). D'ou un defaut bas.
    bridge_contour_ratio_threshold: float = 0.05
    # parois supplementaires resserrees sur le perimetre d'un contour pont,
    # en plus de wall_count, avant de generer l'infill balaye - lecture
    # directe de CLAUDE.md "augmente la densite d'ancrage sur les bords"
    # (renfort de perimetre, pas de densite d'infill - cf. plan step6).
    bridge_anchor_extra_walls: int = 2
    # aire residuelle (relative a l'aire du contour d'origine) sous laquelle
    # un offset est considere comme ayant fait disparaitre la zone.
    degenerate_area_ratio: float = 1e-6
    # Exposant de la ponderation inverse-distance utilisee pour donner un
    # offset hors-plan aux points GENERES par shapely (parois n>=1, infill),
    # interpole depuis les offsets REELS connus sur les boucles de contour.
    # None = DESACTIVE (offset nul pour ces points, comportement historique) -
    # meme discipline opt-in que direction_blend_radius=0.0 (step4) ou
    # max_region_merge_angle_deg=None (step3).
    #
    # Sans ce mecanisme, ces points restent colles au plan local MOYEN de la
    # couche alors que la couche reelle s'en ecarte (mesure sur le modele en
    # T : jusqu'a 0.845 mm d'ecart), et peuvent donc tomber hors du solide -
    # residu que la limite #31 de CLAUDE.md documentait comme hors de son
    # perimetre. Le contour n=0 n'est JAMAIS concerne : il garde ses offsets
    # exacts, mesures, jamais interpoles.
    out_of_plane_interpolation_power: float | None = None

    # Decoupe chaque ANNEAU FERME (paroi, contour original) en troncons
    # monotones, chacun emis comme une POLYLIGNE SEPAREE et parcouru dans le
    # sens montant - voir Decision 5 plus bas. False = DESACTIVE (anneau
    # entier en une polyligne, comportement historique), meme discipline
    # opt-in que out_of_plane_interpolation_power ci-dessus. Necessite ce
    # dernier actif pour avoir un effet reel : sans lui l'offset hors-plan
    # des points generes est nul, chaque anneau est deja plat et ne produit
    # qu'un seul troncon.
    split_rings_into_runs: bool = False

    # Etend le decoupage en troncons monotones (ci-dessus) aux polylignes
    # OUVERTES - lignes d'infill balayees, cf. `_generate_infill` - et plus
    # seulement aux anneaux fermes. False = DESACTIVE (comportement
    # historique), meme discipline opt-in que les champs ci-dessus.
    #
    # Mesure sur death.stl (donnees capturees, 5 couches les plus touchees) :
    # les polylignes ouvertes portent 211 des 243 paires descendantes d'une
    # couche (87%) - le decoupage historique, limite aux anneaux, ne couvrait
    # que ~5% des paires de la trajectoire. Voir `_split_rings_in_toolpath`.
    split_open_polylines: bool = False

    # Continuite des parois (2026-10-01). Deux segments n'etaient JAMAIS
    # extrudes, mesure sur le pipeline par defaut :
    #   - la jonction entre deux troncons monotones (`_split_ring_into_runs`) :
    #     troncons a indices DISJOINTS, donc le segment qui les relie devient
    #     un trajet a vide - 7,7 % de la longueur deposee perdue sur death.stl,
    #     15,7 % sur l_shape_overhang.stl (des trous dans les parois) ;
    #   - la fermeture du contour d'origine n=0 et de ses trous, emis OUVERTS
    #     (0,1 % et 0,85 %).
    # True : chaque troncon reprend le dernier sommet du precedent (la
    # jonction reste monotone, elle prend le sens du troncon suivant), et les
    # contours n=0 recoivent leur point de fermeture. False = historique.
    continuous_rings: bool = False

    # Amplitude minimale (mm) d'un changement de sens pour COUPER un troncon
    # monotone. Le defaut 1e-12 reproduit l'ancien comportement (toute
    # oscillation, meme numerique, coupe). Mesure sur un cylindre plein
    # (couche plate) : l'anneau exterieur a 257 points et une etendue d'offset
    # de 3,8 um, soit 122 changements de sens -> 160 polylignes par couche au
    # lieu de 4, chacune avec son trajet a vide. Une oscillation sous cette
    # amplitude est du bruit du marching tetrahedra, bien en dessous de
    # l'epaisseur d'un cordon : elle ne cree aucun appui manquant.
    run_split_min_amplitude: float = 1e-12

    # --- OSCILLATING_SPIRAL (opt-in) --------------------------------------
    # Point de l'axe du plateau (x, y, z, repere machine) : centre de la
    # spirale. Source de verite unique = MachineProfile.platform_origin
    # (renseigne par main) ; None = centre de la zone (la monotonie de theta
    # n'est alors garantie qu'autour de ce centre, pas autour de l'axe).
    spiral_axis_origin: tuple[float, float, float] | None = None
    # "Seuil d'Alternance" (mm) : ecart radial des pics. None = 0,6 x la
    # limite (pas - buse) / 2. Toujours plafonne a (pas - buse) / 2 : deux
    # tours voisins ne se touchent jamais.
    spiral_alternation_threshold: float | None = None
    # Distance (mm, le long de la trajectoire) entre deux pics. None = 4 buses.
    spiral_peak_spacing: float | None = None
    # "Seuil de Randomisation" dans [0, 1] : 0 = oscillation reguliere,
    # 1 = amplitude, espacement, forme (pic net / courbe) et pics purement
    # radiaux tires au hasard. Graine fixe : resultat reproductible.
    spiral_randomization_threshold: float = 0.5
    spiral_seed: int = 0
