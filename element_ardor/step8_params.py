"""Etape 8 : paramètres et erreur de résolution de trajectoire.

Découpé de `step8_machine_toolpath.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step8_machine_toolpath`.
"""

from __future__ import annotations

from dataclasses import dataclass


class ToolpathResolutionError(Exception):
    """
    Levee uniquement pour une entree structurellement invalide (toolpath_2d
    manquant, schedule pointant vers un segment/couche inexistant) - JAMAIS
    pour une orientation inatteignable, qui est un cas attendu couvert par
    Step8Diagnostic, pas une erreur bloquante.
    """


@dataclass
class Step8Params:
    # Repli par grille porte en natif (kinematics.cpp, GridSolver). False rend le repli Python,
    # point par point (mesure d'equivalence, tests). R-030 point 3 : plus de variable d'environnement.
    native_grid_fallback: bool = True
    # Nombre d'echantillons par axe pour la recherche du plus proche voisin
    # (composed_rotation_angles_and_directions) - cout total = ~samples**N.
    ik_samples_per_axis: int = 10
    ik_refine_iterations: int = 8
    # Multiplicateur du pas initial d'affinement, PAS une constante fixe en
    # degres : le pas de depart = ik_refine_step_multiplier * (360 /
    # ik_samples_per_axis), a l'echelle de l'ecart entre deux points voisins
    # de la grille NN (jusqu'a ~36deg pour samples_per_axis=10). Un pas fixe
    # trop petit laisse l'affinement piege pres du candidat NN au lieu de
    # converger vers le vrai optimum (mesure : residu ~15-20deg avec un pas
    # fixe de 2deg ; <0.1deg avec ce pas derive - voir tests).
    ik_refine_step_multiplier: float = 0.5
    # Distance en-deca de laquelle la position courante est consideree SUR
    # l'axe de rotation (singularite : angle indetermine, herite du point
    # precedent plutot que recalcule).
    singularity_radius_epsilon: float = 1e-3
    # Feedrate des trajets a vide inseres aux frontieres de polyligne. 0 =>
    # reutilise le feedrate du point suivant (pas de vitesse magique par
    # defaut non documentee).
    travel_feedrate: float = 0.0
    default_feedrate: float = 1000.0
    # Au-dela de ce seuil, l'erreur residuelle d'un point est jugee assez
    # grande pour meriter un Step8Diagnostic (signal vers l'etape 10).
    max_unreachable_angular_error_deg: float = 5.0
    # Troncature des angles par les limites de vitesse/acceleration du profil
    # (_clamp_rate). DESACTIVEE par defaut depuis le 2026-09-28 : le firmware
    # de la machine COORDONNE les axes (confirme par l'utilisateur) -- quand un
    # axe atteint sa limite, c'est tout le mouvement qui ralentit, l'angle
    # reste exact. Tronquer l'angle fait prendre a la tete un retard qu'elle
    # n'aura pas (mesure en recherche : jusqu'a 62 deg d'ecart a la normale sur
    # la coque du Benchy). Le bon modele est un ralentissement de l'avance, a
    # faire a l'export G-code (feuille de route). Les limites restent dans le
    # profil machine ; True retrouve l'ancien comportement.
    clamp_angle_rates: bool = False
    # Filtre MEDIAN (en nombre de points, 0 = aucun) des axes d'orientation
    # NON couples a la position, par polyligne (2026-09-30, mesure en
    # recherche). La normale de couche est constante par tetraedre : pres de
    # l'axe du plateau, C etant impose, B oscille 0 <-> ~14 deg d'un point a
    # l'autre (5,9 % des segments du corps de Squirtle), et le firmware, qui
    # coordonne les axes, ralentit tout le mouvement a chaque inversion
    # (x3,7 sur le temps du corps). 9 points : -41 % de temps. Opt-in : a
    # appliquer AVANT tout controle de collision (step9 valide les angles
    # filtres). L'ecart residuel a la direction visee est recalcule.
    # ACTIVE PAR DEFAUT a 9 (2026-10-06) : mesure sur les 4 modeles de res/test_models/ avec la
    # normale par point active (limite #69/#70) : temps d'impression -19 / -20 / -31 / -21 %,
    # collisions de la vraie tete en baisse (death 137 -> 118, turkish 1613 -> 1406), matiere
    # sans appui inchangee ; cout : ecart median tete/normale +0,2 a +0,4 deg. 0 = aucun filtre.
    free_axis_median_window: int = 9
    # Normale PAR POINT (phase 3.1, session du 2026-10-01). Par defaut, la cible
    # de l'IK est la normale du plan MOYEN de la couche, une seule pour toute la
    # couche. Sur la machine actuelle (C impose par la position, limite #59), une
    # normale inclinee fixe n'est dans le plan radial que d'une poignee de
    # points : mesure, 95 % des points de l_shape en `Step8Diagnostic`. Vrai =
    # la cible devient le gradient LOCAL du champ en chaque point, fourni par
    # l'orchestrateur (`run(..., direction_field=...)`) ; sans champ fourni,
    # comportement historique. Meme choix que la recherche (full_step89 :
    # « normale du corps = gradient de phi »).
    # ACTIVE PAR DEFAUT (2026-10-06, decision de l'utilisateur) et reglable : mesure
    # sur les 4 modeles, ecart median tete/normale divise par 2 a 3, mais temps
    # d'impression +65 a +130 % (voir CLAUDE.md, limite #69). False = normale du
    # plan moyen de la couche (historique).
    per_point_normal: bool = True
    # Variante continue de la normale par point (gradient moyenne aux noeuds,
    # interpole en barycentrique) ; sans effet si per_point_normal est faux.
    per_point_normal_nodal: bool = False
