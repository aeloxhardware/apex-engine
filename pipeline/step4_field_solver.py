"""
pipeline/step4_field_solver.py

Etape 4 du pipeline : assemblage et resolution du champ scalaire continu sur
le TetMesh, dont les isosurfaces (etape 5) definiront les couches courbes.

Depend de TetMesh + MeshData (step1), OverhangMap (step2),
RegionDependencyGraph (step3), et des Constraint[] envoyees dynamiquement
par l'UI (points d'attache, contraintes esthetiques/structurelles).

Le systeme assemble est A.phi = b, symetrique defini positif, somme de :
  - un terme de lissage (rigidite FEM P1) ;
  - un terme de surplomb (le gradient doit suivre la direction d'impression
    de chaque region, pres des faces a risque) ;
  - un ancrage automatique du plateau (leve la singularite du systeme) ;
  - les contraintes utilisateur, assemblees par core/constraints.py.

Structure en trois fonctions (cf. plan, decision 6) : `build_static_system()`
assemble tout ce qui NE depend PAS des contraintes,
`solve_with_constraints()` ajoute les contraintes et resout, et `run()` est
l'enveloppe des deux pour un slicing complet en une passe. Ce decoupage
permet a l'orchestrateur de re-resoudre rapidement apres l'edition d'un
point d'attache, sans reassembler la partie structurelle ni relancer
step1/2/3.

Sortie : ScalarField
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.csgraph as spgraph
import scipy.sparse.linalg as spla
import scipy.spatial

from core import _cpp
from core.constraints import ConstraintAssemblyContext, assemble_all_constraints
from core.mesh_types import (
    Constraint,
    UpwardGrowthDiagnostic,
    MeshData,
    OverhangMap,
    RegionDependencyGraph,
    ScalarField,
    SmoothnessModulationDiagnostic,
    SupportFeedbackDiagnostic,
    TetMesh,
)
from core.spatial import KDTree


class FieldSolverError(Exception):
    """Levee si le systeme assemble n'est pas exploitable ou si le solveur ne converge pas."""


@dataclass
class FieldSolverParams:
    laplacian_smoothness_weight: float = 1.0
    gradient_overhang_weight: float = 5.0
    resolution: float = 1.5
    max_solver_iterations: int = 2000
    convergence_tolerance: float = 1e-8

    # Poids de base du terme "le gradient suit la direction d'impression de
    # la region", applique a TOUS les tetraedres. Sans lui, le champ ne
    # serait pilote que pres des surplombs et resterait indetermine (donc
    # nul, l'ancrage plateau valant 0) dans le reste du volume : une piece
    # sans surplomb n'aurait aucune couche. C'est ce terme qui produit le
    # slicing plan classique par defaut ; gradient_overhang_weight ne fait
    # que le renforcer la ou le risque est eleve.
    gradient_baseline_weight: float = 1.0

    # Nombre de sauts d'adjacence tetraedre-tetraedre sur lesquels propager
    # l'influence d'une face a risque. 1 = seulement les tetraedres au
    # contact direct de la face frontiere concernee.
    overhang_influence_depth: int = 1

    # Ancrage automatique des regions racines (contact plateau) a
    # platform_anchor_value. Necessaire pour lever la singularite du
    # systeme (un Laplacien seul a un noyau constant : phi = c est toujours
    # solution) quand aucune Constraint[] utilisateur n'est fournie.
    enable_platform_anchor: bool = True
    platform_anchor_weight: float = 1e6
    platform_anchor_value: float = 0.0

    # Verticale machine globale, utilisee UNIQUEMENT par
    # _blend_tet_directions (le plateau y est traite comme une source
    # d'influence directionnelle "verticale locale", au meme titre qu'une
    # region). Meme valeur par defaut que SegmentationParams.up_vector
    # (step3) - deux dataclasses distinctes par construction du projet
    # (chaque etape porte ses propres parametres), mais meme convention de
    # repere, a fournir coherente par l'appelant (main.py) si une autre
    # valeur est utilisee en step3.
    up_vector: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        if self.up_vector is None:
            self.up_vector = np.array([0.0, 0.0, 1.0])

    # Seuil (nombre de noeuds du TetMesh) en-dessous duquel un solveur
    # direct (spsolve) est utilise plutot qu'un solveur iteratif (cg).
    # Le remplissage d'une factorisation directe 3D croit vite, d'ou un
    # seuil bien inferieur a ce qu'on tolererait en 2D.
    direct_solve_node_threshold: int = 30000
    use_ilu_preconditioner: bool = False

    # Rayon d'influence du melange spatial des directions de gradient (voir
    # _blend_tet_directions). 0.0 = DESACTIVE, comportement historique de
    # plus-proche-voisin strict conserve (winner-take-all, cf. CLAUDE.md) -
    # changement de comportement uniquement sur opt-in explicite, meme
    # discipline que direction_deviation_weight de step3 (limites #31/#32).
    # Au-dela de ce rayon, une source (plateau ou region) n'a plus aucune
    # influence sur un tetraedre - borne le cout de calcul, pas juste le
    # poids (une coupure nette, pas une decroissance infinie a calculer
    # pour toutes les sources sur tous les tetraedres).
    direction_blend_radius: float = 0.0

    # Exposant de la ponderation inverse-distance (w_i = 1/(d_i+eps)^p).
    # Plus eleve = transition plus abrupte pres de chaque source ; plus bas
    # = melange plus etale sur tout le rayon d'influence.
    direction_blend_power: float = 2.0

    # Seuil de degenerescence du melange : si la norme de la moyenne
    # ponderee des directions (avant renormalisation) tombe sous ce seuil
    # (sources a directions quasi-opposees et poids proches, annulation
    # vectorielle), _blend_tet_directions NE normalise PAS un vecteur
    # instable - elle replie sur la direction de la source au poids
    # maximal pour ce tetraedre precis (equivalent local au comportement
    # winner-take-all, uniquement la ou le melange continu serait
    # mathematiquement degenere). Une moyenne de vecteurs unitaires a une
    # norme maximale de 1.0 (directions alignees) et tend vers 0 quand
    # elles s'annulent - 0.1 est un seuil de depart, a calibrer.
    direction_blend_min_norm: float = 0.1

    # --- Redressement des directions de croissance descendantes -----------
    # Un gradient de champ descendant empile les couches du HAUT VERS LE BAS :
    # chaque couche atterrit dans le vide. Sur une machine dont le plateau ne
    # tourne qu'autour de l'axe de la gravite, la piece n'est JAMAIS retournee
    # - "le bas" dans le repere piece est toujours -up_vector, donc une
    # direction de croissance descendante est physiquement irrealisable, pas
    # seulement sous-optimale.
    #
    # Mecanisme : les directions cibles par tetraedre (sortie de
    # _blend_tet_directions) dont la composante selon up est negative sont
    # RETOURNEES (d -> -d) avant l'assemblage du terme de direction.
    #
    # Pourquoi retourner la cible et non penaliser le gradient obtenu : le
    # terme de direction pousse grad(phi) VERS tet_directions. Majorer le
    # poids d'un tetraedre descendant renforce donc sa cible - et si cette
    # cible pointe elle-meme vers le bas, on amplifie exactement la descente
    # qu'on veut supprimer. Mesure directe sur death.stl : une penalite
    # ponderee iterative fait passer le volume descendant de 67.6% a 70.3%
    # (elle AGGRAVE), la ou le redressement de cible le fait tomber a 0.1%.
    #
    # Pourquoi c'est geometriquement sans perte : phi et -phi ont exactement
    # les MEMES isosurfaces, donc les memes couches. Retourner une direction
    # cible ne deforme aucune couche - cela n'inverse que le SENS de parcours
    # de la pile, c'est-a-dire precisement le defaut corrige.
    #
    # False = DESACTIVE, comportement historique strictement inchange - meme
    # discipline que direction_blend_radius=0.0.
    enforce_upward_growth: bool = False

    # Tolerance au bruit numerique : une composante a -1e-12 n'est pas une
    # vraie descente. Jamais un test contre 0.0 exact.
    upward_growth_epsilon: float = 1e-9

    # --- BORNE DE PLONGEE, alternative a l'interdit binaire ci-dessus -----
    # `enforce_upward_growth` repose sur une hypothese FAUSSE pour cette
    # machine : que l'on ne peut pas deposer vers le bas. La tete couvre
    # B ∈ [0°, 185°], donc elle le peut (voir MachineProfile). L'interdit
    # empechait d'exprimer le seul ordre correct pour une piece suspendue -
    # la faire croitre vers le bas depuis son attache.
    #
    # Le vrai besoin n'est pas d'interdire la descente mais de ne pas passer
    # trop de temps a l'extreme de la course de la tete. Cette borne fixe donc
    # l'angle maximal SOUS L'HORIZONTALE de la direction de croissance :
    #   0°  -> rien ne plonge sous l'horizontale
    #   45° -> plongee jusqu'a 45° (B ≈ 135°)
    #   90° -> descente verticale libre (B ≈ 180°)
    # La direction bornee conserve son AZIMUT : on l'aplatit, on ne la
    # retourne pas.
    #
    # PRECEDENCE : quand ce parametre est renseigne, il REMPLACE
    # `enforce_upward_growth` (les deux repondent a la meme question et se
    # contrediraient). None (defaut) = borne inactive, comportement
    # historique strictement inchange.
    max_downward_growth_angle_deg: Optional[float] = None

    # --- Direction cible issue de la DISTANCE GEODESIQUE au plateau -------
    # Remplace les directions par REGION (step3, melangees spatialement par
    # _blend_tet_directions) par le gradient normalise de la distance
    # geodesique "plus court chemin jusqu'au plateau EN RESTANT DANS LA
    # MATIERE", calculee par tetraedre.
    #
    # Pourquoi c'est mieux fonde qu'une direction par region : une region
    # recoit UNE direction pour l'ensemble de ses faces, choisie en
    # minimisant le surplomb residuel - critere qui, pour un plafond de
    # cavite, elit une direction descendante physiquement absurde (il
    # faudrait commencer a imprimer en l'air, cf. limite #51). La geodesique
    # repond a la question physique directe : "par ou la matiere rejoint-elle
    # le plateau, et dans quel sens s'en eloigne-t-elle ICI" - une reponse
    # par tetraedre, jamais un compromis moyenne sur toute une region.
    #
    # Le champ n'est PAS remplace par la geodesique brute : celle-ci est trop
    # irreguliere pour servir d'isosurfaces (Dijkstra suit les aretes en
    # zigzag - mesure sur death.stl : ||grad|| median 4.81 au lieu de 1.0
    # attendu pour une distance, dispersion 0.87 contre 0.29 pour le champ
    # FEM). Elle sert donc uniquement de DIRECTION CIBLE ; le lissage FEM et
    # l'ancrage plateau, inchanges, apportent la regularite.
    #
    # Mesure qui a motive ce mecanisme (death.stl, tetraedres du dome) : le
    # champ FEM actuel ne fait JAMAIS pointer le depot vers le bas (0% des
    # tetraedres, composante verticale mediane +0.955) la ou le capot
    # redescend pourtant ; la geodesique y pointe vers le bas sur 28% d'entre
    # eux (mediane +0.299).
    #
    # False (defaut) = DESACTIVE, comportement historique strictement
    # inchange - meme discipline que direction_blend_radius=0.0.
    use_geodesic_directions: bool = False

    # --- Terme geodesique COMME SOURCE, pas comme substitution ------------
    # `use_geodesic_directions` REMPLACE tet_directions par le gradient brut
    # de la distance geodesique. Mesure sur real_overhang.stl (resolution
    # 0.4) : ce gradient a des POINTS DE SELLE - deux chemins geodesiques
    # quasi equidistants font basculer sa direction de ~177 deg d'un
    # tetraedre a son voisin, alors que ||grad(distance)|| reste 0.999 des
    # deux cotes (donc CE N'EST PAS un artefact de maillage, ni un bruit
    # haute-frequence lissable : un lissage spatial par sauts topologiques,
    # teste jusqu'a 8 sauts, laisse le MAXIMUM inchange a ~170 deg). Impose
    # comme cible RIGIDE, ce bruit survit au Laplacien et froisse les
    # couches : residu de planeite mesure x12-19 pire que la reference
    # (0.050mm median -> 0.646mm, 1.1% -> 20.7% relatif a l'etendue de la
    # couche).
    #
    # Ce parametre ajoute a la place un terme SOURCE : la meme energie de
    # correspondance de gradient que `_assemble_gradient_direction_term`
    # (E_e = w_e * V_e * ||grad(phi) - d||^2), avec d = direction geodesique
    # et un poids INDEPENDANT du poids de direction principal. Le systeme
    # arbitre alors lui-meme entre lissage (Laplacien, deja dans A) et
    # fidelite a la geodesique - au lieu de figer une cible bruitee puis
    # d'esperer que le lissage la corrige apres coup, ce qui a ete mesure
    # comme ne fonctionnant pas.
    #
    # Le poids porte sur la MEME variable que `gradient_overhang_weight`
    # (tet_weights, poids du terme de direction PRINCIPAL) : les deux sont
    # additifs dans le systeme, jamais exclusifs. 0.0 (defaut) = terme
    # absent, comportement historique strictement inchange - meme discipline
    # que direction_blend_radius=0.0.
    #
    # NON MONOTONE en ce poids jusqu'a preuve du contraire (piege deja paye
    # avec anchor_growth_weight) : ne jamais interpoler entre deux valeurs
    # mesurees.
    geodesic_source_weight: float = 0.0

    # Seuil (en degres) au-dela duquel un ecart angulaire entre deux
    # tetraedres voisins est considere comme un POINT DE SELLE plutot qu'une
    # variation normale du champ (voir `_geodesic_saddle_mask`). Les
    # tetraedres marques sont EXCLUS du terme source geodesique (poids nul
    # localement), le systeme retombant sur le terme de direction principal
    # deja lisse a cet endroit.
    #
    # 100.0 par defaut : mesure sur real_overhang.stl, seules 2.6% des paires
    # de voisins depassent ce seuil, et leur qualite de maillage n'est pas
    # significativement pire que la mediane (r=+0.36 avec l'ecart angulaire :
    # correlation faible, ce n'est pas un defaut de maillage a corriger en
    # amont). Une variation "normale" mesuree meme aux abords d'une jonction
    # reelle reste sous ce seuil (jonction pilier/bras : p90 a 56 deg).
    geodesic_saddle_threshold_deg: float = 100.0

    # --- Lissage VECTORIEL iteratif du champ geodesique --------------------
    # `use_geodesic_directions` impose le gradient geodesique BRUT comme
    # cible rigide - mesure comme degradant (limite CLAUDE.md #52, x2.2 sur
    # death.stl) precisement a cause des points de selle ci-dessus : un
    # lissage LINEAIRE composante par composante (`_diffuse_tet_directions`,
    # utilise pour `transition_distance_mm`) resout un systeme (I+t*L)x=x0
    # SEPAREMENT sur x/y/z puis renormalise UNE SEULE FOIS a la fin - a un
    # point de selle (deux voisins a ~180 deg), cela moyenne deux vecteurs
    # quasi opposes en un vecteur proche de zero, renormalise ensuite
    # n'importe comment. C'est exactement le mecanisme qui laisse le lissage
    # a 8 sauts topologiques deja teste bloque a ~170 deg de max (meme
    # limite #52).
    #
    # Ce parametre active un lissage NON-LINEAIRE ITERATIF a la place :
    # `_iterative_vector_smoothing`, qui renormalise a CHAQUE iteration (pas
    # une fois a la fin). Mesure sur death.stl (resolution 0.8, seuil de
    # discordance = fraction de paires de voisins a plus de 100 deg) :
    #   brut ..................... 3.79%
    #   10 iterations, pas 0.5 ... 0.03%
    #   50 iterations, pas 0.5 ... 0.00%
    # Le champ se debloque progressivement parce qu'a chaque passe, la
    # moyenne locale est deja legerement infléchie par les voisins plus
    # loin du point de selle (qui n'ont pas ce probleme), au lieu de figer
    # l'artefact de moyenne lineaire en un seul calcul.
    #
    # S'applique UNIQUEMENT quand `use_geodesic_directions=True` (jamais de
    # sens seul) : lisse le champ geodesique AVANT qu'il ne remplace
    # tet_directions, jamais apres.
    #
    # 0 (defaut) = DESACTIVE, comportement historique inchange - meme
    # discipline que use_geodesic_directions=False. NON MESURE sur step9b -
    # bon indicateur intermediaire (discordance angulaire) necessaire mais
    # pas suffisant, piege deja paye 5 fois dans ce projet (voir CLAUDE.md,
    # "aucun indicateur intermediaire ne predit step9b").
    geodesic_smoothing_iterations: int = 0

    # Pas de la relaxation a chaque iteration (0 = fige, 1 = remplace
    # entierement par la moyenne des voisins). 0.5 mesure comme un bon
    # compromis convergence/stabilite sur death.stl - voir docstring
    # ci-dessus. NON MONOTONE en general (piege anchor_growth_weight) : ne
    # jamais interpoler entre deux valeurs mesurees.
    geodesic_smoothing_step: float = 0.5

    # --- BOUCLE DE RETROACTION DE SURPLOMB --------------------------------
    # Tout le pipeline evalue le surplomb en BOUCLE OUVERTE : step2 mesure
    # l'angle des normales contre la verticale MACHINE globale, une seule
    # fois, AVANT que le champ n'existe ; ce risque fige pondere ensuite le
    # terme de gradient ci-dessus. Le solveur n'est donc jamais penalise pour
    # les surplombs que ses PROPRES couches creent - step9b les mesure, mais
    # apres coup et sans rien reinjecter.
    #
    # Cette boucle ferme le circuit : resoudre phi, re-evaluer le surplomb
    # contre `grad(phi)` (la direction de depot REELLEMENT produite, pas la
    # verticale machine), faire pivoter les directions cibles fautives juste
    # assez pour sortir du seuil, re-resoudre. Le critere est exactement celui
    # de step2 (`angle(normale, -direction) < seuil`), applique avec la
    # verticale LOCALE du champ au lieu de la verticale globale.
    #
    # PLAFOND d'iterations, pas un nombre fixe : la boucle s'arrete d'elle-meme
    # des que le gain d'une iteration a la suivante devient negligeable (voir
    # `support_feedback_min_improvement`) ou que le surplomb est resorbe. Borne
    # dur a SUPPORT_FEEDBACK_ITERATION_CAP quelle que soit la valeur demandee.
    #
    # 0 (defaut) = DESACTIVE, comportement historique strictement inchange.
    # Chaque iteration coute une resolution complete du systeme.
    support_feedback_iterations: int = 0

    # Directions cibles ATTEIGNABLES par la machine actuelle (CLAUDE.md
    # limite #59, phase 3.1 de la session du 2026-10-01). C est impose par la
    # position du point : la tete ne peut incliner la normale de couche QUE
    # dans le plan (radial, vertical) du point, et seulement VERS L'EXTERIEUR
    # (+R, B >= 0). Toute composante tangentielle est inatteignable sans axe D.
    # Mesure avant ce mecanisme : `step8_diagnostics` = 95 % des points sur
    # l_shape, ~ tous sur death.
    # Renseigne (origine du plateau), chaque direction cible est projetee sur
    # ce demi-plan : composante tangentielle retiree, composante radiale
    # negative ramenee a 0. None (defaut) = DESACTIVE, historique inchange.
    radial_tilt_origin: Optional[np.ndarray] = None

    # Critere d'arret dynamique : gain RELATIF minimal d'une iteration a la
    # suivante sur la fraction de volume en surplomb auto-infligee. En dessous,
    # la boucle s'arrete et se declare `stalled` (jamais `converged` - voir
    # SupportFeedbackDiagnostic, s'arreter faute de progres n'est pas avoir
    # resolu le probleme).
    #
    # Relatif et non absolu : le residu part de ~26% du volume et descend d'un
    # ordre de grandeur ; un seuil absolu serait soit inatteignable au debut,
    # soit atteint immediatement a la fin.
    support_feedback_min_improvement: float = 0.02

    # Seuil de surplomb de la boucle, en degres. DOIT refleter
    # `OverhangParams.overhang_angle_threshold_deg` (meme critere geometrique,
    # applique a une verticale differente) : step4 ne recoit pas OverhangParams,
    # donc la valeur est duplique ici faute de source unique. Si le seuil de
    # step2 est change, changer celui-ci - piege deja rencontre avec
    # `platform_z` (deux sources de verite qui divergent en silence).
    support_feedback_threshold_deg: float = 45.0

    # Sur-ponderation appliquee aux tetraedres dont la cible a ete corrigee :
    # sans elle, la cible corrigee serait noyee par le lissage et la
    # correction n'aurait aucun effet mesurable.
    support_feedback_weight: float = 4.0

    # Cible corrigee = direction GEODESIQUE LOCALE du tetraedre, au lieu de la
    # rotation minimale par defaut.
    #
    # PUREMENT PAR TETRAEDRE, jamais par zone ni par region - c'est la
    # condition meme pour que ce soit utile. Une direction moyennee sur toute
    # une region ne peut pas epouser sa geometrie : mesure sur
    # `real_overhang.stl`, la direction unique attribuee au bras pointe a
    # CONTRESENS de son axe (alignement -0.45), ce qui fait qu'une seule couche
    # enjambe le bras entier et se retrouve deposee a 10mm de toute matiere.
    # La geodesique du MEME endroit pointe dans le bon sens (+0.39).
    #
    # Les deux quantites croisees ici sont locales et recalculees a chaque
    # passe : le surplomb auto-inflige (par tetraedre, dependant du champ
    # courant) et le gradient geodesique (par tetraedre, "par ou la matiere
    # rejoint le plateau ICI"). Aucune moyenne globale n'intervient.
    #
    # La geodesique n'est retenue que si elle RESOUT effectivement le surplomb
    # a ce tetraedre ; sinon on retombe sur la rotation minimale. Elle n'est
    # donc jamais imposee sur la foi de sa seule provenance - ce qui evite de
    # reproduire l'echec de `use_geodesic_directions`, qui l'appliquait PARTOUT
    # et ecrasait des directions deja correctes (mesure : x3,7 de degradation).
    #
    # False (defaut) = rotation minimale, comportement inchange.
    support_feedback_use_geodesic_target: bool = False

    # --- MODULATION SPATIALE DU LISSAGE -----------------------------------
    # `laplacian_smoothness_weight` ci-dessus est UN scalaire, le meme partout.
    # C'est la derniere echelle a laquelle une moyenne globale subsiste dans le
    # pipeline, et elle a un cout geometrique precis : sur une piece qui se
    # ramifie, le lissage uniforme force le champ a ETALER de part et d'autre
    # du branchement le compromis entre deux directions de depot incompatibles
    # (monter dans le pilier / partir a l'horizontale dans le bras). D'ou des
    # couches qui se courbent AVANT l'angle, au lieu de rester planes puis de
    # casser franchement.
    #
    # Ce parametre abaisse LOCALEMENT le poids de lissage la ou la piece se
    # ramifie, par tetraedre :
    #
    #     poids_t = max(1 / (1 + strength * indicateur_t), floor)
    #
    # Le champ reste UNIQUE - c'est ce qui distingue ce mecanisme d'un
    # decoupage en branches avec un champ par branche. Consequences directes :
    # les couches restent les lignes de niveau d'une seule fonction, donc ne
    # peuvent toujours pas se croiser ni laisser de vide a une interface ; les
    # contraintes utilisateur continuent de diffuser dans toute la piece ; et
    # step5 a step8 continuent de recevoir des isovaleurs globalement
    # monotones, sans aucune modification.
    #
    # 0.0 (defaut) = DESACTIVE. Le poids vaut alors exactement 1.0 partout et
    # `_assemble_stiffness_matrix` est appele sans ponderation du tout : le
    # comportement historique est retrouve a l'identique, pas a une tolerance
    # pres.
    #
    # NON MONOTONE en `strength` jusqu'a preuve du contraire : ne jamais
    # interpoler entre deux valeurs mesurees (piege deja paye avec
    # `anchor_growth_weight`, ou w=0.25 AGGRAVAIT avant que w>=0.5 n'ameliore).
    smoothness_modulation_strength: float = 0.0

    # Plancher du poids de lissage module. Necessaire parce que l'indicateur
    # |div(X)| n'est PAS borne : mesure jusqu'a 1.8 en mediane de tranche au
    # bout du pilier suspendu de real_overhang.stl. Sans plancher, le contraste
    # entre zone raide et zone souple serait dicte par la valeur extreme de
    # l'indicateur sur le modele, donc imprevisible d'un modele a l'autre.
    #
    # Il borne aussi le CONDITIONNEMENT : le rapport de raideur entre elements
    # ne peut pas depasser 1/floor. C'est ce qui rend le surcout du solveur
    # iteratif previsible sur les gros modeles (death.stl, 76220 noeuds, seul
    # des quatre modeles de reference a passer au-dessus de
    # `direct_solve_node_threshold`).
    #
    # 0.02 = contraste maximal de 50x.
    smoothness_modulation_floor: float = 0.02

    # --- DISTANCE DE TRANSITION entre branches (en mm) --------------------
    # Repond a une question qu'aucun autre parametre ne sait poser : sur
    # QUELLE DISTANCE les couches doivent-elles commencer a se deformer pour
    # prendre l'orientation de la branche suivante ?
    #
    # Agit par DIFFUSION du champ de DIRECTIONS CIBLES (`_diffuse_tet_directions`),
    # pas sur le lissage du champ scalaire. Ce choix vient d'un audit apres
    # une premiere version MESUREE SANS EFFET, et le resultat negatif vaut
    # d'etre garde : moduler le poids du terme laplacien en fonction de la
    # distance a la jonction change bien le champ (nombre de couches
    # 284 -> 333 sur real_overhang.stl) mais laisse le profil d'inclinaison
    # des couches RIGOUREUSEMENT inchange (ecarts de 1-2 deg, du bruit, de
    # 2 mm a 40 mm de transition). Raison : le terme de lissage pilote
    # l'ESPACEMENT des isosurfaces, tandis que leur ORIENTATION est dictee
    # par le terme de direction - dont la cible etait constante par branche,
    # donc insensible a toute modulation du lissage.
    #
    # La diffusion agit au bon endroit : a l'interieur d'une branche les
    # directions cibles sont constantes, donc moyenner un tetraedre avec ses
    # voisins ne le change pas. Le mecanisme est donc AUTO-LIMITE aux
    # frontieres entre branches, la ou les directions different reellement -
    # sans avoir besoin d'identifier explicitement "la branche d'en face" ni
    # de masque geometrique.
    #
    #   PETIT (ex: 2 mm)  -> cassure tardive et nette, les couches gardent
    #       l'orientation de leur branche jusqu'au dernier moment ;
    #   GRAND (ex: 20 mm) -> transition precoce et progressive, les couches
    #       commencent a tourner bien avant la jonction.
    #
    # C'est exactement le controle demande : "plus le seuil est eleve, plus
    # les couches doivent commencer a changer d'orientation tot, sur une
    # longue distance".
    #
    # None (defaut) = DESACTIVE, comportement historique strictement
    # inchange - meme discipline que direction_blend_radius=0.0.
    transition_distance_mm: Optional[float] = None

    # Temps de diffusion de la methode de la chaleur, en multiples du CARRE de
    # la longueur d'arete mediane du maillage (la formulation est dimensionnee
    # en longueur au carre, une constante absolue n'aurait aucun sens d'un
    # modele a l'autre).
    #
    # Trop petit : la chaleur ne diffuse pas assez et la direction du flux est
    # dominee par le bruit numerique. Trop grand : elle lisse par-dessus la
    # geometrie et la distance perd ses details.
    #
    # MESURE (real_overhang.stl, resolution 0.4, 26590 tetraedres, modele haut
    # de 49 aretes) - l'ancien defaut de 1.0, repris de la litterature et
    # explicitement laisse "A MESURER", etait DEGENERE :
    #
    #   facteur   ||grad d|| p50   tetraedres a gradient nul   r vs Dijkstra
    #      1.0        0.000              11166 / 26590 (42%)       +0.846
    #     10.0        0.998                  0 / 26590             +0.986
    #    100.0        0.999                  0 / 26590             +0.985
    #   1000.0        1.000                  0 / 26590             +0.985
    #
    # A 1.0, t = 0.56 pour une arete mediane de 0.75 : la chaleur ne diffuse
    # que sur UNE arete alors que le modele en fait 49 de haut. Le champ
    # SATURE (valeur constante 22.174 partout au-dessus de z~9.7) et la
    # direction de son gradient n'est plus que du bruit numerique au-dela du
    # front de saturation. La recommandation "de l'ordre de 1" de la
    # litterature suppose un maillage bien plus fin rapporte a l'etendue du
    # domaine.
    #
    # 100.0 retenu : les valeurs sont convergees de 10 a 1000 (55.30 contre
    # 55.36, soit 0.1% d'ecart), donc la marge est prise sans rien couter.
    heat_method_time_factor: float = 100.0

    # Le champ scalaire EST la distance geodesique lisse au plateau, au lieu
    # d'etre la solution du systeme lissage/surplomb/ancrage habituel.
    #
    # Consequence recherchee : les couches deviennent des coquilles
    # EQUIDISTANTES a travers la matiere. Leur etalement geodesique tombe a
    # zero par construction - il est aujourd'hui de 11,70mm en mediane pour
    # une epaisseur de 0,26mm, ce qui signifie qu'un meme contour relie des
    # points tres proches du plateau et d'autres tres lointains a travers le
    # solide. C'est ce melange, et non l'ordonnancement, qui empeche depuis le
    # debut de deposer la matiere dans un ordre coherent.
    #
    # LIMITES ASSUMEES, signalees plutot que masquees :
    #  - les Constraint[] utilisateur ne sont PAS honorees dans ce mode (le
    #    champ ne resout plus leur systeme) ;
    #  - le terme de surplomb et l'orientation par region de step3 n'ont plus
    #    d'effet : l'orientation des couches devient purement geometrique.
    # Ces deux points sont exactement ce qui reste a arbitrer avant d'en faire
    # autre chose qu'un mode opt-in.
    #
    # False (defaut) = DESACTIVE, comportement historique inchange.
    use_heat_geodesic_field: bool = False


@dataclass
class FieldSolverStaticSystem:
    """
    Partie du systeme lineaire qui NE depend PAS des Constraint[] : elle
    reste valide tant que tet_mesh, mesh_data, overhang_map,
    region_dependency_graph et les poids/seuils de FieldSolverParams sont
    inchanges.

    L'orchestrateur la construit une fois apres step3, la met en cache, et
    la reutilise a chaque edition de point d'attache via
    `solve_with_constraints()`. Si un parametre non lie aux contraintes
    change (poids de lissage, de surplomb, d'ancrage, seuils solveur), le
    cache doit etre invalide et `build_static_system()` rappele.
    """
    A_structural: sp.csr_matrix
    b_structural: np.ndarray
    ctx: ConstraintAssemblyContext
    # Residu du redressement des directions descendantes, mesure pendant
    # l'assemblage structurel (avant toute Constraint[] utilisateur).
    # None = mecanisme desactive (enforce_upward_growth=False).
    upward_growth_diagnostic: Optional[UpwardGrowthDiagnostic] = None
    # Residu de la boucle de retroaction de surplomb. None = mecanisme
    # desactive (support_feedback_iterations=0).
    support_feedback_diagnostic: Optional[SupportFeedbackDiagnostic] = None
    # Ce qu'a fait la modulation spatiale du lissage. None = mecanisme
    # desactive (smoothness_modulation_strength=0.0).
    smoothness_modulation_diagnostic: Optional[SmoothnessModulationDiagnostic] = None
    # Solution du systeme STRUCTUREL SEUL (sans aucune Constraint[]), quand
    # elle a deja ete calculee pour le diagnostic de croissance montante.
    # `solve_with_constraints` la reutilise quand les contraintes n'apportent
    # rien - voir la note de reutilisation dans cette fonction. None = pas de
    # sonde calculee (diagnostic desactive).
    probe_values: Optional[np.ndarray] = None


# ---------------------------------------------------------------------------
# Decision 1 - Matrice de rigidite FEM P1 (Laplacien discret geometrique)
# ---------------------------------------------------------------------------

def _compute_shape_gradients(tet_mesh: TetMesh) -> tuple[np.ndarray, np.ndarray]:
    """
    Pour chaque tetraedre, calcule G (3,4) - le gradient constant de chacune
    des 4 fonctions de forme barycentriques - et son volume.

    Methode standard P1 : les colonnes de l'inverse transpose de la matrice
    des aretes (x1-x0, x2-x0, x3-x0) donnent les gradients des fonctions de
    forme 1..3 ; le gradient de la fonction de forme 0 est l'oppose de leur
    somme (partition de l'unite : la somme des 4 fonctions vaut 1 partout,
    donc la somme de leurs gradients est nulle).

    Retourne (shape_gradients (T,3,4), volumes (T,)). Les tetraedres
    degeneres (volume ~ 0) ont un gradient nul et un volume nul : leur
    contribution disparait naturellement des assemblages.
    """
    tet_count = tet_mesh.tet_count

    if _cpp.enabled():
        # Formule des cofacteurs en C++ plutot que `np.linalg.inv` (LAPACK) :
        # deux algorithmes differents pour le meme inverse 3x3, donc des
        # arrondis differents. Ecart RELATIF mesure sur death.stl : 7,2e-14,
        # soit la precision machine - voir la note d'equivalence de step4,
        # ou la tolerance remplace l'egalite bit a bit (limite #56, pas 7).
        gradients, volumes = _cpp.module().compute_shape_gradients(
            np.ascontiguousarray(tet_mesh.tet_vertices, dtype=np.float64),
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64),
        )
        return gradients, volumes

    shape_gradients = np.zeros((tet_count, 3, 4))
    volumes = np.zeros(tet_count)

    vertices_per_tet = tet_mesh.tet_vertices[tet_mesh.tetrahedra]  # (T,4,3)
    edge_matrices = np.transpose(
        vertices_per_tet[:, 1:, :] - vertices_per_tet[:, 0:1, :], (0, 2, 1)
    )  # (T,3,3), colonnes = x1-x0, x2-x0, x3-x0

    determinants = np.linalg.det(edge_matrices)
    signed_volumes = determinants / 6.0
    valid = np.abs(signed_volumes) > 1e-14

    if np.any(valid):
        inverses = np.linalg.inv(edge_matrices[valid])          # (V,3,3)
        grad_1_3 = np.transpose(inverses, (0, 2, 1))            # (V,3,3)
        grad_0 = -grad_1_3.sum(axis=2, keepdims=True)           # (V,3,1)
        shape_gradients[valid] = np.concatenate([grad_0, grad_1_3], axis=2)
        volumes[valid] = np.abs(signed_volumes[valid])

    return shape_gradients, volumes


class _FemPattern:
    """
    Motif de sparsite CSR de l'assemblage FEM, partageable entre matrices.

    Les deux matrices de step4 (rigidite et terme de direction) ont
    EXACTEMENT le meme motif : celui de l'adjacence noeud-noeud du maillage,
    qui ne depend que de la topologie, jamais des poids. Le chemin scipy le
    reconstruisait pourtant deux fois, en triant a chaque fois les 6,2 M
    entrees COO (16 par tetraedre) puis en sommant les doublons.

    Construit explicitement par l'appelant et passe aux deux assemblages -
    jamais memorise sur le TetMesh : un cache invisible survivrait a une
    modification de `tetrahedra` et produirait une matrice silencieusement
    fausse.
    """

    __slots__ = ("_native", "indptr", "indices", "node_count")

    def __init__(self, tet_mesh: TetMesh):
        self.node_count = tet_mesh.node_count
        self._native = _cpp.module().FemSparsity(
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64), tet_mesh.node_count
        )
        self.indptr = self._native.indptr()
        self.indices = self._native.indices()

    def gram(
        self,
        tet_mesh: TetMesh,
        shape_gradients: np.ndarray,
        scale: np.ndarray,
        active: Optional[np.ndarray],
    ) -> sp.csr_matrix:
        """Matrice `sum_e scale_e * G_e^T G_e`, assemblee DIRECTEMENT en CSR."""
        data = self._native.assemble_gram(
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64),
            np.ascontiguousarray(shape_gradients, dtype=np.float64),
            np.ascontiguousarray(scale, dtype=np.float64),
            None if active is None else np.ascontiguousarray(active, dtype=np.int32),
        )
        return sp.csr_matrix(
            (data, self.indices, self.indptr), shape=(self.node_count, self.node_count)
        )


def _assemble_stiffness_matrix(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    tet_weights: Optional[np.ndarray] = None,
    pattern: Optional["_FemPattern"] = None,
) -> sp.csr_matrix:
    """
    Assemble la matrice de rigidite FEM P1 globale (Laplacien discret
    geometriquement correct) : K_e = w_e * V_e * G_e^T G_e par element, sommee
    aux indices de sommets globaux. Symetrique semi-definie positive, de
    noyau les champs constants.

    `tet_weights` (T,) module le poids de lissage PAR ELEMENT (voir
    `FieldSolverParams.smoothness_modulation_strength`). `None` = poids
    uniforme, strictement identique au comportement historique - et pas
    seulement a une tolerance pres : c'est la MEME expression, sans
    multiplication supplementaire.

    Le MOTIF DE SPARSITE ne depend pas des poids, seulement de la topologie
    du maillage. Une factorisation directe (`spsolve`) coute donc exactement
    la meme chose avec ou sans modulation ; seul un solveur iteratif peut
    voir son conditionnement change.
    """
    scale = volumes if tet_weights is None else volumes * tet_weights

    if _cpp.enabled():
        # Assemblage DIRECT en CSR : ni matrices locales materialisees, ni
        # COO intermediaire de 6,2 M entrees a trier puis dedupliquer.
        pattern = pattern if pattern is not None else _FemPattern(tet_mesh)
        return pattern.gram(tet_mesh, shape_gradients, scale, None)

    local_matrices = scale[:, None, None] * np.einsum(
        "tki,tkj->tij", shape_gradients, shape_gradients
    )  # (T,4,4)
    return _scatter_local_matrices(tet_mesh, np.arange(tet_mesh.tet_count), local_matrices)


def _scatter_local_matrices(
    tet_mesh: TetMesh, tet_indices: np.ndarray, local_matrices: np.ndarray
) -> sp.csr_matrix:
    """
    Disperse des matrices locales (n,4,4) vers la matrice creuse globale,
    aux indices de sommets des tetraedres `tet_indices`. Les doublons
    d'indices sont sommes par coo_matrix, ce qui realise l'assemblage FEM.
    """
    node_count = tet_mesh.node_count
    if len(tet_indices) == 0:
        return sp.csr_matrix((node_count, node_count))

    tets = tet_mesh.tetrahedra[tet_indices]                     # (n,4)
    rows = np.repeat(tets, 4, axis=1).ravel()                   # (n*16,)
    cols = np.tile(tets, (1, 4)).ravel()                        # (n*16,)
    data = local_matrices.reshape(-1)

    return sp.coo_matrix(
        (data, (rows, cols)), shape=(node_count, node_count)
    ).tocsr()


def _compute_per_tet_gradient(
    tet_mesh: TetMesh, shape_gradients: np.ndarray, values: np.ndarray
) -> np.ndarray:
    """
    Gradient constant par tetraedre du champ P1 defini par `values` (une
    valeur par noeud) : grad(phi)|_e = G_e @ phi_e. Representation exacte
    d'un champ lineaire par element, pas une approximation.
    """
    values_per_tet = values[tet_mesh.tetrahedra]                # (T,4)
    return np.einsum("tij,tj->ti", shape_gradients, values_per_tet)


# ---------------------------------------------------------------------------
# Decision 2 - Pont face de surface -> tetraedre volumique
# ---------------------------------------------------------------------------

def _tet_face_table(tet_mesh: TetMesh) -> dict[tuple[int, int, int], list[int]]:
    """
    Table face (triplet trie de noeuds) -> liste des tetraedres proprietaires.
    Une face interne a 2 proprietaires, une face frontiere un seul.
    """
    face_to_tets: dict[tuple[int, int, int], list[int]] = {}
    local_faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))

    for t_idx, tet in enumerate(tet_mesh.tetrahedra):
        for local_face in local_faces:
            key = tuple(sorted(int(tet[i]) for i in local_face))
            face_to_tets.setdefault(key, []).append(t_idx)

    return face_to_tets


def _map_boundary_faces_to_surface(
    tet_mesh: TetMesh, mesh_data: MeshData
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Identifie les faces frontiere du TetMesh (celles n'appartenant qu'a un
    seul tetraedre) et, pour chacune, retrouve la face du mesh de surface
    d'origine la plus proche via MeshData.kdtree (construit a l'etape 1 sur
    les centres de faces).

    Un mapping par index direct serait fragile : tetgen est appele sans le
    switch `Y`, donc il peut inserer des points de Steiner sur la frontiere
    et celle-ci n'est pas identique noeud-a-noeud au mesh source.

    Retourne (boundary_face_nodes (K,3), boundary_face_tets (K,),
    nearest_surface_face_id (K,), face_pairs (M,2) | None).

    `face_pairs` est un SOUS-PRODUIT gratuit du chemin natif : la meme passe
    topologique produit les faces frontiere ET les paires de tetraedres
    partageant une face interne. Le rendre evite a `_build_tet_tet_adjacency`
    de refaire integralement ce calcul (0,2 s mesurees sur death.stl @ 0.8).
    `None` sur le chemin Python, ou l'appelant retombe sur le calcul separe.
    """
    face_pairs = None
    if _cpp.enabled():
        nodes, owners, pairs = _cpp.module().tet_face_topology(
            np.asarray(tet_mesh.tetrahedra, dtype=np.int64)
        )
        boundary_face_nodes = nodes.astype(np.int64)
        boundary_face_tets = owners.astype(np.int64)
        face_pairs = pairs
        if len(boundary_face_tets) == 0:
            empty_int = np.zeros((0, 3), dtype=np.int64)
            return empty_int, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), face_pairs
    else:
        face_to_tets = _tet_face_table(tet_mesh)

        boundary_nodes = []
        boundary_tets = []
        for key, owners in face_to_tets.items():
            if len(owners) == 1:
                boundary_nodes.append(key)
                boundary_tets.append(owners[0])

        if not boundary_nodes:
            empty_int = np.zeros((0, 3), dtype=np.int64)
            return empty_int, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), None

        boundary_face_nodes = np.array(boundary_nodes, dtype=np.int64)
        boundary_face_tets = np.array(boundary_tets, dtype=np.int64)

    centroids = tet_mesh.tet_vertices[boundary_face_nodes].mean(axis=1)
    _, nearest_surface_face_id = mesh_data.kdtree.query(centroids)

    return (
        boundary_face_nodes,
        boundary_face_tets,
        np.asarray(nearest_surface_face_id, dtype=np.int64),
        face_pairs,
    )


def _build_tet_tet_adjacency(
    tet_mesh: TetMesh, face_pairs: Optional[np.ndarray] = None
) -> sp.csr_matrix:
    """
    Matrice d'adjacence tetraedre-tetraedre (deux tetraedres sont adjacents
    s'ils partagent une face). Distincte de TetMesh.tet_adjacency, qui est
    une adjacence noeud-noeud : ici on propage un masque defini PAR
    TETRAEDRE, pas par noeud.

    `face_pairs` evite de RECALCULER la topologie des faces : elle est deja
    produite par `_map_boundary_faces_to_surface`, qui n'en garde que les
    faces frontiere et jette les paires internes - exactement ce dont cette
    fonction a besoin. Les deux appels coutaient 0,2 s chacun sur death.stl
    @ 0.8 pour le meme calcul, chacun jetant la moitie du resultat.
    """
    tet_count = tet_mesh.tet_count

    if face_pairs is not None:
        pairs = face_pairs
        if len(pairs) == 0:
            return sp.csr_matrix((tet_count, tet_count))
        rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
        cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
        return sp.csr_matrix(
            (np.ones(len(rows)), (rows, cols)), shape=(tet_count, tet_count)
        )

    if _cpp.enabled():
        _nodes, _owners, pairs = _cpp.module().tet_face_topology(
            np.asarray(tet_mesh.tetrahedra, dtype=np.int64)
        )
        if len(pairs) == 0:
            return sp.csr_matrix((tet_count, tet_count))
        # Symetrisation explicite, comme le Python : chaque paire donne les
        # deux sens.
        rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
        cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
        data = np.ones(len(rows))
        return sp.csr_matrix((data, (rows, cols)), shape=(tet_count, tet_count))

    face_to_tets = _tet_face_table(tet_mesh)

    rows, cols = [], []
    for owners in face_to_tets.values():
        if len(owners) == 2:
            t1, t2 = owners
            rows.extend([t1, t2])
            cols.extend([t2, t1])

    if not rows:
        return sp.csr_matrix((tet_count, tet_count))

    data = np.ones(len(rows))
    return sp.csr_matrix((data, (rows, cols)), shape=(tet_count, tet_count))


# ---------------------------------------------------------------------------
# Decision 2bis - Melange spatial des directions de gradient (remplace le
# plus-proche-voisin strict par une moyenne ponderee par distance inverse,
# sur le plateau + chaque region comme sources d'influence de meme nature -
# resout a la fois la transition plateau->region ET les jonctions entre
# regions avec un seul mecanisme. Voir CLAUDE.md pour le diagnostic complet
# (transition brusque mesuree sur l_shape_overhang.stl, fragmentation de
# contour 1->4 boucles mesuree sur un modele tripode a 3 branches).
# ---------------------------------------------------------------------------

def _blend_tet_directions(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    region_dependency_graph: RegionDependencyGraph,
    up_vector: np.ndarray,
    params: FieldSolverParams,
) -> np.ndarray:
    """
    Direction de gradient cible par tetraedre (T,3), par moyenne ponderee
    par distance inverse sur les sources d'influence proches (le plateau,
    direction `up_vector`, + chaque region, sa `target_direction`) - plutot
    que le plus-proche-voisin strict (une seule source "gagne" par
    tetraedre, transition discontinue a la frontiere).

    `params.direction_blend_radius <= 0.0` : DESACTIVE, retombe sur le
    comportement historique (plus-proche-voisin global sur TOUTES les
    faces de surface, une seule requete KDTree) - aucun changement de
    comportement par defaut.

    Sinon, pour chaque tetraedre, chaque source dans un rayon
    `direction_blend_radius` recoit un poids `1/(distance+eps)^power`
    (`params.direction_blend_power`) ; la direction resultante est la
    moyenne ponderee normalisee des directions sources. Cas degenere
    (sources quasi-opposees a poids proches, moyenne de norme quasi nulle)
    : repli sur la source au poids maximal plutot que de normaliser un
    vecteur instable (`params.direction_blend_min_norm`).
    """
    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    tet_count = tet_mesh.tet_count

    if params.direction_blend_radius <= 0.0:
        _, nearest_face_per_tet = mesh_data.kdtree.query(tet_centers)
        nearest_face_per_tet = np.asarray(nearest_face_per_tet, dtype=np.int64)

        face_direction = np.zeros((mesh_data.face_count, 3))
        face_direction[:, 2] = 1.0  # defaut : verticale, si une face n'a pas de region
        for region in region_dependency_graph.regions.values():
            face_direction[region.face_ids] = region.target_direction

        return face_direction[nearest_face_per_tet]

    # Sources : une par region (KDTree sur les centres de SES SEULES faces,
    # pas toute la surface - distance a CETTE region precisement), plus le
    # plateau (plan, distance verticale directe - pas besoin de KDTree).
    source_directions: list[np.ndarray] = []
    source_distances: list[np.ndarray] = []

    face_centers = mesh_data.vertices[mesh_data.faces].mean(axis=1)
    region_points: list[np.ndarray] = []
    for region in region_dependency_graph.regions.values():
        if len(region.face_ids) == 0:
            continue
        region_points.append(face_centers[region.face_ids])
        source_directions.append(region.target_direction)

    platform_level = _compute_platform_level(mesh_data, region_dependency_graph)
    if platform_level is not None:
        norm = np.linalg.norm(up_vector)
        source_directions.append(up_vector / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0]))

    if source_directions and _cpp.enabled():
        # Chemin natif : meme calcul, mais chaque source ECARTE d'emblee les
        # tetraedres hors de son AABB elargie du rayon. Mesure sur death.stl :
        # 69,2 % des paires (tetraedre, region) sont ainsi evitees, et seules
        # 21,4 % etaient reellement dans le rayon - le Python calculait 8,5 M
        # de distances dont la fenetre de Hann annulait la majorite.
        # L'elagage est EXACT : une source hors boite est a plus de `radius`,
        # donc de poids nul ; les tetraedres qu'aucune source n'atteint sont
        # repris sans elagage pour retrouver l'argmin du repli.
        scale = float(np.ptp(tet_mesh.tet_vertices, axis=0).max()) or 1.0
        offsets = np.zeros(len(region_points) + 1, dtype=np.int32)
        if region_points:
            offsets[1:] = np.cumsum([len(p) for p in region_points])
        packed = (
            np.ascontiguousarray(np.vstack(region_points), dtype=np.float64)
            if region_points
            else np.zeros((0, 3), dtype=np.float64)
        )
        return _cpp.module().blend_tet_directions(
            np.ascontiguousarray(tet_centers, dtype=np.float64),
            packed,
            offsets,
            np.ascontiguousarray(np.array(source_directions, dtype=np.float64)),
            platform_level,
            params.direction_blend_radius,
            params.direction_blend_power,
            params.direction_blend_min_norm,
            1e-6 * scale,
        )

    for points in region_points:
        distances, _ = KDTree(points).query(tet_centers)
        source_distances.append(np.asarray(distances, dtype=np.float64))
    if platform_level is not None:
        source_distances.append(np.abs(tet_centers[:, 2] - platform_level))

    if not source_directions:
        # aucune source (ni region, ni plateau identifiable) : comportement
        # de repli identique au defaut historique (verticale partout).
        return np.tile(np.array([0.0, 0.0, 1.0]), (tet_count, 1))

    distances = np.stack(source_distances, axis=1)  # (T, S)
    directions = np.stack(source_directions, axis=0)  # (S, 3)

    scale = float(np.ptp(tet_mesh.tet_vertices, axis=0).max()) or 1.0
    epsilon = 1e-6 * scale

    # Fenetrage LISSE jusqu'a 0 au bord du rayon (pas une coupure nette) :
    # une coupure binaire (distance <= radius) fait basculer l'ENSEMBLE des
    # sources actives de facon discontinue des qu'une source franchit le
    # rayon, produisant un effet d'escalier mesure visuellement (CLAUDE.md)
    # meme si chaque poids individuel decroit "en douceur" par ailleurs -
    # le probleme n'est pas la magnitude du poids, c'est le changement brusque
    # du NOMBRE de sources sommees. Fenetre cosinus (Hann) : vaut 1 au centre
    # (distance=0), decroit continument a 0 pile a `radius` (valeur ET pente
    # nulles au bord - C1, pas seulement C0), et est nulle au-dela - jamais
    # de saut de valeur NI de pente quand une source entre/sort du rayon.
    within_radius = distances <= params.direction_blend_radius
    # tetraedre sans aucune source dans le rayon : elargit implicitement a
    # la source la plus proche seule, pour ne jamais laisser une direction
    # indeterminee (meme esprit que le defaut "verticale" du chemin
    # winner-take-all quand une face n'a pas de region).
    no_source_in_radius = ~within_radius.any(axis=1)
    nearest_source = np.argmin(distances, axis=1)
    within_radius[no_source_in_radius, nearest_source[no_source_in_radius]] = True
    # rayon effectif par tetraedre : le rayon nominal, sauf si la source la
    # plus proche est deja hors du rayon nominal (repli ci-dessus) - alors
    # la fenetre est centree sur cette distance precise pour ne pas creer
    # une fenetre nulle qui diviserait par zero.
    effective_radius = np.maximum(
        params.direction_blend_radius, distances[np.arange(tet_count), nearest_source] + epsilon
    )

    inverse_distance_weight = 1.0 / np.power(distances + epsilon, params.direction_blend_power)
    smooth_window = np.where(
        within_radius,
        0.5 * (1.0 + np.cos(np.pi * np.minimum(distances / effective_radius[:, None], 1.0))),
        0.0,
    )
    weights = inverse_distance_weight * smooth_window
    weights_sum = weights.sum(axis=1, keepdims=True)
    weights_normalized = weights / weights_sum  # weights_sum > 0 garanti par within_radius ci-dessus

    blended = weights_normalized @ directions  # (T,3)
    blended_norm = np.linalg.norm(blended, axis=1)

    degenerate = blended_norm < params.direction_blend_min_norm
    result = np.zeros((tet_count, 3))
    result[~degenerate] = blended[~degenerate] / blended_norm[~degenerate, None]

    if np.any(degenerate):
        fallback_source = np.argmax(weights[degenerate], axis=1)
        result[degenerate] = directions[fallback_source]

    return result


# ---------------------------------------------------------------------------
# Decision 3 - Terme de direction de gradient (base + sur-ponderation surplomb)
# ---------------------------------------------------------------------------

def _spread_risk_influence(
    tet_tet_adjacency: sp.csr_matrix, seed_tets: np.ndarray, hops: int
) -> np.ndarray:
    """
    Propage un masque booleen de tetraedres "a risque" sur `hops` sauts
    d'adjacence tetraedre-tetraedre, en partant de `seed_tets`.
    """
    tet_count = tet_tet_adjacency.shape[0]
    mask = np.zeros(tet_count, dtype=bool)
    if len(seed_tets) == 0:
        return mask
    mask[seed_tets] = True

    frontier = mask.copy()
    for _ in range(max(0, hops)):
        reached = (tet_tet_adjacency @ frontier.astype(np.float64)) > 0
        new_frontier = reached & (~mask)
        if not np.any(new_frontier):
            break
        mask |= new_frontier
        frontier = new_frontier

    return mask


def _assemble_gradient_direction_term(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    active_tets: np.ndarray,
    tet_directions: np.ndarray,
    tet_weights: np.ndarray,
    pattern: Optional["_FemPattern"] = None,
) -> tuple[sp.csr_matrix, np.ndarray]:
    """
    Aligne le gradient du champ sur la direction d'impression de la region,
    via l'energie de correspondance de gradient
    E_e = w_e * V_e * ||grad(phi) - d||^2. Sa minimisation donne :
        A_local = w_e * V_e * G_e^T G_e
        b_local = w_e * V_e * G_e^T d

    Le terme en b est la force motrice du champ : une penalite purement
    homogene (sur la seule composante perpendiculaire, G_e^T (I - d d^T) G_e)
    n'injecterait rien dans le second membre et laisserait le systeme
    s'effondrer sur phi = 0 partout en l'absence de contrainte de valeur.

    Applique a tous les tetraedres de `active_tets`, avec un poids par
    tetraedre : poids de base partout (le champ doit definir des couches
    dans tout le volume, pas seulement pres des surplombs), majore autour
    des zones a risque.
    """
    node_count = tet_mesh.node_count
    active = np.where(active_tets & (tet_weights > 0.0) & (volumes > 0.0))[0]
    if len(active) == 0:
        return sp.csr_matrix((node_count, node_count)), np.zeros(node_count)

    gradients = shape_gradients[active]                          # (n,3,4)
    directions = tet_directions[active]                          # (n,3)
    norms = np.linalg.norm(directions, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    directions = directions / norms

    scale = tet_weights[active] * volumes[active]                # (n,)

    if _cpp.enabled():
        pattern = pattern if pattern is not None else _FemPattern(tet_mesh)
        A = pattern.gram(tet_mesh, gradients, scale, active)
        b = _cpp.module().assemble_gradient_rhs(
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64),
            np.ascontiguousarray(gradients, dtype=np.float64),
            np.ascontiguousarray(scale, dtype=np.float64),
            np.ascontiguousarray(directions, dtype=np.float64),
            np.ascontiguousarray(active, dtype=np.int32),
            node_count,
        )
        return A, b

    local_matrices = scale[:, None, None] * np.einsum(
        "tki,tkj->tij", gradients, gradients
    )                                                            # (n,4,4)
    local_rhs = scale[:, None] * np.einsum("tki,tk->ti", gradients, directions)  # (n,4)

    A = _scatter_local_matrices(tet_mesh, active, local_matrices)

    b = np.zeros(node_count)
    np.add.at(b, tet_mesh.tetrahedra[active].ravel(), local_rhs.ravel())

    return A, b


# ---------------------------------------------------------------------------
# Decision 3bis - Redressement des directions de croissance descendantes
#
# Diagnostic a l'origine de ce mecanisme (mesure sur res/test_models/death.stl) :
# 67.6% du VOLUME du TetMesh avait grad(phi).up < 0, ce qui empile les couches
# du haut vers le bas. Consequence mesuree sur la trajectoire machine finale :
# 11.25% de la matiere extrudee etait deposee sans aucune matiere deja
# imprimee en dessous. Cause : _blend_tet_directions melange les directions
# cibles de regions voisines par proximite spatiale sans garantir que le
# melange reste montant, et step3 choisit parfois lui-meme une direction de
# region descendante pour une region qui n'est pas une racine plateau
# (mesure : 4 regions sur 22, dont la plus grande du modele).
#
# MECANISME ECARTE, mesure avant de l'ecarter : une penalite quadratique
# reponderee iterativement (majorer le poids des tetraedres descendants, puis
# re-resoudre) AGGRAVE le probleme - 67.6% -> 70.3% sur death.stl. La raison
# est structurelle, pas un defaut de reglage : _assemble_gradient_direction_term
# pousse grad(phi) VERS tet_directions, donc majorer le poids d'un tetraedre
# descendant RENFORCE sa cible. Or 61.6% du volume a une cible qui pointe
# elle-meme vers le bas : la penalite amplifiait exactement ce qu'elle devait
# supprimer. Le levier correct est la CIBLE, pas le poids.
#
# MECANISME RETENU : retourner (d -> -d) les directions cibles descendantes
# avant l'assemblage. Sans perte geometrique - phi et -phi ont exactement les
# memes isosurfaces, donc les memes couches ; seul le SENS de parcours de la
# pile change, c'est-a-dire precisement le defaut corrige. Mesure sur
# death.stl : 67.6% -> 0.1% du volume, a cout de calcul inchange (aucune
# resolution supplementaire, c'est une transformation de tableau).
# ---------------------------------------------------------------------------

def _assemble_lumped_mass_matrix(
    tet_mesh: TetMesh, volumes: np.ndarray
) -> sp.csr_matrix:
    """
    Matrice de MASSE P1 condensee (lumped) : chaque tetraedre distribue son
    volume a parts egales sur ses 4 sommets, et la matrice est diagonale.

    Condensee et non consistante : la methode de la chaleur ne demande qu'un
    terme de masse bien conditionne, et une diagonale garantit que
    `M + t*L` reste definie positive sans hypothese sur la qualite du
    maillage - point qui compte ici, `tetgen` ne garantissant aucune
    regularite (meme raison que le `abs(volume)` deja utilise ailleurs dans
    ce fichier).
    """
    node_count = tet_mesh.node_count
    contributions = np.repeat(np.abs(volumes) / 4.0, 4)
    indices = tet_mesh.tetrahedra.ravel()
    diagonal = np.zeros(node_count)
    np.add.at(diagonal, indices, contributions)
    return sp.diags(diagonal, format="csr")


def _assemble_divergence_rhs(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    tet_field: np.ndarray,
) -> np.ndarray:
    """
    Second membre de l'equation de Poisson `L*phi = div(X)` pour un champ
    vectoriel `X` constant par tetraedre.

    En P1, `integral(grad(w_i) . X)` sur un element vaut `V_e * (G_e^T X_e)_i` :
    c'est exactement la meme contraction que
    `_assemble_gradient_direction_term` utilise pour son second membre, avec
    le volume comme poids au lieu du poids de surplomb. Aucune nouvelle
    formule d'integration n'est introduite.
    """
    local = volumes[:, None] * np.einsum("tki,tk->ti", shape_gradients, tet_field)
    rhs = np.zeros(tet_mesh.node_count)
    np.add.at(rhs, tet_mesh.tetrahedra.ravel(), local.ravel())
    return rhs


def _heat_method_distance(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    source_nodes: np.ndarray,
    params: FieldSolverParams,
) -> np.ndarray:
    """
    Distance geodesique LISSE au plateau, par la methode de la chaleur
    (Crane, Weischedel, Wardetzky).

    Pourquoi elle plutot que Dijkstra sur les aretes : Dijkstra suit les
    aretes du maillage en zigzag, ce qui donne une distance dont la norme du
    gradient vaut 4,81 en mediane sur death.stl au lieu de 1,0 attendu pour
    une distance. Cette irregularite la rend inutilisable comme champ
    generateur de couches. La methode de la chaleur produit une distance
    CONTINUE dont `||grad|| ~ 1` par construction, sur le meme domaine et
    avec les memes matrices FEM deja assemblees ici.

    Trois etapes :
      1. diffusion `(M + t*L) u = source` sur un temps court `t` ;
      2. direction du flux `X = -grad(u)/||grad(u)||` par tetraedre - seule
         la DIRECTION est retenue, c'est ce qui elimine l'irregularite
         d'amplitude ;
      3. integration `L*phi = div(X)`, puis recalage a 0 sur la source.

    `t` est derive de la taille caracteristique du maillage (carre de la
    longueur d'arete mediane), jamais une constante absolue : la formulation
    de la chaleur est dimensionnee en longueur au carre, une valeur fixe
    n'aurait pas de sens d'un modele a l'autre.

    Retourne un tableau par NOEUD. Les deux systemes sont singuliers du seul
    champ constant (Laplacien pur), ce que la source pour le premier et le
    recalage final pour le second suffisent a lever.
    """
    node_count = tet_mesh.node_count
    if len(source_nodes) == 0:
        return np.full(node_count, np.inf)

    L = _assemble_stiffness_matrix(tet_mesh, shape_gradients, volumes)
    M = _assemble_lumped_mass_matrix(tet_mesh, volumes)

    # Temps de diffusion, homogene a une longueur au carre.
    edges = tet_mesh.tet_vertices[tet_mesh.tetrahedra[:, [0, 1]]]
    h = float(np.median(np.linalg.norm(edges[:, 0] - edges[:, 1], axis=1))) or 1.0
    t = params.heat_method_time_factor * h * h

    source = np.zeros(node_count)
    source[source_nodes] = 1.0
    u = _solve_linear_system((M + t * L).tocsr(), source, params)

    gradients = _compute_per_tet_gradient(tet_mesh, shape_gradients, u)
    norms = np.linalg.norm(gradients, axis=1)
    flow = np.zeros_like(gradients)
    usable = norms > 1e-12
    flow[usable] = -gradients[usable] / norms[usable][:, None]

    rhs = _assemble_divergence_rhs(tet_mesh, shape_gradients, volumes, flow)

    # Le Laplacien pur est singulier (noyau = constantes). On ancre les
    # noeuds source a 0 par penalite forte plutot que d'eliminer des lignes :
    # meme mecanisme que l'ancrage plateau deja en place, et le systeme reste
    # symetrique defini positif.
    anchor = sp.diags(
        np.where(np.isin(np.arange(node_count), source_nodes),
                 params.platform_anchor_weight, 0.0),
        format="csr",
    )
    phi = _solve_linear_system((L + anchor).tocsr(), rhs, params)

    # Une distance est positive et nulle a la source ; le signe de la
    # solution de Poisson n'est pas garanti par le systeme lui-meme.
    if phi[source_nodes].mean() > phi.mean():
        phi = -phi
    return phi - phi[source_nodes].min()


def _geodesic_distance_to_platform(
    tet_mesh: TetMesh, platform_nodes: np.ndarray
) -> np.ndarray:
    """
    Distance GEODESIQUE de chaque noeud au plateau, mesuree A TRAVERS LA
    MATIERE : plus court chemin sur le graphe des aretes du TetMesh, pondere
    par les longueurs euclidiennes reelles, depuis l'ensemble des noeuds
    plateau (Dijkstra multi-source, `min_only=True`).

    Suivre les aretes du maillage garantit que le chemin reste dans le
    solide - c'est precisement ce qui distingue cette distance de la
    hauteur : une piece suspendue est proche du plateau A VOL D'OISEAU mais
    geodesiquement lointaine, puisque sa matiere ne le rejoint qu'en faisant
    le tour par ce qui la porte.

    Un noeud inatteignable (composante de matiere deconnectee du plateau)
    recoit `inf` ; l'appelant decide quoi en faire plutot que de se voir
    imposer une valeur invente ici.
    """
    vertices = tet_mesh.tet_vertices
    tets = tet_mesh.tetrahedra

    # Les 6 aretes de chaque tetraedre.
    pairs = np.vstack([
        tets[:, [0, 1]], tets[:, [0, 2]], tets[:, [0, 3]],
        tets[:, [1, 2]], tets[:, [1, 3]], tets[:, [2, 3]],
    ])
    lengths = np.linalg.norm(vertices[pairs[:, 0]] - vertices[pairs[:, 1]], axis=1)

    node_count = tet_mesh.node_count
    graph = sp.coo_matrix(
        (
            np.concatenate([lengths, lengths]),
            (
                np.concatenate([pairs[:, 0], pairs[:, 1]]),
                np.concatenate([pairs[:, 1], pairs[:, 0]]),
            ),
        ),
        shape=(node_count, node_count),
    ).tocsr()

    if len(platform_nodes) == 0:
        return np.full(node_count, np.inf)

    return spgraph.dijkstra(
        graph, directed=False, indices=platform_nodes, min_only=True
    )


def _geodesic_tet_directions(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    platform_nodes: np.ndarray,
    fallback_directions: np.ndarray,
) -> tuple[np.ndarray, int]:
    """
    Direction cible par tetraedre = gradient NORMALISE de la distance
    geodesique au plateau : il pointe, en chaque point, dans le sens ou la
    matiere s'eloigne du plateau a travers le solide.

    Seule la DIRECTION est retenue, jamais la magnitude : la distance issue
    de Dijkstra est trop irreguliere pour ca (voir
    `FieldSolverParams.use_geodesic_directions`).

    Repli sur `fallback_directions` (les directions par region) la ou le
    gradient geodesique est inexploitable - tetraedre dont un noeud est
    inatteignable, ou gradient de norme quasi nulle (plateau de distance,
    typiquement a egale distance de deux chemins concurrents). Retourne
    aussi le nombre de tetraedres ainsi repliees, pour que l'appelant puisse
    le rapporter plutot que de le masquer.
    """
    geodesic = _geodesic_distance_to_platform(tet_mesh, platform_nodes)

    finite = np.isfinite(geodesic)
    if not finite.any():
        return fallback_directions, tet_mesh.tet_count

    # Les noeuds inatteignables sont mis au maximum fini AVANT le calcul du
    # gradient : un `inf` contaminerait tout tetraedre les touchant.
    filled = np.where(finite, geodesic, geodesic[finite].max())
    gradients = _compute_per_tet_gradient(tet_mesh, shape_gradients, filled)

    norms = np.linalg.norm(gradients, axis=1)
    usable = norms > 1e-9
    # Un tetraedre touchant un noeud inatteignable a un gradient sans
    # signification physique : replier plutot que de l'utiliser.
    usable &= finite[tet_mesh.tetrahedra].all(axis=1)

    directions = fallback_directions.copy()
    directions[usable] = gradients[usable] / norms[usable][:, None]
    return directions, int((~usable).sum())


def _geodesic_saddle_mask(
    tet_mesh: TetMesh,
    directions: np.ndarray,
    threshold_deg: float,
) -> np.ndarray:
    """
    Detecte les tetraedres situes a un POINT DE SELLE du gradient geodesique -
    la ou deux chemins vers le plateau quasi equidistants se disputent la
    direction locale.

    Mesure a l'origine de cette fonction (real_overhang.stl, resolution 0.4) :
    ce n'est PAS un bruit haute-frequence diffus, ni un artefact de mauvaise
    qualite de maillage (correlation qualite/discordance faible, r=+0.36 ;
    ||grad(distance)|| reste 0.999 des DEUX cotes d'une paire discordante,
    identique au reste du modele). C'est un phenomene LOCALISE : 2.6% des
    paires de tetraedres voisins depassent 100 deg d'ecart, avec des pics
    a 177 deg, en plein corps d'un pilier par ailleurs parfaitement droit.
    Un lissage spatial par sauts topologiques, teste jusqu'a 8 sauts,
    n'attaque pas cette echelle : le MAXIMUM reste ~170 deg quel que soit le
    rayon de lissage, parce que moyenner a travers un point de selle donne un
    vecteur proche de zero puis renormalise n'importe comment, pas une
    direction corrigee.

    Un tetraedre est marque des qu'un SEUL de ses voisins directs s'ecarte de
    plus de `threshold_deg` - la discordance est un phenomene de PAIRE, donc
    les deux membres d'une paire fautive sont marques, jamais un seul.

    Retourne un masque booleen (T,) : True = tetraedre a exclure du terme
    source geodesique (voir `geodesic_source_weight`), le repli se faisant
    sur le terme de direction principal (tet_directions) qui, lui, reste
    lisse par construction (melange spatial de `_blend_tet_directions`).
    """
    adjacency = _build_tet_tet_adjacency(tet_mesh).tocoo()
    dots = np.einsum("ij,ij->i", directions[adjacency.row], directions[adjacency.col])
    angles_deg = np.degrees(np.arccos(np.clip(dots, -1.0, 1.0)))

    discordant_pairs = angles_deg > threshold_deg
    saddle = np.zeros(tet_mesh.tet_count, dtype=bool)
    saddle[adjacency.row[discordant_pairs]] = True
    saddle[adjacency.col[discordant_pairs]] = True
    return saddle


def _iterative_vector_smoothing(
    tet_mesh: TetMesh,
    directions: np.ndarray,
    iterations: int,
    step: float,
) -> np.ndarray:
    """
    Lissage NON-LINEAIRE ITERATIF d'un champ de vecteurs unitaires par
    tetraedre - voir `FieldSolverParams.geodesic_smoothing_iterations` pour
    pourquoi ceci reussit la ou `_diffuse_tet_directions` (solve lineaire en
    une seule passe, composante par composante) echoue aux points de selle.

    A CHAQUE iteration : chaque tetraedre se deplace d'une fraction `step`
    vers la moyenne (non normalisee) de ses voisins directs, PUIS le
    resultat est renormalise - jamais un seul calcul suivi d'une
    renormalisation finale. C'est cette renormalisation REPETEE qui permet
    au champ de "basculer" progressivement hors d'un point de selle plutot
    que de rester fige au vecteur quasi nul de la moyenne lineaire.

    Implementation vectorisee (`np.add.at`, pas de boucle Python sur les
    tetraedres) : le cout par iteration est O(aretes du graphe tet-tet),
    donc lineaire dans la taille du maillage, independant du nombre
    d'iterations demande cote complexite asymptotique.
    """
    adjacency = _build_tet_tet_adjacency(tet_mesh).tocoo()
    row, col = adjacency.row, adjacency.col
    tet_count = tet_mesh.tet_count

    current = directions.copy()
    for _ in range(iterations):
        neighbor_sum = np.zeros_like(current)
        neighbor_count = np.zeros(tet_count)
        np.add.at(neighbor_sum, row, current[col])
        np.add.at(neighbor_count, row, 1)
        neighbor_count[neighbor_count == 0] = 1
        neighbor_avg = neighbor_sum / neighbor_count[:, None]

        blended = (1.0 - step) * current + step * neighbor_avg
        norms = np.linalg.norm(blended, axis=1)
        usable = norms > 1e-9
        current[usable] = blended[usable] / norms[usable][:, None]
    return current


def _branch_indicator(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    platform_nodes: np.ndarray,
    params: FieldSolverParams,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Indicateur de RAMIFICATION par tetraedre : |div(X)|, ou X est le flux
    geodesique normalise (gradient de la distance lisse au plateau).

    Pourquoi la divergence et pas la courbure du flux : dans un tube de
    section constante, le front de distance avance sans changer d'aire, donc
    div(X) ~ 0 ; a un branchement, l'aire du front change d'un coup (elle
    augmente quand une branche part, diminue quand le front entre dans un
    tube plus etroit) et la divergence s'eleve. Le SIGNE distingue les deux
    cas, la VALEUR ABSOLUE les traite comme ce qu'ils sont tous les deux :
    un endroit ou la direction de depot doit tourner vite.

    PISTE ECARTEE PAR LA MESURE - la divergence ANGULAIRE du flux entre
    tetraedres voisins, `max(1 - |g_i . g_j|)`, qui detecte le lieu de
    coupure de la distance. Elle ne s'allume PAS a un branchement, et la
    mesure a explique pourquoi : sur une piece ARBORESCENTE les geodesiques
    montent le tronc puis entrent dans la branche sans jamais entrer en
    concurrence, donc la distance y est parfaitement lisse. Un lieu de
    coupure n'existe qu'aux EXTREMITES, la ou deux fronts se rencontrent -
    ce que cet indicateur trouvait effectivement, et exclusivement. Ce qui
    change a un branchement est la TOPOLOGIE de l'ensemble de niveau, pas la
    regularite de la distance.

    Mesure sur `real_overhang.stl` (rapport a la mediane du corps du pilier) :
        corps du pilier, 20mm de long ......   1.0x  (plat)
        milieu du bras .....................   2.0x
        approche de la jonction ............   8x -> 29x -> 56x
        jonction pilier/bras ...............  103x
        jonction bras/pilier suspendu ......   66x

    FAUX POSITIFS CONNUS, mesures et rendus explicites plutot que masques :
    l'indicateur s'eleve aussi aux FRONTIERES du champ de distance, qui ne
    sont pas des branchements - la source plateau (26x) et les extremites
    libres (257x au sommet du pilier, 805x au bout du pilier suspendu). Ils
    sont comptes dans `boundary_tet_count` et bornes par le plancher de
    poids, jamais supprimes : les retirer demanderait un seuil arbitraire, et
    un tetraedre d'extremite est de toute facon un endroit ou le champ est
    peu contraint.

    Retourne `(indicateur (T,), masque des tetraedres de frontiere (T,))`.
    """
    distance = _heat_method_distance(
        tet_mesh, shape_gradients, volumes, platform_nodes, params
    )

    gradients = _compute_per_tet_gradient(tet_mesh, shape_gradients, distance)
    norms = np.linalg.norm(gradients, axis=1)
    usable = norms > 1e-12
    flow = np.zeros_like(gradients)
    flow[usable] = gradients[usable] / norms[usable][:, None]

    # div(X) par noeud = (integrale de grad(phi_i).X) / masse condensee. Meme
    # assemblage que l'etape 3 de la methode de la chaleur, reutilise tel quel.
    rhs = _assemble_divergence_rhs(tet_mesh, shape_gradients, volumes, flow)
    mass = _assemble_lumped_mass_matrix(tet_mesh, volumes).diagonal()
    mass_ok = mass > 1e-14
    div_node = np.zeros(tet_mesh.node_count)
    div_node[mass_ok] = rhs[mass_ok] / mass[mass_ok]

    indicator = np.abs(div_node[tet_mesh.tetrahedra].mean(axis=1))
    # Un tetraedre dont le flux n'est pas exploitable ne dit rien : lui
    # attribuer 0 le laisse au lissage nominal, plutot que de lui inventer une
    # ramification.
    indicator[~usable] = 0.0

    boundary_nodes = np.zeros(tet_mesh.node_count, dtype=bool)
    boundary_nodes[platform_nodes] = True
    for key, owners in _tet_face_table(tet_mesh).items():
        if len(owners) == 1:                     # face sans voisin = frontiere
            boundary_nodes[list(key)] = True
    boundary_tets = boundary_nodes[tet_mesh.tetrahedra].any(axis=1)

    return indicator, boundary_tets


def _smoothness_modulation_weights(
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    platform_nodes: np.ndarray,
    params: FieldSolverParams,
) -> tuple[np.ndarray, SmoothnessModulationDiagnostic]:
    """
    Poids de lissage par tetraedre, derives de `_branch_indicator` :

        poids_t = max(1 / (1 + strength * indicateur_t), floor)

    Forme choisie pour trois proprietes simultanees :
      - `strength = 0` donne exactement 1.0 partout (l'appelant n'appelle
        alors meme pas cette fonction, mais la formule reste coherente) ;
      - elle est CONTINUE, jamais un seuil binaire - piege deja paye avec la
        penalite d'ancrage de step3, ou un `max(0, ...)` rendait toutes les
        directions non-opposees strictement equivalentes et laissait l'argmin
        trancher au hasard. Une consequence utile : sur une piece SANS
        ramification l'indicateur est plat, les poids sont uniformes, et le
        mecanisme est inerte par construction plutot que par cas particulier ;
      - le plancher borne le contraste de raideur a `1/floor`, donc le
        conditionnement du systeme, donc le surcout du solveur iteratif.
    """
    indicator, boundary_tets = _branch_indicator(
        tet_mesh, shape_gradients, volumes, platform_nodes, params
    )

    weights = 1.0 / (1.0 + params.smoothness_modulation_strength * indicator)
    floored = weights < params.smoothness_modulation_floor
    weights = np.maximum(weights, params.smoothness_modulation_floor)

    levels = [10.0, 50.0, 90.0, 99.0, 100.0]
    diagnostic = SmoothnessModulationDiagnostic(
        indicator_quantile_levels=levels,
        indicator_quantiles=[float(v) for v in np.percentile(indicator, levels)],
        min_weight=float(weights.min()),
        median_weight=float(np.median(weights)),
        floored_tet_count=int(floored.sum()),
        boundary_tet_count=int(boundary_tets.sum()),
    )
    return weights, diagnostic


def _diffuse_tet_directions(
    tet_mesh: TetMesh,
    tet_directions: np.ndarray,
    transition_distance_mm: float,
) -> tuple[np.ndarray, int]:
    """
    Diffuse le champ de DIRECTIONS CIBLES sur le graphe tetraedre-tetraedre,
    pour que la transition entre deux branches d'orientation differente
    s'etale sur `transition_distance_mm` au lieu d'etre un saut net a la
    frontiere (voir `FieldSolverParams.transition_distance_mm`).

    AUTO-LIMITE aux frontieres par construction : a l'interieur d'une
    branche, les directions cibles sont identiques d'un tetraedre a l'autre,
    donc les moyenner ne change rien. Seules les zones ou les directions
    different reellement - c'est-a-dire les jonctions - sont affectees.
    Aucun masque geometrique n'est donc necessaire, et aucune identification
    explicite de "la branche d'en face" non plus.

    Resolution IMPLICITE `(I + t*L) x = x0` par composante, avec `L` le
    Laplacien du graphe pondere par les distances reelles entre centres de
    tetraedres, plutot qu'un moyennage iteratif explicite : le nombre
    d'iterations d'un schema explicite croitrait en `(distance/arete)^2`
    (~2800 iterations pour 40 mm sur un maillage a 0.75 mm d'arete),
    alors que le schema implicite coute trois resolutions creuses quelle que
    soit la distance demandee.

    `t = transition^2 / 2` : pour une diffusion, la longueur caracteristique
    est `sqrt(2*t)`, donc ce choix fait de `transition_distance_mm` la
    largeur reelle de la zone de transition, et non un reglage abstrait.

    Retourne `(directions diffusees, nombre de tetraedres reellement
    modifies au-dela du bruit numerique)` - ce compte est la donnee qui
    permet de verifier que le mecanisme a agi, plutot que de le supposer.
    """
    adjacency = _build_tet_tet_adjacency(tet_mesh).tocoo()
    tet_count = tet_mesh.tet_count
    if adjacency.nnz == 0:
        return tet_directions, 0

    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    lengths = np.linalg.norm(
        tet_centers[adjacency.row] - tet_centers[adjacency.col], axis=1
    )
    # Poids de conduction 1/distance : deux tetraedres proches echangent plus
    # que deux tetraedres eloignes, ce qui rend la diffusion coherente sur un
    # maillage non uniforme (le raffinement local pres du plateau en produit
    # un, voir step1).
    conduction = 1.0 / np.maximum(lengths, 1e-9)

    weights = sp.coo_matrix(
        (conduction, (adjacency.row, adjacency.col)), shape=(tet_count, tet_count)
    ).tocsr()
    degree = np.asarray(weights.sum(axis=1)).ravel()
    laplacian = sp.diags(degree) - weights

    t = 0.5 * float(transition_distance_mm) ** 2
    system = (sp.identity(tet_count, format="csr") + t * laplacian).tocsc()

    diffused = np.empty_like(tet_directions)
    for axis in range(3):
        diffused[:, axis] = spla.spsolve(system, tet_directions[:, axis])

    norms = np.linalg.norm(diffused, axis=1)
    usable = norms > 1e-9
    result = tet_directions.copy()
    result[usable] = diffused[usable] / norms[usable][:, None]

    changed = int(np.count_nonzero(
        np.linalg.norm(result - tet_directions, axis=1) > 1e-6
    ))
    return result, changed


# Plafond DUR du nombre d'iterations de la boucle de retroaction, quel que
# soit `FieldSolverParams.support_feedback_iterations`. Chaque iteration coute
# une resolution complete du systeme (~15s sur death.stl) : sans borne, un
# parametrage distrait rendrait le pipeline inutilisable sans jamais lever
# d'erreur. Le critere d'arret dynamique termine de toute facon bien avant sur
# les modeles mesures.
SUPPORT_FEEDBACK_ITERATION_CAP = 50


def _self_overhang_angles(
    tet_gradients: np.ndarray, tet_normals: np.ndarray
) -> np.ndarray:
    """
    Angle (radians) entre la normale de surface heritee par chaque tetraedre
    et la verticale descendante LOCALE du champ, `-grad(phi)`.

    C'est exactement la quantite calculee par
    `step2._compute_overhang_angles`, a une seule difference - et c'est tout
    l'objet du mecanisme : step2 la mesure contre `-up_vector` (la verticale
    MACHINE, figee avant que le champ n'existe), alors qu'ici la reference est
    la direction de depot que le champ produit REELLEMENT a cet endroit. Un
    angle petit signifie que la surface regarde dans le sens oppose au depot :
    la matiere y serait posee dans le vide.

    Un tetraedre dont le gradient est quasi nul (champ localement plat) n'a
    pas de direction de depot definie : il recoit `pi` (jamais un surplomb),
    plutot qu'un angle calcule sur un vecteur de bruit numerique.
    """
    norms = np.linalg.norm(tet_gradients, axis=1)
    usable = norms > 1e-12
    down = np.zeros_like(tet_gradients)
    down[usable] = -tet_gradients[usable] / norms[usable][:, None]

    dots = np.clip(np.einsum("ij,ij->i", tet_normals, down), -1.0, 1.0)
    angles = np.arccos(dots)
    angles[~usable] = np.pi
    return angles


def _rotate_directions_out_of_overhang(
    directions: np.ndarray,
    normals: np.ndarray,
    angles: np.ndarray,
    mask: np.ndarray,
    threshold_rad: float,
) -> np.ndarray:
    """
    Fait pivoter les directions cibles marquees par `mask` du MINIMUM
    necessaire pour sortir du seuil de surplomb, et rien de plus.

    Geometrie : en posant `u = -d`, le surplomb s'ecrit `angle(n, u) <
    seuil`. Tourner `u` autour de `n x u` d'un angle `beta` augmente
    exactement `angle(n, u)` de `beta` (rotation dans le plan (n, u), sens
    verifie analytiquement). On applique donc `beta = seuil - angle`, la
    correction minimale : une rotation plus ample eloignerait la cible de ce
    que la geometrie de la region demande, sans rien resoudre de plus.

    Cas degenere `n x u ~ 0` - la normale est colineaire a la verticale
    locale, c'est-a-dire un plafond parfaitement plat, precisement le pire
    cas : aucune direction de rotation n'est privilegiee par la geometrie. Un
    axe perpendiculaire stable est alors construit (meme repli que
    `_perpendicular_basis` : reference [1,0,0], bascule sur [0,1,0] si
    quasi-colineaire) plutot que de laisser passer une rotation indefinie.
    """
    if not mask.any():
        return directions

    corrected = directions.copy()
    idx = np.where(mask)[0]

    d = directions[idx]
    norms = np.linalg.norm(d, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    u = -d / norms
    n = normals[idx]

    axis = np.cross(n, u)
    axis_norm = np.linalg.norm(axis, axis=1)

    degenerate = axis_norm < 1e-9
    if degenerate.any():
        reference = np.tile(np.array([1.0, 0.0, 0.0]), (degenerate.sum(), 1))
        colinear = np.abs(u[degenerate] @ np.array([1.0, 0.0, 0.0])) > 0.9
        reference[colinear] = np.array([0.0, 1.0, 0.0])
        fallback = np.cross(u[degenerate], reference)
        axis[degenerate] = fallback
        axis_norm[degenerate] = np.linalg.norm(fallback, axis=1)

    axis_norm[axis_norm < 1e-12] = 1.0
    axis = axis / axis_norm[:, None]

    beta = threshold_rad - angles[idx]

    # Rodrigues, applique par paire (chaque vecteur a son propre axe et son
    # propre angle) - pas `_rodrigues_rotate`, qui fait tourner un lot de
    # vecteurs autour d'UN axe commun.
    cos_b = np.cos(beta)[:, None]
    sin_b = np.sin(beta)[:, None]
    dot = np.einsum("ij,ij->i", axis, u)[:, None]
    rotated = (
        u * cos_b
        + np.cross(axis, u) * sin_b
        + axis * dot * (1.0 - cos_b)
    )

    corrected[idx] = -rotated
    return corrected


def _apply_geodesic_targets(
    directions: np.ndarray,
    geodesic_directions: np.ndarray,
    normals: np.ndarray,
    mask: np.ndarray,
    threshold_rad: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Pour les tetraedres de `mask`, substitue la direction GEODESIQUE LOCALE -
    mais UNIQUEMENT la ou elle resout effectivement le surplomb a ce
    tetraedre-la.

    Decision prise tetraedre par tetraedre, jamais par zone : c'est une
    direction moyennee sur toute une region qui cree le defaut qu'on corrige
    ici (voir `support_feedback_use_geodesic_target`). Substituer par zone
    reproduirait le probleme a une autre echelle.

    Le test d'acceptation est le MEME critere que partout ailleurs
    (`angle(normale, -direction) >= seuil`), evalue sur la direction candidate.
    Une geodesique qui serait elle-meme en surplomb ici est donc refusee, et
    l'appelant retombe sur la rotation minimale - la provenance geodesique
    n'accorde aucun passe-droit.

    Retourne `(directions, masque des tetraedres NON traites)`, ce second
    masque etant exactement ceux qu'il reste a corriger par rotation.
    """
    result = directions.copy()
    if not mask.any():
        return result, mask

    idx = np.where(mask)[0]
    candidate = geodesic_directions[idx]
    norms = np.linalg.norm(candidate, axis=1)
    usable = norms > 1e-12

    unit = np.zeros_like(candidate)
    unit[usable] = candidate[usable] / norms[usable][:, None]
    dots = np.clip(np.einsum("ij,ij->i", normals[idx], -unit), -1.0, 1.0)
    accepted = usable & (np.arccos(dots) >= threshold_rad)

    result[idx[accepted]] = unit[accepted]

    remaining = mask.copy()
    remaining[idx[accepted]] = False
    return result, remaining


def _apply_support_feedback(
    A_without_gradient: sp.csr_matrix,
    b_without_gradient: np.ndarray,
    tet_mesh: TetMesh,
    shape_gradients: np.ndarray,
    volumes: np.ndarray,
    tet_directions: np.ndarray,
    tet_weights: np.ndarray,
    tet_normals: np.ndarray,
    params: FieldSolverParams,
    geodesic_directions: Optional[np.ndarray] = None,
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray, SupportFeedbackDiagnostic]:
    """
    Boucle de point fixe : resoudre, re-evaluer le surplomb contre le champ
    REELLEMENT obtenu, corriger les cibles fautives, recommencer.

    Retourne le systeme (A, b) assemble avec les cibles corrigees, les
    directions finales, et le diagnostic. Le systeme rendu est toujours
    coherent avec la derniere correction appliquee, y compris a la
    convergence.

    Ne resout jamais les Constraint[] utilisateur : celles-ci sont ajoutees
    plus tard par `solve_with_constraints`, et les inclure ici ferait diverger
    la boucle d'un probleme que l'appelant ne resoudra pas ainsi. La boucle
    travaille donc sur le seul systeme structurel - meme raison qui fait
    assembler l'ancrage plateau AVANT le terme de direction.
    """
    threshold_rad = np.radians(params.support_feedback_threshold_deg)
    total_volume = volumes.sum()

    diagnostic = SupportFeedbackDiagnostic()
    directions = tet_directions
    weights = tet_weights.copy()
    corrected_ever = np.zeros(tet_mesh.tet_count, dtype=bool)

    max_iterations = min(
        params.support_feedback_iterations, SUPPORT_FEEDBACK_ITERATION_CAP
    )

    A = b = None
    previous_fraction = None
    for iteration in range(max_iterations):
        gradient_A, gradient_b = _assemble_gradient_direction_term(
            tet_mesh,
            shape_gradients,
            volumes,
            np.ones(tet_mesh.tet_count, dtype=bool),
            directions,
            weights,
        )
        A = (A_without_gradient + gradient_A).tocsr()
        b = b_without_gradient + gradient_b

        values = _solve_linear_system(A, b, params)
        field_gradients = _compute_per_tet_gradient(tet_mesh, shape_gradients, values)
        angles = _self_overhang_angles(field_gradients, tet_normals)
        offending = angles < threshold_rad

        count = int(offending.sum())
        fraction = (
            float(volumes[offending].sum() / total_volume) if total_volume > 0.0 else 0.0
        )
        if iteration == 0:
            diagnostic.initial_self_overhang_tet_count = count
            diagnostic.initial_self_overhang_volume_fraction = fraction
        diagnostic.residual_self_overhang_tet_count = count
        diagnostic.residual_self_overhang_volume_fraction = fraction
        diagnostic.iterations_used = iteration + 1
        diagnostic.volume_fraction_history.append(fraction)

        if count == 0:
            diagnostic.converged = True
            break

        # Arret dynamique : le gain d'une iteration a la suivante est devenu
        # negligeable. Le test porte sur le gain RELATIF, et un residu qui
        # REMONTE (gain negatif) arrete aussi la boucle - continuer a corriger
        # un champ qui se degrade n'a aucune raison de le redresser.
        if previous_fraction is not None and previous_fraction > 0.0:
            improvement = (previous_fraction - fraction) / previous_fraction
            if improvement < params.support_feedback_min_improvement:
                diagnostic.stalled = True
                break
        previous_fraction = fraction

        # Derniere iteration autorisee : corriger maintenant serait sans
        # effet, la correction ne pouvant plus etre ni assemblee ni resolue.
        # S'en abstenir garde (A, b), `directions` et le diagnostic tous
        # coherents avec le MEME etat.
        if iteration == max_iterations - 1:
            break

        # Cible corrigee, tetraedre par tetraedre : la geodesique LOCALE quand
        # elle resout le surplomb ici, la rotation minimale pour tout le reste.
        to_rotate = offending
        if geodesic_directions is not None:
            directions, to_rotate = _apply_geodesic_targets(
                directions, geodesic_directions, tet_normals, offending, threshold_rad
            )
            diagnostic.geodesic_target_count += int(
                offending.sum() - to_rotate.sum()
            )

        directions = _rotate_directions_out_of_overhang(
            directions, tet_normals, angles, to_rotate, threshold_rad
        )
        # La cible corrigee doit peser plus lourd que le lissage, sinon la
        # correction est absorbee sans effet mesurable. Le poids est pose
        # (pas cumule) : un tetraedre corrige a plusieurs iterations ne doit
        # pas voir son poids exploser au fil de la boucle.
        weights[offending] = np.maximum(
            weights[offending],
            params.gradient_baseline_weight * params.support_feedback_weight,
        )
        corrected_ever |= offending

    diagnostic.corrected_tet_count = int(corrected_ever.sum())

    if A is None:  # support_feedback_iterations <= 0, jamais atteint par run()
        gradient_A, gradient_b = _assemble_gradient_direction_term(
            tet_mesh, shape_gradients, volumes,
            np.ones(tet_mesh.tet_count, dtype=bool), directions, weights,
        )
        A = (A_without_gradient + gradient_A).tocsr()
        b = b_without_gradient + gradient_b

    return A, b, directions, diagnostic


def _build_gradient_up_operator(
    tet_mesh: TetMesh, shape_gradients: np.ndarray, up_vector: np.ndarray
) -> sp.csr_matrix:
    """
    Operateur creux G_up de forme (T, P) tel que
    `(G_up @ phi)[e] == grad(phi)|_e . up_vector` pour chaque tetraedre e.

    C'est exactement la meme algebre que `_compute_per_tet_gradient` (le champ
    P1 a un gradient constant par element, `grad|_e = G_e @ phi_e`), projetee
    sur `up_vector` et exprimee comme une matrice creuse globale plutot que
    comme un einsum dense par element. L'interet de la forme matricielle :
    elle ne depend QUE de la geometrie, donc elle est construite une seule
    fois et chaque passe de la boucle ne paie plus qu'un produit
    matrice-vecteur creux au lieu de re-derouler l'einsum.

    Ligne `e` : 4 entrees non nulles, aux colonnes `tetrahedra[e]`, de valeurs
    `up_vector @ shape_gradients[e]` (shape (4,)).

    Les tetraedres degeneres ont un `shape_gradients` nul par construction
    (voir `_compute_shape_gradients`), donc une ligne nulle ici - aucun
    filtrage supplementaire n'est necessaire.
    """
    tet_count = tet_mesh.tet_count
    node_count = tet_mesh.node_count
    if tet_count == 0:
        return sp.csr_matrix((0, node_count))

    up = np.asarray(up_vector, dtype=np.float64)
    norm = float(np.linalg.norm(up))
    up = up / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0])

    # (T,4) : contribution de chacun des 4 noeuds du tetraedre a grad.up
    data = np.einsum("i,tik->tk", up, shape_gradients)
    rows = np.repeat(np.arange(tet_count, dtype=np.int64), 4)
    cols = tet_mesh.tetrahedra.ravel()

    return sp.coo_matrix(
        (data.ravel(), (rows, cols)), shape=(tet_count, node_count)
    ).tocsr()


def _descending_volume_fraction(
    grad_dot_up: np.ndarray, volumes: np.ndarray, epsilon: float
) -> tuple[float, int]:
    """
    Fraction du VOLUME total dont le gradient descend (`grad.up < -epsilon`),
    et compte brut de tetraedres concernes.

    Pondere par volume et non par compte de tetraedres : un maillage
    tetraedrique n'est pas uniforme (le raffinement local pres du plateau en
    est la preuve directe, cf. `step1_ingest.near_platform_resolution`), donc
    un compte brut surponderait mecaniquement les zones finement maillees.
    Meme discipline que la mediane ponderee par volume de
    `step5._representative_gradient_norm`.
    """
    descending = (grad_dot_up < -epsilon) & (volumes > 0.0)
    total_volume = float(volumes.sum())
    if total_volume <= 0.0:
        return 0.0, int(np.count_nonzero(descending))
    return float(volumes[descending].sum() / total_volume), int(np.count_nonzero(descending))


def _clamp_downward_directions(
    tet_directions: np.ndarray,
    up_vector: np.ndarray,
    max_downward_deg: float,
    epsilon: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    BORNE l'inclinaison descendante des directions de croissance au lieu de
    l'interdire.

    `_enforce_upward_directions` retourne toute direction descendante, sur
    l'hypothese qu'on ne peut pas deposer vers le bas. **Cette hypothese est
    fausse sur cette machine** : la tete couvre B ∈ [0°, 185°], donc elle PEUT
    deposer vers le bas (voir MachineProfile). L'interdit binaire empechait
    donc structurellement d'exprimer le seul ordre correct pour une piece
    suspendue - la faire croitre vers le bas depuis son attache.

    Ce que l'on veut vraiment n'est pas d'interdire la descente, mais d'eviter
    d'y passer trop de temps a l'extreme de la course (B proche de 185°). On
    BORNE donc l'angle sous l'horizontale a `max_downward_deg` :

      - 0°  : rien ne descend sous l'horizontale (le plus proche de l'ancien
              comportement, sans jamais retourner une direction) ;
      - 45° : la croissance peut plonger jusqu'a 45° sous l'horizontale,
              soit B ≈ 135° ;
      - 90° : descente verticale autorisee (B ≈ 180°), aucune borne effective.

    La direction bornee **conserve son azimut** (sa composante horizontale) :
    seule son inclinaison est relevee. La retourner, comme le fait le flip,
    changeait le sens du depot ; ici on ne fait que l'aplatir.

    Cas degenere - direction exactement verticale descendante, sans composante
    horizontale : aucun azimut n'est privilegie par la geometrie. Un azimut
    stable est construit a partir d'une reference fixe plutot que de laisser
    passer une direction indefinie.

    Retourne `(directions bornees, masque des directions modifiees)`. L'entree
    n'est jamais mutee.
    """
    up = np.asarray(up_vector, dtype=np.float64)
    norm = float(np.linalg.norm(up))
    up = up / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0])

    directions = np.array(tet_directions, dtype=np.float64, copy=True)
    lengths = np.linalg.norm(directions, axis=1)
    usable = lengths > 1e-12
    unit = np.zeros_like(directions)
    unit[usable] = directions[usable] / lengths[usable][:, None]

    vertical = unit @ up                                   # = -sin(plongee)
    limit = -np.sin(np.radians(max_downward_deg))
    too_low = usable & (vertical < limit - abs(epsilon))
    if not too_low.any():
        return directions, too_low

    horizontal = unit[too_low] - np.outer(unit[too_low] @ up, up)
    h_norm = np.linalg.norm(horizontal, axis=1)

    # Descente strictement verticale : aucun azimut naturel, on en construit
    # un stable (meme repli que `_perpendicular_basis`).
    degenerate = h_norm < 1e-9
    if degenerate.any():
        reference = np.tile(np.array([1.0, 0.0, 0.0]), (int(degenerate.sum()), 1))
        colinear = np.abs(reference @ up) > 0.9
        reference[colinear] = np.array([0.0, 1.0, 0.0])
        fallback = reference - np.outer(reference @ up, up)
        horizontal[degenerate] = fallback
        h_norm[degenerate] = np.linalg.norm(fallback, axis=1)

    h_norm[h_norm < 1e-12] = 1.0
    horizontal = horizontal / h_norm[:, None]

    clamped = (
        horizontal * np.cos(np.radians(max_downward_deg))
        + up * limit
    )
    directions[too_low] = clamped * lengths[too_low][:, None]
    return directions, too_low


def _enforce_upward_directions(
    tet_directions: np.ndarray, up_vector: np.ndarray, epsilon: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Retourne les directions cibles dont la composante selon `up_vector` est
    negative, pour qu'aucune direction de croissance ne pointe vers le bas.

    Retourne `(directions redressees, masque des directions retournees)`.
    L'entree n'est jamais mutee (copie) - `_blend_tet_directions` reste une
    fonction pure de son cote.

    Une direction exactement perpendiculaire a `up` (|d.up| <= epsilon) est
    laissee TELLE QUELLE : elle ne descend pas, et la retourner ne changerait
    rien a la geometrie tout en rendant le resultat sensible au bruit
    numerique du melange spatial.
    """
    up = np.asarray(up_vector, dtype=np.float64)
    norm = float(np.linalg.norm(up))
    up = up / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0])

    directions = np.array(tet_directions, dtype=np.float64, copy=True)
    descending = (directions @ up) < -abs(epsilon)
    directions[descending] = -directions[descending]
    return directions, descending


# ---------------------------------------------------------------------------
# Decision 4 - Ancrage automatique du plateau
# ---------------------------------------------------------------------------

def _compute_platform_level(
    mesh_data: MeshData, region_dependency_graph: RegionDependencyGraph
) -> Optional[float]:
    """
    Niveau geometrique (z) du plateau, deduit de
    `RegionDependencyGraph.platform_contact_faces` (step3, seule etape qui
    connait la position reelle du plateau) - sommet le plus haut des faces
    en contact, la bande de tolerance de step3 est donc deja incluse.
    Retourne None si aucune face de contact n'est disponible (segmentation
    sans plateau identifie) - source de verite UNIQUE, reutilisee a la fois
    par `_identify_platform_nodes` (ancrage de phi) et
    `_blend_tet_directions` (source d'influence "plateau" du melange de
    direction) pour ne jamais faire diverger silencieusement deux calculs
    du meme niveau.
    """
    platform_faces = region_dependency_graph.platform_contact_faces
    if platform_faces is None:
        return None
    platform_faces = np.asarray(platform_faces, dtype=bool)
    if not np.any(platform_faces):
        return None
    platform_vertex_z = mesh_data.vertices[mesh_data.faces[platform_faces]][:, :, 2]
    return float(platform_vertex_z.max())


def _identify_platform_nodes(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    region_dependency_graph: RegionDependencyGraph,
) -> np.ndarray:
    """
    Retrouve les noeuds du TetMesh reposant sur le plateau.

    `RegionDependencyGraph.platform_contact_faces` (calcule par step3, seule
    etape qui connait la position du plateau) designe QUELLES faces touchent
    le plateau ; le niveau geometrique du plateau en est deduit directement
    (z maximal de leurs sommets), puis les noeuds du TetMesh sont testes
    contre ce plan. Aucun parametre `platform_z` n'est redefini ici : la
    source de verite reste step3.

    Deux pieges evites :
      - NE PAS deduire le plateau de la topologie du graphe ("regions sans
        arete entrante") : des que la segmentation ne produit qu'une seule
        region - cas de toute forme simple, cube compris - cette region est
        a la fois racine et couvre tout le modele, ce qui ancrerait toute la
        surface a 0 et empecherait le champ de croitre.
      - NE PAS passer par la face de surface la plus proche de chaque face
        frontiere du TetMesh : ce mapping (correct pour heriter d'une
        direction de region, ou une petite erreur est lissee) est trop
        grossier pour un predicat binaire. Mesure sur un cube : il ancrait
        des noeuds de paroi laterale jusqu'a z=-1.82 tout en manquant de
        vrais noeuds de base a z=-5.
    """
    platform_faces = region_dependency_graph.platform_contact_faces
    if platform_faces is None:
        raise FieldSolverError(
            "RegionDependencyGraph.platform_contact_faces est absent : "
            "impossible d'ancrer le champ au plateau. Utilisez la sortie de "
            "pipeline/step3_segmentation.run(), ou desactivez explicitement "
            "l'ancrage via FieldSolverParams(enable_platform_anchor=False)."
        )

    platform_level = _compute_platform_level(mesh_data, region_dependency_graph)
    if platform_level is None:
        return np.zeros(0, dtype=np.int64)

    scale = float(np.ptp(tet_mesh.tet_vertices[:, 2])) or 1.0
    epsilon = 1e-9 * scale

    return np.where(tet_mesh.tet_vertices[:, 2] <= platform_level + epsilon)[0]


def _assemble_platform_anchor(
    node_count: int, platform_nodes: np.ndarray, weight: float, anchor_value: float
) -> tuple[sp.csr_matrix, np.ndarray]:
    """
    Penalite molle forcant field(node) ~= anchor_value sur les noeuds du
    plateau. Meme mecanique que les contraintes utilisateur (un seul
    mecanisme d'assemblage dans tout le fichier), avec un poids eleve
    dedie plutot qu'une elimination Dirichlet dure.
    """
    A = sp.csr_matrix((node_count, node_count))
    b = np.zeros(node_count)
    if len(platform_nodes) == 0:
        return A, b

    A = sp.csr_matrix(
        (np.full(len(platform_nodes), weight), (platform_nodes, platform_nodes)),
        shape=(node_count, node_count),
    )
    b[platform_nodes] = weight * anchor_value
    return A, b


# ---------------------------------------------------------------------------
# Decision 5 - Resolution du systeme lineaire
# ---------------------------------------------------------------------------

def _solve_linear_system(
    A: sp.csr_matrix,
    b: np.ndarray,
    params: FieldSolverParams,
    initial_guess: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Resout A.phi = b. Le systeme est symetrique defini positif par
    construction (rigidite FEM + penalites quadratiques), donc CG est le
    solveur iteratif approprie - pas GMRES, reserve aux matrices non
    symetriques.
    """
    node_count = A.shape[0]

    if node_count < params.direct_solve_node_threshold:
        try:
            return spla.spsolve(A.tocsc(), b)
        except Exception as exc:
            raise FieldSolverError(f"Le solveur direct (spsolve) a echoue: {exc}") from exc

    if params.use_ilu_preconditioner:
        try:
            ilu = spla.spilu(A.tocsc())
        except Exception as exc:
            raise FieldSolverError(f"La factorisation ILU a echoue: {exc}") from exc
        preconditioner = spla.LinearOperator(A.shape, ilu.solve)
    else:
        diagonal = A.diagonal().copy()
        diagonal[diagonal == 0.0] = 1.0
        preconditioner = spla.LinearOperator(A.shape, lambda x: x / diagonal)

    values, info = spla.cg(
        A,
        b,
        x0=initial_guess,
        rtol=params.convergence_tolerance,
        maxiter=params.max_solver_iterations,
        M=preconditioner,
    )

    if info > 0:
        raise FieldSolverError(
            f"Le solveur iteratif (cg) n'a pas converge apres "
            f"{params.max_solver_iterations} iterations (info={info}). Augmentez "
            f"max_solver_iterations ou convergence_tolerance, ou activez "
            f"use_ilu_preconditioner."
        )
    if info < 0:
        raise FieldSolverError(f"Le solveur iteratif (cg) a recu une entree illegale (info={info}).")

    return values


# ---------------------------------------------------------------------------
# Decision 6 - Systeme statique (invariant aux contraintes) / resolution
# ---------------------------------------------------------------------------


def _project_to_reachable_radial_tilt(
    tet_mesh: TetMesh, directions: np.ndarray, origin: np.ndarray, up_vector: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Projette chaque direction sur le demi-plan (+radial, vertical) de son
    tetraedre : seul degre d'orientation de la machine actuelle quand C est
    impose par la position (limite #59). Rend (directions, masque des
    tetraedres modifies de plus de 1e-6). Pres de l'axe (rayon < 1e-6), la
    radiale est indefinie et C y est libre : direction laissee intacte."""
    up = np.asarray(up_vector, dtype=float)
    up = up / np.linalg.norm(up)
    centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1) - np.asarray(origin, dtype=float)
    radial = centers - np.outer(centers @ up, up)
    rn = np.linalg.norm(radial, axis=1)
    ok = rn > 1e-6
    rhat = np.zeros_like(radial)
    rhat[ok] = radial[ok] / rn[ok, None]
    d_r = np.einsum("ij,ij->i", directions, rhat)
    d_z = directions @ up
    d_r = np.maximum(d_r, 0.0)
    projected = d_r[:, None] * rhat + d_z[:, None] * up
    norm = np.linalg.norm(projected, axis=1)
    degenerate = norm < 1e-9          # direction purement tangentielle : verticale
    projected[degenerate] = up
    norm[degenerate] = 1.0
    projected /= norm[:, None]
    result = np.where(ok[:, None], projected, directions)
    changed = np.linalg.norm(result - directions, axis=1) > 1e-6
    return result, changed

def build_static_system(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    region_dependency_graph: RegionDependencyGraph,
    params: FieldSolverParams,
) -> FieldSolverStaticSystem:
    """
    Assemble la partie du systeme qui NE depend PAS des Constraint[] :
    lissage (rigidite FEM), terme de surplomb, ancrage plateau, plus les
    structures partagees (KDTree des noeuds, gradients de fonctions de
    forme, volumes).

    A appeler UNE FOIS apres step3. Le resultat reste valide tant que la
    geometrie et les poids/seuils ne changent pas ; seule l'edition des
    Constraint[] utilisateur passe ensuite par `solve_with_constraints()`.
    """
    shape_gradients, volumes = _compute_shape_gradients(tet_mesh)

    # Motif de sparsite construit UNE FOIS et partage par les deux matrices
    # (rigidite et terme de direction) : il ne depend que de la topologie du
    # maillage, jamais des poids. Passe explicitement plutot que memorise -
    # voir `_FemPattern`.
    fem_pattern = _FemPattern(tet_mesh) if _cpp.enabled() else None

    # Modulation spatiale du lissage : le poids scalaire global ci-dessous
    # devient un poids PAR TETRAEDRE, abaisse la ou la piece se ramifie (voir
    # `smoothness_modulation_strength`). Desactive par defaut, auquel cas
    # `_assemble_stiffness_matrix` est appele sans ponderation du tout.
    smoothness_weights = None
    smoothness_modulation_diagnostic = None
    if params.smoothness_modulation_strength > 0.0:
        smoothness_weights, smoothness_modulation_diagnostic = (
            _smoothness_modulation_weights(
                tet_mesh,
                shape_gradients,
                volumes,
                _identify_platform_nodes(tet_mesh, mesh_data, region_dependency_graph),
                params,
            )
        )

    A = params.laplacian_smoothness_weight * _assemble_stiffness_matrix(
        tet_mesh, shape_gradients, volumes, smoothness_weights, pattern=fem_pattern
    )
    b = np.zeros(tet_mesh.node_count)

    boundary_face_nodes, boundary_face_tets, nearest_surface_face_id, face_pairs = (
        _map_boundary_faces_to_surface(tet_mesh, mesh_data)
    )

    # Risque par TETRAEDRE (pas seulement les tetraedres frontiere) :
    # chaque tetraedre herite du risque de la face de surface la plus
    # proche de son centroide, via le kdtree deja construit a l'etape 1 sur
    # les centres de faces. Reste un plus-proche-voisin strict (le melange
    # spatial de _blend_tet_directions ne s'applique qu'a la DIRECTION,
    # pas au risque - hors perimetre de ce fix, cf. CLAUDE.md).
    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    _, nearest_face_per_tet = mesh_data.kdtree.query(tet_centers)
    nearest_face_per_tet = np.asarray(nearest_face_per_tet, dtype=np.int64)

    # Direction d'impression par tetraedre : melange spatial ponderee par
    # distance inverse sur le plateau + chaque region (voir
    # _blend_tet_directions) - remplace le plus-proche-voisin strict qui
    # produisait une transition brusque a la frontiere plateau/region ou
    # entre regions (CLAUDE.md, diagnostic complet). Desactive par defaut
    # (direction_blend_radius=0.0), comportement historique inchange.
    tet_directions = _blend_tet_directions(
        tet_mesh, mesh_data, region_dependency_graph, params.up_vector, params
    )

    # Direction cible issue de la distance geodesique au plateau, quand elle
    # est demandee : elle REMPLACE les directions par region ci-dessus (qui
    # restent le repli la ou la geodesique est inexploitable). Calculee ici,
    # apres le melange spatial, pour que le repli beneficie de ce melange.
    geodesic_fallback_count = None
    geodesic_directions_for_source = None
    if params.use_geodesic_directions or params.geodesic_source_weight > 0.0:
        platform_nodes_for_geodesic = _identify_platform_nodes(
            tet_mesh, mesh_data, region_dependency_graph
        )
        geodesic_directions_for_source, geodesic_fallback_count = _geodesic_tet_directions(
            tet_mesh, shape_gradients, platform_nodes_for_geodesic, tet_directions
        )
        # Lissage VECTORIEL iteratif du champ geodesique, AVANT qu'il ne
        # remplace tet_directions - corrige les points de selle qui rendaient
        # `use_geodesic_directions` seul degradant (limite CLAUDE.md #52),
        # sans passer par le lissage lineaire deja mesure insuffisant (voir
        # `geodesic_smoothing_iterations`). N'affecte QUE le champ utilise
        # comme cible principale ici - `geodesic_targets` (boucle de
        # retroaction, plus bas) reste calcule separement sur le champ brut.
        if params.geodesic_smoothing_iterations > 0:
            geodesic_directions_for_source = _iterative_vector_smoothing(
                tet_mesh,
                geodesic_directions_for_source,
                params.geodesic_smoothing_iterations,
                params.geodesic_smoothing_step,
            )
        if params.use_geodesic_directions:
            tet_directions = geodesic_directions_for_source

    # DIFFUSION des directions cibles : etale la transition entre deux
    # branches d'orientation differente sur `transition_distance_mm` au lieu
    # d'un saut net a la frontiere (voir ce parametre). Appliquee ICI, apres
    # que les directions cibles soient definitives (region ou geodesique) et
    # AVANT le redressement/la borne de plongee, qui doivent voir les
    # directions reellement utilisees. Auto-limitee aux jonctions : a
    # l'interieur d'une branche, moyenner des directions identiques ne change
    # rien.
    transition_diffused_tet_count = None
    if params.transition_distance_mm is not None:
        tet_directions, transition_diffused_tet_count = _diffuse_tet_directions(
            tet_mesh, tet_directions, params.transition_distance_mm
        )

    # Directions geodesiques tenues A PART, comme cible de correction possible
    # de la boucle de retroaction (voir
    # `support_feedback_use_geodesic_target`). Elles ne remplacent RIEN ici :
    # seuls les tetraedres que la boucle mesurera comme fautifs pourront les
    # adopter, un par un. Calculees une seule fois - la geodesique ne depend
    # que de la geometrie, jamais du champ, donc elle est constante d'une
    # iteration a l'autre.
    geodesic_targets = None
    if params.support_feedback_iterations > 0 and params.support_feedback_use_geodesic_target:
        geodesic_targets, _ = _geodesic_tet_directions(
            tet_mesh,
            shape_gradients,
            _identify_platform_nodes(tet_mesh, mesh_data, region_dependency_graph),
            np.zeros((tet_mesh.tet_count, 3), dtype=float),  # repli NUL = refuse
        )

    # Risque effectif : une face posee sur le plateau n'est pas un surplomb,
    # meme si sa normale pointe vers le bas (limite de conception #1 de
    # CLAUDE.md). Meme raisonnement que step3._effective_risk_score().
    effective_risk = np.asarray(overhang_map.risk_score, dtype=np.float64).copy()
    if region_dependency_graph.platform_contact_faces is not None:
        effective_risk[np.asarray(region_dependency_graph.platform_contact_faces, dtype=bool)] = 0.0

    tet_risk = effective_risk[nearest_face_per_tet]

    # Sur-ponderation autour des surplombs : les tetraedres frontiere a
    # risque servent de germes, propages sur overhang_influence_depth sauts.
    risky_tets = np.zeros(tet_mesh.tet_count, dtype=bool)
    if len(boundary_face_tets) > 0:
        risky_boundary = effective_risk[nearest_surface_face_id] > 0.0
        seed_tets = np.unique(boundary_face_tets[risky_boundary])
        if len(seed_tets) > 0:
            risky_tets = _spread_risk_influence(
                _build_tet_tet_adjacency(tet_mesh, face_pairs),
                seed_tets,
                params.overhang_influence_depth,
            )

    # Tous les tetraedres recoivent le poids de base ; ceux marques a risque
    # recoivent en plus la sur-ponderation proportionnelle au risque.
    tet_weights = np.full(tet_mesh.tet_count, params.gradient_baseline_weight, dtype=np.float64)
    tet_weights[risky_tets] += params.gradient_overhang_weight * tet_risk[risky_tets]

    # L'ancrage plateau est assemble AVANT le terme de direction : la
    # penalite iteree anti-descente re-resout le systeme complet a chaque
    # passe, et un systeme sans ancrage serait singulier (un Laplacien seul a
    # un noyau constant) - la boucle resoudrait alors un probleme different
    # de celui que resout finalement solve_with_constraints.
    if params.enable_platform_anchor:
        platform_nodes = _identify_platform_nodes(
            tet_mesh, mesh_data, region_dependency_graph
        )
        anchor_A, anchor_b = _assemble_platform_anchor(
            tet_mesh.node_count,
            platform_nodes,
            params.platform_anchor_weight,
            params.platform_anchor_value,
        )
        A = A + anchor_A
        b = b + anchor_b

    # Redressement des directions descendantes (voir Decision 3bis) : une
    # direction de croissance qui pointe vers le bas empile les couches du
    # haut vers le bas. Desactive par defaut.
    # La BORNE de plongee remplace l'interdit binaire quand elle est
    # renseignee : les deux repondent a la meme question (que faire d'une
    # direction descendante) et se contrediraient s'ils s'appliquaient tous
    # les deux. Voir `max_downward_growth_angle_deg`.
    upward_flip_mask = None
    if params.max_downward_growth_angle_deg is not None:
        tet_directions, upward_flip_mask = _clamp_downward_directions(
            tet_directions,
            params.up_vector,
            params.max_downward_growth_angle_deg,
            params.upward_growth_epsilon,
        )
    elif params.enforce_upward_growth:
        tet_directions, upward_flip_mask = _enforce_upward_directions(
            tet_directions, params.up_vector, params.upward_growth_epsilon
        )

    if params.radial_tilt_origin is not None:
        tet_directions, _radial_projected = _project_to_reachable_radial_tilt(
            tet_mesh, tet_directions, params.radial_tilt_origin,
            params.up_vector if params.up_vector is not None else np.array([0.0, 0.0, 1.0]),
        )

    # Boucle de retroaction de surplomb : ferme le circuit ouvert decrit sur
    # `FieldSolverParams.support_feedback_iterations`. Elle assemble
    # elle-meme le terme de direction a chaque iteration, d'ou le branchement
    # exclusif ci-dessous plutot qu'un post-traitement.
    support_feedback_diagnostic = None
    if params.support_feedback_iterations > 0:
        A, b, tet_directions, support_feedback_diagnostic = _apply_support_feedback(
            A.tocsr(),
            b,
            tet_mesh,
            shape_gradients,
            volumes,
            tet_directions,
            tet_weights,
            mesh_data.face_normals[nearest_face_per_tet],
            params,
            geodesic_targets,
        )
    else:
        gradient_A, gradient_b = _assemble_gradient_direction_term(
            tet_mesh,
            shape_gradients,
            volumes,
            np.ones(tet_mesh.tet_count, dtype=bool),
            tet_directions,
            tet_weights,
            pattern=fem_pattern,
        )
        A = A + gradient_A
        b = b + gradient_b

    # Terme SOURCE geodesique (voir `geodesic_source_weight`) : ADDITIF au
    # terme de direction principal ci-dessus, jamais une substitution. Meme
    # energie de correspondance de gradient (_assemble_gradient_direction_term,
    # reutilisee telle quelle), avec la direction geodesique comme cible et un
    # poids INDEPENDANT de `tet_weights`. Assemble APRES la boucle de
    # retroaction (qui peut avoir modifie A/b et tet_directions) plutot que
    # dans l'un ou l'autre de ses branches : ce terme est orthogonal aux deux,
    # il s'ajoute au systeme final quel que soit le chemin qui l'a produit.
    geodesic_saddle_tet_count = None
    if params.geodesic_source_weight > 0.0 and geodesic_directions_for_source is not None:
        # Exclusion des points de selle : voir `_geodesic_saddle_mask`. Sans
        # elle, mesure sur real_overhang.stl, aucun poids testable ne donne a
        # la fois des couches plates dans le corps du pilier ET un residu de
        # planeite faible (le compromis global degrade l'un pour ameliorer
        # l'autre, sans jamais les ameliorer ensemble).
        saddle_mask = _geodesic_saddle_mask(
            tet_mesh, geodesic_directions_for_source, params.geodesic_saddle_threshold_deg
        )
        geodesic_saddle_tet_count = int(saddle_mask.sum())
        geo_weights = np.where(
            saddle_mask, 0.0, params.geodesic_source_weight
        ).astype(np.float64)
        geo_A, geo_b = _assemble_gradient_direction_term(
            tet_mesh,
            shape_gradients,
            volumes,
            np.ones(tet_mesh.tet_count, dtype=bool),
            geodesic_directions_for_source,
            geo_weights,
        )
        A = A + geo_A
        b = b + geo_b

    # Diagnostic mesure sur le champ REELLEMENT resolu (pas sur les cibles) :
    # redresser les cibles ne garantit pas mathematiquement que le gradient
    # resultant soit montant partout - le lissage FEM peut encore produire un
    # gradient localement descendant. Le residu est donc mesure, jamais
    # suppose, et rapporte explicitement.
    upward_growth_diagnostic = None
    probe_values = None
    if params.enforce_upward_growth or params.max_downward_growth_angle_deg is not None:
        probe_values = _solve_linear_system(A.tocsr(), b, params)
        gradient_up = _build_gradient_up_operator(
            tet_mesh, shape_gradients, params.up_vector
        )
        fraction, count = _descending_volume_fraction(
            gradient_up @ probe_values, volumes, params.upward_growth_epsilon
        )
        upward_growth_diagnostic = UpwardGrowthDiagnostic(
            flipped_tet_count=int(np.count_nonzero(upward_flip_mask)),
            flipped_volume_fraction=float(
                volumes[upward_flip_mask].sum() / volumes.sum()
            ) if volumes.sum() > 0.0 else 0.0,
            residual_descending_volume_fraction=fraction,
            residual_descending_tet_count=count,
        )

    ctx = ConstraintAssemblyContext(
        tet_mesh=tet_mesh,
        mesh_data=mesh_data,
        overhang_map=overhang_map,
        region_dependency_graph=region_dependency_graph,
        node_kdtree=KDTree(tet_mesh.tet_vertices),
        shape_gradients=shape_gradients,
        tet_volumes=volumes,
        params=params,
    )

    return FieldSolverStaticSystem(
        A_structural=A.tocsr(),
        b_structural=b,
        ctx=ctx,
        upward_growth_diagnostic=upward_growth_diagnostic,
        support_feedback_diagnostic=support_feedback_diagnostic,
        smoothness_modulation_diagnostic=smoothness_modulation_diagnostic,
        probe_values=probe_values,
    )


def solve_with_constraints(
    static_system: FieldSolverStaticSystem,
    constraints: Optional[list[Constraint]],
    params: FieldSolverParams,
    initial_guess: Optional[np.ndarray] = None,
) -> ScalarField:
    """
    Ajoute les contraintes utilisateur au systeme statique et resout.
    C'est le chemin rapide apres l'edition d'un point d'attache : seules
    les contraintes sont reassemblees, la partie structurelle est reutilisee
    telle quelle.

    `initial_guess` (typiquement le champ precedent) accelere la convergence
    du solveur iteratif quand l'edition est une petite perturbation ; sans
    effet sur le chemin direct.
    """
    constraint_A, constraint_b = assemble_all_constraints(constraints, static_system.ctx)

    A = (static_system.A_structural + constraint_A).tocsr()
    b = static_system.b_structural + constraint_b

    # Reutilisation de la sonde du diagnostic de croissance montante.
    #
    # Quand `enforce_upward_growth` (ou la borne de plongee) est actif,
    # `build_static_system` a DEJA resolu le systeme structurel pour mesurer
    # le residu descendant sur le champ reellement resolu. Si les contraintes
    # utilisateur n'apportent rien - cas par defaut du pipeline, `main.py`
    # passant `constraints=None` -, le systeme a resoudre ici est le MEME,
    # au bit pres : meme A, meme b, meme solveur deterministe. Le resoudre une
    # seconde fois coutait 2,0 s sur les 9,45 s de step4 (mesure sur
    # death.stl @ resolution 0.8) pour reproduire un resultat deja disponible.
    #
    # Les trois conditions sont necessaires :
    #   - une sonde existe (le diagnostic etait actif) ;
    #   - les contraintes sont vides OU de contribution strictement nulle,
    #     donc `A` et `b` sont inchanges (teste, jamais suppose) ;
    #   - aucun `initial_guess` n'est fourni : le solveur iteratif partant
    #     d'un autre point ne rend pas le meme vecteur au bit pres, et le
    #     contrat de cet argument est justement d'influencer la resolution.
    reuse_probe = (
        static_system.probe_values is not None
        and initial_guess is None
        and constraint_A.nnz == 0
        and not constraint_b.any()
    )
    if reuse_probe:
        # COPIE, jamais l'alias : `solve_with_constraints` est appelable
        # plusieurs fois sur le meme systeme statique (c'est son role apres
        # edition d'un point d'attache). Rendre le tableau cache tel quel
        # exposerait la sonde a une mutation par l'appelant, qui reviendrait
        # silencieusement dans le resultat de l'appel suivant.
        values = static_system.probe_values.copy()
    else:
        values = _solve_linear_system(A, b, params, initial_guess=initial_guess)
    gradient = _compute_per_tet_gradient(
        static_system.ctx.tet_mesh, static_system.ctx.shape_gradients, values
    )

    return ScalarField(
        values=values,
        gradient=gradient,
        upward_growth_diagnostic=static_system.upward_growth_diagnostic,
        support_feedback_diagnostic=static_system.support_feedback_diagnostic,
        smoothness_modulation_diagnostic=static_system.smoothness_modulation_diagnostic,
    )


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    region_dependency_graph: RegionDependencyGraph,
    constraints: Optional[list[Constraint]],
    params: FieldSolverParams,
) -> ScalarField:
    """
    Assemble et resout le champ scalaire continu sur le TetMesh, en une
    passe. Enveloppe de `build_static_system()` + `solve_with_constraints()`
    pour un slicing complet ; pour une re-resolution rapide apres edition
    d'un point d'attache, appeler ces deux fonctions separement et mettre
    le systeme statique en cache.

    Args:
        tet_mesh: sortie de pipeline/step1_ingest.py.
        mesh_data: sortie de pipeline/step1_ingest.py (pour le mapping
            surface -> volume via son kdtree).
        overhang_map: sortie de pipeline/step2_overhang.py.
        region_dependency_graph: sortie de pipeline/step3_segmentation.py.
        constraints: contraintes utilisateur, potentiellement vide/None.
        params: poids et parametres de resolution.

    Returns:
        ScalarField : une valeur par noeud du TetMesh, un gradient constant
        par tetraedre.

    Raises:
        FieldSolverError: si le solveur ne converge pas ou si le systeme
            assemble n'est pas exploitable.
    """
    # Le champ EST la distance geodesique lisse au plateau : ses isosurfaces
    # sont des coquilles equidistantes A TRAVERS LA MATIERE, donc chaque couche
    # est geodesiquement homogene par construction (etalement nul), et l'ordre
    # d'impression par isovaleur croissante devient la distance elle-meme -
    # une piece suspendue vient apres ce qui la porte parce qu'elle en est plus
    # loin, sans mecanisme d'ordonnancement dedie.
    #
    # Court-circuite deliberement tout l'assemblage de step4 (surplomb,
    # ancrage, contraintes utilisateur) : ce n'est PAS un terme de plus dans le
    # meme systeme, c'est un autre champ. Les Constraint[] utilisateur ne sont
    # donc pas honorees dans ce mode - limite assumee, signalee ici plutot que
    # silencieuse, et raison pour laquelle il reste opt-in.
    if params.use_heat_geodesic_field:
        shape_gradients, volumes = _compute_shape_gradients(tet_mesh)
        platform_nodes = _identify_platform_nodes(
            tet_mesh, mesh_data, region_dependency_graph
        )
        values = _heat_method_distance(
            tet_mesh, shape_gradients, volumes, platform_nodes, params
        )
        return ScalarField(
            values=values,
            gradient=_compute_per_tet_gradient(tet_mesh, shape_gradients, values),
        )

    static_system = build_static_system(
        tet_mesh, mesh_data, overhang_map, region_dependency_graph, params
    )
    return solve_with_constraints(static_system, constraints, params)
