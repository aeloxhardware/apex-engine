"""
main.py

Orchestrateur : enchaine step1_ingest a step10_adjustment sur un fichier
STL/3MF reel, avec des parametres par defaut raisonnables pour chaque
etape (les memes que ceux deja actes dans chaque `stepN_*.Params`, jamais
reinventes ici).

Ne contient AUCUNE logique metier propre - un simple relais entre les 10
etapes, chacune restant l'unique source de verite pour son propre
comportement. Les exceptions metier (MeshIngestError, SegmentationError,
etc.) ne sont jamais interceptees : elles remontent telles quelles.

Usage CLI :
    python main.py res/test_models/l_shape_overhang.stl --output-glb out.glb

Le `resegmentation_hook` de l'etape 10 est laisse a None ici : l'implementer
reellement demanderait de re-jouer step3->step8 avec exclusion d'angles
(voir CLAUDE.md, "Detail de pipeline/step10_adjustment.py", limite #23) -
hors perimetre d'un premier passage dont l'objectif est de VOIR le
resultat. Sans hook, toute collision non resolue par reorientation locale
bascule directement en `fallback_support`, un comportement deja explicite
et teste, pas un trou dans le pipeline.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np

from core.mesh_types import (
    CollisionReport,
    Constraint,
    ColorMap,
    Layer,
    MachineProfile,
    MeshData,
    OverhangMap,
    PrintSchedule,
    RegionDependencyGraph,
    RotationAxisSpec,
    ScalarField,
    Step8Diagnostic,
    SupportCheckReport,
    TetMesh,
    ToolpathMachine,
    TravelRoutingDiagnostic,
)
from pipeline.step1_ingest import check_fits_platform
from pipeline.step1_ingest import run as ingest_run
from pipeline.step1b_junction_detection import detect_junction_zones
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import FieldSolverParams, run as field_solver_run
from pipeline.step5_layers import LayerParams, run as layers_run
from pipeline.step5b_collision_order import CollisionOrderParams, run as collision_order_run
from pipeline.step6_toolpath2d import Step6Params, run as toolpath2d_run
from pipeline.step7_scheduling import Step7Params, run as scheduling_run
from pipeline.step8_machine_toolpath import gradient_direction_field
from pipeline.step8_machine_toolpath import Step8Params, run as machine_toolpath_run
from pipeline.step8b_travel_routing import TravelRoutingParams, route_long_travels
from pipeline.step9b_support_check import SupportCheckParams, run as support_check_run
from pipeline.step9_collision_sim import CollisionSimParams, run as collision_run
from pipeline.step10_adjustment import (AdjustmentReport, Step10Params, clamp_head_height, reorient_coupled_tiers,
                                        run as adjustment_run)




def _effective_field_solver_params(params: "PipelineParams"):
    """FieldSolverParams effectifs : `radial_tilt_origin` derive de la machine
    quand `reachable_tilt_from_machine` est actif et que le premier axe est
    couple a la position. Copie : `params` n'est jamais mute."""
    from dataclasses import replace as _replace
    fs = params.field_solver
    m = params.machine_profile
    if (params.reachable_tilt_from_machine and fs.radial_tilt_origin is None
            and list(getattr(m, "position_coupled_axes", None) or []) == [0]
            and m.platform_origin is not None):
        return _replace(fs, radial_tilt_origin=np.asarray(m.platform_origin, dtype=float))
    return fs

def _effective_toolpath2d_params(params: "PipelineParams"):
    """Step6Params effectifs : `spiral_axis_origin` (centre de la spirale
    OSCILLATING_SPIRAL) derive de `MachineProfile.platform_origin`, source de
    verite unique de l'axe du plateau. Copie : `params` n'est jamais mute."""
    from dataclasses import replace as _replace
    tp = params.toolpath2d
    origin = params.machine_profile.platform_origin
    if tp.spiral_axis_origin is None and origin is not None:
        return _replace(tp, spiral_axis_origin=tuple(float(x) for x in origin))
    return tp


def check_model_fits_machine(mesh_path: str, machine_profile: MachineProfile) -> None:
    """Refus DES L'ENTREE d'un modele qui ne tient sur le plateau dans aucune
    position (ModelTooLargeError), avant toute tetraedrisation. Placement-
    independant : la position sur le plateau est decidee par l'interface."""
    radius = getattr(machine_profile, "platform_radius", None)
    bv = machine_profile.build_volume
    max_height = None if bv is None else float(np.asarray(bv.max_corner)[2] - np.asarray(bv.min_corner)[2])
    if radius is None and max_height is None:
        return
    import trimesh as _trimesh
    loaded = _trimesh.load(mesh_path, force="mesh")
    check_fits_platform(np.asarray(loaded.vertices), radius, max_height)

def _default_machine_profile() -> MachineProfile:
    """
    Cinematique par defaut = la machine actuelle documentee dans CLAUDE.md
    (plateau C tour complet + tete B in [0,180], perpendiculaires) - UN
    exemple d'usage concret parmi d'autres cinematiques possibles (limite
    de conception #18 de CLAUDE.md), pas la reference figee du type
    MachineProfile lui-meme.

    Charge aetheris-4.glb comme nozzle_model si present ; repli sur None
    (log, pas de crash) sinon - step9 degrade gracieusement dans ce cas
    (`nozzle_model_checked=False`), voir CLAUDE.md.
    """
    import os

    platform_axis = np.array([0.0, 0.0, 1.0])
    head_axis = np.array([1.0, 0.0, 0.0])

    nozzle_model = None
    glb_path = os.path.join(
        os.path.dirname(__file__), "res", "3d", "toolheads", "aetheris-4", "aetheris-4.glb"
    )
    if os.path.isfile(glb_path):
        try:
            from file_io.mesh_loader import load_toolhead_mesh

            nozzle_model = load_toolhead_mesh(glb_path)
        except Exception as exc:  # pragma: no cover - chemin de repli, pas une logique testee
            print(f"[main] Avertissement : echec du chargement de {glb_path} ({exc}) - "
                  f"simulation de collision desactivee (nozzle_model=None).", file=sys.stderr)

    return MachineProfile(
        rotation_axes=[
            RotationAxisSpec(axis=platform_axis, angle_range_deg=(0.0, 360.0)),
            RotationAxisSpec(axis=head_axis, angle_range_deg=(0.0, 185.0)),
        ],
        platform_origin=np.array([0.0, 0.0, 0.0]),
        nozzle_model=nozzle_model,
        nozzle_diameter=0.4,
        # Limites de la tete (axe B). Fichier RRF ACTIF : config.g (confirme par
        # l'utilisateur, 2026-10-06), M203 B100000 (/min) = 1667 deg/s et M201 B1300
        # = 1300 deg/s2 (unites supposees en degres). 4axisinput.g (360 deg/s,
        # 2000 deg/s2) et les valeurs de 2026-09-28 (750, 300) ne sont pas actifs.
        # Le profil n'a qu'une valeur pour tous les axes de rotation ; C, couple a la
        # position, n'est jamais lisse (limite #59) : elles ne s'appliquent qu'a B.
        # Non modelises : Y (600 deg/s, 4000 deg/s2) et Z (6,7 mm/s, 70 mm/s2).
        max_angular_velocity=float(np.radians(1666.7)),
        max_angular_acceleration=float(np.radians(1300.0)),
        # Machine ACTUELLE : R C Z B, pas d'axe Y -- C est impose par la
        # position du point (limite #59). +R = -Y a angles nuls : B > 0 (axe
        # X) penche la normale de couche vers -Y, c'est-a-dire vers +R (le
        # surplomb s'etend vers l'exterieur, convention verifiee avec
        # l'utilisateur, cf. scripts/gen_elbow_tube_gcode.py). Apres l'upgrade
        # (axe lineaire de plus) : position_coupled_axes=[].
        position_coupled_axes=[0],
        radial_direction=np.array([0.0, -1.0, 0.0]),
        # Plateau de 300 mm de diametre (utilisateur, 2026-09-29).
        platform_radius=150.0,
    )


@dataclass
class PipelineParams:
    """
    Un champ par etape. Chaque valeur par defaut est celle DEJA ACTEE dans
    le `Params` de l'etape correspondante - ce fichier ne reinvente aucun
    defaut, il les expose juste au niveau orchestrateur/CLI.
    """
    resolution: float = 1.5  # partage step1 (densite de tetraedrisation) et FieldSolverParams.resolution
    allow_multi_mesh: bool = False
    # near_platform_resolution/margin calibres sur res/test_models/l_shape_overhang.stl
    # (CLAUDE.md) : l'erreur angulaire pres du plateau (voir field_solver
    # ci-dessous) est une erreur de discretisation FEM, pas de ponderation -
    # un tetraedre P1 assez grand pour s'etendre du contact plateau jusqu'a
    # une zone ou la direction cible a deja beaucoup tourne ne peut
    # representer qu'un seul gradient constant. Raffiner LOCALEMENT
    # (near_platform_resolution=0.5, margin=5.0, via TetGen.add_region -
    # bgmesh/metric teste et ecarte, crash natif reproductible) donne une
    # amelioration angulaire quasi identique a un raffinement global complet
    # (resolution=0.5 : 4.77 deg moyenne) - 4.91 deg moyenne mesuree, pour
    # ~5.4x moins de tetraedres (1677 contre 9088).
    near_platform_resolution: Optional[float] = 0.5
    near_platform_margin: float = 5.0
    # Raffinement local de maillage aux jonctions multi-region (CLAUDE.md
    # limite #47/#48) : generalise near_platform_resolution ci-dessus
    # au-dela du seul cas plateau, corrige la non-monotonie du champ
    # scalaire dans les zones ou plusieurs regions a directions divergentes
    # sont adjacentes (motif festonne dense, mesure sur death.stl). Opt-in,
    # None = DESACTIVE par defaut, meme discipline que near_platform_resolution
    # (qui lui est actif par defaut - la difference est deliberee : ce
    # mecanisme ajoute une passe legere supplementaire de step1+2+3 avant la
    # vraie tetraedrisation, cout non nul meme minime, jamais impose par
    # defaut sans mesure prealable sur le modele de l'appelant).
    junction_refinement_resolution: Optional[float] = None
    junction_refinement_angle_deg: float = 30.0
    junction_refinement_margin: float = 3.0
    # Pose le modele sur le plateau (translation en z seule, min(z) -> 0) a
    # l'ingestion : le plateau de la machine est a z = 0 et l'export refuse
    # une pointe sous 0. Voir step1_ingest.run(place_on_platform=...).
    # ACTIVE par defaut (2026-10-13) : mesure sur les 4 modeles de
    # res/test_models/ (out/measure_placed_2026-10-13.jsonl) - l_shape et
    # l-triple deja a z = 0 (identiques au bit), death (min z = 51,4 mm)
    # 0,22 -> 0,21 % de matiere sans appui, turkish (30,3 mm) 0,17 -> 0,19 %,
    # collisions de la tete 118 -> 116 et 1406 -> 1431, 0 collision plateau.
    # Sans pose, death et turkish flottaient au-dessus du plateau : collision
    # plateau jamais exercee, et un G-code qui demarre en l'air.
    place_on_platform: bool = True
    overhang: OverhangParams = field(default_factory=OverhangParams)
    # anchor_growth_weight VOLONTAIREMENT NON ACTIVE ici (defaut de la
    # dataclass = None). Il l'a ete brievement, puis retire APRES MESURE sur
    # un 4e modele - le detail vaut d'etre garde, c'est une erreur de methode
    # facile a refaire :
    #   - death.stl                      : 0.53% -> 0.47% de matiere sans
    #                                      appui (-11% relatif) ;
    #   - l_shape_overhang.stl           : strictement neutre ;
    #   - l-triple-axis-overhanged-shape : strictement neutre ;
    #   - turkish.stl                    : 2.78% -> 4.89% (**+76%**),
    #                                      defauts d'ordre 5475 -> 7715.
    # Trois modeles sur quatre disaient "ameliorant ou neutre" : c'est ce
    # qui a fait croire a un reglage universel. Le quatrieme montre que non.
    # Le mecanisme reste donc OPT-IN, a mesurer modele par modele - meme
    # discipline que `support_driven_order`, et exactement la mise en garde
    # de la limite #50 ("l'activer partout repeterait l'erreur des cinq
    # tentatives precedentes"). Voir limite #51 pour ce qu'il corrige
    # reellement et pourquoi il n'est pas generalisable.
    segmentation: SegmentationParams = field(default_factory=SegmentationParams)
    # direction_blend_radius/power calibres sur res/test_models/l_shape_overhang.stl
    # (CLAUDE.md) : sans melange (0.0), l'angle du plan local d'une couche
    # proche du plateau atteignait 33 deg en moyenne (jusqu'a 42.7 deg) -
    # avec radius=10.0/power=2.0, le GRADIENT REEL du champ (pas
    # l'approximation de plan local de step5) est mesure parfaitement
    # vertical (0.0 deg) sur les tetraedres les plus proches du plateau ;
    # l'angle residuel mesure sur le plan local de step5 (~10 deg) est un
    # artefact de sa propre moyenne sur toute une couche, pas du melange
    # lui-meme (voir CLAUDE.md, hors perimetre de ce fix).
    # enforce_upward_growth active ici (defaut de la dataclass = False), meme
    # precedent que direction_blend_radius : une direction de croissance
    # descendante empile les couches du haut vers le bas et chaque couche
    # atterrit dans le vide. Mesure sur death.stl : 67.6% -> 0.06% du volume
    # dont le gradient descend, a cout de calcul quasi inchange.
    #
    # DESACTIVATION TENTEE PUIS ANNULEE (session de debug de l'ordre des
    # couches), sur l'argument - correct en soi - que la tete de la machine
    # reelle couvre 0-185 deg et PEUT donc deposer vers le bas, rendant
    # caduque la justification d'origine ("physiquement irrealisable").
    # Mesure du desactiver sur death.stl, qui invalide la tentative :
    #   - le gradient suit alors bien la cible de step3 (ecart 158.5 -> 29.4
    #     deg) : le solveur n'a jamais eu de bug, il obeissait au flip ;
    #   - MAIS la couche la plus basse remonte de z=51.49 a z=53.19 pour un
    #     plateau a z=51.40 : plus AUCUNE couche ne touche le plateau, et
    #     100% d'entre elles (contre 86.8%) sont `support_unresolved`.
    # Cause identifiee : la cible descendante que step3 assigne au dome
    # (region 3, [-0.26, 0.08, -0.96]) est elle-meme physiquement absurde
    # POUR CE MODELE - la suivre voudrait dire commencer a imprimer en l'air
    # au sommet du capot (z~100) puis descendre, alors que le capot ne tient
    # au pilier que par sa BASE (z~59, couches 34-36 mesurees a 0.12-0.36 mm
    # du pilier, toutes les suivantes s'en ecartant jusqu'a 36 mm). Sa
    # direction de croissance correcte pointe vers l'EXTERIEUR depuis cette
    # jonction, pas vers le bas. `enforce_upward_growth` reste donc actif
    # comme garde-fou contre une cible de region mal choisie ; le vrai fix
    # est dans step3 (choix de direction tenant compte de l'ancrage, pas
    # seulement du surplomb local) - chantier distinct, non traite ici.
    field_solver: FieldSolverParams = field(
        default_factory=lambda: FieldSolverParams(
            direction_blend_radius=10.0,
            direction_blend_power=2.0,
            # VOIE PRINCIPALE (decision de l'utilisateur, 2026-10-25) : champ = distance au
            # plateau par la methode de la chaleur (limite #79), flux projete sur les
            # inclinaisons atteignables, couches ordonnees par isovaleur (variante
            # `heat_iso_reach` de scripts/measure_matrix.py). FEM directionnel et Bird retires.
            enforce_upward_growth=False,
            use_heat_geodesic_field=True,
            heat_project_reachable=True,
        )
    )
    constraints: Optional[list[Constraint]] = None
    layers: LayerParams = field(default_factory=lambda: LayerParams(isovalue_order=True))
    # out_of_plane_interpolation_power active par defaut ici (defaut de la
    # dataclass = None/desactive), meme precedent que field_solver ci-dessus
    # qui force direction_blend_radius=10.0 : sans lui, les points generes par
    # shapely (parois n>=1, infill) restent colles au plan local MOYEN de la
    # couche et peuvent tomber HORS du solide la ou la couche est fortement
    # non plane. Mesure sur res/test_models/l-triple-axis-overhanged-shape.stl :
    # 74 points debordants (max 0.390 mm) sans interpolation, contre 0
    # (max 0.002 mm) avec power=2.0 - power=1.0 laisse 20 points, power=3.0
    # en laisse 1, d'ou ce choix. Le contour n=0 n'est jamais concerne : ses
    # offsets restent exacts (limite #31).
    # split_rings_into_runs actif par defaut EN PLUS de
    # out_of_plane_interpolation_power (dont il depend - sans interpolation,
    # les anneaux generes sont deja plats et ne produisent qu'un troncon).
    # Corrige la cause racine mesuree du residu du mecanisme 2 (step5,
    # reordonnancement de COUCHES ENTIERES) : 70% des defauts d'ordre
    # restants sur death.stl etaient INTRA-anneau (deux sommets consecutifs
    # d'une meme paroi sur une couche courbe), hors de portee de step5 - voir
    # pipeline/step6_toolpath2d.py Decision 5. L'ORDRE des troncons entre eux
    # n'est PAS decide (aucun mecanisme actif actuellement) : une tentative de
    # tri par altitude locale a ete mesuree et s'est averee une regression
    # (voir historique de session) - retiree, pas reintroduite sans nouvelle
    # mesure positive.
    # continuous_rings (2026-10-01) : le decoupage en troncons laissait le
    # segment de jonction NON extrude (7,7 % de la longueur deposee sur death,
    # 15,7 % sur l_shape : des trous dans les parois), et le contour n=0 sans
    # son segment de fermeture. Mesure step9b sur les 4 modeles : death -15 %,
    # turkish -2,6 %, l-triple +0,8 %, l_shape +3 % ; active comme correction
    # d'un defaut physique (CLAUDE.md, limite #65).
    toolpath2d: Step6Params = field(
        default_factory=lambda: Step6Params(
            out_of_plane_interpolation_power=2.0, split_rings_into_runs=True, continuous_rings=True,
            run_split_min_amplitude=0.05,
        )
    )
    scheduling: Step7Params = field(default_factory=Step7Params)
    color_map: Optional[ColorMap] = None
    machine_toolpath: Step8Params = field(default_factory=Step8Params)
    # path_sampling_step/safety_margin n'ont PAS de defaut dans
    # CollisionSimParams (step9) - l'orchestrateur doit choisir une valeur
    # raisonnable plutot que de laisser l'appelant en decouvrir l'absence.
    collision: CollisionSimParams = field(
        # fast_kdtree (limite #62) : mesure contre l'exact sur les 4 modeles le
        # 2026-10-01 (l-triple et l_shape entiers, turkish et death par
        # echantillon) : 0 faux positif, 0 faux negatif, ecart <= 0,087 mm ;
        # death 16 s au lieu de plusieurs heures.
        default_factory=lambda: CollisionSimParams(path_sampling_step=2.0, safety_margin=1.0, fast_kdtree=True)
    )
    adjustment: Step10Params = field(default_factory=Step10Params)
    # Borne de l'angle de tete selon la hauteur au-dessus du plateau, AVANT
    # step9 (CLAUDE.md limite #62, step10.clamp_head_height). Opt-in, non
    # mesure sur les quatre modeles de res/test_models/.
    height_clamp: bool = False
    # Reorientation par paliers des points en contact, machine a axe couple
    # (limite #62, step10.reorient_coupled_tiers), APRES step9 et avant
    # step10.run ; step9 est rejoue apres chaque tour. Opt-in (None).
    coupled_reorientation_tiers: Optional[tuple[float, ...]] = None      # ex. (45.0, 90.0, 181.0)
    coupled_reorientation_rounds: int = 3
    machine_profile: MachineProfile = field(default_factory=_default_machine_profile)
    # Routage des trajets a vide longs (CLAUDE.md limite #46) : palliatif
    # lift-and-traverse valide par collision, PAS un vrai path-planning 3D
    # (voir pipeline/step8b_travel_routing.py). Opt-in, None = DESACTIVE
    # par defaut - meme discipline que junction_refinement_resolution :
    # ce mecanisme ajoute un cout de calcul non nul (re-teste chaque trajet
    # a vide long contre l'index de matiere imprimee), jamais impose sans
    # activation explicite de l'appelant.
    travel_routing: Optional[TravelRoutingParams] = None
    # step5b (2026-10-01) : ordre des couches contraint par la tete - la
    # couche sur laquelle la tete toucherait une autre passe avant elle.
    # ACTIVE par defaut (CLAUDE.md, limite #67) : collisions de la vraie tete
    # death 41,6 -> 1,2 %, turkish 53,6 -> 0,6 % ; matiere sans appui turkish
    # 0,70 -> 0,25 % mais death 0,227 -> 0,311 % (A VALIDER par l'utilisateur) ;
    # l_shape et l-triple inchanges. None = ordre de step5 seul.
    collision_order: Optional[CollisionOrderParams] = None  # step5b inadapte au champ de la chaleur (limite #79)
    # Verification du support de depot sur la trajectoire machine reelle
    # (voir pipeline/step9b_support_check.py) : quelle matiere est extrudee
    # sans rien de deja imprime en dessous, et s'agit-il d'un defaut d'ORDRE
    # (le support existe, depose trop tard) ou d'un surplomb REEL.
    # Opt-in, None = DESACTIVE par defaut - meme discipline que
    # travel_routing : le mecanisme ajoute un cout non nul (une passe sur
    # tous les points en extrusion), jamais impose sans activation explicite.
    support_check: Optional[SupportCheckParams] = None
    # Directions cibles de step4 projetees sur les inclinaisons ATTEIGNABLES
    # quand le premier axe est couple a la position (machine actuelle, limite
    # #59) : voir FieldSolverParams.radial_tilt_origin et limite #64. ACTIVE
    # par defaut apres mesure sur les QUATRE modeles de res/test_models/
    # (matiere sans appui : l_shape -9 %, l-triple -45 %, death -23 %,
    # turkish -30 % ; ecart tete/normale median 16-49 deg -> 3-6 deg).
    # Sans effet sur une machine a C libre (position_coupled_axes vide) ou si
    # radial_tilt_origin est deja renseigne.
    reachable_tilt_from_machine: bool = True


@dataclass
class PipelineResult:
    """
    Tout ce que le pipeline a produit, etape par etape - pour qu'un export
    de visualisation (ou un futur appelant UI) n'ait pas besoin de
    re-derouler le pipeline pour retrouver une sortie intermediaire.
    """
    mesh_data: MeshData
    tet_mesh: TetMesh
    overhang_map: OverhangMap
    region_dependency_graph: RegionDependencyGraph
    scalar_field: ScalarField
    layers: list[Layer]                       # sortie step6 (avant recoloration/decoupage step7)
    colored_layers: list[Layer]                # sortie step7 (ce que step8 a reellement consomme)
    schedule: PrintSchedule
    toolpath: ToolpathMachine                  # AVANT ajustement step10 (mais APRES routage des trajets a vide, si actif)
    step8_diagnostics: list[Step8Diagnostic]
    travel_routing_diagnostics: list[TravelRoutingDiagnostic]
    support_check_report: SupportCheckReport
    collision_report: CollisionReport
    adjusted_toolpath: ToolpathMachine
    adjustment_report: AdjustmentReport


def export_gcode_checked(toolpath: ToolpathMachine, params: PipelineParams, out_path: str,
                         pivot_to_tip: Optional[float] = None, check_travel: bool = True):
    """Export G-code controle (comme scripts/export_model_gcode.py) : pointe jamais
    sous le plateau (`strict`), hauteur mini de la VRAIE tete selon B pour les
    degagements, L = derniere valeur de config/pivot_to_tip.json sinon 46 mm.
    `check_travel` : degagements controles contre la matiere deja deposee
    (`material_travel_check`, comme la chaine de recherche), en natif si actif."""
    from file_io.gcode_export import (GcodeParams, export_gcode, head_zlow_from_table, material_travel_check,
                                      toolpath_arrays)
    from file_io.pivot_length import current_pivot_to_tip
    from pipeline import step9_collision_sim as s9
    machine, sim = params.machine_profile, params.collision
    if pivot_to_tip is None:
        pivot_to_tip, _origin = current_pivot_to_tip()
    head_zlow = None
    if machine.nozzle_model is not None:
        head = s9._filter_nozzle_contact_zone(s9._surface_sample_points(machine.nozzle_model, sim.toolhead_sample_spacing),
                                              sim.nozzle_contact_clearance)
        lo_b, hi_b = machine.rotation_axes[1].angle_range_deg
        grid_b = np.arange(lo_b, hi_b + 1e-9, 0.5)
        zlow = np.array([s9._positioned_toolhead_vertices(head, np.zeros(3), [0.0, float(b)], machine.rotation_axes)
                         [:, 2].min() for b in grid_b])
        head_zlow = head_zlow_from_table(grid_b, zlow)
    gp = GcodeParams(pivot_to_tip=pivot_to_tip, strict=True)
    arrays = toolpath_arrays(toolpath)
    travel_check = material_travel_check(arrays, machine, gp) if check_travel else None
    return export_gcode(arrays, machine, gp, out_path, head_zlow=head_zlow, travel_check=travel_check)


def run_pipeline(mesh_path: str, params: PipelineParams) -> PipelineResult:
    """
    Enchaine step1 a step10 dans l'ordre exact du pipeline. Trois
    dependances multi-etapes a respecter (deja verifiees, pas a
    re-decouvrir) :
      - step4/step5/step7 ont besoin de `mesh_data`/`tet_mesh` de step1
        DIRECTEMENT, pas seulement de la sortie de l'etape precedente.
      - step8 doit consommer `colored_layers` (1er element du tuple
        retourne par step7), PAS les `layers` d'entree de step7.
      - step9 et step10 partagent la MEME instance de `machine_profile` et
        de `CollisionSimParams` (`params.collision`).

    Ne catch aucune exception metier : MeshIngestError, SegmentationError,
    LayerExtractionError, FieldSolverError, SchedulingError,
    ToolpathResolutionError remontent telles quelles a l'appelant.
    """
    check_model_fits_machine(mesh_path, params.machine_profile)
    junction_zones = None
    if params.junction_refinement_resolution is not None:
        # Passe LEGERE, jetable : mesh + step2 (overhang) + step3
        # (segmentation) SEULEMENT, jamais la tetraedrisation fine ni
        # step4-10 - juste pour localiser les frontieres inter-region a
        # forte divergence directionnelle AVANT la vraie tetraedrisation
        # (CLAUDE.md limite #47/#48). Chicken-and-egg resolu ainsi : step2/
        # step3 ne dependent que de MeshData (le mesh de SURFACE), jamais du
        # TetMesh volumique - `_light_tet_mesh` ci-dessous n'est jamais
        # utilise. Chaque step reste pur (run(input, params) -> output) ;
        # c'est UNIQUEMENT cette fonction d'orchestration qui sait qu'elle
        # rejoue une partie du pipeline, jamais un step lui-meme.
        print(f"[main] Passe legere - detection des zones de jonction : {mesh_path}")
        light_mesh_data, light_tet_mesh, light_ingest_report = ingest_run(
            mesh_path, resolution=params.resolution, allow_multi_mesh=params.allow_multi_mesh,
            place_on_platform=params.place_on_platform,
        )
        light_overhang_map = overhang_run(light_mesh_data, params.overhang)
        light_segmentation_params = (
            replace(params.segmentation, platform_z=light_ingest_report.platform_z_used)
            if params.segmentation.platform_z is None
            else params.segmentation
        )
        # light_tet_mesh EST utilise si branch_cut_threshold est actif (la
        # segmentation par branche a besoin du volume) - inchange sinon.
        light_region_graph = segmentation_run(
            light_mesh_data, light_overhang_map, light_segmentation_params, light_tet_mesh
        )
        junction_zones = detect_junction_zones(
            light_mesh_data, light_region_graph, params.junction_refinement_angle_deg,
            params.junction_refinement_margin,
        )
        print(f"[main]   {len(junction_zones)} zone(s) de jonction detectee(s)")

    print(f"[main] Etape 1/10 - ingestion : {mesh_path}")
    mesh_data, tet_mesh, ingest_report = ingest_run(
        mesh_path,
        resolution=params.resolution,
        allow_multi_mesh=params.allow_multi_mesh,
        near_platform_resolution=params.near_platform_resolution,
        near_platform_margin=params.near_platform_margin,
        junction_zones=junction_zones,
        junction_refinement_resolution=params.junction_refinement_resolution,
        place_on_platform=params.place_on_platform,
    )
    print(
        f"[main]   watertight avant/apres reparation : "
        f"{ingest_report.was_watertight_initially}/{ingest_report.is_watertight_after_repair}, "
        f"{ingest_report.tetrahedra_generated} tetraedres"
    )

    print("[main] Etape 2/10 - detection des surplombs")
    overhang_map = overhang_run(mesh_data, params.overhang)
    print(f"[main]   {int(np.sum(overhang_map.is_overhang))} faces en surplomb sur {mesh_data.face_count}")

    print("[main] Etape 3/10 - segmentation en regions")
    # step1 a deja resolu platform_z (fourni ou detecte automatiquement,
    # IngestReport.platform_z_used) - propage cette MEME valeur a step3
    # plutot que de la laisser re-deriver independamment (source de
    # divergence silencieuse potentielle si les deux calculs different),
    # sauf si l'appelant a explicitement fixe params.segmentation.platform_z.
    segmentation_params = (
        replace(params.segmentation, platform_z=ingest_report.platform_z_used)
        if params.segmentation.platform_z is None
        else params.segmentation
    )
    region_dependency_graph = segmentation_run(mesh_data, overhang_map, segmentation_params, tet_mesh)
    print(f"[main]   {len(region_dependency_graph.regions)} region(s)")

    print("[main] Etape 4/10 - resolution du champ scalaire")
    _fs_params = _effective_field_solver_params(params)
    if _fs_params.planar_parts:
        from pipeline.step4c_planar_parts import build_planar_parts
        tet_mesh, scalar_field = build_planar_parts(tet_mesh, mesh_data, region_dependency_graph, _fs_params,
                                                    overhang_map)
    else:
        scalar_field = field_solver_run(
            tet_mesh, mesh_data, overhang_map, region_dependency_graph, params.constraints, _fs_params
        )

    print("[main] Etape 5/10 - extraction des couches courbes")
    layers = layers_run(tet_mesh, mesh_data, scalar_field, region_dependency_graph, params.layers)
    if params.collision_order is not None:
        layers, order_report = collision_order_run(layers, params.machine_profile, params.collision_order)
        print(f"[main]   ordre contraint par la tete : {order_report.constraints} contrainte(s), "
              f"{order_report.moved_layers} couche(s) deplacee(s), {order_report.cycle_breaks} violee(s)")
    print(f"[main]   {len(layers)} couche(s)")

    print("[main] Etape 6/10 - remplissage 2D par couche")
    layers = toolpath2d_run(layers, mesh_data, overhang_map, _effective_toolpath2d_params(params))

    print("[main] Etape 7/10 - couleur + ordonnancement")
    colored_layers, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, params.color_map, region_dependency_graph, scalar_field, params.scheduling
    )
    print(f"[main]   {len(schedule.passes)} passe(s) d'impression")

    print("[main] Etape 8/10 - trajectoire machine (resolution IK)")
    toolpath, step8_diagnostics = machine_toolpath_run(
        colored_layers, schedule, params.machine_profile, params.machine_toolpath,
        direction_field=(gradient_direction_field(tet_mesh, scalar_field,
                                                  nodal=params.machine_toolpath.per_point_normal_nodal)
                         if params.machine_toolpath.per_point_normal else None),
    )
    print(
        f"[main]   {len(toolpath.points)} point(s) de trajectoire, "
        f"{len(step8_diagnostics)} orientation(s) inatteignable(s)"
    )

    travel_routing_diagnostics: list[TravelRoutingDiagnostic] = []
    if params.travel_routing is not None:
        # 🔴 Palliatif pour CLAUDE.md limite #46 (trajets a vide sans
        # routage) - lift-and-traverse valide par collision, PAS un vrai
        # path-planning 3D. Voir pipeline/step8b_travel_routing.py.
        print("[main] Etape 8b/10 - routage des trajets a vide longs")
        toolpath, travel_routing_diagnostics = route_long_travels(
            toolpath, params.machine_profile, params.collision, params.travel_routing
        )
        print(
            f"[main]   {len(toolpath.points)} point(s) de trajectoire apres routage, "
            f"{len(travel_routing_diagnostics)} trajet(s) encore en collision apres routage"
        )

    support_check_report = SupportCheckReport(checked=False)
    if params.support_check is not None:
        # Place APRES step8b (la trajectoire routee est celle reellement
        # executee) et AVANT step9 : ne depend ni du CollisionReport ni de
        # step10, qui ne modifie jamais `position` - seulement `axis_angles`.
        print("[main] Etape 8c/10 - verification du support de depot")
        support_check_report = support_check_run(
            toolpath,
            params.machine_profile,
            params.layers.layer_height_base,
            ingest_report.platform_z_used,
            params.support_check,
        )
        print(
            f"[main]   {len(support_check_report.events)} depot(s) sans support "
            f"({support_check_report.unsupported_fraction:.2%} de la matiere extrudee) : "
            f"{support_check_report.ordering_defect_count} defaut(s) d'ordre, "
            f"{support_check_report.true_overhang_count} surplomb(s) reel(s)"
        )

    if params.height_clamp:
        toolpath, clamped, _deviations = clamp_head_height(toolpath, params.machine_profile, params.collision)
        print(f"[main]   borne de l'angle de tete selon la hauteur : {clamped} point(s) ramene(s)")

    print("[main] Etape 9/10 - simulation de collision")
    collision_report = collision_run(toolpath, params.machine_profile, params.collision)
    print(
        f"[main]   {len(collision_report.events)} collision(s) detectee(s) "
        f"(nozzle_model_checked={collision_report.nozzle_model_checked}, "
        f"chuck_checked={collision_report.chuck_checked})"
    )
    if params.coupled_reorientation_tiers and collision_report.events:
        for _round in range(params.coupled_reorientation_rounds):
            contacts = sorted({e.point_index for e in collision_report.events})
            toolpath, reorient_log = reorient_coupled_tiers(
                toolpath, contacts, params.machine_profile, params.collision,
                tiers_deg=tuple(params.coupled_reorientation_tiers),
            )
            collision_report = collision_run(toolpath, params.machine_profile, params.collision)
            print(f"[main]   reorientation par paliers : {sum(1 for _i, b, _t in reorient_log if b is not None)}"
                  f"/{len(contacts)} point(s) resolu(s) ; {len(collision_report.events)} collision(s) restante(s)")
            if not collision_report.events:
                break
            if not any(found is not None for _i, found, _t in reorient_log):
                break   # aucun point resolu : trajectoire inchangee, le tour suivant serait identique

    print("[main] Etape 10/10 - ajustement")
    adjusted_toolpath, adjustment_report = adjustment_run(
        toolpath,
        collision_report,
        params.machine_profile,
        params.collision,
        params.adjustment,
        step8_diagnostics=step8_diagnostics,
        resegmentation_hook=None,
    )
    print(
        f"[main]   converged={adjustment_report.converged}, "
        f"{len(adjustment_report.entries)} ajustement(s), "
        f"{len(adjustment_report.unresolved_collisions)} collision(s) non resolue(s)"
    )

    return PipelineResult(
        mesh_data=mesh_data,
        tet_mesh=tet_mesh,
        overhang_map=overhang_map,
        region_dependency_graph=region_dependency_graph,
        scalar_field=scalar_field,
        layers=layers,
        colored_layers=colored_layers,
        schedule=schedule,
        toolpath=toolpath,
        step8_diagnostics=step8_diagnostics,
        travel_routing_diagnostics=travel_routing_diagnostics,
        support_check_report=support_check_report,
        collision_report=collision_report,
        adjusted_toolpath=adjusted_toolpath,
        adjustment_report=adjustment_report,
    )


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Orchestrateur du pipeline de slicing 4 axes (step1 -> step10)."
    )
    parser.add_argument("mesh_path", help="Fichier STL/3MF source.")
    parser.add_argument(
        "--resolution", type=float, default=1.5, help="Taille cible de tetraedre (mm). Defaut : 1.5."
    )
    parser.add_argument(
        "--layer-height", type=float, default=0.2, help="Hauteur de couche cible (mm). Defaut : 0.2."
    )
    parser.add_argument(
        "--output-glb", default=None,
        help="Chemin du fichier .glb de visualisation a produire (mesh + couches + trajectoire).",
    )
    parser.add_argument(
        "--no-color-by-print-order", action="store_true",
        help="Desactive le degrade par ordre d'impression sur le mesh source "
             "(bleu #0000ff = imprime en premier, jaune #fff700 = en dernier) "
             "et revient au gris semi-transparent.",
    )
    parser.add_argument(
        "--show-toolpath", action="store_true",
        help="Reaffiche la trajectoire finale ajustee dans le GLB (desactivee par defaut).",
    )
    parser.add_argument(
        "--output-gcode", default=None,
        help="G-code a produire (trajectoire ajustee, export controle : pointe jamais sous le plateau).",
    )
    args = parser.parse_args(argv)

    params = PipelineParams(resolution=args.resolution)
    params.layers = LayerParams(layer_height_base=args.layer_height)

    result = run_pipeline(args.mesh_path, params)
    if args.output_gcode:
        _lines, report = export_gcode_checked(result.adjusted_toolpath, params, args.output_gcode)
        print(f"[main] G-code ecrit : {args.output_gcode} ({report.moves} mouvements)")

    if args.output_glb:
        from file_io.visualization_export import export_visualization_glb

        export_visualization_glb(
            result,
            args.output_glb,
            color_by_print_order=not args.no_color_by_print_order,
            show_toolpath=args.show_toolpath,
        )
        print(f"[main] Visualisation ecrite : {args.output_glb}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
