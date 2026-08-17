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
from dataclasses import dataclass, field
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
    TetMesh,
    ToolpathMachine,
)
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run
from pipeline.step4_field_solver import FieldSolverParams, run as field_solver_run
from pipeline.step5_layers import LayerParams, run as layers_run
from pipeline.step6_toolpath2d import Step6Params, run as toolpath2d_run
from pipeline.step7_scheduling import Step7Params, run as scheduling_run
from pipeline.step8_machine_toolpath import Step8Params, run as machine_toolpath_run
from pipeline.step9_collision_sim import CollisionSimParams, run as collision_run
from pipeline.step10_adjustment import AdjustmentReport, Step10Params, run as adjustment_run


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
            RotationAxisSpec(axis=head_axis, angle_range_deg=(0.0, 180.0)),
        ],
        platform_origin=np.array([0.0, 0.0, 0.0]),
        nozzle_model=nozzle_model,
        nozzle_diameter=0.4,
        max_angular_velocity=0.0,
        max_angular_acceleration=0.0,
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
    overhang: OverhangParams = field(default_factory=OverhangParams)
    # direction_deviation_weight active explicitement ici (defaut 0.0 dans
    # SegmentationParams lui-meme, jamais change silencieusement pour les
    # appelants existants) : sans lui, une seule direction compatible avec
    # TOUTES les faces d'une piece (ex: un L) suffit a fusionner toute la
    # piece en une region unique, meme inclinee a plus de 50 deg de la
    # verticale la ou aucun surplomb ne le justifie (mesure concretement
    # sur res/test_models/l_shape_overhang.stl : base verticale a 54.8 deg
    # sans raison avant ce fix). Valeur 0.2 calibree par balayage direct
    # (0.1 a 0.5 tous robustes sur ce cas), avec region_count_penalty=0.005
    # (recalibre dans SegmentationParams meme, voir CLAUDE.md) pour que la
    # scission base/bras se maintienne. Pipeline reel = premier appelant a
    # activer ce mecanisme, apres validation par
    # tests/test_step3_segmentation.py::test_l_shape_with_deviation_weight_splits_base_from_arm.
    segmentation: SegmentationParams = field(
        default_factory=lambda: SegmentationParams(direction_deviation_weight=0.2)
    )
    field_solver: FieldSolverParams = field(default_factory=FieldSolverParams)
    constraints: Optional[list[Constraint]] = None
    layers: LayerParams = field(default_factory=LayerParams)
    toolpath2d: Step6Params = field(default_factory=Step6Params)
    scheduling: Step7Params = field(default_factory=Step7Params)
    color_map: Optional[ColorMap] = None
    machine_toolpath: Step8Params = field(default_factory=Step8Params)
    # path_sampling_step/safety_margin n'ont PAS de defaut dans
    # CollisionSimParams (step9) - l'orchestrateur doit choisir une valeur
    # raisonnable plutot que de laisser l'appelant en decouvrir l'absence.
    collision: CollisionSimParams = field(
        default_factory=lambda: CollisionSimParams(path_sampling_step=2.0, safety_margin=1.0)
    )
    adjustment: Step10Params = field(default_factory=Step10Params)
    machine_profile: MachineProfile = field(default_factory=_default_machine_profile)


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
    toolpath: ToolpathMachine                  # AVANT ajustement step10
    step8_diagnostics: list[Step8Diagnostic]
    collision_report: CollisionReport
    adjusted_toolpath: ToolpathMachine
    adjustment_report: AdjustmentReport


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
    print(f"[main] Etape 1/10 - ingestion : {mesh_path}")
    mesh_data, tet_mesh, ingest_report = ingest_run(
        mesh_path, resolution=params.resolution, allow_multi_mesh=params.allow_multi_mesh
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
    region_dependency_graph = segmentation_run(mesh_data, overhang_map, params.segmentation)
    print(f"[main]   {len(region_dependency_graph.regions)} region(s)")

    print("[main] Etape 4/10 - resolution du champ scalaire")
    scalar_field = field_solver_run(
        tet_mesh, mesh_data, overhang_map, region_dependency_graph, params.constraints, params.field_solver
    )

    print("[main] Etape 5/10 - extraction des couches courbes")
    layers = layers_run(tet_mesh, mesh_data, scalar_field, region_dependency_graph, params.layers)
    print(f"[main]   {len(layers)} couche(s)")

    print("[main] Etape 6/10 - remplissage 2D par couche")
    layers = toolpath2d_run(layers, mesh_data, overhang_map, params.toolpath2d)

    print("[main] Etape 7/10 - couleur + ordonnancement")
    colored_layers, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, params.color_map, region_dependency_graph, scalar_field, params.scheduling
    )
    print(f"[main]   {len(schedule.passes)} passe(s) d'impression")

    print("[main] Etape 8/10 - trajectoire machine (resolution IK)")
    toolpath, step8_diagnostics = machine_toolpath_run(
        colored_layers, schedule, params.machine_profile, params.machine_toolpath
    )
    print(
        f"[main]   {len(toolpath.points)} point(s) de trajectoire, "
        f"{len(step8_diagnostics)} orientation(s) inatteignable(s)"
    )

    print("[main] Etape 9/10 - simulation de collision")
    collision_report = collision_run(toolpath, params.machine_profile, params.collision)
    print(
        f"[main]   {len(collision_report.events)} collision(s) detectee(s) "
        f"(nozzle_model_checked={collision_report.nozzle_model_checked}, "
        f"chuck_checked={collision_report.chuck_checked})"
    )

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
    args = parser.parse_args(argv)

    params = PipelineParams(resolution=args.resolution)
    params.layers = LayerParams(layer_height_base=args.layer_height)

    result = run_pipeline(args.mesh_path, params)

    if args.output_glb:
        from file_io.visualization_export import export_visualization_glb

        export_visualization_glb(result, args.output_glb)
        print(f"[main] Visualisation ecrite : {args.output_glb}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
