"""
scripts/export_slicing_preview.py

Export GLB du resultat de SLICING SEUL (step1 -> step8), en court-circuitant
volontairement step9 (collision) et step10 (ajustement).

Pourquoi ce script existe a cote de main.py : la detection de collision
produit actuellement un faux positif structurel a CHAQUE point de
trajectoire (la pointe de buse, origine du modele de tete, touche par
construction la matiere deposee et le plateau - voir CLAUDE.md, point de
reprise "zone de garde nozzle_contact_clearance"). Tant que ce point n'est
pas traite, step9/step10 rendent le pipeline complet inexploitable pour un
simple retour visuel, alors que les couches et la trajectoire, elles, sont
deja calculees correctement en une trentaine de secondes.

Ce script ne remplace pas main.py et ne duplique aucune logique de pipeline :
il appelle les memes fonctions run() dans le meme ordre, et reutilise
file_io.visualization_export tel quel.

Usage :
    python scripts/export_slicing_preview.py <mesh.stl> [-o sortie.glb]
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.mesh_types import (  # noqa: E402
    AdjustmentReport,
    CollisionReport,
    SupportCheckReport,
)
from file_io.visualization_export import export_visualization_glb  # noqa: E402
from main import PipelineParams, PipelineResult  # noqa: E402
from pipeline.step1_ingest import run as ingest_run  # noqa: E402
from pipeline.step1b_junction_detection import detect_junction_zones  # noqa: E402
from pipeline.step2_overhang import run as overhang_run  # noqa: E402
from pipeline.step3_segmentation import run as segmentation_run  # noqa: E402
from pipeline.step4_field_solver import run as field_solver_run  # noqa: E402
from pipeline.step5_layers import run as layers_run  # noqa: E402
from pipeline.step5b_collision_order import run as collision_order_run  # noqa: E402
from pipeline.step6_toolpath2d import run as toolpath2d_run  # noqa: E402
from pipeline.step7_scheduling import run as scheduling_run  # noqa: E402
from pipeline.step8_machine_toolpath import gradient_direction_field
from pipeline.step8_machine_toolpath import run as machine_toolpath_run  # noqa: E402
from pipeline.step8b_travel_routing import route_long_travels  # noqa: E402
from pipeline.step9b_support_check import run as support_check_run  # noqa: E402


def run_slicing_only(mesh_path: str, params: PipelineParams) -> PipelineResult:
    """
    step1 -> step8, puis un CollisionReport/AdjustmentReport vides explicites.

    Les deux rapports vides ne pretendent pas qu'aucune collision n'existe :
    ils disent qu'aucune verification n'a eu lieu (`nozzle_model_checked` et
    `chuck_checked` a False, `converged` a False) - meme discipline que le
    reste du projet, une donnee explicite plutot qu'un silence.

    Reproduit la MEME orchestration de double passe legere que
    main.run_pipeline pour params.junction_refinement_resolution
    (CLAUDE.md limite #47/#48) - sans elle, ce script ignorait
    silencieusement ce parametre (bug reel mesure : deux runs consecutifs
    avec/sans junction_refinement_resolution actif produisaient un
    tet_count et un nombre d'anomalies IDENTIQUES, alors que le mecanisme
    fonctionne correctement quand il est reellement exerce).
    """
    t0 = time.time()
    from main import _effective_field_solver_params, _effective_toolpath2d_params, check_model_fits_machine
    check_model_fits_machine(mesh_path, params.machine_profile)

    junction_zones = None
    if params.junction_refinement_resolution is not None:
        print(f"[preview] Passe legere - detection des zones de jonction : {mesh_path}")
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
        light_region_graph = segmentation_run(
            light_mesh_data, light_overhang_map, light_segmentation_params, light_tet_mesh
        )
        junction_zones = detect_junction_zones(
            light_mesh_data, light_region_graph, params.junction_refinement_angle_deg,
            params.junction_refinement_margin,
        )
        print(f"[preview]   {len(junction_zones)} zone(s) de jonction detectee(s)")

    print(f"[preview] Etape 1/8 - ingestion : {mesh_path}")
    mesh_data, tet_mesh, ingest_report = ingest_run(
        mesh_path,
        params.resolution,
        near_platform_resolution=params.near_platform_resolution,
        near_platform_margin=params.near_platform_margin,
        junction_zones=junction_zones,
        junction_refinement_resolution=params.junction_refinement_resolution,
        place_on_platform=params.place_on_platform,
    )

    print("[preview] Etape 2/8 - detection des surplombs")
    overhang_map = overhang_run(mesh_data, params.overhang)

    print("[preview] Etape 3/8 - segmentation en regions")
    # step1 a deja resolu platform_z (fourni ou detecte automatiquement) -
    # propage cette MEME valeur a step3, meme correctif que main.py.
    segmentation_params = (
        replace(params.segmentation, platform_z=ingest_report.platform_z_used)
        if params.segmentation.platform_z is None
        else params.segmentation
    )
    region_graph = segmentation_run(mesh_data, overhang_map, segmentation_params)

    print("[preview] Etape 4/8 - resolution du champ scalaire")
    _fs_params = _effective_field_solver_params(params)
    scalar_field = field_solver_run(tet_mesh, mesh_data, overhang_map, region_graph, _fs_params)

    print("[preview] Etape 5/8 - extraction des couches courbes")
    layers = layers_run(tet_mesh, mesh_data, scalar_field, region_graph, params.layers)
    if params.collision_order is not None:
        layers, order_report = collision_order_run(layers, params.machine_profile, params.collision_order)
        print(f"[preview]   ordre contraint par la tete : {order_report.constraints} contrainte(s), "
              f"{order_report.moved_layers} couche(s) deplacee(s), {order_report.cycle_breaks} violee(s)")
    print(f"[preview]   {len(layers)} couche(s)")

    print("[preview] Etape 6/8 - remplissage 2D par couche")
    layers = toolpath2d_run(layers, mesh_data, overhang_map, _effective_toolpath2d_params(params))

    print("[preview] Etape 7/8 - couleur + ordonnancement")
    colored_layers, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, params.color_map, region_graph, scalar_field, params.scheduling
    )
    print(f"[preview]   {len(schedule.passes)} passe(s) d'impression")

    print("[preview] Etape 8/8 - trajectoire machine (resolution IK)")
    toolpath, diagnostics = machine_toolpath_run(
        colored_layers, schedule, params.machine_profile, params.machine_toolpath,
        direction_field=(gradient_direction_field(tet_mesh, scalar_field,
                                                  nodal=params.machine_toolpath.per_point_normal_nodal)
                         if params.machine_toolpath.per_point_normal else None),
    )
    print(
        f"[preview]   {len(toolpath.points)} point(s) de trajectoire, "
        f"{len(diagnostics)} orientation(s) inatteignable(s)"
    )

    travel_routing_diagnostics = []
    if params.travel_routing is not None:
        # 🔴 Palliatif pour CLAUDE.md limite #46 - meme mecanisme que
        # main.run_pipeline, voir pipeline/step8b_travel_routing.py.
        print("[preview] Etape 8b/8 - routage des trajets a vide longs")
        toolpath, travel_routing_diagnostics = route_long_travels(
            toolpath, params.machine_profile, params.collision, params.travel_routing
        )
        print(
            f"[preview]   {len(toolpath.points)} point(s) de trajectoire apres routage, "
            f"{len(travel_routing_diagnostics)} trajet(s) encore en collision apres routage"
        )

    support_check_report = SupportCheckReport(checked=False)
    if params.support_check is not None:
        # Meme mecanisme et meme position relative que main.run_pipeline
        # (apres step8b, avant step9) : step9b ne depend ni de step9 ni de
        # step10, donc il tourne REELLEMENT ici, contrairement au
        # CollisionReport/AdjustmentReport vides retournes plus bas.
        print("[preview] Etape 8c/8 - verification du support de depot")
        support_check_report = support_check_run(
            toolpath,
            params.machine_profile,
            params.layers.layer_height_base,
            ingest_report.platform_z_used,
            params.support_check,
        )
        print(
            f"[preview]   {len(support_check_report.events)} depot(s) sans support "
            f"({support_check_report.unsupported_fraction:.2%} de la matiere extrudee) : "
            f"{support_check_report.ordering_defect_count} defaut(s) d'ordre, "
            f"{support_check_report.true_overhang_count} surplomb(s) reel(s)"
        )

    print(f"[preview] step1->step8 termine en {time.time() - t0:.1f}s")
    print("[preview] step9/step10 VOLONTAIREMENT SAUTES (voir docstring)")

    return PipelineResult(
        mesh_data=mesh_data,
        tet_mesh=tet_mesh,
        overhang_map=overhang_map,
        region_dependency_graph=region_graph,
        scalar_field=scalar_field,
        layers=layers,
        colored_layers=colored_layers,
        schedule=schedule,
        toolpath=toolpath,
        step8_diagnostics=diagnostics,
        travel_routing_diagnostics=travel_routing_diagnostics,
        support_check_report=support_check_report,
        collision_report=CollisionReport(
            events=[], chuck_checked=False, nozzle_model_checked=False
        ),
        adjusted_toolpath=toolpath,  # aucun ajustement applique
        adjustment_report=AdjustmentReport(
            entries=[], converged=False, unresolved_collisions=[], iterations_used=0
        ),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh_path")
    parser.add_argument("-o", "--output-glb", default="slicing_preview.glb")
    parser.add_argument("--resolution", type=float, default=None)
    parser.add_argument(
        "--layer-stride", type=int, default=1, help="1 = toutes les couches, N = une sur N"
    )
    parser.add_argument(
        "--no-color-by-print-order", action="store_true",
        help="Desactive le degrade par ordre d'impression sur le mesh source "
             "(bleu #0000ff = premier, jaune #fff700 = dernier)",
    )
    parser.add_argument(
        "--show-toolpath", action="store_true",
        help="Reaffiche la trajectoire dans le GLB (desactivee par defaut).",
    )
    args = parser.parse_args(argv)

    params = PipelineParams()
    if args.resolution is not None:
        params.resolution = args.resolution

    result = run_slicing_only(args.mesh_path, params)

    print(f"[preview] export GLB -> {args.output_glb}")
    export_visualization_glb(
        result,
        args.output_glb,
        layer_stride=args.layer_stride,
        color_by_print_order=not args.no_color_by_print_order,
        show_toolpath=args.show_toolpath,
    )
    size = os.path.getsize(args.output_glb)
    print(f"[preview] Visualisation ecrite : {os.path.abspath(args.output_glb)} ({size} octets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
