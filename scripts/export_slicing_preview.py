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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.mesh_types import AdjustmentReport, CollisionReport  # noqa: E402
from file_io.visualization_export import export_visualization_glb  # noqa: E402
from main import PipelineParams, PipelineResult  # noqa: E402
from pipeline.step1_ingest import run as ingest_run  # noqa: E402
from pipeline.step2_overhang import run as overhang_run  # noqa: E402
from pipeline.step3_segmentation import run as segmentation_run  # noqa: E402
from pipeline.step4_field_solver import run as field_solver_run  # noqa: E402
from pipeline.step5_layers import run as layers_run  # noqa: E402
from pipeline.step6_toolpath2d import run as toolpath2d_run  # noqa: E402
from pipeline.step7_scheduling import run as scheduling_run  # noqa: E402
from pipeline.step8_machine_toolpath import run as machine_toolpath_run  # noqa: E402


def run_slicing_only(mesh_path: str, params: PipelineParams) -> PipelineResult:
    """
    step1 -> step8, puis un CollisionReport/AdjustmentReport vides explicites.

    Les deux rapports vides ne pretendent pas qu'aucune collision n'existe :
    ils disent qu'aucune verification n'a eu lieu (`nozzle_model_checked` et
    `chuck_checked` a False, `converged` a False) - meme discipline que le
    reste du projet, une donnee explicite plutot qu'un silence.
    """
    t0 = time.time()

    print(f"[preview] Etape 1/8 - ingestion : {mesh_path}")
    mesh_data, tet_mesh, _report = ingest_run(mesh_path, params.resolution)

    print("[preview] Etape 2/8 - detection des surplombs")
    overhang_map = overhang_run(mesh_data, params.overhang)

    print("[preview] Etape 3/8 - segmentation en regions")
    region_graph = segmentation_run(mesh_data, overhang_map, params.segmentation)

    print("[preview] Etape 4/8 - resolution du champ scalaire")
    scalar_field = field_solver_run(
        tet_mesh, mesh_data, overhang_map, region_graph, params.constraints, params.field_solver
    )

    print("[preview] Etape 5/8 - extraction des couches courbes")
    layers = layers_run(tet_mesh, mesh_data, scalar_field, region_graph, params.layers)
    print(f"[preview]   {len(layers)} couche(s)")

    print("[preview] Etape 6/8 - remplissage 2D par couche")
    layers = toolpath2d_run(layers, mesh_data, overhang_map, params.toolpath2d)

    print("[preview] Etape 7/8 - couleur + ordonnancement")
    colored_layers, schedule = scheduling_run(
        layers, mesh_data, tet_mesh, params.color_map, region_graph, scalar_field, params.scheduling
    )
    print(f"[preview]   {len(schedule.passes)} passe(s) d'impression")

    print("[preview] Etape 8/8 - trajectoire machine (resolution IK)")
    toolpath, diagnostics = machine_toolpath_run(
        colored_layers, schedule, params.machine_profile, params.machine_toolpath
    )
    print(
        f"[preview]   {len(toolpath.points)} point(s) de trajectoire, "
        f"{len(diagnostics)} orientation(s) inatteignable(s)"
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
    args = parser.parse_args(argv)

    params = PipelineParams()
    if args.resolution is not None:
        params.resolution = args.resolution

    result = run_slicing_only(args.mesh_path, params)

    print(f"[preview] export GLB -> {args.output_glb}")
    export_visualization_glb(result, args.output_glb, layer_stride=args.layer_stride)
    size = os.path.getsize(args.output_glb)
    print(f"[preview] Visualisation ecrite : {os.path.abspath(args.output_glb)} ({size} octets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
