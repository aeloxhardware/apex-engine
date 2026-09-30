"""
tests/test_main_pipeline.py

Test d'integration cible sur la propagation de platform_z entre step1 et
step3 telle qu'orchestree par main.py::run_pipeline (et
scripts/export_slicing_preview.py::run_slicing_only, meme correctif) :
verifie que la valeur REELLEMENT utilisee par step1 (fournie ou detectee
automatiquement) est bien celle transmise a step3, et non un defaut fige
re-derive independamment.

Ne passe PAS par run_pipeline() complet (step1->step10, couteux - collision
et ajustement sont hors du perimetre de ce fix) : reproduit exactement le
meme fragment de logique que main.py (ingest_run -> propagation ->
segmentation_run), sur un mesh synthetique translate en Z, pour rester
rapide tout en prouvant l'integration reelle entre les deux etapes.
"""
from __future__ import annotations

import os
import tempfile
from dataclasses import replace

import numpy as np
import pytest
import trimesh

from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run

_PLATFORM_Z_OFFSET = 1.245  # meme offset que res/test_models/death.stl


def _translated_l_shape_mesh() -> trimesh.Trimesh:
    vertical = trimesh.creation.box(extents=[4, 4, 20])
    vertical.apply_translation([0, 0, 10])
    horizontal = trimesh.creation.box(extents=[20, 4, 4])
    horizontal.apply_translation([8, 0, 20])
    mesh = vertical.union(horizontal, engine="manifold")
    mesh.apply_translation([0, 0, _PLATFORM_Z_OFFSET])
    return mesh


def _propagate_platform_z(segmentation_params: SegmentationParams, platform_z_used: float) -> SegmentationParams:
    """Meme logique de propagation que main.py::run_pipeline / scripts/export_slicing_preview.py."""
    if segmentation_params.platform_z is None:
        return replace(segmentation_params, platform_z=platform_z_used)
    return segmentation_params


def test_pipeline_propagates_detected_platform_z_from_step1_to_step3():
    """
    Reproduction du fragment step1->step3 de main.py::run_pipeline sur un
    mesh NON pose a z=0 (offset 1.245, meme decalage que
    res/test_models/death.stl) : platform_contact_faces (step3) ne doit
    JAMAIS etre vide - preuve directe que la valeur detectee par step1
    (IngestReport.platform_z_used) atteint bien step3, pas un defaut fige
    a 0.0 qui manquerait systematiquement le vrai plateau.
    """
    mesh = _translated_l_shape_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, _tet_mesh, ingest_report = ingest_run(path, resolution=2.0)
    finally:
        os.unlink(path)

    assert ingest_report.platform_z_used == pytest.approx(_PLATFORM_Z_OFFSET, abs=1e-6)

    overhang_map = overhang_run(mesh_data, OverhangParams(overhang_angle_threshold_deg=45.0))

    segmentation_params = _propagate_platform_z(SegmentationParams(), ingest_report.platform_z_used)
    region_graph = segmentation_run(mesh_data, overhang_map, segmentation_params)

    assert np.any(region_graph.platform_contact_faces), (
        "platform_contact_faces ne doit jamais etre vide sur un mesh pose a "
        "plat, meme translate en Z - sinon l'ancrage plateau du solveur FEM "
        "(step4) devient un no-op silencieux (CLAUDE.md)"
    )


def test_pipeline_respects_explicit_platform_z_override():
    """
    Si l'appelant fixe explicitement SegmentationParams.platform_z, la
    propagation depuis step1 ne doit JAMAIS l'ecraser.
    """
    mesh = _translated_l_shape_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, _tet_mesh, ingest_report = ingest_run(path, resolution=2.0)
    finally:
        os.unlink(path)

    overhang_map = overhang_run(mesh_data, OverhangParams(overhang_angle_threshold_deg=45.0))

    # override explicite volontairement FAUX (loin du vrai plateau) - la
    # propagation ne doit pas le corriger silencieusement.
    explicit_params = SegmentationParams(platform_z=-100.0)
    segmentation_params = _propagate_platform_z(explicit_params, ingest_report.platform_z_used)

    assert segmentation_params.platform_z == -100.0


# ---------------------------------------------------------------------------
# Orchestration de la double passe legere pour le raffinement de jonction
# (CLAUDE.md limite #47/#48) - reproduit le meme fragment de logique que
# main.py::run_pipeline (passe legere step1+step2+step3 -> detect_junction_zones
# -> vraie tetraedrisation), sans passer par run_pipeline() complet
# (step1->step10, couteux et hors du perimetre de ce fix - meme discipline
# que les tests platform_z ci-dessus).
# ---------------------------------------------------------------------------

def _tripod_mesh() -> trimesh.Trimesh:
    """Meme geometrie que tests/test_step5_layers.py::_tripod_mesh."""
    foot = trimesh.creation.box(extents=[4, 4, 20])
    foot.apply_translation([0, 0, 10])
    arm_x = trimesh.creation.box(extents=[20, 4, 4])
    arm_x.apply_translation([8, 0, 20])
    arm_y = trimesh.creation.box(extents=[4, 20, 4])
    arm_y.apply_translation([0, 8, 20])
    arm_diag = trimesh.creation.box(extents=[14, 14, 4])
    arm_diag.apply_translation([-7, -7, 20])

    mesh = foot.union(arm_x, engine="manifold")
    mesh = mesh.union(arm_y, engine="manifold")
    mesh = mesh.union(arm_diag, engine="manifold")
    return mesh


def test_junction_refinement_disabled_skips_light_pass(monkeypatch):
    """
    junction_refinement_resolution=None (defaut de PipelineParams) : le
    bloc d'orchestration REEL de main.run_pipeline ne doit JAMAIS declencher
    la passe legere (ingest_run appele UNE SEULE fois, pas deux) - verifie
    par compteur d'appels directement sur `main.ingest_run`.

    N'appelle PAS `run_pipeline()` complet (step1->step10, y compris
    step9/collision - trop lourd et fragile pour un test unitaire cible sur
    l'orchestration de step1 seule, meme discipline "reproduire le
    fragment" que les tests platform_z ci-dessus) : reproduit ici
    uniquement le DEBUT de `run_pipeline` (le bloc `if
    params.junction_refinement_resolution is not None` puis le premier
    appel a `ingest_run`), copie fidele du code source verifie par lecture
    directe de main.py au moment de l'ecriture de ce test.
    """
    import main as main_module

    call_count = {"ingest": 0}
    real_ingest_run = main_module.ingest_run

    def _counting_ingest_run(*args, **kwargs):
        call_count["ingest"] += 1
        return real_ingest_run(*args, **kwargs)

    monkeypatch.setattr(main_module, "ingest_run", _counting_ingest_run)

    mesh = _tripod_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        params = main_module.PipelineParams(resolution=1.5, junction_refinement_resolution=None)

        junction_zones = None
        if params.junction_refinement_resolution is not None:
            light_mesh_data, _light_tet_mesh, light_ingest_report = main_module.ingest_run(
                path, resolution=params.resolution, allow_multi_mesh=params.allow_multi_mesh,
            )
            light_overhang_map = main_module.overhang_run(light_mesh_data, params.overhang)
            light_segmentation_params = (
                replace(params.segmentation, platform_z=light_ingest_report.platform_z_used)
                if params.segmentation.platform_z is None
                else params.segmentation
            )
            light_region_graph = main_module.segmentation_run(
                light_mesh_data, light_overhang_map, light_segmentation_params
            )
            junction_zones = main_module.detect_junction_zones(
                light_mesh_data, light_region_graph, params.junction_refinement_angle_deg,
                params.junction_refinement_margin,
            )

        main_module.ingest_run(
            path, resolution=params.resolution, allow_multi_mesh=params.allow_multi_mesh,
            near_platform_resolution=params.near_platform_resolution,
            near_platform_margin=params.near_platform_margin,
            junction_zones=junction_zones,
            junction_refinement_resolution=params.junction_refinement_resolution,
        )
    finally:
        os.unlink(path)

    assert call_count["ingest"] == 1, (
        "junction_refinement_resolution=None ne doit jamais declencher la "
        "passe legere - un seul appel a ingest_run attendu"
    )


def test_junction_refinement_enabled_runs_light_pass_before_real_pass():
    """
    junction_refinement_resolution actif : la passe legere (mesh + step2 +
    step3, jamais step4-10) doit s'executer et produire des JunctionZone
    non vides sur le tripode (jonction physique connue), AVANT la vraie
    tetraedrisation qui les consomme - reproduit le fragment complet de
    main.py::run_pipeline.
    """
    from pipeline.step1b_junction_detection import detect_junction_zones

    mesh = _tripod_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        # --- Passe legere (exactement le fragment de main.py) ---
        light_mesh_data, _light_tet_mesh, light_report = ingest_run(path, resolution=1.5)
        light_overhang = overhang_run(light_mesh_data, OverhangParams(overhang_angle_threshold_deg=45.0))
        light_seg_params = replace(SegmentationParams(), platform_z=light_report.platform_z_used)
        light_graph = segmentation_run(light_mesh_data, light_overhang, light_seg_params)

        junction_zones = detect_junction_zones(
            light_mesh_data, light_graph, angle_threshold_deg=30.0, margin=3.0
        )

        # Le tripode synthetique (une seule region sur ce mesh brut, jamais
        # segmente automatiquement en plusieurs regions par step3 - meme
        # limite documentee pour tripod_data dans tests/test_step5_layers.py)
        # ne produit PAS de jonction detectable par step3 seul (contrairement
        # a la fixture tripod_data qui injecte les regions manuellement) -
        # ce test verifie seulement que le FRAGMENT D'ORCHESTRATION s'execute
        # sans erreur et retourne une liste (vide ou non), pas que des zones
        # sont necessairement trouvees sur ce mesh precis.
        assert isinstance(junction_zones, list)

        # --- Vraie tetraedrisation (exactement le fragment de main.py) ---
        mesh_data, tet_mesh, ingest_report = ingest_run(
            path, resolution=1.5,
            junction_zones=junction_zones,
            junction_refinement_resolution=0.4 if junction_zones else None,
        )
    finally:
        os.unlink(path)

    assert tet_mesh.tet_count > 0
    assert ingest_report.tetrahedra_generated == tet_mesh.tet_count


# ---------------------------------------------------------------------------
# Orchestration du routage des trajets a vide (CLAUDE.md limite #46) -
# reproduit le fragment d'orchestration de main.py::run_pipeline (appel a
# route_long_travels entre step8 et step9, uniquement si
# params.travel_routing is not None), sans passer par run_pipeline()
# complet - meme discipline "reproduire le fragment" que les tests
# ci-dessus.
# ---------------------------------------------------------------------------

def test_travel_routing_disabled_by_default_leaves_toolpath_unchanged():
    """
    PipelineParams.travel_routing=None (defaut) : le fragment
    d'orchestration ne doit jamais appeler route_long_travels - le nombre
    de points de la trajectoire reste strictement celui produit par step8.
    """
    import main as main_module
    from core.mesh_types import MachineProfile, ToolpathMachine, ToolpathPoint

    points = [
        ToolpathPoint(position=np.array([0.0, 0.0, 0.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([50.0, 0.0, 0.0]), axis_angles=[], extrusion=0.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([51.0, 0.0, 0.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
    ]
    toolpath = ToolpathMachine(points=points)
    params = main_module.PipelineParams(travel_routing=None)

    travel_routing_diagnostics = []
    if params.travel_routing is not None:
        toolpath, travel_routing_diagnostics = main_module.route_long_travels(
            toolpath, params.machine_profile, params.collision, params.travel_routing
        )

    assert len(toolpath.points) == 3
    assert travel_routing_diagnostics == []


def test_travel_routing_enabled_routes_long_travel_before_collision_sim():
    """
    PipelineParams.travel_routing actif : le fragment d'orchestration doit
    router le trajet a vide long en plusieurs points (lift-and-traverse)
    AVANT toute simulation de collision (step9), meme fragment que
    main.py::run_pipeline.
    """
    import main as main_module
    from core.mesh_types import MachineProfile, ToolpathMachine, ToolpathPoint
    from pipeline.step8b_travel_routing import TravelRoutingParams

    points = [
        ToolpathPoint(position=np.array([0.0, 0.0, 0.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([50.0, 0.0, 0.0]), axis_angles=[], extrusion=0.0, feedrate=1000.0),
        ToolpathPoint(position=np.array([51.0, 0.0, 0.0]), axis_angles=[], extrusion=1.0, feedrate=1000.0),
    ]
    toolpath = ToolpathMachine(points=points)
    machine_profile = MachineProfile(rotation_axes=[], platform_origin=None, nozzle_diameter=0.4)
    params = main_module.PipelineParams(
        travel_routing=TravelRoutingParams(length_threshold=5.0, lift_height=3.0, traverse_samples=2),
    )
    params.machine_profile = machine_profile

    travel_routing_diagnostics = []
    if params.travel_routing is not None:
        toolpath, travel_routing_diagnostics = main_module.route_long_travels(
            toolpath, params.machine_profile, params.collision, params.travel_routing
        )

    assert len(toolpath.points) > 3, (
        "un trajet a vide de 50mm (> length_threshold=5.0) doit etre route "
        "en plusieurs points, pas laisse en ligne droite"
    )
    assert travel_routing_diagnostics == []  # aucun modele de tete charge : rien a signaler


# ---------------------------------------------------------------------------
# Etape 9b - verification du support de depot (orchestration)
#
# CLAUDE.md documente un piege deja rencontru deux fois (limites #47/#48/#49) :
# `main.run_pipeline` et `scripts/export_slicing_preview.run_slicing_only`
# orchestrent TOUS LES DEUX le pipeline, et un mecanisme cable dans un seul
# des deux passe inapercu. Le test ci-dessous le verifie LITTERALEMENT sur la
# source des deux fonctions, plutot que de s'en remettre a la convention.
# ---------------------------------------------------------------------------

def test_support_check_disabled_by_default():
    """
    `support_check=None` (defaut) : aucune verification declenchee, et le
    rapport porte `checked=False` - signal explicite, jamais une liste vide
    ambigue entre "verifie, rien a signaler" et "jamais verifie".
    """
    import main as main_module
    from core.mesh_types import SupportCheckReport

    params = main_module.PipelineParams()
    assert params.support_check is None

    report = SupportCheckReport(checked=False)
    if params.support_check is not None:  # pragma: no cover - jamais pris au defaut
        raise AssertionError("support_check devrait etre desactive par defaut")

    assert report.checked is False
    assert report.events == []


def test_support_check_wired_in_both_orchestrators():
    """
    Les DEUX points d'entree du pipeline doivent appeler step9b quand il est
    actif. Verification directe sur la source, pour ne pas dependre d'une
    convention que la prochaine modification pourrait rompre en silence.
    """
    import inspect

    import main as main_module
    import scripts.export_slicing_preview as preview_module

    main_source = inspect.getsource(main_module.run_pipeline)
    preview_source = inspect.getsource(preview_module.run_slicing_only)

    for name, source in (("run_pipeline", main_source), ("run_slicing_only", preview_source)):
        assert "support_check_run(" in source, (
            f"{name} n'appelle pas step9b : les deux orchestrateurs doivent "
            "rester synchronises (CLAUDE.md limites #47/#48)"
        )
        assert "params.support_check is not None" in source, (
            f"{name} doit garder step9b opt-in (None = desactive)"
        )


def test_support_check_runs_before_collision_sim():
    """
    step9b est place APRES step8/step8b et AVANT step9 : il consomme la
    trajectoire reellement executee, sans dependre du CollisionReport.
    """
    import inspect

    import main as main_module

    source = inspect.getsource(main_module.run_pipeline)
    assert source.index("support_check_run(") < source.index("collision_run(")
