"""
element_ardor/step1_ingest.py

Etape 1 du pipeline : ingestion du mesh source (STL/3MF) et preparation de
toutes les structures dont les etapes suivantes ont besoin.

Cette etape est la seule du pipeline qui touche a l'I/O disque cote geometrie.
Elle ne depend d'aucune autre etape - c'est le point d'entree du pipeline et
le fichier le plus simple a tester en isolation (aucune donnee amont requise,
juste un chemin de fichier).

Sortie : (MeshData, TetMesh, IngestReport)
"""

from __future__ import annotations

import dataclasses
import os
import pickle
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Optional

import math

import numpy as np
import scipy.sparse as sp
import scipy.spatial
import trimesh
import trimesh.transformations as trimesh_tf
from trimesh import geometry as trimesh_geometry
from trimesh import grouping as trimesh_grouping
from trimesh.creation import triangulate_polygon
from trimesh.path import polygons as trimesh_polygons

from core import progress as _progress
from core import _cpp, runtime
from core.mesh_types import JunctionZone, MeshData, TetMesh
from core.spatial import KDTree, RayIntersector
from element_ardor._self_intersection_resolve import resolve_self_intersections

from element_ardor.step1_types import (  # noqa: F401  (réexport : imports existants)
    IngestReport,
    MeshIngestError,
    ModelTooLargeError,
    SurfaceInconsistentError,
    _load_raw_mesh,
    _min_enclosing_circle,
    check_fits_platform,
)
from element_ardor.step1_repair import (  # noqa: F401  (réexport : imports existants)
    _COPLANAR_EPS,
    _DOMINANT_SHELL_VOLUME_RATIO,
    _KEPT_PIECE_VOLUME_RATIO,
    _LOCAL_REPAIR_MAX_GROWTH_ITERATIONS,
    _LOCAL_REPAIR_MAX_REMOVED_RATIO,
    _LOCAL_REPAIR_MIN_TRIANGLE_AREA,
    _boundary_loops_oriented,
    _build_half_edge_map,
    _clean_and_repair,
    _count_self_intersecting_faces,
    _extract_dominant_shell,
    _fill_loop_min_weight,
    _find_self_intersecting_faces,
    _repair_local_self_intersections,
    _repair_with_pymeshfix,
    _triangles_cross,
)
from element_ardor.step1_tetgen_seeds import (  # noqa: F401  (réexport : imports existants)
    _INTERIOR_POINT_INITIAL_SAMPLES,
    _INTERIOR_POINT_MAX_SAMPLES,
    _INTERIOR_POINT_SEED,
    _find_guaranteed_interior_point,
    _find_internal_cavities,
    _find_rest_seed_outside_junction_zones,
    _merged_junction_seed,
    _shared_cut_cap,
    _split_for_platform_refinement,
)
from element_ardor.step1_tetgen_run import (  # noqa: F401  (réexport : imports existants)
    _GEOMETRIC_RESOLUTION_MAX_TETGEN_PROBES,
    _NativeTetgen,
    _RECOVERABLE_TETGEN_ERROR_SUBSTRINGS,
    _TETGEN_BRACKETED_INDICES,
    _TetgenCrashError,
    _TetgenNamedConflictError,
    _build_tet_adjacency,
    _is_recoverable_tetgen_error,
    _native_tetgen_worker_path,
    _parse_tetgen_conflict_vertices,
    _probe_tetgen_conflicts,
    _tetrahedralize_once,
    _tetrahedralize_once_isolated,
    _TetOptions,
)
from element_ardor.step1_build import (  # noqa: F401  (réexport : imports existants)
    TETRANET_CELL_PER_RESOLUTION,
    TETRANET_STATUS_INCONSISTENT_SURFACE,
    _build_acceleration_structures,
    _build_mesh_data,
    _ingest_with_tetranet,
    tetranet_available,
)


# ---------------------------------------------------------------------------
# Etape 2bis - Auto-intersections
#
# TetGen REFUSE ("The input surface mesh contain self-intersections") tout
# maillage dont des triangles se traversent, meme parfaitement watertight,
# a winding coherent, d'un seul tenant et sans face degeneree - c'est
# exactement le cas de res/test_models/death.stl (mesure : 5024 faces
# fautives sur 12662, soit 39.7%, reparties sur tout le modele).
#
# Aucune des reparations deja en place (merge_vertices, nondegenerate_faces,
# unique_faces, fill_holes, fix_normals, fix_winding) ne traite ce defaut :
# elles corrigent la TOPOLOGIE, pas la geometrie qui s'interpenetre. Une
# union booleenne manifold3d ne le traite pas non plus (mesure : rend le
# maillage a l'identique, meme nombre de faces, meme volume).
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Etape 6 - Tetraedrisation volumique (tetgen)
# ---------------------------------------------------------------------------


def _finish_geometric_resolution(
    resolved: trimesh.Trimesh,
    resolution: float,
) -> trimesh.Trimesh:
    """
    Affine iterativement le resultat de `resolve_self_intersections` avec TetGen LUI-MEME comme oracle
    (limite #43). A chaque iteration : vraie tetraedrisation (`_tetrahedralize_once_isolated`) ; si TetGen la
    refuse pour une raison RECUPERABLE et nomme des sommets en conflit (`_probe_tetgen_conflicts`), cible une
    reparation locale (`_repair_local_self_intersections`, `extra_seed_faces`) sur EXACTEMENT ce qu'il refuse.
    Mesure sur turkish.stl : 8 auto-intersections locales invisibles au test SAT global ; apres reparation
    ciblee, TetGen aboutit sur la geometrie COMPLETE (330306 tetraedres, z=[30.26,90.78]).

    Deliberement SANS raffinement plateau/jonction : le PLC combine renumerote les sommets
    (`_split_for_platform_refinement`), un conflit nomme ne se traduirait donc pas en `extra_seed_faces` sur le
    mesh de surface repare ici (bug mesure : aucun sommet exploitable). La question est « ce mesh de SURFACE
    est-il tetraedrisable EN SOI » ; la compatibilite avec le raffinement plateau reste a la cascade de
    `_tetrahedralize` (repli « sans raffinement plateau », limite #44).

    Best-effort, ne leve jamais : rend `(mesh, TetMesh de la sonde ou None)` ; le meilleur mesh obtenu des
    qu'un probe reussit, qu'une reparation echoue, qu'aucun sommet n'est nommable ou que le budget est epuise.
    """
    current = resolved
    for _probe in range(_GEOMETRIC_RESOLUTION_MAX_TETGEN_PROBES):
        try:
            probe = _tetrahedralize_once_isolated(
                current, resolution, None, 5.0, 0.0, None, None,
            )
            # deja verifie tetraedrisable tel quel (sans raffinement) ; le
            # TetMesh de la sonde est RENDU : la tentative « sans raffinement
            # plateau » de `_tetrahedralize` refait exactement ce calcul
            # (meme maillage, meme resolution ; marge et platform_z ignores
            # sans raffinement). Mesure turkish.stl : 7,26 s economisees.
            return current, probe
        except _TetgenCrashError:
            return current, None  # rien a cibler depuis un crash natif
        except RuntimeError as exc:
            if not _is_recoverable_tetgen_error(exc):
                return current, None
            # `exc.tetgen_output` (capture depuis la VRAIE tentative) s'est
            # mesuree peu fiable - tronquee avant les avertissements
            # nommes eux-memes sur ce chemin de code precis (voir
            # `_probe_tetgen_conflicts`). Un sondage DEDIE en mode
            # detection, plus leger, recupere le texte complet de facon
            # fiable.
            tetgen_output = _probe_tetgen_conflicts(current)
            conflict_vertices = _parse_tetgen_conflict_vertices(tetgen_output)
            if not conflict_vertices:
                return current, None
            faces = np.asarray(current.faces)
            seed_faces = [
                i for i, face in enumerate(faces)
                if len(conflict_vertices.intersection(int(v) for v in face)) >= 2
            ]
            if not seed_faces:
                return current, None
            repaired = _repair_local_self_intersections(
                current, extra_seed_faces=seed_faces, require_euler_preserved=False,
            )
            if repaired is None:
                return current, None
            current = repaired
    return current, None


def _tetrahedralize(
    mesh: trimesh.Trimesh,
    resolution: float,
    near_platform_resolution: Optional[float] = None,
    near_platform_margin: float = 5.0,
    platform_z: float = 0.0,
    junction_zones: Optional[list[JunctionZone]] = None,
    junction_refinement_resolution: Optional[float] = None,
) -> TetMesh:
    """
    Genere un maillage tetraedrique volumique respectant la frontiere du mesh de surface, via TetGen.
    `resolution` donne le volume maximal par tetraedre (`-a`), plafond global.

    `near_platform_resolution` (opt-in) raffine LOCALEMENT pres du plateau (z < platform_z +
    near_platform_margin) via deux regions TetGen (`_split_for_platform_refinement` + `add_region`) ; le
    mecanisme `bgmesh`/`metric` est ecarte (ignore par tetgen 0.8.4 avec `switches=`, segfault avec `-m`,
    limite #35). Corrige l'artefact d'« escalier » des couches pres du plateau pour ~5.4x moins de tetraedres
    qu'un raffinement global. Limite #42 : `_split_for_platform_refinement` garantit un seed valide pour les
    deux regions (`_find_guaranteed_interior_point`).

    `junction_zones`/`junction_refinement_resolution` (opt-in) raffinent les jonctions multi-region a forte
    divergence directionnelle (limite #47/#48) : un seed interne supplementaire, PLUSIEURS zones FUSIONNEES en
    un seul seed (`_merged_junction_seed`), jamais une par zone (le bug #42 se generalise a N regions).

    Cascade reactive (TetGen est le seul juge fiable des auto-intersections ; un maillage accepte n'est jamais
    modifie) : tel quel, apres resolution geometrique, apres pymeshfix, puis sans raffinement plateau.
    """
    if not mesh.is_watertight:
        raise MeshIngestError(
            "Impossible de tetraedriser un mesh non etanche apres reparation. "
            "Verifiez la geometrie source."
        )
    opts = _TetOptions(resolution, near_platform_margin, platform_z, junction_zones, junction_refinement_resolution)
    attempts: list[tuple[str, trimesh.Trimesh, Optional[float]]] = [("tel quel", mesh, near_platform_resolution)]
    speculative_attempts: dict = {}
    probe_results: dict = {}
    resolved, resolved_collapsed = _geometric_resolution_attempt(
        mesh, opts, near_platform_resolution, speculative_attempts, probe_results)
    if resolved is not None:
        attempts.append(("apres resolution geometrique", resolved, near_platform_resolution))
    repaired, pymeshfix_collapsed = _pymeshfix_attempt(mesh)
    if repaired is not None:
        attempts.append(("apres reparation pymeshfix", repaired, near_platform_resolution))
    # Dernier recours : abandonner le raffinement plateau (faces de bouchage coplanaires que TetGen peut
    # refuser). Preference : `resolved` > `repaired` > `mesh` (sur turkish.stl, retomber sur `mesh` brut
    # donnait un echec garanti pour une raison sans rapport avec celle deja reglee).
    if near_platform_resolution is not None:
        source = resolved if resolved is not None else (repaired if repaired is not None else mesh)
        attempts.append(("sans raffinement plateau", source, None))

    tet_mesh, last_error = _run_tetgen_attempts(attempts, opts, speculative_attempts, probe_results)
    if tet_mesh is not None:
        return tet_mesh
    hint = _tetgen_failure_hint(repaired, pymeshfix_collapsed, resolved_collapsed)
    raise MeshIngestError(
        f"Le maillage contient des auto-intersections (ou un autre defaut) que "
        f"TetGen refuse ou sur lequel il plante, et que la reparation "
        f"automatique n'a pas suffi a corriger : {last_error}.{hint}"
    )


# Repli par TRANSLATION (2026-10-06). La cascade de reparation des auto-intersections
# depend de l'arrondi flottant des coordonnees ABSOLUES : turkish.stl passe a sa
# position d'origine mais echoue decale de 20, 40 ou 70 mm en x (resolution
# geometrique laissant 20 auto-intersections residuelles au lieu de 8, que la
# reparation locale ne sait plus traiter). Rejouer la cascade sur le maillage
# translate EN MEMOIRE (centre de la boite en x, y = 0) puis ramener les noeuds
# a leur place fait passer les trois cas. Uniquement apres un echec : un maillage
# deja accepte n'est jamais modifie. Decalages essayes dans l'ordre (mm, en x, y).
_TRANSLATION_RETRY_OFFSETS = ((0.0, 0.0), (0.37, 0.21))


def _tetrahedralize_with_translation_retry(
    mesh, resolution, near_platform_resolution, near_platform_margin, platform_z,
    junction_zones, junction_refinement_resolution,
) -> TetMesh:
    try:
        return _tetrahedralize(
            mesh, resolution, near_platform_resolution, near_platform_margin, platform_z,
            junction_zones, junction_refinement_resolution,
        )
    except MeshIngestError as first_error:
        if junction_zones:      # zones de jonction en coordonnees absolues : pas de repli
            raise
        center = (np.asarray(mesh.bounds[0]) + np.asarray(mesh.bounds[1])) / 2.0
        for dx, dy in _TRANSLATION_RETRY_OFFSETS:
            offset = np.array([dx - center[0], dy - center[1], 0.0])
            moved = mesh.copy()
            moved.apply_translation(offset)
            try:
                tm = _tetrahedralize(
                    moved, resolution, near_platform_resolution, near_platform_margin, platform_z,
                    None, junction_refinement_resolution,
                )
            except MeshIngestError:
                continue
            print(f"[step1] tetraedrisation reussie apres translation de {offset[:2].round(3).tolist()} mm "
                  f"(repli d'arrondi flottant), noeuds ramenes a leur place")
            return TetMesh(tet_vertices=np.asarray(tm.tet_vertices, dtype=float) - offset,
                           tetrahedra=tm.tetrahedra, tet_adjacency=tm.tet_adjacency)
        raise first_error


def _tet_mesh_unusable(tm: TetMesh) -> Optional[str]:
    """Raison pour laquelle un maillage tetraedrique est inexploitable (valences
    aberrantes : millions de tetraedres de 1e-18 mm3 autour de details quasi coincidents,
    Fighter jet F18), ou None. Meme critere que l autodocteur."""
    from element_ardor.mesh_doctor import SurfaceTet, tet_mesh_problem

    return tet_mesh_problem(SurfaceTet(points=np.asarray(tm.tet_vertices), cells=np.asarray(tm.tetrahedra)))


def _repair_and_tetrahedralize(
    mesh: trimesh.Trimesh, first_error: Exception, resolution: float, near_platform_resolution,
    near_platform_margin, platform_z, junction_zones, junction_refinement_resolution,
) -> tuple[trimesh.Trimesh, TetMesh, str]:
    """REPARATION SYSTEMATIQUE de la surface (2026-10-22), appelee quand l'autodocteur trouve un
    defaut AVANT TetGen, ou quand TetGen refuse / produit un maillage inexploitable. Cascade
    la reparation d'Apex (`element_ardor.mesh_doctor.repair_candidates` : soudure des aretes courtes,
    union des coques qui se traversent, resolution geometrique, pymeshfix), chaque candidat passe
    par la tetraedrisation habituelle, rejete s'il ampute la piece (etendue) ou donne un
    maillage inexploitable ; en dernier recours TetGen sans point de Steiner en surface (`Y`).
    Rend (surface retenue, maillage tetraedrique, libelle) ; leve `MeshIngestError` avec TOUTES
    les raisons si rien ne passe - jamais un maillage faux en silence."""
    from element_ardor.mesh_doctor import repair_candidates, tetrahedralize as doctor_tetrahedralize

    cleaned, labelled = repair_candidates(mesh)
    reference = np.asarray(mesh.extents, dtype=float)
    errors = [f"premier refus : {str(first_error)[-300:]}"]
    usable: list[tuple[str, trimesh.Trimesh]] = []
    for label, candidate in labelled:
        if np.any(np.asarray(candidate.extents) < 0.99 * reference):
            errors.append(f"{label} : rejetee (etendue {np.round(candidate.extents, 1)} contre {np.round(reference, 1)})")
            continue
        usable.append((label, candidate))
        try:
            tm = _tetrahedralize_with_translation_retry(
                candidate, resolution, near_platform_resolution, near_platform_margin, platform_z,
                junction_zones, junction_refinement_resolution,
            )
        except (MeshIngestError, RuntimeError) as exc:
            errors.append(f"{label} : TetGen refuse ({str(exc)[-200:]})")
            continue
        problem = _tet_mesh_unusable(tm)
        if problem is not None:
            errors.append(f"{label} : maillage tetraedrique inexploitable ({problem})")
            continue
        return candidate, tm, label
    for label, candidate in usable:
        try:
            doctor_mesh = doctor_tetrahedralize(np.asarray(candidate.vertices), np.asarray(candidate.faces),
                                       resolution ** 3, keep_boundary=True)
        except RuntimeError as exc:
            errors.append(f"{label} + Y : TetGen refuse ({str(exc)[-200:]})")
            continue
        tm = TetMesh(tet_vertices=np.asarray(doctor_mesh.points, dtype=float), tetrahedra=np.asarray(doctor_mesh.cells),
                     tet_adjacency=_build_tet_adjacency(np.asarray(doctor_mesh.cells), len(doctor_mesh.points)))
        problem = _tet_mesh_unusable(tm)
        if problem is not None:
            errors.append(f"{label} + Y : maillage tetraedrique inexploitable ({problem})")
            continue
        return candidate, tm, f"{label} + TetGen sans point de Steiner en surface"
    raise MeshIngestError("surface irreparable : " + " | ".join(errors))


def run(
    mesh_path: str,
    resolution: float,
    options: Optional[IngestOptions] = None,
    **overrides,
) -> tuple[MeshData, TetMesh, IngestReport]:
    """
    Charge, nettoie et prepare le mesh source pour le reste du pipeline.

    `resolution` : taille cible de tetraedre. `options` (`IngestOptions`, `overrides` = champs nommes) :
    `allow_multi_mesh` fusionne les geometries d'une Scene au lieu de lever une erreur ;
    `near_platform_resolution`/`near_platform_margin` (opt-in) raffinement pres du plateau ; `platform_z` :
    `None` = detection automatique (`min(z)` du mesh nettoye), sinon override total (valeur utilisee dans
    `IngestReport.platform_z_used`) ; `junction_zones` / `junction_refinement_resolution` (opt-in) zones de
    jonction a raffiner (limite #47/#48) ; `check_surface_coverage` (opt-in) part de la surface hors du volume
    tetraedrise ; `place_on_platform` (opt-in) translation en z SEULE (decalage dans
    `IngestReport.placement_offset_z`).

    Returns: (MeshData, TetMesh, IngestReport).
    Raises: MeshIngestError si le fichier est illisible, vide, ou non tetraedrisable meme apres reparation.
    """
    opts = dataclasses.replace(options or IngestOptions(), **overrides)
    allow_multi_mesh, near_platform_resolution = opts.allow_multi_mesh, opts.near_platform_resolution
    near_platform_margin, platform_z, junction_zones = opts.near_platform_margin, opts.platform_z, opts.junction_zones
    junction_refinement_resolution = opts.junction_refinement_resolution
    check_surface_coverage, place_on_platform = opts.check_surface_coverage, opts.place_on_platform
    mesher, tetranet_refuse_ratio = opts.mesher, opts.tetranet_refuse_ratio
    raw_mesh = _load_raw_mesh(mesh_path, allow_multi_mesh)
    _progress.milestone(0.5)                 # R-017 : fichier lu (etape « import »)
    if mesher == "tetranet":
        if junction_zones:
            print("[step1] TetraNet : zones de jonction non prises en charge, retour a TetGen")
        elif not tetranet_available():
            print("[step1] TetraNet : noyau natif absent, retour a TetGen")
        else:
            if near_platform_resolution is not None:
                print("[step1] TetraNet : raffinement pres du plateau pas encore porte, ignore")
            _progress.milestone(0.9)
            _progress.switch_stage("tetrahedralization")   # TetraNet nettoie et maille en un appel (R-022)
            return _ingest_with_tetranet(raw_mesh, resolution, platform_z, place_on_platform,
                                         check_surface_coverage, tetranet_refuse_ratio,
                                         opts.tetranet_feature_angle_deg)
    cleaned_mesh, report = _clean_and_repair(raw_mesh)
    _progress.milestone(0.9)                 # surface nettoyee et reparee : fin de l'« import » (R-022)
    _progress.switch_stage("tetrahedralization")   # reste le maillage volumique (TetGen, indivisible)

    if place_on_platform and len(cleaned_mesh.vertices):
        _place_mesh_on_platform(cleaned_mesh, report, platform_z)

    resolved_platform_z = (
        float(cleaned_mesh.vertices[:, 2].min()) if platform_z is None else platform_z
    )
    report.platform_z_used = resolved_platform_z
    cleaned_mesh, tet_mesh = _tetrahedralize_with_autodoctor(
        cleaned_mesh, report, resolution, near_platform_resolution, near_platform_margin, resolved_platform_z,
        junction_zones, junction_refinement_resolution,
    )
    _progress.milestone(0.9)                 # volume tetraedrise
    mesh_data = _build_mesh_data(cleaned_mesh)

    report.tetrahedra_generated = tet_mesh.tet_count
    if check_surface_coverage:
        _surface_coverage_into_report(report, cleaned_mesh, tet_mesh)
    return mesh_data, tet_mesh, report


def _geometric_resolution_attempt(
    mesh: trimesh.Trimesh, opts: _TetOptions, near_platform_resolution: Optional[float],
    speculative_attempts: dict, probe_results: dict,
):
    """Resolution GEOMETRIQUE des auto-intersections, tentee EN PREMIER (avant pymeshfix) : elle preserve la
    topologie et l'aire de la geometrie non affectee. Best-effort : ne leve jamais. Retourne
    `(resolved, collapsed)` ; `resolved` est None si le defaut n'est pas de cette nature ou si le resultat
    s'effondre (moins de la moitie des faces : garde-fou limite #51/#52). Alimente `speculative_attempts` et
    `probe_results` (indexes par `id(resolved)`)."""
    resolved = resolve_self_intersections(mesh)
    if resolved is not None and len(resolved.faces) < 0.5 * len(mesh.faces):
        return None, True
    if resolved is None:
        return None, False
    # Auto-intersections RESIDUELLES introduites par le soudage de trous minuscules : nettoyees sans la garde
    # d'Euler stricte (le soudage a deja change l'euler). turkish.stl : 8 -> 0, ecart de volume 2.9e-6 %.
    if _find_self_intersecting_faces(resolved):
        cleaned_resolved = _repair_local_self_intersections(resolved, require_euler_preserved=False)
        if cleaned_resolved is not None:
            resolved = cleaned_resolved
    # Execution SPECULATIVE (2026-10-01) : la sonde de `_finish_geometric_resolution` est une vraie
    # tetraedrisation SANS raffinement, jetee ; la tentative suivante (AVEC raffinement) est un autre calcul.
    # Quand la sonde reussit, elle rend `resolved` INTACT : la tentative finale est lancee EN PARALLELE et son
    # issue n'est gardee que si la sonde n'a rien modifie (resultat identique, temps = max au lieu de somme).
    speculative = None
    if runtime.options().speculative_tetgen and not opts.junction_zones:
        import concurrent.futures as _cf
        _executor = _cf.ThreadPoolExecutor(max_workers=1)
        speculative = _executor.submit(
            runtime.bind(_tetrahedralize_once_isolated), resolved, opts.resolution, near_platform_resolution,
            opts.near_platform_margin, opts.platform_z, opts.junction_zones, opts.junction_refinement_resolution,
        )
        _executor.shutdown(wait=False)
    before_finish = resolved
    resolved, probe_tet_mesh = _finish_geometric_resolution(resolved, opts.resolution)
    if probe_tet_mesh is not None and not opts.junction_zones:
        probe_results[id(resolved)] = probe_tet_mesh
    if speculative is not None and resolved is before_finish:
        speculative_attempts[id(resolved)] = speculative
    return resolved, False


def _pymeshfix_attempt(mesh: trimesh.Trimesh):
    """Reparation pymeshfix avec le garde-fou anti-effondrement de `_clean_and_repair` (limite #51) : MeshFix
    peut rendre un mesh watertight en supprimant une portion entiere de la geometrie (turkish.stl : 10104 ->
    1488 faces, z tronque), que TetGen tetraedrise sans erreur. Retourne `(repaired, collapsed)`."""
    repaired = _repair_with_pymeshfix(mesh)
    if repaired is not None and len(repaired.faces) < 0.5 * len(mesh.faces):
        return None, True
    return repaired, False


def _run_tetgen_attempts(attempts: list, opts: _TetOptions, speculative_attempts: dict, probe_results: dict):
    """Essaie chaque tentative dans l'ordre. Retourne `(tet_mesh, last_error)` ; `tet_mesh` est None si toutes
    ont echoue sur une erreur recuperable (crash natif ou auto-intersection)."""
    last_error: Optional[Exception] = None
    for _label, candidate, npr in attempts:
        try:
            future = speculative_attempts.get(id(candidate)) if npr == attempts[0][2] else None
            if future is not None:
                return future.result(), None
            if npr is None and id(candidate) in probe_results:
                return probe_results[id(candidate)], None
            return _tetrahedralize_once_isolated(
                candidate, opts.resolution, npr, opts.near_platform_margin, opts.platform_z,
                opts.junction_zones, opts.junction_refinement_resolution,
            ), None
        except _TetgenCrashError as exc:
            # Crash NATIF dans le sous-processus isole : meme signal d'echec qu'une auto-intersection, le repli
            # suivant a la meme chance de reussir. DOIT rester AVANT `except RuntimeError` (sous-classe).
            last_error = exc
        except RuntimeError as exc:
            if not _is_recoverable_tetgen_error(exc):
                raise
            last_error = exc
    return None, last_error


def _tetgen_failure_hint(repaired, pymeshfix_collapsed: bool, resolved_collapsed: bool) -> str:
    if repaired is not None:
        hint = ""
    elif pymeshfix_collapsed:
        hint = (
            " La reparation `pymeshfix` a produit un resultat watertight "
            "mais radicalement plus petit que le mesh source (effondrement "
            "geometrique, rejete par garde-fou) - le defaut source n'a pas "
            "de reparation automatique fiable sur ce mesh."
        )
    else:
        hint = (
            " La reparation automatique n'a pas pu s'appliquer (dependance "
            "optionnelle `pymeshfix` absente ou en echec) - l'installer peut suffire."
        )
    if resolved_collapsed:
        hint += (
            " La resolution geometrique des auto-intersections a produit un "
            "resultat radicalement plus petit que le mesh source (rejete par "
            "garde-fou)."
        )
    return hint


@dataclass(frozen=True)
class IngestOptions:
    """Options de `run` (memes noms et defauts que les anciens parametres nommes)."""

    allow_multi_mesh: bool = False
    near_platform_resolution: Optional[float] = None
    near_platform_margin: float = 5.0
    platform_z: Optional[float] = None
    junction_zones: Optional[list[JunctionZone]] = None
    junction_refinement_resolution: Optional[float] = None
    check_surface_coverage: bool = False
    place_on_platform: bool = False
    mesher: str = "tetgen"
    tetranet_refuse_ratio: float = 0.05
    tetranet_feature_angle_deg: float = 0.0     # > 0 : arêtes vives protégées par TetraNet (R-030 point 8, piste C)


def _surface_coverage_into_report(report: IngestReport, mesh: trimesh.Trimesh, tet_mesh: TetMesh) -> None:
    from core.tet_coverage import surface_outside_volume

    coverage = surface_outside_volume(mesh.vertices, mesh.faces, tet_mesh.tet_vertices, tet_mesh.tetrahedra)
    report.surface_outside_volume_fraction = coverage.far_fraction
    report.surface_outside_volume_area = coverage.far_area
    report.surface_max_distance_to_tets = coverage.max_distance
    report.tet_volume_ratio = coverage.volume_ratio
    report.surface_outside_volume_clusters = [(c.area, c.center) for c in coverage.clusters[:10]]


def _tetrahedralize_with_autodoctor(
    cleaned_mesh: trimesh.Trimesh, report: IngestReport, resolution: float,
    near_platform_resolution: Optional[float], near_platform_margin: float, resolved_platform_z: float,
    junction_zones: Optional[list[JunctionZone]], junction_refinement_resolution: Optional[float],
):
    """Autodocteur : la surface est TOUJOURS analysee AVANT TetGen. Reparation immediate quand des coques se
    traversent, que la surface n'est pas etanche ou que son orientation est incoherente (TetGen peut alors
    ACCEPTER la surface et produire des tetraedres degeneres : Fighter jet F18) ; reparation aussi apres un
    refus de TetGen ou un maillage inexploitable. Une surface saine n'est jamais modifiee. Une
    auto-intersection dans UNE coque garde la cascade historique de `_tetrahedralize` (limite #43).
    Retourne `(cleaned_mesh, tet_mesh)` (le maillage de surface peut avoir ete repare)."""
    from element_ardor.mesh_doctor import diagnose_surface

    diagnosis = diagnose_surface(cleaned_mesh)
    report.surface_problems = list(diagnosis.problems)
    pre_repair = (diagnosis.shells > 1 and diagnosis.self_intersecting_faces > 0) or not diagnosis.watertight
    tet_mesh = None
    first_error: Optional[Exception] = None
    if pre_repair:
        first_error = MeshIngestError("defauts de surface detectes avant TetGen : " + ", ".join(diagnosis.problems))
    else:
        try:
            tet_mesh = _tetrahedralize_with_translation_retry(
                cleaned_mesh, resolution, near_platform_resolution, near_platform_margin, resolved_platform_z,
                junction_zones, junction_refinement_resolution,
            )
            problem = _tet_mesh_unusable(tet_mesh)
            if problem is not None:
                first_error, tet_mesh = MeshIngestError(f"maillage tetraedrique inexploitable ({problem})"), None
        except (MeshIngestError, RuntimeError) as exc:
            if isinstance(exc, ModelTooLargeError):
                raise
            first_error = exc
    if tet_mesh is None:
        print(f"[step1] reparation de la surface : {first_error}")
        final_mesh, tet_mesh, label = _repair_and_tetrahedralize(
            cleaned_mesh, first_error, resolution, near_platform_resolution, near_platform_margin,
            resolved_platform_z, junction_zones, junction_refinement_resolution,
        )
        report.surface_repair = label
        report.is_watertight_after_repair = bool(final_mesh.is_watertight)
        print(f"[step1] surface reparee : {label} ({tet_mesh.tet_count} tetraedres)")
        cleaned_mesh = final_mesh
    return cleaned_mesh, tet_mesh


def _place_mesh_on_platform(mesh: trimesh.Trimesh, report: IngestReport, platform_z: Optional[float]) -> None:
    """Translation en z SEULE : min(z) -> `platform_z` (0.0 si None)."""
    target_z = 0.0 if platform_z is None else float(platform_z)
    offset_z = target_z - float(mesh.vertices[:, 2].min())
    if offset_z != 0.0:
        mesh.apply_translation([0.0, 0.0, offset_z])
    report.placement_offset_z = offset_z

from element_ardor.step1_repair import (  # noqa: E402,F401  (reexport : imports et patchs existants)
    _boundary_successors,
    _boundary_cycles,
    _split_cycles_into_simple_loops,
    _LoopTriangleWeights,
    _fill_all_loops,
    _validated_repair,
    _grow_removed_around,
)
from element_ardor.step1_tetgen_run import (  # noqa: E402,F401  (reexport : imports et patchs existants)
    _run_worker_process,
    _tetgen_with_platform_refinement,
    _tetgen_with_junction_refinement,
    _check_seeds_outside_junctions,
)
