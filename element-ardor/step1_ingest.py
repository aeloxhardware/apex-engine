"""
pipeline/step1_ingest.py

Etape 1 du pipeline : ingestion du mesh source (STL/3MF) et preparation de
toutes les structures dont les etapes suivantes ont besoin.

Cette etape est la seule du pipeline qui touche a l'I/O disque cote geometrie.
Elle ne depend d'aucune autre etape - c'est le point d'entree du pipeline et
le fichier le plus simple a tester en isolation (aucune donnee amont requise,
juste un chemin de fichier).

Sortie : (MeshData, TetMesh, IngestReport)
"""

from __future__ import annotations

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
from core import _cpp
from core.mesh_types import JunctionZone, MeshData, TetMesh
from core.spatial import KDTree, RayIntersector
from pipeline._self_intersection_resolve import resolve_self_intersections

from pipeline.step1_types import (  # noqa: F401  (réexport : imports existants)
    IngestReport,
    MeshIngestError,
    ModelTooLargeError,
    SurfaceInconsistentError,
    _load_raw_mesh,
    _min_enclosing_circle,
    check_fits_platform,
)
from pipeline.step1_repair import (  # noqa: F401  (réexport : imports existants)
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
from pipeline.step1_tetgen_seeds import (  # noqa: F401  (réexport : imports existants)
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
from pipeline.step1_tetgen_run import (  # noqa: F401  (réexport : imports existants)
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
)
from pipeline.step1_build import (  # noqa: F401  (réexport : imports existants)
    TETRANET_CELL_PER_RESOLUTION,
    TETRANET_STATUS_INCONSISTENT_SURFACE,
    _build_acceleration_structures,
    _build_mesh_data,
    _ingest_with_tetranet,
    tetranet_available,
)


# Tentative finale de TetGen lancee en parallele de la sonde de
# `_finish_geometric_resolution` (voir `_tetrahedralize`). Interrupteur de
# module pour les tests d'equivalence (resultat identique dans les deux cas).
_SPECULATIVE_TETGEN = True


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
    Affine iterativement le resultat de `resolve_self_intersections` en
    utilisant TetGen LUI-MEME comme oracle (principe deja acte, limite
    #43 : "TetGen est le seul juge fiable de ce qui le bloque").

    A chaque iteration : tente une VRAIE tetraedrisation (subprocess isole,
    `_tetrahedralize_once_isolated`) ; si TetGen la refuse pour une raison
    RECUPERABLE et nomme des sommets en conflit dans sa sortie textuelle
    (via `_probe_tetgen_conflicts`), cible une reparation locale
    (`_repair_local_self_intersections`, `extra_seed_faces` = les faces
    portant ces sommets) sur EXACTEMENT ce que TetGen vient de refuser, au
    lieu d'un repli generique. Mesure sur turkish.stl : 8 auto-
    intersections locales (issues du soudage de trous residuels de
    `resolve_self_intersections`) rendaient TetGen incapable d'aboutir
    malgre un SAT=0 global (contact EXACT entre deux segments de faces
    partageant un sommet, donc invisible au test SAT qui exclut ces
    paires) ; ce mecanisme reduit le residu jusqu'a un dernier segment
    nomme, cible et repare en 1 iteration supplementaire - TetGen aboutit
    alors sur la geometrie COMPLETE (330306 tetraedres, z=[30.26,90.78],
    contre l'echec total precedent), pour un ecart de volume negligeable.

    Deliberement SANS raffinement plateau/jonction ici (jamais
    `near_platform_resolution`/`junction_refinement_resolution`), meme si
    l'appelant final en utilise un : ce raffinement construit un PLC
    COMBINE avec de nouveaux sommets injectes au plan de coupe
    (`_split_for_platform_refinement`), dont la numerotation ne correspond
    PAS a celle de `current` - un conflit nomme par TetGen sur CE PLC
    combine ne peut donc pas etre traduit fidelement en `extra_seed_faces`
    sur le mesh de surface qu'on repare ici. Bug mesure en integrant ce
    mecanisme : avec raffinement plateau actif, le sondage ne trouvait
    jamais aucun sommet exploitable (mauvais referentiel), la fonction
    rendait la main sans le moindre progres, et la tentative reelle de
    l'appelant (elle, avec raffinement) echouait alors identiquement a
    chaque fois. La question posee ici est volontairement plus etroite :
    « ce mesh de SURFACE est-il tetraedrisable EN SOI » - la compatibilite
    avec le raffinement plateau reste de la responsabilite de la cascade
    normale de `_tetrahedralize`, qui sait deja s'en degrader proprement
    (repli "sans raffinement plateau", limite #44).

    Best-effort, ne leve jamais : rend le meilleur mesh obtenu (le dernier
    `resolved` valide) des qu'un probe reussit, qu'une reparation echoue,
    qu'aucun sommet n'est nommable, ou que le budget d'iterations est
    epuise - l'appelant (`_tetrahedralize`) refera de toute facon une VRAIE
    tentative sur ce resultat dans sa cascade normale, avec son cortege
    habituel de replis (pymeshfix, sans raffinement plateau).
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
    Genere un maillage tetraedrique volumique respectant la frontiere du mesh
    de surface, via TetGen. `resolution` est converti en volume maximal
    autorise par tetraedre (parametre -a de tetgen), qui reste le plafond
    global.

    `near_platform_resolution` (opt-in, None = desactive par defaut - meme
    discipline que `direction_blend_radius=0.0` en step4) raffine LOCALEMENT
    la densite pres du plateau (z < platform_z + near_platform_margin), via
    deux regions TetGen (`_split_for_platform_refinement` + `add_region`,
    une pour la tranche pres du plateau a `near_platform_resolution`, une
    pour le reste a `resolution`) - PAS via un maillage de fond
    (`bgmesh`/`metric`), mecanisme teste et ecarte : silencieusement ignore
    par le binding tetgen installe (0.8.4) des qu'un `switches=` non vide
    est passe, et provoque un segfault natif reproductible si on force `-m`
    malgre tout (teste sous plusieurs formes, voir CLAUDE.md). Corrige une
    erreur de discretisation FEM mesuree pres du plateau (CLAUDE.md) : un
    tetraedre P1 assez grand pour s'etendre du contact plateau jusqu'a une
    zone ou la direction cible a deja beaucoup tourne ne peut representer
    qu'un seul gradient constant sur toute son etendue, produisant un
    artefact d'« escalier » visible dans l'inclinaison des couches - affine
    localement plutot que globalement pour ne pas payer le cout d'un
    maillage fin partout ou ce n'est pas necessaire (mesure sur
    l_shape_overhang.stl : amelioration angulaire quasi identique a un
    raffinement global complet, pour ~5.4x moins de tetraedres).

    Bug corrige (CLAUDE.md limite #42) : `add_region` ignorait
    silencieusement le volume de la region "rest" (donc `resolution`
    perdait tout effet hors de la tranche pres du plateau) quand le seed de
    cette region tombait hors de sa propre coque sur une geometrie non
    convexe - `_split_for_platform_refinement` garantit desormais un seed
    valide pour les deux regions (`_find_guaranteed_interior_point`).

    `junction_zones`/`junction_refinement_resolution` (opt-in, None =
    desactive par defaut - meme discipline) raffinent LOCALEMENT les
    jonctions multi-region a forte divergence directionnelle (CLAUDE.md
    limite #47/#48), generalisant le mecanisme ci-dessus au-dela du seul
    cas plateau. Contrairement au plateau, une jonction n'a pas de plan de
    coupe unique : pas de decoupe de surface, juste un seed interne
    supplementaire (voir `_merged_junction_seed`) - PLUSIEURS zones
    detectees sont FUSIONNEES en un seul seed avant d'etre passees a
    `add_region`, jamais une par zone : mesure directe
    (scripts/measure_add_region_n_zones.py) que le bug limite #42 se
    generalise a N regions actives simultanement (aucune des zones fines
    n'atteignait le facteur de densite attendu, une meme MOINS dense que
    la region censee etre grossiere) - le fallback fusionne fonctionne
    correctement (scripts/measure_add_region_merged_zone.py, facteur ~48x).
    """
    if not mesh.is_watertight:
        raise MeshIngestError(
            "Impossible de tetraedriser un mesh non etanche apres reparation. "
            "Verifiez la geometrie source."
        )

    # Strategie de repli, du moins invasif au plus invasif. TetGen est le
    # seul juge fiable des auto-intersections qui le bloquent : ce sont
    # souvent des contacts EXACTS ("Two segments exactly intersect", "A
    # vertex lies on a facet") que ni trimesh (watertight/winding
    # coherents), ni manifold3d (rend le maillage a l'identique), ni un test
    # geometrique de traversee ne signalent - mesure sur
    # res/test_models/death.stl : 3 triangles fautifs sur 12662, zero
    # detecte par le test SAT de _count_self_intersecting_faces. D'ou une
    # reparation REACTIVE (on repare parce que TetGen a refuse), jamais
    # preventive : un maillage deja accepte n'est jamais modifie.
    attempts: list[tuple[str, trimesh.Trimesh, Optional[float]]] = [
        ("tel quel", mesh, near_platform_resolution),
    ]

    # Resolution GEOMETRIQUE des auto-intersections (decoupe le long de la
    # courbe d'intersection + extraction de l'enveloppe exterieure par
    # nombre d'enroulement) - tentee EN PREMIER, avant pymeshfix :
    # contrairement a pymeshfix (qui reecrit tout le maillage), elle
    # preserve la topologie et l'aire de la geometrie non affectee, donc
    # strictement moins invasive. Best-effort : ne leve jamais, retourne
    # None si le defaut n'est pas de cette nature ou si le resultat ne
    # satisfait pas ses propres post-conditions (voir docstring du module).
    # Meme garde-fou anti-effondrement que pour pymeshfix ci-dessous (limite
    # #51/#52) : un resultat radicalement plus petit que l'entree est rejete
    # plutot qu'accepte silencieusement. Mesure sur turkish.stl (session
    # ayant introduit ce mecanisme) : reduit les auto-intersections de 50 a
    # 0 faces SAT et de 51 a 2 segments nommes par TetGen, mais ne suffit
    # pas encore a lui seul a faire accepter ce modele precis (recouvrements
    # coplanaires non traites) - reste neanmoins utile en cascade, avant le
    # repli plus destructeur pymeshfix.
    speculative_attempts: dict = {}
    probe_results: dict = {}
    resolved = resolve_self_intersections(mesh)
    resolved_collapsed = False
    if resolved is not None and len(resolved.faces) < 0.5 * len(mesh.faces):
        resolved = None
        resolved_collapsed = True
    if resolved is not None:
        # `resolve_self_intersections` peut rendre un mesh etanche et
        # connexe mais avec quelques auto-intersections RESIDUELLES,
        # introduites par son propre soudage de trous minuscules
        # (`_close_small_boundary_loops` - fusionner des sommets proches
        # peut faire chevaucher deux triangles fins voisins). Ce n'est plus
        # le defaut d'origine (recouvrement CSG transversal), c'est un
        # artefact local du soudage - nettoye ici avec le meme mecanisme
        # que ci-dessous (retrait + retriangulation Klincsek), sans la
        # garde de caracteristique d'Euler stricte (inadaptee : le soudage
        # de trous a deja legitimement change l'euler avant cet appel).
        # Mesure sur turkish.stl : 8 auto-intersections residuelles -> 0,
        # corps unique preserve, ecart de volume 2.9e-6%.
        residual_sat = _find_self_intersecting_faces(resolved)
        if residual_sat:
            cleaned_resolved = _repair_local_self_intersections(
                resolved, require_euler_preserved=False,
            )
            if cleaned_resolved is not None:
                resolved = cleaned_resolved
        # Le test SAT ci-dessus ignore volontairement toute paire de faces
        # partageant un sommet (contact par construction) - or un contact
        # EXACT entre deux faces qui partagent DEJA un sommet ailleurs
        # (ex: un eventail de triangles fins autour d'un sommet issu du
        # soudage) lui est invisible, alors que TetGen le refuse.
        # `_finish_geometric_resolution` boucle sur TetGen lui-meme comme
        # oracle final pour ce residu : mesure sur turkish.stl, un dernier
        # segment nomme (2 faces partageant le sommet 890) restait apres
        # le nettoyage SAT ci-dessus ; cible et repare en 1 iteration
        # supplementaire, TetGen aboutit alors sur la geometrie COMPLETE.
        # Execution SPECULATIVE (session du 2026-10-01) : la sonde de
        # `_finish_geometric_resolution` est une vraie tetraedrisation SANS
        # raffinement, jetee ; la tentative suivante (AVEC raffinement) est
        # un autre calcul. Mesure sur death.stl : 6,55 s + 6,85 s en serie.
        # Quand la sonde reussit (cas nominal), elle rend `resolved` INTACT :
        # la tentative finale porte alors sur ce meme objet. On la lance donc
        # EN PARALLELE de la sonde (deux sous-processus TetGen) et on ne
        # garde son issue que si la sonde n'a rien modifie - resultat
        # identique dans tous les cas, temps = max au lieu de somme.
        speculative = None
        if _SPECULATIVE_TETGEN and not junction_zones:
            import concurrent.futures as _cf
            _executor = _cf.ThreadPoolExecutor(max_workers=1)
            speculative = _executor.submit(
                _tetrahedralize_once_isolated, resolved, resolution, near_platform_resolution,
                near_platform_margin, platform_z, junction_zones, junction_refinement_resolution,
            )
            _executor.shutdown(wait=False)
        before_finish = resolved
        resolved, probe_tet_mesh = _finish_geometric_resolution(resolved, resolution)
        if probe_tet_mesh is not None and not junction_zones:
            probe_results[id(resolved)] = probe_tet_mesh
        if speculative is not None and resolved is before_finish:
            speculative_attempts[id(resolved)] = speculative
        attempts.append(("apres resolution geometrique", resolved, near_platform_resolution))

    # Meme garde-fou anti-effondrement que `_clean_and_repair` (limite #51) :
    # MeshFix peut rendre un mesh watertight en supprimant une portion
    # entiere de la geometrie plutot qu'en la reparant. Bug reel mesure sur
    # `turkish.stl` : ce deuxieme appel a `_repair_with_pymeshfix`, propre a
    # la cascade de tetraedrisation et jusqu'ici SANS garde-fou (seul
    # `_clean_and_repair` en avait un), acceptait un resultat passant de
    # 10104 a 1488 faces - watertight, donc silencieusement accepte par la
    # cascade, mais amputant tout le haut du modele (z tronque de 90.78 a
    # 63.26mm). TetGen tetraedrise alors ce mesh mutile SANS AUCUNE erreur
    # ni warning (il est parfaitement valide, juste faux) : le pipeline
    # entier (step2 a step8) traite ensuite fidelement un modele auquel il
    # manque les deux tiers de sa hauteur, sans jamais lever d'exception -
    # silence total jusqu'a l'inspection visuelle du GLB final.
    repaired = _repair_with_pymeshfix(mesh)
    pymeshfix_collapsed = False
    if repaired is not None and len(repaired.faces) < 0.5 * len(mesh.faces):
        # Effondrement detecte : ce resultat n'est utilisable NULLE PART
        # dans la cascade, ni comme tentative dediee ni comme source de
        # repli pour "sans raffinement plateau" ci-dessous (meme mesh
        # mutile, memes deux tiers de hauteur disparus).
        repaired = None
        pymeshfix_collapsed = True
    if repaired is not None:
        attempts.append(("apres reparation pymeshfix", repaired, near_platform_resolution))

    # Dernier recours : abandonner le raffinement local pres du plateau.
    # `_split_for_platform_refinement` recolle deux coques obtenues par
    # slice_plane(cap=True) ; leurs faces de bouchage sont coplanaires et
    # superposees le long de la frontiere de coupe, ce que TetGen peut
    # refuser sur une geometrie dense qui touche ce plan en de nombreux
    # points - mesure sur death.stl repare : 107037 tetraedres SANS
    # raffinement, refus net AVEC. Mieux vaut un maillage un peu plus
    # grossier pres du plateau qu'aucun maillage du tout.
    #
    # Ordre de preference : `resolved` (resolution geometrique, la moins
    # invasive) > `repaired` (pymeshfix) > `mesh` (brut). Bug reel mesure
    # en integrant ce repli avec `resolved` : la source ne considerait
    # QUE `repaired`/`mesh`, jamais `resolved` - sur turkish.stl, avec
    # `near_platform_resolution` actif, l'attempt "apres resolution
    # geometrique" echouait pour une raison PROPRE au raffinement plateau
    # ("split_segment", nouvelle variante limite #44 sur le PLC combine
    # avec ce mesh precis), mais ce repli retombait alors sur `mesh` NON
    # REPARE (51 auto-intersections brutes) au lieu du mesh deja nettoye
    # par la resolution geometrique - echec garanti pour une raison sans
    # aucun rapport avec celle qui vient d'etre reglee.
    if near_platform_resolution is not None:
        source = resolved if resolved is not None else (
            repaired if repaired is not None else mesh
        )
        attempts.append(("sans raffinement plateau", source, None))

    last_error: Optional[Exception] = None
    for _label, candidate, npr in attempts:
        try:
            future = speculative_attempts.get(id(candidate)) if npr == near_platform_resolution else None
            if future is not None:
                return future.result()
            if npr is None and id(candidate) in probe_results:
                return probe_results[id(candidate)]
            return _tetrahedralize_once_isolated(
                candidate, resolution, npr, near_platform_margin, platform_z,
                junction_zones, junction_refinement_resolution,
            )
        except _TetgenCrashError as exc:
            # Crash NATIF (segfault/access violation) dans le sous-processus
            # isole, jamais une RuntimeError Python normale - traite comme
            # le meme signal d'echec que "self-intersection" ci-dessous : le
            # repli suivant de la cascade (reparation pymeshfix, puis sans
            # raffinement plateau) a la meme chance de reussir, un crash sur
            # auto-intersection etant la cause la plus probable mesuree
            # (voir CLAUDE.md limite #43 et le worker isole plus bas).
            # DOIT rester AVANT `except RuntimeError` ci-dessous :
            # `_TetgenCrashError` en est une sous-classe, et Python retient
            # le PREMIER except applicable - inverser l'ordre ferait
            # tomber tout crash natif dans la branche "self-intersection"
            # (dont le message ne contient jamais ce texte), le relevant
            # immediatement au lieu de tenter le repli suivant. Bug reel
            # mesure sur `death.stl` modifie par l'utilisateur (2340
            # triangles) : le crash etait bien detecte (code de sortie
            # 3221225477 = 0xC0000005, access violation Windows) mais
            # remontait quand meme jusqu'a l'appelant au lieu de declencher
            # la reparation pymeshfix, a cause de cet ordre incorrect.
            last_error = exc
        except RuntimeError as exc:
            if not _is_recoverable_tetgen_error(exc):
                raise
            last_error = exc

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
    from pipeline.mesh_doctor import SurfaceTet, tet_mesh_problem

    return tet_mesh_problem(SurfaceTet(points=np.asarray(tm.tet_vertices), cells=np.asarray(tm.tetrahedra)))


def _repair_and_tetrahedralize(
    mesh: trimesh.Trimesh, first_error: Exception, resolution: float, near_platform_resolution,
    near_platform_margin, platform_z, junction_zones, junction_refinement_resolution,
) -> tuple[trimesh.Trimesh, TetMesh, str]:
    """REPARATION SYSTEMATIQUE de la surface (2026-10-22), appelee quand l'autodocteur trouve un
    defaut AVANT TetGen, ou quand TetGen refuse / produit un maillage inexploitable. Cascade
    la reparation d'Apex (`pipeline.mesh_doctor.repair_candidates` : soudure des aretes courtes,
    union des coques qui se traversent, resolution geometrique, pymeshfix), chaque candidat passe
    par la tetraedrisation habituelle, rejete s'il ampute la piece (etendue) ou donne un
    maillage inexploitable ; en dernier recours TetGen sans point de Steiner en surface (`Y`).
    Rend (surface retenue, maillage tetraedrique, libelle) ; leve `MeshIngestError` avec TOUTES
    les raisons si rien ne passe - jamais un maillage faux en silence."""
    from pipeline.mesh_doctor import repair_candidates, tetrahedralize as doctor_tetrahedralize

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
    allow_multi_mesh: bool = False,
    near_platform_resolution: Optional[float] = None,
    near_platform_margin: float = 5.0,
    platform_z: Optional[float] = None,
    junction_zones: Optional[list[JunctionZone]] = None,
    junction_refinement_resolution: Optional[float] = None,
    check_surface_coverage: bool = False,
    place_on_platform: bool = False,
    mesher: str = "tetgen",
    tetranet_refuse_ratio: float = 0.05,
) -> tuple[MeshData, TetMesh, IngestReport]:
    """
    Charge, nettoie et prepare le mesh source pour le reste du pipeline.

    Args:
        mesh_path: chemin vers le fichier STL/3MF source.
        resolution: taille cible de tetraedre (FieldSolverParams.resolution),
            utilisee pour controler la densite de la tetraedrisation.
        allow_multi_mesh: si True, fusionne les geometries multiples d'une
            Scene plutot que de lever une erreur.
        near_platform_resolution: taille de tetraedre cible pres du plateau,
            plus fine que `resolution` (opt-in, None = desactive par defaut,
            comportement historique de resolution globale unique inchange).
            Corrige une erreur de discretisation FEM mesuree pres du plateau
            (CLAUDE.md) sans payer le cout d'un maillage fin partout.
        near_platform_margin: hauteur (mm) au-dessus de `platform_z` sur
            laquelle le raffinement s'applique, avec une transition lisse
            vers `resolution` au-dela. Ignore si `near_platform_resolution`
            est None.
        platform_z: niveau du plateau. `None` (defaut) = DETECTION
            AUTOMATIQUE depuis le mesh nettoye (`min(z)` des sommets) -
            meme discipline de derivation geometrique que
            `_find_guaranteed_interior_point`/`_find_internal_cavities`
            plutot que d'exiger un parametre externe suppose correct. Une
            valeur explicite reste un override total. La valeur REELLEMENT
            utilisee (fournie ou deduite) est toujours remontee via
            `IngestReport.platform_z_used` - jamais un silence sur une
            valeur derivee, meme discipline que
            `self_intersecting_faces_before`. Corrige un bug reel mesure
            sur `res/test_models/death.stl` (min(z)=1.245, pas 0.0) : avec
            l'ancien defaut `platform_z=0.0` fige, `_compute_platform_contact`
            (step3) ne detectait AUCUNE face de contact plateau sur ce
            modele, cassant en cascade l'ancrage plateau du solveur FEM
            (step4) - voir CLAUDE.md.
        junction_zones: zones de jonction multi-region a raffiner localement
            (opt-in, None = desactive par defaut), typiquement produites par
            `pipeline.step1b_junction_detection.detect_junction_zones` sur
            une premiere passe legere (mesh + step2 + step3) orchestree par
            `main.run_pipeline` - `step1_ingest.run()` lui-meme ne les
            calcule jamais, il les recoit deja prealablement calculees.
            Corrige la non-monotonie du champ scalaire aux jonctions
            (CLAUDE.md limite #47/#48), portee laissee explicitement ouverte
            par le mecanisme `near_platform_resolution` ci-dessus (limite
            #36 : "traiter le cas jonction... non fait").
        junction_refinement_resolution: taille de tetraedre cible aux zones
            de jonction, plus fine que `resolution` (opt-in, None =
            desactive par defaut). Ignore si `junction_zones` est vide ou
            None.
        check_surface_coverage: mesure la part de la surface hors du volume
            tetraedrise et la remonte dans l'`IngestReport` (opt-in, defaut
            False : calcul supplementaire, aucun effet sur le maillage).
        place_on_platform: pose le modele sur le plateau par une translation
            en z SEULE (min(z) -> `platform_z`, ou 0.0 si `platform_z` est
            None) : la machine place son plateau a z = 0 et l'export refuse
            une pointe sous 0, or un STL peut avoir son minimum ailleurs
            (death.stl : +1.245 mm, real_overhang.stl : -15 mm). Opt-in
            (defaut False : coordonnees du fichier inchangees). Le decalage
            applique est remonte via `IngestReport.placement_offset_z`.

    Returns:
        Tuple (MeshData, TetMesh, IngestReport).

    Raises:
        MeshIngestError: si le fichier est illisible, vide, ou ne peut pas
            etre reduit a un volume tetraedrisable meme apres reparation.
    """
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
                                         check_surface_coverage, tetranet_refuse_ratio)
    cleaned_mesh, report = _clean_and_repair(raw_mesh)
    _progress.milestone(0.9)                 # surface nettoyee et reparee : fin de l'« import » (R-022)
    _progress.switch_stage("tetrahedralization")   # reste le maillage volumique (TetGen, indivisible)

    if place_on_platform and len(cleaned_mesh.vertices):
        target_z = 0.0 if platform_z is None else float(platform_z)
        offset_z = target_z - float(cleaned_mesh.vertices[:, 2].min())
        if offset_z != 0.0:
            cleaned_mesh.apply_translation([0.0, 0.0, offset_z])
        report.placement_offset_z = offset_z

    resolved_platform_z = (
        float(cleaned_mesh.vertices[:, 2].min()) if platform_z is None else platform_z
    )
    report.platform_z_used = resolved_platform_z

    # Autodocteur : la surface est TOUJOURS analysee AVANT TetGen. Reparation immediate quand
    # des coques se traversent, que la surface n'est pas etanche ou que son orientation est
    # incoherente (TetGen peut alors ACCEPTER la surface et produire des tetraedres
    # degeneres : Fighter jet F18) ; reparation aussi apres un refus de TetGen ou un maillage
    # inexploitable. Une surface saine n'est jamais modifiee. Une auto-intersection dans UNE
    # coque garde la cascade historique de `_tetrahedralize` (TetGen pour juge, limite #43).
    from pipeline.mesh_doctor import diagnose_surface

    diagnosis = diagnose_surface(cleaned_mesh)
    report.surface_problems = list(diagnosis.problems)
    pre_repair = (diagnosis.shells > 1 and diagnosis.self_intersecting_faces > 0) or not diagnosis.watertight
    final_mesh = cleaned_mesh
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
    _progress.milestone(0.9)                 # volume tetraedrise
    mesh_data = _build_mesh_data(cleaned_mesh)

    report.tetrahedra_generated = tet_mesh.tet_count

    if check_surface_coverage:
        from core.tet_coverage import surface_outside_volume

        coverage = surface_outside_volume(
            cleaned_mesh.vertices, cleaned_mesh.faces, tet_mesh.tet_vertices, tet_mesh.tetrahedra
        )
        report.surface_outside_volume_fraction = coverage.far_fraction
        report.surface_outside_volume_area = coverage.far_area
        report.surface_max_distance_to_tets = coverage.max_distance
        report.tet_volume_ratio = coverage.volume_ratio
        report.surface_outside_volume_clusters = [(c.area, c.center) for c in coverage.clusters[:10]]

    return mesh_data, tet_mesh, report
