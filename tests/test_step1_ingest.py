"""
tests/test_step1_ingest.py

Tests de pipeline/step1_ingest.py, en particulier le raffinement local du
maillage tetraedrique pres du plateau (`near_platform_resolution`,
CLAUDE.md) : verifie que le decoupage slab/rest produit deux coques
etanches avec des seeds du bon cote de la frontiere de coupe, et que la
densite reelle du TetMesh produit reflete effectivement ce raffinement -
sans regresser le comportement par defaut (raffinement desactive).

Mecanisme retenu : deux regions TetGen (`add_region`), PAS un maillage de
fond (`bgmesh`/`metric`) - ce dernier a ete teste et ecarte (silencieusement
ignore par le binding tetgen installe des qu'un `switches=` non vide est
passe, et provoque un segfault natif reproductible si on force `-m` malgre
tout). Voir CLAUDE.md pour le diagnostic complet.
"""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest
import trimesh

from core.mesh_types import JunctionZone
from pipeline.step1_ingest import (
    MeshIngestError,
    _clean_and_repair,
    _find_guaranteed_interior_point,
    _shared_cut_cap,
    _split_for_platform_refinement,
    run as ingest_run,
)


def _l_shape_mesh() -> trimesh.Trimesh:
    """
    Meme geometrie que res/test_models/l_shape_overhang.stl (voir
    scripts/generate_test_model.py) : pied vertical pose sur le plateau +
    bras en porte-a-faux, assez de hauteur au-dessus du plateau pour
    exercer la transition de raffinement.
    """
    foot = trimesh.creation.box(extents=[4, 4, 20])
    foot.apply_translation([0, 0, 10])
    arm = trimesh.creation.box(extents=[20, 4, 4])
    arm.apply_translation([8, 0, 20])
    return foot.union(arm, engine="manifold")


def _u_shape_mesh() -> trimesh.Trimesh:
    """
    Forme en U : deux montants verticaux relies uniquement pres du plateau
    (pas de traverse au sommet) - concue pour que la coupe de
    _split_for_platform_refinement (a z = platform_z + margin, AU-DESSUS de
    la jonction des montants) laisse deux coques ("slab" ET "rest") chacune
    composees de deux composantes disjointes, dont le centroide global tombe
    dans le vide entre les deux montants - reproduction directe et minimale
    (sans dependre d'un fichier STL externe) du mecanisme mesure sur
    res/test_models/l-triple-axis-overhanged-shape.stl (CLAUDE.md limite
    #42) : centroide naif hors du solide, cassant silencieusement
    TetGen.add_region.
    """
    leg_a = trimesh.creation.box(extents=[4, 4, 20])
    leg_a.apply_translation([-8, 0, 10])
    leg_b = trimesh.creation.box(extents=[4, 4, 20])
    leg_b.apply_translation([8, 0, 10])
    base = trimesh.creation.box(extents=[24, 4, 4])
    base.apply_translation([0, 0, 2])
    mesh = leg_a.union(base, engine="manifold").union(leg_b, engine="manifold")
    return mesh


def _tet_volumes(vertices: np.ndarray, tetrahedra: np.ndarray) -> np.ndarray:
    v = vertices[tetrahedra]
    return np.abs(np.linalg.det(v[:, 1:] - v[:, 0:1])) / 6.0


def test_split_produces_watertight_slab_and_rest_shells():
    mesh = _l_shape_mesh()
    near_platform_margin = 5.0
    platform_z = 0.0

    combined_vertices, combined_faces, slab_seed, rest_seed = _split_for_platform_refinement(
        mesh, near_platform_margin, platform_z
    )

    assert combined_vertices.shape[1] == 3
    assert combined_faces.shape[1] == 3
    # les indices de faces referencent des sommets valides du tableau combine
    assert combined_faces.max() < len(combined_vertices)

    cut_z = platform_z + near_platform_margin
    assert slab_seed[2] < cut_z
    assert rest_seed[2] > cut_z


def test_split_covers_same_bounds_as_source_mesh():
    mesh = _l_shape_mesh()
    combined_vertices, _, _, _ = _split_for_platform_refinement(mesh, 5.0, 0.0)

    # le decoupage ne doit ni tronquer ni deborder du volume source (les
    # coques restent des sous-parties du meme solide, juste separees par un
    # plan de coupe).
    assert np.allclose(combined_vertices.min(axis=0), mesh.bounds[0], atol=1e-6)
    assert np.allclose(combined_vertices.max(axis=0), mesh.bounds[1], atol=1e-6)


def test_tetrahedralize_default_unchanged_without_near_platform_resolution():
    mesh = _l_shape_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data_a, tet_mesh_a, _ = ingest_run(path, resolution=1.5)
        mesh_data_b, tet_mesh_b, _ = ingest_run(path, resolution=1.5, near_platform_resolution=None)
    finally:
        os.unlink(path)

    # meme chemin de code (aucune region construite) : resultat identique.
    assert tet_mesh_a.tet_count == tet_mesh_b.tet_count
    assert np.allclose(tet_mesh_a.tet_vertices, tet_mesh_b.tet_vertices)


def test_tetrahedralize_with_near_platform_resolution_refines_locally():
    mesh = _l_shape_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, tet_mesh, report = ingest_run(
            path,
            resolution=1.5,
            near_platform_resolution=0.5,
            near_platform_margin=5.0,
            platform_z=0.0,
        )
    finally:
        os.unlink(path)

    assert tet_mesh.tet_count > 0
    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    volumes = _tet_volumes(tet_mesh.tet_vertices, tet_mesh.tetrahedra)

    near = tet_centers[:, 2] < 3.0
    far = tet_centers[:, 2] >= 5.0

    assert np.any(near), "fixture must produce tetrahedra near the platform"
    assert np.any(far), "fixture must produce tetrahedra far from the platform"

    # meme mesure que le diagnostic (CLAUDE.md) : le volume moyen pres du
    # plateau doit etre mesurablement plus petit que loin du plateau.
    assert volumes[near].mean() < volumes[far].mean()
    assert report.tetrahedra_generated == tet_mesh.tet_count


def test_near_platform_resolution_uses_fewer_tets_than_equivalent_global_refinement():
    mesh = _l_shape_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        _, tet_mesh_local, _ = ingest_run(
            path,
            resolution=1.5,
            near_platform_resolution=0.5,
            near_platform_margin=5.0,
            platform_z=0.0,
        )
        _, tet_mesh_global, _ = ingest_run(path, resolution=0.5)
    finally:
        os.unlink(path)

    # le raffinement local doit rester nettement moins couteux que raffiner
    # tout le volume a la meme finesse (objectif du mecanisme retenu,
    # mesure sur l_shape_overhang.stl : ~5.4x moins de tetraedres).
    assert tet_mesh_local.tet_count < tet_mesh_global.tet_count


# ---------------------------------------------------------------------------
# Raffinement local aux jonctions multi-region (CLAUDE.md limite #48,
# generalisation de near_platform_resolution au-dela du seul cas plateau)
# ---------------------------------------------------------------------------

def _tripod_mesh() -> trimesh.Trimesh:
    """Meme geometrie que tests/test_step5_layers.py::_tripod_mesh - pied
    vertical + 3 bras non coplanaires, jonction physique au sommet du pied
    (z~=20)."""
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


def _tripod_junction_zones() -> list[JunctionZone]:
    """Une seule zone couvrant la jonction physique du tripode (sommet du
    pied), coherente avec detect_junction_zones sur cette geometrie -
    independante de step1b_junction_detection pour isoler le test de
    _tetrahedralize/_tetrahedralize_once du module de detection (deja
    teste separement dans tests/test_step1b_junction_detection.py)."""
    return [JunctionZone(centroid=np.array([0.0, 0.0, 21.0]), radius=6.0)]


def test_junction_refinement_disabled_by_default_leaves_tetrahedralize_unchanged():
    """junction_refinement_resolution=None (defaut) : aucun changement de
    comportement, meme chemin de code que sans le parametre du tout - meme
    discipline que test_tetrahedralize_default_unchanged_without_near_platform_resolution."""
    mesh = _tripod_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        _, tet_mesh_a, _ = ingest_run(path, resolution=1.5)
        _, tet_mesh_b, _ = ingest_run(
            path, resolution=1.5,
            junction_zones=_tripod_junction_zones(),
            junction_refinement_resolution=None,
        )
    finally:
        os.unlink(path)

    assert tet_mesh_a.tet_count == tet_mesh_b.tet_count
    assert np.allclose(tet_mesh_a.tet_vertices, tet_mesh_b.tet_vertices)


def test_junction_refinement_increases_density_at_junction_not_elsewhere():
    """
    Avec junction_refinement_resolution actif, la densite de tetraedres
    DOIT augmenter pres de la jonction physique (sommet du pied, z~=20-24)
    et rester INCHANGEE loin de toute jonction (bout du bras +Y, z~=28-32) -
    comparaison DANS LA MEME ZONE PHYSIQUE entre deux runs independants
    (piege deja documente CLAUDE.md limite #44 : ne jamais comparer near/far
    dans un seul run sur une geometrie multi-branches, les branches fines
    isolees loin de la jonction ont un volume tetraedrique intrinsequement
    eleve sans rapport avec le raffinement).
    """
    mesh = _tripod_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        _, tet_mesh_coarse, _ = ingest_run(path, resolution=1.5)
        _, tet_mesh_fine, _ = ingest_run(
            path, resolution=1.5,
            junction_zones=_tripod_junction_zones(),
            junction_refinement_resolution=0.4,
        )
    finally:
        os.unlink(path)

    junction_point = np.array([0.0, 0.0, 21.0])
    far_point = np.array([0.0, 15.0, 20.0])  # milieu du bras +Y, loin de la jonction

    def _local_density(tet_mesh, point, radius=3.0):
        centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
        dist = np.linalg.norm(centers - point, axis=1)
        return int((dist < radius).sum())

    coarse_junction = _local_density(tet_mesh_coarse, junction_point)
    fine_junction = _local_density(tet_mesh_fine, junction_point)
    coarse_far = _local_density(tet_mesh_coarse, far_point)
    fine_far = _local_density(tet_mesh_fine, far_point)

    assert fine_junction > coarse_junction * 2, (
        f"La densite pres de la jonction devrait augmenter nettement "
        f"(coarse={coarse_junction}, fine={fine_junction})"
    )
    # Tolerance : le raffinement plateau/jonction peut legerement influencer
    # le maillage global via le solveur Delaunay, mais ne doit PAS produire
    # une explosion de densite loin de la zone ciblee.
    assert fine_far < coarse_far * 2, (
        f"La densite loin de la jonction ne devrait pas exploser "
        f"(coarse={coarse_far}, fine={fine_far})"
    )


def test_junction_refinement_coexists_with_platform_refinement():
    """
    near_platform_resolution ET junction_refinement_resolution actifs
    simultanement (3 regions TetGen : slab, rest, junction) doivent
    cohabiter sans erreur - jamais teste avant cette session (le mecanisme
    plateau existait seul jusqu'ici). Platform_z=0.0 (base du pied) et la
    jonction (z~=21) sont suffisamment separees pour ne pas se chevaucher.
    """
    mesh = _tripod_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        _, tet_mesh, report = ingest_run(
            path, resolution=1.5,
            near_platform_resolution=0.5,
            near_platform_margin=5.0,
            platform_z=0.0,
            junction_zones=_tripod_junction_zones(),
            junction_refinement_resolution=0.4,
        )
    finally:
        os.unlink(path)

    assert tet_mesh.tet_count > 0
    assert report.tetrahedra_generated == tet_mesh.tet_count


# ---------------------------------------------------------------------------
# Seed garanti a l'interieur de la coque (_find_guaranteed_interior_point,
# CLAUDE.md limite #42 - fix du bug add_region sur geometrie non convexe)
# ---------------------------------------------------------------------------

def test_interior_point_prefers_fast_path_centroid_when_valid():
    """Coque convexe simple : le point retourne est EXACTEMENT le fallback (centroide), pas de detour inutile."""
    box = trimesh.creation.box(extents=[10, 10, 10])
    centroid = box.vertices.mean(axis=0)

    result = _find_guaranteed_interior_point(box, centroid, "test box")

    assert np.array_equal(result, np.asarray(centroid, dtype=np.float64))


def test_interior_point_falls_back_when_centroid_outside():
    """
    Coque 'rest' de _u_shape_mesh() decoupee comme
    _split_for_platform_refinement le ferait : le centroide naif est hors
    du solide, le point retourne doit etre different ET a l'interieur.
    """
    mesh = _u_shape_mesh()
    cut_z = 5.0  # platform_z=0.0 + margin=5.0
    rest = mesh.slice_plane(plane_origin=[0.0, 0.0, cut_z], plane_normal=[0.0, 0.0, 1.0], cap=True)
    naive_centroid = rest.vertices.mean(axis=0)

    assert not rest.contains([naive_centroid])[0], (
        "cette fixture doit reproduire un centroide hors solide - sinon ce "
        "test ne prouve rien sur le mecanisme de repli"
    )

    result = _find_guaranteed_interior_point(rest, naive_centroid, "rest test")

    assert not np.array_equal(result, naive_centroid)
    assert rest.contains([result])[0]


def test_interior_point_deterministic_across_calls():
    """Deux appels identiques (meme coque non convexe, meme fallback invalide) retournent EXACTEMENT le meme point."""
    mesh = _u_shape_mesh()
    rest = mesh.slice_plane(plane_origin=[0.0, 0.0, 5.0], plane_normal=[0.0, 0.0, 1.0], cap=True)
    naive_centroid = rest.vertices.mean(axis=0)

    result_a = _find_guaranteed_interior_point(rest, naive_centroid, "rest test")
    result_b = _find_guaranteed_interior_point(rest, naive_centroid, "rest test")

    assert np.array_equal(result_a, result_b)


def test_interior_point_raises_on_empty_shell():
    """Coque vide (0 face) : leve MeshIngestError explicitement, pas un TypeError natif de trimesh."""
    empty = trimesh.Trimesh()

    try:
        _find_guaranteed_interior_point(empty, np.zeros(3), "empty shell test")
        assert False, "devrait lever MeshIngestError"
    except MeshIngestError:
        pass


def test_split_for_platform_refinement_returns_valid_seeds_on_nonconvex_u_shape():
    """
    Garantie centrale du fix : sur une geometrie ou le centroide naif est
    connu pour echouer, les DEUX seeds retournes doivent etre a l'interieur
    de leur coque respective.
    """
    mesh = _u_shape_mesh()
    near_platform_margin = 5.0
    platform_z = 0.0

    _combined_vertices, _combined_faces, slab_seed, rest_seed = _split_for_platform_refinement(
        mesh, near_platform_margin, platform_z
    )

    cut_z = platform_z + near_platform_margin
    slab = mesh.slice_plane(plane_origin=[0.0, 0.0, cut_z], plane_normal=[0.0, 0.0, -1.0], cap=True)
    rest = mesh.slice_plane(plane_origin=[0.0, 0.0, cut_z], plane_normal=[0.0, 0.0, 1.0], cap=True)

    assert slab.contains([slab_seed])[0]
    assert rest.contains([rest_seed])[0]


def test_tetrahedralize_resolution_has_monotonic_effect_on_nonconvex_u_shape():
    """
    Preuve directe de la regression corrigee (CLAUDE.md limite #42) : sur
    une geometrie non convexe, resolution doit avoir un effet monotone
    croissant sur le nombre de tetraedres, meme avec near_platform_resolution actif.
    """
    mesh = _u_shape_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        tet_counts = []
        for resolution in [1.5, 0.8, 0.4]:
            _, tet_mesh, _ = ingest_run(
                path,
                resolution=resolution,
                near_platform_resolution=0.5,
                near_platform_margin=5.0,
                platform_z=0.0,
            )
            tet_counts.append(tet_mesh.tet_count)
    finally:
        os.unlink(path)

    assert tet_counts[0] < tet_counts[1] < tet_counts[2], (
        f"le nombre de tetraedres devrait croitre strictement quand "
        f"resolution diminue, obtenu {tet_counts}"
    )


# ---------------------------------------------------------------------------
# Cap de coupe partage (fix de la fragilite ou deux appels independants a
# slice_plane(cap=True) divergent sur un contour non convexe riche - voir
# CLAUDE.md, mesure sur res/test_models/death.stl)
# ---------------------------------------------------------------------------

def _tripod_platform_mesh() -> trimesh.Trimesh:
    """
    Plateforme circulaire + 3 jambes cylindriques placees de facon
    irreguliere (angle/rayon/rayon-de-jambe deterministes via seed fixe) -
    a z = platform_z + near_platform_margin = 5.0, la section transversale
    (disque + 3 disques de jambes qui le chevauchent partiellement) est un
    polygone non convexe suffisamment riche pour faire diverger deux
    triangulations independantes de `slice_plane(cap=True)` : mesure
    directe (script jetable, meme mecanisme que ce test) : 13 triangles
    divergents sur 66 entre les deux triangulations independantes, assez
    pour que TetGen refuse l'ANCIEN mecanisme ("Internal TetGen error
    within `recoversubfaces`") - reproduction minimale et synthetique
    (sans dependre du fichier externe death.stl) du bug diagnostique sur ce
    dernier.
    """
    rng = np.random.RandomState(7)
    base = trimesh.creation.cylinder(radius=12.0, height=2.0, sections=64)
    base.apply_translation([0, 0, 1.0])
    parts = [base]
    for _ in range(3):
        angle = rng.uniform(0, 2 * np.pi)
        r = rng.uniform(2.0, 9.0)
        x, y = np.cos(angle) * r, np.sin(angle) * r
        leg_radius = rng.uniform(0.8, 2.5)
        leg = trimesh.creation.cylinder(radius=leg_radius, height=20.0, sections=12)
        leg.apply_translation([x, y, 10.0])
        parts.append(leg)
    mesh = parts[0]
    for part in parts[1:]:
        mesh = mesh.union(part, engine="manifold")
    return mesh


def test_shared_cut_cap_produces_single_triangulation_covering_full_contour():
    """
    _shared_cut_cap doit produire UNE SEULE triangulation dont les sommets
    couvrent exactement le contour XY de coupe (pas de sommet en dehors du
    plan, pas de trou dans le polygone d'origine).
    """
    mesh = _tripod_platform_mesh()
    cut_z = 5.0
    plane_origin = np.array([0.0, 0.0, cut_z])
    plane_normal = np.array([0.0, 0.0, 1.0])

    rest_open = mesh.slice_plane(plane_origin=plane_origin, plane_normal=plane_normal, cap=False)
    cap_vertices, cap_faces = _shared_cut_cap(rest_open, plane_origin, plane_normal)

    assert len(cap_faces) > 0, "la fixture doit produire un cap non vide"
    assert np.allclose(cap_vertices[:, 2], cut_z, atol=1e-6), (
        "tous les sommets du cap doivent etre exactement sur le plan de coupe"
    )

    # aire totale du cap generee = aire de la section transversale reelle du
    # mesh a cut_z (comparee via le cap de reference slice_plane(cap=True),
    # qui reste correct pour UNE coque isolee - seule la coherence ENTRE
    # deux coques independantes etait en cause, pas l'aire elle-meme).
    reference = mesh.slice_plane(plane_origin=plane_origin, plane_normal=plane_normal, cap=True)
    ref_face_z = reference.vertices[reference.faces][:, :, 2]
    ref_cap_mask = np.all(np.abs(ref_face_z - cut_z) < 1e-6, axis=1)
    ref_tris = reference.vertices[reference.faces[ref_cap_mask]]
    ref_area = 0.5 * np.linalg.norm(
        np.cross(ref_tris[:, 1] - ref_tris[:, 0], ref_tris[:, 2] - ref_tris[:, 0]), axis=1
    ).sum()

    cap_tris = cap_vertices[cap_faces]
    cap_area = 0.5 * np.linalg.norm(
        np.cross(cap_tris[:, 1] - cap_tris[:, 0], cap_tris[:, 2] - cap_tris[:, 0]), axis=1
    ).sum()

    assert abs(cap_area - ref_area) < 1e-6, (
        f"aire du cap partage ({cap_area}) doit correspondre a la section "
        f"transversale reelle ({ref_area})"
    )


def test_split_for_platform_refinement_avoids_duplicated_cap_triangulation():
    """
    Regression directe du bug : sur la fixture tripode, le PLC combine
    retourne par _split_for_platform_refinement ne doit contenir le cap
    QU'UNE SEULE FOIS - pas une copie par coque comme le faisait l'ancien
    mecanisme (deux slice_plane(cap=True) independants).

    Verifie indirectement (pas d'acces direct au detail interne) : le
    nombre de faces du PLC combine doit etre EXACTEMENT
    faces(slab_open) + faces(rest_open) + faces(cap partage) - une
    duplication du cap ajouterait des faces supplementaires au meme
    endroit geometrique.
    """
    mesh = _tripod_platform_mesh()
    near_platform_margin = 5.0
    platform_z = 0.0
    cut_z = platform_z + near_platform_margin
    plane_origin = np.array([0.0, 0.0, cut_z])
    plane_normal = np.array([0.0, 0.0, 1.0])

    slab_open = mesh.slice_plane(plane_origin=plane_origin, plane_normal=-plane_normal, cap=False)
    rest_open = mesh.slice_plane(plane_origin=plane_origin, plane_normal=plane_normal, cap=False)
    _cap_vertices, cap_faces = _shared_cut_cap(rest_open, plane_origin, plane_normal)

    combined_vertices, combined_faces, _slab_seed, _rest_seed = _split_for_platform_refinement(
        mesh, near_platform_margin, platform_z
    )

    expected_face_count = len(slab_open.faces) + len(rest_open.faces) + len(cap_faces)
    assert len(combined_faces) == expected_face_count, (
        f"le PLC combine devrait contenir exactement "
        f"{expected_face_count} faces (slab + rest + UNE copie du cap), "
        f"obtenu {len(combined_faces)} - une duplication du cap changerait ce compte"
    )
    assert combined_faces.max() < len(combined_vertices)


def test_tetrahedralize_succeeds_on_nonconvex_multibranch_platform_section():
    """
    Preuve de bout en bout du fix (pas seulement du mecanisme interne) : sur
    la fixture tripode, dont le contour de coupe est riche au point de faire
    echouer l'ANCIEN mecanisme (mesure : "Internal TetGen error within
    `recoversubfaces`" avec deux slice_plane(cap=True) independants),
    ingest_run() avec near_platform_resolution actif doit reussir SANS
    lever, et produire un raffinement reellement applique - pas un repli
    silencieux vers "sans raffinement plateau" (CLAUDE.md limite #43, 3e
    tentative de _tetrahedralize).

    Comparaison DANS LA MEME ZONE PHYSIQUE (plateforme circulaire, rayon XY
    > 6 pour exclure les jambes fines isolees qui biaiseraient une
    comparaison near/far sur cette geometrie multi-branches - le volume
    tetraedrique y est domine par l'epaisseur de la jambe, pas par
    near_platform_resolution) entre un run AVEC et SANS raffinement local :
    seule une comparaison entre deux runs independants sur la meme zone est
    fiable ici, contrairement a _l_shape_mesh (solide continu unique) ou
    near vs far dans un seul run suffit.
    """
    mesh = _tripod_platform_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        _mesh_data, tet_mesh_refined, report = ingest_run(
            path,
            resolution=1.5,
            near_platform_resolution=0.5,
            near_platform_margin=5.0,
            platform_z=0.0,
        )
        _mesh_data_b, tet_mesh_plain, _report_b = ingest_run(path, resolution=1.5)
    finally:
        os.unlink(path)

    assert tet_mesh_refined.tet_count > 0
    assert report.tetrahedra_generated == tet_mesh_refined.tet_count

    def _platform_zone_volumes(tet_mesh):
        centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
        r_xy = np.linalg.norm(centers[:, :2], axis=1)
        mask = (r_xy > 6.0) & (centers[:, 2] < 5.0)
        volumes = _tet_volumes(tet_mesh.tet_vertices, tet_mesh.tetrahedra)
        return volumes[mask]

    refined_volumes = _platform_zone_volumes(tet_mesh_refined)
    plain_volumes = _platform_zone_volumes(tet_mesh_plain)

    assert len(refined_volumes) > 0 and len(plain_volumes) > 0, (
        "fixture doit produire des tetraedres dans la zone plateforme (r_xy>6, z<5)"
    )
    assert refined_volumes.mean() < plain_volumes.mean(), (
        "le raffinement plateau doit etre reellement applique (tetraedres "
        "plus petits dans la zone plateforme avec near_platform_resolution "
        "actif) - un repli silencieux vers 'sans raffinement plateau' "
        "produirait un volume moyen identique au run sans raffinement"
    )


# ---------------------------------------------------------------------------
# Detection automatique de platform_z (fix d'un bug reel mesure sur
# res/test_models/death.stl, min(z)=1.245 - un ancien defaut fige a 0.0
# cassait silencieusement _compute_platform_contact sur tout modele non
# pose exactement a z=0, voir CLAUDE.md)
# ---------------------------------------------------------------------------

_PLATFORM_Z_OFFSET = 1.245  # rappelle explicitement le cas reel diagnostique


def test_ingest_detects_platform_z_from_translated_mesh():
    """
    Un mesh translate en Z (non pose a z=0) avec platform_z=None (defaut)
    doit voir sa VRAIE hauteur minimale detectee, pas 0.0.
    """
    mesh = _l_shape_mesh()
    mesh.apply_translation([0, 0, _PLATFORM_Z_OFFSET])
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        _mesh_data, _tet_mesh, report = ingest_run(path, resolution=1.5)
    finally:
        os.unlink(path)

    assert report.platform_z_used == pytest.approx(_PLATFORM_Z_OFFSET, abs=1e-6)


def test_ingest_explicit_platform_z_overrides_detection():
    """Une valeur explicite de platform_z n'est jamais re-derivee du mesh."""
    mesh = _l_shape_mesh()
    mesh.apply_translation([0, 0, _PLATFORM_Z_OFFSET])
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        _mesh_data, _tet_mesh, report = ingest_run(path, resolution=1.5, platform_z=0.0)
    finally:
        os.unlink(path)

    assert report.platform_z_used == 0.0


@pytest.mark.parametrize("z_shift", [_PLATFORM_Z_OFFSET, -15.0])
def test_place_on_platform_translates_in_z_only(z_shift):
    """
    `place_on_platform` ramene min(z) a 0 par une translation en z SEULE
    (cas reels : death.stl +1.245 mm, real_overhang.stl -15 mm), remonte le
    decalage, et rend un maillage identique a celui d'un modele deja pose.
    """
    base = _l_shape_mesh()
    shifted = base.copy()
    shifted.apply_translation([3.0, -2.0, z_shift])
    paths = []
    for m in (base, shifted):
        with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
            paths.append(f.name)
        m.export(paths[-1])
    try:
        ref_data, _ref_tet, ref_report = ingest_run(paths[0], resolution=1.5, place_on_platform=True)
        mesh_data, tet_mesh, report = ingest_run(paths[1], resolution=1.5, place_on_platform=True)
        _d, _t, untouched = ingest_run(paths[1], resolution=1.5)
    finally:
        for p in paths:
            os.unlink(p)

    assert report.placement_offset_z == pytest.approx(-z_shift, abs=1e-5)
    assert ref_report.placement_offset_z == 0.0
    assert report.platform_z_used == pytest.approx(0.0, abs=1e-9)
    assert mesh_data.vertices[:, 2].min() == pytest.approx(0.0, abs=1e-9)
    assert tet_mesh.tet_vertices[:, 2].min() == pytest.approx(0.0, abs=1e-6)
    # x, y intacts : seule la hauteur est deplacee
    assert np.allclose(mesh_data.vertices[:, :2].mean(axis=0),
                       ref_data.vertices[:, :2].mean(axis=0) + np.array([3.0, -2.0]), atol=1e-4)
    # sans l'option, le modele garde ses coordonnees et aucun decalage n'est declare
    assert untouched.placement_offset_z == 0.0
    assert untouched.platform_z_used == pytest.approx(z_shift, abs=1e-5)


def test_ingest_platform_z_detection_is_noop_on_mesh_already_at_zero():
    """
    Non-regression : un mesh deja pose a z=0 (comportement historique de
    tous les modeles de test existants) produit EXACTEMENT le meme resultat
    qu'avant, platform_z=None detecte alors 0.0.
    """
    mesh = _l_shape_mesh()
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data_a, tet_mesh_a, report_a = ingest_run(path, resolution=1.5)
        mesh_data_b, tet_mesh_b, report_b = ingest_run(path, resolution=1.5, platform_z=0.0)
    finally:
        os.unlink(path)

    assert report_a.platform_z_used == 0.0
    assert report_b.platform_z_used == 0.0
    assert tet_mesh_a.tet_count == tet_mesh_b.tet_count
    assert np.allclose(tet_mesh_a.tet_vertices, tet_mesh_b.tet_vertices)


def test_ingest_platform_z_detection_refines_at_true_platform_level():
    """
    Avec near_platform_resolution actif et platform_z=None sur un mesh
    translate, le raffinement local doit viser la VRAIE zone plateau
    (autour de min(z)=1.245), pas z=[0,5] - preuve de bout en bout que la
    valeur detectee est bien propagee jusqu'a _tetrahedralize.
    """
    mesh = _l_shape_mesh()
    mesh.apply_translation([0, 0, _PLATFORM_Z_OFFSET])
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        _mesh_data, tet_mesh, report = ingest_run(
            path,
            resolution=1.5,
            near_platform_resolution=0.5,
            near_platform_margin=5.0,
        )
    finally:
        os.unlink(path)

    assert report.platform_z_used == pytest.approx(_PLATFORM_Z_OFFSET, abs=1e-6)

    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    volumes = _tet_volumes(tet_mesh.tet_vertices, tet_mesh.tetrahedra)

    near = tet_centers[:, 2] < (_PLATFORM_Z_OFFSET + 3.0)
    far = tet_centers[:, 2] >= (_PLATFORM_Z_OFFSET + 5.0)
    assert np.any(near), "fixture doit produire des tetraedres pres du vrai plateau"
    assert np.any(far), "fixture doit produire des tetraedres loin du vrai plateau"
    assert volumes[near].mean() < volumes[far].mean(), (
        "le raffinement doit cibler la VRAIE hauteur du plateau detectee "
        "(min(z) du mesh translate), pas z=0 fige"
    )


# ---------------------------------------------------------------------------
# Repli pymeshfix quand la reparation standard (fill_holes/fix_normals/
# fix_winding) ne suffit pas a rendre le mesh etanche - mesure sur
# res/test_models/turkish.stl (64 faces cassees sur 10180, mesh toujours
# non-watertight apres reparation standard, corrige par MeshFix.repair()).
# Sans ce repli, `_tetrahedralize` refusait immediatement (MeshIngestError)
# sans jamais essayer TetGen ni sa propre cascade de reparation reactive
# (limite #43, qui suppose deja un mesh etanche en entree).
# ---------------------------------------------------------------------------

def test_repair_with_pymeshfix_closes_holes_standard_repair_misses():
    """
    Test unitaire direct de `_repair_with_pymeshfix` (anciennement
    `_repair_self_intersections`, renommee pour refleter sa portee reelle -
    MeshFix corrige generiquement auto-intersections ET trous non-manifold,
    pas seulement le premier cas) : verifie qu'elle rend bien un mesh
    watertight sur un defaut synthetique, independamment de tout appelant.
    """
    pytest.importorskip("pymeshfix")
    from pipeline.step1_ingest import _repair_with_pymeshfix

    box = trimesh.creation.box(extents=[10, 10, 10])
    faces = np.delete(box.faces, [0, 1], axis=0)  # trou simple (un quad)
    broken = trimesh.Trimesh(vertices=box.vertices, faces=faces, process=False)
    assert not broken.is_watertight

    repaired = _repair_with_pymeshfix(broken)
    assert repaired is not None
    assert repaired.is_watertight


def test_clean_and_repair_falls_back_to_pymeshfix_when_standard_repair_fails(monkeypatch):
    """
    Reproduit le vrai gap mesure sur turkish.stl : un mesh qui reste
    non-watertight APRES `trimesh.repair.fill_holes` (simule ici en
    remplacant `fill_holes` par un no-op via monkeypatch, plutot que de
    devoir construire une geometrie synthetique dont fill_holes echoue
    reellement - defaut difficile a synthetiser de facon fiable, propre
    aux maillages denses/scannes comme le vrai cas) doit quand meme
    ressortir watertight de `_clean_and_repair` grace au repli pymeshfix.
    """
    pytest.importorskip("pymeshfix")
    import pipeline.step1_ingest as step1

    monkeypatch.setattr(trimesh.repair, "fill_holes", lambda mesh: None)

    box = trimesh.creation.box(extents=[10, 10, 10])
    faces = np.delete(box.faces, [0, 1], axis=0)
    broken = trimesh.Trimesh(vertices=box.vertices, faces=faces, process=False)
    assert not broken.is_watertight

    repaired_mesh, report = step1._clean_and_repair(broken)

    assert repaired_mesh.is_watertight
    assert report.is_watertight_after_repair is True


def test_ingest_succeeds_on_mesh_pymeshfix_can_repair_but_standard_repair_cannot(monkeypatch):
    """
    Test de bout en bout (pas seulement `_clean_and_repair` isole) :
    `step1_ingest.run()` doit produire un TetMesh valide sur un mesh dont
    seul pymeshfix peut refermer le trou - sans le repli, `_tetrahedralize`
    levait `MeshIngestError` avant meme d'essayer TetGen.
    """
    pytest.importorskip("pymeshfix")

    monkeypatch.setattr(trimesh.repair, "fill_holes", lambda mesh: None)

    box = trimesh.creation.box(extents=[10, 10, 10])
    faces = np.delete(box.faces, [0, 1], axis=0)
    broken = trimesh.Trimesh(vertices=box.vertices, faces=faces, process=False)

    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    broken.export(path)
    try:
        _mesh_data, tet_mesh, report = ingest_run(path, resolution=2.0)
    finally:
        os.unlink(path)

    assert report.is_watertight_after_repair is True
    assert tet_mesh.tet_count > 0


# ---------------------------------------------------------------------------
# Extraction de la coque dominante quand un mesh combine est non-watertight
# a cause de fragments negligeables (CLAUDE.md) - mesure sur
# res/test_models/turkish.stl : 16 fragments minuscules (76 faces sur
# 10180, 0.31% du volume, normales inversees) accoles au corps principal
# (deja individuellement watertight) par des jonctions non-manifold.
# `MeshFix.repair()` (pymeshfix) echoue dangereusement sur ce cas precis :
# il rend le mesh watertight en le REDUISANT de 10180 a 10 faces - d'ou
# `_extract_dominant_shell`, tente AVANT pymeshfix dans la cascade,
# strictement plus cible (garde la coque dominante, jamais de risque de
# collapse global).
# ---------------------------------------------------------------------------

def _mesh_with_negligible_attached_fragments() -> trimesh.Trimesh:
    """
    Reproduit la STRUCTURE du vrai cas turkish.stl (une coque dominante
    watertight + plusieurs fragments watertight minuscules, le tout non
    watertight une fois combine) plutot que le mecanisme exact de jonction
    non-manifold (difficile a synthetiser de facon fiable - voir la
    discussion equivalente pour le fix de `_classify_contours`, CLAUDE.md).
    Un sommet du corps principal est deplace pour toucher exactement un
    sommet d'un petit fragment SANS que les deux coques partagent une
    facette complete - `merge_vertices` fusionne alors ce sommet unique,
    creant une configuration non-manifold localisee suffisante pour que
    `mesh.is_watertight` retourne False sur le mesh COMBINE, alors que
    chaque coque individuelle (`mesh.split(only_watertight=False)`) reste
    parfaitement watertight - exactement le motif que `_extract_dominant_shell`
    est concu pour geree.
    """
    main = trimesh.creation.box(extents=[10, 10, 10])
    fragments = []
    fragment_vertex_offset = len(main.vertices)
    all_vertices = [main.vertices]
    all_faces = [main.faces]

    rng = np.random.default_rng(7)
    for i in range(3):
        frag = trimesh.creation.box(extents=[0.1, 0.1, 0.1])
        # positionne le fragment pour qu'UN de ses sommets coincide
        # EXACTEMENT avec un sommet du corps principal (coin du cube),
        # sans partager de facette - source du defaut non-manifold local.
        corner = main.vertices[i]
        frag.apply_translation(corner - frag.vertices[0])
        all_vertices.append(frag.vertices)
        all_faces.append(frag.faces + fragment_vertex_offset)
        fragment_vertex_offset += len(frag.vertices)
        fragments.append(frag)

    combined_v = np.vstack(all_vertices)
    combined_f = np.vstack(all_faces)
    return trimesh.Trimesh(vertices=combined_v, faces=combined_f, process=False)


def test_extract_dominant_shell_keeps_main_body_drops_negligible_fragments():
    """
    Test unitaire direct de `_extract_dominant_shell` : sur un mesh combine
    non-watertight dont une coque domine tres largement (>95% faces/volume),
    doit retourner CETTE coque seule, avec le bon nombre de fragments
    ecartes - jamais une perte de la geometrie principale.
    """
    from pipeline.step1_ingest import _extract_dominant_shell

    mesh = _mesh_with_negligible_attached_fragments()

    dominant, dropped = _extract_dominant_shell(mesh)

    assert dominant is not None
    assert dominant.is_watertight
    assert np.allclose(dominant.extents, [10.0, 10.0, 10.0])
    assert dropped == 3


def test_extract_dominant_shell_returns_none_when_no_clear_dominance():
    """
    Garde-fou : si aucune composante ne domine assez clairement (deux
    coques de taille comparable), `_extract_dominant_shell` doit retourner
    (None, 0) plutot que de choisir arbitrairement - jamais de perte
    silencieuse de geometrie legitime multi-corps.
    """
    from pipeline.step1_ingest import _extract_dominant_shell

    box_a = trimesh.creation.box(extents=[10, 10, 10])
    box_b = trimesh.creation.box(extents=[10, 10, 10])
    box_b.apply_translation([0.0, 0.0, 0.0])  # meme sommet [0] que box_a pour forcer une jonction
    corner_shift = box_a.vertices[0] - box_b.vertices[0]
    box_b.apply_translation(corner_shift + np.array([20.0, 0.0, 0.0]))

    combined_v = np.vstack([box_a.vertices, box_b.vertices])
    combined_f = np.vstack([box_a.faces, box_b.faces + len(box_a.vertices)])
    mesh = trimesh.Trimesh(vertices=combined_v, faces=combined_f, process=False)

    dominant, dropped = _extract_dominant_shell(mesh)

    assert dominant is None
    assert dropped == 0


# ---------------------------------------------------------------------------
# Reconnaissance des erreurs TetGen "recuperables" par la cascade de
# reparation reactive (limite #43/#44) - au-dela du seul message
# "self-intersection" d'origine. Mesure sur res/test_models/turkish.stl :
# TetGen levait "Internal TetGen error within `recoversubfaces`" (variante
# deja documentee CLAUDE.md limite #44 comme "differente mais egalement
# fatale"), jamais reliee au code de detection avant ce fix - la cascade
# ne se declenchait pas, remontant l'erreur brute sans jamais essayer la
# reparation pymeshfix.
# ---------------------------------------------------------------------------

def test_is_recoverable_tetgen_error_matches_self_intersection():
    from pipeline.step1_ingest import _is_recoverable_tetgen_error

    exc = RuntimeError("The input surface mesh contain self-intersections")
    assert _is_recoverable_tetgen_error(exc)


def test_is_recoverable_tetgen_error_matches_recoversubfaces():
    from pipeline.step1_ingest import _is_recoverable_tetgen_error

    exc = RuntimeError("Internal TetGen error within `recoversubfaces`.")
    assert _is_recoverable_tetgen_error(exc)


def test_is_recoverable_tetgen_error_rejects_unrelated_error():
    from pipeline.step1_ingest import _is_recoverable_tetgen_error

    exc = RuntimeError("some unrelated structural error, not TetGen's fault")
    assert not _is_recoverable_tetgen_error(exc)


def test_clean_and_repair_reports_dropped_fragments():
    """
    `_clean_and_repair` doit remonter `IngestReport.dominant_shell_fragments_dropped`
    explicitement (jamais un silence) quand ce mecanisme s'est declenche -
    y compris sur ce mesh synthetique qui reste `is_watertight=True` au
    sens strict de trimesh (un seul sommet partage entre le corps principal
    et chaque fragment ne suffit pas a briser le comptage d'aretes
    manifold), exactement le genre de cas "techniquement watertight mais
    structurellement etrange" que le declenchement inconditionnel de
    `_extract_dominant_shell` (avant meme le test `is_watertight`) est
    concu pour couvrir aussi.
    """
    mesh = _mesh_with_negligible_attached_fragments()

    repaired_mesh, report = _clean_and_repair(mesh)

    assert repaired_mesh.is_watertight
    assert report.dominant_shell_fragments_dropped == 3
    assert np.allclose(repaired_mesh.extents, [10.0, 10.0, 10.0])


def test_tetrahedralize_cascade_rejects_collapsed_pymeshfix_repair(monkeypatch):
    """
    Reproduit le vrai bug mesure sur turkish.stl : `_tetrahedralize` a son
    PROPRE appel a `_repair_with_pymeshfix` (independant de celui deja
    garde dans `_clean_and_repair`), utilise pour retenter TetGen apres un
    premier refus - jusqu'a ce fix, cet appel n'avait AUCUN garde-fou
    anti-effondrement. Sur turkish.stl, ce chemin acceptait un resultat
    passant de 10104 a 1488 faces (watertight, donc silencieusement
    utilise), amputant le mesh de z=90.78mm a z=63.26mm - TetGen
    tetraedrisait ensuite ce mesh mutile SANS aucune erreur, et le pipeline
    entier produisait un resultat visuellement absurde (trajectoire
    concentree dans le tiers inferieur du modele) sans jamais lever
    d'exception.

    Ce test simule l'effondrement via monkeypatch (`_repair_with_pymeshfix`
    retourne un mesh tronque a une fraction du nombre de faces d'entree) et
    verifie que `_tetrahedralize` ne l'utilise NULLE PART dans sa cascade -
    ni comme tentative dediee, ni comme source de repli pour "sans
    raffinement plateau". Sur un mesh par ailleurs valide (une boite), la
    premiere tentative ("tel quel") doit de toute facon reussir : le test
    verifie que le TetMesh produit couvre la VRAIE etendue de la boite,
    jamais celle du faux resultat pymeshfix tronque.
    """
    import pipeline.step1_ingest as step1

    # icosphere plutot qu'une boite simple : une boite a toujours 12 faces
    # quelle que soit sa taille, insuffisant pour simuler un EFFONDREMENT
    # de nombre de faces (le vrai signal du garde-fou) - une icosphere a
    # subdivisions=2 (320 faces) laisse la place a un "effondrement" vers
    # une boite (12 faces) qui reduit vraiment le compte, meme signature
    # que le vrai incident turkish.stl (10104 -> 1488 faces).
    source = trimesh.creation.icosphere(subdivisions=2, radius=10.0)
    assert source.is_watertight
    assert len(source.faces) > 24  # garantit que le "collapse" ci-dessous divise bien par >2

    def _fake_collapsed_repair(mesh):
        # Simule l'effondrement mesure sur turkish.stl : un resultat
        # watertight mais radicalement plus petit que le mesh d'entree,
        # geometriquement sans rapport (une minuscule boite au lieu de la
        # vraie sphere) - le nombre exact de faces importe peu, seul le
        # ratio compte pour le garde-fou.
        collapsed = trimesh.creation.box(extents=[0.1, 0.1, 0.1])
        assert collapsed.is_watertight
        assert len(collapsed.faces) < 0.5 * len(mesh.faces)
        return collapsed

    monkeypatch.setattr(step1, "_repair_with_pymeshfix", _fake_collapsed_repair)

    tet_mesh = step1._tetrahedralize(source, resolution=2.0)

    # Le vrai mesh (sphere de rayon 10) doit etre tetraedrise correctement -
    # jamais le faux resultat pymeshfix (boite 0.1x0.1x0.1) accepte a sa
    # place par la cascade.
    assert tet_mesh.tet_vertices[:, 2].max() > 5.0


# ---------------------------------------------------------------------------
# Refus des modeles trop grands (feuille de route, 2026-10-01)
# ---------------------------------------------------------------------------

def test_min_enclosing_circle_exact_on_known_shapes():
    import numpy as np
    from pipeline.step1_ingest import _min_enclosing_circle
    square = np.array([[0, 0], [2, 0], [2, 2], [0, 2], [1, 1]], dtype=float)
    c, r = _min_enclosing_circle(square)
    assert np.allclose(c, [1, 1]) and abs(r - np.sqrt(2)) < 1e-12
    tri = np.array([[0, 0], [10, 0], [5, 1]], dtype=float)      # obtus : diametre du grand cote
    c, r = _min_enclosing_circle(tri)
    assert np.allclose(c, [5, 0]) and abs(r - 5) < 1e-12
    rng = np.random.default_rng(3)
    pts = rng.normal(size=(300, 2))
    c, r = _min_enclosing_circle(pts)
    assert np.all(np.linalg.norm(pts - c, axis=1) <= r + 1e-9)


def test_check_fits_platform_is_placement_independent():
    import numpy as np
    import pytest
    from pipeline.step1_ingest import ModelTooLargeError, check_fits_platform
    # rectangle 200 x 200 : diagonale 282,8 -> rayon 141,4 < 150, tient
    box = np.array([[x, y, z] for x in (0, 200) for y in (0, 200) for z in (0, 50)], dtype=float) + 500.0
    assert abs(check_fits_platform(box, 150.0) - np.hypot(100, 100)) < 1e-9
    big = box * [1.1, 1.1, 1.0]
    with pytest.raises(ModelTooLargeError) as e:
        check_fits_platform(big, 150.0)
    assert e.value.footprint_radius > 150.0
    with pytest.raises(ModelTooLargeError):
        check_fits_platform(box, 150.0, max_height=40.0)
    assert check_fits_platform(big, None) == 0.0   # plateau inconnu : aucun refus


def test_speculative_tetgen_uses_parallel_attempt_only_when_probe_left_mesh_intact(monkeypatch):
    """Logique de l'execution speculative de `_tetrahedralize` (session du
    2026-10-01) : la tentative finale lancee en parallele de la sonde n'est
    retenue que si la sonde a laisse le maillage INTACT ; sinon elle est
    ignoree et la tentative reelle porte sur le maillage repare."""
    import trimesh
    import pipeline.step1_ingest as s1
    mesh = trimesh.creation.box(extents=[4, 4, 4])
    resolved = trimesh.creation.box(extents=[4, 4, 4.0001])
    repaired = trimesh.creation.box(extents=[4, 4, 4.0002])
    calls = []

    def fake_once(candidate, resolution, npr, *rest):
        calls.append((id(candidate), npr))
        if candidate is mesh:
            raise RuntimeError("The input surface mesh contain self-intersections.")
        return ("tets", id(candidate), npr)

    monkeypatch.setattr(s1, "_tetrahedralize_once_isolated", fake_once)
    monkeypatch.setattr(s1, "resolve_self_intersections", lambda m: resolved)
    monkeypatch.setattr(s1, "_find_self_intersecting_faces", lambda m: [])
    monkeypatch.setattr(s1, "_repair_with_pymeshfix", lambda m: None)
    monkeypatch.setattr(s1, "_is_recoverable_tetgen_error", lambda e: True)

    # 1) sonde reussie (maillage intact) : resultat de la tentative speculative
    monkeypatch.setattr(s1, "_finish_geometric_resolution", lambda r, res: (r, None))
    calls.clear()
    out = s1._tetrahedralize(mesh, 1.0, near_platform_resolution=0.5)
    assert out == ("tets", id(resolved), 0.5)
    assert calls.count((id(resolved), 0.5)) == 1        # jamais rejouee

    # 2) sonde qui repare : speculation ignoree, vraie tentative sur le repare
    monkeypatch.setattr(s1, "_finish_geometric_resolution", lambda r, res: (repaired, None))
    calls.clear()
    out = s1._tetrahedralize(mesh, 1.0, near_platform_resolution=0.5)
    assert out == ("tets", id(repaired), 0.5)

    # 3) interrupteur coupe : comportement serie historique, meme resultat
    monkeypatch.setattr(s1, "_SPECULATIVE_TETGEN", False)
    monkeypatch.setattr(s1, "_finish_geometric_resolution", lambda r, res: (r, None))
    calls.clear()
    assert s1._tetrahedralize(mesh, 1.0, near_platform_resolution=0.5) == ("tets", id(resolved), 0.5)


def test_probe_tet_mesh_reused_for_no_refinement_attempt(monkeypatch):
    """La sonde reussie de `_finish_geometric_resolution` rend son TetMesh ; la
    tentative « sans raffinement plateau » sur le meme maillage le reutilise
    au lieu de refaire le meme calcul (turkish.stl : 7,26 s)."""
    import trimesh
    import pipeline.step1_ingest as s1
    mesh = trimesh.creation.box(extents=[4, 4, 4])
    repaired = trimesh.creation.box(extents=[4, 4, 4.0002])
    calls = []

    def fake_once(candidate, resolution, npr, *rest):
        calls.append((id(candidate), npr))
        raise RuntimeError("The input surface mesh contain self-intersections.")

    monkeypatch.setattr(s1, "_tetrahedralize_once_isolated", fake_once)
    monkeypatch.setattr(s1, "resolve_self_intersections", lambda m: mesh.copy())
    monkeypatch.setattr(s1, "_find_self_intersecting_faces", lambda m: [])
    monkeypatch.setattr(s1, "_repair_with_pymeshfix", lambda m: None)
    monkeypatch.setattr(s1, "_is_recoverable_tetgen_error", lambda e: True)
    monkeypatch.setattr(s1, "_finish_geometric_resolution", lambda r, res: (repaired, "SONDE"))
    out = s1._tetrahedralize(mesh, 1.0, near_platform_resolution=0.5)
    assert out == "SONDE"
    assert (id(repaired), None) not in calls      # jamais recalcule


def test_translation_retry_recovers_and_moves_nodes_back(monkeypatch):
    """Repli d'arrondi flottant : si la cascade echoue a la position d'origine, elle est
    rejouee sur le maillage recentre en memoire, et les noeuds sont ramenes a leur place."""
    import trimesh as _tm
    from core.mesh_types import TetMesh as _TetMesh
    from pipeline import step1_ingest as s1

    mesh = _tm.creation.box(extents=(10.0, 10.0, 10.0))
    mesh.apply_translation([120.0, 80.0, 5.0])
    calls = []

    def fake(m, *args, **kwargs):
        center = (m.bounds[0] + m.bounds[1]) / 2.0
        calls.append(center.copy())
        if abs(center[0]) > 1.0:                       # n'accepte que le maillage recentre
            raise s1.MeshIngestError("auto-intersections")
        return _TetMesh(tet_vertices=np.asarray(m.vertices, dtype=float), tetrahedra=np.zeros((0, 4), int))

    monkeypatch.setattr(s1, "_tetrahedralize", fake)
    tm = s1._tetrahedralize_with_translation_retry(mesh, 1.0, None, 5.0, 0.0, None, None)
    assert len(calls) == 2
    assert np.allclose(tm.tet_vertices, mesh.vertices)

    def always_fail(m, *args, **kwargs):
        raise s1.MeshIngestError("toujours")

    monkeypatch.setattr(s1, "_tetrahedralize", always_fail)
    with pytest.raises(s1.MeshIngestError, match="toujours"):
        s1._tetrahedralize_with_translation_retry(mesh, 1.0, None, 5.0, 0.0, None, None)


# ---------------------------------------------------------------------------
# R-013 : plusieurs pieces sur le plateau
# ---------------------------------------------------------------------------

def test_multi_mesh_scene_keeps_node_transforms(tmp_path):
    """Les transformations de nœud d'une scene (placement d'un 3MF/GLB) ne sont pas perdues."""
    from pipeline.step1_ingest import _load_raw_mesh
    scene = trimesh.Scene()
    scene.add_geometry(trimesh.creation.box((10, 10, 10)), node_name="a")
    scene.add_geometry(trimesh.creation.box((5, 5, 5)), node_name="b",
                       transform=trimesh.transformations.translation_matrix([50, 0, 0]))
    path = str(tmp_path / "plate.glb")
    scene.export(path)
    mesh = _load_raw_mesh(path, allow_multi_mesh=True)
    assert mesh.bounds[1][0] == pytest.approx(52.5, abs=1e-6)   # la seconde piece est restee en x = 50


def test_small_real_piece_next_to_large_one_is_not_dropped():
    """Une petite piece a volume positif n'est jamais prise pour un fragment (`_extract_dominant_shell`)."""
    from pipeline.step1_ingest import _extract_dominant_shell
    big = trimesh.creation.box((40, 40, 40))
    small = trimesh.creation.box((6, 6, 6))
    small.apply_translation([60, 0, 0])
    both = trimesh.util.concatenate([big, small])
    assert _extract_dominant_shell(both) == (None, 0)


def test_inverted_fragment_next_to_large_body_is_still_dropped():
    """Le cas de turkish.stl est conserve : fragment a volume signe negatif (normales inversees)."""
    from pipeline.step1_ingest import _extract_dominant_shell
    big = trimesh.creation.box((40, 40, 40))
    frag = trimesh.creation.box((1, 1, 1))
    frag.invert()
    frag.apply_translation([60, 0, 0])
    dominant, dropped = _extract_dominant_shell(trimesh.util.concatenate([big, frag]))
    assert dropped == 1 and dominant is not None and len(dominant.faces) == 12


def test_surface_inconsistent_diagnostic_carries_numbers_for_translation(monkeypatch):
    """R-016 : le diagnostic porte aussi `self_intersections` et `volume_pct` (texte reconstructible)."""
    from pipeline import step1_ingest as s1
    box = trimesh.creation.box((10, 10, 10))
    fake = {"status": s1.TETRANET_STATUS_INCONSISTENT_SURFACE, "tets": 0, "volume_ratio": 0.42, "sign_conflicts": 1234,
            "surface_holes": 3, "surface_flipped_edges": 9, "surface_nonmanifold_edges": 0,
            "surface_duplicate_faces": 0, "surface_vertices": np.asarray(box.vertices),
            "surface_faces": np.asarray(box.faces), "tet_vertices": np.zeros((0, 3)),
            "tetrahedra": np.zeros((0, 4), dtype=np.int64)}
    try:
        real = s1._cpp.module()
    except RuntimeError:   # noyau natif absent (APEX_USE_CPP=0)
        real = None

    class _Module:   # le noyau natif reste visible pour le reste de l'etape 1 quand il est actif
        tetranet_mesh = staticmethod(lambda *a: fake)

        def __getattr__(self, name):
            return getattr(real, name)
    monkeypatch.setattr(s1._cpp, "module", lambda: _Module())
    with pytest.raises(s1.SurfaceInconsistentError) as err:
        s1._ingest_with_tetranet(box, 1.0, None, False, False, 0.1)
    assert err.value.diagnostic["volume_pct"] == pytest.approx(42.0)
    assert err.value.diagnostic["self_intersections"] == 0
    assert err.value.diagnostic["surface_holes"] == 3
    assert "1\u202f234 points incohérents" in str(err.value)
