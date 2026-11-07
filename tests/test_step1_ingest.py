"""
tests/test_step1_ingest.py

Tests de element_ardor/step1_ingest.py, en particulier le raffinement local du
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
from element_ardor.step1_ingest import (
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


# ---------------------------------------------------------------------------
# Repli pymeshfix quand la reparation standard (fill_holes/fix_normals/
# fix_winding) ne suffit pas a rendre le mesh etanche - mesure sur
# res/test_models/turkish.stl (64 faces cassees sur 10180, mesh toujours
# non-watertight apres reparation standard, corrige par MeshFix.repair()).
# Sans ce repli, `_tetrahedralize` refusait immediatement (MeshIngestError)
# sans jamais essayer TetGen ni sa propre cascade de reparation reactive
# (limite #43, qui suppose deja un mesh etanche en entree).
# ---------------------------------------------------------------------------


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
