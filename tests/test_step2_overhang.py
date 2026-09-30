"""
tests/test_step2_overhang.py

Tests pour pipeline/step2_overhang.py, sur des mesh synthetiques generes a
la volee (aucun fichier externe requis).

Execution : pytest tests/test_step2_overhang.py -v
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import OverhangParams, run as overhang_run


# ---------------------------------------------------------------------------
# Fixtures : mesh synthetiques
# ---------------------------------------------------------------------------

def _mesh_data_from_trimesh(mesh: trimesh.Trimesh, resolution: float = 1.5):
    """Passe par step1_ingest.run() via un fichier temporaire, pour tester
    step2 exactement dans les conditions du pipeline reel (MeshData avec
    bvh/kdtree deja construits)."""
    import tempfile
    import os

    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, _, _ = ingest_run(path, resolution=resolution)
    finally:
        os.unlink(path)
    return mesh_data


@pytest.fixture
def flat_box_mesh_data():
    """Boite simple posee a plat : aucune face ne devrait etre en surplomb ni pont."""
    mesh = trimesh.creation.box(extents=[10, 10, 10])
    return _mesh_data_from_trimesh(mesh)


@pytest.fixture
def overhang_l_shape_mesh_data():
    """
    Forme en L couche (cf. schema discute dans la conversation) : une base
    verticale + un bras horizontal en porte-a-faux. Le dessous du bras
    horizontal doit etre detecte comme surplomb.
    """
    vertical = trimesh.creation.box(extents=[4, 4, 20])
    vertical.apply_translation([0, 0, 10])

    horizontal = trimesh.creation.box(extents=[20, 4, 4])
    horizontal.apply_translation([8, 0, 20])

    combined = vertical.union(horizontal, engine="manifold")
    return _mesh_data_from_trimesh(combined, resolution=2.0)


@pytest.fixture
def hollow_box_mesh_data():
    """
    Boite creuse avec cavite interne : le "plafond" de la cavite est
    horizontal (donc pas is_overhang) mais n'a rien en dessous (donc
    is_bridge attendu).
    """
    outer = trimesh.creation.box(extents=[10, 10, 10])
    inner = trimesh.creation.box(extents=[4, 4, 4])
    hollow = outer.difference(inner, engine="manifold")
    return _mesh_data_from_trimesh(hollow, resolution=1.5)


@pytest.fixture
def sphere_mesh_data():
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=5.0)
    return _mesh_data_from_trimesh(mesh)


# ---------------------------------------------------------------------------
# Cas 1 - Cube pose a plat : aucun surplomb, aucun pont
# ---------------------------------------------------------------------------

def test_flat_box_sides_and_top_have_no_overhang(flat_box_mesh_data):
    """
    Note : la face du DESSOUS d'un cube pose a plat pointe geometriquement
    a angle 0 de la verticale descendante et n'a, du point de vue du
    ray-casting, rien en dessous d'elle dans le mesh (le plateau
    d'impression n'est pas un objet du mesh). Cette etape est donc
    incapable, par construction, de distinguer "repose sur le plateau" de
    "en surplomb dans le vide" - elle marque les deux is_overhang=True et
    potentiellement is_bridge=True. C'est une limite connue de cette etape
    seule : la connaissance du plateau d'impression (sa position Z) devra
    etre injectee soit ici via un parametre explicite, soit filtree a une
    etape ulterieure du pipeline (ex: step3_segmentation, qui connait la
    position du plateau via MachineProfile). Ce test verifie donc les
    faces laterales et superieure uniquement, qui ne devraient jamais etre
    en surplomb ni en pont sur un cube droit plein.
    """
    params = OverhangParams(overhang_angle_threshold_deg=45.0)
    result = overhang_run(flat_box_mesh_data, params)

    face_centers = flat_box_mesh_data.vertices[flat_box_mesh_data.faces].mean(axis=1)
    is_bottom_face = face_centers[:, 2] < face_centers[:, 2].min() + 1e-3

    non_bottom_overhang = result.is_overhang & (~is_bottom_face)
    assert not np.any(non_bottom_overhang), (
        f"Seule la face du dessous devrait etre marquee surplomb sur un "
        f"cube droit, mais {non_bottom_overhang.sum()} face(s) laterale/"
        f"superieure sont aussi marquees"
    )
    non_bottom_bridge = result.is_bridge & (~is_bottom_face)
    assert not np.any(non_bottom_bridge), (
        f"Un cube plein pose a plat ne devrait avoir aucun pont interne "
        f"hors face du dessous (limite connue, cf. docstring), trouve "
        f"{non_bottom_bridge.sum()} face(s) laterale/superieure marquees"
    )


# ---------------------------------------------------------------------------
# Cas 2 - Forme en L : le bras horizontal en porte-a-faux doit etre detecte
# ---------------------------------------------------------------------------

def test_l_shape_detects_overhang_underside(overhang_l_shape_mesh_data):
    params = OverhangParams(overhang_angle_threshold_deg=45.0)
    result = overhang_run(overhang_l_shape_mesh_data, params)

    assert np.any(result.is_overhang), (
        "La face inferieure du bras horizontal en porte-a-faux devrait etre "
        "detectee comme surplomb, mais aucune face n'est marquee"
    )


# ---------------------------------------------------------------------------
# Cas 3 - Boite creuse : le plafond interne doit etre un pont, pas un surplomb
# ---------------------------------------------------------------------------

def test_hollow_box_ceiling_is_bridge(hollow_box_mesh_data):
    """
    Note de conception : is_bridge et is_overhang ne sont PAS mutuellement
    exclusifs. Un plafond de cavite quasi-horizontal (normale pointant vers
    le bas, angle proche de 0) est a la fois "en surplomb" par un critere
    d'angle simple ET testable comme "pont" par le mecanisme de ray-cast -
    les deux detections capturent des aspects differents du meme risque.
    Ce qui distingue vraiment un plafond de cavite d'un simple surplomb en
    pente, c'est la PRESENCE d'un vide immediat en dessous (teste ici),
    pas l'exclusivite avec is_overhang.
    """
    params = OverhangParams(
        overhang_angle_threshold_deg=45.0,
        bridge_detection_ray_density=6,
        bridge_support_ratio_threshold=0.5,
        bridge_max_gap_distance=1.0,  # cavite vide sur 4 unites : tout gap > 1 indique une absence de support immediat
    )
    result = overhang_run(hollow_box_mesh_data, params)

    assert np.any(result.is_bridge), (
        "Le plafond de la cavite interne devrait etre detecte comme pont, "
        "mais aucune face n'est marquee is_bridge"
    )

    # Specifiquement, les faces du plafond de la cavite interne (proches du
    # centre en Z, pas le fond exterieur de la boite qui repose sur le sol)
    # doivent etre marquees is_bridge.
    face_centers = hollow_box_mesh_data.vertices[hollow_box_mesh_data.faces].mean(axis=1)
    is_near_center_z = np.abs(face_centers[:, 2]) < 3.0  # la cavite est a mi-hauteur, le fond exterieur est a z=-5
    ceiling_candidates = result.is_bridge & is_near_center_z
    assert np.any(ceiling_candidates), (
        "Le plafond de la cavite interne (proche du centre en Z) devrait "
        "specifiquement etre marque is_bridge"
    )

    # bridge_span_direction doit etre defini (non-NaN) pour les faces pont
    bridge_indices = np.where(result.is_bridge)[0]
    directions = result.bridge_span_direction[bridge_indices]
    assert not np.any(np.isnan(directions)), (
        "bridge_span_direction devrait etre calculee pour toutes les faces "
        "marquees is_bridge"
    )
    # les directions doivent etre normalisees (issues d'une PCA/SVD)
    norms = np.linalg.norm(directions, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-6), "Les directions de portee devraient etre des vecteurs unitaires"


# ---------------------------------------------------------------------------
# Cas 4 - Sphere : le score de risque doit croitre progressivement vers le pole bas
# ---------------------------------------------------------------------------

def test_sphere_risk_score_increases_toward_bottom(sphere_mesh_data):
    params = OverhangParams(overhang_angle_threshold_deg=45.0)
    result = overhang_run(sphere_mesh_data, params)

    face_centers = sphere_mesh_data.vertices[sphere_mesh_data.faces].mean(axis=1)
    z_coords = face_centers[:, 2]

    # correlation entre hauteur (z) et risk_score : plus on descend, plus le
    # risque devrait etre eleve -> correlation negative attendue
    correlation = np.corrcoef(z_coords, result.risk_score)[0, 1]
    assert correlation < -0.5, (
        f"Le score de risque devrait etre fortement anti-correle avec la "
        f"hauteur z sur une sphere (pole bas = risque max), "
        f"correlation obtenue: {correlation:.3f}"
    )

    # aucune coupure brutale : le score doit rester dans [0, 1]
    assert np.all(result.risk_score >= 0.0) and np.all(result.risk_score <= 1.0)

    # le point le plus bas de la sphere doit avoir un risque proche de 1
    bottom_face_idx = np.argmin(z_coords)
    assert result.risk_score[bottom_face_idx] > 0.9, (
        "La face au pole inferieur de la sphere devrait avoir un score de "
        "risque proche de 1 (pointe droit vers le bas)"
    )


# ---------------------------------------------------------------------------
# Cas 5 - Sensibilite au seuil : changer le seuil doit changer le compte de faces
# ---------------------------------------------------------------------------

def test_overhang_count_increases_with_lower_threshold(sphere_mesh_data):
    strict_params = OverhangParams(overhang_angle_threshold_deg=20.0)
    lenient_params = OverhangParams(overhang_angle_threshold_deg=70.0)

    strict_result = overhang_run(sphere_mesh_data, strict_params)
    lenient_result = overhang_run(sphere_mesh_data, lenient_params)

    strict_count = strict_result.is_overhang.sum()
    lenient_count = lenient_result.is_overhang.sum()

    assert lenient_count > strict_count, (
        f"Un seuil plus permissif (70 deg) devrait marquer plus de faces en "
        f"surplomb qu'un seuil strict (20 deg). Obtenu: strict={strict_count}, "
        f"lenient={lenient_count}"
    )


# ---------------------------------------------------------------------------
# Test de coherence supplementaire : les shapes des tableaux de sortie
# ---------------------------------------------------------------------------

def test_output_shapes_match_face_count(flat_box_mesh_data):
    params = OverhangParams()
    result = overhang_run(flat_box_mesh_data, params)

    n = flat_box_mesh_data.face_count
    assert result.risk_score.shape == (n,)
    assert result.is_overhang.shape == (n,)
    assert result.is_bridge.shape == (n,)
    assert result.bridge_span_direction.shape == (n, 3)
