"""Suite découpée de `test_step1_ingest.py` (R-030 point 1, limite de 800 lignes par fichier) : tests déplacés tels quels."""

from __future__ import annotations

import os
import tempfile
import numpy as np
import pytest
import trimesh
from pipeline.step1_ingest import _clean_and_repair, run as ingest_run
from test_step1_ingest import (  # noqa: F401
    _l_shape_mesh,
    _tet_volumes,
)


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
