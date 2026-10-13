"""
tests/test_cpp_equivalence.py

ORACLE D'EQUIVALENCE du portage C++ : chaque fonction portee doit produire
exactement ce que produisait sa version Python. C'est le garde-fou central de
toute la migration — le Python reste la reference tant que le C++ n'a pas
prouve qu'il la reproduit.

Ces tests se SAUTENT proprement si le module natif n'est pas compile : lancer
la suite ne doit jamais exiger un compilateur C++ (voir core/_cpp.py).

Pourquoi comparer des ENSEMBLES et non des listes : `_PrintedGeometryIndex.
candidates_near` rend un sur-ensemble de candidats dont l'ordre est, cote
Python, celui d'iteration d'un `set` — donc deja non significatif. Les cinq
sites d'appel ont ete relus un par un et sont tous insensibles a l'ordre (ils
cherchent un minimum de distance ou l'existence d'un support). La comparaison
porte donc sur le contenu, jamais sur l'ordre.
"""

from __future__ import annotations

import numpy as np
import pytest

from core import _cpp
from core.mesh_types import AABB
from pipeline.step9_collision_sim import _PrintedGeometryIndex

pytestmark = pytest.mark.skipif(
    not _cpp.available(),
    reason=f"noyau natif apex_core non compile ({_cpp.load_error()})",
)

NOZZLE_DIAMETER = 0.4


@pytest.fixture
def restore_switch():
    """Restaure l'interrupteur apres le test - un reglage ne doit jamais fuiter."""
    previous = _cpp.set_enabled(None)
    yield
    _cpp.set_enabled(previous)


def _random_segments(count: int, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(seed)
    starts = rng.uniform(-20.0, 20.0, size=(count, 3))
    # Cordons courts, comme de vrais segments d'extrusion consecutifs.
    ends = starts + rng.uniform(-1.5, 1.5, size=(count, 3))
    return [(starts[i], ends[i]) for i in range(count)]


def _build(segments, *, native: bool, mode: int | None = None) -> _PrintedGeometryIndex:
    previous = _cpp.set_enabled(native)
    try:
        index = _PrintedGeometryIndex(NOZZLE_DIAMETER, cell_size_multiplier=4.0)
        if native and mode is not None:
            index._native.mode = mode
        for start, end in segments:
            index.append_segment(start, end)
        return index
    finally:
        _cpp.set_enabled(previous)


def _sorted_endpoints(index: _PrintedGeometryIndex, aabb: AABB, margin: float) -> np.ndarray:
    """Candidats sous forme canonique (C,6), tries — comparables sans ordre."""
    candidates = index.candidates_near(aabb, margin)
    if not candidates:
        return np.zeros((0, 6), dtype=float)
    rows = np.array(
        [np.concatenate([segment.start, segment.end]) for segment in candidates], dtype=float
    )
    order = np.lexsort(rows.T[::-1])
    return rows[order]


def _query_boxes(count: int, seed: int) -> list[AABB]:
    rng = np.random.default_rng(seed)
    centers = rng.uniform(-22.0, 22.0, size=(count, 3))
    extents = rng.uniform(0.2, 6.0, size=(count, 3))
    return [
        AABB(min_corner=centers[i] - extents[i], max_corner=centers[i] + extents[i])
        for i in range(count)
    ]


def test_selftest_roundtrip_is_bit_exact():
    """La chaine numpy -> C++ -> numpy ne doit alterer aucun bit."""
    values = np.array([0.0, -1.5, 3.25, 1e-17, -2.0e18, np.pi])
    assert np.array_equal(_cpp.module().selftest_roundtrip(values), values)


def test_candidates_identical_to_python(restore_switch):
    """Coeur de l'oracle : memes candidats que le Python, requete par requete."""
    segments = _random_segments(400, seed=17)
    python_index = _build(segments, native=False)
    native_index = _build(segments, native=True, mode=_cpp.GRID_MODE_STAMPED)

    compared_non_empty = 0
    for aabb in _query_boxes(150, seed=23):
        for margin in (0.0, 0.75):
            expected = _sorted_endpoints(python_index, aabb, margin)
            actual = _sorted_endpoints(native_index, aabb, margin)
            assert np.array_equal(expected, actual)
            if len(expected):
                compared_non_empty += 1

    # Sans cette garde, un bug rendant TOUT vide des deux cotes passerait
    # inapercu : le test comparerait des ensembles vides et serait vert.
    assert compared_non_empty > 50


def test_both_grid_modes_agree(restore_switch):
    """
    Le mode `faithful` (transcription du `set` Python) et le mode `stamped`
    (dedoublonnage par estampilles) doivent rendre le meme ensemble.

    Les deux coexistent pour mesurer separement le gain du LANGAGE et celui de
    l'ALGORITHME ; s'ils divergeaient, cette mesure ne voudrait rien dire.
    """
    segments = _random_segments(300, seed=5)
    faithful = _build(segments, native=True, mode=_cpp.GRID_MODE_FAITHFUL)
    stamped = _build(segments, native=True, mode=_cpp.GRID_MODE_STAMPED)

    for aabb in _query_boxes(80, seed=11):
        assert np.array_equal(
            _sorted_endpoints(faithful, aabb, 0.0), _sorted_endpoints(stamped, aabb, 0.0)
        )


@pytest.mark.parametrize("native", [False, True])
def test_candidates_endpoints_matches_candidates_near(native, restore_switch):
    """
    `candidates_endpoints` est le chemin rapide des trois consommateurs qui
    reconstruisaient ces tableaux a la main. Il doit rendre exactement le meme
    contenu que `candidates_near`, DANS LES DEUX implementations - sinon le
    chemin rapide serait une seconde source de verite.
    """
    segments = _random_segments(250, seed=31)
    index = _build(segments, native=native)

    for aabb in _query_boxes(60, seed=37):
        starts, ends = index.candidates_endpoints(aabb, margin=0.0)
        assert starts.shape == ends.shape
        assert starts.ndim == 2 and starts.shape[1] == 3

        rows = (
            np.concatenate([starts, ends], axis=1)
            if len(starts)
            else np.zeros((0, 6), dtype=float)
        )
        order = np.lexsort(rows.T[::-1]) if len(rows) else np.array([], dtype=int)
        assert np.array_equal(rows[order], _sorted_endpoints(index, aabb, 0.0))


@pytest.mark.parametrize("native", [False, True])
def test_empty_index_yields_empty_results(native, restore_switch):
    """Un index vide ne doit ni lever, ni rendre une forme inattendue."""
    index = _build([], native=native)
    aabb = AABB(min_corner=np.array([-1.0, -1.0, -1.0]), max_corner=np.array([1.0, 1.0, 1.0]))

    assert index.candidates_near(aabb, margin=0.0) == []
    starts, ends = index.candidates_endpoints(aabb, margin=0.0)
    assert starts.shape == (0, 3)
    assert ends.shape == (0, 3)


# ---------------------------------------------------------------------------
# Pas 2 : cinematique inverse N axes (step8)
# ---------------------------------------------------------------------------

def _resolve_both_ways(raw_points, axes, reference_direction, platform_origin, params):
    """Resout la meme trajectoire par les deux chemins, dans le meme processus."""
    from pipeline.step8_machine_toolpath import _resolve_all_axis_angles

    import pipeline.step8_machine_toolpath as step8

    results = {}
    for native in (False, True):
        previous = _cpp.set_enabled(native)
        try:
            results[native] = _resolve_all_axis_angles(
                raw_points, axes, reference_direction, platform_origin, params
            )
        finally:
            _cpp.set_enabled(previous)
        if native:
            # SANS cette verification, un repli silencieux du natif ferait
            # comparer Python a lui-meme et le test serait vert sans avoir
            # rien teste. C'est exactement le trou qui a masque, un temps,
            # le fait que le chemin natif n'etait pas emprunte du tout.
            assert step8._LAST_NATIVE_FALLBACK_COUNT is not None, (
                "le chemin natif n'a pas ete emprunte du tout"
            )
    return results[False], results[True]


def _assert_ik_equivalent(py_result, cpp_result, max_error_deg=5.0):
    """
    Equivalence de la resolution IK, avec DEUX niveaux d'exigence distincts.

    **Angles : 1e-9 degre.** C'est la sortie qui pilote la machine ; elle doit
    etre identique a la precision de calcul pres. La tolerance n'est pas nulle
    parce que le bracketing Python evalue ses 64 echantillons par un lot numpy
    la ou le C++ boucle en scalaire : memes operations, ordonnancement
    different, donc dernier bit potentiellement different.

    **Residuelles : 1e-6 radian**, tolerance BEAUCOUP plus lache et il faut
    dire pourquoi plutot que de la masquer. La residuelle vaut
    `acos(dot(achieved, target))`, et l'IK etant exacte ce produit scalaire
    vaut ~1-1e-15. Or `acos` est mal conditionnee en 1 : sa derivee y diverge
    (~1e7 pour cet ecart), si bien qu'un ecart d'UN ULP sur le produit
    scalaire produit ~1e-8 rad sur la residuelle. L'ecart observe est donc du
    bruit de conditionnement, pas une divergence de calcul.

    **Et surtout** : la seule chose que la residuelle DECIDE est l'emission
    d'un `Step8Diagnostic` au-dela de `max_unreachable_angular_error_deg`
    (5 degres par defaut, soit ~0,087 rad - quatre ordres de grandeur au-dessus
    du bruit). C'est cette decision qui est comparee exactement.
    """
    py_angles, py_residuals = py_result
    cpp_angles, cpp_residuals = cpp_result

    assert len(cpp_angles) == len(py_angles)
    np.testing.assert_allclose(np.array(cpp_angles), np.array(py_angles), rtol=0.0, atol=1e-9)
    np.testing.assert_allclose(
        np.array(cpp_residuals), np.array(py_residuals), rtol=0.0, atol=1e-6
    )

    py_flags = np.degrees(np.array(py_residuals)) > max_error_deg
    cpp_flags = np.degrees(np.array(cpp_residuals)) > max_error_deg
    assert np.array_equal(py_flags, cpp_flags), (
        "les deux implementations doivent emettre exactement les memes Step8Diagnostic"
    )


def _ik_fixture(axes, seed=101, count=250, with_origin=True):
    from pipeline.step8_machine_toolpath import Step8Params, _RawPoint

    rng = np.random.default_rng(seed)
    directions = rng.normal(size=(count, 3))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    positions = rng.uniform(-30.0, 30.0, size=(count, 3))
    # Quelques points EXACTEMENT sur l'axe de rotation : c'est la singularite
    # (atan2(0,0)) ou l'heritage d'angle du point precedent se declenche, donc
    # precisement le chemin ou une divergence Python/C++ passerait inapercue
    # sur des points quelconques.
    positions[::40] = np.array([0.0, 0.0, 5.0])

    raw_points = [
        _RawPoint(
            position=positions[i], desired_direction=directions[i], extrusion=1.0,
            feedrate=60.0, color_id=None, layer_id=0, segment_index=0, is_travel=False,
        )
        for i in range(count)
    ]
    origin = np.array([0.0, 0.0, 0.0]) if with_origin else None
    return raw_points, axes, np.array([0.0, 0.0, -1.0]), origin, Step8Params()


def _current_machine_axes():
    """Cinematique de la machine reelle : C (tour complet) + B [0,185]."""
    from core.mesh_types import RotationAxisSpec

    return [
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0)),
    ]


def test_ik_matches_python_on_real_machine_kinematics(restore_switch):
    """
    Coeur de l'oracle du pas 2, sur la cinematique reellement utilisee
    (C sur un tour complet + B dans [0,185]). Les niveaux de tolerance et
    leur justification sont dans `_assert_ik_equivalent`.
    """
    fixture = _ik_fixture(_current_machine_axes())
    py_result, cpp_result = _resolve_both_ways(*fixture)
    _assert_ik_equivalent(py_result, cpp_result)


def test_ik_resumes_after_grid_fallback(restore_switch):
    """
    LE test du mecanisme de reprise, et le plus important du pas 2.

    Une cible colineaire a l'axe externe (+Z, un depot vertical) met en echec
    la resolution analytique de la paire : le residu du bracketing y est une
    TANGENCE a zero, jamais un changement de signe. Le natif s'arrete alors et
    rend la main pour ce point.

    On en place plusieurs, dont le TOUT PREMIER - le cas reellement rencontre
    sur `death.stl`, ou le premier point d'une trajectoire a fond plat depose
    verticalement. Une version anterieure abandonnait le lot entier au premier
    repli : un seul point sur 297 994 annulait alors tout le gain, sans qu'aucun
    test ne s'en apercoive.

    Verifie donc deux choses : que des replis ont bien eu lieu (sans quoi le
    test n'exerce pas la reprise), et que le resultat reste identique au Python.
    """
    import pipeline.step8_machine_toolpath as step8

    raw_points, axes, reference, origin, params = _ik_fixture(
        _current_machine_axes(), seed=71, count=150
    )
    up = np.array([0.0, 0.0, 1.0])
    for index in (0, 37, 149):
        raw_points[index].desired_direction = up

    py_result, cpp_result = _resolve_both_ways(raw_points, axes, reference, origin, params)

    assert step8._LAST_NATIVE_FALLBACK_COUNT >= 3, (
        "les cibles +Z doivent declencher le repli grille ; sans repli, ce test "
        f"n'exerce pas la reprise (count={step8._LAST_NATIVE_FALLBACK_COUNT})"
    )
    _assert_ik_equivalent(py_result, cpp_result)


def test_ik_matches_python_for_one_axis(restore_switch):
    """Cas N=1 : chemin `_resolve_single_axis_analytic`, distinct de la paire."""
    from core.mesh_types import RotationAxisSpec

    axes = [RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0))]
    _assert_ik_equivalent(*_resolve_both_ways(*_ik_fixture(axes, seed=7)))


def test_ik_matches_python_for_locked_and_parallel_axes(restore_switch):
    """
    Deux degenerescences que la transcription doit reproduire : un axe
    VERROUILLE (lo == hi) et deux axes PARALLELES (un seul DOF effectif,
    plage combinee). Chacune emprunte une branche distincte de
    `_resolve_axis_pair_analytic` que les cas generiques n'atteignent jamais.
    """
    from core.mesh_types import RotationAxisSpec

    locked = [
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
        RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(45.0, 45.0)),
    ]
    parallel = [
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 180.0)),
        RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 180.0)),
    ]
    for axes in (locked, parallel):
        _assert_ik_equivalent(*_resolve_both_ways(*_ik_fixture(axes, seed=13, count=120)))


def test_ik_matches_python_without_platform_origin(restore_switch):
    """`platform_origin is None` : aucune detection de singularite ne doit courir."""
    fixture = _ik_fixture(_current_machine_axes(), seed=23, count=120, with_origin=False)
    _assert_ik_equivalent(*_resolve_both_ways(*fixture))


# ---------------------------------------------------------------------------
# Pas 3 : test de support fusionne (requete spatiale + geometrie)
# ---------------------------------------------------------------------------

def _support_directions(count: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    directions = rng.normal(size=(count, 3))
    return directions / np.linalg.norm(directions, axis=1)[:, None]


@pytest.mark.parametrize("direction_count", [1, 5])
def test_has_support_matches_python(direction_count, restore_switch):
    """
    Oracle du pas 3. Le chemin natif FUSIONNE la requete spatiale et le test
    geometrique, avec ARRET ANTICIPE au premier appui trouve ; le chemin
    Python materialise les candidats puis vectorise sur tous. Les deux doivent
    rendre le meme booleen — l'arret anticipe ne change rien, puisque
    `np.any(...)` est vrai des qu'un seul candidat convient.
    """
    segments = _random_segments(500, seed=53)
    python_index = _build(segments, native=False)
    native_index = _build(segments, native=True, mode=_cpp.GRID_MODE_STAMPED)

    directions = _support_directions(direction_count, seed=59)
    rng = np.random.default_rng(61)
    positions = rng.uniform(-20.0, 20.0, size=(300, 3))

    true_count = 0
    for position in positions:
        for radius, height in ((0.6, 0.6), (2.5, 2.5)):
            aabb = AABB(
                min_corner=position - max(radius, height),
                max_corner=position + max(radius, height),
            )
            expected = python_index.has_support(aabb, position, directions, radius, height)
            actual = native_index.has_support(aabb, position, directions, radius, height)
            assert expected == actual
            true_count += int(expected)

    # Sans cette garde, un bug rendant TOUT False des deux cotes passerait
    # inapercu : le test comparerait des faux a des faux et serait vert.
    assert true_count > 20, f"le jeu de test n'exerce pas assez de cas positifs ({true_count})"


def test_has_support_respects_the_strictly_below_condition(restore_switch):
    """
    Le critere exige un appui STRICTEMENT en dessous (`delta_up > 1e-9`) : une
    capsule exactement a la hauteur du point ne le supporte pas. Cas limite
    facile a perdre en transcrivant, et invisible sur des positions aleatoires.
    """
    segments = [(np.array([-1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]))]
    position = np.array([0.0, 0.0, 0.0])  # exactement SUR la capsule
    directions = np.array([[0.0, 0.0, 1.0]])
    aabb = AABB(min_corner=position - 2.0, max_corner=position + 2.0)

    for native in (False, True):
        index = _build(segments, native=native)
        assert index.has_support(aabb, position, directions, 1.0, 1.0) is False

    # Décalé vers le haut : la capsule passe alors strictement en dessous.
    above = np.array([0.0, 0.0, 0.5])
    for native in (False, True):
        index = _build(segments, native=native)
        assert index.has_support(aabb, above, directions, 1.0, 1.0) is True


# ---------------------------------------------------------------------------
# Pas 4 : tables topologiques tetraedriques
# ---------------------------------------------------------------------------

def _random_tetrahedra(tet_count: int, node_count: int, seed: int) -> np.ndarray:
    """
    Tetraedres synthetiques partageant beaucoup de faces — indispensable pour
    que le test exerce reellement les faces INTERNES (deux proprietaires) et
    pas seulement des faces frontiere isolees.
    """
    rng = np.random.default_rng(seed)
    tets = rng.integers(0, node_count, size=(tet_count, 4), dtype=np.int64)
    # Elimine les tetraedres a sommets repetes : le Python les traiterait, mais
    # ils n'ont pas de sens geometrique et brouillent la comparaison.
    valid = np.array([len(set(t.tolist())) == 4 for t in tets])
    return tets[valid]


def test_tet_face_topology_matches_python(restore_switch):
    """
    L'ORDRE des faces frontiere doit etre celui d'INSERTION, pas l'ordre trie.

    Ce n'est pas cosmetique : `_map_boundary_faces_to_surface` interroge
    ensuite le kdtree avec ces faces dans cet ordre, et l'ordre du resultat
    s'y propage. Une permutation produirait des tableaux decales en aval, sans
    qu'aucune exception ne soit levee.
    """
    from pipeline.step4_field_solver import _tet_face_table

    class _FakeTetMesh:
        def __init__(self, tetrahedra):
            self.tetrahedra = tetrahedra
            self.tet_count = len(tetrahedra)

    tetrahedra = _random_tetrahedra(4000, 220, seed=91)
    tet_mesh = _FakeTetMesh(tetrahedra)

    face_to_tets = _tet_face_table(tet_mesh)
    expected_nodes = [key for key, owners in face_to_tets.items() if len(owners) == 1]
    expected_tets = [owners[0] for owners in face_to_tets.values() if len(owners) == 1]
    expected_pairs = sorted(
        tuple(sorted(owners)) for owners in face_to_tets.values() if len(owners) == 2
    )

    nodes, owners, pairs = _cpp.module().tet_face_topology(tetrahedra)

    assert [tuple(row) for row in nodes.tolist()] == expected_nodes
    assert owners.tolist() == expected_tets
    assert sorted(tuple(sorted(p)) for p in pairs.tolist()) == expected_pairs

    # Gardes : sans elles, un jeu degenere (aucune face interne, ou aucune
    # frontiere) rendrait le test vert sans rien exercer.
    assert len(expected_nodes) > 50
    assert len(expected_pairs) > 50


def test_tet_unique_edges_matches_python(restore_switch):
    """Memes aretes uniques que la boucle a `set` de `_build_tet_adjacency`."""
    tetrahedra = _random_tetrahedra(3000, 180, seed=97)

    expected = set()
    for tet in tetrahedra:
        for i in range(4):
            for j in range(i + 1, 4):
                a, b = int(tet[i]), int(tet[j])
                expected.add((min(a, b), max(a, b)))

    edges = _cpp.module().tet_unique_edges(tetrahedra)
    assert {tuple(e) for e in edges.tolist()} == expected
    assert len(expected) > 100


def test_tet_adjacency_matrix_identical_both_paths(restore_switch):
    """
    Bout en bout : la matrice creuse produite doit etre identique, ce qui
    verifie aussi que l'ordre des aretes est bien sans importance (csr trie
    ses indices a la construction).
    """
    from pipeline.step1_ingest import _build_tet_adjacency

    tetrahedra = _random_tetrahedra(2500, 150, seed=101)
    node_count = 150

    previous = _cpp.set_enabled(False)
    try:
        expected = _build_tet_adjacency(tetrahedra, node_count)
    finally:
        _cpp.set_enabled(previous)

    previous = _cpp.set_enabled(True)
    try:
        actual = _build_tet_adjacency(tetrahedra, node_count)
    finally:
        _cpp.set_enabled(previous)

    assert (expected != actual).nnz == 0
    assert expected.nnz > 100


@pytest.mark.parametrize("native", [False, True])
def test_append_segments_batch_matches_one_by_one(native, restore_switch):
    """
    L'insertion en lot doit rendre exactement le meme index que l'insertion
    capsule par capsule — mêmes candidats pour les mêmes requetes.
    """
    segments = _random_segments(300, seed=107)
    starts = np.array([s for s, _ in segments])
    ends = np.array([e for _, e in segments])

    one_by_one = _build(segments, native=native)

    previous = _cpp.set_enabled(native)
    try:
        batched = _PrintedGeometryIndex(NOZZLE_DIAMETER, cell_size_multiplier=4.0)
        batched.append_segments(starts, ends)
    finally:
        _cpp.set_enabled(previous)

    for aabb in _query_boxes(60, seed=109):
        assert np.array_equal(
            _sorted_endpoints(one_by_one, aabb, 0.0), _sorted_endpoints(batched, aabb, 0.0)
        )


def test_native_is_not_used_when_switch_is_off(restore_switch):
    """
    Non-regression du defaut : interrupteur eteint => aucun objet natif, donc
    le chemin Python historique est bien celui emprunte. C'est ce qui garantit
    qu'un environnement sans compilateur se comporte exactement comme avant.
    """
    index = _build(_random_segments(10, seed=3), native=False)
    assert index._native is None
    assert index._grid, "le chemin Python doit avoir rempli sa propre grille"


# ---------------------------------------------------------------------------
# Pas 5 : arbre k-d spatial (core/spatial.py, cpp/src/core/kdtree.cpp)
# ---------------------------------------------------------------------------
#
# Portage ARCHITECTURAL, pas une optimisation : `scipy.spatial.cKDTree` est
# deja natif. L'oracle verifie donc surtout la NON-DIVERGENCE - memes indices
# de plus proche voisin, memes distances (au bit), memes ensembles de rayon -
# sur des donnees non degenerees (aucune egalite exacte de distance, cas ou
# scipy ne specifie pas quel indice il rend). Voir cpp/include/apex/kdtree.hpp.

from core.spatial import KDTree as _WrappedKDTree  # noqa: E402


def _random_cloud(count: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).uniform(-40.0, 40.0, size=(count, 3))


def test_kdtree_query_nn_matches_scipy_bit_exact(restore_switch):
    """Indices identiques et distances IDENTIQUES AU BIT face a cKDTree."""
    import scipy.spatial

    points = _random_cloud(6000, seed=41)
    queries = np.random.default_rng(43).uniform(-45.0, 45.0, size=(2500, 3))

    reference = scipy.spatial.cKDTree(points)
    ref_dist, ref_idx = reference.query(queries)

    _cpp.set_enabled(True)
    native = _cpp.module().KDTree(np.ascontiguousarray(points))
    nat_dist, nat_idx = native.query_nn(np.ascontiguousarray(queries))

    assert np.array_equal(nat_idx, ref_idx)
    assert np.array_equal(nat_dist, ref_dist)  # bit-exact, pas allclose

    # Garde anti-faux-vert : le jeu doit exercer beaucoup de voisins distincts,
    # pas retomber toujours sur le meme point.
    assert len(np.unique(ref_idx)) > 200


def test_kdtree_query_radius_matches_scipy(restore_switch):
    """Memes indices dans un rayon que `cKDTree.query_ball_point` (comme ensembles)."""
    import scipy.spatial

    points = _random_cloud(4000, seed=51)
    reference = scipy.spatial.cKDTree(points)

    _cpp.set_enabled(True)
    native = _cpp.module().KDTree(np.ascontiguousarray(points))

    centers = np.random.default_rng(53).uniform(-20.0, 20.0, size=(40, 3))
    non_empty = 0
    for center in centers:
        for radius in (3.0, 9.0, 22.0):
            expected = set(reference.query_ball_point(center, radius))
            actual = set(native.query_radius(np.ascontiguousarray(center), radius))
            assert expected == actual
            non_empty += int(bool(expected))
    assert non_empty > 30


def test_wrapped_kdtree_dispatches_and_agrees(restore_switch):
    """
    `core.spatial.KDTree` : `._native is None` interrupteur eteint, un objet
    natif allume, et les deux chemins rendent le meme `.query()`.
    """
    points = _random_cloud(3000, seed=61)
    queries = np.random.default_rng(63).uniform(-45.0, 45.0, size=(900, 3))

    _cpp.set_enabled(False)
    python_tree = _WrappedKDTree(points)
    assert python_tree._native is None
    py_dist, py_idx = python_tree.query(queries)

    _cpp.set_enabled(True)
    native_tree = _WrappedKDTree(points)
    assert native_tree._native is not None
    nat_dist, nat_idx = native_tree.query(queries)

    assert np.array_equal(py_idx, nat_idx)
    assert np.array_equal(py_dist, nat_dist)

    # `.data` preserve, et une requete 1-D rend bien des scalaires des deux cotes.
    assert np.array_equal(python_tree.data, native_tree.data)
    one_py = python_tree.query(queries[0])
    one_nat = native_tree.query(queries[0])
    assert isinstance(one_nat[0], float) and isinstance(one_nat[1], int)
    assert one_py[1] == one_nat[1]
    assert one_py[0] == one_nat[0]


def test_wrapped_kdtree_query_ball_point_agrees(restore_switch):
    """`query_ball_point` du wrapper : meme ensemble dans les deux positions."""
    points = _random_cloud(2500, seed=71)

    _cpp.set_enabled(False)
    python_tree = _WrappedKDTree(points)
    _cpp.set_enabled(True)
    native_tree = _WrappedKDTree(points)

    non_empty = 0
    for center in np.random.default_rng(73).uniform(-15.0, 15.0, size=(30, 3)):
        for radius in (4.0, 12.0):
            expected = set(python_tree.query_ball_point(center, radius))
            actual = set(native_tree.query_ball_point(center, radius))
            assert expected == actual
            non_empty += int(bool(expected))
    assert non_empty > 20


def test_wrapped_kdtree_native_absent_when_switch_off(restore_switch):
    """Non-regression : interrupteur eteint => strictement le chemin scipy."""
    import scipy.spatial

    _cpp.set_enabled(False)
    tree = _WrappedKDTree(_random_cloud(50, seed=3))
    assert tree._native is None
    assert isinstance(tree._scipy, scipy.spatial.cKDTree)


# ---------------------------------------------------------------------------
# Pas 6 : BVH de triangles (core/spatial.RayIntersector, cpp/src/core/bvh.cpp)
# ---------------------------------------------------------------------------
#
# Contrairement au k-d tree, la reference Python n'est PAS native ici
# (`pyembree` absent du venv) : le gain est reel, mesure ×32 sur death.stl et
# ×139 sur death_hi-res.stl. L'oracle porte sur l'ENSEMBLE des couples
# (index_ray, location) dedoublonnes — les deux appelants du pipeline
# (step2, step7) ignorent `index_tri` et reduisent par `.min()`, donc l'ordre
# n'est pas significatif. Voir cpp/include/apex/bvh.hpp.

import trimesh  # noqa: E402
import trimesh.ray.ray_triangle  # noqa: E402

from core.spatial import RayIntersector as _WrappedRayIntersector  # noqa: E402


def _hit_rows(locations: np.ndarray, index_ray: np.ndarray) -> np.ndarray:
    """Impacts sous forme canonique (H,4) triee — comparable sans ordre."""
    if len(locations) == 0:
        return np.zeros((0, 4), dtype=float)
    rows = np.column_stack([np.asarray(locations, dtype=float),
                            np.asarray(index_ray, dtype=float)])
    return rows[np.lexsort(rows.T[::-1])]


def _ray_test_mesh() -> trimesh.Trimesh:
    """Maillage synthetique non trivial : une sphere moins une boite, donc des
    faces internes ET externes, et de vraies aretes partagees (le cas que le
    dedoublonnage de `unique_rows` doit collapser)."""
    sphere = trimesh.creation.icosphere(subdivisions=2, radius=10.0)
    box = trimesh.creation.box(extents=(6.0, 6.0, 30.0))
    return sphere.difference(box, engine="manifold")


@pytest.mark.parametrize(
    "direction_mode", ["vertical", "random_unit", "unnormalized"]
)
def test_bvh_intersects_location_matches_trimesh(direction_mode, restore_switch):
    """
    Memes impacts que `RayMeshIntersector.intersects_location`, AU BIT.

    Les trois modes de direction ne sont pas cosmetiques : `vertical` exerce
    le cas `1/0 = inf` du test de tranches, et `unnormalized` verifie que le
    filtre "forward" utilise bien le produit scalaire avec la direction
    BRUTE (comme `util.diagonal_dot`) et non une direction renormalisee.
    """
    mesh = _ray_test_mesh()
    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)

    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )

    rng = np.random.default_rng(1301)
    lo, hi = mesh.bounds
    count = 800
    origins = rng.uniform(lo, hi, size=(count, 3))
    if direction_mode == "vertical":
        directions = np.tile([0.0, 0.0, -1.0], (count, 1))
    elif direction_mode == "random_unit":
        directions = rng.normal(size=(count, 3))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    else:
        directions = rng.normal(size=(count, 3)) * rng.uniform(0.1, 5.0, size=(count, 1))

    ref_loc, ref_ray, _ref_tri = reference.intersects_location(origins, directions)
    nat_loc, nat_ray, _nat_tri = native.intersects_location(
        np.ascontiguousarray(origins), np.ascontiguousarray(directions)
    )

    assert np.array_equal(_hit_rows(ref_loc, ref_ray), _hit_rows(nat_loc, nat_ray))

    # Garde anti-faux-vert : sans elle, deux resultats VIDES seraient egaux.
    assert len(ref_loc) > 100, f"jeu de test trop pauvre ({len(ref_loc)} impacts)"


def test_bvh_handles_misses_and_empty_input(restore_switch):
    """Rayons qui ratent tout, et lot de rayons vide : ni exception, ni forme
    inattendue, des deux cotes."""
    mesh = _ray_test_mesh()
    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )

    far = np.tile(mesh.bounds[1] + 1000.0, (10, 1))
    away = np.tile([0.0, 0.0, 1.0], (10, 1))
    for loc, ray, _tri in (
        reference.intersects_location(far, away),
        native.intersects_location(np.ascontiguousarray(far), np.ascontiguousarray(away)),
    ):
        assert len(loc) == 0 and len(ray) == 0

    empty = np.zeros((0, 3))
    loc, ray, _tri = native.intersects_location(empty, empty)
    assert len(loc) == 0 and len(ray) == 0


def test_bvh_deduplicates_edge_hits_like_trimesh(restore_switch):
    """
    Un impact PILE sur l'arete partagee par deux triangles ne compte qu'UNE
    fois — c'est le role de `grouping.unique_rows` cote trimesh, et la raison
    pour laquelle un comptage par triangle differe d'un comptage par point.

    Cas construit explicitement : deux triangles coplanaires partageant
    l'arete x=0, et un rayon vertical qui tombe exactement dessus. Un
    dedoublonnage manquant donnerait 2 impacts au lieu d'1, et aucune donnee
    aleatoire ne tomberait sur ce cas.
    """
    vertices = np.array(
        [[-1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 1.0, 0.0],
         [1.0, -1.0, 0.0], [1.0, 1.0, 0.0]]
    )
    faces = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]])
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    origin = np.array([[0.0, 0.0, 5.0]])          # pile au-dessus de l'arete x=0
    direction = np.array([[0.0, 0.0, -1.0]])

    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    ref_loc, ref_ray, _ = reference.intersects_location(origin, direction)

    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )
    nat_loc, nat_ray, _ = native.intersects_location(origin, direction)

    assert np.array_equal(_hit_rows(ref_loc, ref_ray), _hit_rows(nat_loc, nat_ray))
    # Le cas doit VRAIMENT toucher l'arete, sinon le test ne verifie rien.
    assert len(ref_loc) == 1, f"le rayon devrait donner 1 impact dedoublonne, pas {len(ref_loc)}"


def test_wrapped_ray_intersector_dispatches_and_agrees(restore_switch):
    """`core.spatial.RayIntersector` : dispatch correct et memes impacts."""
    mesh = _ray_test_mesh()
    rng = np.random.default_rng(1307)
    lo, hi = mesh.bounds
    origins = rng.uniform(lo, hi, size=(400, 3))
    directions = np.tile([0.0, 0.0, -1.0], (400, 1))

    _cpp.set_enabled(False)
    python_side = _WrappedRayIntersector(mesh)
    assert python_side._native is None
    py_loc, py_ray, _ = python_side.intersects_location(origins, directions)

    _cpp.set_enabled(True)
    native_side = _WrappedRayIntersector(mesh)
    assert native_side._native is not None
    nat_loc, nat_ray, _ = native_side.intersects_location(origins, directions)

    assert np.array_equal(_hit_rows(py_loc, py_ray), _hit_rows(nat_loc, nat_ray))
    assert len(py_loc) > 50


def test_wrapped_ray_intersector_rejects_unsupported_options(restore_switch):
    """
    Une option passee au chemin natif doit LEVER, jamais etre ignoree en
    silence. Aucun appelant du pipeline n'en passe aujourd'hui ; si l'un s'y
    met, l'ecart de comportement doit etre bruyant.
    """
    _cpp.set_enabled(True)
    intersector = _WrappedRayIntersector(_ray_test_mesh())
    with pytest.raises(TypeError, match="non supportees"):
        intersector.intersects_location(
            np.zeros((1, 3)), np.tile([0.0, 0.0, -1.0], (1, 1)), multiple_hits=False
        )


def test_bvh_excludes_the_originating_face(restore_switch):
    """
    NON-REGRESSION D'UN BUG REEL, trouve par la mesure sur le pipeline et non
    par des rayons aleatoires.

    `step2` part d'un point SUR une face, decale de 1e-6 vers le bas, et tire
    vers le bas. L'impact avec sa PROPRE face tombe a `dot = -1e-6`, soit
    exactement le seuil du filtre "vers l'avant" — et le calcul flottant reel
    donne -9.999999994736442e-07, strictement SUPERIEUR au seuil. La narrow
    phase de trimesh l'accepterait donc ; ce qui l'ecarte pour une face
    HORIZONTALE, c'est son BROAD PHASE, qui clippe `t` a
    `buffer_dist = 1e-5` : la boite du rayon s'arrete alors juste sous la
    face, dont l'AABB plate ne la recoupe plus.

    Une premiere version du BVH descendait par un test rayon/AABB classique,
    plus permissif : elle rendait cet impact parasite a distance ~1e-6, et
    `step2` tombait de 49 ponts detectes a 0 sur death.stl. Aucun des 3 000
    rayons aleatoires des tests ci-dessus ne l'avait vu, parce qu'aucun ne
    partait exactement d'une face.

    La fixture est une boite creuse : son plafond de cavite est plat et
    horizontal, exactement la geometrie ou le broad phase tranche.
    """
    outer = trimesh.creation.box(extents=(10.0, 10.0, 10.0))
    inner = trimesh.creation.box(extents=(4.0, 4.0, 4.0))
    mesh = outer.difference(inner, engine="manifold")

    # Point sur le plafond de la cavite (z = +2), decale comme le fait step2.
    origins = np.array([[0.0, 0.0, 2.0 - 1e-6]])
    directions = np.array([[0.0, 0.0, -1.0]])

    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    ref_loc, ref_ray, _ = reference.intersects_location(origins, directions)

    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )
    nat_loc, nat_ray, _ = native.intersects_location(
        np.ascontiguousarray(origins), np.ascontiguousarray(directions)
    )

    assert np.array_equal(_hit_rows(ref_loc, ref_ray), _hit_rows(nat_loc, nat_ray))

    # L'assertion qui verrouille le bug : aucun impact a distance quasi nulle.
    # Sans elle, deux implementations egalement fautives passeraient.
    for loc, ray in ((ref_loc, ref_ray), (nat_loc, nat_ray)):
        distances = np.linalg.norm(loc - origins[np.asarray(ray, dtype=np.int64)], axis=1)
        assert (distances > 1e-4).all(), "impact parasite sur la face d'origine"
    # Garde : le rayon doit VRAIMENT traverser la cavite puis le fond, sinon
    # le test ne verifie rien (il passerait sur un resultat vide).
    assert len(ref_loc) == 2, f"attendu 2 impacts (fond de cavite + fond de boite), vu {len(ref_loc)}"


def test_bvh_matches_trimesh_for_rays_starting_on_faces(restore_switch):
    """
    Cas general du precedent : des centaines d'origines posees SUR les faces,
    decalees comme step2. Ici trimesh accepte parfois l'auto-impact (face non
    plate, dont l'AABB recoupe encore la boite du rayon) — l'oracle exige
    seulement que le natif fasse EXACTEMENT le meme choix, pas qu'il applique
    une regle plus propre.
    """
    mesh = _ray_test_mesh()
    origins = mesh.triangles.mean(axis=1)[:300] + np.array([0.0, 0.0, -1e-6])
    directions = np.tile([0.0, 0.0, -1.0], (len(origins), 1))

    reference = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    ref_loc, ref_ray, _ = reference.intersects_location(origins, directions)

    _cpp.set_enabled(True)
    native = _cpp.module().TriangleBVH(
        np.ascontiguousarray(mesh.triangles), np.ascontiguousarray(mesh.face_normals)
    )
    nat_loc, nat_ray, _ = native.intersects_location(
        np.ascontiguousarray(origins), np.ascontiguousarray(directions)
    )

    assert np.array_equal(_hit_rows(ref_loc, ref_ray), _hit_rows(nat_loc, nat_ray))
    assert len(ref_loc) > 50, f"jeu de test trop pauvre ({len(ref_loc)} impacts)"


# ---------------------------------------------------------------------------
# Pas 7 : noyau numerique de step4 (assemblage FEM + solveur)
# ---------------------------------------------------------------------------
#
# PREMIER PORTAGE DU PROJET DONT L'EQUIVALENCE N'EST PAS BIT-A-BIT, et c'est
# une decision explicite, pas un accident. Deux causes, toutes deux
# irreductibles :
#   - `compute_shape_gradients` inverse une 3x3 par la formule des
#     cofacteurs, la ou numpy appelle LAPACK (LU a pivot partiel) : meme
#     inverse, autres arrondis. Ecart RELATIF mesure : 7,2e-14.
#   - le solveur : `scipy.sparse.linalg.cg` et le gradient conjugue d'Eigen
#     appliquent le meme algorithme avec le meme preconditionneur diagonal,
#     mais leurs produits scalaires n'accumulent pas dans le meme ordre.
#
# Les tolerances ci-dessous sont donc DECLAREES, jamais subies : chacune est
# justifiee par la mesure qui l'accompagne.

def _step4_inputs(resolution: float = 3.0):
    """Pipeline step1->step3 sur une forme en L, assez fine pour exercer
    plusieurs regions (donc plusieurs sources de melange)."""
    import tempfile
    import os

    import trimesh

    from pipeline.step1_ingest import run as ingest_run
    from pipeline.step2_overhang import OverhangParams, run as overhang_run
    from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run

    vertical = trimesh.creation.box(extents=[6, 6, 24])
    vertical.apply_translation([0, 0, 12])
    horizontal = trimesh.creation.box(extents=[24, 6, 6])
    horizontal.apply_translation([9, 0, 24])
    mesh = vertical.union(horizontal, engine="manifold")

    handle, path = tempfile.mkstemp(suffix=".stl")
    os.close(handle)
    try:
        mesh.export(path)
        mesh_data, tet_mesh, *_ = ingest_run(path, resolution=resolution)
    finally:
        os.unlink(path)

    overhang = overhang_run(mesh_data, OverhangParams())
    regions = segmentation_run(mesh_data, overhang, SegmentationParams(platform_z=0.0))
    return mesh_data, tet_mesh, overhang, regions


def test_shape_gradients_match_numpy_to_machine_precision(restore_switch):
    """
    Cofacteurs (C++) contre `np.linalg.inv` (LAPACK) : meme inverse 3x3 par
    deux algorithmes differents. La tolerance porte sur l'ecart RELATIF - les
    gradients atteignent des magnitudes de plusieurs milliers sur un maillage
    fin, et un seuil absolu y serait ininterpretable.
    """
    from pipeline.step4_field_solver import _compute_shape_gradients

    _mesh_data, tet_mesh, _overhang, _regions = _step4_inputs()

    _cpp.set_enabled(False)
    gradients_py, volumes_py = _compute_shape_gradients(tet_mesh)
    _cpp.set_enabled(True)
    gradients_native, volumes_native = _compute_shape_gradients(tet_mesh)

    scale = max(float(np.abs(gradients_py).max()), 1e-30)
    assert np.abs(gradients_py - gradients_native).max() / scale < 1e-12
    assert np.allclose(volumes_py, volumes_native, rtol=1e-12, atol=0.0)

    # Garde anti-faux-vert : un maillage vide comparerait deux tableaux vides.
    assert tet_mesh.tet_count > 100
    assert volumes_py.sum() > 0.0


def test_blend_tet_directions_matches_python(restore_switch):
    """
    Le chemin natif ECARTE par AABB les sources hors rayon, que le Python
    calcule puis annule par la fenetre de Hann. L'elagage est exact - une
    source ecartee a un poids strictement nul - donc les directions doivent
    coincider a la precision machine, pas seulement "a peu pres".
    """
    from pipeline.step4_field_solver import FieldSolverParams, _blend_tet_directions

    mesh_data, tet_mesh, _overhang, regions = _step4_inputs()
    params = FieldSolverParams(direction_blend_radius=8.0, direction_blend_power=2.0)

    _cpp.set_enabled(False)
    expected = _blend_tet_directions(tet_mesh, mesh_data, regions, params.up_vector, params)
    _cpp.set_enabled(True)
    actual = _blend_tet_directions(tet_mesh, mesh_data, regions, params.up_vector, params)

    assert expected.shape == actual.shape
    assert np.abs(expected - actual).max() < 1e-12

    # Les directions doivent rester unitaires des deux cotes.
    assert np.allclose(np.linalg.norm(actual, axis=1), 1.0, atol=1e-9)

    # Garde anti-faux-vert, sur ce qui compte VRAIMENT : le melange doit
    # produire des directions variees. Compter les regions ne suffirait pas -
    # le plateau est lui aussi une source, et une seule region melangee avec
    # lui exerce deja pleinement la ponderation.
    assert float(np.ptp(actual, axis=0).max()) > 1e-3, (
        "le melange rend une direction quasi constante : il n'est pas exerce"
    )


def test_blend_pruning_agrees_on_a_radius_that_excludes_almost_everything(restore_switch):
    """
    Cas limite de l'elagage : un rayon minuscule ecarte presque toutes les
    paires, si bien que la quasi-totalite des tetraedres empruntent le REPLI
    `nearest_source`. C'est precisement le chemin que la seconde passe sans
    elagage doit rattraper - un elagage naif y divergerait.
    """
    from pipeline.step4_field_solver import FieldSolverParams, _blend_tet_directions

    mesh_data, tet_mesh, _overhang, regions = _step4_inputs()
    params = FieldSolverParams(direction_blend_radius=1e-3, direction_blend_power=2.0)

    _cpp.set_enabled(False)
    expected = _blend_tet_directions(tet_mesh, mesh_data, regions, params.up_vector, params)
    _cpp.set_enabled(True)
    actual = _blend_tet_directions(tet_mesh, mesh_data, regions, params.up_vector, params)

    assert np.abs(expected - actual).max() < 1e-12


def test_eigen_solver_agrees_with_scipy_within_solver_tolerance(restore_switch):
    """
    Eigen et scipy resolvent le MEME systeme au MEME critere d'arret relatif,
    par le meme algorithme et le meme preconditionneur - mais pas dans le meme
    ordre d'accumulation. On verifie donc ce qui a un sens : les deux champs
    concordent a l'echelle de la tolerance du solveur, et surtout leurs
    RESIDUS sont tous deux sous cette tolerance. Comparer les vecteurs entre
    eux sans regarder le residu ne dirait pas lequel resout reellement.
    """
    import scipy.sparse as sp

    from pipeline.step4_field_solver import (
        FieldSolverParams,
        _solve_linear_system,
        build_static_system,
    )

    mesh_data, tet_mesh, overhang, regions = _step4_inputs()
    params = FieldSolverParams(direction_blend_radius=8.0)

    _cpp.set_enabled(False)
    system = build_static_system(tet_mesh, mesh_data, overhang, regions, params)
    A = system.A_structural.tocsr()
    A.sort_indices()
    b = system.b_structural

    reference = _solve_linear_system(A, b, params)

    _cpp.set_enabled(True)
    direct = A.shape[0] < params.direct_solve_node_threshold
    values, info = _cpp.module().solve_linear_system(
        np.ascontiguousarray(A.indptr, dtype=np.int32),
        np.ascontiguousarray(A.indices, dtype=np.int32),
        np.ascontiguousarray(A.data),
        np.ascontiguousarray(b),
        direct,
        params.convergence_tolerance,
        params.max_solver_iterations,
        None,
    )
    assert info["converged"], "le solveur natif doit converger sur ce systeme"

    norm_b = np.linalg.norm(b)
    assert norm_b > 0.0, "un second membre nul rendrait le test vide de sens"
    residual_native = np.linalg.norm(A @ values - b) / norm_b
    residual_reference = np.linalg.norm(A @ reference - b) / norm_b

    # Le critere qui compte : les deux resolvent reellement le systeme.
    assert residual_native <= max(params.convergence_tolerance * 10.0, 1e-7)
    assert residual_reference <= max(params.convergence_tolerance * 10.0, 1e-7)

    # Et les deux solutions decrivent bien le meme champ.
    scale = max(float(np.abs(reference).max()), 1e-30)
    assert np.abs(reference - values).max() / scale < 1e-5


def test_fem_assembly_matches_scipy_coo_path(restore_switch):
    """
    Assemblage DIRECT en CSR contre le chemin `coo_matrix(...).tocsr()`.

    Deux choses doivent coincider, et la seconde est la plus fragile :
      - les COEFFICIENTS, a la precision machine (l'ordre d'accumulation
        differe de `sum_duplicates`, donc pas au bit) ;
      - le MOTIF de sparsite, exactement. Un motif plus large passerait
        inapercu sur les valeurs (les entrees surnumeraires valant zero) tout
        en changeant le cout de toute factorisation en aval.
    """
    from pipeline.step4_field_solver import (
        _assemble_stiffness_matrix,
        _compute_shape_gradients,
    )

    _mesh_data, tet_mesh, _overhang, _regions = _step4_inputs()

    _cpp.set_enabled(False)
    gradients, volumes = _compute_shape_gradients(tet_mesh)
    expected = _assemble_stiffness_matrix(tet_mesh, gradients, volumes)

    _cpp.set_enabled(True)
    actual = _assemble_stiffness_matrix(tet_mesh, gradients, volumes)

    assert expected.shape == actual.shape
    assert expected.nnz == actual.nnz, "le motif de sparsite doit etre identique"
    assert np.array_equal(expected.indptr, actual.indptr)
    assert np.array_equal(expected.indices, actual.indices)

    scale = max(float(abs(expected).max()), 1e-30)
    assert abs(expected - actual).max() / scale < 1e-13

    assert expected.nnz > 400, f"jeu de test trop pauvre ({expected.nnz} entrees)"


def test_gradient_direction_term_matches_python(restore_switch):
    """
    Terme de direction : la matrice ET le second membre. Le second membre est
    le plus revelateur - c'est lui la force motrice du champ, et une erreur
    d'indexation entre tetraedres actifs et globaux y serait invisible sur la
    seule matrice.
    """
    from pipeline.step4_field_solver import (
        _assemble_gradient_direction_term,
        _compute_shape_gradients,
    )

    _mesh_data, tet_mesh, _overhang, _regions = _step4_inputs()

    _cpp.set_enabled(False)
    gradients, volumes = _compute_shape_gradients(tet_mesh)

    rng = np.random.default_rng(4201)
    directions = rng.normal(size=(tet_mesh.tet_count, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    weights = rng.uniform(0.5, 3.0, size=tet_mesh.tet_count)
    # Un sous-ensemble ACTIF strict : le chemin natif indexe alors
    # `shape_gradients` et `scale` par rang dans `active`, pas par tetraedre.
    active = np.zeros(tet_mesh.tet_count, dtype=bool)
    active[::3] = True

    expected_A, expected_b = _assemble_gradient_direction_term(
        tet_mesh, gradients, volumes, active, directions, weights
    )
    _cpp.set_enabled(True)
    actual_A, actual_b = _assemble_gradient_direction_term(
        tet_mesh, gradients, volumes, active, directions, weights
    )

    # Le motif natif est un SUR-ENSEMBLE assume : il couvre tous les
    # tetraedres du maillage, la ou le chemin scipy ne cree d'entrees que
    # pour les tetraedres ACTIFS. Les entrees surnumeraires valent exactement
    # zero, et elles ne coutent rien en pratique : `build_static_system`
    # additionne ce terme a la matrice de rigidite, dont le motif est deja
    # complet. Ce qui doit etre verifie est donc l'egalite MATHEMATIQUE, pas
    # l'egalite des motifs.
    assert actual_A.nnz >= expected_A.nnz
    difference = (actual_A - expected_A).tocoo()
    scale_A = max(float(abs(expected_A).max()), 1e-30)
    assert (np.abs(difference.data).max() if difference.nnz else 0.0) / scale_A < 1e-13

    scale_b = max(float(np.abs(expected_b).max()), 1e-30)
    assert np.abs(expected_b - actual_b).max() / scale_b < 1e-13

    assert int(active.sum()) < tet_mesh.tet_count, "le masque actif doit etre strict"
    assert scale_b > 0.0, "un second membre nul ne testerait rien"


def test_shared_face_topology_matches_separate_computation(restore_switch):
    """
    `_build_tet_tet_adjacency` accepte desormais les paires de faces deja
    calculees par `_map_boundary_faces_to_surface`, pour ne pas refaire la
    meme passe topologique. La matrice produite doit etre identique a celle
    du calcul separe - sinon l'optimisation changerait la propagation du
    risque de surplomb en aval.
    """
    from pipeline.step4_field_solver import (
        _build_tet_tet_adjacency,
        _map_boundary_faces_to_surface,
    )

    mesh_data, tet_mesh, _overhang, _regions = _step4_inputs()

    _cpp.set_enabled(True)
    *_head, face_pairs = _map_boundary_faces_to_surface(tet_mesh, mesh_data)
    assert face_pairs is not None, "le chemin natif doit fournir les paires internes"

    shared = _build_tet_tet_adjacency(tet_mesh, face_pairs)
    separate = _build_tet_tet_adjacency(tet_mesh)

    assert (shared != separate).nnz == 0
    assert shared.nnz > 100, "un maillage sans face interne ne testerait rien"


# ---------------------------------------------------------------------------
# Pas 8 : SAT triangle-triangle (step1_ingest, cpp/src/core/mesh_repair.cpp)
# ---------------------------------------------------------------------------
#
# step1_ingest.py est domine (~90%) par de la logique de reparation de mesh
# qui appelle trimesh/pymeshfix directement (slice_plane, booleans,
# triangulation, MeshFix) - pas du code maison portable. Le seul point chaud
# reel identifie par lecture est le test SAT triangle-triangle
# (`_triangles_cross`) et la boucle qui l'appelle sur toutes les paires
# candidates (`_find_self_intersecting_faces`). Transcription bit-exacte
# (memes 11 axes, meme tolerance _COPLANAR_EPS=1e-9) - voir
# cpp/include/apex/mesh_repair.hpp. Le broad-phase natif (BVH maison) n'a
# PAS besoin d'etre identique a bounds_tree (rtree, cote Python) : un
# sur-ensemble conservateur d'AABB exactes suffit, la narrow-phase SAT est
# deterministe et bit-exacte des deux cotes - contrairement au BVH de rayons
# (pas 6), ou le broad-phase de trimesh fait partie du comportement observable.

from pipeline.step1_ingest import (  # noqa: E402
    _find_self_intersecting_faces as _py_find_self_intersecting_faces,
)
from pipeline.step1_ingest import _triangles_cross as _py_triangles_cross  # noqa: E402


def test_triangles_cross_matches_python_on_crossing_pair(restore_switch):
    """Deux triangles qui se traversent reellement (pas juste un contact)."""
    t0 = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    t1 = np.array([[1.0, 1.0, -1.0], [1.0, -1.0, -1.0], [1.0, 0.0, 1.0]])

    _cpp.set_enabled(False)
    expected = _py_triangles_cross(t0, t1)
    assert expected is True, "fixture invalide : ces deux triangles doivent se croiser"

    _cpp.set_enabled(True)
    actual = bool(_cpp.module().triangles_cross(
        np.ascontiguousarray(t0), np.ascontiguousarray(t1)
    ))
    assert actual == expected


def test_triangles_cross_matches_python_on_edge_contact(restore_switch):
    """Deux triangles qui se touchent SANS se traverser (bord commun, coplanaires)."""
    t0 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    t1 = np.array([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])

    _cpp.set_enabled(False)
    expected = _py_triangles_cross(t0, t1)
    assert expected is False, "fixture invalide : contact simple, pas une intersection"

    _cpp.set_enabled(True)
    actual = bool(_cpp.module().triangles_cross(
        np.ascontiguousarray(t0), np.ascontiguousarray(t1)
    ))
    assert actual == expected


def test_triangles_cross_matches_python_random_pairs(restore_switch):
    """Balayage aleatoire, dont une fraction reelle de paires qui se croisent."""
    rng = np.random.default_rng(4242)
    mismatches = 0
    crossing_count = 0
    for _ in range(400):
        t0 = rng.uniform(-2.0, 2.0, size=(3, 3))
        t1 = rng.uniform(-2.0, 2.0, size=(3, 3))

        _cpp.set_enabled(False)
        expected = _py_triangles_cross(t0, t1)
        _cpp.set_enabled(True)
        actual = bool(_cpp.module().triangles_cross(
            np.ascontiguousarray(t0), np.ascontiguousarray(t1)
        ))

        if expected:
            crossing_count += 1
        if actual != expected:
            mismatches += 1

    assert crossing_count > 0, "aucune paire croisee tiree : jeu d'essai peu discriminant"
    assert mismatches == 0


def _self_intersecting_test_mesh() -> "trimesh.Trimesh":
    """
    Maillage synthetique watertight portant de VRAIES auto-intersections :
    une boite dont une face a ete remplacee par deux triangles fins tordus
    de facon a traverser une paroi voisine (motif similaire au recouvrement
    CSG mesure sur turkish.stl/death.stl - voir CLAUDE.md).
    """
    import trimesh as _trimesh

    box = _trimesh.creation.box(extents=(4.0, 4.0, 4.0))
    spike = _trimesh.Trimesh(
        vertices=np.array([
            [-3.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 3.0],
        ]),
        faces=np.array([[0, 1, 2]]),
        process=False,
    )
    return _trimesh.util.concatenate([box, spike])


def test_find_self_intersecting_faces_matches_python(restore_switch):
    mesh = _self_intersecting_test_mesh()
    triangles = np.ascontiguousarray(mesh.triangles, dtype=np.float64)
    faces = np.ascontiguousarray(mesh.faces, dtype=np.int64)

    _cpp.set_enabled(False)
    expected = _py_find_self_intersecting_faces(mesh)
    assert expected, "fixture invalide : doit contenir de vraies auto-intersections"

    _cpp.set_enabled(True)
    actual_array = _cpp.module().find_self_intersecting_faces(triangles, faces)
    actual = {int(i) for i in actual_array}

    assert actual == expected


def test_find_self_intersecting_faces_matches_python_on_clean_mesh(restore_switch):
    """Un maillage propre (sans auto-intersection) doit rendre un ensemble vide des deux cotes."""
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=5.0)
    triangles = np.ascontiguousarray(mesh.triangles, dtype=np.float64)
    faces = np.ascontiguousarray(mesh.faces, dtype=np.int64)

    _cpp.set_enabled(False)
    expected = _py_find_self_intersecting_faces(mesh)

    _cpp.set_enabled(True)
    actual_array = _cpp.module().find_self_intersecting_faces(triangles, faces)
    actual = {int(i) for i in actual_array}

    assert expected == set()
    assert actual == expected


def test_step1_dispatches_to_native_when_enabled(restore_switch):
    """
    `step1_ingest._find_self_intersecting_faces` (la fonction du pipeline,
    pas la primitive testee ci-dessus) doit rendre le meme resultat que le
    Python pur, que l'interrupteur soit actif ou non - non-regression du
    branchement dans le fichier du pipeline lui-meme.
    """
    from pipeline import step1_ingest

    mesh = _self_intersecting_test_mesh()

    _cpp.set_enabled(False)
    expected = step1_ingest._find_self_intersecting_faces(mesh)
    assert expected

    _cpp.set_enabled(True)
    actual = step1_ingest._find_self_intersecting_faces(mesh)

    assert actual == expected


# ---------------------------------------------------------------------------
# Pas 9 : segments d'intersection + winding number (_self_intersection_resolve,
# cpp/src/core/mesh_intersect.cpp)
# ---------------------------------------------------------------------------
#
# Mesure au profil (pipeline step1->step8, death.stl, resolution 0.8) :
# resolve_self_intersections pesait ~13.5s/40.4s, le plus gros poste Python
# restant apres le pas 8. Deux primitives PORTEES A TOLERANCE (pas bit-exact,
# meme discipline que step4/pas 7) : triangle_triangle_segment (resout un
# systeme 3x3 par cofacteurs plutot que LAPACK, ~1e-14 relatif mesure) et
# generalized_winding_number (formule fermee, ~1e-15 mesure - reduction
# numpy non identique terme a terme). Voir cpp/include/apex/mesh_intersect.hpp
# pour le detail et pour la limite topologique ASSUMEE de
# find_intersection_segments (l'ordre des segments par face n'est PAS
# reproductible face a rtree, divergence mineure mesuree sur turkish.stl :
# euler_number -34 vs -35, volume identique a 1.13e-5 relatif pres,
# vitesse x11.5 - decision explicite de l'utilisateur de garder ce gain).

from pipeline._self_intersection_resolve import (  # noqa: E402
    _intersection_segments as _py_intersection_segments,
)
from pipeline._self_intersection_resolve import (  # noqa: E402
    generalized_winding_number as _py_generalized_winding_number,
)
from pipeline._self_intersection_resolve import (  # noqa: E402
    triangle_triangle_segment as _py_triangle_triangle_segment,
)

_SEGMENT_TOLERANCE = 1e-9  # voir _NATIVE_SEGMENT_TOLERANCE, pipeline/_self_intersection_resolve.py


def test_triangle_triangle_segment_matches_python_random_pairs(restore_switch):
    """
    Balayage aleatoire : la primitive isolee (pas le pipeline complet) doit
    concorder a tolerance serree, avec une garde anti-faux-vert sur le
    nombre de paires reellement croisees tirees.
    """
    rng = np.random.default_rng(777)
    crossing_count = 0
    max_err = 0.0
    for _ in range(500):
        t0 = rng.uniform(-5.0, 5.0, size=(3, 3))
        t1 = rng.uniform(-5.0, 5.0, size=(3, 3))

        _cpp.set_enabled(False)
        expected = _py_triangle_triangle_segment(t0, t1)
        _cpp.set_enabled(True)
        actual = _cpp.module().triangle_triangle_segment(
            np.ascontiguousarray(t0), np.ascontiguousarray(t1)
        )

        assert (expected is None) == (actual is None)
        if expected is None:
            continue
        crossing_count += 1
        err = float(np.max(np.abs(np.asarray(expected) - actual)))
        max_err = max(max_err, err)
        assert err < _SEGMENT_TOLERANCE

    assert crossing_count > 0, "aucune paire croisee tiree : jeu d'essai peu discriminant"


def test_triangle_triangle_segment_degenerate_and_parallel_cases(restore_switch):
    """Cas limites explicites : triangle degenere, plans paralleles, contact ponctuel."""
    degenerate = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    normal_tri = np.array([[0.0, 0.0, -1.0], [1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]])

    _cpp.set_enabled(False)
    expected = _py_triangle_triangle_segment(degenerate, normal_tri)
    assert expected is None, "fixture invalide : triangle degenere attendu"
    _cpp.set_enabled(True)
    actual = _cpp.module().triangle_triangle_segment(
        np.ascontiguousarray(degenerate), np.ascontiguousarray(normal_tri)
    )
    assert actual is None

    coplanar_a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    coplanar_b = np.array([[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]])
    _cpp.set_enabled(False)
    expected = _py_triangle_triangle_segment(coplanar_a, coplanar_b)
    assert expected is None, "fixture invalide : plans paralleles attendus"
    _cpp.set_enabled(True)
    actual = _cpp.module().triangle_triangle_segment(
        np.ascontiguousarray(coplanar_a), np.ascontiguousarray(coplanar_b)
    )
    assert actual is None


def test_generalized_winding_number_matches_python(restore_switch):
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=5.0)
    triangles = np.ascontiguousarray(mesh.triangles, dtype=np.float64)

    rng = np.random.default_rng(42)
    points = rng.uniform(-8.0, 8.0, size=(300, 3))

    _cpp.set_enabled(False)
    expected = _py_generalized_winding_number(points, triangles)
    _cpp.set_enabled(True)
    actual = _cpp.module().generalized_winding_number(
        np.ascontiguousarray(points, dtype=np.float64), triangles
    )

    max_err = float(np.max(np.abs(expected - actual)))
    assert max_err < 1e-10, f"ecart inattendu : {max_err}"
    # Au moins un point clairement interieur (winding ~1) et un exterieur
    # (winding ~0) : garde anti-faux-vert contre un jeu d'essai degenere.
    assert np.any(expected > 0.9)
    assert np.any(expected < 0.1)


def _self_intersecting_turkish_style_mesh() -> "trimesh.Trimesh":
    """Deux nervures qui se traversent - meme motif que turkish.stl (docstring du module)."""
    box = trimesh.creation.box(extents=(4.0, 4.0, 4.0))
    spike = trimesh.Trimesh(
        vertices=np.array([
            [-3.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 3.0],
        ]),
        faces=np.array([[0, 1, 2]]),
        process=False,
    )
    return trimesh.util.concatenate([box, spike])


def test_intersection_segments_same_faces_and_counts(restore_switch):
    """
    `find_intersection_segments` (broad+narrow phase complet) doit affecter
    exactement les MEMES faces, avec le MEME nombre de segments par face,
    que la version Python - seul l'ORDRE a l'interieur d'une face peut
    differer (limite documentee, cpp/include/apex/mesh_intersect.hpp).
    """
    mesh = _self_intersecting_turkish_style_mesh()
    scale = float(mesh.scale) or 1.0
    min_length = 1e-7 * scale

    _cpp.set_enabled(False)
    expected = _py_intersection_segments(mesh, min_length)
    assert expected, "fixture invalide : doit produire des segments d'intersection"

    _cpp.set_enabled(True)
    actual = _py_intersection_segments(mesh, min_length)

    assert set(expected.keys()) == set(actual.keys())
    for face_index in expected:
        assert len(expected[face_index]) == len(actual[face_index])


def test_resolve_self_intersections_agrees_within_documented_tolerance(restore_switch):
    """
    Oracle bout-en-bout sur le vrai cas d'usage (motif turkish.stl) :
    resultat topologiquement PROCHE mais pas garanti identique - la
    divergence assumee est documentee dans mesh_intersect.hpp. On verifie
    ici que l'ecart mesure reste dans l'ordre de grandeur observe (volume
    relatif < 1e-3), pas une exactitude bit-a-bit.
    """
    from pipeline._self_intersection_resolve import resolve_self_intersections

    mesh = _self_intersecting_turkish_style_mesh()

    _cpp.set_enabled(False)
    expected = resolve_self_intersections(mesh)
    _cpp.set_enabled(True)
    actual = resolve_self_intersections(mesh)

    assert (expected is None) == (actual is None)
    if expected is None:
        return
    relative_volume_diff = abs(float(expected.volume) - float(actual.volume)) / max(
        abs(float(expected.volume)), 1e-9
    )
    assert relative_volume_diff < 1e-3
    assert expected.is_watertight == actual.is_watertight


# ---------------------------------------------------------------------------
# Pas 10 : interpolation inverse-distance (step6, cpp/src/core/interpolate.cpp)
# ---------------------------------------------------------------------------
#
# DERNIER portage cible par le profil : apres le pas 9, plus aucun poste
# Python isole ne depasse 2,1 s hors TetGen (22,4 s / 52 %, plancher dur), et
# `_interpolate_out_of_plane` etait ce dernier poste (2,10 s propre, 2853
# appels). Tolerance serree et non bit-exact - MESURE : ecart max 8,88e-16
# sur 200 tirages, 36/200 bit-identiques, du a l'ordre de reduction de numpy
# (meme phenomene qu'au pas 9). Voir cpp/include/apex/interpolate.hpp.

from pipeline.step6_toolpath2d import (  # noqa: E402
    _interpolate_out_of_plane as _py_interpolate_out_of_plane,
)


def test_interpolate_out_of_plane_matches_python(restore_switch):
    """Balayage aleatoire, tailles variees, dont des points coincidents forces."""
    rng = np.random.default_rng(11)
    max_err = 0.0
    coincident_trials = 0
    for trial in range(120):
        k = int(rng.integers(1, 60))
        m = int(rng.integers(1, 40))
        query = rng.uniform(-10.0, 10.0, size=(k, 2))
        sources = rng.uniform(-10.0, 10.0, size=(m, 2))
        offsets = rng.uniform(-2.0, 2.0, size=m)
        if trial % 7 == 0:
            # Coincidence EXACTE : exerce la branche `_COINCIDENT_EPS`, qui
            # doit rendre la valeur source telle quelle (jamais 0**-power).
            query[0] = sources[int(rng.integers(0, m))]
            coincident_trials += 1

        _cpp.set_enabled(False)
        expected = _py_interpolate_out_of_plane(query, sources, offsets, 2.0)
        _cpp.set_enabled(True)
        actual = _py_interpolate_out_of_plane(query, sources, offsets, 2.0)

        max_err = max(max_err, float(np.max(np.abs(expected - actual))))

    assert coincident_trials > 0, "aucun cas coincident tire : branche non exercee"
    assert max_err < 1e-12, f"ecart inattendu : {max_err}"


def test_interpolate_out_of_plane_empty_inputs(restore_switch):
    """Les deux cas vides du Python (aucune requete / aucune source)."""
    sources = np.array([[0.0, 0.0], [1.0, 1.0]])
    offsets = np.array([0.5, -0.5])

    for native in (False, True):
        _cpp.set_enabled(native)
        assert len(_py_interpolate_out_of_plane(np.zeros((0, 2)), sources, offsets, 2.0)) == 0
        empty_source = _py_interpolate_out_of_plane(
            np.array([[3.0, 3.0]]), np.zeros((0, 2)), np.zeros(0), 2.0
        )
        assert np.array_equal(empty_source, np.zeros(1))


def test_interpolate_out_of_plane_returns_source_value_on_exact_hit(restore_switch):
    """Une requete posee EXACTEMENT sur une source rend son offset, pas une moyenne."""
    sources = np.array([[0.0, 0.0], [10.0, 0.0]])
    offsets = np.array([1.0, -1.0])
    query = np.array([[0.0, 0.0]])

    for native in (False, True):
        _cpp.set_enabled(native)
        result = _py_interpolate_out_of_plane(query, sources, offsets, 2.0)
        assert result[0] == pytest.approx(1.0, abs=1e-12)


# ---------------------------------------------------------------------------
# T9 (2026-09-30) : distance exacte point -> maillage (MeshDistance)
# ---------------------------------------------------------------------------

def _mesh_distance_or_skip():
    from core import _cpp
    mod = _cpp.module() if _cpp.available() else None
    if mod is None or not hasattr(mod, "MeshDistance"):
        pytest.skip("apex_core sans MeshDistance (module non reconstruit)")
    return mod


@pytest.mark.parametrize("backend", [0, 1])
def test_mesh_distance_matches_trimesh(backend):
    mod = _mesh_distance_or_skip()
    mesh = trimesh.creation.icosphere(subdivisions=3, radius=10.0)
    mesh.apply_transform(trimesh.transformations.scale_and_translate([1.0, 0.6, 1.3], [1.0, -2.0, 0.5]))
    rng = np.random.default_rng(3)
    pts = rng.uniform(-20, 20, (5000, 3))
    pts[:200] = mesh.vertices[:200] + 1e-3 * rng.normal(size=(200, 3))     # tres pres de la surface
    ref_q, ref_d, _ = trimesh.proximity.closest_point(mesh, pts)
    index = mod.MeshDistance(np.asarray(mesh.vertices, dtype=float), np.asarray(mesh.faces, dtype=np.int64))
    d, q = index.query(pts, backend=backend, threads=4)
    # trimesh n'est PAS exact tout pres des sommets (mesure : jusqu'a 1,9e-5 mm
    # d'erreur sur 15 points sur 5000, toujours par exces). La ou les deux
    # different, l'arbitre est la force brute sur tous les triangles.
    assert np.max(np.abs(d - ref_d)) < 1e-4
    assert np.all(d <= ref_d + 1e-12)
    bad = np.flatnonzero(np.abs(d - ref_d) > 1e-9)
    tri = mesh.triangles
    for i in bad:
        brute = np.linalg.norm(trimesh.triangles.closest_point(tri, np.repeat(pts[i][None], len(tri), 0)) - pts[i],
                               axis=1).min()
        assert abs(d[i] - brute) < 1e-12
    np.testing.assert_allclose(np.linalg.norm(q - pts, axis=1), d, atol=1e-12)
    good = np.abs(d - ref_d) <= 1e-9
    assert np.max(np.linalg.norm(q - ref_q, axis=1)[good & (ref_d > 1e-6)]) < 1e-6


def test_mesh_distance_vulkan_falls_back_to_cpu():
    mod = _mesh_distance_or_skip()
    backend, reason = mod.mesh_distance_backend(2)
    assert backend in (1, 2)
    if backend == 1:
        assert reason


# ---------------------------------------------------------------------------
# Pas 13 - marching tetrahedra + soudure (step5), equivalence AU BIT
# ---------------------------------------------------------------------------

def _random_tet_field(seed: int, n: int = 600):
    from scipy.spatial import Delaunay
    rng = np.random.default_rng(seed)
    pts = rng.uniform(-10.0, 10.0, size=(n, 3))
    tets = Delaunay(pts).simplices.astype(np.int64)
    values = 0.3 * pts[:, 0] + np.sin(pts[:, 1] * 0.4) + 0.05 * pts[:, 2] ** 2 + rng.normal(0, 0.05, n)
    return pts, tets, values


def _python_extract(pts, tets, values, iso, tol):
    from core.mesh_types import TetMesh
    from pipeline.step5_layers import _extract_isosurface_triangles, _weld_vertices
    tm = TetMesh(tet_vertices=pts, tetrahedra=tets, tet_adjacency=None)
    per = values[tets]
    v, f, src = _extract_isosurface_triangles(tm, values, iso, per.min(axis=1), per.max(axis=1))
    wv, wf = _weld_vertices(v, f, tol)
    return wv, wf, src


def test_isosurface_extractor_bit_exact_vs_python():
    core = _cpp.module()
    if not hasattr(core, "IsosurfaceExtractor"):
        pytest.skip("module natif anterieur au pas 13")
    total_caps = total_quads = merged = 0
    for seed in range(4):
        pts, tets, values = _random_tet_field(seed)
        ex = core.IsosurfaceExtractor(pts, tets, values)
        isos = list(np.quantile(values, [0.1, 0.3, 0.5, 0.7, 0.9])) + [float(values[7])]  # iso sur un noeud exact
        for iso in isos:
            for tol in (0.0, 1e-6, 0.05):
                wv, wf, src = _python_extract(pts, tets, values, float(iso), tol)
                nv, nf, nsrc = ex.extract(float(iso), tol)
                assert np.array_equal(src, nsrc)
                assert np.array_equal(wf, nf)
                assert np.array_equal(wv, nv)          # au bit : aucune reduction
                if tol > 0:
                    merged += 3 * len(wf) - len(wv)
            per = values[tets]
            above = (per >= iso).sum(axis=1)
            cand = (per.min(axis=1) <= iso) & (iso <= per.max(axis=1))
            total_caps += int(((above == 1) | (above == 3))[cand].sum())
            total_quads += int((above == 2)[cand].sum())
    # garde anti-faux-vert : les deux cas geometriques et la soudure ont bien ete exerces
    assert total_caps > 100 and total_quads > 100 and merged > 100


def test_isosurface_extractor_empty_and_validation():
    core = _cpp.module()
    if not hasattr(core, "IsosurfaceExtractor"):
        pytest.skip("module natif anterieur au pas 13")
    pts, tets, values = _random_tet_field(9, 200)
    ex = core.IsosurfaceExtractor(pts, tets, values)
    v, f, s = ex.extract(float(values.max()) + 1.0, 1e-6)
    assert v.shape == (0, 3) and f.shape == (0, 3) and s.shape == (0,)
    with pytest.raises(ValueError):
        core.IsosurfaceExtractor(pts, np.array([[0, 1, 2, 10**9]]), values)


# ---------------------------------------------------------------------------
# Pas 14 - trace de bord natif (outline), identique a trimesh sur boucles simples
# ---------------------------------------------------------------------------

def _outline_python(vertices, faces):
    previous = _cpp.set_enabled(False)
    try:
        from pipeline.step5_layers import _extract_closed_contours
        return _extract_closed_contours(vertices, faces)
    finally:
        _cpp.set_enabled(previous)


def test_boundary_loops_match_trimesh_outline(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "boundary_loops") or not hasattr(core, "IsosurfaceExtractor"):
        pytest.skip("module natif anterieur au pas 14")
    import pipeline.step5_layers as s5
    loops_seen = native_used = 0
    for seed in range(6):
        pts, tets, values = _random_tet_field(seed, 800)
        ex = core.IsosurfaceExtractor(pts, tets, values)
        for iso in np.quantile(values, [0.15, 0.35, 0.5, 0.65, 0.85]):
            v, f, _ = ex.extract(float(iso), 1e-6)
            expected = _outline_python(v, f)
            _cpp.set_enabled(True)
            before = s5._NATIVE_OUTLINE_FALLBACKS
            got = s5._extract_closed_contours(v, f)
            native_used += int(s5._NATIVE_OUTLINE_FALLBACKS == before)
            assert len(got) == len(expected)
            for a, b in zip(got, expected):
                assert np.array_equal(a, b)
            loops_seen += len(expected)
    assert loops_seen > 20 and native_used > 20   # garde anti-faux-vert


def test_boundary_loops_fall_back_on_non_manifold_border():
    core = _cpp.module()
    if not hasattr(core, "boundary_loops"):
        pytest.skip("module natif anterieur au pas 14")
    bowtie = np.array([[0, 1, 2], [0, 3, 4]], dtype=np.int64)   # sommet 0 de degre 4
    assert core.boundary_loops(bowtie, 5) is None
    square = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    loops = core.boundary_loops(square, 4)
    assert [list(l) for l in loops] == [[0, 1, 2, 3, 0]]
    degenerate = np.array([[0, 0, 1]], dtype=np.int64)
    assert core.boundary_loops(degenerate, 2) == []


# ---------------------------------------------------------------------------
# Pas 17 - IK couplee a la position (step8, machine actuelle), tolerance declaree
# ---------------------------------------------------------------------------

def test_coupled_ik_native_matches_python(restore_switch):
    import pipeline.step8_machine_toolpath as s8
    from core.mesh_types import MachineProfile, RotationAxisSpec
    machine = MachineProfile(
        rotation_axes=[RotationAxisSpec(axis=np.array([0.0, 0.0, 1.0]), angle_range_deg=(0.0, 360.0)),
                       RotationAxisSpec(axis=np.array([1.0, 0.0, 0.0]), angle_range_deg=(0.0, 185.0))],
        platform_origin=np.zeros(3), position_coupled_axes=[0], radial_direction=np.array([0.0, -1.0, 0.0]))
    rng = np.random.default_rng(17)
    pos = rng.uniform(-40, 40, (600, 3))
    pos[::50, :2] = 0.0                                     # points SUR l'axe : heritage de C
    dirs = rng.normal(size=(600, 3))
    dirs[:, 2] = np.abs(dirs[:, 2]) + 0.3
    dirs[::7] = [0.0, 0.0, 1.0]                             # depots verticaux : repli grille du natif
    raw = [s8._RawPoint(position=p, desired_direction=d / np.linalg.norm(d), layer_id=0, segment_index=0,
                        color_id=None, feedrate=1.0, extrusion=1.0, is_travel=False) for p, d in zip(pos, dirs)]
    ref = np.array([0.0, 0.0, 1.0])
    params = s8.Step8Params()
    _cpp.set_enabled(False)
    a_py, r_py = s8._resolve_all_axis_angles_coupled(raw, machine.rotation_axes, ref, machine, params)
    _cpp.set_enabled(True)
    a_cc, r_cc = s8._resolve_all_axis_angles_coupled(raw, machine.rotation_axes, ref, machine, params)
    a_py, a_cc = np.array(a_py), np.array(a_cc)
    assert np.abs(a_py - a_cc).max() < 1e-7            # degres
    assert np.abs(np.array(r_py) - np.array(r_cc)).max() < 1e-6   # radians (acos mal conditionnee pres de 0)
    # memes diagnostics a 5 deg : l'assertion qui compte reellement (pas 2)
    lim = np.radians(params.max_unreachable_angular_error_deg)
    assert np.array_equal(np.array(r_py) > lim, np.array(r_cc) > lim)
    assert (np.array(r_py) > lim).sum() > 50          # garde : des cibles inatteignables existent


# ---------------------------------------------------------------------------
# Pas 12 - point_in_polygon_2d natif, au bit
# ---------------------------------------------------------------------------

def test_point_in_polygon_2d_native_matches_python(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "point_in_polygon_2d"):
        pytest.skip("module natif anterieur au pas 12")
    from core.geometry_utils import point_in_polygon_2d
    rng = np.random.default_rng(12)
    inside_count = total = 0
    for _ in range(60):
        n = int(rng.integers(6, 40))
        ang = np.sort(rng.uniform(0, 2 * np.pi, n))
        rad = rng.uniform(1.0, 5.0, n)
        poly = np.column_stack([rad * np.cos(ang), rad * np.sin(ang)])
        pts = rng.uniform(-6, 6, (80, 2))
        pts[:5] = poly[:5]                          # sommets exacts : cas limites
        pts[5:8, 1] = poly[:3, 1]                   # meme ordonnee qu'un sommet
        for p in pts:
            _cpp.set_enabled(False)
            a = point_in_polygon_2d(p, poly)
            _cpp.set_enabled(True)
            b = point_in_polygon_2d(p, poly)
            assert a == b
            inside_count += a
            total += 1
    assert 0.1 * total < inside_count < 0.9 * total   # garde : les deux issues exercees


# ---------------------------------------------------------------------------
# Pas 19 - matiere deja deposee (step9 rapide), minimum exact au bit
# ---------------------------------------------------------------------------

def test_past_material_index_matches_python_fast_checker(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "PastMaterialIndex"):
        pytest.skip("module natif anterieur au pas 19")
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(19)
    dense = rng.uniform(0, 30, (20000, 3))
    owner = rng.integers(0, 5000, 20000).astype(np.int64)
    idx = core.PastMaterialIndex(dense, owner)
    tree = cKDTree(dense)
    owner_pad = np.r_[owner, np.iinfo(np.int64).max]
    reach = 1.2
    hits = 0
    for trial in range(300):
        verts = rng.uniform(0, 30, (60, 3))
        limit = int(rng.integers(0, 5000))
        # reference : algorithme Python de FastCollisionChecker
        dist, ii = tree.query(verts, k=16, distance_upper_bound=reach)
        past = owner_pad[ii] < limit
        d_past = np.where(past, dist, np.inf).min(axis=1)
        sat = np.isfinite(dist[:, -1]) & ~np.isfinite(d_past)
        for vi in np.flatnonzero(sat):
            near = np.asarray(tree.query_ball_point(verts[vi], reach), dtype=int)
            near = near[owner[near] < limit]
            if len(near):
                d_past[vi] = float(np.linalg.norm(dense[near] - verts[vi], axis=1).min())
        ref = float(d_past.min())
        got = idx.min_past_distance(verts, limit, reach)
        assert (ref == got) or (not np.isfinite(ref) and not np.isfinite(got))
        hits += np.isfinite(ref)
    assert 50 < hits < 300          # garde : cas avec et sans matiere passee


# ---------------------------------------------------------------------------
# Pas 24 - export G-code natif, texte identique ligne a ligne
# ---------------------------------------------------------------------------

def _random_machine_toolpath(seed, n=3000):
    rng = np.random.default_rng(seed)
    steps = rng.normal(0, 0.4, (n, 3))
    steps[:, 2] = np.abs(steps[:, 2]) * 0.05
    pos = np.cumsum(steps, axis=0) + [20.0, 5.0, 0.1]
    pos[100:110, :2] = 0.0                                  # passages au centre (Y herite)
    jumps = rng.choice(n, 60, replace=False)
    pos[jumps] += rng.uniform(-15, 15, (60, 3)) * [1, 1, 0.1]   # trajets longs
    pos[:, 2] = np.maximum(pos[:, 2], 0.05)
    c = (np.degrees(np.arctan2(pos[:, 1], pos[:, 0])) + 90.0) % 360.0
    b = np.clip(rng.normal(20, 30, n), 0, 185)
    ext = (rng.uniform(size=n) > 0.08).astype(float)
    ext[0] = 0
    skin = rng.uniform(size=n) > 0.7
    return {"position": pos, "axis_angles": np.column_stack([c, b]), "extrusion": ext, "is_skin": skin}


@pytest.mark.parametrize("variant", ["defaut", "pivot", "sans_table", "decoupage"])
def test_gcode_export_native_identical_text(restore_switch, variant):
    core = _cpp.module()
    if not hasattr(core, "export_gcode_body"):
        pytest.skip("module natif anterieur au pas 24")
    import main
    from file_io.gcode_export import GcodeParams, export_gcode, head_zlow_from_table
    machine = main._default_machine_profile()
    machine.nozzle_model = None
    if variant == "decoupage":       # bisection max_tip_deviation en natif (2026-10-07)
        params = GcodeParams()
    elif variant == "pivot":
        params = GcodeParams(z_reference="pivot", tilt_sign=-1.0, max_tip_deviation=None)
    else:
        params = GcodeParams(max_tip_deviation=None)
    hz = None if variant == "sans_table" else head_zlow_from_table(np.arange(0, 185.5, 0.5), -np.linspace(1, 30, 371))
    for seed in range(3):
        tp = _random_machine_toolpath(seed)
        _cpp.set_enabled(False)
        a, ra = export_gcode(tp, machine, params, head_zlow=hz)
        _cpp.set_enabled(True)
        b, rb = export_gcode(tp, machine, params, head_zlow=hz)
        assert a == b
        for k in ("moves", "travels", "retractions", "clearance_axes", "filament_mm", "c_mismatch_deg",
                  "tip_z_min", "head_z_min", "tip_below_plate", "head_below_margin", "fw_x_range",
                  "fw_y_range", "fw_z_range", "print_time_s", "subdivided_moves", "max_tip_deviation_before"):
            assert getattr(ra, k) == getattr(rb, k), k
        if variant == "decoupage":
            assert ra.subdivided_moves > 50                                  # garde : le decoupage agit
        assert ra.retractions > 5 and ra.clearance_axes["R"] + ra.clearance_axes["Z"] > 5   # garde


def test_feed_plan_native_bit_exact(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "feed_plan_caps"):
        pytest.skip("module natif anterieur au pas 24")
    from file_io.gcode_export import feed_plan
    tp = _random_machine_toolpath(7, 20000)
    _cpp.set_enabled(False)
    a = feed_plan(tp["position"], tp["axis_angles"], 750.0, 300.0, [0], 25.0)
    _cpp.set_enabled(True)
    b = feed_plan(tp["position"], tp["axis_angles"], 750.0, 300.0, [0], 25.0)
    assert a[0] == b[0] and a[1] == b[1] and np.array_equal(a[2], b[2])
    assert (a[2] < 1.0).sum() > 100     # garde : le ralentissement est bien exerce


# ---------------------------------------------------------------------------
# Pas 15 : remplissage 2D d'une zone (step6), Clipper2 contre shapely/GEOS.
# PAS au bit (entiers 1e-6 mm, autre depart d'anneau) : aires des parois et
# longueurs de remplissage a tolerance declaree, sur des zones sans angle
# aigu (la ou les onglets de GEOS et Clipper2 different, voir CLAUDE.md #56).
# ---------------------------------------------------------------------------

def _star_region(rng):
    # contour lisse (pas d'angle aigu) : rayon modulé par une sinusoïde
    n = int(rng.integers(40, 120))
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    rad = 10.0 + rng.uniform(0.5, 2.0) * np.sin(int(rng.integers(2, 6)) * ang + rng.uniform(0, 6))
    outer = np.column_stack([rad * np.cos(ang), rad * np.sin(ang)])
    holes = []
    if rng.random() < 0.5:
        t = np.linspace(0, 2 * np.pi, 24, endpoint=False)
        holes.append(np.column_stack([2.5 * np.cos(t), 2.5 * np.sin(t)]))
    return outer, holes


def test_fill_region_native_matches_shapely(restore_switch):
    core = _cpp.module()
    if not hasattr(core, "fill_region"):
        pytest.skip("module natif anterieur au pas 15")
    from shapely.geometry import Polygon
    from pipeline import step6_toolpath2d as s6

    def length(polys):
        return sum(float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum()) for p in polys)

    rng = np.random.default_rng(15)
    codes = {s6.InfillPattern.LINES: 0, s6.InfillPattern.GRID: 1, s6.InfillPattern.CROSSED_ZIGZAG: 2}
    checked = 0
    for k in range(40):
        outer, holes = _star_region(rng)
        pattern = list(codes)[k % 3]
        angle = float(rng.uniform(0, np.pi))
        walls, innermost = s6._offset_walls(outer, holes, 0.4, 2, 1e-6)
        infill = s6._generate_infill(innermost, 0.2, pattern, 0.4, angle)
        valid, n_infill, n_walls = core.fill_region(outer, holes, 0.4, 2, 1e-6, 0.2, codes[pattern], angle)
        assert valid
        py_walls = walls[: len(walls) - 1 - len(holes)]
        assert len(n_walls) == len(py_walls)
        area_py = sum(abs(Polygon(w).area) for w in py_walls)
        area_c = sum(abs(Polygon(w).area) for w in n_walls)
        assert abs(area_c - area_py) <= 1e-3 * area_py
        for w in n_walls:
            assert np.array_equal(w[0], w[-1])          # anneaux fermes, comme shapely
        assert abs(length(n_infill) - length(infill)) <= 1e-3 * max(length(infill), 1.0)
        checked += length(infill) > 0
    assert checked >= 30, "le jeu d'essai doit produire du remplissage"


def test_fill_region_orientation_and_degenerate():
    core = _cpp.module()
    if not hasattr(core, "fill_region"):
        pytest.skip("module natif anterieur au pas 15")
    from shapely.geometry import LinearRing
    square = np.array([[0, 0], [10, 0], [10, 10], [0, 10]], float)
    hole = np.array([[4, 4], [6, 4], [6, 6], [4, 6]], float)
    valid, _infill, walls = core.fill_region(square, [hole], 0.4, 1, 0.0, 0.2, 2, 0.0)
    assert valid and len(walls) == 2
    # orientation de GEOS : exterieur horaire, trou anti-horaire
    assert not LinearRing(walls[0]).is_ccw and LinearRing(walls[1]).is_ccw
    # zone d'aire nulle : ignoree entierement (comme `_offset_walls` -> [])
    flat = np.array([[0, 0], [5, 0], [10, 0]], float)
    valid, infill, walls = core.fill_region(flat, [], 0.4, 2, 0.0, 0.2, 2, 0.0)
    assert not valid and not infill and not walls


# ---------------------------------------------------------------------------
# Pas 21a : etape 2 entiere en natif (surplombs, risque, ponts, portee).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("shape", ["hollow", "l_shape", "sphere"])
def test_overhang_map_native_matches_python(restore_switch, shape):
    core = _cpp.module()
    if not hasattr(core, "overhang_map"):
        pytest.skip("module natif anterieur au pas 21a")
    import os
    import tempfile
    import trimesh
    from pipeline import step2_overhang as s2
    from pipeline.step1_ingest import run as ingest_run

    if shape == "hollow":
        mesh = trimesh.creation.box(extents=[10, 10, 10]).difference(
            trimesh.creation.box(extents=[4, 4, 4]), engine="manifold")
    elif shape == "l_shape":
        v = trimesh.creation.box(extents=[4, 4, 20]); v.apply_translation([0, 0, 10])
        h = trimesh.creation.box(extents=[20, 4, 4]); h.apply_translation([8, 0, 20])
        mesh = v.union(h, engine="manifold")
    else:
        mesh = trimesh.creation.icosphere(subdivisions=2, radius=5.0)
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, _, _ = ingest_run(path, resolution=2.0)
    finally:
        os.unlink(path)
    params = s2.OverhangParams()
    _cpp.set_enabled(False)
    ref = s2.run(mesh_data, params)
    _cpp.set_enabled(True)
    assert s2._native_overhang_map(mesh_data, params) is not None, "le natif doit avoir tourne"
    nat = s2.run(mesh_data, params)
    assert np.array_equal(ref.is_overhang, nat.is_overhang)
    assert np.array_equal(ref.is_bridge, nat.is_bridge)
    assert np.allclose(ref.risk_score, nat.risk_score, rtol=0, atol=1e-12)
    assert np.array_equal(np.isnan(ref.bridge_span_direction), np.isnan(nat.bridge_span_direction))
    # La portee est un axe principal : unique (au signe canonique pres) sauf
    # quand deux variances sont egales (plafond carre). On verifie donc ce qui
    # est garanti : meme variance le long des deux directions, par composante.
    bridges = np.flatnonzero(ref.is_bridge)
    for row in np.unique(ref.bridge_span_direction[bridges], axis=0):
        faces = bridges[np.all(ref.bridge_span_direction[bridges] == row, axis=1)]
        pts = mesh_data.vertices[np.unique(mesh_data.faces[faces])]
        pts = pts - pts.mean(axis=0)
        native_rows = np.unique(nat.bridge_span_direction[faces], axis=0)
        assert len(native_rows) == 1
        var_ref, var_nat = np.var(pts @ row), np.var(pts @ native_rows[0])
        assert abs(var_ref - var_nat) <= 1e-9 * max(var_ref, 1e-12)
    if shape == "hollow":
        assert ref.is_bridge.any(), "la cavite doit produire un pont, sinon le test ne verifie rien"


# ---------------------------------------------------------------------------
# Pas 21b : etape 3, chemin par defaut, IDENTIQUE (ordre des dependances compris).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("merge", [None, 30.0])
def test_segmentation_native_matches_python(restore_switch, merge):
    core = _cpp.module()
    if not hasattr(core, "segment_default"):
        pytest.skip("module natif anterieur au pas 21b")
    import os
    import tempfile
    from dataclasses import replace
    import trimesh
    from pipeline import step2_overhang as s2, step3_segmentation as s3
    from pipeline.step1_ingest import run as ingest_run

    # T couche + bras pendant : plusieurs regions et dependances
    v = trimesh.creation.box(extents=[4, 4, 20]); v.apply_translation([0, 0, 10])
    h = trimesh.creation.box(extents=[24, 4, 4]); h.apply_translation([0, 0, 22])
    d = trimesh.creation.box(extents=[4, 4, 10]); d.apply_translation([10, 0, 15])
    mesh = v.union(h, engine="manifold").union(d, engine="manifold")
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        path = f.name
    mesh.export(path)
    try:
        mesh_data, tet_mesh, report = ingest_run(path, resolution=2.0)
    finally:
        os.unlink(path)
    om = s2.run(mesh_data, s2.OverhangParams())
    params = replace(s3.SegmentationParams(), platform_z=report.platform_z_used, max_region_merge_angle_deg=merge)
    _cpp.set_enabled(False)
    a = s3.run(mesh_data, om, params, tet_mesh)
    _cpp.set_enabled(True)
    b = s3.run(mesh_data, om, params, tet_mesh)
    if merge is not None:
        assert len(a.regions) > 1, "la fixture doit produire plusieurs regions"
    assert a.regions.keys() == b.regions.keys()
    assert a.dependencies == b.dependencies
    assert a.region_print_depth == b.region_print_depth
    for k in a.regions:
        assert np.array_equal(a.regions[k].face_ids, b.regions[k].face_ids)
        assert np.array_equal(a.regions[k].target_direction, b.regions[k].target_direction)


# ---------------------------------------------------------------------------
# Pas 21c : etape 7 entiere (couleur, lissage, trous, ordonnancement).
# ---------------------------------------------------------------------------

def _schedule_signature(colored, schedule):
    layers = [([np.asarray(p).tolist() for p in (l.toolpath_2d or [])], list(l.segment_colors or []), l.color_id)
              for l in colored]
    passes = [(list(p.segments), p.color_id, p.is_hole_fill, tuple(p.layer_range), list(p.deferred_regions))
              for p in schedule.passes]
    return layers, passes


@pytest.mark.parametrize("case", ["none", "two_tone", "lateral", "lateral_merge", "cheap_holes"])
def test_schedule_native_matches_python(restore_switch, case):
    core = _cpp.module()
    if not hasattr(core, "schedule_layers"):
        pytest.skip("module natif anterieur au pas 21c")
    import trimesh
    from pipeline import step7_scheduling as s7
    from tests.test_step7_scheduling import (
        _lateral_color_map, _pipeline_through_step6, _two_tone_color_map,
    )
    mesh = trimesh.creation.box(extents=[10, 10, 10]).subdivide().subdivide()
    mesh_data, tet_mesh, regions, field, layers = _pipeline_through_step6(mesh, platform_z=-5.0)
    params = s7.Step7Params()
    color_map = None
    if case == "two_tone":
        color_map = _two_tone_color_map(mesh_data)
    elif case == "lateral":
        color_map = _lateral_color_map(mesh_data, split_x=0.7)
    elif case == "lateral_merge":
        color_map = _lateral_color_map(mesh_data, split_x=0.7, merge_tolerance=3.0)
    elif case == "cheap_holes":
        color_map = _lateral_color_map(mesh_data, split_x=0.7)
        params = s7.Step7Params(
            hole_cost_model=s7.HoleCostModel(reheat_penalty=0.0, adhesion_risk_weight=0.0),
            purge_cost_model=s7.PurgeCostModel(volume_per_purge=1000.0, time_per_purge=1000.0),
            max_open_holes=100,
        )
    _cpp.set_enabled(False)
    ref = s7.run(layers, mesh_data, tet_mesh, color_map, regions, field, params)
    _cpp.set_enabled(True)
    assert s7._native_schedule(layers, mesh_data, tet_mesh, color_map, field, params) is not None
    nat = s7.run(layers, mesh_data, tet_mesh, color_map, regions, field, params)
    sig_ref, sig_nat = _schedule_signature(*ref), _schedule_signature(*nat)
    assert sig_ref[1] == sig_nat[1]
    assert len(sig_ref[0]) == len(sig_nat[0])
    for (pa, ca, da), (pb, cb, db) in zip(sig_ref[0], sig_nat[0]):
        assert ca == cb and da == db and len(pa) == len(pb)
        for a, b in zip(pa, pb):
            assert np.array_equal(np.asarray(a), np.asarray(b))
    if case == "cheap_holes":
        assert any(p.is_hole_fill for p in nat[1].passes), "la branche trou doit etre exercee"


def test_past_material_batch_matches_single_calls():
    """Requete par lots, en parallele : identique AU BIT aux appels unitaires."""
    core = _cpp.module()
    if not hasattr(core.PastMaterialIndex, "min_past_distance_batch"):
        pytest.skip("module natif anterieur au lot parallele")
    rng = np.random.default_rng(19)
    pts = rng.uniform(0, 20, (5000, 3))
    owners = rng.integers(0, 400, 5000).astype(np.int64)
    index = core.PastMaterialIndex(pts, owners)
    poses = rng.uniform(0, 20, (300, 40, 3))
    limits = rng.integers(0, 400, 300).astype(np.int64)
    batch = index.min_past_distance_batch(poses, limits, 1.5, 4)
    single = np.array([index.min_past_distance(np.ascontiguousarray(poses[k]), int(limits[k]), 1.5)
                       for k in range(300)])
    assert np.array_equal(batch, single)
    assert np.isfinite(batch).sum() > 30, "des poses doivent toucher de la matiere"
