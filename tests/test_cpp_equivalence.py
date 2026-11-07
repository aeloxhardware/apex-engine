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

import dataclasses

import numpy as np
import pytest

from core import _cpp, runtime
from core.mesh_types import AABB
from element_ardor.step9_collision_sim import _PrintedGeometryIndex

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
    from element_ardor.step8_machine_toolpath import _resolve_all_axis_angles

    import element_ardor.step8_machine_toolpath as step8

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
            assert runtime.state().native_fallbacks is not None, (
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
    from element_ardor.step8_machine_toolpath import Step8Params, _RawPoint

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
    import element_ardor.step8_machine_toolpath as step8

    raw_points, axes, reference, origin, params = _ik_fixture(
        _current_machine_axes(), seed=71, count=150
    )
    up = np.array([0.0, 0.0, 1.0])
    for index in (0, 37, 149):
        raw_points[index].desired_direction = up

    py_result, cpp_result = _resolve_both_ways(raw_points, axes, reference, origin, params)

    assert runtime.state().native_fallbacks >= 3, (
        "les cibles +Z doivent declencher le repli grille ; sans repli, ce test "
        f"n'exerce pas la reprise (count={runtime.state().native_fallbacks})"
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
    from element_ardor.step1_ingest import _build_tet_adjacency

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

from element_ardor.step1_ingest import (  # noqa: E402
    _find_self_intersecting_faces as _py_find_self_intersecting_faces,
)
from element_ardor.step1_ingest import _triangles_cross as _py_triangles_cross  # noqa: E402


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

from element_ardor._self_intersection_resolve import (  # noqa: E402
    _intersection_segments as _py_intersection_segments,
)
from element_ardor._self_intersection_resolve import (  # noqa: E402
    generalized_winding_number as _py_generalized_winding_number,
)
from element_ardor._self_intersection_resolve import (  # noqa: E402
    triangle_triangle_segment as _py_triangle_triangle_segment,
)


#
# DERNIER portage cible par le profil : apres le pas 9, plus aucun poste
# Python isole ne depasse 2,1 s hors TetGen (22,4 s / 52 %, plancher dur), et
# `_interpolate_out_of_plane` etait ce dernier poste (2,10 s propre, 2853
# appels). Tolerance serree et non bit-exact - MESURE : ecart max 8,88e-16
# sur 200 tirages, 36/200 bit-identiques, du a l'ordre de reduction de numpy
# (meme phenomene qu'au pas 9). Voir cpp/include/apex/interpolate.hpp.

from element_ardor.step6_toolpath2d import (  # noqa: E402
    _interpolate_out_of_plane as _py_interpolate_out_of_plane,
)


# ---------------------------------------------------------------------------
# Pas 15 : remplissage 2D d'une zone (step6), Clipper2 contre shapely/GEOS.
# PAS au bit (entiers 1e-6 mm, autre depart d'anneau) : aires des parois et
# longueurs de remplissage a tolerance declaree, sur des zones sans angle
# aigu (la ou les onglets de GEOS et Clipper2 different, voir CLAUDE.md #56).
# ---------------------------------------------------------------------------


# --- Worker TetGen C++ (executable separe, limite #56) -----------------------
