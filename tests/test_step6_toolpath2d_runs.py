"""Suite découpée de `test_step6_toolpath2d.py` (R-030 point 1, limite de 800 lignes par fichier) : tests déplacés tels quels."""

from __future__ import annotations

from dataclasses import replace
import numpy as np
import pytest
from core.mesh_types import Layer
from pipeline.step6_toolpath2d import InfillPattern, Step6Params, _generate_infill, _generate_spiral_infill, _offset_walls, _lift_to_3d, _oscillating_spiral_polyline, _SpiralShape, _plate_axis_in_local_2d, _merge_short_runs, _monotonic_runs, _project_to_local_2d, _split_ring_into_runs, _split_rings_in_toolpath, run as toolpath_run
from test_step6_toolpath2d import (  # noqa: F401
    _wavy_square_layer,
    flat_box_layers,
)


def _unsupported_indices(ring: np.ndarray, horizontal_radius: float, vertical_height: float) -> list[int]:
    """
    Meme geometrie que le test de support de step9b/step5 (segments
    degeneres a un point), appliquee directement a un anneau (K,3) deja
    depose dans l'ordre - sert a MESURER l'effet du reordonnancement dans
    ces tests, independamment de l'implementation interne.
    """
    unsupported = []
    for i in range(1, len(ring)):
        pos = ring[i]
        prior = ring[:i]
        dz = pos[2] - prior[:, 2]
        horiz = np.hypot(pos[0] - prior[:, 0], pos[1] - prior[:, 1])
        if not np.any((dz > 1e-9) & (dz <= vertical_height) & (horiz <= horizontal_radius)):
            unsupported.append(i)
    return unsupported


def _wavy_ring(n: int, radius: float, amplitude: float, cycles: int) -> np.ndarray:
    """Anneau ferme (K+1,3) dont l'offset (3e colonne) ondule sinusoidalement -
    pas de plateaux (jamais deux offsets consecutifs egaux), point d'espacement
    dense (nozzle-scale) pour que le mecanisme ait une chance geometrique
    reelle de fonctionner (cf. diagnostic : sur un anneau grossier, l'espacement
    entre sommets depasse deja le rayon de recherche independamment de l'ordre)."""
    theta = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    u = radius * np.cos(theta)
    v = radius * np.sin(theta)
    offset = amplitude * np.sin(theta * cycles)
    ring = np.stack([u, v, offset], axis=1)
    return np.vstack([ring, ring[0]])


def test_monotonic_runs_single_run_for_flat_ring():
    """Un anneau plat (comportement historique sans interpolation, offset
    partout nul) ne produit qu'UN SEUL troncon - rien a reordonner."""
    ring = _wavy_ring(n=20, radius=5.0, amplitude=0.0, cycles=3)
    runs = _monotonic_runs(ring[:, 2])
    assert len(runs) == 1
    assert runs[0] == (0, len(ring))


def test_monotonic_runs_partition_without_overlap():
    """Les troncons sont des intervalles [start, stop) DISJOINTS qui
    partitionnent exactement l'anneau - aucun sommet partage ni oublie (bug
    reel trouve et corrige : une version precedente partageait un sommet a
    chaque frontiere, dupliquant des points une fois les troncons reordonnes)."""
    ring = _wavy_ring(n=64, radius=10.0, amplitude=2.0, cycles=6)
    runs = _monotonic_runs(ring[:, 2])
    assert runs[0][0] == 0
    assert runs[-1][1] == len(ring)
    for (_, stop_a), (start_b, _) in zip(runs, runs[1:]):
        assert stop_a == start_b
    total = sum(stop - start for start, stop in runs)
    assert total == len(ring)


def test_split_ring_preserves_every_point():
    """Le decoupage ne peut ni ajouter ni supprimer de points : la
    concatenation des troncons contient exactement le meme multi-ensemble
    d'offsets que l'anneau d'origine (bug reel deja mesure : une version
    precedente partageait un sommet a chaque frontiere et ressortait avec
    PLUS de points que l'entree)."""
    ring = _wavy_ring(n=200, radius=25.0, amplitude=2.0, cycles=6)
    runs = _split_ring_into_runs(ring)

    total = sum(len(r) for r in runs)
    assert total == len(ring)
    concatenated = np.concatenate(runs, axis=0)
    assert np.allclose(sorted(ring[:, 2].tolist()), sorted(concatenated[:, 2].tolist()))


def test_split_ring_runs_are_all_ascending():
    """
    Coeur du fix (Decision 5) : chaque troncon ressort parcouru dans le sens
    MONTANT (un troncon descendant est inverse). C'est ce qui le rend
    auto-portant - chaque point s'appuie sur celui qui vient d'etre depose
    juste avant lui, plus bas et proche horizontalement.
    """
    ring = _wavy_ring(n=400, radius=25.0, amplitude=2.0, cycles=6)
    runs = _split_ring_into_runs(ring)

    assert len(runs) > 1, "la fixture doit vraiment onduler, sinon le test ne verifie rien"
    for run in runs:
        if len(run) >= 2:
            assert run[-1, 2] >= run[0, 2]


def test_split_ring_runs_have_at_least_two_points():
    """
    Un troncon d'UN SEUL point deviendrait une polyligne d'un point : step8 y
    insererait un trajet a vide puis un point isole, et
    step9._grow_printed_index_to (qui exige DEUX extremites en extrusion) n'en
    ferait aucun cordon - ce point cesserait d'etre imprime. `_merge_short_runs`
    doit donc les fusionner dans leur voisin, sans jamais perdre de point.
    """
    # ondulation a haute frequence : force des troncons tres courts
    ring = _wavy_ring(n=40, radius=5.0, amplitude=1.0, cycles=9)
    runs = _split_ring_into_runs(ring)

    assert sum(len(r) for r in runs) == len(ring)
    for run in runs:
        assert len(run) >= 2


def test_split_ring_noop_on_flat_ring():
    """Un anneau plat (un seul troncon) ressort en UNE polyligne identique -
    non-regression stricte du comportement historique (offset nul)."""
    ring = _wavy_ring(n=30, radius=5.0, amplitude=0.0, cycles=4)
    runs = _split_ring_into_runs(ring)
    assert len(runs) == 1
    assert np.allclose(runs[0], ring)


def test_split_rings_in_toolpath_skips_open_polylines():
    """Une ligne d'infill OUVERTE n'est jamais un anneau et ressort
    strictement inchangee : une droite clippee n'a pas de "retour au point de
    depart" a gerer."""
    open_line = np.array([[0.0, 0.0, 2.0], [1.0, 0.0, 1.0], [2.0, 0.0, -2.0], [3.0, 0.0, 0.0]])
    closed_ring = _wavy_ring(n=64, radius=10.0, amplitude=2.0, cycles=6)

    result = _split_rings_in_toolpath([open_line, closed_ring])

    assert np.allclose(result[0], open_line)
    assert len(result) > 2, "l'anneau ondule doit avoir ete decoupe en plusieurs troncons"
    rebuilt = np.concatenate(result[1:], axis=0)
    assert np.allclose(sorted(rebuilt[:, 2].tolist()), sorted(closed_ring[:, 2].tolist()))


def test_split_disabled_by_default_is_bit_identical(flat_box_layers):
    """Comportement historique inchange quand `split_rings_into_runs` reste a
    son defaut (False) - meme discipline opt-in que les autres mecanismes de
    ce fichier (out_of_plane_interpolation_power, etc.)."""
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _wavy_square_layer()
    params = Step6Params(
        nozzle_width=0.4, wall_count=2, infill_density=0.3,
        out_of_plane_interpolation_power=2.0,
    )

    default = toolpath_run([layer], mesh_data, overhang_map, params)[0].toolpath_2d
    explicit_false = toolpath_run(
        [layer], mesh_data, overhang_map, replace(params, split_rings_into_runs=False)
    )[0].toolpath_2d

    assert len(default) == len(explicit_false)
    for a, b in zip(default, explicit_false):
        assert np.allclose(a, b)


def _dense_wavy_square_layer() -> Layer:
    """
    Variante DENSE de `_wavy_square_layer` (48 points, espacement 1.67mm -
    superieur a horizontal_radius=1.0, donc structurellement incapable de
    demontrer l'effet du mecanisme). Meme forme (carre, ondulation
    sinusoidale d'amplitude 1.0 autour du plan moyen), mais espacement proche
    de nozzle_width, pour que le test d'integration porte reellement sur
    l'effet du decoupage et non sur une resolution de fixture insuffisante.
    """
    n_per_side = 60
    t = np.linspace(0.0, 1.0, n_per_side, endpoint=False)
    side = 20.0
    half = side / 2.0
    xs = np.concatenate([-half + t * side, np.full(n_per_side, half),
                         half - t * side, np.full(n_per_side, -half)])
    ys = np.concatenate([np.full(n_per_side, -half), -half + t * side,
                         np.full(n_per_side, half), half - t * side])
    theta = np.arctan2(ys, xs)
    zs = 1.0 * np.sin(2.0 * theta)

    loop = np.stack([xs, ys, zs], axis=1)
    loop = np.vstack([loop, loop[0]])

    local_frame = (np.zeros(3), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    return Layer(id=0, contour_3d=[loop], region_ids=[0], local_frame=local_frame, isovalue=0.0)


def test_split_increases_polyline_count_without_losing_points(flat_box_layers):
    """
    Test d'integration : sur une couche NON PLANE dense, le decoupage doit
    produire STRICTEMENT PLUS de polylignes (chaque frontiere de troncon
    devient une vraie frontiere de segment, donc un vrai trajet a vide insere
    par step8._iter_schedule_points au lieu d'un cordon parasite extrude en
    travers de la piece) sans jamais perdre ni dupliquer un point.
    """
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _dense_wavy_square_layer()

    base_params = Step6Params(
        nozzle_width=0.4, wall_count=2, infill_density=0.3,
        out_of_plane_interpolation_power=2.0, split_rings_into_runs=False,
    )
    split_params = replace(base_params, split_rings_into_runs=True)

    before = toolpath_run([layer], mesh_data, overhang_map, base_params)[0].toolpath_2d
    after = toolpath_run([layer], mesh_data, overhang_map, split_params)[0].toolpath_2d

    assert len(after) > len(before)
    assert sum(len(p) for p in after) == sum(len(p) for p in before)


def test_split_reduces_unsupported_points_end_to_end(flat_box_layers):
    """
    Meme fixture, mesure du VRAI effet : en comptant les points deposes sans
    rien en dessous (avec un trajet a vide a chaque frontiere de polyligne -
    la geometrie que step8 produira reellement), le decoupage en troncons
    montants doit faire baisser le compte.
    """
    mesh_data, overhang_map, _layers = flat_box_layers
    layer = _dense_wavy_square_layer()
    horizontal_radius, vertical_height = 1.0, 1.5

    base_params = Step6Params(
        nozzle_width=0.4, wall_count=2, infill_density=0.3,
        out_of_plane_interpolation_power=2.0, split_rings_into_runs=False,
    )
    split_params = replace(base_params, split_rings_into_runs=True)

    def _count_unsupported(toolpath_2d):
        positions, extrusions = [], []
        for seg_idx, poly in enumerate(toolpath_2d):
            pts_3d = _lift_to_3d(poly, layer.local_frame)
            if seg_idx > 0:
                positions.append(pts_3d[0])
                extrusions.append(0.0)   # trajet a vide, comme step8
            for p in pts_3d:
                positions.append(p)
                extrusions.append(1.0)
        positions = np.array(positions)
        extrusions = np.array(extrusions)

        unsupported = 0
        for i in range(1, len(positions)):
            if extrusions[i] <= 0.0:
                continue
            pos = positions[i]
            prior = positions[:i][extrusions[:i] > 0.0]
            if len(prior) == 0:
                unsupported += 1
                continue
            dz = pos[2] - prior[:, 2]
            horiz = np.hypot(pos[0] - prior[:, 0], pos[1] - prior[:, 1])
            if not np.any((dz > 1e-9) & (dz <= vertical_height) & (horiz <= horizontal_radius)):
                unsupported += 1
        return unsupported

    before = _count_unsupported(toolpath_run([layer], mesh_data, overhang_map, base_params)[0].toolpath_2d)
    after = _count_unsupported(toolpath_run([layer], mesh_data, overhang_map, split_params)[0].toolpath_2d)

    assert after < before


def test_merge_short_runs_carries_dip_forward_not_backward():
    """Un point isole PLUS BAS que la fin du troncon precedent doit demarrer le
    troncon SUIVANT (il en devient le point le plus bas), jamais etre colle a la
    fin du precedent (il y serait depose apres plus haut que lui)."""
    #        idx: 0    1    2    3     4    5
    offsets = np.array([0.0, 1.0, 2.0, 1.9, 3.0, 4.0])
    runs = [(0, 3), (3, 4), (4, 6)]     # (3,4) = le micro-creux, 1 seul point

    merged = _merge_short_runs(runs, offsets)

    assert (0, 3) in merged, "le troncon montant ne doit pas absorber le creux"
    assert (3, 6) in merged, "le creux doit demarrer le troncon suivant"


def test_merge_short_runs_still_merges_backward_when_monotony_holds():
    """Si le point isole est PLUS HAUT que la fin du troncon precedent, la
    fusion arriere preserve la montee - comportement conserve."""
    offsets = np.array([0.0, 1.0, 2.0, 2.5, 1.0, 0.5])
    runs = [(0, 3), (3, 4), (4, 6)]

    merged = _merge_short_runs(runs, offsets)

    assert (0, 4) in merged, "2.5 >= 2.0 : la fusion arriere reste monotone"


def test_merge_short_runs_partitions_exactly_with_offsets():
    """Aucun point perdu ni duplique, quelle que soit la direction choisie."""
    offsets = np.array([0.0, 1.0, 0.9, 2.0, 3.0, 2.9, 4.0])
    runs = [(0, 2), (2, 3), (3, 5), (5, 6), (6, 7)]

    merged = _merge_short_runs(runs, offsets)

    assert merged[0][0] == 0
    assert merged[-1][1] == len(offsets)
    for (_, stop_a), (start_b, _) in zip(merged, merged[1:]):
        assert stop_a == start_b
    assert sum(stop - start for start, stop in merged) == len(offsets)


def test_merge_short_runs_legacy_path_unchanged():
    """`offsets=None` conserve strictement l'ancien comportement (fusion
    arriere systematique) - non-regression pour tout appelant sans altitude."""
    runs = [(0, 3), (3, 4), (4, 6)]
    assert _merge_short_runs(runs) == [(0, 4), (4, 6)]


def _square(side):
    from shapely.geometry import box
    return box(0.0, 0.0, side, side)


def test_infill_spacing_follows_nozzle_and_density():
    """Carre de 40 mm, buse 0,4, densite 0,2 : une ligne tous les 2 mm, soit
    environ 20 lignes. L'ancienne formule (fondee sur la diagonale de la
    couche) en donnait ~5 : c'est le defaut qui laissait les grandes couches
    presque vides (mesure sur les couches du Benchy)."""
    lines = _generate_infill(_square(40.0), 0.2, InfillPattern.LINES, nozzle_width=0.4)
    assert 18 <= len(lines) <= 22
    ys = sorted(float(line[0][1]) for line in lines)
    assert np.allclose(np.diff(ys), 2.0, atol=1e-6)


def test_infill_spacing_does_not_depend_on_layer_size():
    """Meme densite, couches de tailles differentes : meme espacement, donc
    un nombre de lignes proportionnel a la taille."""
    small = _generate_infill(_square(10.0), 0.2, InfillPattern.LINES, nozzle_width=0.4)
    large = _generate_infill(_square(80.0), 0.2, InfillPattern.LINES, nozzle_width=0.4)
    assert abs(len(large) / len(small) - 8.0) < 1.0


def test_infill_full_density_lines_are_one_nozzle_apart():
    lines = _generate_infill(_square(10.0), 1.0, InfillPattern.LINES, nozzle_width=0.4)
    ys = sorted(float(line[0][1]) for line in lines)
    assert np.allclose(np.diff(ys), 0.4, atol=1e-6)


def _ring(xmin, ymin, xmax, ymax):
    return np.array([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]], dtype=float)


def _ring_bounds(ring):
    return ring[:, 0].min(), ring[:, 1].min(), ring[:, 0].max(), ring[:, 1].max()


def test_hole_receives_its_contour_and_walls():
    """Carre de 20 mm perce d'un trou de 8 mm, 2 parois : le trou doit etre
    borde par son contour (n=0) ET par une paroi a chaque niveau, placee
    cote matiere (trou elargi de 0,4 puis 0,8 mm). Avant la correction, seul
    le bord exterieur de chaque offset etait trace : le trou n'avait rien."""
    outer, hole = _ring(0, 0, 20, 20), _ring(6, 6, 14, 14)
    walls, _ = _offset_walls(outer, [hole], nozzle_width=0.4, wall_count=2, degenerate_area_ratio=0.01)

    around_hole = [w for w in walls if _ring_bounds(w)[2] - _ring_bounds(w)[0] < 12.0]
    widths = sorted(round(_ring_bounds(w)[2] - _ring_bounds(w)[0], 3) for w in around_hole)
    assert widths == [8.0, 8.8, 9.6], f"bords du trou attendus a 8 / 8,8 / 9,6 mm, obtenus {widths}"


def test_split_offset_keeps_every_piece_with_walls_and_infill():
    """Haltere : deux carres de 10 mm relies par un col de 0,6 mm. Des la
    premiere paroi (0,4 mm), le col disparait et la zone se coupe en deux.
    Chaque moitie doit garder ses parois et son remplissage ; avant la
    correction, seul le plus grand morceau etait garde."""
    outer = np.array([[0, 0], [10, 0], [10, 4.7], [14, 4.7], [14, 0], [24, 0],
                      [24, 10], [14, 10], [14, 5.3], [10, 5.3], [10, 10], [0, 10]], dtype=float)
    walls, innermost = _offset_walls(outer, [], nozzle_width=0.4, wall_count=2, degenerate_area_ratio=0.01)

    inner_walls = walls[:-1]              # le dernier est le contour d'origine
    left = [w for w in inner_walls if _ring_bounds(w)[2] < 12.0]
    right = [w for w in inner_walls if _ring_bounds(w)[0] > 12.0]
    assert len(left) == 2 and len(right) == 2, "2 parois attendues dans chaque moitie"

    infill = _generate_infill(innermost, 0.2, InfillPattern.LINES, nozzle_width=0.4)
    xs = np.concatenate([line[:, 0] for line in infill])
    assert (xs < 12.0).any() and (xs > 12.0).any(), "remplissage attendu dans les deux moities"


def _segments(polyline):
    return np.diff(np.asarray(polyline)[:, :2], axis=0)


def test_zigzag_links_a_convex_zone_into_one_continuous_path():
    """Carre de 40 mm : un seul zigzag, qui contient toutes les lignes de
    balayage (~20, une tous les 2 mm) reliees par des raccords courts."""
    chains = _generate_infill(_square(40.0), 0.2, InfillPattern.CROSSED_ZIGZAG, nozzle_width=0.4)
    assert len(chains) == 1
    seg = _segments(chains[0])
    sweeps = np.abs(seg[:, 0]) > 10.0          # lignes de balayage (selon u)
    assert 18 <= sweeps.sum() <= 22
    assert np.all(np.linalg.norm(seg[~sweeps], axis=1) <= 2.5 * 2.0 + 1e-9)
    # le zigzag change de sens a chaque ligne
    assert np.all(np.sign(seg[sweeps, 0])[1:] != np.sign(seg[sweeps, 0])[:-1])


def test_zigzag_never_jumps_across_a_hole_or_a_gap():
    """Carre perce d'un grand trou, et haltere coupee en deux par la paroi :
    aucun raccord du zigzag ne traverse le vide - chaque segment reste dans
    la zone a remplir."""
    from shapely.geometry import LineString, box

    holed = box(0, 0, 30, 30).difference(box(8, 8, 22, 22))
    split = box(0, 0, 10, 10).union(box(14, 0, 24, 10))
    for zone in (holed, split):
        chains = _generate_infill(zone, 0.2, InfillPattern.CROSSED_ZIGZAG, nozzle_width=0.4)
        assert len(chains) >= 2, "une zone non convexe doit donner plusieurs zigzags"
        inside = zone.buffer(1e-6)
        for chain in chains:
            for a, b in zip(chain[:-1], chain[1:]):
                assert inside.covers(LineString([a, b]))


def test_grid_crosses_two_sweeps_in_the_same_layer():
    """GRID etait recu et ignore (lignes simples). Il produit maintenant deux
    balayages a 90 deg, chacun a la moitie de la densite."""
    lines = _generate_infill(_square(40.0), 0.2, InfillPattern.GRID, nozzle_width=0.4)
    dirs = np.array([line[-1] - line[0] for line in lines])
    along_u = np.abs(dirs[:, 0]) > np.abs(dirs[:, 1])
    assert 8 <= along_u.sum() <= 11 and 8 <= (~along_u).sum() <= 11


def _flat_square_layer(layer_id, isovalue, side=20.0):
    half = side / 2.0
    loop = np.array([[-half, -half, 0.0], [half, -half, 0.0], [half, half, 0.0],
                     [-half, half, 0.0], [-half, -half, 0.0]])
    frame = (np.zeros(3), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    return Layer(id=layer_id, contour_3d=[loop], region_ids=[0], local_frame=frame, isovalue=isovalue)


def _zigzag_axis(layer):
    """0 si le zigzag de la couche balaie selon u, 1 selon v."""
    chain = max((p for p in layer.toolpath_2d if not np.allclose(p[0], p[-1])), key=len)
    seg = np.abs(_segments(chain))
    return int(seg[:, 1].sum() > seg[:, 0].sum())


def test_zigzag_alternates_between_stacked_layers(flat_box_layers):
    """Lignes croisees-alternees : deux couches EMPILEES (isovaleurs
    consecutives) se croisent a 90 deg. La parite suit le rang de
    l'isovaleur, pas la position dans la liste (l'ordre d'impression en
    arbre fait se succeder des couches de branches differentes), et deux
    morceaux d'une meme isovaleur gardent la meme orientation."""
    mesh_data, overhang_map, _layers = flat_box_layers
    layers = [_flat_square_layer(0, 0.0), _flat_square_layer(1, 0.4),
              _flat_square_layer(2, 0.2), _flat_square_layer(3, 0.2)]
    params = Step6Params(nozzle_width=0.4, wall_count=2, infill_density=0.2)
    assert params.infill_pattern is InfillPattern.CROSSED_ZIGZAG, "motif par defaut"
    result = toolpath_run(layers, mesh_data, overhang_map, params)
    axes = [_zigzag_axis(layer) for layer in result]
    assert axes == [0, 0, 1, 1]


def test_zigzag_does_not_retrace_the_innermost_wall():
    """La zone passee a _generate_infill est le polygone de la paroi la plus
    interieure : son bord EST l'axe de cette paroi. Le zigzag doit s'en tenir
    a une demi-buse, sinon ses raccords, qui suivent le bord, repasseraient
    sur la paroi."""
    from shapely.geometry import Point

    zone = _square(20.0)
    chains = _generate_infill(zone, 0.2, InfillPattern.CROSSED_ZIGZAG, nozzle_width=0.4)
    points = np.vstack(chains)
    distances = [zone.exterior.distance(Point(p)) for p in points]
    assert min(distances) >= 0.2 - 1e-9


def test_split_ring_shared_vertices_extrude_every_segment():
    """
    `continuous_rings` (2026-10-01) : avec des troncons a indices DISJOINTS, le
    segment qui relie deux troncons n'est jamais extrude (il devient un trajet
    a vide) - 7,7 % de la longueur perdue sur death.stl. Avec le sommet de
    jonction partage, chaque segment de l'anneau appartient a un troncon, et
    chaque troncon reste montant.
    """
    ring = _wavy_ring(n=400, radius=25.0, amplitude=2.0, cycles=6)
    runs = _split_ring_into_runs(ring, share_vertices=True)
    assert len(runs) > 1

    def edges(poly):
        keys = set()
        for a, b in zip(poly[:-1], poly[1:]):
            keys.add(frozenset((tuple(np.round(a, 9)), tuple(np.round(b, 9)))))
        return keys

    covered = set().union(*(edges(r) for r in runs))
    assert edges(ring) <= covered
    for run in runs:
        assert run[-1, 2] >= run[0, 2]
    # sans partage, des segments manquent : le test verifie bien quelque chose
    disjoint = set().union(*(edges(r) for r in _split_ring_into_runs(ring)))
    assert not edges(ring) <= disjoint


def test_continuous_rings_closes_original_contour(flat_box_layers):
    """Le contour d'origine n=0 sort OUVERT de la classification : sans
    point de fermeture, son dernier segment n'est jamais trace. Avec
    `continuous_rings`, il est ferme (premier point repete)."""
    mesh_data, overhang_map, layers = flat_box_layers
    base = Step6Params(nozzle_width=0.4, wall_count=1, infill_density=0.2)
    closed = Step6Params(nozzle_width=0.4, wall_count=1, infill_density=0.2, continuous_rings=True)

    def original(result):
        found = []
        for layer in result:
            contour_2d = _project_to_local_2d(layer.contour_3d[0][:-1], layer.local_frame)
            for polyline in layer.toolpath_2d or []:
                if len(polyline) >= len(contour_2d) and np.allclose(polyline[: len(contour_2d), :2], contour_2d):
                    found.append((polyline, len(contour_2d)))
        return found

    open_ones = original(toolpath_run(layers, mesh_data, overhang_map, base))
    closed_ones = original(toolpath_run(layers, mesh_data, overhang_map, closed))
    assert open_ones and closed_ones
    assert all(len(p) == n for p, n in open_ones)
    assert all(len(p) == n + 1 and np.allclose(p[0], p[-1]) for p, n in closed_ones)


def test_monotonic_runs_ignores_oscillation_below_min_amplitude():
    # Bruit de 4 um autour d'une couche plate : une seule course, pas une coupe par oscillation.
    offsets = 0.002 * np.sin(np.arange(200) * 1.7)
    assert len(_monotonic_runs(offsets)) > 50                      # defaut historique : coupe partout
    assert _monotonic_runs(offsets, min_amplitude=0.05) == [(0, 200)]


def test_monotonic_runs_still_cuts_a_real_reversal_with_min_amplitude():
    offsets = np.concatenate([np.linspace(0, 1, 50), np.linspace(1, 0, 50)[1:]])
    assert _monotonic_runs(offsets, min_amplitude=0.05) == [(0, 50), (50, 99)]


def _spiral_params(**overrides):
    base = dict(nozzle_width=0.4, wall_count=2, infill_density=0.2,
                infill_pattern=InfillPattern.OSCILLATING_SPIRAL)
    base.update(overrides)
    return Step6Params(**base)


def _disk(radius=13.0):
    from shapely.geometry import Point
    return Point(0.0, 0.0).buffer(radius, resolution=64)


def _raw_spiral(rho, pitch=3.1, amplitude=0.48, seed=0, rank=0):
    rng = np.random.default_rng([seed, rank])
    return _oscillating_spiral_polyline(
        _SpiralShape(np.zeros(2), 0.0, 13.0, pitch, amplitude, 1.6, rho, 0.4), rng,
        first_sign=1.0, angle_offset=0.0,
    )


@pytest.mark.parametrize("rho", [0.0, 0.5, 1.0])
def test_oscillating_spiral_never_goes_backward_and_never_crosses_itself(rho):
    """Le plateau ne recule jamais (theta croissant) et la trajectoire ne se
    recoupe pas : deux pics purement radiaux consecutifs repasseraient sur la
    meme ligne (defaut trouve a la mise au point)."""
    from shapely.geometry import LineString
    line = _raw_spiral(rho)
    angles = np.unwrap(np.arctan2(line[:, 1], line[:, 0]))
    assert np.all(np.diff(angles) >= -1e-9)
    assert LineString(line).is_simple


def test_oscillating_spiral_has_radial_peaks_that_do_not_advance_only_when_randomized():
    def stationary(line):
        radii = np.hypot(line[:, 0], line[:, 1])
        angles = np.unwrap(np.arctan2(line[:, 1], line[:, 0]))
        return int(np.sum((np.abs(np.diff(angles)) < 1e-9) & (np.abs(np.diff(radii)) > 1e-3)))
    assert stationary(_raw_spiral(0.0)) == 0
    assert stationary(_raw_spiral(1.0)) > 0


def test_oscillating_spiral_amplitude_is_capped_so_neighbouring_turns_never_touch():
    pitch, width = 3.1, 0.4
    line = _raw_spiral(1.0, pitch=pitch, amplitude=50.0)          # seuil demande absurde
    angles = np.unwrap(np.arctan2(line[:, 1], line[:, 0]))
    deviation = np.hypot(line[:, 0], line[:, 1]) - pitch * angles / (2 * np.pi)
    assert np.abs(deviation).max() <= 0.5 * (pitch - width) + 1e-9


@pytest.mark.parametrize("rho", [0.0, 0.5, 1.0])
def test_oscillating_spiral_infill_hits_the_requested_density_inside_the_zone(rho):
    from shapely.geometry import LineString
    params = _spiral_params(spiral_randomization_threshold=rho)
    chains = _generate_spiral_infill(_disk(), params, 0, np.zeros(2))
    zone = _disk().buffer(-0.2, join_style="mitre")
    length = sum(float(np.linalg.norm(np.diff(c, axis=0), axis=1).sum()) for c in chains)
    assert length * 0.4 / zone.area == pytest.approx(0.2, rel=0.12)
    assert all(zone.buffer(1e-3).covers(LineString(c)) for c in chains)


def test_oscillating_spiral_is_reproducible_and_alternates_between_layers():
    params = _spiral_params()
    a = _generate_spiral_infill(_disk(), params, 0, np.zeros(2))
    b = _generate_spiral_infill(_disk(), params, 0, np.zeros(2))
    assert len(a) == len(b) and all(np.array_equal(x, y) for x, y in zip(a, b))
    other_seed = _generate_spiral_infill(_disk(), _spiral_params(spiral_seed=7), 0, np.zeros(2))
    assert not np.array_equal(a[0], other_seed[0])
    odd = _generate_spiral_infill(_disk(), params, 1, np.zeros(2))
    assert not np.array_equal(a[0][: len(odd[0])], odd[0][: len(a[0])])


def test_oscillating_spiral_layers_share_the_same_median_line():
    """Les couches se superposent sur la meme ligne mediane : l'ecart radial
    de chaque couche a la spirale de base reste borne (plus de decalage d'un
    demi-pas entre couches)."""
    pitch = 2.0
    params = _spiral_params(spiral_randomization_threshold=0.7)
    for rank in (0, 1):
        for chain in _generate_spiral_infill(_disk(), params, rank, np.zeros(2)):
            angles = np.unwrap(np.arctan2(chain[:, 1], chain[:, 0]))
            radius = np.hypot(chain[:, 0], chain[:, 1])
            offset = np.mod(radius - pitch * angles / (2 * np.pi) + pitch / 2, pitch) - pitch / 2
            assert np.abs(offset).max() < 0.5 * pitch


def test_oscillating_spiral_connectors_never_make_the_plate_go_backward():
    from pipeline.step6_toolpath2d import _angle_nondecreasing
    for rank in (0, 1):
        for chain in _generate_spiral_infill(_disk(), _spiral_params(spiral_randomization_threshold=1.0),
                                             rank, np.zeros(2)):
            assert _angle_nondecreasing(chain, np.zeros(2))


def test_plate_axis_in_local_2d_follows_the_layer_plane():
    frame = (np.array([5.0, 2.0, 1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    assert np.allclose(_plate_axis_in_local_2d(frame, (0.0, 0.0, 0.0)), [-5.0, -2.0])
    # plan incline de 45 deg autour de x, passant par l'origine : l'axe vertical
    # x = y = 0 le coupe a l'origine, donc (0, 0) en local
    tilt = (np.zeros(3), np.array([1.0, 0.0, 0.0]), np.array([0.0, np.cos(np.pi / 4), np.sin(np.pi / 4)]))
    assert np.allclose(_plate_axis_in_local_2d(tilt, (0.0, 0.0, 3.0)), [0.0, 0.0], atol=1e-9)
    # un axe decale de 1 mm en x se projette a 1 mm sur u
    assert _plate_axis_in_local_2d(tilt, (1.0, 0.0, 0.0))[0] == pytest.approx(1.0)


def test_oscillating_spiral_through_run_is_monotone_around_the_axis_and_alternates(flat_box_layers):
    mesh_data, overhang_map, _layers = flat_box_layers
    layers = [_flat_square_layer(0, 0.0), _flat_square_layer(1, 0.2)]
    params = _spiral_params(spiral_axis_origin=(0.0, 0.0, 0.0), spiral_randomization_threshold=0.5)
    result = toolpath_run(layers, mesh_data, overhang_map, params)
    infills = []
    for layer in result:
        chains = [p for p in layer.toolpath_2d if not np.allclose(p[0, :2], p[-1, :2])]
        assert chains, "le remplissage en spirale doit exister"
        for chain in chains:
            angles = np.unwrap(np.arctan2(chain[:, 1], chain[:, 0]))
            assert np.all(np.diff(angles) >= -1e-9)
        infills.append(max(chains, key=len))
    n = min(len(infills[0]), len(infills[1]))
    assert not np.allclose(infills[0][:n, :2], infills[1][:n, :2])
