"""Etape 6 : motifs de remplissage (lignes, zigzag, spirale oscillante).

Découpé de `step6_toolpath2d.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step6_toolpath2d`.
"""

from __future__ import annotations

import numpy as np
import shapely
import shapely.affinity
from shapely.geometry import LineString, Polygon
from core.mesh_types import Layer
from pipeline.step6_geometry import (  # noqa: F401
    _polygon_pieces,
)
from pipeline.step6_params import (  # noqa: F401
    InfillPattern,
)


_MAX_INFILL_LINES = 5000


# longueur maximale d'un raccord de zigzag, en espacements de balayage
_ZIGZAG_MAX_CONNECTOR_SPACINGS = 2.5


def _generate_infill(
    boundary: Polygon,
    density: float,
    pattern: InfillPattern,
    nozzle_width: float,
    sweep_angle_rad: float = 0.0,
) -> list[np.ndarray]:
    """
    Genere l'infill a l'interieur de `boundary` par balayage de lignes
    paralleles clippees par le polygone (shapely.intersection), selon
    `pattern` :
      - LINES : lignes independantes (une polyligne de 2 points chacune),
        orientees selon `sweep_angle_rad` ;
      - GRID : deux balayages croises a 90 deg, chacun a densite / 2 ;
      - CROSSED_ZIGZAG : les lignes d'un meme balayage sont reliees bout a
        bout en zigzag (`_connect_zigzag`) ; l'alternance entre couches est
        portee par l'appelant via `sweep_angle_rad`
        (`_infill_angle_for_layer`).

    L'espacement entre lignes d'un balayage vaut nozzle_width / density :
    1.0 = lignes jointives (une largeur de buse), 0.2 = une ligne toutes les
    cinq largeurs de buse.

    Correction (2026-09-26) : l'espacement etait derive de la TAILLE de la
    couche (`diagonal / density * 0.05`), jamais de la buse - si bien qu'une
    couche de 35 mm a 20 % recevait une ligne tous les 8,75 mm au lieu de
    2 mm, et que plus une couche etait grande, plus elle etait creuse.
    Mesure sur les couches du Benchy (43 mm de remplissage au lieu
    d'environ 160 sur une couche de 424 mm2). Et `pattern` etait recu sans
    jamais etre lu : GRID produisait des lignes simples.

    sweep_angle_rad force l'orientation du balayage (ponts : traverser
    perpendiculairement a la portee, cf. Decision 2 du plan ; couches
    normales en CROSSED_ZIGZAG : 0 ou pi/2 selon la couche).
    """
    if boundary is None or boundary.is_empty or density <= 0:
        return []

    if pattern is InfillPattern.GRID:
        half = min(density, 1.0) / 2.0
        return (
            _generate_infill(boundary, half, InfillPattern.LINES, nozzle_width, sweep_angle_rad)
            + _generate_infill(
                boundary, half, InfillPattern.LINES, nozzle_width, sweep_angle_rad + np.pi / 2.0
            )
        )

    if pattern is InfillPattern.CROSSED_ZIGZAG:
        # `boundary` est le polygone de la paroi la plus interieure : son bord
        # EST l'axe de cette paroi. Les raccords du zigzag suivant ce bord, ils
        # repasseraient sur la paroi (double depot) : on les en ecarte d'une
        # demi-largeur de buse, cordons jointifs. Une zone plus etroite
        # qu'une buse ne recoit pas de remplissage (les parois s'y touchent).
        zone = boundary.buffer(-nozzle_width / 2.0, join_style="mitre")
        if zone.is_empty:
            return []
        scans, spacing = _clipped_scanlines(zone, density, nozzle_width, sweep_angle_rad)
        return _connect_zigzag(scans, zone, spacing)

    scans, spacing = _clipped_scanlines(boundary, density, nozzle_width, sweep_angle_rad)
    return [segment for scan in scans for segment in scan]


def _clipped_scanlines(
    boundary: Polygon, density: float, nozzle_width: float, sweep_angle_rad: float
) -> tuple[list[list[np.ndarray]], float]:
    """
    Lignes de balayage clippees par `boundary`, regroupees PAR ligne de
    balayage (dans l'ordre de balayage) et, au sein d'une ligne, triees et
    orientees dans le sens du balayage. Rend aussi l'espacement retenu.
    """
    minx, miny, maxx, maxy = boundary.bounds
    diagonal = float(np.hypot(maxx - minx, maxy - miny))
    if diagonal <= 0:
        return [], 0.0

    spacing = nozzle_width / min(density, 1.0)
    # borne basse pour eviter un nombre de lignes explosif (buse ou densite
    # aberrante) : le balayage couvre 2 x diagonal, donc au plus
    # _MAX_INFILL_LINES lignes par couche.
    spacing = max(spacing, 2.0 * diagonal / _MAX_INFILL_LINES, 1e-6)

    center = np.array([(minx + maxx) / 2.0, (miny + maxy) / 2.0])
    half_extent = diagonal  # sur-couvre largement pour rester valide apres rotation

    offsets = np.arange(-half_extent, half_extent, spacing)
    cos_a, sin_a = np.cos(sweep_angle_rad), np.sin(sweep_angle_rad)
    direction = np.array([cos_a, sin_a])
    perpendicular = np.array([-sin_a, cos_a])

    scans: list[list[np.ndarray]] = []
    for offset in offsets:
        base = center + perpendicular * offset
        p0 = base - direction * half_extent
        p1 = base + direction * half_extent
        clipped = LineString([p0, p1]).intersection(boundary)
        if clipped.is_empty:
            continue
        parts = clipped.geoms if hasattr(clipped, "geoms") else [clipped]
        segments = []
        for part in parts:
            if part.geom_type != "LineString":
                continue
            coords = np.asarray(part.coords)
            if len(coords) < 2:
                continue
            if np.dot(coords[-1] - coords[0], direction) < 0:
                coords = coords[::-1]
            segments.append(coords)
        if segments:
            segments.sort(key=lambda c: float(np.dot(c[0], direction)))
            scans.append(segments)
    return scans, spacing


def _connect_zigzag(
    scans: list[list[np.ndarray]], boundary: Polygon, spacing: float
) -> list[np.ndarray]:
    """
    Relie les segments de lignes de balayage SUCCESSIVES en zigzag : depuis
    la fin d'un segment, on enchaine sur l'extremite la plus proche d'un
    segment libre de la ligne suivante, que l'on parcourt alors en sens
    inverse. Le raccord SUIT LE BORD de la zone entre les deux extremites
    (`_boundary_connector`), comme dans un slicer classique ; a defaut (les
    deux extremites ne sont pas sur un meme bord), un raccord droit est
    accepte s'il reste dans `boundary`. Dans les deux cas il mesure au plus
    `_ZIGZAG_MAX_CONNECTOR_SPACINGS` espacements - au-dela il longerait un
    bord presque parallele au balayage et deposerait une ligne en trop.
    Sinon le zigzag s'arrete et un nouveau commence (une zone concave ou
    trouee donne donc plusieurs zigzags, relies par des trajets a vide).

    Un raccord droit seul ne suffit pas : mesure sur les couches Bird du
    Benchy, dont les bords sont des polygones irreguliers issus du maillage,
    le raccord droit sortait legerement de la zone a presque chaque ligne et
    la mediane tombait a UNE ligne par zigzag.
    """
    from shapely.prepared import prep

    inside = prep(boundary.buffer(max(1e-6, spacing * 1e-4)))
    max_connector = _ZIGZAG_MAX_CONNECTOR_SPACINGS * spacing
    rings = [
        LineString(ring.coords)
        for piece in _polygon_pieces(boundary)
        for ring in (piece.exterior, *piece.interiors)
    ]
    used = [[False] * len(scan) for scan in scans]
    chains: list[np.ndarray] = []

    for k, scan in enumerate(scans):
        for j, segment in enumerate(scan):
            if used[k][j]:
                continue
            used[k][j] = True
            parts = [segment]
            end = segment[-1]
            current = k
            while current + 1 < len(scans):
                best = None
                for jj, candidate in enumerate(scans[current + 1]):
                    if used[current + 1][jj]:
                        continue
                    for reverse, start in ((False, candidate[0]), (True, candidate[-1])):
                        distance = float(np.linalg.norm(start - end))
                        if best is None or distance < best[0]:
                            best = (distance, jj, reverse)
                if best is None or best[0] > max_connector:
                    break
                distance, jj, reverse = best
                candidate = scans[current + 1][jj]
                start = candidate[-1] if reverse else candidate[0]
                connector = _boundary_connector(end, start, rings, max_connector)
                if connector is None:
                    if distance > 0 and not inside.covers(LineString([end, start])):
                        break
                    connector = np.empty((0, 2))
                used[current + 1][jj] = True
                oriented = candidate[::-1] if reverse else candidate
                parts.append(connector)
                parts.append(oriented)
                end = oriented[-1]
                current += 1
            chains.append(np.vstack(parts))
    return chains


def _boundary_connector(
    a: np.ndarray, b: np.ndarray, rings: list[LineString], max_length: float
) -> np.ndarray | None:
    """
    Points intermediaires du plus court chemin LE LONG D'UN BORD (anneau
    exterieur ou trou) entre deux points situes sur ce bord, extremites
    exclues. None si aucun anneau ne porte les deux points, ou si le chemin
    depasse `max_length`.
    """
    from shapely.geometry import Point
    from shapely.ops import substring

    pa, pb = Point(a), Point(b)
    tolerance = 1e-6 * max(1.0, max_length)
    for ring in rings:
        if ring.distance(pa) > tolerance or ring.distance(pb) > tolerance:
            continue
        length = ring.length
        sa, sb = ring.project(pa), ring.project(pb)
        forward = (sb - sa) % length
        backward = (sa - sb) % length
        if min(forward, backward) > max_length:
            return None

        def along(s0, s1):
            if s1 >= s0:
                return np.asarray(substring(ring, s0, s1).coords)
            return np.vstack([np.asarray(substring(ring, s0, length).coords),
                              np.asarray(substring(ring, 0.0, s1).coords)])

        path = along(sa, sb) if forward <= backward else along(sb, sa)[::-1]
        return path[1:-1] if len(path) > 2 else np.empty((0, 2))
    return None


_SPIRAL_SAMPLE_STEP_FACTOR = 0.5       # pas d'echantillonnage = 0,5 buse


_SPIRAL_MIN_RADIUS_PITCHES = 0.5       # rayon de depart : pres de l'axe, C tourne tres vite


_SPIRAL_DENSITY_TOLERANCE = 0.05


_SPIRAL_DENSITY_ITERATIONS = 6


_SPIRAL_STOP_PROBABILITY = 0.3         # part (x randomisation) de pics purement radiaux


def _layer_rank(layer: Layer, isovalue_rank: dict) -> int:
    if layer.isovalue is not None and layer.isovalue in isovalue_rank:
        return int(isovalue_rank[layer.isovalue])
    return int(layer.id)


def _plate_axis_in_local_2d(
    local_frame: tuple[np.ndarray, np.ndarray, np.ndarray], axis_origin: tuple[float, float, float]
) -> np.ndarray:
    """
    Intersection de l'axe du plateau (vertical, passant par `axis_origin`) avec
    le plan de la couche, en coordonnees locales (u, v). Couche quasi
    verticale (pas d'intersection stable) : projection du point lui-meme.
    """
    origin, axis_u, axis_v = (np.asarray(a, dtype=float) for a in local_frame)
    normal = np.cross(axis_u, axis_v)
    normal /= max(float(np.linalg.norm(normal)), 1e-12)
    a0 = np.asarray(axis_origin, dtype=float)
    denom = float(normal[2])
    point = a0 if abs(denom) < 0.1 else a0 + np.array([0.0, 0.0, float(normal @ (origin - a0)) / denom])
    d = point - origin
    return np.array([float(d @ axis_u), float(d @ axis_v)])


def _oscillating_spiral_polyline(
    center: np.ndarray, r_min_zone: float, r_max: float, pitch: float, amplitude: float,
    wavelength: float, randomization: float, nozzle_width: float, rng: np.random.Generator,
    first_sign: float, angle_offset: float,
) -> np.ndarray:
    """
    Spirale r = pitch * theta / 2pi, oscillante, en coordonnees polaires
    (theta, r) autour de `center` : une suite de points de controle dont theta
    ne decroit JAMAIS (deux points au meme theta = pic qui n'avance pas,
    mouvement purement radial). Entre deux points, la deviation radiale suit
    une droite (pic net) ou un cosinus (courbe), toujours monotone en theta.
    `angle_offset` pivote toute la spirale (pi = decalage d'un demi-pas).
    """
    rho = float(np.clip(randomization, 0.0, 1.0))
    cap = max(0.0, 0.5 * (pitch - nozzle_width))   # deux tours voisins restent a >= 1 buse
    amp = min(amplitude, cap)
    step = nozzle_width * _SPIRAL_SAMPLE_STEP_FACTOR
    r_start = max(_SPIRAL_MIN_RADIUS_PITCHES * pitch, r_min_zone - 2.0 * amp)
    theta = 2.0 * np.pi * r_start / pitch
    theta_end = 2.0 * np.pi * (r_max + 2.0 * amp) / pitch
    base = lambda th: pitch * th / (2.0 * np.pi)

    ctrl_theta = [theta]
    ctrl_dev = [0.0]
    ctrl_curve = [False]
    sign = first_sign
    min_half = nozzle_width
    k = 0
    prev_stop = False
    while theta < theta_end and amp > 0.0:
        r_eff = max(base(theta), r_start)
        half = max(min_half, 0.5 * wavelength * (1.0 + rho * 0.6 * rng.uniform(-1.0, 1.0)))
        # jamais deux pics radiaux de suite : la tete repasserait sur la meme ligne
        stop = k > 0 and not prev_stop and rng.random() < _SPIRAL_STOP_PROBABILITY * rho
        prev_stop = stop
        a_k = float(np.clip(amp * (1.0 + rho * rng.uniform(-0.6, 0.0)), 0.25 * amp, amp))
        curve = bool(rho > 0.0 and rng.random() < 0.5)
        theta = theta if stop else theta + half / r_eff
        ctrl_theta.append(theta)
        ctrl_dev.append(sign * a_k)
        ctrl_curve.append(curve)
        sign = -sign
        k += 1
    if amp <= 0.0:
        ctrl_theta.append(theta_end)
        ctrl_dev.append(0.0)
        ctrl_curve.append(False)

    thetas: list[np.ndarray] = [np.array([ctrl_theta[0]])]
    devs: list[np.ndarray] = [np.array([ctrl_dev[0]])]
    for i in range(1, len(ctrl_theta)):
        t0, t1 = ctrl_theta[i - 1], ctrl_theta[i]
        d0, d1 = ctrl_dev[i - 1], ctrl_dev[i]
        r_eff = max(base(t1), r_start)
        length = max((t1 - t0) * r_eff, abs(d1 - d0), 1e-9)
        n = max(1, int(np.ceil(length / step)))
        s = np.linspace(0.0, 1.0, n + 1)[1:]
        thetas.append(t0 + (t1 - t0) * s)
        ease = (1.0 - np.cos(np.pi * s)) / 2.0 if ctrl_curve[i] else s
        devs.append(d0 + (d1 - d0) * ease)
    th = np.concatenate(thetas)
    r = base(th) + np.concatenate(devs)
    ang = th + angle_offset
    return np.column_stack([center[0] + r * np.cos(ang), center[1] + r * np.sin(ang)])


def _angle_nondecreasing(points: np.ndarray, center: np.ndarray, tolerance: float = 1e-9) -> bool:
    """Angle polaire autour de `center` (deroule) jamais decroissant le long de `points`."""
    d = points - center
    if np.any(np.hypot(d[:, 0], d[:, 1]) < 1e-9):
        return True            # passage par l'axe : angle indetermine, rien a verifier
    angles = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
    return bool(np.all(np.diff(angles) >= -tolerance))


def _generate_spiral_infill(
    boundary: Polygon, params: "Step6Params", rank: int, axis_local: np.ndarray | None
) -> list[np.ndarray]:
    """
    Remplissage OSCILLATING_SPIRAL d'une zone. Toutes les couches suivent la
    MEME spirale mediane (elles se superposent) ; seules les asperites
    different : oscillation en opposition de phase une couche sur deux et
    tirages differents (graine + rang). Le pas est corrige
    d'apres la longueur REELLEMENT deposee pour tenir `infill_density`.
    Les raccords entre morceaux suivent le bord (comme le zigzag) et ne sont
    pas garantis monotones en theta.
    """
    from shapely.prepared import prep

    if boundary is None or boundary.is_empty or params.infill_density <= 0:
        return []
    width = params.nozzle_width
    zone = boundary.buffer(-width / 2.0, join_style="mitre")
    if zone.is_empty or zone.area <= 0:
        return []
    density = min(params.infill_density, 1.0)
    coords_all = np.vstack([np.asarray(piece.exterior.coords) for piece in _polygon_pieces(zone)])
    center = np.asarray(axis_local, dtype=float) if axis_local is not None else np.array(zone.centroid.coords[0])
    radii = np.linalg.norm(coords_all - center, axis=1)
    r_max = float(radii.max())
    r_min_zone = 0.0 if zone.contains(shapely.Point(center)) else float(zone.distance(shapely.Point(center)))
    minx, miny, maxx, maxy = zone.bounds
    diagonal = float(np.hypot(maxx - minx, maxy - miny))
    min_pitch = max(2.0 * diagonal / _MAX_INFILL_LINES, 1e-6)

    pitch = max(width / density, min_pitch)
    odd = rank % 2 == 1
    best: LineString | None = None
    # amplitude et longueur d'onde fixees UNE fois (a partir du pas nominal) :
    # sinon le pas corrige les ferait varier et la correction de densite ne
    # convergerait pas.
    amplitude = (
        params.spiral_alternation_threshold
        if params.spiral_alternation_threshold is not None
        else 0.6 * max(0.0, 0.5 * (pitch - width))
    )
    wavelength = params.spiral_peak_spacing if params.spiral_peak_spacing is not None else 4.0 * width
    for _ in range(_SPIRAL_DENSITY_ITERATIONS):
        rng = np.random.default_rng([int(params.spiral_seed) & 0xFFFFFFFF, int(rank)])
        line = _oscillating_spiral_polyline(
            center, r_min_zone, r_max, pitch, amplitude, wavelength,
            params.spiral_randomization_threshold, width, rng,
            first_sign=-1.0 if odd else 1.0, angle_offset=0.0,
        )
        clipped = LineString(line).intersection(zone)
        length = float(clipped.length)
        best = clipped
        achieved = length * width / zone.area
        ratio = achieved / density if density > 0 else 1.0
        if abs(ratio - 1.0) <= _SPIRAL_DENSITY_TOLERANCE or length <= 0:
            break
        pitch = max(pitch * float(np.clip(ratio, 0.5, 2.0)), min_pitch)
    if best is None or best.is_empty:
        return []

    pieces = [np.asarray(g.coords) for g in (best.geoms if hasattr(best, "geoms") else [best])
              if g.geom_type == "LineString" and len(g.coords) >= 2]
    if not pieces:
        return []
    # ordre le long de la spirale : position du premier point sur la ligne
    line_geom = LineString(line)
    order = shapely.line_locate_point(line_geom, shapely.points(np.array([p[0] for p in pieces])))
    pieces = [pieces[i] for i in np.argsort(order, kind="stable")]

    inside = prep(zone.buffer(max(1e-6, pitch * 1e-4)))
    max_connector = _ZIGZAG_MAX_CONNECTOR_SPACINGS * pitch
    rings = [LineString(ring.coords) for piece in _polygon_pieces(zone) for ring in (piece.exterior, *piece.interiors)]
    chains: list[np.ndarray] = []
    current = pieces[0]
    for nxt in pieces[1:]:
        end, start = current[-1], nxt[0]
        gap = float(np.linalg.norm(start - end))
        connector = _boundary_connector(end, start, rings, max_connector) if gap <= max_connector else None
        if connector is None and gap > 0 and (gap > max_connector or not inside.covers(LineString([end, start]))):
            chains.append(current)
            current = nxt
            continue
        joint = np.vstack([end[None, :], connector if connector is not None else np.empty((0, 2)), start[None, :]])
        if not _angle_nondecreasing(joint, center):
            # un raccord qui ferait reculer le plateau : trajet a vide a la place
            chains.append(current)
            current = nxt
            continue
        parts = [current] + ([connector] if connector is not None and len(connector) else []) + [nxt]
        current = np.vstack(parts)
    chains.append(current)
    return chains


def _infill_angle_for_layer(layer: Layer, isovalue_rank: dict) -> float:
    """
    Orientation du balayage d'une couche normale (hors pont) en
    CROSSED_ZIGZAG : axe u du repere local (0) une couche sur deux, axe v
    (pi/2) l'autre. La parite est prise sur le RANG de l'isovaleur, pas sur
    la position dans la liste : l'ordre d'impression (step5, en arbre) fait
    se succeder des couches de branches differentes, alors que deux couches
    EMPILEES ont des isovaleurs consecutives - ce sont elles qui doivent se
    croiser. Les morceaux d'une meme isovaleur (clusters) partagent donc la
    meme orientation. Sans isovaleur, repli sur l'identifiant de couche.
    """
    if layer.isovalue is not None and layer.isovalue in isovalue_rank:
        rank = isovalue_rank[layer.isovalue]
    else:
        rank = layer.id
    return 0.0 if rank % 2 == 0 else np.pi / 2.0
