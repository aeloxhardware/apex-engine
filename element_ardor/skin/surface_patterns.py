"""Lissage non planaire : motifs de base (lignes, concentriques, spirales) et outils de chemin.

Découpé de `surface_paths.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `surface_paths`.
"""

import numpy as np
from element_ardor.skin.surface_mesh import (  # noqa: F401
    NOZZLE_WIDTH,
    PatternResult,
    SEED_JITTER,
    ZoneMesh,
)


def clip_keep_above(lines, field, threshold):
    """Coupe les polylignes la ou `field` < threshold (morceaux gardes >= 2
    points). Accepte des Polyline ou des (k, Polyline) ; rend le meme type."""
    out = []
    for item in lines:
        k, pl = item if isinstance(item, tuple) else (None, item)
        keep = pl.value(field) >= threshold - 1e-9
        start = None
        for q in range(len(pl) + 1):
            if q < len(pl) and keep[q]:
                start = q if start is None else start
            elif start is not None:
                if q - start >= 2:
                    piece = pl.subset(slice(start, q))
                    out.append(piece if k is None else (k, piece))
                start = None
    return out


def first_level(width, overlap):
    """Axe du premier cordon : a w/2 - recouvrement du bord (>= 0,01 mm)."""
    return max(width / 2.0 - overlap, 0.01)


def _nested(zm, field, start, width, overlap, order, d_edge, clip=True):
    """Anneaux = lignes de niveau de `field` a start, start + w, ... ; rognes
    au bord ; ordre par niveau croissant (outside_in) ou decroissant."""
    levels = np.arange(start, field.max(), width)
    lines = zm.lines_multi(field, levels)
    if clip:
        lines = clip_keep_above(lines, d_edge, first_level(width, overlap))
    lines.sort(key=lambda kp: kp[0])
    if order != "outside_in":
        lines = lines[::-1]
    return PatternResult(order_by_level([(k, pl.points) for k, pl in lines]),
                         [(k, pl.points) for k, pl in lines], d_edge)


def order_by_level(level_pieces):
    """Ordre de depot des morceaux, niveau par niveau (ordre des niveaux tel
    que recu) ; dans un niveau, plus proche voisin depuis la position
    courante : un morceau ouvert peut etre parcouru a l'envers, un anneau
    ferme commence a son point le plus proche. Reduit les trajets a vide
    quand un anneau est coupe en arcs par le bord."""
    groups = []
    for k, pts in level_pieces:
        if not groups or groups[-1][0] != k:
            groups.append((k, []))
        groups[-1][1].append(pts)
    out, pos = [], None
    for _k, pieces in groups:
        left = list(pieces)
        while left:
            if pos is None:
                out.append(left.pop(0))
                pos = out[-1][-1]
                continue
            best = None
            for q, pts in enumerate(left):
                closed = len(pts) > 2 and np.linalg.norm(pts[0] - pts[-1]) < 1e-6
                if closed:
                    i = int(np.argmin(np.linalg.norm(pts - pos, axis=1)))
                    cand, dist = np.vstack([pts[i:], pts[1:i + 1]]), np.linalg.norm(pts[i] - pos)
                else:
                    da, db = np.linalg.norm(pts[0] - pos), np.linalg.norm(pts[-1] - pos)
                    cand, dist = (pts, da) if da <= db else (pts[::-1], db)
                if best is None or dist < best[0]:
                    best = (dist, q, cand)
            left.pop(best[1])
            out.append(best[2])
            pos = out[-1][-1]
    return out


def concentric_conform(zm, width=NOZZLE_WIDTH, overlap=0.2, order="outside_in", d_edge=None, **_):
    """Anneaux a distance constante du bord (mesuree sur la surface)."""
    d = zm.distance_to_boundary() if d_edge is None else d_edge
    return _nested(zm, d, first_level(width, overlap), width, overlap, order, d, clip=False)


def zone_center(zm, d_edge):
    """Centre d'une zone : sommet le plus loin du bord."""
    return int(np.argmax(d_edge))


def distance_from_center(zm, center, seed_radius=2 * NOZZLE_WIDTH):
    """Distance geodesique au centre = rayon de graine + distance (chaleur
    signee, precise a 0,002 mm) au petit cercle |x - c| = seed_radius (3D ;
    pour un rayon de 0,8 mm, cercle 3D et cercle geodesique coincident a
    mieux que le micron sur les courbures de ces pieces). Remplace le rayon du
    log map, qui se degradait loin du centre sur les zones allongees (cercles
    rapproches a 0,03 mm sur la coque du Benchy) et coutait 14 s."""
    e = np.linalg.norm(zm.V - zm.V[center], axis=1)
    seeds = zm.raw_isolines(e - SEED_JITTER, seed_radius)
    if not seeds:
        return e
    d = zm.distance_to_curves(seeds)
    return np.where(e < seed_radius, e, seed_radius + d)


def concentric_circles(zm, width=NOZZLE_WIDTH, overlap=0.2, order="outside_in", d_edge=None, **_):
    """Cercles geodesiques autour du centre de la zone, rognes au bord. Le
    rayon croit vers l'EXTERIEUR : sens de tri inverse de celui des anneaux."""
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    r = distance_from_center(zm, zone_center(zm, d_edge))
    return _nested(zm, r, width, width, overlap, "inside_out" if order == "outside_in" else "outside_in", d_edge)


def line_direction(zm, angle_deg):
    """Direction des lignes : direction machine (angle depuis X, dans XY)
    projetee dans le plan moyen de la zone ; si elle y est presque normale,
    direction de Z projetee, tournee de l'angle."""
    normal, _c = zm.outward_frame()
    a = np.radians(angle_deg)
    d = np.array([np.cos(a), np.sin(a), 0.0])
    t = d - normal * (d @ normal)
    if np.linalg.norm(t) < 0.2:
        z = np.array([0.0, 0.0, 1.0]) - normal * normal[2]
        z /= max(np.linalg.norm(z), 1e-12)
        t = np.cos(a) * z + np.sin(a) * np.cross(normal, z)
    return t / np.linalg.norm(t), normal


def parallel_lines(zm, width, overlap, angle_deg, aligned, d_edge):
    """Lignes paralleles espacees de `width` MESURE SUR LA SURFACE : niveaux
    de la distance signee a une graine = coupe de la zone par un plan
    contenant la direction des lignes et la normale moyenne ; `aligned` :
    plan par l'origine machine (meme phase d'une zone a l'autre). Rend
    [(k, Polyline)] rognees, orientees selon +t, triees par (k, position)."""
    t, normal = line_direction(zm, angle_deg)
    m = np.cross(t, normal)
    m /= np.linalg.norm(m)
    _n, center = zm.outward_frame()
    # decalage infime de la graine : un plan passant EXACTEMENT par une rangee
    # de sommets (demi-cylindre, z = 0) fragmentait la graine en 605 morceaux
    # degeneres et faisait echouer la chaleur signee
    # aligne : plan a un multiple entier de w de l'origine machine, le plus
    # proche du centre (le plan par l'origine elle-meme ne coupe presque
    # jamais la zone : le repli au centre effacait l'alignement)
    offset = width * np.round((center @ m) / width) if aligned else center @ m
    h = zm.V @ m - offset - SEED_JITTER
    seeds = zm.raw_isolines(h, 0.0)
    if not seeds:                                   # le plan ne coupe pas la zone : graine au centre
        h = (zm.V - center) @ m - SEED_JITTER
        seeds = zm.raw_isolines(h, 0.0)
    u = zm.distance_to_curves(seeds)
    signed = np.where(h >= 0, u, -u)
    k0, k1 = int(np.floor(signed.min() / width)), int(np.ceil(signed.max() / width))
    ks = np.arange(k0, k1 + 1)
    lines = clip_keep_above(zm.lines_multi(signed, ks * width), d_edge, first_level(width, overlap))
    out = []
    for q, pl in sorted(lines, key=lambda kp: (kp[0], float(kp[1].points.mean(axis=0) @ t))):
        out.append((int(ks[q]), pl if (pl.points[-1] - pl.points[0]) @ t >= 0 else pl.subset(slice(None, None, -1))))
    return out


def rectilinear(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=45.0, aligned=False, d_edge=None, **options):
    """Rectiligne : lignes paralleles en zigzag, extremites voisines
    raccordees le long du premier anneau si le raccord fait moins de
    max_connector largeurs (sinon trajet a vide). Options nommees : `zigzag`
    (defaut True), `max_connector` (defaut 2.5) ; les autres sont ignorees."""
    zigzag, max_connector = options.get("zigzag", True), options.get("max_connector", 2.5)
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    lines = parallel_lines(zm, width, overlap, angle_deg, aligned, d_edge)
    level_lines = [(k, pl.points) for k, pl in lines]
    if not zigzag:
        return PatternResult([pl.points for _k, pl in lines], level_lines, d_edge)
    rings = zm.isolines(d_edge, first_level(width, overlap))
    return PatternResult(zigzag_strips([(k, pl.points) for k, pl in lines], rings, max_connector * width),
                         level_lines, d_edge)


def zigzag_strips(pieces, rings, max_len):
    """Zigzags par BANDES : depuis le morceau restant de plus petit niveau,
    on enchaine vers un morceau du niveau suivant joignable par un raccord le
    long du premier anneau (le plus court), en alternant le sens ; sans
    raccord possible, la bande se ferme et la suivante commence. Une zone
    coupee en deux par une encoche ou un trou est ainsi couverte bande par
    bande, au lieu d'alterner gauche/droite a chaque ligne (un trajet a vide
    par ligne). `pieces` : [(k, points orientes +t)] tries par (k, position)."""
    by_k = {}
    for q, (k, _pts) in enumerate(pieces):
        by_k.setdefault(k, []).append(q)
    left = set(range(len(pieces)))
    paths = []
    while left:
        q = min(left)
        left.discard(q)
        k, pts = pieces[q]
        current, flip = [pts], False
        while True:
            end = current[-1][-1]
            best = None
            for r in by_k.get(k + 1, []):
                if r not in left:
                    continue
                cand = pieces[r][1] if flip else pieces[r][1][::-1]
                if np.linalg.norm(cand[0] - end) > max_len:
                    continue
                link = _connector(end, cand[0], rings, max_len)
                if link is None:
                    continue
                length = np.linalg.norm(np.diff(link, axis=0), axis=1).sum() if len(link) > 1 else 0.0
                if best is None or length < best[0]:
                    best = (length, r, link, cand)
            if best is None:
                break
            _l, r, link, cand = best
            left.discard(r)
            current.extend([link, cand])
            k, flip = k + 1, not flip
        paths.append(np.vstack(current))
    return paths


def aligned_rectilinear(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=45.0, d_edge=None, **kw):
    """Rectiligne dont la phase est celle d'un plan passant par l'origine
    machine : les lignes de zones voisines se prolongent."""
    kw.pop("aligned", None)
    return rectilinear(zm, width, overlap, angle_deg, aligned=True, d_edge=d_edge, **kw)


def monotonic_line(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=45.0, d_edge=None, **_):
    """Ligne monotone : lignes paralleles toutes parcourues dans le MEME sens
    (+t), sans raccord (trajet a vide entre deux lignes) : aspect uniforme."""
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    lines = parallel_lines(zm, width, overlap, angle_deg, False, d_edge)
    # ordre ligne apres ligne sur toute la zone (pas par bandes : mesure -3 a -12 % de vide
    # seulement, et l'ordre de recouvrement uniforme est la raison d'etre du motif)
    return PatternResult([pl.points for _k, pl in lines],
                         [(k, pl.points) for k, pl in lines], d_edge)


def monotonic_strips(pieces, max_len):
    """Ordre des lignes monotones par BANDES : apres une ligne du niveau k, la
    ligne du niveau k+1 dont le DEBUT est le plus proche du debut de celle-ci
    (a moins de max_len) ; sinon, nouvelle bande depuis le plus petit niveau
    restant. Meme principe que zigzag_strips, sans raccord (le trajet a vide
    ramene la buse au debut de la ligne suivante)."""
    by_k = {}
    for q, (k, _pts) in enumerate(pieces):
        by_k.setdefault(k, []).append(q)
    left, out = set(range(len(pieces))), []
    while left:
        q = min(left)
        while q is not None:
            left.discard(q)
            k, pts = pieces[q]
            out.append(pts)
            cands = [(np.linalg.norm(pieces[r][1][0] - pts[0]), r) for r in by_k.get(k + 1, []) if r in left]
            cands = [c for c in cands if c[0] <= max_len]
            q = min(cands)[1] if cands else None
    return out


# --- motifs de forme et spirales (etape 4c) -----------------------------------
def signed_line_distance(zm, direction, point):
    """Distance signee MESUREE SUR LA SURFACE a la coupe de la zone par le plan
    passant par `point`, de normale `direction` (positive du cote de
    `direction`). Meme graine que le rectiligne : ses lignes de niveau sont des
    paralleles exactement espacees, meme sur surface courbe."""
    h = (zm.V - point) @ direction - SEED_JITTER
    seeds = zm.raw_isolines(h, 0.0)
    if not seeds:
        return h
    u = zm.distance_to_curves(seeds)
    return np.where(h >= 0, u, -u)


def shape_field(zm, shape, angle_deg, center):
    """Champ dont les lignes de niveau sont des formes imbriquees espacees de
    leur valeur : chaque cote appartient a une famille de paralleles exactes
    (signed_line_distance), donc l'espacement des cotes est exact et les coins
    restent vifs. Carre = max(|s0|, |s90|) ; triangle equilateral = max des
    trois distances signees a 120 deg ; octagramme (etoile a 8 branches) =
    union de deux carres tournes de 45 deg = min des deux champs carres."""
    t, normal = line_direction(zm, angle_deg)
    m = np.cross(t, normal)
    m /= np.linalg.norm(m)
    p = zm.V[center]

    def s(phi_deg):
        phi = np.radians(phi_deg)
        return signed_line_distance(zm, np.cos(phi) * m + np.sin(phi) * t, p)

    if shape == "square":
        return np.maximum(np.abs(s(0)), np.abs(s(90)))
    if shape == "triangle":
        return np.maximum.reduce([s(90), s(210), s(330)])
    if shape == "octagram":
        return np.minimum(np.maximum(np.abs(s(0)), np.abs(s(90))),
                          np.maximum(np.abs(s(45)), np.abs(s(135))))
    raise ValueError(shape)


def concentric_shape(zm, shape, width=NOZZLE_WIDTH, overlap=0.2, order="outside_in", angle_deg=45.0,
                     d_edge=None):
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    f = shape_field(zm, shape, angle_deg, zone_center(zm, d_edge))
    return _nested(zm, f, width, width, overlap, "inside_out" if order == "outside_in" else "outside_in", d_edge)


def concentric_square(zm, width=NOZZLE_WIDTH, overlap=0.2, order="outside_in", angle_deg=45.0, d_edge=None, **_):
    """Carres imbriques (cotes a angle_deg), espaces de w, rognes au bord."""
    return concentric_shape(zm, "square", width, overlap, order, angle_deg, d_edge)


def concentric_triangle(zm, width=NOZZLE_WIDTH, overlap=0.2, order="outside_in", angle_deg=45.0, d_edge=None, **_):
    """Triangles equilateraux imbriques, espaces de w, rognes au bord."""
    return concentric_shape(zm, "triangle", width, overlap, order, angle_deg, d_edge)


def merge_polylines(pieces, tol=1e-6):
    """Recolle les morceaux dont une extremite coincide (a tol pres), en
    retournant au besoin : les tours d'une spirale coupes par la couture."""
    pieces = [np.asarray(p) for p in pieces]
    merged = True
    while merged and len(pieces) > 1:
        merged = False
        ends = np.array([[p[0], p[-1]] for p in pieces])            # (P,2,3)
        for a in range(len(pieces)):
            d = np.linalg.norm(ends[:, :, None, :] - ends[a][None, None, :, :], axis=3)   # (P,2,2)
            d[a] = np.inf
            b, eb, ea = np.unravel_index(np.argmin(d), d.shape)
            if d[b, eb, ea] > tol:
                continue
            pa = pieces[a] if ea == 1 else pieces[a][::-1]          # pa se termine au raccord
            pb = pieces[b] if eb == 0 else pieces[b][::-1]          # pb commence au raccord
            new = np.vstack([pa, pb[1:]])
            pieces = [p for q, p in enumerate(pieces) if q not in (a, b)] + [new]
            merged = True
            break
    return pieces


def _nearest_value(V, field):
    from scipy.spatial import cKDTree
    tree = cKDTree(V)
    return lambda x: float(field[tree.query(x)[1]])


def spiral(zm, radial, center, width, overlap, order, d_edge):
    """Spirale = ligne de niveau de g = radial - w.theta/2pi, theta = angle
    polaire geodesique (log map) autour du centre. A theta fixe, deux tours
    successifs sont a w l'un de l'autre en `radial` (distance exacte) :
    l'espacement vient de `radial`, l'angle ne fait que la forme. g saute de w
    a la coupure theta = +-pi : le maillage y est OUVERT (sommets du cote
    theta < 0 dupliques avec theta + 2pi dans les triangles a cheval), les
    niveaux k*w s'y enchainent (fin du niveau k = debut du niveau k+1, meme
    point 3D) et sont recolles en une seule courbe."""
    from types import SimpleNamespace
    lm = zm.log_map(center)
    theta = np.arctan2(lm[:, 1], lm[:, 0])
    F = zm.F.copy()
    th = theta[F]
    straddle = (th.max(axis=1) - th.min(axis=1) > np.pi) & ~(F == center).any(axis=1)
    low = straddle[:, None] & (th < 0)
    dup = np.unique(F[low])
    remap = np.full(len(zm.V), -1)
    remap[dup] = len(zm.V) + np.arange(len(dup))
    F[low] = remap[F[low]]
    V2 = np.vstack([zm.V, zm.V[dup]])
    theta2 = np.r_[theta, theta[dup] + 2 * np.pi]
    radial2 = np.r_[radial, radial[dup]]
    edge2 = np.r_[d_edge, d_edge[dup]]
    g = radial2 - width * theta2 / (2 * np.pi)
    ks = np.arange(np.floor(g.min() / width), np.ceil(g.max() / width) + 1)
    lines = ZoneMesh.lines_multi(SimpleNamespace(V=V2, F=F), g, ks * width)
    # rognage : au bord (premier niveau) et au centre (rayon < w, angle mal defini)
    keep = np.minimum(edge2 - first_level(width, overlap), radial2 - width)
    lines = clip_keep_above(lines, keep, 0.0)
    level_lines = [(int(k), pl.points) for k, pl in lines]
    # valeur EXACTE (interpolee sur l'arete) du rayon aux extremites : deux arcs
    # successifs de l'etoile ne different que de w/8 = 0,05 mm, moins que
    # l'erreur de la valeur au sommet le plus proche (~0,1 mm), qui melangeait
    # l'ordre (remarque de l'utilisateur, octagramme sur le disque)
    exact = {}
    for _k, pl in lines:
        r = pl.value(radial2)
        exact[tuple(np.round(pl.points[0], 9))] = r[0]
        exact[tuple(np.round(pl.points[-1], 9))] = r[-1]
    fallback = _nearest_value(zm.V, radial)

    def value(x):
        key = tuple(np.round(x, 9))
        return exact[key] if key in exact else fallback(x)

    paths = merge_polylines([pl.points for _k, pl in lines])
    oriented = []
    for p in paths:
        inward = value(p[0]) >= value(p[-1])
        oriented.append(p if inward == (order == "outside_in") else p[::-1])
    sign = -1.0 if order == "outside_in" else 1.0
    # le rayon croit de facon monotone le long de la spirale : trier par la
    # valeur au DEBUT de chaque morceau (deja oriente) suit l'ordre de la
    # spirale elle-meme
    oriented.sort(key=lambda p: sign * value(p[0]))
    return PatternResult(oriented, level_lines, d_edge)


def archimedean(zm, width=NOZZLE_WIDTH, overlap=0.2, order="outside_in", d_edge=None, **_):
    """Spirale d'Archimede geodesique autour du centre : un tour = w."""
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    c = zone_center(zm, d_edge)
    return spiral(zm, distance_from_center(zm, c), c, width, overlap, order, d_edge)


def octagram_spiral(zm, width=NOZZLE_WIDTH, overlap=0.2, order="outside_in", angle_deg=45.0, d_edge=None, **_):
    """Spirale en etoile a 8 branches : meme construction, rayon = champ
    octagramme (distance exacte depuis les cotes de l'etoile)."""
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    c = zone_center(zm, d_edge)
    return spiral(zm, shape_field(zm, "octagram", angle_deg, c), c, width, overlap, order, d_edge)


def path_neighbor_gaps(paths, width, skip_arc=1.5):
    """Ecart entre passes voisines pour un motif sans lignes de niveau
    (Hilbert) : pour chaque point, distance 3D au point le plus proche d'une
    AUTRE portion de trajectoire (a plus de skip_arc * w le long de la courbe,
    ou d'une autre passe)."""
    from scipy.spatial import cKDTree
    pts, arc, pid = [], [], []
    for q, p in enumerate(paths):
        d = densify(p, width / 8)
        pts.append(d)
        arc.append(np.r_[0, np.cumsum(np.linalg.norm(np.diff(d, axis=0), axis=1))])
        pid.append(np.full(len(d), q))
    pts, arc, pid = np.vstack(pts), np.concatenate(arc), np.concatenate(pid)
    tree = cKDTree(pts)
    dist, idx = tree.query(pts, k=64, distance_upper_bound=2 * width)
    out = np.full(len(pts), np.nan)
    for c in range(1, idx.shape[1]):
        j = idx[:, c]
        valid = j < len(pts)
        jj = np.where(valid, j, 0)
        other = valid & ((pid[jj] != pid) | (np.abs(arc[jj] - arc) > skip_arc * width))
        new = np.isnan(out) & other
        out[new] = dist[new, c]
    return out, pts                                   # NaN : aucun voisin a moins de 2 w


def _connector(a, b, rings, max_len):
    """Chemin le long d'un anneau entre ses points les plus proches de a et
    de b (arc le plus court), ou None si trop long."""
    best = None
    for ring in rings:
        ia = int(np.argmin(np.linalg.norm(ring - a, axis=1)))
        ib = int(np.argmin(np.linalg.norm(ring - b, axis=1)))
        if np.linalg.norm(ring[ia] - a) > max_len or np.linalg.norm(ring[ib] - b) > max_len:
            continue
        closed = np.linalg.norm(ring[0] - ring[-1]) < 1e-6
        cum = np.r_[0, np.cumsum(np.linalg.norm(np.diff(ring, axis=0), axis=1))]
        lo, hi = sorted((ia, ib))
        direct = cum[hi] - cum[lo]
        if closed and cum[-1] - direct < direct:
            arc = np.vstack([ring[hi:], ring[1:lo + 1]])
            arc = arc if ia == hi else arc[::-1]
            length = cum[-1] - direct
        else:
            arc = ring[lo:hi + 1] if ia == lo else ring[lo:hi + 1][::-1]
            length = direct
        if length <= max_len and (best is None or length < best[0]):
            best = (length, arc)
    return None if best is None else best[1]


# --- mesures de controle (independantes des champs) ---------------------------
def densify(polyline, step=0.02):
    out = [polyline[:1]]
    for a, b in zip(polyline[:-1], polyline[1:]):
        n = max(1, int(np.ceil(np.linalg.norm(b - a) / step)))
        out.append(a + np.linspace(0, 1, n + 1)[1:, None] * (b - a))
    return np.vstack(out)


def gaps_between(lines_a, lines_b):
    """Distance 3D de chaque sommet des lignes A a la ligne B la plus proche
    (B densifiee a 0,02 mm : erreur de mesure < 0,001 mm)."""
    from scipy.spatial import cKDTree
    if not lines_a or not lines_b:
        return np.empty(0)
    tree = cKDTree(np.vstack([densify(b) for b in lines_b]))
    return tree.query(np.vstack(lines_a))[0]
