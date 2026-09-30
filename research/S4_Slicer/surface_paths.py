"""
Lissage non-planaire -- P4 : trajectoires des couches de lissage, tracees
DIRECTEMENT sur la surface (decision de l'utilisateur, 2026-09-27).

Principe : 9 motifs sur 10 sont des LIGNES DE NIVEAU d'une distance mesuree
sur la surface depuis une graine (bord, centre, forme englobante, ligne de
coupe), espacees de la largeur de buse. Une distance geodesique a un gradient
de norme 1 sur la surface : l'espacement est exact par construction, meme sur
une forme bombee (seule la forme du motif se deforme, jamais l'espacement).

Etape 4a (ce fichier, pour l'instant) : outils de surface.
  - ZoneMesh : sous-maillage d'une zone de P2-bis, raffine (aretes <= REFINE)
    -- mesure : a 0,38 mm d'arete, espacement 0,401-0,425 mm pour 0,4 vise ;
    a 0,19 mm, 0,401-0,408 (demi-sphere R = 10) ;
  - distances geodesiques par la methode de la CHALEUR SIGNEE (potpourri3d,
    MeshSignedHeatSolver) depuis des courbes (bord, graines) ou des sommets.
    Mesure de la POSITION des lignes (ecart a l'exact median / 95e / max) :
      demi-sphere : signee 0,002/0,016/0,031 ; fast marching 0,007/0,035/0,066 ;
                    chaleur simple 0,045/0,148/0,213 (lignes toutes GLISSEES :
                    espacement juste, premier anneau et recouvrement faux) ;
      disque en triangles-aiguilles : signee 0,059/0,228/0,357 ; chaleur
                    0,152/0,435/0,720 ; fast marching 0,348/1,100/1,278 (le fast
                    marching se degrade sur les triangles obtus) ;
    raccourcir la diffusion de la chaleur simple EMPIRE (t = 0,1 : x6) ;
    courbe OUVERTE traversant la zone de bord a bord (graine du rectiligne) :
    distance SIGNEE puis valeur absolue, 0,001 mm (demi-cylindre) -- en mode
    non signe, la chaleur signee donne n'importe quoi (3,5 mm d'ecart median) ;
    depuis un POINT, la chaleur signee est mediocre (0,13-0,17 mm) : la
    distance au centre est le RAYON des coordonnees polaires (0,000-0,011) ;
  - coordonnees polaires geodesiques (log map, methode de la chaleur
    vectorielle) ;
  - lignes de niveau (marching triangles) rendues en 3D avec leur face et
    leurs coordonnees barycentriques (pour la normale, en P4 plus loin).

Dependance : potpourri3d (MIT ; numpy, scipy), venv d'Apex.
"""
import os
import sys

import numpy as np
import potpourri3d as pp3d
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

NOZZLE_WIDTH = 0.4
SEED_JITTER = 1.234567e-4        # mm : ecarte une graine plane des rangees de sommets
REFINE = NOZZLE_WIDTH / 2.0
# Qualite mediane des triangles au-dessus de laquelle les distances passent par
# le FAST MARCHING : mesure sur les vraies zones (qualite 0,83), anneaux du
# conforme 0,386-0,411 mm contre 0,357-0,490 en chaleur signee, et la chaleur
# signee perd sa precision loin de la source (plateau de distance sur 13 % d'une
# zone du dino a 26 mm) ; en dessous (triangles-aiguilles), le fast marching se
# degrade (1,1 mm) et la chaleur signee reste la plus sure.
GOOD_MESH_QUALITY = 0.5      # arete maximale du maillage de calcul (mesure : +-2 % d'espacement)


def boundary_edges(faces):
    """Aretes n'appartenant qu'a une face (bord), orientees comme dans leur face."""
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    key = np.sort(e, axis=1)
    _u, inv, count = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    return e[count[inv.ravel()] == 1]


def boundary_loops(faces):
    """Boucles de bord (listes ordonnees de sommets)."""
    nxt = {int(a): int(b) for a, b in boundary_edges(faces)}
    loops, seen = [], set()
    for start in list(nxt):
        if start in seen:
            continue
        loop, v = [], start
        while v not in seen and v in nxt:
            seen.add(v)
            loop.append(v)
            v = nxt[v]
        loops.append(np.array(loop))
    return loops


def split_bowtie_vertices(vertices, faces):
    """Duplique chaque sommet ou la zone se touche elle-meme en un seul point
    (plusieurs eventails de faces sans arete commune autour du sommet) : une
    copie par eventail. potpourri3d refuse ces maillages (zone 19 du Benchy,
    100 mm2). La geometrie est inchangee ; seule la topologie est ouverte."""
    faces = faces.copy()
    vertices = list(vertices)
    incident = {}
    for fi, tri in enumerate(faces):
        for v in tri:
            incident.setdefault(int(v), []).append(fi)
    for v, fs in incident.items():
        if len(fs) < 2:
            continue
        parent = {x: x for x in fs}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for a_i, a in enumerate(fs):
            for b in fs[a_i + 1:]:
                if len(set(faces[a]) & set(faces[b])) == 2:
                    parent[find(a)] = find(b)
        fans = {}
        for x in fs:
            fans.setdefault(find(x), []).append(x)
        for extra in list(fans.values())[1:]:
            vertices.append(vertices[v])
            nv = len(vertices) - 1
            for fi in extra:
                faces[fi][faces[fi] == v] = nv
    return np.asarray(vertices), faces


class ZoneMesh:
    """Maillage de calcul d'une zone : faces de la zone, raffinees."""

    def __init__(self, vertices, faces, refine=REFINE):
        used, local = np.unique(faces, return_inverse=True)
        v, f = split_bowtie_vertices(vertices[used], local.reshape(-1, 3))
        v, f, parent = trimesh.remesh.subdivide_to_size(v, f, refine, max_iter=20, return_index=True)
        self.V, self.F = np.asarray(v, float), np.asarray(f, int)
        self.parent = np.asarray(parent)                # face de zone d'origine de chaque face raffinee
        self.E = pp3d.edges(self.V, self.F)
        self.boundary = np.unique(boundary_edges(self.F))
        # qualite des triangles (1 = equilateral) : decide la methode de distance
        tri = self.V[self.F]
        area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
        e2 = (np.linalg.norm(tri - np.roll(tri, 1, axis=1), axis=2) ** 2).sum(axis=1)
        self.quality = float(np.median(4 * np.sqrt(3) * area / np.maximum(e2, 1e-30)))
        self.method = "fast_marching" if self.quality >= GOOD_MESH_QUALITY else "signed_heat"
        self._heat = self._vector = self._march = self._fm = self._signed = None

    # --- solveurs (construits a la demande, reutilises) ----------------------
    @property
    def heat(self):
        if self._heat is None:
            self._heat = pp3d.MeshHeatMethodDistanceSolver(self.V, self.F)
        return self._heat

    @property
    def vector(self):
        if self._vector is None:
            self._vector = pp3d.MeshVectorHeatSolver(self.V, self.F)
        return self._vector

    @property
    def fm(self):
        if self._fm is None:
            self._fm = pp3d.MeshFastMarchingDistanceSolver(self.V, self.F)
        return self._fm

    @property
    def signed(self):
        if self._signed is None:
            self._signed = pp3d.MeshSignedHeatSolver(self.V, self.F)
        return self._signed

    def distance_to_curves(self, curves):
        """Distance a des courbes FERMEES ou traversant la zone de bord a bord,
        donnees en points barycentriques potpourri3d [[(element, [coords]),
        ...], ...] : chaleur signee, puis valeur absolue (voir en tete)."""
        curves = [c for c in (self.clean_curve(c) for c in curves) if len(c) >= 2]
        if self.method == "fast_marching":
            return np.abs(np.asarray(self.fm.compute_distance(list(curves), [[] for _ in curves], False)).ravel())
        try:
            d = self.signed.compute_distance(list(curves), [True] * len(curves), [])
            if not np.all(np.isfinite(d)):
                raise ValueError("distance non finie")
        except (ValueError, RuntimeError):          # solveur en echec : repli fast marching
            self.fallbacks = getattr(self, "fallbacks", 0) + 1
            d = self.fm.compute_distance(list(curves), [[] for _ in curves], False)
        return np.abs(np.asarray(d).ravel())

    def clean_curve(self, curve, snap=1e-3):
        """Graine sans segment quasi nul : un point d'arete a moins de `snap`
        mm d'une extremite est REMPLACE par ce sommet (qui appartient a
        toutes les faces voisines : chaque segment reste dans une face), puis
        les doublons consecutifs sont fusionnes. Un segment de 1,4 micron
        (graine circulaire sur la coque du Benchy) donnait des NaN ; retirer
        simplement le point faisait enjamber deux faces (refuse par
        potpourri3d)."""
        out = []
        for el, bary in curve:
            if len(bary) == 1:
                a, b = self.E[el]
                length = np.linalg.norm(self.V[b] - self.V[a])
                if bary[0] * length < snap:
                    el, bary = int(a), []
                elif (1 - bary[0]) * length < snap:
                    el, bary = int(b), []
            item = (int(el), list(bary))
            if out and not item[1] and not out[-1][1] and out[-1][0] == item[0]:
                continue                                   # meme sommet deux fois de suite
            out.append(item)
        return out

    def distance_to_curve(self, curve_points):
        return self.distance_to_curves([curve_points])

    def distance_to_center(self, center_vertex):
        """Distance au centre = rayon des coordonnees polaires geodesiques."""
        return np.linalg.norm(self.log_map(center_vertex), axis=1)

    def distance_to_boundary(self):
        return self.distance_to_curves([[(int(i), []) for i in np.r_[lp, lp[:1]]]
                                        for lp in boundary_loops(self.F)])

    def log_map(self, center_vertex):
        """Coordonnees polaires geodesiques (V,2) autour d'un sommet."""
        return np.asarray(self.vector.compute_log_map(int(center_vertex)))

    # --- lignes de niveau ----------------------------------------------------
    def raw_isolines(self, field, level):
        """Lignes de niveau au format barycentrique potpourri3d (utilisable
        directement comme graine d'une distance)."""
        if self._march is None:
            self._march = pp3d.MeshMarchingTrianglesSolver(self.V, self.F)
        return self._march.marching_triangles(np.asarray(field, float), float(level))

    def lines(self, field, level):
        """Composantes de la ligne de niveau : liste de Polyline (points 3D +
        gabarit (i, j, t) pour y evaluer tout champ aux sommets)."""
        return [self.polyline(curve) for curve in self.raw_isolines(field, level)]

    def isolines(self, field, level):
        """Composantes de la ligne de niveau : liste de tableaux de points (n,3)."""
        return [pl.points for pl in self.lines(field, level)]

    def lines_multi(self, field, levels):
        """Lignes de niveau pour TOUS les niveaux en un passage vectorise (le
        marching de potpourri3d parcourt tout le maillage a chaque niveau :
        40 ms x 708 niveaux = 28 s mesures sur une zone du dino). Rend une
        liste de (indice de niveau, Polyline)."""
        field = np.asarray(field, float)
        levels = np.asarray(levels, float)
        if not len(levels):
            return []
        # valeurs exactement sur un niveau : decalees d'un epsilon (evite les
        # segments degeneres passant par un sommet)
        f = field.copy()
        step = np.min(np.diff(levels)) if len(levels) > 1 else 1.0
        eps = 1e-9 * max(step, 1e-6)
        k_near = np.clip(np.searchsorted(levels, f), 0, len(levels) - 1)
        on = np.abs(f - levels[k_near]) < eps
        f[on] += 2 * eps
        fe = f[self.F]                                        # (T,3)
        lo, hi = fe.min(axis=1), fe.max(axis=1)
        k0 = np.searchsorted(levels, lo, side="left")        # premier niveau > min
        k1 = np.searchsorted(levels, hi, side="left")        # niveaux < max : [k0, k1)
        count = np.maximum(k1 - k0, 0)
        if not count.any():
            return []
        tri = np.repeat(np.arange(len(self.F)), count)
        kk = np.concatenate([np.arange(a, b) for a, b in zip(k0[count > 0], k1[count > 0])])
        lv = levels[kk]
        v = self.F[tri]                                        # (S,3)
        s = f[v] - lv[:, None]
        pos = s > 0
        alone = np.where(pos.sum(axis=1) == 1, np.argmax(pos, axis=1), np.argmin(pos, axis=1))
        a = np.take_along_axis(v, alone[:, None], 1)[:, 0]
        b = np.take_along_axis(v, ((alone + 1) % 3)[:, None], 1)[:, 0]
        c = np.take_along_axis(v, ((alone + 2) % 3)[:, None], 1)[:, 0]
        # deux points du segment : sur les aretes (a,b) et (a,c), identifies par (arete, niveau)
        def pt_key(p, q):
            lo_, hi_ = np.minimum(p, q), np.maximum(p, q)
            return lo_, hi_

        (i1, j1), (i2, j2) = pt_key(a, b), pt_key(a, c)
        nv = len(self.V)
        key1 = (i1.astype(np.int64) * nv + j1) * len(levels) + kk
        key2 = (i2.astype(np.int64) * nv + j2) * len(levels) + kk
        keys, inv = np.unique(np.concatenate([key1, key2]), return_inverse=True)
        inv = inv.ravel()
        p1, p2 = inv[:len(tri)], inv[len(tri):]
        # geometrie des points
        kidx = keys % len(levels)
        ej = (keys // len(levels)) % nv
        ei = (keys // len(levels)) // nv
        t = (levels[kidx] - f[ei]) / np.where(f[ej] != f[ei], f[ej] - f[ei], 1.0)
        t = np.clip(t, 0.0, 1.0)
        # chainage : chaque point a au plus deux segments
        adj = [[] for _ in range(len(keys))]
        for sidx, (u, w) in enumerate(zip(p1.tolist(), p2.tolist())):
            adj[u].append(w)
            adj[w].append(u)
        seen = np.zeros(len(keys), bool)
        out = []
        order = sorted(range(len(keys)), key=lambda q: len(adj[q]))   # extremites (degre 1) d'abord
        for start in order:
            if seen[start]:
                continue
            chain, prev, cur = [start], -1, start
            seen[start] = True
            while True:
                nxt = [q for q in adj[cur] if q != prev and not seen[q]]
                if not nxt:
                    if len(adj[cur]) == 2 and start in adj[cur] and len(chain) > 2:
                        chain.append(start)                  # boucle fermee
                    break
                prev, cur = cur, nxt[0]
                seen[cur] = True
                chain.append(cur)
            if len(chain) >= 2:
                ch = np.array(chain)
                i, j, tt = ei[ch], ej[ch], t[ch]
                out.append((int(kidx[ch[0]]), Polyline(self.V[i] + tt[:, None] * (self.V[j] - self.V[i]), i, j, tt)))
        return out

    def polyline(self, curve):
        n = len(curve)
        i, j, t = np.empty(n, int), np.empty(n, int), np.zeros(n)
        for k, (el, bary) in enumerate(curve):
            if len(bary) == 1:                          # point sur une arete
                i[k], j[k] = self.E[el]
                t[k] = bary[0]
            elif len(bary) == 0:                        # sommet
                i[k] = j[k] = el
            else:                                       # dans une face : ramene sur l'arete la plus proche
                f = self.F[el]
                b = list(bary) + ([1.0 - sum(bary)] if len(bary) == 2 else [])
                a0, a1 = np.argsort(b)[::-1][:2]
                i[k], j[k] = f[a0], f[a1]
                t[k] = b[a1] / max(b[a0] + b[a1], 1e-12)
        return Polyline(self.V[i] + t[:, None] * (self.V[j] - self.V[i]), i, j, t)

    def outward_frame(self):
        """Normale moyenne (ponderee par l'aire) et centre de la zone."""
        tri = self.V[self.F]
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        normal = n.sum(axis=0)
        return normal / max(np.linalg.norm(normal), 1e-12), tri.mean(axis=1).mean(axis=0)


class Polyline:
    """Points 3D d'une ligne tracee sur la surface, avec le gabarit
    d'interpolation (point = V[i] + t (V[j] - V[i]))."""

    def __init__(self, points, i, j, t):
        self.points, self.i, self.j, self.t = points, i, j, t

    def value(self, field):
        return field[self.i] + self.t * (field[self.j] - field[self.i])

    def subset(self, sl):
        return Polyline(self.points[sl], self.i[sl], self.j[sl], self.t[sl])

    def __len__(self):
        return len(self.points)


# --- motifs (etape 4b) --------------------------------------------------------
class PatternResult:
    """Sortie commune des motifs : `paths` = trajectoires (N,3) dans l'ordre de
    depot, un trajet a vide entre deux paths ; `level_lines` = [(k, points)]
    lignes brutes numerotees par niveau (controle d'espacement entre niveaux
    VOISINS) ; `d_edge` = distance au bord par sommet."""

    def __init__(self, paths, level_lines, d_edge):
        self.paths, self.level_lines, self.d_edge = paths, level_lines, d_edge


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


def rectilinear(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=45.0, aligned=False, d_edge=None,
                zigzag=True, max_connector=2.5, **_):
    """Rectiligne : lignes paralleles en zigzag, extremites voisines
    raccordees le long du premier anneau si le raccord fait moins de
    max_connector largeurs (sinon trajet a vide)."""
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


# --- Hilbert : depliage ARAP (etape 4d) ----------------------------------------
def arap_unfold(zm, iterations=30):
    """Depliage au plus rigide (ARAP, libigl, MPL-2.0) initialise par LSCM
    (deux sommets de bord les plus eloignes epingles). ARAP vise des
    triangles plans ISOMETRIQUES aux triangles 3D : l'echelle est conservee,
    l'ecart restant est mesure par les valeurs singulieres (arap_stretch)."""
    import igl
    V, F = np.ascontiguousarray(zm.V, float), np.ascontiguousarray(zm.F, np.int64)
    bnd = np.unique(boundary_edges(zm.F))
    i0 = int(bnd[np.argmax(np.linalg.norm(zm.V[bnd] - zm.V[bnd].mean(0), axis=1))])
    i1 = int(bnd[np.argmax(np.linalg.norm(zm.V[bnd] - zm.V[i0], axis=1))])
    b = np.array([i0, i1], np.int64)
    bc = np.array([[0.0, 0.0], [np.linalg.norm(zm.V[i1] - zm.V[i0]), 0.0]])
    uv0, _q = igl.lscm(V, F, b, bc)
    data = igl.ARAPData()
    data.max_iter = iterations
    # un sommet epingle (sa position LSCM) : sans aucune contrainte, la
    # precalculation echoue sur le demi-cylindre (systeme singulier) ; deux
    # epingles imposeraient la corde 3D au lieu de la distance sur la surface
    igl.arap_precomputation(V, F, 2, np.array([i0], np.int32), data)
    return np.asarray(igl.arap_solve(np.ascontiguousarray(uv0[[i0]]), data, np.ascontiguousarray(uv0)))


def arap_stretch(zm, uv):
    """Valeurs singulieres (T,2) du passage triangle 3D -> plan (1 = isometrie)
    et aire 3D par triangle."""
    tri = zm.V[zm.F]
    e1, e2 = tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]
    n = np.cross(e1, e2)
    area = 0.5 * np.linalg.norm(n, axis=1)
    x = e1 / np.linalg.norm(e1, axis=1)[:, None]
    y = np.cross(n / (2 * area)[:, None], x)
    dq = np.stack([np.stack([(e1 * x).sum(1), (e1 * y).sum(1)], 1),
                   np.stack([(e2 * x).sum(1), (e2 * y).sum(1)], 1)], axis=2)       # (T,2,2)
    t = uv[zm.F]
    du = np.stack([t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]], axis=2)
    return np.linalg.svd(du @ np.linalg.inv(dq), compute_uv=False), area


def locate_uv(zm, uv, pts2d, ks=(16, 64, 256, 1024)):
    """Triangle et barycentriques de chaque point 2D dans le depliage ; -1 si
    hors. Candidats = triangles aux centres les plus proches, recherche
    elargie (ks) pour les points non trouves : sur un maillage d'aiguilles
    (demi-cylindre), 12 candidats laissaient 8 % de points non localises en
    PLEIN milieu de zone, donc des trous dans la courbe."""
    from scipy.spatial import cKDTree
    tri = uv[zm.F]
    cen = tri.mean(axis=1)
    tree = cKDTree(cen)
    face = np.full(len(pts2d), -1)
    bary = np.zeros((len(pts2d), 3))
    # hors du depliage : plus loin de tout centre que le plus grand rayon de triangle
    rmax = np.linalg.norm(tri - cen[:, None, :], axis=2).max()
    outside = tree.query(pts2d)[0] > rmax + 1e-9
    done = 0
    for k in ks:
        todo_idx = np.flatnonzero((face < 0) & ~outside)
        if not len(todo_idx):
            break
        k = min(k, len(zm.F))
        _d, cand = tree.query(pts2d[todo_idx], k=k)
        cand = cand.reshape(len(todo_idx), -1)
        for c in range(done, cand.shape[1]):
            sub = face[todo_idx] < 0
            if not sub.any():
                break
            ii = todo_idx[sub]
            f = cand[sub, c]
            a, b_, c_ = tri[f, 0], tri[f, 1], tri[f, 2]
            v0, v1, v2 = b_ - a, c_ - a, pts2d[ii] - a
            den = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
            den = np.where(np.abs(den) < 1e-30, 1e-30, den)
            l1 = (v2[:, 0] * v1[:, 1] - v1[:, 0] * v2[:, 1]) / den
            l2 = (v0[:, 0] * v2[:, 1] - v2[:, 0] * v0[:, 1]) / den
            inside = (l1 >= -1e-9) & (l2 >= -1e-9) & (l1 + l2 <= 1 + 1e-9)
            face[ii[inside]] = f[inside]
            bary[ii[inside]] = np.column_stack([1 - l1[inside] - l2[inside], l1[inside], l2[inside]])
        done = cand.shape[1]
        if k == len(zm.F):
            break
    return face, bary


def hilbert(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=0.0, d_edge=None, uv=None, **_):
    """Courbe de Hilbert au pas w tracee dans le depliage ARAP de la zone,
    ramenee sur la surface par barycentriques (points densifies a w/4 : la
    courbe suit la surface entre deux noeuds de grille), rognee au premier
    niveau. Le pas reel sur la surface = w x etirement local du depliage (mesure,
    pas corrige : decision de l'utilisateur apres mesure)."""
    import math
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    uv = arap_unfold(zm) if uv is None else uv
    a = np.radians(angle_deg)
    rot = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    q = uv @ rot                                            # grille tournee de angle_deg
    lo, hi = q.min(axis=0), q.max(axis=0)
    extent = float(max(hi - lo))
    order = max(1, math.ceil(math.log2(max(extent / width, 1.0))))
    n = 2 ** order
    xs, ys = _hilbert_d2xy(n, np.arange(n * n))
    grid = lo + (np.column_stack([xs, ys]) + 0.5) * width   # pas exact w (grille plus grande que la zone)
    # densification a w/4 le long de la courbe
    sub = 4
    tt = np.arange(sub) / sub
    dense = (grid[:-1, None, :] + tt[None, :, None] * (grid[1:] - grid[:-1])[:, None, :]).reshape(-1, 2)
    dense = np.vstack([dense, grid[-1:]])
    face, bary = locate_uv(zm, q, dense)
    ok = face >= 0
    pts = np.zeros((len(dense), 3))
    de = np.full(len(dense), -1.0)
    fv = zm.F[face[ok]]
    pts[ok] = (bary[ok, :, None] * zm.V[fv]).sum(axis=1)
    de[ok] = (bary[ok] * d_edge[fv]).sum(axis=1)
    keep = de >= first_level(width, overlap)
    paths, start = [], None
    for i in range(len(dense) + 1):
        if i < len(dense) and keep[i]:
            start = i if start is None else start
        elif start is not None:
            if i - start >= 2:
                paths.append(pts[start:i])
            start = None
    return PatternResult(paths, [], d_edge)


def _hilbert_d2xy(n, d):
    """Indice le long de la courbe -> (x, y) sur une grille n x n (vectorise)."""
    x = np.zeros_like(d)
    y = np.zeros_like(d)
    t = d.copy()
    s = 1
    while s < n:
        rx = 1 & (t // 2)
        ry = 1 & (t ^ rx)
        swap = ry == 0
        flip = swap & (rx == 1)
        x = np.where(flip, s - 1 - x, x)
        y = np.where(flip, s - 1 - y, y)
        x, y = np.where(swap, y, x), np.where(swap, x, y)
        x = x + s * rx
        y = y + s * ry
        t //= 4
        s *= 2
    return x, y


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


# --- cordon de bord (decision utilisateur 2026-09-28) --------------------------
EDGE_OVERLAP = 0.25        # recouvrement motif / cordon de bord, en largeurs (usage des slicers plans ;
                           # mesure : 0,15 -> -0,5 a -2,6 pts de couverture ; a fixer apres essai)


def with_edge_ring(zm, pattern, width=NOZZLE_WIDTH, overlap=0.2, edge_rings=1, edge_overlap=EDGE_OVERLAP,
                   d_edge=None, **kw):
    """Motif precede de `edge_rings` cordons qui suivent le bord de la zone :
    lignes de niveau de la distance au bord MESUREE SUR LA SURFACE (le cordon
    epouse la courbe de la piece), le premier a w/2 - recouvrement du bord,
    les suivants tous les w. Le motif commence un cordon plus loin, moins
    `edge_overlap` largeur (le conforme continue exactement au pas w, ses
    anneaux etant paralleles au cordon). Le motif recoit un recouvrement
    EFFECTIF tel que son premier niveau tombe la : ses raccords (zigzag)
    longent sa propre limite, jamais le cordon deja depose."""
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    fn = pattern if callable(pattern) else PATTERNS[pattern]
    if edge_rings <= 0:
        return fn(zm, width, overlap, d_edge=d_edge, **kw)
    base = first_level(width, overlap)
    rings = [(q, pl.points) for q, pl in zm.lines_multi(d_edge, base + width * np.arange(edge_rings))]
    rings.sort(key=lambda kp: kp[0])
    inner = base + width * (edge_rings - (0.0 if pattern == "concentric_conform" else edge_overlap))
    res = fn(zm, width, width / 2.0 - inner, d_edge=d_edge, **kw)
    ring_paths = order_by_level(rings)
    # dernier cordon ferme : il commence et finit au point le plus proche du
    # debut du motif (sinon un trajet a vide traversait la zone)
    if ring_paths and res.paths:
        last = ring_paths[-1]
        if len(last) > 2 and np.linalg.norm(last[0] - last[-1]) < 1e-6:
            i = int(np.argmin(np.linalg.norm(last - res.paths[0][0], axis=1)))
            ring_paths[-1] = np.vstack([last[i:], last[1:i + 1]])
    return PatternResult(ring_paths + list(res.paths), res.level_lines, d_edge)


# --- Hilbert : avertissement et remedes (decision utilisateur 2026-09-28) ------
HILBERT_STRETCH_TOLERANCE = 0.05     # etirement local admis (5 % = tolerance d'espacement)
HILBERT_WARN_FRACTION = 0.10         # avertissement si plus de 10 % de la surface depasse

HILBERT_OPTIONS = {
    "A": ("Garder tel quel",
          "Le motif reste parfaitement regulier, mais sur les parties bombees les lignes "
          "s'ecartent ou se serrent un peu : de petits jours ou de legeres surepaisseurs "
          "peuvent apparaitre a ces endroits."),
    "B": ("Decouper la zone",
          "La zone est partagee en morceaux assez plats pour que les lignes restent bien "
          "espacees. Le motif reste regulier dans chaque morceau, mais on verra des raccords "
          "entre les morceaux."),
    "C": ("Resserrer la grille",
          "L'ecart entre les lignes redevient juste en moyenne, mais les carres du motif se "
          "deforment en losanges. Sur une forme tres bombee, les pointes des losanges se "
          "touchent presque : les lignes s'y serrent et de petits manques peuvent apparaitre."),
}
# mesure (test du 2026-09-28, part des lignes a w +- 5 %, reference plane 94 %) :
#   demi-sphere A 27 %  B 90 %  C 74 % ; coque Benchy A 88 %  B 93 %  C 88 % ;
#   dino A 52 %  B 92 %  C 74 %. B coute 7 a 18 s par grande zone (ARAP par morceau).


def hilbert_stretch(zm, uv=None):
    """Depliage ARAP et part de surface etiree au-dela de la tolerance."""
    uv = arap_unfold(zm) if uv is None else uv
    sv, area = arap_stretch(zm, uv)
    over = np.abs(sv - 1).max(axis=1) > HILBERT_STRETCH_TOLERANCE
    return uv, float(area[over].sum() / area.sum())


def geodesic_grid_coords(zm, angle_deg, center):
    """Remede C : coordonnees (U, V) = distances signees MESUREES SUR LA
    SURFACE a deux lignes perpendiculaires passant par le centre (meme graine
    que le rectiligne). Chaque famille de lignes de grille est alors exactement
    espacee de w ; les deux familles ne sont plus perpendiculaires sur une
    surface bombee : les carres du motif se deforment."""
    t, normal = line_direction(zm, angle_deg)
    m = np.cross(t, normal)
    m /= np.linalg.norm(m)
    p = zm.V[center]
    return np.column_stack([signed_line_distance(zm, t, p), signed_line_distance(zm, m, p)])


def hilbert_geodesic(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=0.0, d_edge=None, **_):
    """Hilbert, remede C. Les triangles retournes dans (U, V) (les deux champs
    de distance deviennent paralleles loin du centre) sont ecartes de la
    localisation : aucun point n'y est pose deux fois."""
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    uv = geodesic_grid_coords(zm, angle_deg, zone_center(zm, d_edge))
    t = uv[zm.F]
    a, b = t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]
    det = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    good = det * np.sign(np.median(det)) > 0
    sub = ZoneMesh.__new__(ZoneMesh)
    sub.V, sub.F = zm.V, zm.F[good]
    return hilbert(sub, width, overlap, 0.0, d_edge=d_edge, uv=uv)


def split_zone(zm, min_area=25 * NOZZLE_WIDTH ** 2):
    """Remede B : coupe recursive en deux (plan par le barycentre, normal a
    l'axe principal de la zone) jusqu'a ce que chaque morceau se deplie avec
    moins de HILBERT_WARN_FRACTION de surface etiree, ou devienne trop petit.
    Rend [(ZoneMesh du morceau, uv ARAP, part etiree)]. Coupe au niveau des
    triangles du maillage de calcul (0,2 mm) : limite en escalier fin, lisse
    ensuite par le cordon de bord."""
    uv, frac = hilbert_stretch(zm)
    tri = zm.V[zm.F]
    area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    if frac <= HILBERT_WARN_FRACTION or area.sum() < 2 * min_area:
        return [(zm, uv, frac)]
    c = tri.mean(axis=1)
    g = (c * area[:, None]).sum(0) / area.sum()
    _u, _s, vt = np.linalg.svd((c - g) * np.sqrt(area)[:, None], full_matrices=False)
    side = (c - g) @ vt[0] >= 0
    out = []
    for keep in (side, ~side):
        piece = ZoneMesh(zm.V, zm.F[keep], refine=np.inf)
        piece.parent = zm.parent[keep][piece.parent]
        out.extend(split_zone(piece, min_area))
    return out


def hilbert_split(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=0.0, d_edge=None, piece_rings=0,
                  piece_overlap=None, **_):
    """Hilbert, remede B : un Hilbert par morceau, chacun avec sa propre
    distance au bord et, s'il est demande, son propre cordon de bord : le
    long d'une coupe, les deux cordons voisins sont a w l'un de l'autre."""
    piece_overlap = EDGE_OVERLAP if piece_overlap is None else piece_overlap
    paths = []
    for piece, uv, _frac in split_zone(zm):
        paths.extend(with_edge_ring(piece, hilbert, width, overlap, piece_rings, piece_overlap,
                                    angle_deg=angle_deg, uv=uv).paths)
    return PatternResult(paths, [], zm.distance_to_boundary() if d_edge is None else d_edge)


# --- micro-zones (decision utilisateur 2026-09-28) -----------------------------
def merge_micro_zones(vertices, faces, areas, zone, zone_moment, last_access, micro_inradius):
    """Rattache chaque micro-zone (plus grande distance au bord, mesuree sur
    la surface, < micro_inradius : trop etroite pour un cordon de bord et un
    motif) a une zone VOISINE (arete commune) ou la fusion reste imprimable :
    la zone fusionnee s'imprime au plus tard des deux moments, qui doit rester
    <= au dernier moment accessible de TOUTES ses faces. Parmi les voisines
    admissibles : la plus longue frontiere commune. Repete jusqu'a stabilite.
    Rend (zone, moments {z: m}, journal [(z, cible, m_z, m_cible, m, frontiere)],
    restantes {z: "isolee" | "voisines inaccessibles"}) ; les restantes ne
    sont PAS lissees (decision utilisateur) : c'est a l'appelant de les retirer."""
    zone = zone.copy()
    moment = {z: int(m) for z, m in enumerate(zone_moment)}
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    fid = np.tile(np.arange(len(faces)), 3)
    key = np.sort(e, axis=1)
    o = np.lexsort((key[:, 1], key[:, 0]))
    key, fid = key[o], fid[o]
    same = np.all(key[1:] == key[:-1], axis=1)
    pa, pb = fid[:-1][same], fid[1:][same]                     # paires de faces voisines
    elen = np.linalg.norm(vertices[key[:-1][same, 0]] - vertices[key[:-1][same, 1]], axis=1)
    radius = {}

    def inradius(z):
        if z not in radius:
            try:
                radius[z] = float(ZoneMesh(vertices, faces[zone == z]).distance_to_boundary().max())
            except Exception:                                  # zone degeneree : traitee comme micro
                radius[z] = 0.0
        return radius[z]

    log = []
    changed = True
    while changed:
        changed = False
        live = [z for z in moment if (zone == z).any()]
        for z in sorted(live, key=lambda q: areas[zone == q].sum()):
            if not (zone == z).any() or inradius(z) >= micro_inradius:
                continue
            za, zb = zone[pa], zone[pb]
            touch = (za == z) ^ (zb == z)
            border = {}
            for q, w in zip(np.where(za == z, zb, za)[touch], elen[touch]):
                if q >= 0:
                    border[int(q)] = border.get(int(q), 0.0) + float(w)
            best = None
            for q, length in border.items():
                m = max(moment[z], moment[q])
                if m <= last_access[(zone == z) | (zone == q)].min() and (best is None or length > best[1]):
                    best = (q, length, m)
            if best is None:
                continue
            q, length, m = best
            zone[zone == z] = q
            log.append((z, q, moment[z], moment[q], m, length))
            moment[q] = m
            radius.pop(q, None)                                # la cible a grandi : rayon a recalculer
            changed = True
    remaining = {}
    for z in [z for z in moment if (zone == z).any()]:
        if inradius(z) < micro_inradius:
            za, zb = zone[pa], zone[pb]
            touch = (za == z) ^ (zb == z)
            has_neighbour = (np.where(za == z, zb, za)[touch] >= 0).any()
            remaining[z] = "voisines inaccessibles" if has_neighbour else "isolee"
    return zone, moment, log, remaining


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


# --- registre des motifs (memes noms que surface_patterns.PATTERNS) ----------
PATTERNS = {
    "concentric_conform": concentric_conform,
    "concentric": concentric_circles,
    "rectilinear": rectilinear,
    "aligned_rectilinear": aligned_rectilinear,
    "monotonic_line": monotonic_line,
    "concentric_square": concentric_square,
    "concentric_triangle": concentric_triangle,
    "archimedean": archimedean,
    "octagram_spiral": octagram_spiral,
    "hilbert": hilbert,
}
