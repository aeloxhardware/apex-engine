"""Lissage non planaire : sous-maillage d'une zone (ZoneMesh), bords, polylignes et résultat de motif.

Découpé de `surface_paths.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `surface_paths`.
"""

import numpy as np
import potpourri3d as pp3d
import trimesh


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
