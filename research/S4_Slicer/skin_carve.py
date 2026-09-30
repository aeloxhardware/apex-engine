"""
Lissage non-planaire -- P3 : creuser le corps de la piece sous les zones qui
recoivent une couche de surface, et poser la premiere couche avec le motif
choisi.

Champ psi aux noeuds du maillage = distance REELLE (piece d'origine, repere
Bird) a la plus proche face retenue par surface_skin.py, moins l'epaisseur de
lissage. psi < 0 : la matiere appartient a la couche de surface, pas au corps.

Une couche est un plan z = h dans l'espace deforme. Sa section par un
tetraedre est un polygone convexe (triangle ou quadrilatere) sur lequel psi
est LINEAIRE : la decoupe psi >= 0 y est donc exacte. La couche creusee est
l'union de ces morceaux. Seules les couches TOUCHEES (un tetraedre coupe a
un noeud a psi < 0) sont recalculees ; les autres gardent exactement leur
contour de l'etape 5.

Venv d'Apex.
"""
import os
import sys

import numpy as np
import shapely
import trimesh
from shapely.geometry import MultiPolygon, Polygon

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from pipeline.step6_toolpath2d import _offset_walls, _polygon_pieces  # noqa: E402

_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

# Garde-fous memoire (un calcul non borne a sature 47 Go et fait planter le PC) :
DISTANCE_BATCH = 50_000       # points par appel a trimesh.proximity.closest_point
MAX_REFINE_TRIANGLES = 400_000  # triangles en attente d'une meme profondeur, par couche
MAX_PROCESS_BYTES = 12 * 2**30  # memoire du processus au-dela de laquelle on s'arrete
# threads de la distance native par appel : 1 dans les processus du creusement
# parallele (deja un processus par coeur), 0 (= tous les coeurs) sinon
_NATIVE_THREADS = int(os.environ.get("APEX_P3_NATIVE_THREADS", "0"))


def process_bytes():
    """Memoire de travail du processus (Windows, sans psutil ; 0 ailleurs)."""
    if os.name != "nt":
        return 0
    import ctypes
    from ctypes import wintypes

    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
    c = Counters()
    c.cb = ctypes.sizeof(c)
    k32 = ctypes.windll.kernel32
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
    k32.K32GetProcessMemoryInfo.restype = wintypes.BOOL
    if not k32.K32GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(c), c.cb):
        raise OSError("GetProcessMemoryInfo a echoue : garde-fou memoire inoperant")
    return c.WorkingSetSize


def check_memory(where):
    used = process_bytes()
    if used > MAX_PROCESS_BYTES:
        raise MemoryError(f"{where} : {used / 2**30:.1f} Go utilises (> {MAX_PROCESS_BYTES / 2**30:.1f}) -- arret")


def skin_psi(original_nodes, skin, thickness):
    """psi par noeud. `skin` = contenu de loop_*_skin_pass<k>.npz."""
    keep = (skin["status"] == 0) & (skin["zone"] >= 0)
    psi = np.full(len(original_nodes), np.inf)
    if not keep.any():
        return psi
    retained = trimesh.Trimesh(skin["vertices"], skin["faces"][keep], process=False)
    _closest, distance, _tri = trimesh.proximity.closest_point(retained, original_nodes)
    return distance - thickness


def _sections(cells, deformed, psi, h):
    """Sections des tetraedres coupes par z = h : liste de (xy (m,2), psi (m,), tetraedre)."""
    z = deformed[cells][:, :, 2]
    cut = np.flatnonzero((z.min(axis=1) < h) & (z.max(axis=1) > h))
    out = []
    if not len(cut):
        return out
    nodes = deformed[cells[cut]]
    values = psi[cells[cut]]
    zc = nodes[:, :, 2]
    points = np.full((len(cut), 6, 2), np.nan)
    pvals = np.full((len(cut), 6), np.nan)
    ids = cells[cut]
    rows = np.arange(len(cut))
    for e, (a0, b0) in enumerate(_EDGES):
        # sens CANONIQUE (du noeud de plus petit numero vers l'autre) : deux
        # tetraedres voisins calculent alors le point de coupe de leur arete
        # commune AU BIT PRES identique. Sinon, ~1e-16 d'ecart laissait des
        # micro-fentes entre sections, donc de faux trous interieurs dans
        # l'union (contours mesures jusqu'a 12 mm a l'interieur de dino).
        swap = ids[:, a0] > ids[:, b0]
        a = np.where(swap, b0, a0)
        b = np.where(swap, a0, b0)
        za, zb = zc[rows, a], zc[rows, b]
        crossing = (za - h) * (zb - h) < 0
        t = np.where(crossing, (h - za) / np.where(crossing, zb - za, 1.0), np.nan)
        pa, pb = nodes[rows, a, :2], nodes[rows, b, :2]
        points[:, e] = pa + t[:, None] * (pb - pa)
        va, vb = values[rows, a], values[rows, b]
        pvals[:, e] = va + t * (vb - va)
    for i in range(len(cut)):
        ok = ~np.isnan(pvals[i])
        xy, pv = points[i, ok], pvals[i, ok]
        if len(xy) < 3:
            continue
        center = xy.mean(axis=0)
        order = np.argsort(np.arctan2(xy[:, 1] - center[1], xy[:, 0] - center[0]))
        out.append((xy[order], pv[order], int(cut[i])))
    return out


def _clip(xy, pv):
    """Partie psi >= 0 d'un polygone convexe (psi lineaire) -- Sutherland-Hodgman
    sur un demi-plan."""
    kept = []
    n = len(xy)
    for i in range(n):
        a, b = xy[i], xy[(i + 1) % n]
        va, vb = pv[i], pv[(i + 1) % n]
        if va >= 0:
            kept.append(a)
        if (va >= 0) != (vb >= 0):
            kept.append(a + (va / (va - vb)) * (b - a))
    return np.array(kept) if len(kept) >= 3 else None


def layer_is_touched(cells, deformed, psi, h):
    z = deformed[cells][:, :, 2]
    cut = (z.min(axis=1) < h) & (z.max(axis=1) > h)
    return bool((psi[cells[cut]] < 0).any())


class ExactDistance:
    """Distance REELLE exacte a la peau en des points de l'espace deforme :
    chaque point est ramene dans la piece d'origine par la transformation
    (affine, donc exacte) de son tetraedre, puis sa distance aux faces
    retenues est calculee (trimesh)."""

    def __init__(self, skin, cells, deformed, original, thickness):
        keep = (skin["status"] == 0) & (skin["zone"] >= 0)
        self.mesh = trimesh.Trimesh(skin["vertices"], skin["faces"][keep], process=False) if keep.any() else None
        self.cells, self.deformed, self.original, self.thickness = cells, deformed, original, thickness
        # Distance native (T9 du 2026-09-30, apex_core.MeshDistance, BVH + threads)
        # par defaut (APEX_P3_NATIVE_DISTANCE=0 pour la voie trimesh) : exacte (trimesh se trompe jusqu'a 2e-5 mm
        # pres des sommets), donc resultats a ~1e-5 mm de la voie trimesh.
        self.native = None
        if self.mesh is not None and os.environ.get("APEX_P3_NATIVE_DISTANCE", "1") == "1":
            try:
                from core import _cpp  # noqa: PLC0415
                self.native = _cpp.module().MeshDistance(np.asarray(self.mesh.vertices, dtype=float),
                                                     np.asarray(self.mesh.faces, dtype=np.int64))
            except Exception:  # noqa: BLE001 -- module absent ou ancien : voie trimesh
                self.native = None
        edges = original[cells][:, [0, 0, 0, 1, 1, 2]] - original[cells][:, [1, 2, 3, 2, 3, 3]]
        self.tet_size = np.linalg.norm(edges, axis=2).max(axis=1)   # plus longue arete, par tetraedre

    def psi(self, xy, h, tets):
        from bird_to_apex import barycentric
        points = np.column_stack([xy, np.full(len(xy), h)])
        with np.errstate(all="ignore"):
            bary = barycentric(points, self.deformed[self.cells[tets]])
        real = np.einsum("nv,nvk->nk", bary, self.original[self.cells[tets]])
        if self.native is not None:
            return self.native.query(np.ascontiguousarray(real), 1, _NATIVE_THREADS)[0] - self.thickness, real
        distance = np.empty(len(real))
        for start in range(0, len(real), DISTANCE_BATCH):
            chunk = real[start:start + DISTANCE_BATCH]
            _c, distance[start:start + DISTANCE_BATCH], _t = trimesh.proximity.closest_point(self.mesh, chunk)
        return distance - self.thickness, real


def _split4(tri):
    """Chaque triangle (n,3,2) coupe en 4 par les milieux de ses cotes."""
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    return np.concatenate([np.stack([a, ab, ca], 1), np.stack([ab, b, bc], 1),
                           np.stack([ca, bc, c], 1), np.stack([ab, bc, ca], 1)])


def carved_region(cells, deformed, psi, h, exact=None, step=None, close=True):
    """Region de matiere du CORPS dans le plan z = h (espace deforme).

    Sans `exact` : psi interpole lineairement depuis les noeuds (exact dans
    chaque tetraedre, mais une frontiere a 0,2 mm de profondeur est mal
    suivie par des tetraedres de 1 a 5 mm, et l'ame d'une paroi mince lissee
    des deux cotes disparait : tous ses noeuds sont "dans la peau").
    Avec `exact` (ExactDistance) : pres de la peau, chaque section est coupee
    en triangles d'arete <= step et psi y est la VRAIE distance, calculee en
    leurs sommets ramenes dans la piece."""
    pieces, fine_tri, fine_tet = [], [], []
    for xy, pv, tet in _sections(cells, deformed, psi, h):
        near = exact is not None and pv.min() < exact.tet_size[tet]
        if near:
            # section convexe -> eventail de triangles, affine ensuite par la
            # boucle adaptative (jamais de subdivision uniforme a priori)
            fan = np.stack([np.repeat(xy[:1], len(xy) - 2, 0), xy[1:-1], xy[2:]], axis=1)
            fine_tri.append(fan)
            fine_tet.append(np.full(len(fan), tet))
        elif (pv >= 0).all():
            pieces.append(Polygon(xy))
        elif (pv < 0).all():
            continue
        else:
            clipped = _clip(xy, pv)
            if clipped is not None:
                pieces.append(Polygon(clipped))
    if fine_tri:
        tris = np.concatenate(fine_tri)
        tets = np.concatenate(fine_tet)
        # ADAPTATIF : la distance est 1-lipschitzienne. Un triangle dont les
        # trois sommets sont du meme cote de la frontiere, et plus loin d'elle
        # que sa plus grande arete (dans la piece REELLE, ou la distance est
        # mesuree), ne peut pas etre traverse par elle : il est garde ou ecarte
        # d'un bloc. Seuls les autres sont redivises -- le travail suit la
        # LONGUEUR de la frontiere, pas l'aire (la subdivision uniforme a 0,2 mm
        # depassait 10 min sur dino).
        for _depth in range(12):
            if not len(tris):
                break
            check_memory(f"couche z={h:.3f}, profondeur {_depth}")
            if len(tris) > MAX_REFINE_TRIANGLES:
                raise MemoryError(f"couche z={h:.3f} : {len(tris)} triangles a raffiner "
                                  f"(> {MAX_REFINE_TRIANGLES}) -- arret plutot que saturer la memoire")
            values, real = exact.psi(tris.reshape(-1, 2), h, np.repeat(tets, 3))
            values, real = values.reshape(-1, 3), real.reshape(-1, 3, 3)
            longest = np.linalg.norm(real - np.roll(real, 1, axis=1), axis=2).max(axis=1)
            inside = values.min(axis=1) > longest
            outside = values.max(axis=1) < -longest
            small = np.linalg.norm(tris - np.roll(tris, 1, axis=1), axis=2).max(axis=1) <= step
            pieces.extend(shapely.polygons(tris[inside]))
            final = small & ~inside & ~outside
            for tri, pv in zip(tris[final], values[final]):
                if (pv >= 0).all():
                    pieces.append(Polygon(tri))
                elif (pv >= 0).any():
                    clipped = _clip(tri, pv)
                    if clipped is not None:
                        pieces.append(Polygon(clipped))
            refine = ~inside & ~outside & ~small
            tris = _split4(tris[refine])
            tets = np.tile(tets[refine], 4)
    if not pieces:
        return MultiPolygon()
    region = shapely.union_all(pieces)
    if not close:
        return region
    # micro-fentes entre triangles de deux tetraedres voisins (subdivisions
    # differentes, raccords pas au bit pres) : resorbees en arrondissant les
    # coordonnees a 1e-6 mm. PAS de fermeture buffer(+e).buffer(-e) en onglet :
    # mesuree sur tilted benchy (couche 157), elle AJOUTAIT 4,7 mm2 de matiere
    # et bouchait un vrai trou de 3 mm2 (angles aigus des morceaux).
    return shapely.set_precision(region, 1e-6)


def region_loops(region, h):
    """Boucles fermees (x, y, h), convention de step5 (point de fermeture duplique)."""
    loops = []
    for piece in _polygon_pieces(region):
        for ring in (piece.exterior, *piece.interiors):
            xy = np.asarray(ring.coords)
            if len(xy) >= 4:
                loops.append(np.column_stack([xy, np.full(len(xy), h)]))
    return loops


def first_layer_toolpath(loops_2d, pattern, nozzle_width, wall_count, degenerate_area_ratio):
    """Premiere couche PLEINE : parois de l'etape 6, puis le motif choisi
    (surface_patterns) dans la zone interieure, a une buse d'ecart."""
    import surface_patterns
    from pipeline.step6_toolpath2d import _classify_contours

    paths = []
    for outer, holes in _classify_contours(loops_2d):
        walls, innermost = _offset_walls(outer, holes, nozzle_width, wall_count, degenerate_area_ratio)
        zone = innermost.buffer(-nozzle_width / 2.0, join_style="mitre") if innermost is not None else None
        fill = surface_patterns.generate(pattern, zone, nozzle_width) if zone is not None else []
        paths.extend(fill)
        paths.extend(walls)
    return paths


# --- creusement parallele par couche (refonte P3 du 2026-09-30) ---------------
_WORKER = {}


def carve_worker_init(skin, cells, deformed, original, psi, thickness):
    """Initialisation d'un processus : une ExactDistance par processus, et un
    seul thread pour la distance native (deja un processus par coeur)."""
    global _NATIVE_THREADS
    _NATIVE_THREADS = 1
    _WORKER.update(cells=cells, deformed=deformed, psi=psi, thickness=thickness,
                   exact=ExactDistance(skin, cells, deformed, original, thickness))


def carve_worker_layer(isovalue):
    w = _WORKER
    region = carved_region(w["cells"], w["deformed"], w["psi"], isovalue, exact=w["exact"], step=w["thickness"])
    return region_loops(region, isovalue)
