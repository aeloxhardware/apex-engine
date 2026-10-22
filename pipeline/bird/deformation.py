"""
Deformation de Bird : champ de rotation par tetraedre, puis deformation du
maillage qui le realise au mieux (moindres carres).

Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer),
cellules 2, 3, 4 et 7 du notebook d'origine, via research/S4_Slicer/bird_core.py.

Portage dans Apex (2026-10-07) : SANS pyvista. Seule la TOPOLOGIE change
d'implementation (voisinages, surface exterieure, centres de cellules, calcules
avec numpy) ; les problemes d'optimisation (`least_squares`, residus,
jacobiennes) sont repris a l'identique. Particularites conservees VOLONTAIREMENT,
parce qu'elles definissent le comportement mesure et valide par l'utilisateur :
  - les residus sont des CARRES (le solveur minimise donc des puissances 4) ;
  - les jacobiennes sont calculees en float32.

Representation : `BirdTet` (points (N,3), cells (T,4)) ; `BirdTopology`, calcule
une fois par maillage ; `BirdAttributes`, recalcule apres chaque deformation.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field, replace
from typing import Callable, Optional

import networkx as nx
import numpy as np
import scipy.sparse as sp
import trimesh
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

UP_VECTOR = np.array([0, 0, 1])
IN_AIR_THRESHOLD = 1
BOTTOM_CELL_TOLERANCE = 0.3
_N = np.eye(4) - 1 / 4 * np.ones((4, 4))
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ProgressCallback = Callable[[str, int, int, float], None]


# ---------------------------------------------------------------------------
# Donnees
# ---------------------------------------------------------------------------

@dataclass
class BirdTet:
    points: np.ndarray          # (N,3)
    cells: np.ndarray           # (T,4)

    @property
    def number_of_cells(self) -> int:
        return len(self.cells)

    @property
    def number_of_points(self) -> int:
        return len(self.points)

    def cell_centers(self) -> np.ndarray:
        return self.points[self.cells].mean(axis=1)

    def with_points(self, points: np.ndarray) -> "BirdTet":
        return BirdTet(points=np.asarray(points, dtype=float), cells=self.cells)


@dataclass
class BirdTopology:
    """Tout ce qui ne change pas d'une passe a l'autre (equivalent de ModelState
    et des field_data topologiques de la recherche)."""
    cell_neighbours: dict            # {"point"|"edge"|"face": (P,2) paires i<j}
    cell_neighbour_dict: dict        # {"point"|"edge"|"face": {cellule: [voisins]}}
    cell_neighbour_graph: nx.Graph   # aretes "point", poids = distance entre centres d'ORIGINE
    surface_faces: np.ndarray        # (F,3) triangles exterieurs, orientes vers l'exterieur
    cell_to_face: dict               # {cellule: [indices de surface_faces]}
    bottom_cells: Optional[np.ndarray] = None
    bottom_cells_mask: Optional[np.ndarray] = None


@dataclass
class BirdAttributes:
    """Attributs par cellule recalcules a partir de la geometrie COURANTE
    (update_tet_attributes de la recherche)."""
    face_normal: np.ndarray
    face_center: np.ndarray
    cell_center: np.ndarray
    is_bottom: np.ndarray
    overhang_angle: np.ndarray
    overhang_direction: np.ndarray
    in_air: np.ndarray
    extra: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Ingestion
# ---------------------------------------------------------------------------

def tetrahedralize(vertices: np.ndarray, faces: np.ndarray, max_tet_volume: Optional[float] = None,
                   keep_boundary: bool = False) -> BirdTet:
    """TetGen en sous-processus. Memes switches que la recherche
    (`pq1.414a<v>`, ou les defauts du binding). `keep_boundary` ajoute `Y` : aucun point
    de Steiner sur la surface. Remede des surfaces a details quasi coincidents, ou TetGen,
    en decoupant les segments de bord, fabriquait des millions de tetraedres de 1e-18 mm3
    (Fighter jet F18, 2026-10-18) ; `max_tet_volume` est alors obligatoire."""
    switches = f"pq1.414a{max_tet_volume}" if max_tet_volume is not None else ""
    if keep_boundary:
        switches = f"pYq1.414a{max_tet_volume}"
    with tempfile.TemporaryDirectory() as tmp:
        src, dst = os.path.join(tmp, "in.npz"), os.path.join(tmp, "out.npz")
        np.savez(src, vertices=np.asarray(vertices, dtype=float), faces=np.asarray(faces, dtype=np.int64),
                 switches=np.array(switches))
        try:
            proc = subprocess.run([sys.executable, "-m", "pipeline.bird._tetgen_worker", src, dst],
                                  cwd=_ROOT, capture_output=True, text=True, timeout=TETGEN_TIMEOUT_S)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"TetGen (voie Bird) depasse {TETGEN_TIMEOUT_S} s (maillage degenere probable)")
        if proc.returncode != 0 or not os.path.exists(dst):
            raise RuntimeError(f"TetGen (voie Bird) a echoue, code {proc.returncode} : {proc.stderr[-2000:]}")
        with np.load(dst) as out:       # fermeture explicite : sous Windows, un .npz ouvert bloque la suppression
            return BirdTet(points=np.array(out["node"], dtype=float), cells=np.array(out["elem"], dtype=np.int64))


# ---------------------------------------------------------------------------
# Autodocteur : diagnostic systematique AVANT TetGen, controle du maillage APRES
# ---------------------------------------------------------------------------

# Un nœud partage par des dizaines de milliers de tetraedres rend `neighbour_pairs`
# quadratique en memoire (Fighter jet F18, 2026-10-18 : valence max 88 568, 17 milliards
# d'entrees, 125 Go demandes) alors que TetGen avait ACCEPTE la surface (146 faces
# auto-intersectantes, 2 coques qui se traversent). Seuils : max_valence absolu, et somme
# des valences^2 (cout memoire du produit d'incidence) bornee par cellule.
MAX_NODE_VALENCE = 2000
MAX_VALENCE_SQUARES_PER_CELL = 300.0


@dataclass
class MeshDiagnosis:
    faces: int
    watertight: bool
    winding_consistent: bool
    shells: int
    self_intersecting_faces: int
    degenerate_faces: int

    @property
    def problems(self) -> list[str]:
        found = []
        if not self.watertight:
            found.append("non etanche")
        if not self.winding_consistent:
            found.append("orientation incoherente")
        if self.self_intersecting_faces:
            found.append(f"{self.self_intersecting_faces} faces auto-intersectantes")
        if self.degenerate_faces:
            found.append(f"{self.degenerate_faces} faces degenerees")
        if self.shells > 1 and self.self_intersecting_faces:
            found.append(f"{self.shells} coques qui se traversent")
        return found

    @property
    def needs_repair(self) -> bool:
        return bool(self.problems)


def diagnose_surface(mesh: trimesh.Trimesh) -> MeshDiagnosis:
    """Analyse la surface AVANT toute tetraedrisation (jamais seulement apres un refus)."""
    from pipeline import step1_ingest as s1

    merged = mesh.copy()
    merged.merge_vertices()
    return MeshDiagnosis(
        faces=len(merged.faces),
        watertight=bool(merged.is_watertight),
        winding_consistent=bool(merged.is_winding_consistent),
        shells=len(merged.split(only_watertight=False)),
        self_intersecting_faces=len(s1._find_self_intersecting_faces(merged)),
        degenerate_faces=int((merged.area_faces < 1e-12).sum()),
    )


WELD_RELATIVE_TOLERANCE = 5e-4      # fraction de la plus grande dimension (Fighter : ~0,08 mm)
WELD_MAX_VOLUME_CHANGE = 0.05       # une soudure qui change le volume de plus de 5 % est rejetee
TETGEN_TIMEOUT_S = 300


def weld_short_edges(mesh: trimesh.Trimesh, rel_tol: float = WELD_RELATIVE_TOLERANCE,
                     max_volume_change: float = WELD_MAX_VOLUME_CHANGE) -> Optional[trimesh.Trimesh]:
    """Nettoyage de surface : fusionne les sommets distants de moins de `rel_tol` x la plus
    grande dimension (aretes quasi nulles, tôles d'epaisseur ~0, bords de fuite), retire les
    faces devenues degenerees ou en double. TetGen produisait des millions de tetraedres de
    1e-18 mm3 autour de ces details (Fighter jet F18, 2026-10-18). None si le resultat
    change le volume de plus de `max_volume_change`, rapetisse l'etendue, ou n'est plus un
    solide valide : le modele n'est alors jamais deforme en silence."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    vertices = np.asarray(mesh.vertices, dtype=float)
    tolerance = rel_tol * float(np.max(mesh.extents))
    pairs = cKDTree(vertices).query_pairs(tolerance, output_type="ndarray")
    n = len(vertices)
    graph = coo_matrix((np.ones(len(pairs), dtype=np.int8), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    count, labels = connected_components(graph, directed=False)
    if count == n:
        return mesh                                   # rien a souder
    sums = np.zeros((count, 3))
    np.add.at(sums, labels, vertices)
    merged = sums / np.bincount(labels, minlength=count)[:, None]
    faces = labels[np.asarray(mesh.faces)]
    faces = faces[(faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & (faces[:, 0] != faces[:, 2])]
    _, first = np.unique(np.sort(faces, axis=1), axis=0, return_index=True)   # faces en double
    faces = faces[np.sort(first)]
    welded = trimesh.Trimesh(vertices=merged, faces=faces, process=False)
    welded.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(welded)
    if not welded.is_watertight:
        trimesh.repair.fill_holes(welded)
    if not (welded.is_watertight and welded.is_winding_consistent and welded.volume > 0):
        return None
    if abs(welded.volume - mesh.volume) > max_volume_change * abs(mesh.volume):
        return None
    if np.any(np.asarray(welded.extents) < 0.99 * np.asarray(mesh.extents)):
        return None
    return welded


def tet_mesh_problem(tet: BirdTet) -> Optional[str]:
    """None si le maillage tetraedrique est exploitable, sinon la raison."""
    valence = np.bincount(np.asarray(tet.cells).ravel(), minlength=len(tet.points)).astype(float)
    squares = float((valence ** 2).sum())
    if valence.max() > MAX_NODE_VALENCE:
        return f"valence maximale {int(valence.max())} (> {MAX_NODE_VALENCE})"
    if squares > MAX_VALENCE_SQUARES_PER_CELL * len(tet.cells):
        return f"somme des valences^2 = {squares:.3g} (> {MAX_VALENCE_SQUARES_PER_CELL:g} par cellule)"
    return None


def repair_candidates(mesh: trimesh.Trimesh) -> tuple[trimesh.Trimesh, list[tuple[str, trimesh.Trimesh]]]:
    """Surfaces de remplacement candidates, dans l'ordre de preference, chacune avec son
    libelle. Rend `(maillage nettoye, candidats)` ; le garde-fou d'etendue (une reparation
    qui ampute la piece est rejetee) est applique par l'appelant. Partagee par la voie Bird
    et par l'ingestion FEM (`step1_ingest`) : une seule cascade de reparation. L'ordre est
    celui de la voie Bird d'origine (identique au bit)."""
    from pipeline import step1_ingest as s1

    candidates: list[tuple[str, trimesh.Trimesh]] = []
    cleaned, _report = s1._clean_and_repair(mesh.copy())

    # Nettoyage de surface d'abord : souder les details quasi coincidents, puis l'union des
    # coques qui se traversent. Un candidat soude est ensuite essaye avant les autres.
    welded = weld_short_edges(cleaned)
    if welded is not None and welded is not cleaned:
        w_bodies = welded.split(only_watertight=True)
        if len(w_bodies) > 1:
            try:
                candidates.append(("soudure des aretes courtes + union des coques",
                                   trimesh.boolean.union(list(w_bodies), engine="manifold")))
            except Exception:      # noqa: BLE001 -- best-effort
                pass
        candidates.append(("soudure des aretes courtes", welded))
    candidates.append(("nettoyage seul", cleaned))
    # Coques qui se traversent (corps concatenes sans union) : l'union booleenne les
    # fusionne sans rien perdre, la ou pymeshfix garde UNE coque (mesure : deux boites
    # croisees -> une seule, a exactement 50 % des faces, sous le garde-fou de step1).
    bodies = cleaned.split(only_watertight=True)
    if len(bodies) > 1:
        try:
            candidates.append(("union des coques qui se traversent",
                               trimesh.boolean.union(list(bodies), engine="manifold")))
        except Exception:      # noqa: BLE001 -- best-effort, la cascade continue
            pass
    resolved = s1.resolve_self_intersections(cleaned)
    if resolved is not None and len(resolved.faces) >= 0.5 * len(cleaned.faces):
        if s1._find_self_intersecting_faces(resolved):
            local = s1._repair_local_self_intersections(resolved, require_euler_preserved=False)
            if local is not None:
                resolved = local
        resolved, _probe = s1._finish_geometric_resolution(resolved, 1.0)
        candidates.append(("resolution geometrique des auto-intersections", resolved))
    fixed = s1._repair_with_pymeshfix(cleaned)
    if fixed is not None and len(fixed.faces) >= 0.5 * len(cleaned.faces):
        candidates.append(("pymeshfix", fixed))
    return cleaned, candidates


def _tetrahedralize_repaired(mesh: trimesh.Trimesh, max_tet_volume: Optional[float],
                             first_error: Exception) -> BirdTet:
    """Reparation de la surface (autodocteur) : la cascade de reparation de step1, dans le
    meme ordre de preference (nettoyage, resolution geometrique des auto-intersections avec
    TetGen pour oracle, puis pymeshfix), chaque candidat rejete s'il perd plus de la moitie
    des faces, ou si son maillage tetraedrique est inexploitable (`tet_mesh_problem`).
    Appelee quand `diagnose_surface` trouve un defaut, ou quand TetGen refuse ou produit
    un maillage inexploitable. Un maillage sain n'est jamais modifie (Benchy inchange)."""
    cleaned, labelled = repair_candidates(mesh)
    candidates = [candidate for _label, candidate in labelled]
    errors = [first_error]
    reference = cleaned.extents
    for candidate in candidates:
        # Garde-fou d'effondrement plus strict que le nombre de faces : une reparation
        # qui ampute une partie de la piece (turkish, limite step1 ; corps perdu par
        # pymeshfix) en reduit l'etendue.
        if np.any(np.asarray(candidate.extents) < 0.99 * reference):
            errors.append(RuntimeError(f"reparation rejetee : etendue {np.round(candidate.extents, 2)} "
                                       f"contre {np.round(reference, 2)}"))
            continue
        # Deux essais par candidat : TetGen normal, puis sans point de Steiner sur la surface
        # (`Y`, volume impose : ~30 000 tetraedres a defaut de `max_tet_volume`).
        keep_volume = max_tet_volume if max_tet_volume is not None else abs(float(candidate.volume)) / 30000.0
        for keep_boundary in (False, True):
            try:
                tet = tetrahedralize(np.asarray(candidate.vertices), np.asarray(candidate.faces),
                                     keep_volume if keep_boundary else max_tet_volume, keep_boundary)
            except RuntimeError as exc:
                errors.append(exc)
                continue
            problem = tet_mesh_problem(tet)
            if problem is not None:
                errors.append(RuntimeError(f"candidat rejete{' (Y)' if keep_boundary else ''} : "
                                           f"maillage tetraedrique inexploitable ({problem})"))
                continue
            if keep_boundary:
                print(f"[bird] autodocteur : TetGen sans point de Steiner en surface (Y), "
                      f"{len(tet.cells)} tetraedres", flush=True)
            return tet
    raise RuntimeError("TetGen (voie Bird) refuse le maillage meme apres reparation : "
                       + " | ".join(str(e)[-300:] for e in errors))


def load_and_tetrahedralize(stl_path: str, max_tet_volume: Optional[float] = None, surface_subdivisions: int = 0,
                            flatten_bottom: Optional[float] = None) -> BirdTet:
    """Charge, aplatit le dessous (option), tetraedrise, recentre en XY et pose
    au sol en Z (bounds du maillage tetraedrique, comme la recherche)."""
    mesh = trimesh.load(stl_path)
    if flatten_bottom is not None:
        from pipeline.bird.mesh_flatten import flatten_bottom as _flatten_bottom
        mesh, _displacement = _flatten_bottom(mesh, flatten_bottom)
    for _ in range(surface_subdivisions):
        mesh = mesh.subdivide()
    # Autodocteur : la surface est TOUJOURS analysee ; un defaut declenche la reparation
    # avant TetGen. Apres TetGen, le maillage est controle (valences) : un maillage
    # inexploitable declenche la meme reparation au lieu d'une erreur memoire plus loin.
    diagnosis = diagnose_surface(mesh)
    if diagnosis.needs_repair:
        print(f"[bird] autodocteur : {', '.join(diagnosis.problems)} -> reparation avant TetGen", flush=True)
        tet = _tetrahedralize_repaired(mesh, max_tet_volume, RuntimeError("defaut de surface detecte"))
    else:
        try:
            tet = tetrahedralize(np.asarray(mesh.vertices), np.asarray(mesh.faces), max_tet_volume)
            problem = tet_mesh_problem(tet)
            if problem is not None:
                print(f"[bird] autodocteur : maillage tetraedrique inexploitable ({problem}) -> reparation",
                      flush=True)
                tet = _tetrahedralize_repaired(mesh, max_tet_volume, RuntimeError(problem))
        except RuntimeError as first_error:
            tet = _tetrahedralize_repaired(mesh, max_tet_volume, first_error)
    lo, hi = tet.points.min(axis=0), tet.points.max(axis=0)
    tet.points = tet.points - np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, lo[2]])
    return tet


# ---------------------------------------------------------------------------
# Topologie (remplace pyvista : cell_neighbors, extract_surface, extract_points)
# ---------------------------------------------------------------------------

_TET_FACES = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
_TET_OPPOSITE = np.array([3, 2, 1, 0])


def _shared_node_counts(cells: np.ndarray, n_points: int) -> sp.csr_matrix:
    t = len(cells)
    inc = sp.csr_matrix((np.ones(t * 4, dtype=np.int32), (np.repeat(np.arange(t), 4), cells.ravel())),
                        shape=(t, n_points))
    return (inc @ inc.T).tocsr()


def neighbour_pairs(cells: np.ndarray, n_points: int) -> dict:
    """Paires (i<j) de cellules partageant au moins 1 noeud ("point"), 2 noeuds
    (une arete, "edge") ou 3 noeuds (une face, "face") -- meme definition que
    `cell_neighbors(i, "<type>s")` de VTK. Triees par (i, j)."""
    counts = sp.triu(_shared_node_counts(cells, n_points), k=1).tocoo()
    order = np.lexsort((counts.col, counts.row))
    rows, cols, shared = counts.row[order], counts.col[order], counts.data[order]
    out = {}
    for name, minimum in (("point", 1), ("edge", 2), ("face", 3)):
        keep = shared >= minimum
        out[name] = np.column_stack([rows[keep], cols[keep]]).astype(np.int64)
    return out


def surface_faces_and_owners(tet: BirdTet) -> tuple[np.ndarray, np.ndarray]:
    """Triangles exterieurs (faces appartenant a une seule cellule), orientes
    vers l'exterieur sur la geometrie de `tet`, et leur cellule proprietaire.
    L'orientation est FIXEE ici, sur la geometrie d'origine, puis reutilisee
    telle quelle apres deformation (comme VTK, dont l'orientation vient de la
    connectivite et non de la geometrie courante)."""
    cells = tet.cells
    tri = cells[:, _TET_FACES].reshape(-1, 3)
    opposite = cells[:, _TET_OPPOSITE].ravel()
    owner = np.repeat(np.arange(len(cells)), 4)
    key = np.sort(tri, axis=1)
    _, inverse, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    boundary = counts[inverse.ravel()] == 1
    tri, opposite, owner = tri[boundary], opposite[boundary], owner[boundary]
    p = tet.points
    normal = np.cross(p[tri[:, 1]] - p[tri[:, 0]], p[tri[:, 2]] - p[tri[:, 0]])
    inward = np.einsum("ij,ij->i", normal, p[opposite] - p[tri[:, 0]]) > 0
    tri[inward] = tri[inward][:, [0, 2, 1]]
    order = np.argsort(owner, kind="stable")
    return tri[order], owner[order]


def build_topology(tet: BirdTet) -> BirdTopology:
    pairs = neighbour_pairs(tet.cells, tet.number_of_points)
    neighbour_dict = {}
    for name, p in pairs.items():
        d = {c: [] for c in range(tet.number_of_cells)}
        for a, b in p:
            d[int(a)].append(int(b))
            d[int(b)].append(int(a))
        neighbour_dict[name] = d
    graph = nx.Graph()
    centers = tet.cell_centers()
    pts = pairs["point"]
    distances = np.linalg.norm(centers[pts[:, 0]] - centers[pts[:, 1]], axis=1)
    graph.add_weighted_edges_from((int(a), int(b), float(w)) for (a, b), w in zip(pts, distances))
    faces, owners = surface_faces_and_owners(tet)
    cell_to_face: dict = {}
    for face_index, cell in enumerate(owners):
        cell_to_face.setdefault(int(cell), []).append(face_index)
    return BirdTopology(cell_neighbours=pairs, cell_neighbour_dict=neighbour_dict, cell_neighbour_graph=graph,
                        surface_faces=faces, cell_to_face=cell_to_face)


def surface_normals_and_centers(tet: BirdTet, topo: BirdTopology) -> tuple[np.ndarray, np.ndarray]:
    tri = tet.points[topo.surface_faces]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.linalg.norm(n, axis=1)[:, None]
    return n, tri.mean(axis=1)


# ---------------------------------------------------------------------------
# CELLULE 2 -- attributs (update_tet_attributes / calculate_tet_attributes)
# ---------------------------------------------------------------------------

def update_tet_attributes(tet: BirdTet, topo: BirdTopology) -> BirdAttributes:
    face_normals_all, face_centers_all = surface_normals_and_centers(tet, topo)
    n_cells = tet.number_of_cells
    face_normal = np.full((n_cells, 3), np.nan)
    face_center = np.full((n_cells, 3), np.nan)
    for cell_index, face_indices in topo.cell_to_face.items():
        fn = face_normals_all[face_indices]
        face_normal[cell_index] = fn[np.argmin(fn[:, 2])]
        fc = face_centers_all[face_indices]
        face_center[cell_index] = fc[np.argmin(fc[:, 2])]
    face_normal = face_normal / np.linalg.norm(face_normal, axis=1)[:, None]
    cell_center = tet.cell_centers()

    bottom_cell_threshold = np.nanmin(face_center[:, 2]) + BOTTOM_CELL_TOLERANCE
    is_bottom = face_center[:, 2] < bottom_cell_threshold
    bottom_cells_local = np.where(is_bottom)[0]

    normals = face_normal.copy()
    normals[is_bottom] = np.nan
    with np.errstate(invalid="ignore"):
        overhang_angle = np.arccos(np.dot(normals, UP_VECTOR))
        overhang_direction = normals[:, :2].copy()
        overhang_direction /= np.linalg.norm(overhang_direction, axis=1)[:, None]

    in_air = np.full(n_cells, False)
    _, paths_to_bottom = nx.multi_source_dijkstra(topo.cell_neighbour_graph, set(bottom_cells_local.tolist()))
    for cell_index in range(n_cells):
        path = paths_to_bottom[cell_index]
        if len(path) > 1:
            if np.any(cell_center[path, 2] > cell_center[cell_index, 2] + IN_AIR_THRESHOLD):
                in_air[cell_index] = True

    return BirdAttributes(face_normal=face_normal, face_center=face_center, cell_center=cell_center,
                          is_bottom=is_bottom, overhang_angle=overhang_angle,
                          overhang_direction=overhang_direction, in_air=in_air)


def calculate_tet_attributes(tet: BirdTet) -> tuple[BirdTopology, BirdAttributes]:
    """Topologie (une fois) + attributs de la geometrie d'origine ; fixe les
    cellules du bas, dont le surplomb est neutralise (NaN)."""
    topo = build_topology(tet)
    attrs = update_tet_attributes(tet, topo)
    topo.bottom_cells_mask = attrs.is_bottom.copy()
    topo.bottom_cells = np.where(topo.bottom_cells_mask)[0]
    attrs.overhang_angle[topo.bottom_cells] = np.nan
    return topo, attrs


# ---------------------------------------------------------------------------
# CELLULE 3 -- gradient de la longueur de chemin vers la base (inchange)
# ---------------------------------------------------------------------------

def plane_fit(points):
    points = np.reshape(points, (np.shape(points)[0], -1))
    ctr = points.mean(axis=1)
    x = points - ctr[:, np.newaxis]
    m = np.dot(x, x.T)
    return ctr, np.linalg.svd(m)[0][:, -1]


def calculate_path_length_to_base_gradient(attrs: BirdAttributes, topo: BirdTopology, MAX_OVERHANG,
                                           INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO):
    cell_neighbour_dict = topo.cell_neighbour_dict
    bottom_cells = topo.bottom_cells
    bottom_set = set(bottom_cells.tolist())
    n_cells = len(attrs.cell_center)
    path_length_to_base_gradient = np.zeros(n_cells)

    cell_distance_to_bottom = np.full(n_cells, np.nan)
    distances_to_bottom, paths_to_bottom = nx.multi_source_dijkstra(topo.cell_neighbour_graph, bottom_set)
    closest_bottom_cell_indices = np.zeros(n_cells, dtype=int)
    with np.errstate(invalid="ignore"):
        for cell_index in range(n_cells):
            face_normal = attrs.face_normal[cell_index]
            cell_is_overhang = np.arccos(np.dot(face_normal, [0, 0, 1])) > np.deg2rad(90 + MAX_OVERHANG)
            if cell_is_overhang and cell_index not in bottom_set:
                closest_bottom_cell_indices[cell_index] = paths_to_bottom[cell_index][0]
                cell_distance_to_bottom[cell_index] = distances_to_bottom[cell_index]

    attrs.extra["cell_distance_to_bottom"] = cell_distance_to_bottom
    centers = attrs.cell_center

    with np.errstate(invalid="ignore", divide="ignore"):
        for cell_index in range(n_cells):
            if np.isnan(cell_distance_to_bottom[cell_index]):
                continue
            local_cells = cell_neighbour_dict["edge"][cell_index]
            local_cells = np.hstack((local_cells, cell_index)).astype(int)
            local_lengths = cell_distance_to_bottom[local_cells]
            local_cells = local_cells[~np.isnan(local_lengths)]
            local_lengths = local_lengths[~np.isnan(local_lengths)]

            if len(local_lengths) < 3:
                location_to_roll_to = centers[closest_bottom_cell_indices[cell_index], :2]
                direction_to_bottom = location_to_roll_to - centers[cell_index, :2]
                direction_to_bottom /= np.linalg.norm(direction_to_bottom)
                cell_center = centers[cell_index, :2].copy()
                cell_center /= np.linalg.norm(cell_center)
                optimal = np.dot(cell_center, direction_to_bottom) / np.abs(np.dot(cell_center, direction_to_bottom))
                if np.isnan(optimal):
                    optimal = 0
                path_length_to_base_gradient[cell_index] = optimal
            else:
                points = np.hstack((centers[local_cells, :2], local_lengths[:, None]))
                _, plane_normal = plane_fit(points.T)
                radial = centers[cell_index, :2] / np.linalg.norm(centers[cell_index, :2])
                gradient = np.dot(radial, plane_normal[:2])
                if np.isnan(gradient):
                    neighbours = path_length_to_base_gradient[local_cells]
                    gradient = np.mean(neighbours[~np.isnan(neighbours)])
                    if np.isnan(gradient):
                        gradient = 0
                path_length_to_base_gradient[cell_index] = gradient

    if INITIAL_ROTATION_FIELD_SMOOTHING != 0:
        # Dans l'original, chaque tour repart du champ NON lisse (la variable
        # lissee n'est reaffectee qu'apres la boucle) : les N tours rendent
        # exactement le meme resultat. Un seul calcul, resultat identique.
        for _ in range(1):
            smoothed = np.zeros(n_cells)
            for cell_index in range(n_cells):
                if path_length_to_base_gradient[cell_index] != 0:
                    neighbours = cell_neighbour_dict["point"][cell_index]
                    local_cells = neighbours.copy()
                    for neighbour in neighbours:
                        local_cells.extend(cell_neighbour_dict["point"][neighbour])
                    local_cells = np.array(list(set(local_cells)))
                    local_cells = local_cells[path_length_to_base_gradient[local_cells] != 0]
                    smoothed[cell_index] = np.mean(path_length_to_base_gradient[local_cells])
        path_length_to_base_gradient = smoothed

    if not SET_INITIAL_ROTATION_TO_ZERO:
        path_length_to_base_gradient[path_length_to_base_gradient == 0] = np.nan
    attrs.extra["path_length_to_base_gradient"] = path_length_to_base_gradient
    return path_length_to_base_gradient


# ---------------------------------------------------------------------------
# CELLULE 4 -- champ de rotation initial et matrices (inchange)
# ---------------------------------------------------------------------------

def calculate_initial_rotation_field(attrs, topo, MAX_OVERHANG, ROTATION_MULTIPLIER, STEEP_OVERHANG_COMPENSATION,
                                     INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO,
                                     MAX_POS_ROTATION, MAX_NEG_ROTATION):
    initial = np.abs(np.deg2rad(90 + MAX_OVERHANG) - attrs.overhang_angle)
    gradient = calculate_path_length_to_base_gradient(attrs, topo, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING,
                                                      SET_INITIAL_ROTATION_TO_ZERO)
    if STEEP_OVERHANG_COMPENSATION:
        initial[attrs.in_air] += 2 * (np.deg2rad(180) - attrs.overhang_angle[attrs.in_air])
    initial *= gradient
    initial = np.clip(initial * ROTATION_MULTIPLIER, -np.deg2rad(360), np.deg2rad(360))
    initial = np.clip(initial, MAX_NEG_ROTATION, MAX_POS_ROTATION)
    attrs.extra["initial_rotation_field"] = initial
    return initial


def calculate_rotation_matrices(cell_center: np.ndarray, rotation_field: np.ndarray) -> np.ndarray:
    xy_3d = np.hstack([cell_center[:, :2], np.zeros((len(cell_center), 1))])
    tangential = np.cross(np.array([0, 0, 1]), xy_3d)
    with np.errstate(invalid="ignore"):
        tangential /= np.linalg.norm(tangential, axis=1)[:, None]
    tangential[np.isnan(tangential).any(axis=1)] = [1, 0, 0]
    return Rotation.from_rotvec(rotation_field[:, None] * tangential).as_matrix()


def planar_base_mask(tet: BirdTet, topo: BirdTopology, overhang_deg: float, layer_height: float,
                     max_height: Optional[float] = None, min_height: Optional[float] = None):
    """Base imprimee A PLAT : tout le volume sous le surplomb le plus bas (voir
    la recherche, bird_core.planar_base_mask, pour les mesures qui ont fixe ce
    choix). Rend (masque des cellules, masque des noeuds, z du plus bas surplomb)."""
    cells, points = tet.cells, tet.points
    tri = points[topo.surface_faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.maximum(np.linalg.norm(normals, axis=1), 1e-12)[:, None]
    z0 = points[:, 2].min()
    on_plate = (tri[:, :, 2].max(axis=1) < z0 + layer_height) & (normals[:, 2] < 0.0)
    angle = np.degrees(np.arccos(np.clip(normals[:, 2], -1.0, 1.0)))
    overhanging = (angle > 90.0 + overhang_deg) & ~on_plate
    if min_height is not None:
        # Chanfrein / arrondi de pied (Fighter jet F18 : 45 deg a 0,134 mm du plateau) : un
        # « surplomb » de quelques dixiemes de mm ne limite pas la base, il s'imprime sans support.
        overhanging &= tri[:, :, 2].min(axis=1) >= z0 + min_height
    z_first = float(tri[overhanging][:, :, 2].min()) if overhanging.any() else float(points[:, 2].max())
    top = max(z_first, z0 + layer_height, z0 + (min_height or 0.0))
    if max_height is not None:
        top = min(top, z0 + max(max_height, layer_height))
    node_mask = points[:, 2] <= top + 1e-6
    cell_mask = node_mask[cells].all(axis=1)
    return cell_mask, node_mask, z_first


# ---------------------------------------------------------------------------
# Optimisation des rotations (inchange)
# ---------------------------------------------------------------------------

def optimize_rotations(tet, topo, attrs, NEIGHBOUR_LOSS_WEIGHT, MAX_OVERHANG, ROTATION_MULTIPLIER, ITERATIONS,
                       STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO,
                       MAX_POS_ROTATION, MAX_NEG_ROTATION, progress_callback=None, fixed_zero_cells=None):
    initial = calculate_initial_rotation_field(attrs, topo, MAX_OVERHANG, ROTATION_MULTIPLIER,
                                               STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING,
                                               SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION)
    face_pairs = topo.cell_neighbours["face"]
    native = _native()
    if native is not None:
        fixed = (np.asarray(fixed_zero_cells, dtype=bool) if fixed_zero_cells is not None
                 else np.zeros(tet.number_of_cells, dtype=bool))
        free = np.flatnonzero(~fixed)
        column = np.full(tet.number_of_cells, -1, dtype=np.int64)
        column[free] = np.arange(len(free))
        pairs = face_pairs[~(fixed[face_pairs[:, 0]] & fixed[face_pairs[:, 1]])]
        targets = np.flatnonzero(~np.isnan(initial) & ~fixed)
        x, cost, _nfev, _status = native.bird_lsq_rotation(
            np.ascontiguousarray(pairs, dtype=np.int64), targets.astype(np.int64),
            np.ascontiguousarray(np.nan_to_num(initial), dtype=float), column, float(NEIGHBOUR_LOSS_WEIGHT),
            int(ITERATIONS), 1e-6)
        field_ = np.zeros(tet.number_of_cells)
        field_[free] = x
        return field_, cost
    if fixed_zero_cells is not None and np.asarray(fixed_zero_cells).any():
        return _optimize_rotations_with_fixed(tet, face_pairs, initial, NEIGHBOUR_LOSS_WEIGHT, ITERATIONS,
                                              np.asarray(fixed_zero_cells, dtype=bool), progress_callback)
    n_cells = tet.number_of_cells
    n_targets = int(np.sum(~np.isnan(initial)))
    n_pairs = len(face_pairs)
    c1, c2 = face_pairs[:, 0], face_pairs[:, 1]
    valid = np.where(~np.isnan(initial))[0]
    rows = np.concatenate([np.arange(n_pairs), np.arange(n_pairs), n_pairs + np.arange(len(valid))])
    cols = np.concatenate([c1, c2, valid])
    shape = (n_pairs + n_targets, n_cells)

    def objective_function(field_):
        d = field_[c1] - field_[c2]
        return np.concatenate((NEIGHBOUR_LOSS_WEIGHT * d ** 2, (field_[valid] - initial[valid]) ** 2))

    def objective_jacobian(field_):
        d = field_[c1] - field_[c2]
        data = np.concatenate([2 * NEIGHBOUR_LOSS_WEIGHT * d, -2 * NEIGHBOUR_LOSS_WEIGHT * d,
                               2 * (field_[valid] - initial[valid])]).astype(np.float32)
        return sp.coo_matrix((data, (rows, cols)), shape=shape).tocsr()

    sparsity = sp.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=shape).tocsr()
    result = least_squares(objective_function, np.zeros(n_cells), jac=objective_jacobian, max_nfev=ITERATIONS,
                           jac_sparsity=sparsity, verbose=0, method="trf", ftol=1e-6,
                           callback=_progress("rotation", ITERATIONS, progress_callback))
    return result.x, result.cost


def _native():
    """Module natif si `APEX_USE_CPP=1` et compile avec les moindres carres de Bird."""
    from core import _cpp
    module = _cpp.module() if _cpp.enabled() else None
    return module if module is not None and hasattr(module, "bird_lsq_deform") else None


def _progress(phase, iterations, progress_callback):
    if progress_callback is None:
        return None

    def callback(intermediate_result):
        progress_callback(phase, intermediate_result.nfev, iterations, float(intermediate_result.cost))
    return callback


def _optimize_rotations_with_fixed(tet, face_pairs, initial, NEIGHBOUR_LOSS_WEIGHT, ITERATIONS, fixed,
                                   progress_callback=None):
    n_cells = tet.number_of_cells
    free = np.flatnonzero(~fixed)
    column = np.full(n_cells, -1)
    column[free] = np.arange(len(free))
    pairs = face_pairs[~(fixed[face_pairs[:, 0]] & fixed[face_pairs[:, 1]])]
    c1, c2 = pairs[:, 0], pairs[:, 1]
    targets = np.flatnonzero(~np.isnan(initial) & ~fixed)
    n_pairs = len(pairs)

    rows, cols, which = [], [], []
    for side, cells_ in ((0, c1), (1, c2)):
        keep = column[cells_] >= 0
        rows.append(np.flatnonzero(keep))
        cols.append(column[cells_[keep]])
        which.append(np.full(keep.sum(), side))
    rows.append(n_pairs + np.arange(len(targets)))
    cols.append(column[targets])
    rows, cols, which = np.concatenate(rows), np.concatenate(cols), np.concatenate(which)
    shape = (n_pairs + len(targets), len(free))
    pair_rows = rows[: len(rows) - len(targets)]
    pair_which = which[: len(rows) - len(targets)]

    def full(x):
        f = np.zeros(n_cells)
        f[free] = x
        return f

    def objective_function(x):
        f = full(x)
        d = f[c1] - f[c2]
        return np.concatenate((NEIGHBOUR_LOSS_WEIGHT * d ** 2, (f[targets] - initial[targets]) ** 2))

    def objective_jacobian(x):
        f = full(x)
        d = f[c1] - f[c2]
        pair_values = np.where(pair_which == 0, 2 * NEIGHBOUR_LOSS_WEIGHT, -2 * NEIGHBOUR_LOSS_WEIGHT) * d[pair_rows]
        data = np.concatenate([pair_values, 2 * (f[targets] - initial[targets])]).astype(np.float32)
        return sp.coo_matrix((data, (rows, cols)), shape=shape).tocsr()

    sparsity = sp.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=shape).tocsr()
    result = least_squares(objective_function, np.zeros(len(free)), jac=objective_jacobian, max_nfev=ITERATIONS,
                           jac_sparsity=sparsity, verbose=0, method="trf", ftol=1e-6,
                           callback=_progress("rotation", ITERATIONS, progress_callback))
    return full(result.x), result.cost


# ---------------------------------------------------------------------------
# CELLULE 7 -- deformation (inchange)
# ---------------------------------------------------------------------------

def _deformation_native(native, tet, cell_center, rotation_field, ITERATIONS, node_mask, pinned_positions):
    """Meme probleme que _calculate_deformation_pinned (et, sans noeud epingle,
    que calculate_deformation), resolu en C++ (apex_core.bird_lsq_deform)."""
    cells = tet.cells
    base = tet.points.copy()
    base[node_mask] = pinned_positions[node_mask]
    free = np.flatnonzero(~node_mask)
    column = np.full(tet.number_of_points, -1, dtype=np.int64)
    column[free] = np.arange(len(free))
    rotation_matrices = calculate_rotation_matrices(cell_center, rotation_field)
    old = np.einsum("ijk,ikl->ijl", rotation_matrices, (_N @ tet.points[cells]).transpose(0, 2, 1))
    x, cost, _nfev, _status = native.bird_lsq_deform(base, np.ascontiguousarray(cells, dtype=np.int64), column,
                                                     np.ascontiguousarray(old), int(ITERATIONS))
    base[free] = np.asarray(x).reshape(-1, 3)
    return base, cost


def calculate_deformation(tet, cell_center, rotation_field, ITERATIONS, progress_callback=None, pinned=None):
    native = _native()
    if native is not None:
        if pinned is not None and np.asarray(pinned[0]).any():
            return _deformation_native(native, tet, cell_center, rotation_field, ITERATIONS,
                                       np.asarray(pinned[0], dtype=bool), np.asarray(pinned[1], dtype=float))
        return _deformation_native(native, tet, cell_center, rotation_field, ITERATIONS,
                                   np.zeros(tet.number_of_points, dtype=bool), tet.points)
    if pinned is not None and np.asarray(pinned[0]).any():
        return _calculate_deformation_pinned(tet, cell_center, rotation_field, ITERATIONS, pinned, progress_callback)
    cells = tet.cells
    params = tet.points.copy().flatten()
    rotation_matrices = calculate_rotation_matrices(cell_center, rotation_field)
    old_vertices = tet.points[cells]
    old_transformed = np.einsum("ijk,ikl->ijl", rotation_matrices, (_N @ old_vertices).transpose(0, 2, 1))
    cell_idx = np.repeat(np.arange(tet.number_of_cells), 4)
    vertex_idx = np.ravel(cells)
    jac_rows = np.tile(cell_idx, 3)
    jac_cols = np.concatenate([vertex_idx * 3 + dim for dim in range(3)])
    shape = (tet.number_of_cells, len(params))
    n3 = tet.number_of_points * 3

    def objective_function(p):
        v = p[:n3].reshape(-1, 3)
        transformed = (_N @ v[cells]).transpose(0, 2, 1)
        return np.linalg.norm(transformed - old_transformed, axis=(1, 2)) ** 2

    def objective_jacobian(p):
        v = p[:n3].reshape(-1, 3)
        diff = ((_N @ v[cells]).transpose(0, 2, 1) - old_transformed).transpose(0, 2, 1)
        data = np.concatenate([2 * diff[:, :, dim].ravel() for dim in range(3)]).astype(np.float32)
        return sp.coo_matrix((data, (jac_rows, jac_cols)), shape=shape).tocsr()

    sparsity = sp.coo_matrix((np.ones(len(jac_rows), dtype=np.int8), (jac_rows, jac_cols)), shape=shape).tocsr()
    result = least_squares(objective_function, params, max_nfev=ITERATIONS, verbose=0, jac=objective_jacobian,
                           jac_sparsity=sparsity, method="trf", x_scale="jac",
                           callback=_progress("deformation", ITERATIONS, progress_callback))
    return result.x[:n3].reshape(-1, 3), result.cost


def normal_constraint_rows(original, cells):
    p = original[cells]
    m = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=1)
    d = np.array([[-1.0, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]])
    g_op = np.linalg.solve(m, np.broadcast_to(d, (len(cells), 3, 4)))
    c = p.mean(axis=1)[:, :2]
    r = np.column_stack([c, np.zeros(len(c))])
    r /= np.maximum(np.linalg.norm(r, axis=1, keepdims=True), 1e-9)
    t = np.cross([0.0, 0.0, 1.0], r)
    return np.einsum("tk,tkj->tj", r, g_op), np.einsum("tk,tkj->tj", t, g_op)


def _calculate_deformation_pinned(tet, cell_center, rotation_field, ITERATIONS, pinned, progress_callback=None,
                                  normal_weights=None):
    node_mask, pinned_positions = np.asarray(pinned[0], dtype=bool), np.asarray(pinned[1], dtype=float)
    n_points, cells = tet.number_of_points, tet.cells
    free = np.flatnonzero(~node_mask)
    column = np.full(n_points, -1)
    column[free] = np.arange(len(free))
    base_vertices = tet.points.copy()
    base_vertices[node_mask] = pinned_positions[node_mask]
    rotation_matrices = calculate_rotation_matrices(cell_center, rotation_field)
    old_transformed = np.einsum("ijk,ikl->ijl", rotation_matrices, (_N @ tet.points[cells]).transpose(0, 2, 1))

    cell_rows = np.repeat(np.arange(tet.number_of_cells), cells.shape[1])
    vertex_ids = np.ravel(cells)
    keep = column[vertex_ids] >= 0
    jac_rows = np.tile(cell_rows[keep], 3)
    jac_cols = np.concatenate([column[vertex_ids[keep]] * 3 + dim for dim in range(3)])
    shape = (tet.number_of_cells, len(free) * 3)

    def vertices(x):
        v = base_vertices.copy()
        v[free] = x.reshape(-1, 3)
        return v

    def objective_function(x):
        transformed = (_N @ vertices(x)[cells]).transpose(0, 2, 1)
        return np.linalg.norm(transformed - old_transformed, axis=(1, 2)) ** 2

    def objective_jacobian(x):
        diff = ((_N @ vertices(x)[cells]).transpose(0, 2, 1) - old_transformed).transpose(0, 2, 1)
        data = np.concatenate([2 * diff[:, :, dim].ravel()[keep] for dim in range(3)]).astype(np.float32)
        return sp.coo_matrix((data, (jac_rows, jac_cols)), shape=shape).tocsr()

    sparsity = sp.coo_matrix((np.ones(len(jac_rows), dtype=np.int8), (jac_rows, jac_cols)), shape=shape).tocsr()

    if normal_weights is not None:
        w_out, w_tan = normal_weights
        a_r, a_t = normal_constraint_rows(pinned_positions, cells)
        n_cells = tet.number_of_cells
        z_cols = column[cells] * 3 + 2
        free4 = column[cells] >= 0
        nr_rows = np.repeat(np.arange(n_cells), 4).reshape(-1, 4)
        base_objective, base_jacobian = objective_function, objective_jacobian

        def objective_function(x):  # noqa: F811
            z4 = vertices(x)[cells][:, :, 2]
            sr, st = (a_r * z4).sum(axis=1), (a_t * z4).sum(axis=1)
            return np.concatenate([base_objective(x), w_out * np.minimum(0.0, sr), w_tan * st])

        def objective_jacobian(x):  # noqa: F811
            z4 = vertices(x)[cells][:, :, 2]
            active = ((a_r * z4).sum(axis=1) < 0)[:, None]
            rows_ = np.concatenate([nr_rows[free4] + n_cells, nr_rows[free4] + 2 * n_cells])
            cols_ = np.concatenate([z_cols[free4], z_cols[free4]])
            data = np.concatenate([(w_out * a_r * active)[free4], (w_tan * a_t)[free4]])
            extra = sp.coo_matrix((data, (rows_ - n_cells, cols_)), shape=(2 * n_cells, shape[1]))
            return sp.vstack([base_jacobian(x), extra]).tocsr()

        rows_ = np.concatenate([nr_rows[free4], nr_rows[free4] + n_cells])
        cols_ = np.concatenate([z_cols[free4], z_cols[free4]])
        extra_s = sp.coo_matrix((np.ones(len(rows_), dtype=np.int8), (rows_, cols_)), shape=(2 * n_cells, shape[1]))
        sparsity = sp.vstack([sparsity, extra_s]).tocsr()

    result = least_squares(objective_function, base_vertices[free].ravel(), max_nfev=ITERATIONS, verbose=0,
                           jac=objective_jacobian, jac_sparsity=sparsity, method="trf", x_scale="jac",
                           callback=_progress("deformation", ITERATIONS, progress_callback))
    return vertices(result.x), result.cost


# ---------------------------------------------------------------------------
# Une passe complete
# ---------------------------------------------------------------------------

DEFAULT_PARAMS = {
    "NEIGHBOUR_LOSS_WEIGHT": 20,
    "SET_INITIAL_ROTATION_TO_ZERO": False,
    "INITIAL_ROTATION_FIELD_SMOOTHING": 30,
    "MAX_POS_ROTATION_DEG": 360,
    "MAX_NEG_ROTATION_DEG": -360,
    "ROTATE_ITERATIONS": 100,
    "DEFORM_ITERATIONS": 1000,
    "STEEP_OVERHANG_COMPENSATION": True,
}


def run_one_iteration(tet: BirdTet, topo: BirdTopology, attrs: BirdAttributes, params: dict,
                      progress_callback: Optional[ProgressCallback] = None, planar_base=None):
    """Une passe : rotations optimisees, deformation, attributs recalcules.
    Rend (tet deforme, attributs, info). `planar_base` = (masque des cellules,
    masque des noeuds, positions d'ORIGINE) : base figee a chaque passe."""
    cell_mask = planar_base[0] if planar_base is not None else None
    pinned = (planar_base[1], planar_base[2]) if planar_base is not None else None
    rotation_field, rotate_cost = optimize_rotations(
        tet, topo, attrs, params["NEIGHBOUR_LOSS_WEIGHT"], params["MAX_OVERHANG"], params["ROTATION_MULTIPLIER"],
        params["ROTATE_ITERATIONS"], params["STEEP_OVERHANG_COMPENSATION"],
        params["INITIAL_ROTATION_FIELD_SMOOTHING"], params["SET_INITIAL_ROTATION_TO_ZERO"],
        np.deg2rad(params["MAX_POS_ROTATION_DEG"]), np.deg2rad(params["MAX_NEG_ROTATION_DEG"]),
        progress_callback=progress_callback, fixed_zero_cells=cell_mask)

    if params.get("REACHABLE_NORMAL_WEIGHT"):
        if pinned is None:
            raise ValueError("REACHABLE_NORMAL_WEIGHT exige la base figee (positions d'origine)")
        w = float(params["REACHABLE_NORMAL_WEIGHT"])
        new_vertices, deform_cost = _calculate_deformation_pinned(
            tet, attrs.cell_center, rotation_field, params["DEFORM_ITERATIONS"], pinned, progress_callback,
            normal_weights=(w, 0.5 * w))
    else:
        new_vertices, deform_cost = calculate_deformation(
            tet, attrs.cell_center, rotation_field, params["DEFORM_ITERATIONS"], progress_callback, pinned)

    deformed = tet.with_points(new_vertices)
    new_attrs = update_tet_attributes(deformed, topo)
    with np.errstate(invalid="ignore"):
        overhang_fraction = float(np.nanmean(new_attrs.overhang_angle > np.deg2rad(90 + params["MAX_OVERHANG"])))
    info = {"rotate_cost": float(rotate_cost), "deform_cost": float(deform_cost),
            "overhang_fraction": overhang_fraction, "number_of_cells": int(tet.number_of_cells),
            "rotation_field": rotation_field}
    return deformed, new_attrs, info


def surface_mesh(tet: BirdTet, topo: BirdTopology, recentre: bool = True) -> trimesh.Trimesh:
    """Surface exterieure (pour le juge des replis et l'export). `recentre` :
    centre XY et pose au sol, comme export_surface_stl de la recherche."""
    points = tet.points
    if recentre:
        lo, hi = points.min(axis=0), points.max(axis=0)
        points = points - np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, lo[2]])
    used = np.unique(topo.surface_faces)
    remap = np.full(len(points), -1)
    remap[used] = np.arange(len(used))
    return trimesh.Trimesh(points[used], remap[topo.surface_faces], process=False)
