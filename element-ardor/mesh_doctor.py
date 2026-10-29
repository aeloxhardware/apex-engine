"""
Autodocteur de surface : diagnostic AVANT TetGen, reparation (soudure des aretes courtes,
union des coques qui se traversent, resolution geometrique, pymeshfix) et controle du
maillage tetraedrique APRES. Utilise par `step1_ingest`.

Extrait de l'ancien paquet `pipeline/bird` lors de sa suppression (voie Bird retiree,
sauvegardee dans la branche `bird_old`). Cette partie est du code Apex ; le projet reste
inspire par S4 Slicer de Joshua Bird (voir README, Remerciements).
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Optional

import numpy as np
import trimesh

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

@dataclass
class SurfaceTet:
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

    def with_points(self, points: np.ndarray) -> "SurfaceTet":
        return SurfaceTet(points=np.asarray(points, dtype=float), cells=self.cells)



# Ingestion
# ---------------------------------------------------------------------------

def tetrahedralize(vertices: np.ndarray, faces: np.ndarray, max_tet_volume: Optional[float] = None,
                   keep_boundary: bool = False) -> SurfaceTet:
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
            proc = subprocess.run([sys.executable, "-m", "pipeline._tetgen_worker", src, dst],
                                  cwd=_ROOT, capture_output=True, text=True, timeout=TETGEN_TIMEOUT_S)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"TetGen (autodocteur) depasse {TETGEN_TIMEOUT_S} s (maillage degenere probable)")
        if proc.returncode != 0 or not os.path.exists(dst):
            raise RuntimeError(f"TetGen (autodocteur) a echoue, code {proc.returncode} : {proc.stderr[-2000:]}")
        with np.load(dst) as out:       # fermeture explicite : sous Windows, un .npz ouvert bloque la suppression
            return SurfaceTet(points=np.array(out["node"], dtype=float), cells=np.array(out["elem"], dtype=np.int64))


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


def tet_mesh_problem(tet: SurfaceTet) -> Optional[str]:
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
    qui ampute la piece est rejetee) est applique par l'appelant. Partagee par la autodocteur
    et par l'ingestion FEM (`step1_ingest`) : une seule cascade de reparation. L'ordre est
    celui de la autodocteur d'origine (identique au bit)."""
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


