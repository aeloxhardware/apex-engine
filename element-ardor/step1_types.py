"""Etape 1 : rapport d'ingestion, erreurs et contrôle d'emprise sur le plateau.

Découpé de `step1_ingest.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step1_ingest`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import numpy as np
import trimesh


@dataclass
class IngestReport:
    """
    Resume de l'ingestion, remonte a l'UI TypeScript pour donner un retour
    sur la qualite du mesh d'entree avant de lancer le reste du pipeline.
    """
    was_watertight_initially: bool
    is_watertight_after_repair: bool
    faces_removed_degenerate: int
    vertices_merged: int
    tetrahedra_generated: int
    source_face_count: int
    source_vertex_count: int
    # Auto-intersections : donnees EXPLICITES plutot qu'un silence, meme
    # discipline que CollisionReport.chuck_checked (etape 9). Un maillage
    # peut etre watertight, a winding coherent, sans face degeneree, et
    # rester refuse par TetGen parce que des triangles se traversent -
    # `is_watertight_after_repair` seul ne dit rien de ce cas.
    self_intersecting_faces_before: int = 0
    self_intersection_repair_attempted: bool = False
    self_intersection_repair_succeeded: bool = False
    # Niveau de plateau REELLEMENT utilise par cette ingestion - toujours
    # renseigne (jamais None), meme quand `run(platform_z=None)` l'a fait
    # deriver automatiquement du mesh source (min(z) des sommets nettoyes).
    # Donnee explicite plutot qu'un silence, meme discipline que
    # `self_intersecting_faces_before`/`CollisionReport.chuck_checked` :
    # un modele qui n'est pas pose a z=0 (ex: res/test_models/death.stl,
    # min(z)=1.245) casse silencieusement tout ancrage plateau en aval
    # (step3 puis step4) si aucune etape ne remonte la valeur reellement
    # utilisee - voir CLAUDE.md.
    platform_z_used: float = 0.0
    # Translation en z appliquee a l'ingestion par `place_on_platform` (mm),
    # 0.0 si le modele n'a pas ete deplace. Donnee explicite : le maillage
    # remonte par `run` n'est plus aux coordonnees du fichier source.
    placement_offset_z: float = 0.0
    # Fragments negligeables (coques individuellement watertight mais
    # minuscules par rapport a la coque dominante) ecartes pour rendre le
    # mesh watertight - donnee EXPLICITE, jamais une perte silencieuse de
    # geometrie. 0 dans l'immense majorite des cas (mesh deja propre) -
    # mesure sur res/test_models/turkish.stl : 16 fragments (76 faces sur
    # 10180, 0.31% du volume) accoles au corps principal par des jonctions
    # non-manifold, faisant echouer la reparation standard ET pymeshfix
    # (qui reduisait le mesh entier a 10 faces - voir CLAUDE.md).
    dominant_shell_fragments_dropped: int = 0
    # Ecart entre la SURFACE du modele et la frontiere du volume tetraedrise
    # (core/tet_coverage). `None` = non mesure. Un maillage dont la frontiere
    # s'ecarte de la peau laisse des trous dans les couches que step5 ne voit
    # pas : ses contours sont les sections des tetraedres. Donnee explicite.
    surface_outside_volume_fraction: Optional[float] = None   # part de la peau a plus de `tolerance`
    surface_outside_volume_area: float = 0.0                  # mm2
    surface_max_distance_to_tets: float = 0.0                 # mm, pire ecart
    tet_volume_ratio: float = 0.0                             # volume tetra / volume du modele
    surface_outside_volume_clusters: list = field(default_factory=list)  # (aire mm2, centre) par amas
    # Reparation SYSTEMATIQUE de la surface (2026-10-22) : defauts detectes par
    # l'autodocteur (`pipeline.mesh_doctor.diagnose_surface`) et libelle de la
    # reparation retenue (None = surface saine ou acceptee telle quelle par TetGen).
    # Donnee explicite : un modele repare n'est plus tout a fait le fichier source.
    surface_problems: list = field(default_factory=list)
    surface_repair: Optional[str] = None
    # Mailleur utilise ("tetgen" ou "tetranet") et, pour TetraNet, son diagnostic chiffre de la
    # surface d'entree (trous, faces retournees, points incoherents...) : donnee explicite.
    mesher: str = "tetgen"
    mesher_report: Optional[dict] = None


class MeshIngestError(Exception):
    """Levee quand le mesh source ne peut pas etre charge ou exploite."""


class SurfaceInconsistentError(MeshIngestError):
    """Surface trop abimee pour distinguer l'interieur de l'exterieur (R-016) : `str(exc)` est le
    message de l'utilisateur, `exc.diagnostic` le diagnostic chiffre de TetraNet."""

    def __init__(self, message: str, diagnostic: dict):
        super().__init__(message)
        self.diagnostic = diagnostic


class ModelTooLargeError(MeshIngestError):
    """Le modele ne peut tenir sur le plateau DANS AUCUNE position (CLAUDE.md,
    feuille de route « refuser dès l'entrée un modèle trop grand »).
    `footprint_radius` et `platform_radius` en mm, pour l'affichage UI."""

    def __init__(self, footprint_radius: float, platform_radius: float, height: float | None = None,
                 max_height: float | None = None):
        self.footprint_radius = footprint_radius
        self.platform_radius = platform_radius
        self.height = height
        self.max_height = max_height
        parts = []
        if footprint_radius > platform_radius:
            parts.append(f"empreinte au sol de rayon {footprint_radius:.1f} mm > plateau de rayon "
                         f"{platform_radius:.1f} mm")
        if height is not None and max_height is not None and height > max_height:
            parts.append(f"hauteur {height:.1f} mm > hauteur utile {max_height:.1f} mm")
        super().__init__("modele trop grand pour la machine : " + " ; ".join(parts))


def _min_enclosing_circle(points: np.ndarray) -> tuple[np.ndarray, float]:
    """Plus petit cercle englobant de points 2D (Welzl iteratif, ordre
    aleatoire a graine fixe : resultat deterministe). Exact a l'arrondi pres."""
    pts = np.unique(np.asarray(points, dtype=float), axis=0)
    if len(pts) == 0:
        return np.zeros(2), 0.0
    rng = np.random.default_rng(0)
    pts = pts[rng.permutation(len(pts))]

    def circle2(a, b):
        c = (a + b) / 2.0
        return c, float(np.linalg.norm(a - c))

    def circle3(a, b, c):
        d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
        if abs(d) < 1e-15:  # colineaires : le cercle du plus grand segment
            best = max((circle2(a, b), circle2(a, c), circle2(b, c)), key=lambda t: t[1])
            return best
        ux = ((a @ a) * (b[1] - c[1]) + (b @ b) * (c[1] - a[1]) + (c @ c) * (a[1] - b[1])) / d
        uy = ((a @ a) * (c[0] - b[0]) + (b @ b) * (a[0] - c[0]) + (c @ c) * (b[0] - a[0])) / d
        center = np.array([ux, uy])
        return center, float(np.linalg.norm(a - center))

    def inside(c, r, p):
        return np.linalg.norm(p - c) <= r * (1.0 + 1e-12) + 1e-12

    c, r = pts[0].copy(), 0.0
    for i in range(1, len(pts)):
        if inside(c, r, pts[i]):
            continue
        c, r = pts[i].copy(), 0.0
        for j in range(i):
            if inside(c, r, pts[j]):
                continue
            c, r = circle2(pts[i], pts[j])
            for k in range(j):
                if not inside(c, r, pts[k]):
                    c, r = circle3(pts[i], pts[j], pts[k])
    return c, r


def check_fits_platform(vertices: np.ndarray, platform_radius: float | None,
                        max_height: float | None = None) -> float:
    """Leve ModelTooLargeError si le modele ne tient sur le plateau dans aucune
    position : son empreinte XY (plus petit cercle englobant, calcule sur
    l'enveloppe convexe) depasse le rayon du plateau, ou sa hauteur depasse la
    hauteur utile. Independant du PLACEMENT (decide par l'interface) : ne
    refuse que l'impossible. Rend le rayon d'empreinte (mm)."""
    v = np.asarray(vertices, dtype=float)
    if platform_radius is None or len(v) == 0:
        return 0.0
    xy = v[:, :2]
    try:
        from scipy.spatial import ConvexHull
        xy = xy[ConvexHull(xy).vertices]
    except Exception:  # points colineaires/degeneres : cercle sur tous les points
        pass
    _, radius = _min_enclosing_circle(xy)
    height = float(np.ptp(v[:, 2]))
    if radius > platform_radius or (max_height is not None and height > max_height):
        raise ModelTooLargeError(radius, float(platform_radius), height, max_height)
    return radius


def _load_raw_mesh(mesh_path: str, allow_multi_mesh: bool) -> trimesh.Trimesh:
    loaded = trimesh.load(mesh_path, force="mesh" if not allow_multi_mesh else None)

    if isinstance(loaded, trimesh.Scene):
        if not allow_multi_mesh:
            raise MeshIngestError(
                f"'{mesh_path}' contient plusieurs objets (Scene) mais "
                "allow_multi_mesh=False. Fusionnez le fichier ou activez ce parametre."
            )
        if len(loaded.geometry) == 0:
            raise MeshIngestError(f"'{mesh_path}' ne contient aucune geometrie.")
        # `to_geometry()` applique les transformations de nœud de la scene (placement des objets d'un 3MF) ;
        # concatener `geometry.values()` les perdait et empilait les pieces a l'origine (R-013).
        loaded = loaded.to_geometry()

    if not isinstance(loaded, trimesh.Trimesh):
        raise MeshIngestError(
            f"'{mesh_path}' n'a pas produit un Trimesh exploitable "
            f"(type obtenu: {type(loaded)})."
        )

    if loaded.faces.shape[0] == 0:
        raise MeshIngestError(f"'{mesh_path}' est charge mais ne contient aucune face.")

    return loaded
