"""Etape 1 : structures d'accélération du maillage et chemin TetraNet.

Découpé de `step1_ingest.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step1_ingest`.
"""

from __future__ import annotations

from typing import Optional
import numpy as np
import trimesh
from core import _cpp
from core.mesh_types import MeshData, TetMesh
from core.spatial import KDTree, RayIntersector
from pipeline.step1_repair import (  # noqa: F401
    _count_self_intersecting_faces,
)
from pipeline.step1_tetgen_run import (  # noqa: F401
    _build_tet_adjacency,
)
from pipeline.step1_types import (  # noqa: F401
    IngestReport,
    SurfaceInconsistentError,
)


def _build_acceleration_structures(mesh: trimesh.Trimesh):
    # core.spatial : trimesh/scipy ou apex_core selon APEX_USE_CPP, meme
    # interface dans les deux cas (voir core/spatial.py). Le bvh est
    # interroge par step2/step7, le kdtree par step4/5/6 - tous deux partages
    # via MeshData.
    # Cote trimesh, RayMeshIntersector choisit pyembree s'il est installe et
    # retombe sinon sur l'implementation pure numpy ; pyembree est ABSENT de
    # ce venv (verifie), d'ou le gain reel du chemin natif.
    bvh = RayIntersector(mesh)
    kdtree = KDTree(mesh.triangles_center)
    return bvh, kdtree


def _build_mesh_data(mesh: trimesh.Trimesh) -> MeshData:
    bvh, kdtree = _build_acceleration_structures(mesh)
    return MeshData(
        vertices=np.asarray(mesh.vertices, dtype=np.float64),
        faces=np.asarray(mesh.faces, dtype=np.int64),
        face_normals=np.asarray(mesh.face_normals, dtype=np.float64),
        vertex_normals=np.asarray(mesh.vertex_normals, dtype=np.float64),
        is_watertight=bool(mesh.is_watertight),
        bvh=bvh,
        kdtree=kdtree,
    )


# Arete du reseau de TetraNet par unite de `resolution` : TetGen borne le VOLUME d'un tetraedre
# a resolution^3 ; TetraNet fait des tetraedres de volume cell^3 / 12. Le facteur 2/3 donne a
# resolution 1,5 la maille de 1,0 mm validee sur le Benchy (tetranet/docs, CLAUDE.md).
TETRANET_CELL_PER_RESOLUTION = 2.0 / 3.0


TETRANET_STATUS_INCONSISTENT_SURFACE = 4


def tetranet_available() -> bool:
    """Le noyau natif est charge ET contient TetraNet (jamais une erreur s'il manque)."""
    if not _cpp.enabled():
        return False
    module = _cpp.module()
    return module is not None and hasattr(module, "tetranet_mesh")


def _ingest_with_tetranet(
    raw_mesh: trimesh.Trimesh, resolution: float, platform_z: Optional[float], place_on_platform: bool,
    check_surface_coverage: bool, refuse_ratio: float,
) -> tuple[MeshData, TetMesh, IngestReport]:
    """Etape 1 SANS reparation de surface : TetraNet mesure, soude, remet les faces dans le bon sens
    et tetraedrise en un appel natif. Aucune reparation (pymeshfix, union booleenne, TetGen comme
    juge) : le nombre d'enroulement tolere trous et faces retournees. Une surface trop abimee leve
    `SurfaceInconsistentError` avec le message de l'utilisateur (R-016)."""
    vertices = np.ascontiguousarray(raw_mesh.vertices, dtype=np.float64)
    faces = np.ascontiguousarray(raw_mesh.faces, dtype=np.int64)
    placement_offset = 0.0
    if place_on_platform and len(vertices):
        target_z = 0.0 if platform_z is None else float(platform_z)
        placement_offset = target_z - float(vertices[:, 2].min())
        if placement_offset != 0.0:
            vertices = vertices + np.array([0.0, 0.0, placement_offset])
    out = _cpp.module().tetranet_mesh(vertices, faces, float(resolution) * TETRANET_CELL_PER_RESOLUTION,
                                      0.5, 0.05, float(refuse_ratio), 0)
    mesher_report = {k: v for k, v in out.items()
                     if k not in ("tet_vertices", "tetrahedra", "surface_vertices", "surface_faces")}
    if out["tets"] == 0 or out["status"] == TETRANET_STATUS_INCONSISTENT_SURFACE:
        from core.surface_message import surface_error_message
        cleaned = trimesh.Trimesh(out["surface_vertices"], out["surface_faces"], process=False) \
            if len(out["surface_faces"]) else raw_mesh
        self_intersections = int(_count_self_intersecting_faces(cleaned))
        volume_pct = 100.0 * float(out["volume_ratio"])
        # Les deux nombres dont l'interface a besoin pour reconstruire le texte (traduisible) : R-016.
        diagnostic = {**mesher_report, "self_intersections": self_intersections, "volume_pct": volume_pct}
        raise SurfaceInconsistentError(surface_error_message(mesher_report, volume_pct, self_intersections),
                                       diagnostic)
    mesh = trimesh.Trimesh(out["surface_vertices"], out["surface_faces"], process=False)
    tet_mesh = TetMesh(tet_vertices=np.asarray(out["tet_vertices"], dtype=np.float64),
                       tetrahedra=np.asarray(out["tetrahedra"], dtype=np.int64))
    tet_mesh.tet_adjacency = _build_tet_adjacency(tet_mesh.tetrahedra, tet_mesh.tet_vertices.shape[0])
    resolved_platform_z = float(mesh.vertices[:, 2].min()) if platform_z is None else float(platform_z)
    report = IngestReport(
        was_watertight_initially=bool(out["surface_open_edges"] == 0 and out["surface_nonmanifold_edges"] == 0
                                      and out["surface_flipped_edges"] == 0),
        is_watertight_after_repair=bool(mesh.is_watertight),
        faces_removed_degenerate=int(out["surface_degenerate_faces"]),
        vertices_merged=max(0, int(raw_mesh.vertices.shape[0]) - int(mesh.vertices.shape[0])),
        tetrahedra_generated=int(tet_mesh.tet_count),
        source_face_count=int(mesh.faces.shape[0]), source_vertex_count=int(mesh.vertices.shape[0]),
        platform_z_used=resolved_platform_z, placement_offset_z=placement_offset,
        mesher="tetranet", mesher_report=mesher_report,
    )
    mesh_data = _build_mesh_data(mesh)
    if check_surface_coverage:
        from core.tet_coverage import surface_outside_volume
        coverage = surface_outside_volume(mesh.vertices, mesh.faces, tet_mesh.tet_vertices, tet_mesh.tetrahedra)
        report.surface_outside_volume_fraction = coverage.far_fraction
        report.surface_outside_volume_area = coverage.far_area
        report.surface_max_distance_to_tets = coverage.max_distance
        report.tet_volume_ratio = coverage.volume_ratio
        report.surface_outside_volume_clusters = [(c.area, c.center) for c in coverage.clusters[:10]]
    return mesh_data, tet_mesh, report
