"""Étape 1 avec TetraNet (mailleur natif) : volume juste, défauts tolérés, surface trop abîmée refusée."""
import numpy as np
import pytest
import trimesh

from pipeline import step1_ingest

pytestmark = pytest.mark.skipif(not step1_ingest.tetranet_available(), reason="noyau natif sans TetraNet")


def _write(tmp_path, mesh, name="m.stl"):
    path = str(tmp_path / name)
    mesh.export(path)
    return path


def _tet_volume(tm):
    p = tm.tet_vertices[tm.tetrahedra]
    return float(np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]).sum() / 6)


def test_clean_box_volume_and_report(tmp_path):
    box = trimesh.creation.box(extents=(20, 12, 8))
    md, tm, rep = step1_ingest.run(_write(tmp_path, box), resolution=1.5, place_on_platform=True, mesher="tetranet")
    assert rep.mesher == "tetranet" and rep.mesher_report["ambiguous_points"] == 0
    assert abs(_tet_volume(tm) / box.volume - 1.0) < 0.03
    assert md.is_watertight and abs(rep.platform_z_used) < 1e-9          # posé sur le plateau
    assert tm.tetrahedra.min() >= 0 and tm.tetrahedra.max() < tm.tet_vertices.shape[0]


def test_a_missing_face_is_tolerated(tmp_path):
    box = trimesh.creation.box(extents=(20, 12, 8))
    holed = trimesh.Trimesh(box.vertices, box.faces[:-1], process=False)  # un triangle en moins
    md, tm, rep = step1_ingest.run(_write(tmp_path, holed), resolution=1.5, mesher="tetranet")
    assert abs(_tet_volume(tm) / box.volume - 1.0) < 0.05
    assert rep.mesher_report["surface_holes"] >= 1                        # le défaut est compté, pas caché


def test_flipped_faces_are_reoriented(tmp_path):
    sphere = trimesh.creation.icosphere(subdivisions=3, radius=8.0)
    faces = sphere.faces.copy()
    faces[:200] = faces[:200][:, ::-1]                                    # une calotte retournée
    broken = trimesh.Trimesh(sphere.vertices, faces, process=False)
    md, tm, rep = step1_ingest.run(_write(tmp_path, broken), resolution=1.5, mesher="tetranet")
    assert abs(_tet_volume(tm) / sphere.volume - 1.0) < 0.05
    assert rep.mesher_report["surface_reoriented_faces"] > 0


def test_very_broken_surface_is_refused_with_the_agreed_message(tmp_path):
    sphere = trimesh.creation.icosphere(subdivisions=3, radius=8.0)
    keep = sphere.triangles_center[:, 2] < 0                              # une demi-sphère : surface ouverte
    shell = trimesh.Trimesh(sphere.vertices, sphere.faces[keep], process=False)
    with pytest.raises(step1_ingest.SurfaceInconsistentError) as exc:
        step1_ingest.run(_write(tmp_path, shell), resolution=1.5, mesher="tetranet", tetranet_refuse_ratio=0.0001)
    message = str(exc.value)
    assert message.startswith("Erreur de tétraédralisation : la surface du modèle présente un trou,")
    assert "incohérent" in message and exc.value.diagnostic["surface_holes"] == 1


def test_tetgen_stays_the_default_path(tmp_path):
    box = trimesh.creation.box(extents=(10, 10, 10))
    md, tm, rep = step1_ingest.run(_write(tmp_path, box), resolution=2.0)
    assert rep.mesher == "tetgen"


def test_main_exits_with_code_3_and_the_message_on_a_refused_model(tmp_path, capsys):
    import main

    sphere = trimesh.creation.icosphere(subdivisions=3, radius=8.0)
    shell = trimesh.Trimesh(sphere.vertices, sphere.faces[sphere.triangles_center[:, 2] < 0], process=False)
    # seuil de refus par defaut (5 %) : une demi-sphere entiere ouverte le depasse
    code = main.main([_write(tmp_path, shell), "--mesher", "tetranet"])
    err = capsys.readouterr().err
    assert code == 3
    assert "Erreur de tétraédralisation : la surface du modèle présente un trou," in err
