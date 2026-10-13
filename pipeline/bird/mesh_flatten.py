"""
Aplatit le dessous presque plan d'un STL avant tetraedrisation (voie Bird).

Porte tel quel de research/S4_Slicer/mesh_flatten.py (code Apex, pas Bird) :
un dessous etage de quelques centiemes de mm (tilted benchy undertriangulated)
faisait generer a TetGen ~22 000 tetraedres microscopiques ; apres : 2 740.
La piste « aplatir toutes les zones presque planes » a ete mesuree et ecartee
(TetGen refusait 5 modeles sur 9).
"""
import numpy as np
import trimesh
from scipy.sparse.csgraph import connected_components




BOTTOM_FACING_DEG = 60.0


def flatten_bottom(mesh, layer_height, facing_deg=BOTTOM_FACING_DEG):
    """Ramene sur le plateau (z = z_min) les sommets de la premiere couche,
    SAUF ceux dont l'ecrasement serait faux :
      - une face ENTIEREMENT ecrasee qui ne regardait pas vers le bas (a plus
        de `facing_deg` de -Z) deviendrait une face plate chevauchant le
        dessous -- les bords arrondis qui s'enroulent juste au-dessus du
        plateau (pattes de dino, de Squirtle) : TetGen refusait alors ces
        modeles, sans auto-intersection detectee (chevauchement COPLANAIRE) ;
      - une face qui se RETOURNERAIT (orientation projetee inversee).
    Les sommets de ces faces sont retires des candidats, jusqu'a stabilite.
    Une regle plus simple -- n'ecraser un sommet que si toutes ses faces
    regardent vers le bas -- a ete mesuree trop stricte : les sommets etages
    du dessous du Benchy incline sont sur son bord, ou ils touchent aussi les
    flancs de la coque, et rien n'etait aplati."""
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces)
    z_min = vertices[:, 2].min()
    cos_tol = np.cos(np.radians(facing_deg))
    normals = np.asarray(mesh.face_normals)

    def signed_area_xy(v):
        a, b, c = v[faces[:, 0]], v[faces[:, 1]], v[faces[:, 2]]
        return (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])

    before = signed_area_xy(vertices)
    candidate = vertices[:, 2] < z_min + layer_height
    for _ in range(100):
        trial = vertices.copy()
        trial[candidate, 2] = z_min
        flat = candidate[faces].all(axis=1)
        bad = flat & (normals[:, 2] > -cos_tol)
        touched = candidate[faces].any(axis=1)
        after = signed_area_xy(trial)
        flipped = touched & (np.sign(after) != np.sign(before)) & (np.abs(before) > 1e-12)
        drop = np.unique(faces[bad | flipped].ravel())
        drop = drop[candidate[drop]]
        if not len(drop):
            break
        candidate[drop] = False
    displacement = np.zeros(len(vertices))
    displacement[candidate] = vertices[candidate, 2] - z_min
    if not (displacement > 0).any():
        return mesh, displacement             # rien a aplatir : maillage rendu TEL QUEL
    new = vertices.copy()
    new[candidate, 2] = z_min
    return _rebuild(mesh, new), displacement


def _rebuild(mesh, vertices):
    out = trimesh.Trimesh(vertices, np.asarray(mesh.faces), process=True)
    out.update_faces(out.nondegenerate_faces())
    out.remove_unreferenced_vertices()
    return out
