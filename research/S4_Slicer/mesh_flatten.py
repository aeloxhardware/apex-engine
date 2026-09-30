"""
Nettoyage du maillage a l'entree : aplatit le dessous presque plan.

Pourquoi : un STL dont le dessous est en realite etage de quelques
centiemes de millimetre (tilted benchy undertriangulated : 0 ; 0,004 ;
0,01 ; ... 0,13 mm) fait generer a TetGen des milliers de tetraedres
microscopiques pour remplir ce coin extremement mince -- 22 321 sur ce
modele (48 493 tetraedres en tout). Apres flatten_bottom : 2 740 sur 24 345.
Mesure sur les 9 modeles disponibles (test_mesh_flatten.py) : TetGen les
accepte tous, sans auto-intersection ni perte d'etancheite ; les modeles
deja plats dessous sont inchanges.

PISTE MESUREE ET ECARTEE -- aplatir TOUTES les zones presque planes
(demande de l'utilisateur) : zones de faces voisines d'orientation proche
(< 10 deg), sommets projetes sur le plan ajuste s'ils en sont a moins d'une
demi-couche (0,1 mm), sur l'intersection des plans pour un sommet d'arete.
Resultat : TetGen REFUSAIT ensuite 5 modeles sur 9 (Squirtle, les deux
Benchy, death, dino), avec auto-intersections creees (3 a 4) et etancheite
perdue sur deux d'entre eux, et les tetraedres microscopiques AUGMENTAIENT
ailleurs (B axis mount : 598 -> 1 118). Deplacer des sommets de 0,1 mm
suffit, dans les parties fines, a faire se croiser des faces voisines ou a
creer de nouveaux angles aigus. Par ailleurs, la plupart des modeles ont
deja beaucoup de tetraedres microscopiques SANS zone presque plane
(Squirtle 19 410, death 14 024) : ils viennent des petits details du STL,
pas d'etagements. Code retire (aucun cas d'usage).

Tolerances relatives a la hauteur de couche, aucune constante propre a un
modele. Dependances : numpy, trimesh.
"""
import numpy as np
import trimesh
from scipy.sparse.csgraph import connected_components

ANGLE_TOL_DEG = 10.0


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
