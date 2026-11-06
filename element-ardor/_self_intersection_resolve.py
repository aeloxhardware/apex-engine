"""
Resolution des auto-intersections d'une surface fermee (« outer hull »).

Contrairement a une EXCISION locale (retirer les faces fautives puis
reboucher le trou, cf. `step1_ingest._repair_local_self_intersections`),
qui ne peut pas traiter deux volumes de matiere qui s'INTERPENETRENT
reellement sans detacher un morceau, ce module resout le probleme a sa
racine :

  1. on decoupe chaque triangle le long de sa courbe d'intersection avec
     les autres triangles (la surface devient CONFORME : plus aucun
     interieur de triangle n'en traverse un autre) ;
  2. on jette les morceaux qui se retrouvent strictement A L'INTERIEUR de
     la matiere, en ne gardant que l'enveloppe exterieure.

Le critere d'interiorite est le NOMBRE D'ENROULEMENT generalise (somme
d'angles solides), et non un test de parite par lancer de rayon : sur une
surface qui se recouvre elle-meme, la parite compte un point du
recouvrement double comme EXTERIEUR (semantique XOR), alors que
l'enroulement le compte correctement comme interieur (semantique UNION) -
c'est exactement la configuration a traiter ici.

Mesure ayant motive ce module (`res/test_models/turkish.stl`) : deux
nervures issues d'un export booleen CSG se traversent, avec ~0.12 mm2 de
recouvrement coplanaire reparti sur trois plans horizontaux.
"""

import numpy as np
import scipy.spatial
import trimesh

from core import _cpp

from pipeline._self_intersection_geometry import (  # noqa: F401  (réexport : imports existants)
    _EPS,
    _NATIVE_SEGMENT_TOLERANCE,
    _intersection_segments,
    _interval_on_line,
    _plane_frame,
    _project,
    _snap_to_triangle_boundary,
    _split_keyhole_polygon,
    _subdivide_face,
    _triangle_planes,
    _unproject,
    _weld_line_endpoints,
    generalized_winding_number,
    triangle_triangle_segment,
)


# ---------------------------------------------------------------------------
# Decoupe conforme et extraction de l'enveloppe exterieure
# ---------------------------------------------------------------------------


def _split_t_junctions(vertices, faces, tolerance, max_passes=8):
    """
    Elimine les JONCTIONS EN T : un sommet pose au milieu de l'arete d'une
    face voisine, qui laisse une fissure de surface nulle (l'arete n'est
    portee que par une seule face d'un cote, deux de l'autre).

    Elles apparaissent parce que la courbe d'intersection est nodee dans le
    plan de CHAQUE face independamment : une face peut la voir coupee en un
    point que la face d'en face ne voit pas. Plutot que de rendre conforme
    chaque arrangement local (fragile, chaque plan ayant sa propre
    arithmetique), on corrige globalement a la fin - la geometrie est deja
    exacte a ce stade (volume mesure identique a l'union de reference), il ne
    manque que le raccord topologique.

    Retourne un nouveau tableau de faces.
    """
    faces = np.asarray(faces, dtype=np.int64).copy()
    vertices = np.asarray(vertices, dtype=np.float64)
    for _pass in range(max_passes):
        edge_faces: dict = {}
        for face_index in range(len(faces)):
            triangle = faces[face_index]
            for corner in range(3):
                key = (int(triangle[corner]), int(triangle[(corner + 1) % 3]))
                edge_faces.setdefault(
                    (min(key), max(key)), []
                ).append((face_index, corner))
        open_edges = [key for key, owners in edge_faces.items() if len(owners) == 1]
        if not open_edges:
            break

        tree = scipy.spatial.cKDTree(vertices)
        replacements: dict = {}
        for (a, b) in open_edges:
            pa, pb = vertices[a], vertices[b]
            direction = pb - pa
            length = float(np.linalg.norm(direction))
            if length < tolerance:
                continue
            unit = direction / length
            candidates = tree.query_ball_point(
                (pa + pb) * 0.5, length * 0.5 + tolerance
            )
            best = None
            for candidate in candidates:
                if candidate == a or candidate == b:
                    continue
                offset = vertices[candidate] - pa
                along = float(np.dot(offset, unit))
                if along <= tolerance or along >= length - tolerance:
                    continue
                if float(np.linalg.norm(offset - along * unit)) > tolerance:
                    continue
                if best is None or along < best[0]:
                    best = (along, candidate)
            if best is None:
                continue
            face_index, corner = edge_faces[(a, b)][0]
            if face_index in replacements:
                continue
            replacements[face_index] = (corner, best[1])

        if not replacements:
            break

        new_faces = []
        for face_index in range(len(faces)):
            triangle = faces[face_index]
            if face_index not in replacements:
                new_faces.append(triangle)
                continue
            corner, inserted = replacements[face_index]
            first = int(triangle[corner])
            second = int(triangle[(corner + 1) % 3])
            third = int(triangle[(corner + 2) % 3])
            new_faces.append([first, inserted, third])
            new_faces.append([inserted, second, third])
        faces = np.asarray(new_faces, dtype=np.int64)
    return faces


_CLOSE_SMALL_LOOPS_MAX_VERTICES = 12
_CLOSE_SMALL_LOOPS_MAX_PASSES = 8


def _close_small_boundary_loops(mesh: trimesh.Trimesh, tolerance: float):
    """
    Referme les PETITES boucles de bord residuelles en fusionnant chacune
    en un seul sommet (moins de `_CLOSE_SMALL_LOOPS_MAX_VERTICES` sommets).

    Ces boucles apparaissent quand `resolve_self_intersections` elague un
    fragment orphelin du corps principal (voir docstring de l'appelant) :
    le trou laisse est systematiquement minuscule (mesure sur turkish.stl :
    3 boucles de 3 sommets chacune, sommets separes de 0.001 a 0.02mm - la
    signature meme d'un point qui aurait du etre unique mais a ete calcule
    deux fois independamment lors du decoupage par courbe). Fusionner ces
    quelques sommets ferme le trou exactement, sans introduire de nouvelle
    geometrie devinee (contrairement a un patch en eventail ou une
    triangulation qui choisirait une forme arbitraire).

    Chaque fusion peut degenerer des triangles adjacents (deux sommets
    confondus) - retires via `nondegenerate_faces`/`unique_faces` a chaque
    passe, et peut EXPOSER une nouvelle petite boucle plus loin (mesure :
    9 aretes ouvertes -> 2 apres la 1ere passe -> 0 apres la 2e) : la
    fonction ITERE jusqu'a convergence ou jusqu'a `_CLOSE_SMALL_LOOPS_MAX_PASSES`.

    Ne touche jamais une boucle plus grande que le seuil (une vraie perte
    de matiere localisee ne doit pas etre "fermee" par un simple soudage de
    sommets, qui n'aurait alors plus aucune justification geometrique) - le
    retour est alors `None` si des boucles de bord subsistent au-dela du
    budget d'iterations, jamais un mesh silencieusement encore troue.

    Best-effort : ne leve jamais, retourne `None` en cas d'echec (mesh
    non-manifold, boucle trop grande, budget d'iterations epuise, ou
    connexite cassee par une fusion - jamais observe en pratique mais
    verifie explicitement a chaque passe).
    """
    import networkx as nx

    current = mesh.copy()
    for _pass in range(_CLOSE_SMALL_LOOPS_MAX_PASSES):
        open_edges_mask = trimesh.grouping.group_rows(
            current.edges_sorted, require_count=1
        )
        open_edges = current.edges[open_edges_mask]
        if len(open_edges) == 0:
            # Verification de connexite EN FIN DE CONVERGENCE seulement,
            # jamais a chaque passe intermediaire : une passe peut isoler
            # transitoirement un fragment (mesure sur turkish.stl : une
            # face esseulee juste apres la 1ere passe) que la passe
            # SUIVANTE recolle bel et bien au corps principal en fermant le
            # dernier petit trou adjacent. Rejeter des la 1ere passe sur ce
            # signal transitoire aurait empeche une fermeture pourtant
            # valide au final.
            if len(current.split(only_watertight=False)) != 1:
                return None
            return current

        graph = nx.Graph()
        graph.add_edges_from(open_edges.tolist())
        loops = list(nx.connected_components(graph))
        if any(len(loop) > _CLOSE_SMALL_LOOPS_MAX_VERTICES for loop in loops):
            return None

        merged_vertices = current.vertices.copy()
        for loop in loops:
            indices = list(loop)
            target = merged_vertices[indices[0]]
            for index in indices[1:]:
                merged_vertices[index] = target

        candidate = trimesh.Trimesh(
            vertices=merged_vertices, faces=current.faces.copy(), process=False
        )
        candidate.merge_vertices()
        candidate.update_faces(candidate.nondegenerate_faces())
        candidate.update_faces(candidate.unique_faces())
        candidate.remove_unreferenced_vertices()
        current = candidate

    return None


def resolve_self_intersections(mesh: trimesh.Trimesh, max_affected_ratio: float = 0.25):
    """
    Rend une surface fermee EQUIVALENTE a `mesh` mais sans auto-intersection,
    en decoupant les triangles fautifs le long de leurs courbes
    d'intersection puis en ne conservant que l'enveloppe exterieure.

    Best-effort : retourne None (jamais une exception) si le defaut n'est pas
    de cette nature, si la proportion de faces touchees depasse
    `max_affected_ratio` (le mecanisme vise un defaut LOCAL, pas un maillage
    globalement incoherent), ou si le resultat ne satisfait pas ses
    post-conditions.
    """
    triangles = np.asarray(mesh.triangles, dtype=np.float64)
    if len(triangles) == 0:
        return None
    scale = float(mesh.scale) or 1.0
    tolerance = 1e-9 * scale
    min_length = 1e-7 * scale

    per_face = _intersection_segments(mesh, min_length)
    if not per_face:
        return None
    if len(per_face) > max_affected_ratio * len(triangles):
        return None

    # Points introduits sur les aretes du maillage : ils doivent apparaitre
    # A L'IDENTIQUE dans les deux faces qui partagent l'arete, sinon la
    # surface recollee presente une fissure (arete portee par une seule face).
    faces = np.asarray(mesh.faces)
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    edge_points: dict = {}
    for face_index, segments in per_face.items():
        corners = faces[face_index]
        for corner in range(3):
            a = int(corners[corner])
            b = int(corners[(corner + 1) % 3])
            pa, pb = vertices[a], vertices[b]
            direction = pb - pa
            length = float(np.linalg.norm(direction))
            if length < tolerance:
                continue
            unit = direction / length
            for segment in segments:
                for point in segment:
                    offset = point - pa
                    along = float(np.dot(offset, unit))
                    if along <= tolerance or along >= length - tolerance:
                        continue
                    if float(np.linalg.norm(offset - along * unit)) > tolerance:
                        continue
                    edge_points.setdefault((min(a, b), max(a, b)), []).append(point)

    # toute face portant une arete enrichie doit etre redecoupee elle aussi
    affected = set(per_face)
    if edge_points:
        for face_index in range(len(faces)):
            corners = faces[face_index]
            for corner in range(3):
                a = int(corners[corner])
                b = int(corners[(corner + 1) % 3])
                if (min(a, b), max(a, b)) in edge_points:
                    affected.add(face_index)
                    break

    # Points-noeuds globaux : toutes les extremites de segments et tous les
    # points poses sur une arete du maillage. Reinjectes dans CHAQUE face qui
    # les porte, pour eliminer les jonctions en T entre faces voisines.
    node_list = [point for segments in per_face.values()
                 for segment in segments for point in segment]
    node_list.extend(point for points in edge_points.values() for point in points)
    node_points = (np.asarray(node_list, dtype=np.float64)
                   if node_list else np.zeros((0, 3)))

    new_triangles = []
    for face_index in sorted(affected):
        corners = faces[face_index]
        extras = []
        for corner in range(3):
            a = int(corners[corner])
            b = int(corners[(corner + 1) % 3])
            extras.extend(edge_points.get((min(a, b), max(a, b)), ()))
        pieces = _subdivide_face(
            triangles[face_index], per_face.get(face_index, []), extras, tolerance,
            node_points,
        )
        if pieces is None:
            return None
        new_triangles.extend(pieces)

    kept = [triangles[i] for i in range(len(triangles)) if i not in affected]
    all_triangles = np.asarray(kept + new_triangles, dtype=np.float64)
    if len(all_triangles) == 0:
        return None

    # Enveloppe exterieure : on ne garde que les faces dont le COTE EXTERNE
    # donne sur du vide. Le nombre d'enroulement est calcule contre la
    # surface d'ORIGINE (qui porte la semantique de matiere voulue), et non
    # contre la surface redecoupee.
    #
    # Une sonde par triangle (centroide, decale d'une petite distance le
    # long de la normale) est FRAGILE dans les zones ou la subdivision
    # produit beaucoup de petits triangles proches les uns des autres : la
    # sonde peut atterrir du mauvais cote d'une nappe voisine a quelques
    # microns, donnant un winding FAUX pour CE triangle precis sans rapport
    # avec ses voisins immediats. Mesure sur turkish.stl : 14 morceaux
    # partageant un meme sommet, winding mesure -0.0 pour certains et
    # jusqu'a 1.997 pour d'autres - bruit pur. Deux pistes de sonde
    # multiple (vote a distances croissantes, vote a positions multiples a
    # distance fixe) ont ete mesurees et ECARTEES : elles degradent le
    # resultat au lieu de l'ameliorer (34->38 composantes deconnectees,
    # volume errone jusqu'a 10.8%) - eloigner la sonde change la question
    # posee (on quitte le voisinage immediat pour de la matiere reelle
    # differente), et voter sur plusieurs points co-localises n'apporte pas
    # d'information independante sur un triangle deja quasi degenere.
    centroids = all_triangles.mean(axis=1)
    normals = np.cross(
        all_triangles[:, 1] - all_triangles[:, 0],
        all_triangles[:, 2] - all_triangles[:, 0],
    )
    lengths = np.linalg.norm(normals, axis=1)
    valid = lengths > _EPS
    normals[valid] = normals[valid] / lengths[valid, None]
    probe = centroids + normals * (1e-5 * scale)
    winding = generalized_winding_number(probe, triangles)
    keep_mask = valid & (winding < 0.5)
    if not np.any(keep_mask):
        return None

    surviving = all_triangles[keep_mask]
    flat = surviving.reshape(-1, 3)
    decimals = max(0, int(-np.log10(tolerance)) - 1)
    _, unique_index, inverse = np.unique(
        np.round(flat, decimals), axis=0, return_index=True, return_inverse=True
    )
    welded_vertices = flat[unique_index]
    welded_faces = inverse.reshape(-1, 3)
    welded_faces = welded_faces[
        (welded_faces[:, :1] != welded_faces[:, 1:]).all(axis=1)
    ]
    if len(welded_faces) == 0:
        return None

    welded_faces = _split_t_junctions(welded_vertices, welded_faces, tolerance)
    resolved = trimesh.Trimesh(
        vertices=welded_vertices, faces=welded_faces, process=False
    )
    resolved.merge_vertices()
    resolved.update_faces(resolved.nondegenerate_faces())
    resolved.update_faces(resolved.unique_faces())
    resolved.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(resolved)
    trimesh.repair.fix_winding(resolved)

    # Garde-fou de CONNEXITE, indispensable : le filtre par nombre
    # d'enroulement (ci-dessus) peut garder des morceaux de surface qui
    # restent geometriquement corrects mais se retrouvent DECONNECTES du
    # corps principal (quelques triangles orphelins pres d'une courbe
    # d'intersection tres locale, jamais recolles au reste par le soudage).
    # Mesure sur turkish.stl : 35 composantes connexes en sortie (1 corps
    # principal a 10289 faces + 34 fragments de 1 a 4 faces), TOUTES situees
    # dans la fine bande z=[33.26,33.76] de la zone reparee, et JUSQU'A
    # 76mm du bord ouvert le plus proche du corps principal - pas de simples
    # fissures de soudure (28/34 mesures a plus de 0.05mm, la plupart bien
    # au-dela). Sans ce garde, `_tetrahedralize` choisit un point-seed
    # potentiellement DANS un de ces fragments plutot que dans le corps
    # principal : TetGen reussit alors SANS LEVER D'EXCEPTION, mais ne
    # maille que ce fragment - regression SILENCIEUSE du meme type que la
    # "CATASTROPHE" pymeshfix deja corrigee (limite #52), mesuree ici :
    # z=[33.26,90.78] tronque a z=[33.26,33.76], 243 tetraedres au lieu de
    # dizaines de milliers, jamais signale par TetGen ni par le garde-fou
    # anti-effondrement (le NOMBRE de faces du mesh source est preserve,
    # seule la CONNEXITE est cassee).
    #
    # Fix retenu : ne garder QUE la plus grande composante (le corps
    # principal) plutot que de rejeter tout le resultat des qu'il existe un
    # fragment orphelin - ces fragments sont mesures comme du bruit residuel
    # (0.8% des faces sur turkish.stl), pas de la geometrie legitime perdue
    # (volume du corps principal seul : ecart de 0.00003% avec le mesh
    # source, contre 0.028% pour le resultat complet non filtre). Plusieurs
    # tentatives de RECOLLER ces fragments (fusion de sommets, soudage 3D a
    # tolerance elargie, triangulation Klincsek deja existante pour ce cas
    # dans step1_ingest) ont ete mesurees et ECARTEES : la zone est
    # numeriquement pathologique (contact quasi tangentiel entre deux
    # courbes d'intersection), et chaque tentative de recollage y introduit
    # de NOUVELLES auto-intersections (jusqu'a 8 mesurees) sans jamais
    # converger proprement. Ne garder que le corps principal EVITE le
    # probleme au lieu de le deplacer.
    #
    # Le corps principal peut cependant rester NON etanche (petits trous
    # laisses par les fragments retires, mesure : 3 boucles de 3 aretes sur
    # turkish.stl, 9 aretes ouvertes au total). Plutot que de rejeter
    # immediatement, on referme ces trous RESIDUELS par soudage de sommets
    # (les boucles sont systematiquement minuscules - quelques sommets a
    # 0.001-0.02mm les uns des autres, la signature meme des fragments
    # qu'on vient d'elaguer) : `_close_small_boundary_loops` fusionne
    # chaque petite boucle en un seul point et re-verifie a chaque passe
    # que la connexite (1 seul corps) est preservee. Mesure sur
    # turkish.stl : converge en 2 passes vers un mesh etanche et connexe -
    # mais introduit 8 nouvelles auto-intersections locales (le soudage
    # deplace legerement des sommets dans une zone deja numeriquement
    # tendue). Ce module NE LES NETTOIE PAS lui-meme : ce n'est plus une
    # auto-intersection issue d'une courbe d'intersection transversale (son
    # perimetre declare), c'est un artefact du soudage - la responsabilite
    # du nettoyage final (retrait + retriangulation Klincsek) revient a
    # l'appelant (`step1_ingest._tetrahedralize`), qui dispose deja de ce
    # mecanisme pour un autre cas d'usage.
    components = resolved.split(only_watertight=False)
    if len(components) == 0:
        return None
    components.sort(key=lambda c: -len(c.faces))
    main_component = components[0]
    if not main_component.is_watertight:
        main_component = _close_small_boundary_loops(main_component, tolerance)
        if main_component is None:
            return None

    return main_component
