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

from element_ardor._self_intersection_geometry import (  # noqa: F401  (réexport : imports existants)
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

    faces = np.asarray(mesh.faces)
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    edge_points = _edge_points_of(per_face, faces, vertices, tolerance)
    affected = _affected_faces(per_face, faces, edge_points)
    new_triangles = _subdivide_affected(triangles, faces, per_face, edge_points, affected, tolerance)
    if new_triangles is None:
        return None

    kept = [triangles[i] for i in range(len(triangles)) if i not in affected]
    all_triangles = np.asarray(kept + new_triangles, dtype=np.float64)
    if len(all_triangles) == 0:
        return None
    surviving = _exterior_envelope(all_triangles, triangles, scale)
    if surviving is None:
        return None
    resolved = _weld_surface(surviving, tolerance)
    if resolved is None:
        return None
    return _main_component(resolved, tolerance)


def _edge_points_of(per_face: dict, faces, vertices, tolerance: float) -> dict:
    """Points introduits sur les aretes du maillage : ils doivent apparaitre A L'IDENTIQUE dans les deux faces
    qui partagent l'arete, sinon la surface recollee presente une fissure (arete portee par une seule face)."""
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
    return edge_points


def _affected_faces(per_face: dict, faces, edge_points: dict) -> set:
    """Toute face portant une arete enrichie doit etre redecoupee elle aussi."""
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
    return affected


def _subdivide_affected(triangles, faces, per_face: dict, edge_points: dict, affected: set, tolerance: float):
    """Subdivise les faces touchees ; None si l'une d'elles echoue. Les points-noeuds globaux (extremites de
    segments et points d'arete) sont reinjectes dans CHAQUE face qui les porte (pas de jonctions en T)."""
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
    return new_triangles


def _exterior_envelope(all_triangles, triangles, scale: float):
    """Enveloppe exterieure : ne garde que les faces dont le COTE EXTERNE donne sur du vide. Le nombre
    d'enroulement est calcule contre la surface d'ORIGINE (qui porte la semantique de matiere voulue). Une sonde
    par triangle (centroide decale le long de la normale) est fragile dans les zones de petits triangles
    serres (turkish.stl : winding -0.0 a 1.997 pour 14 morceaux d'un meme sommet) ; deux sondes multiples
    (vote a distances croissantes, a positions multiples) ont ete mesurees et ECARTEES (34->38 composantes,
    volume errone jusqu'a 10.8%). Retourne les triangles gardes ou None."""
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
    return all_triangles[keep_mask]


def _weld_surface(surviving, tolerance: float):
    """Soude les triangles gardes en un maillage nettoye ; None si plus aucune face."""
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
    return resolved


def _main_component(resolved, tolerance: float):
    """Garde-fou de CONNEXITE : le filtre par enroulement peut laisser des fragments orphelins (turkish.stl :
    34 fragments de 1 a 4 faces, jusqu'a 76 mm du corps principal) ; TetGen pourrait alors ne mailler QUE un
    fragment sans erreur (regression silencieuse, limite #52). On ne garde donc que la plus grande composante
    (bruit residuel : 0.8% des faces ; recoller les fragments introduit de nouvelles auto-intersections). Le
    corps principal peut rester non etanche (turkish.stl : 3 boucles de 3 aretes) : `_close_small_boundary_loops`
    referme ces trous residuels par soudage (2 passes, connexite verifiee). Le nettoyage des auto-intersections
    locales qui en resultent revient a l'appelant (`step1_ingest._tetrahedralize`)."""
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

from element_ardor._self_intersection_geometry import (  # noqa: E402,F401  (reexport : imports et patchs existants)
    _boundary_edge_lines,
    _intersection_segment_lines,
    _triangulate_polygons,
)
