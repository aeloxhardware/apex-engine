"""Résolution des auto-intersections : primitives géométriques (triangle-triangle, winding number, segments, subdivision de face).

Découpé de `_self_intersection_resolve.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `_self_intersection_resolve`.
"""

import numpy as np
import scipy.spatial
import trimesh
from core import _cpp


_EPS = 1e-10


# Tolerance d'equivalence pas 9 (triangle_triangle_segment) : le C++ resout
# le meme systeme 3x3 par cofacteurs plutot que LAPACK (np.linalg.solve),
# mesure a ~1e-14 relatif d'ecart - jamais le bit, voir
# cpp/include/apex/mesh_intersect.hpp. Tolerance choisie tres au-dessus de
# ce bruit et tres en dessous des tolerances propres du module (_EPS,
# snap_tolerance >= 1e-9 * scale).
_NATIVE_SEGMENT_TOLERANCE = 1e-9


def _triangle_planes(triangles: np.ndarray):
    """Normales (non normalisees) et offsets des plans porteurs."""
    normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    offsets = np.einsum("ij,ij->i", normals, triangles[:, 0])
    return normals, offsets


def _interval_on_line(triangle, distances, direction, scale):
    """
    Intervalle couvert par `triangle` sur la droite d'intersection des deux
    plans, parametree par `direction`.

    `distances` = distances signees des 3 sommets a l'AUTRE plan. On collecte
    directement les points ou le BORD du triangle rencontre ce plan : chaque
    sommet pose dessus, et chaque arete dont les deux extremites sont de part
    et d'autre. Formulation preferee a une logique de « sommet isole » :
    celle-ci degenere des qu'un sommet est exactement sur le plan (mesure :
    intervalle reduit a un point, donc intersection reelle manquee), cas
    frequent sur un maillage issu d'un booleen CSG.
    """
    points = []
    for index in range(3):
        following = (index + 1) % 3
        distance_index = distances[index]
        distance_following = distances[following]
        if abs(distance_index) <= scale:
            points.append(triangle[index])
        crosses = (
            (distance_index > scale and distance_following < -scale)
            or (distance_index < -scale and distance_following > scale)
        )
        if crosses:
            ratio = distance_index / (distance_index - distance_following)
            points.append(
                triangle[index] + ratio * (triangle[following] - triangle[index])
            )
    if len(points) < 2:
        return None  # contact ponctuel, pas une traversee
    projections = [float(np.dot(direction, point)) for point in points]
    return min(projections), max(projections)


def triangle_triangle_segment(triangle_a: np.ndarray, triangle_b: np.ndarray):
    """
    Segment d'intersection de deux triangles NON coplanaires, ou None.

    Methode de Moller : si les deux triangles se coupent, leur intersection
    est un segment porte par la droite d'intersection de leurs deux plans.
    On calcule l'intervalle couvert par chaque triangle sur cette droite,
    et l'intersection des deux intervalles donne le segment cherche.

    Retourne `(point_debut, point_fin)` ou None (pas d'intersection,
    contact ponctuel, ou triangles coplanaires - ce dernier cas est traite
    separement, une intersection coplanaire etant une AIRE et non un
    segment).
    """
    if _cpp.enabled():
        result = _cpp.module().triangle_triangle_segment(
            np.ascontiguousarray(triangle_a, dtype=np.float64),
            np.ascontiguousarray(triangle_b, dtype=np.float64),
        )
        if result is None:
            return None
        return result[0], result[1]

    normal_a = np.cross(triangle_a[1] - triangle_a[0], triangle_a[2] - triangle_a[0])
    normal_b = np.cross(triangle_b[1] - triangle_b[0], triangle_b[2] - triangle_b[0])
    norm_a = np.linalg.norm(normal_a)
    norm_b = np.linalg.norm(normal_b)
    if norm_a < _EPS or norm_b < _EPS:
        return None  # triangle degenere

    offset_a = float(np.dot(normal_a, triangle_a[0]))
    offset_b = float(np.dot(normal_b, triangle_b[0]))

    distances_b = triangle_b @ normal_a - offset_a
    distances_a = triangle_a @ normal_b - offset_b

    scale_a = _EPS * norm_a * float(np.max(np.abs(triangle_a)) + 1.0)
    scale_b = _EPS * norm_b * float(np.max(np.abs(triangle_b)) + 1.0)

    if np.all(distances_b > scale_a) or np.all(distances_b < -scale_a):
        return None
    if np.all(distances_a > scale_b) or np.all(distances_a < -scale_b):
        return None

    direction = np.cross(normal_a, normal_b)
    direction_norm = np.linalg.norm(direction)
    if direction_norm < _EPS:
        return None  # plans paralleles : cas coplanaire, traite ailleurs
    direction = direction / direction_norm

    interval_a = _interval_on_line(triangle_a, distances_a, direction, scale_b)
    interval_b = _interval_on_line(triangle_b, distances_b, direction, scale_a)
    if interval_a is None or interval_b is None:
        return None

    low = max(interval_a[0], interval_b[0])
    high = min(interval_a[1], interval_b[1])
    if high - low <= _EPS:
        return None  # contact ponctuel : pas une vraie traversee

    # un point de la droite : on repart d'un point du plan A projete
    reference = triangle_a[0]
    reference = reference - direction * float(np.dot(direction, reference))
    # recaler le point de reference pour qu'il appartienne aux deux plans
    basis = np.array([normal_a, normal_b, direction])
    target = np.array([offset_a, offset_b, 0.0])
    try:
        reference = np.linalg.solve(basis, target)
    except np.linalg.LinAlgError:
        return None
    return reference + direction * low, reference + direction * high


def generalized_winding_number(points: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """
    Nombre d'enroulement generalise de `points` par rapport a la surface
    fermee `triangles` (K, 3, 3), par somme d'angles solides (Van Oosterom
    et Strackee).

    Vaut ~1 a l'interieur d'une surface fermee simple, ~0 a l'exterieur, et
    ~2 dans une zone couverte deux fois - c'est cette derniere propriete qui
    le rend correct pour une surface qui se recouvre elle-meme, la ou un
    test de parite par lancer de rayon repondrait « exterieur ».

    Vectorise sur les points ET les triangles ; l'appelant doit garder le
    produit `len(points) * len(triangles)` raisonnable (usage vise : quelques
    centaines de points de test contre le maillage entier).
    """
    if _cpp.enabled():
        points_arr = np.ascontiguousarray(
            np.asarray(points, dtype=np.float64).reshape(-1, 3)
        )
        return _cpp.module().generalized_winding_number(
            points_arr, np.ascontiguousarray(triangles, dtype=np.float64)
        )

    points = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    a = triangles[None, :, 0, :] - points[:, None, :]
    b = triangles[None, :, 1, :] - points[:, None, :]
    c = triangles[None, :, 2, :] - points[:, None, :]

    norm_a = np.linalg.norm(a, axis=2)
    norm_b = np.linalg.norm(b, axis=2)
    norm_c = np.linalg.norm(c, axis=2)

    numerator = np.einsum("ijk,ijk->ij", a, np.cross(b, c))
    denominator = (
        norm_a * norm_b * norm_c
        + np.einsum("ijk,ijk->ij", a, b) * norm_c
        + np.einsum("ijk,ijk->ij", b, c) * norm_a
        + np.einsum("ijk,ijk->ij", c, a) * norm_b
    )
    return np.sum(np.arctan2(numerator, denominator), axis=1) / (2.0 * np.pi)


def _plane_frame(triangle):
    """Repere orthonorme du plan porteur d'un triangle : (origine, u, v, n)."""
    origin = triangle[0]
    normal = np.cross(triangle[1] - origin, triangle[2] - origin)
    norm = np.linalg.norm(normal)
    if norm < _EPS:
        return None
    normal = normal / norm
    axis_u = triangle[1] - origin
    axis_u_norm = np.linalg.norm(axis_u)
    if axis_u_norm < _EPS:
        return None
    axis_u = axis_u / axis_u_norm
    axis_v = np.cross(normal, axis_u)
    return origin, axis_u, axis_v, normal


def _project(points, frame):
    origin, axis_u, axis_v, _ = frame
    relative = np.asarray(points, dtype=np.float64).reshape(-1, 3) - origin
    return np.column_stack([relative @ axis_u, relative @ axis_v])


def _unproject(points_2d, frame):
    origin, axis_u, axis_v, _ = frame
    points_2d = np.asarray(points_2d, dtype=np.float64).reshape(-1, 2)
    return origin + points_2d[:, 0:1] * axis_u + points_2d[:, 1:2] * axis_v


def _intersection_segments(mesh, min_length):
    """
    Segments d'intersection par face.

    Les faces ADJACENTES (deux sommets communs, donc une arete partagee) sont
    ignorees : elles se rencontrent le long de cette arete par construction,
    ce n'est pas une auto-intersection. Les faces ne partageant qu'un SEUL
    sommet sont en revanche testees normalement - elles peuvent parfaitement
    se traverser ailleurs qu'en ce sommet, et le filtre de longueur minimale
    ecarte le simple contact ponctuel.
    """
    triangles = mesh.triangles
    faces = np.asarray(mesh.faces)

    if _cpp.enabled():
        # Broad-phase BVH maison + narrow-phase triangle_triangle_segment,
        # transcription de la boucle ci-dessous - voir
        # cpp/include/apex/mesh_intersect.hpp. Meme raisonnement qu'au pas 8
        # (mesh_repair.hpp) : le broad-phase natif n'a pas besoin d'etre
        # identique a bounds_tree (rtree), la narrow-phase est deterministe
        # des deux cotes (a la tolerance _NATIVE_SEGMENT_TOLERANCE pres).
        face_a, face_b, points = _cpp.module().find_intersection_segments(
            np.ascontiguousarray(triangles, dtype=np.float64),
            np.ascontiguousarray(faces, dtype=np.int64),
            float(min_length),
        )
        per_face: dict = {}
        for k in range(len(face_a)):
            segment = (points[k, 0], points[k, 1])
            per_face.setdefault(int(face_a[k]), []).append(segment)
            per_face.setdefault(int(face_b[k]), []).append(segment)
        return per_face

    tree = trimesh.triangles.bounds_tree(triangles)
    per_face = {}
    for i in range(len(triangles)):
        query = np.concatenate([triangles[i].min(axis=0), triangles[i].max(axis=0)])
        for j in tree.intersection(query):
            j = int(j)
            if j <= i:
                continue
            if np.intersect1d(faces[i], faces[j]).size >= 2:
                continue
            segment = triangle_triangle_segment(triangles[i], triangles[j])
            if segment is None:
                continue
            if float(np.linalg.norm(segment[1] - segment[0])) < min_length:
                continue
            per_face.setdefault(i, []).append(segment)
            per_face.setdefault(j, []).append(segment)
    return per_face


def _snap_to_triangle_boundary(point_2d, corners_2d, snap_tolerance):
    """
    Ramene un point sur le BORD du triangle s'il en est distant de moins de
    `snap_tolerance`, sinon le laisse tel quel.

    Indispensable a la robustesse : une extremite de segment d'intersection
    tombe sur le bord du triangle a l'erreur d'arrondi pres (mesure : 1e-16
    a 1e-17 d'ecart). Shapely raisonne sur des predicats EXACTS en double
    precision - un point a 1e-16 du bord n'est PAS sur le bord, la coupe
    reste donc pendante et `polygonize` ne scinde rien du tout (mesure : 1
    seul polygone rendu au lieu de 2, sur 2 des 4 faces testees).

    Priorite absolue aux SOMMETS du triangle (`corners_2d`), retournes tels
    quels (pas recalcules par interpolation) des qu'un sommet est dans
    `snap_tolerance` : sinon, un point proche d'un sommet peut se faire
    "snapper sur l'arete" via une interpolation `a + ratio*edge` qui,
    meme avec `ratio` clampe a 1.0, ne retombe PAS forcement bit-a-bit sur
    ce meme sommet obtenu depuis une AUTRE arete ou directement depuis
    `corners_2d` (deux chemins de calcul flottant differents pour la meme
    valeur mathematique). Ecart mesure sur turkish.stl : 2.8e-9 entre le
    sommet du triangle et le point "snappe sur arete" cense coincider avec
    lui - assez pour que shapely les traite comme deux points distincts et
    produise un polygone en "boucle en epingle" (une arete de longueur
    quasi nulle), source de triangles a rapport d'aspect >1000 en sortie de
    `_subdivide_face` (mesure : 35 slivers, TetGen les rejette comme
    intersections/coincidences "exactes").
    """
    for corner in corners_2d:
        if float(np.linalg.norm(point_2d - corner)) < snap_tolerance:
            return corner
    best_point = point_2d
    best_distance = snap_tolerance
    for index in range(3):
        a = corners_2d[index]
        b = corners_2d[(index + 1) % 3]
        edge = b - a
        length_squared = float(edge @ edge)
        if length_squared < _EPS:
            continue
        ratio = float((point_2d - a) @ edge) / length_squared
        ratio = min(1.0, max(0.0, ratio))
        projected = a + ratio * edge
        distance = float(np.linalg.norm(point_2d - projected))
        if distance < best_distance:
            best_distance = distance
            best_point = projected
    return best_point


def _weld_line_endpoints(lines, snap_tolerance, min_edge_length):
    """
    Fusionne tous les points EXTREMITES des `lines` (LineString) plus
    proches que `snap_tolerance` en un seul representant, et reconstruit les
    lignes avec les coordonnees fusionnees.

    Necessaire au-dela du seul cas "point proche d'un sommet du triangle"
    deja traite par `_snap_to_triangle_boundary` : deux SEGMENTS
    d'intersection distincts peuvent chacun atterrir a quelques 1e-4 unites
    l'un de l'autre pres d'un meme sommet (deux courbes 3D differentes,
    physiquement proches mais pas identiques), sans qu'aucun des deux points
    ne soit assez proche du sommet lui-meme pour declencher ce snap. Mesure
    sur turkish.stl (face 9134) : deux segments a 0.00027 et 0.00042 unites
    du coin, jamais fusionnes entre eux, produisaient un anneau shapely en
    "epingle" (une arete de longueur EXACTEMENT nulle) donc un triangle
    degenere en sortie de `triangulate_polygon` - rapport d'aspect mesure
    jusqu'a 7.6e14 (zero numerique, pas un simple triangle fin). Un
    clustering generique sur TOUS les points de bord (pas seulement ceux
    pres d'un sommet) resout cette classe entiere de defaut d'un coup.
    """
    from shapely.geometry import LineString

    points = []
    for line in lines:
        coords = list(line.coords)
        points.append(coords[0])
        points.append(coords[-1])
    if not points:
        return lines
    points = np.asarray(points, dtype=np.float64)
    tree = scipy.spatial.cKDTree(points)
    pairs = tree.query_pairs(r=snap_tolerance)
    parent = list(range(len(points)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in pairs:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    representative = {}
    for i in range(len(points)):
        root = find(i)
        if root not in representative:
            representative[root] = points[root]
        points[i] = representative[root]

    welded_lines = []
    for index, line in enumerate(lines):
        start = points[2 * index]
        end = points[2 * index + 1]
        if float(np.linalg.norm(start - end)) < min_edge_length:
            continue
        welded_lines.append(LineString([start, end]))
    return welded_lines


def _split_keyhole_polygon(polygon, tolerance):
    """
    Convertit un polygone shapely (avec eventuels trous) en une liste de
    polygones SIMPLES (sans trou), en epissant chaque trou qui TOUCHE le
    contour exterieur en un seul sommet partage - un cas de « trou en trou
    de serrure » degenere pour tout triangulateur base sur des contraintes.

    Mesure sur turkish.stl (face 8883) : `polygonize` peut legitimement
    produire un polygone dont le trou touche le contour exterieur en un
    UNIQUE sommet (deux boucles issues de la meme courbe d'intersection,
    partageant un point de la boucle qui se referme sur elle-meme).
    `triangulate_polygon(..., force_vertices=True OU False)` echoue
    silencieusement sur ce cas precis : il introduit un sommet DUPLIQUE
    (le point de contact compte deux fois) relie par une arete de longueur
    EXACTEMENT nulle - rapport d'aspect mesure jusqu'a 7.6e14 sur le
    triangle resultant. Aucune des deux methodes ne l'evite.

    Fix : repere le(s) sommet(s) de contact (trou et exterieur partageant
    un point a `tolerance` pres), et REINJECTE le trou comme une excursion
    aller-retour dans le contour exterieur a ce point plutot que comme un
    anneau interieur separe - le resultat est un polygone SIMPLE (une seule
    boucle, sans trou), que tout triangulateur traite correctement. Un trou
    qui NE touche PAS le contour exterieur (cas normal, non degenere) est
    laisse tel quel, encode via `holes=` (shapely/trimesh le gerent
    correctement dans ce cas, seul le contact ponctuel est pathologique).
    """
    from shapely.geometry import Polygon

    if not polygon.interiors:
        return [polygon]

    exterior = list(polygon.exterior.coords[:-1])
    remaining_holes = []
    for interior in polygon.interiors:
        hole_coords = list(interior.coords[:-1])
        touch_index = None
        exterior_index = None
        for h_index, hole_point in enumerate(hole_coords):
            for e_index, ext_point in enumerate(exterior):
                if (
                    abs(hole_point[0] - ext_point[0]) < tolerance
                    and abs(hole_point[1] - ext_point[1]) < tolerance
                ):
                    touch_index = h_index
                    exterior_index = e_index
                    break
            if touch_index is not None:
                break
        if touch_index is None:
            remaining_holes.append(hole_coords)
            continue
        # Reordonne la boucle du trou pour commencer au point de contact
        # (exclu, deja porte par `exterior[exterior_index]`), l'insere comme
        # un aller-retour dans le contour exterieur, sans jamais repeter ce
        # sommet de contact plus d'une fois (sinon `Polygon` le voit comme
        # un anneau invalide - mesure : `is_valid=False`, le point y
        # apparaissait 3 fois de suite).
        reordered_hole = hole_coords[touch_index + 1 :] + hole_coords[:touch_index]
        exterior = (
            exterior[: exterior_index + 1]
            + reordered_hole
            + exterior[exterior_index:]
        )

    if not exterior or len(exterior) < 3:
        return [polygon] if not remaining_holes else []
    simple = Polygon(exterior, holes=remaining_holes or None)
    # `simple.is_valid` est souvent FAUX ici au sens strict OGC (l'anneau se
    # retouche lui-meme au point de contact, un "noeud papillon" degenere) -
    # mais `triangulate_polygon` (mapbox_earcut/`triangle`) le traite
    # correctement malgre tout (verifie sur turkish.stl, face 8883 :
    # 6 triangles produits, 0 arete degeneree). Rejeter sur `is_valid`
    # reviendrait a retomber sur le polygone-a-trou d'origine, exactement le
    # cas qui echoue. Seule l'aire (garde-fou contre un anneau effondre) est
    # verifiee.
    if simple.area <= tolerance ** 2:
        return [polygon]
    return [simple]


def _subdivide_face(triangle, segments, extra_edge_points, tolerance, node_points=None):
    """
    Decoupe un triangle par ses segments d'intersection, dans son propre plan (donc sans aucune perte : le
    triangle est plan par construction, tout point 2D se releve exactement en 3D).

    Le BORD n'est PAS passe comme un anneau d'un seul tenant mais comme trois aretes explicitement SUBDIVISEES
    aux points imposes (`extra_edge_points`) : c'est ce qui garantit la CONFORMITE entre deux faces voisines
    (memes sommets sur l'arete partagee, sinon fissure). Un mecanisme approximatif (micro-segment) a ete
    essaye puis abandonne - mesure : maillage final non etanche, euler=1 au lieu de 2.

    Retourne la liste des sous-triangles en 3D, ou None si la decoupe echoue.
    """
    from shapely.ops import polygonize, unary_union

    frame = _plane_frame(triangle)
    if frame is None:
        return None
    corners_2d = _project(triangle, frame)
    snap_tolerance = max(
        tolerance, 1e-9 * float(np.max(np.abs(corners_2d)) + 1.0)
    )
    extras_2d = (
        _project(np.asarray(extra_edge_points), frame)
        if len(extra_edge_points) else np.zeros((0, 2))
    )

    lines = _boundary_edge_lines(corners_2d, extras_2d, snap_tolerance, tolerance)
    nodes_3d = node_points if node_points is not None else np.zeros((0, 3))
    lines.extend(_intersection_segment_lines(segments, nodes_3d, frame, corners_2d, snap_tolerance, tolerance))

    if not lines:
        return None
    lines = _weld_line_endpoints(lines, snap_tolerance, tolerance)
    if not lines:
        return None
    try:
        arrangement = unary_union(lines)
        polygons = [p for p in polygonize(arrangement) if p.area > tolerance ** 2]
    except Exception:
        return None
    if not polygons:
        return None

    output = _triangulate_polygons(polygons, corners_2d, frame, snap_tolerance, tolerance)
    return output or None


def _boundary_edge_lines(corners_2d, extras_2d, snap_tolerance: float, tolerance: float) -> list:
    """Les trois aretes du triangle, explicitement SUBDIVISEES aux points imposes (`extras_2d`)."""
    from shapely.geometry import LineString

    lines = []
    for index in range(3):
        a = corners_2d[index]
        b = corners_2d[(index + 1) % 3]
        edge = b - a
        length_squared = float(edge @ edge)
        if length_squared < _EPS:
            continue
        chain = [(0.0, a), (1.0, b)]
        for point in extras_2d:
            ratio = float((point - a) @ edge) / length_squared
            if ratio <= 0.0 or ratio >= 1.0:
                continue
            projected = a + ratio * edge
            if float(np.linalg.norm(point - projected)) > snap_tolerance:
                continue
            chain.append((ratio, projected))
        # Deduplication des points de la chaine, JAMAIS un saut de sous-segment : sauter un sous-segment trop
        # court TROUE l'anneau de bord, et `polygonize` ne peut plus refermer la grande region restante
        # (mesure : 4 faces sur 12 ne conservaient que 8 a 70 % de leur aire, surface finale non etanche).
        chain.sort(key=lambda item: item[0])
        deduplicated = [chain[0]]
        for item in chain[1:]:
            if float(np.linalg.norm(item[1] - deduplicated[-1][1])) >= tolerance:
                deduplicated.append(item)
        if float(np.linalg.norm(deduplicated[-1][1] - b)) >= tolerance:
            deduplicated.append((1.0, b))
        else:
            deduplicated[-1] = (1.0, b)
        for position in range(len(deduplicated) - 1):
            lines.append(LineString([
                deduplicated[position][1], deduplicated[position + 1][1]
            ]))
    return lines


def _intersection_segment_lines(segments, nodes_3d, frame, corners_2d, snap_tolerance: float,
                                tolerance: float) -> list:
    """Segments d'intersection nodes GLOBALEMENT : un point ou la courbe traverse une arete du maillage scinde
    le segment dans une face mais pas forcement dans la face d'en face (JONCTION EN T, mesure : une arete
    [0.1,-0.63,1]->[0.1,0.94,1] d'un cote, deux aretes coupees en [0.1,0.1,1] de l'autre). On reinjecte donc
    dans chaque segment tous les points-noeuds qui tombent dessus."""
    from shapely.geometry import LineString

    lines = []
    for (start, end) in segments:
        direction = end - start
        length = float(np.linalg.norm(direction))
        if length < tolerance:
            continue
        unit = direction / length
        chain = [0.0, length]
        for node in nodes_3d:
            offset = node - start
            along = float(np.dot(offset, unit))
            if along <= tolerance or along >= length - tolerance:
                continue
            if float(np.linalg.norm(offset - along * unit)) > snap_tolerance:
                continue
            chain.append(along)
        chain.sort()
        deduplicated = [chain[0]]
        for value in chain[1:]:
            if value - deduplicated[-1] >= tolerance:
                deduplicated.append(value)
        deduplicated[-1] = length
        points_3d = [start + unit * value for value in deduplicated]
        projected = _project(np.asarray(points_3d), frame)
        projected[0] = _snap_to_triangle_boundary(projected[0], corners_2d, snap_tolerance)
        projected[-1] = _snap_to_triangle_boundary(projected[-1], corners_2d, snap_tolerance)
        for position in range(len(projected) - 1):
            if float(np.linalg.norm(projected[position + 1] - projected[position])) < tolerance:
                continue
            lines.append(LineString([projected[position], projected[position + 1]]))
    return lines


def _triangulate_polygons(polygons, corners_2d, frame, snap_tolerance: float, tolerance: float):
    """Triangule les polygones dont le point representatif est dans le triangle ; None si l'un echoue."""
    from shapely.geometry import Polygon
    from trimesh.creation import triangulate_polygon

    triangle_polygon = Polygon(corners_2d)
    output = []
    for polygon in polygons:
        if not polygon.representative_point().within(triangle_polygon.buffer(snap_tolerance)):
            continue
        for simple_ring in _split_keyhole_polygon(polygon, tolerance):
            vertices_2d = None
            for force in (True, False):
                try:
                    vertices_2d, faces_2d = triangulate_polygon(simple_ring, force_vertices=force)
                    break
                except Exception:
                    vertices_2d = None
            if vertices_2d is None:
                return None
            vertices_3d = _unproject(vertices_2d, frame)
            for face in faces_2d:
                output.append(vertices_3d[face])
    return output
