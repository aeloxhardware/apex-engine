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


# ---------------------------------------------------------------------------
# Decoupe conforme et extraction de l'enveloppe exterieure
# ---------------------------------------------------------------------------

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
    Decoupe un triangle par ses segments d'intersection, dans son propre
    plan (donc sans aucune perte : le triangle est plan par construction,
    tout point 2D se releve exactement en 3D).

    Le BORD n'est PAS passe comme un anneau d'un seul tenant mais comme
    trois aretes explicitement SUBDIVISEES aux points imposes
    (`extra_edge_points`). C'est ce qui garantit la CONFORMITE entre deux
    faces voisines : chacune doit porter exactement les memes sommets sur
    l'arete qu'elles partagent, sinon la surface recollee presente une
    fissure. Un mecanisme approximatif (materialiser le point par un
    micro-segment) a ete essaye puis abandonne - mesure : maillage final
    non etanche, euler=1 au lieu de 2.

    Retourne la liste des sous-triangles en 3D, ou None si la decoupe echoue.
    """
    from shapely.geometry import LineString, Polygon
    from shapely.ops import polygonize, unary_union
    from trimesh.creation import triangulate_polygon

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
        # Deduplication des points de la chaine, JAMAIS un saut de
        # sous-segment : sauter un sous-segment trop court TROUE l'anneau de
        # bord, et `polygonize` ne peut alors plus refermer la grande region
        # restante. Mesure du defaut : 4 faces sur 12 ne conservaient que 8 a
        # 70 % de leur aire d'origine, d'ou une surface finale non etanche.
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

    # Les courbes d'intersection doivent etre nodees GLOBALEMENT : un point
    # ou la courbe traverse une arete du maillage scinde le segment dans une
    # face mais pas forcement dans la face d'en face, ce qui cree une
    # JONCTION EN T (mesure : une arete [0.1,-0.63,1]->[0.1,0.94,1] d'un
    # cote, deux aretes coupees en [0.1,0.1,1] de l'autre). On reinjecte donc
    # dans chaque segment tous les points-noeuds qui tombent dessus.
    nodes_3d = node_points if node_points is not None else np.zeros((0, 3))
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
    return output or None


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
