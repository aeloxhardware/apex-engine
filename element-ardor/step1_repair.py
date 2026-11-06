"""Etape 1 : réparation du maillage source (auto-intersections, trous, coque dominante).

Découpé de `step1_ingest.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step1_ingest`.
"""

from __future__ import annotations

from typing import Optional
import math
import numpy as np
import trimesh
from core import _cpp
from pipeline.step1_types import (  # noqa: F401
    IngestReport,
)


# Tolerance du test de separation par plan (_triangles_cross) : en dessous,
# un sommet est considere SUR le plan de l'autre triangle, pas d'un cote.
# Deux faces adjacentes coplanaires se toucheraient sinon a tort.
_COPLANAR_EPS = 1e-9


def _count_self_intersecting_faces(mesh: trimesh.Trimesh) -> int:
    """
    Nombre de faces impliquees dans au moins une auto-intersection reelle.
    Enveloppe mince autour de `_find_self_intersecting_faces` : le COMPTE
    seul suffit au rapport d'ingestion, mais la reparation locale
    (`_repair_local_self_intersections`) a besoin des INDICES.
    """
    return len(_find_self_intersecting_faces(mesh))


def _find_self_intersecting_faces(mesh: trimesh.Trimesh) -> set[int]:
    """
    Indices des faces impliquees dans au moins une auto-intersection reelle.

    Broad-phase par arbre d'AABB (`trimesh.triangles.bounds_tree`) puis test
    exact triangle-triangle par separation de plan : deux triangles ne se
    traversent pas si l'un est entierement d'un cote du plan de l'autre.
    Les paires partageant un sommet sont ignorees - elles se touchent par
    construction, ce n'est pas une auto-intersection.

    Compte les FACES et non les paires : c'est ce qui permet de dire si le
    defaut est local (quelques faces, reparable par suppression) ou global
    (une fraction importante du maillage, cf. mesure sur death.stl).
    """
    triangles = mesh.triangles
    faces = mesh.faces
    if len(triangles) == 0:
        return set()

    if _cpp.enabled():
        # Broad-phase BVH maison + narrow-phase SAT, transcription bit-exacte
        # de la boucle ci-dessous - voir cpp/include/apex/mesh_repair.hpp.
        # Le broad-phase natif n'a pas besoin d'etre identique a bounds_tree
        # (rtree) : un sur-ensemble conservateur d'AABB exactes suffit, la
        # narrow-phase (SAT) est deterministe et bit-exacte des deux cotes.
        bad_array = _cpp.module().find_self_intersecting_faces(
            np.ascontiguousarray(triangles, dtype=np.float64),
            np.ascontiguousarray(faces, dtype=np.int64),
        )
        return {int(i) for i in bad_array}

    tree = trimesh.triangles.bounds_tree(triangles)
    bad: set[int] = set()

    for i, triangle in enumerate(triangles):
        query = np.concatenate([triangle.min(axis=0), triangle.max(axis=0)])
        for j in tree.intersection(query):
            j = int(j)
            if j <= i:
                continue
            # sommet partage => contact par construction, pas une intersection
            if np.intersect1d(faces[i], faces[j]).size > 0:
                continue
            if _triangles_cross(triangles[i], triangles[j]):
                bad.add(i)
                bad.add(j)

    return bad


def _triangles_cross(a: np.ndarray, b: np.ndarray) -> bool:
    """
    True si deux triangles s'interpenetrent reellement.

    Theoreme de l'axe separateur (SAT) : deux convexes sont disjoints des
    qu'il existe UN axe sur lequel leurs projections ne se recouvrent pas.
    Pour deux triangles, il suffit de tester les 2 normales de face et les
    9 produits vectoriels arete x arete.

    Le test des seules normales de face (condition NECESSAIRE mais pas
    SUFFISANTE) produisait des faux positifs massifs, mesures : 41 paires
    sur res/test_models/l-triple-axis-overhanged-shape.stl, un modele que
    TetGen accepte pourtant sans broncher - typiquement deux faces qui se
    TOUCHENT le long d'une arete sans se traverser. D'ou les 9 axes croises.

    La comparaison utilise `>= -_COPLANAR_EPS` : un simple contact (arete
    partagee, faces coplanaires jointives) laisse les projections se
    toucher en un point sans recouvrement d'interieur, et n'est donc PAS
    compte comme une auto-intersection.
    """
    if _cpp.enabled():
        return bool(_cpp.module().triangles_cross(
            np.ascontiguousarray(a, dtype=np.float64), np.ascontiguousarray(b, dtype=np.float64)
        ))

    edges_a = [a[1] - a[0], a[2] - a[1], a[0] - a[2]]
    edges_b = [b[1] - b[0], b[2] - b[1], b[0] - b[2]]

    axes = [np.cross(edges_a[0], edges_a[1]), np.cross(edges_b[0], edges_b[1])]
    axes.extend(np.cross(ea, eb) for ea in edges_a for eb in edges_b)

    for axis in axes:
        norm = np.linalg.norm(axis)
        if norm < 1e-12:  # aretes paralleles (ou face degeneree) : axe inutile
            continue
        axis = axis / norm
        pa, pb = a @ axis, b @ axis
        # separes si un intervalle finit avant que l'autre ne commence
        if pa.min() >= pb.max() - _COPLANAR_EPS or pb.min() >= pa.max() - _COPLANAR_EPS:
            return False

    return True


def _repair_with_pymeshfix(mesh: trimesh.Trimesh) -> trimesh.Trimesh | None:
    """
    Reconstruit une surface propre via `pymeshfix` (wrapper de MeshFix) -
    corrige generiquement les auto-intersections ET les trous/defauts de
    non-manifold que `trimesh.repair.fill_holes` ne comble pas (mesure sur
    `res/test_models/turkish.stl` : 64 faces cassees, mesh reste
    non-watertight apres fill_holes/fix_normals/fix_winding standard,
    corrige par MeshFix.repair()). Anciennement nommee
    `_repair_self_intersections` - renommee pour refleter cette portee plus
    large, deja celle du mecanisme sous-jacent (MeshFix ne distingue pas
    ces deux classes de defauts, seul son APPELANT le faisait par le nom de
    la fonction).

    Dependance OPTIONNELLE, importee paresseusement : un environnement sans
    `pymeshfix` continue de fonctionner exactement comme avant sur tous les
    maillages deja valides. Retourne None si la bibliotheque est absente ou
    si la reparation echoue - l'appelant decide alors quoi faire, cette
    fonction ne leve jamais.

    `remove_smallest_components=False` : supprimer les composantes annexes
    reviendrait a jeter silencieusement des morceaux du modele de
    l'utilisateur, ce qui n'est pas a la main de cette etape.
    """
    try:
        from pymeshfix import MeshFix
    except ImportError:
        return None

    try:
        fixer = MeshFix(
            np.asarray(mesh.vertices, dtype=np.float64),
            np.asarray(mesh.faces, dtype=np.int32),
        )
        fixer.repair(remove_smallest_components=False)
        repaired = trimesh.Trimesh(
            vertices=np.asarray(fixer.points), faces=np.asarray(fixer.faces), process=False
        )
    except Exception:
        return None

    if repaired.faces.shape[0] == 0:
        return None
    return repaired


# La croissance de la zone retiree est ADAPTATIVE et PAR TROU : on n'elargit
# que les sites dont le rebouchage a echoue. Mesure directe sur turkish.stl :
# elargir GLOBALEMENT d'un anneau DEGRADE le resultat (152 faces retirees et
# 9 boucles en echec, contre 50 faces et 4 echecs sans elargissement) parce
# que des sites de defaut independants fusionnent en un seul grand trou dont
# le bord contient alors les DEUX nervures qui se croisent - donc forcement
# auto-intersectant. Le retrait minimal n'est pas un compromis de fidelite,
# il est strictement plus correct.
_LOCAL_REPAIR_MAX_GROWTH_ITERATIONS = 8


_LOCAL_REPAIR_MAX_REMOVED_RATIO = 0.05


_LOCAL_REPAIR_MIN_TRIANGLE_AREA = 1e-12


def _build_half_edge_map(faces: np.ndarray):
    """
    `(u, v) -> (indice de face, coin k)` pour chaque demi-arete dirigee.

    Retourne None si une meme arete DIRIGEE apparait deux fois : le maillage
    est alors non-manifold et le pivot autour d'un sommet serait indefini.
    Echec explicite plutot qu'un resultat devine.
    """
    half_edges: dict[tuple[int, int], tuple[int, int]] = {}
    for face_index in range(len(faces)):
        triangle = faces[face_index]
        for corner in range(3):
            key = (int(triangle[corner]), int(triangle[(corner + 1) % 3]))
            if key in half_edges:
                return None
            half_edges[key] = (face_index, corner)
    return half_edges


def _boundary_loops_oriented(faces: np.ndarray):
    """
    Boucles de bord ORIENTEES d'un maillage triangulaire ouvert.

    Chaque arete de bord appartient a exactement UNE face, dont le winding
    lui donne une direction canonique. Le successeur de `(u, v)` s'obtient
    en PIVOTANT autour de `v` a travers l'eventail de faces incidentes
    jusqu'a la prochaine demi-arete sans face opposee. C'est ce pivot qui
    leve l'ambiguite a un sommet de PINCEMENT (ou plusieurs eventails se
    rejoignent) : il apparie les demi-aretes entrantes et sortantes PAR
    EVENTAIL, information qu'un graphe non oriente a deja perdue.

    `trimesh.outline()` ne convient PAS ici (verifie par lecture directe des
    sources, pas suppose) : il passe par `edges_sorted` - le tri detruit la
    direction portee par le winding - puis par un parcours en profondeur sur
    graphe NON ORIENTE (`trimesh.graph.traversals`, `directed=False`), et
    rend des polylignes qui ne sont meme pas garanties fermees.

    Retourne une liste de boucles (listes d'indices de sommets, sans point
    de fermeture duplique), ou None si le maillage est non-manifold ou si un
    invariant de correction est viole.
    """
    half_edges = _build_half_edge_map(faces)
    if half_edges is None:
        return None

    def next_in_face(half_edge):
        face_index, corner = half_edges[half_edge]
        triangle = faces[face_index]
        return (int(triangle[(corner + 1) % 3]), int(triangle[(corner + 2) % 3]))

    boundary = [h for h in half_edges if (h[1], h[0]) not in half_edges]
    if not boundary:
        return []

    successor: dict = {}
    for start in boundary:
        current = start
        guard = 0
        while True:
            guard += 1
            if guard > len(half_edges) + 1:
                return None  # eventail non ferme : geometrie incoherente
            current = next_in_face(current)
            if (current[1], current[0]) not in half_edges:
                successor[start] = current
                break
            current = (current[1], current[0])

    # `successor` doit etre une BIJECTION sur les demi-aretes de bord :
    # invariant de correction bon marche et tres discriminant.
    if len(set(successor.values())) != len(boundary):
        return None

    visited: set = set()
    cycles = []
    for start in sorted(boundary):
        if start in visited:
            continue
        cycle = []
        current = start
        while current not in visited:
            visited.add(current)
            cycle.append(current)
            current = successor[current]
        if current != start:
            return None
        cycles.append(cycle)
    if sum(len(c) for c in cycles) != len(boundary):
        return None

    # Un cycle peut visiter deux fois le meme sommet : c'est CORRECT (le bord
    # se touche en un point) mais la triangulation le refuse. Scission en
    # sous-boucles simples ; l'ordre du cycle encode deja l'appariement par
    # eventail, donc cette scission est univoque - contrairement a celle d'un
    # parcours non oriente, qui devine.
    loops = []
    for cycle in cycles:
        sequence = [h[0] for h in cycle]
        stack: list = []
        position: dict = {}
        for vertex in sequence:
            if vertex in position:
                start_index = position[vertex]
                sub_loop = stack[start_index:]
                for repeated in sub_loop:
                    position.pop(repeated, None)
                del stack[start_index:]
                if len(sub_loop) >= 3:
                    loops.append(sub_loop)
                elif sub_loop:
                    return None  # eperon degenere : jamais ignore en silence
            position[vertex] = len(stack)
            stack.append(vertex)
        if len(stack) >= 3:
            loops.append(stack)
        elif stack:
            return None
    return loops


def _fill_loop_min_weight(
    loop,
    vertices: np.ndarray,
    blocked_edges: set,
    obstacle_faces: np.ndarray,
    obstacle_triangles: np.ndarray,
    obstacle_tree,
    accepted_triangles: list,
):
    """
    Triangulation a POIDS MINIMAL de la boucle `loop` (Klincsek, O(n^3)).

    `best_cost[i][j]` = cout minimal de trianguler le sous-contour `i..j`,
    cout d'un triangle = son aire. C'est l'algorithme classique de bouchage
    de trou (Barequet-Sharir) : il choisit une triangulation respectant la
    forme 3D reelle de la boucle, au lieu de l'aplatir sur un plan moyen.

    Deux contraintes sont integrees DANS le DP sous forme de poids infini,
    et non verifiees apres coup - un patch invalide rend donc la
    triangulation infaisable au lieu d'etre accepte puis regrette :
      1. UNICITE D'ARETE : une corde ne doit pas dupliquer une arete deja
         presente. Sans elle, mesure directe : des aretes se retrouvaient
         portees par 4 faces et le maillage n'etait pas etanche.
      2. NON-INTERSECTION : le triangle candidat ne doit croiser ni la
         surface conservee ni un patch deja accepte (SAT, `_triangles_cross`).
         C'est exactement la verification qui manquait au repli en eventail,
         lequel creait de NOUVELLES auto-intersections (mesure : TetGen
         nommait les sommets d'eventail synthetiques comme fautifs).

    Retourne la liste des triangles `(a, b, c)`, ou None si aucune
    triangulation valide n'existe.
    """
    n = len(loop)
    if n < 3:
        return None

    weight_cache: dict = {}

    def triangle_weight(i: int, k: int, j: int) -> float:
        key = (i, k, j)
        cached = weight_cache.get(key)
        if cached is not None:
            return cached
        a, b, c = loop[i], loop[k], loop[j]
        pa, pb, pc = vertices[a], vertices[b], vertices[c]
        area = 0.5 * float(np.linalg.norm(np.cross(pb - pa, pc - pa)))
        result = area
        if area < _LOCAL_REPAIR_MIN_TRIANGLE_AREA:
            result = math.inf
        else:
            for (x, y) in ((i, k), (k, j), (i, j)):
                separation = abs(x - y)
                if separation == 1 or separation == n - 1:
                    continue  # arete de bord de la boucle, pas une corde
                if tuple(sorted((loop[x], loop[y]))) in blocked_edges:
                    result = math.inf
                    break
        if result != math.inf:
            candidate = np.array([pa, pb, pc])
            own = {a, b, c}
            query = np.concatenate([candidate.min(axis=0), candidate.max(axis=0)])
            for other in obstacle_tree.intersection(query):
                other = int(other)
                if own & {int(v) for v in obstacle_faces[other]}:
                    continue  # sommet partage : contact par construction
                if _triangles_cross(candidate, obstacle_triangles[other]):
                    result = math.inf
                    break
            if result != math.inf:
                for (oa, ob, oc) in accepted_triangles:
                    if own & {oa, ob, oc}:
                        continue
                    if _triangles_cross(candidate, vertices[[oa, ob, oc]]):
                        result = math.inf
                        break
        weight_cache[key] = result
        return result

    best_cost = [[math.inf] * n for _ in range(n)]
    split_at = [[-1] * n for _ in range(n)]
    for i in range(n - 1):
        best_cost[i][i + 1] = 0.0
    for gap in range(2, n):
        for i in range(0, n - gap):
            j = i + gap
            best = math.inf
            best_k = -1
            for k in range(i + 1, j):
                left = best_cost[i][k]
                right = best_cost[k][j]
                if left == math.inf or right == math.inf:
                    continue
                own_cost = triangle_weight(i, k, j)
                if own_cost == math.inf:
                    continue
                total = left + right + own_cost
                if total < best:
                    best = total
                    best_k = k
            best_cost[i][j] = best
            split_at[i][j] = best_k
    if best_cost[0][n - 1] == math.inf:
        return None

    triangles = []
    stack = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if j - i < 2:
            continue
        k = split_at[i][j]
        if k < 0:
            return None
        # Winding INVERSE par rapport au sens de parcours du bord : c'est ce
        # qui rend le maillage recolle coherent (mesure directe - emettre
        # dans le sens de la boucle donne is_winding_consistent=False).
        triangles.append((loop[j], loop[k], loop[i]))
        stack.append((i, k))
        stack.append((k, j))

    # Filet de securite : deux triangles du MEME patch ne doivent pas se
    # croiser non plus (le DP les evalue isolement, jamais deux a deux).
    for x in range(len(triangles)):
        for y in range(x + 1, len(triangles)):
            if set(triangles[x]) & set(triangles[y]):
                continue
            if _triangles_cross(vertices[list(triangles[x])], vertices[list(triangles[y])]):
                return None
    return triangles


def _repair_local_self_intersections(
    mesh: trimesh.Trimesh, extra_seed_faces=None, require_euler_preserved: bool = True,
):
    """
    Retire les faces auto-intersectantes et rebouche les trous par
    triangulation 3D a poids minimal, en laissant TOUT le reste du maillage
    strictement intact.

    `require_euler_preserved=False` desactive UNIQUEMENT le garde-fou de
    caracteristique d'Euler ci-dessous (le garde de CONNEXITE - nombre de
    corps preserve - reste actif dans tous les cas). Necessaire pour un
    appelant dont le mesh d'entree a DEJA subi un changement topologique
    legitime avant cet appel (ex : `resolve_self_intersections` +
    `_close_small_boundary_loops`, qui ferme des trous residuels par
    soudage de sommets - une operation qui change l'euler par construction,
    pas une mutilation). Mesure sur turkish.stl : reparation de 8 faces
    auto-intersectantes en 1 iteration, corps unique preserve, watertight,
    winding coherent, 0 auto-intersection residuelle, ecart de volume
    2.9e-6% - un resultat de tres haute fidelite que la garde stricte
    rejetait a tort (euler -37 au lieu de -36, un decalage de 1 exactement
    attendu pour la fermeture de boucles de tailles differentes de celles
    retirees).

    Contrairement a `_repair_with_pymeshfix` (qui reecrit le maillage entier
    et peut en amputer une grande part - mesure sur turkish.stl : 10104 ->
    1488 faces, hauteur tronquee de 90.78 a 63.26mm), cette reparation est
    CHIRURGICALE : seules les faces effectivement fautives sont retirees, et
    la zone n'est elargie que la ou le rebouchage echoue, trou par trou.

    Mesure de reference sur `res/test_models/turkish.stl` : 50 faces
    auto-intersectantes detectees, 99 faces retirees (0.98% du maillage) et
    93 ajoutees apres convergence en 4 iterations ; etendue Z et
    `euler_number` preserves a l'identique, ecart de volume 3.6e-4 %,
    50 -> 0 auto-intersections.

    Best-effort : retourne None (jamais une exception) si le defaut n'est
    pas de cette nature - maillage non-manifold, boucle non rebouchable,
    budget de retrait depasse, ou post-condition violee.
    """
    faces = np.asarray(mesh.faces, dtype=np.int64)
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    if len(faces) == 0:
        return None

    initial_bad = set(_find_self_intersecting_faces(mesh))
    if extra_seed_faces:
        # Le test SAT ignore volontairement les paires partageant un sommet
        # et ne voit pas les CONTACTS EXACTS (triangles coplanaires
        # superposes, sommet pose sur une facette) - or c'est precisement ce
        # que TetGen refuse. Mesure directe sur turkish.stl : apres avoir
        # repare les 50 faces SAT, TetGen en signalait encore 40. Les faces
        # qu'il nomme lui-meme sont donc indispensables a la graine.
        initial_bad |= {int(i) for i in extra_seed_faces if 0 <= int(i) < len(faces)}
    if not initial_bad:
        return None  # rien a reparer : le maillage n'est jamais modifie

    max_removed = max(1, int(_LOCAL_REPAIR_MAX_REMOVED_RATIO * len(faces)))
    removed = {int(i) for i in initial_bad}
    _body_count_before = len(mesh.split(only_watertight=False))

    edge_to_faces: dict = {}
    for face_index in range(len(faces)):
        triangle = faces[face_index]
        for corner in range(3):
            key = (int(triangle[corner]), int(triangle[(corner + 1) % 3]))
            edge_to_faces.setdefault(tuple(sorted(key)), []).append(face_index)

    for _iteration in range(_LOCAL_REPAIR_MAX_GROWTH_ITERATIONS):
        if len(removed) > max_removed:
            return None

        kept_index = np.array(
            [i for i in range(len(faces)) if i not in removed], dtype=np.int64
        )
        if len(kept_index) == 0:
            return None
        kept_faces = faces[kept_index]

        loops = _boundary_loops_oriented(kept_faces)
        if loops is None:
            return None

        kept_triangles = vertices[kept_faces]
        obstacle_tree = trimesh.triangles.bounds_tree(kept_triangles)
        blocked_edges = set()
        for triangle in kept_faces:
            for corner in range(3):
                blocked_edges.add(tuple(sorted(
                    (int(triangle[corner]), int(triangle[(corner + 1) % 3]))
                )))

        accepted: list = []
        failed = []
        for loop in loops:
            patch = _fill_loop_min_weight(
                loop, vertices, blocked_edges, kept_faces, kept_triangles,
                obstacle_tree, accepted,
            )
            if patch is None:
                failed.append(loop)
                continue
            accepted.extend(patch)
            for (a, b, c) in patch:
                for (x, y) in ((a, b), (b, c), (c, a)):
                    blocked_edges.add(tuple(sorted((int(x), int(y)))))

        if not failed:
            if accepted:
                all_faces = np.vstack([kept_faces, np.array(accepted, dtype=np.int64)])
            else:
                all_faces = kept_faces
            repaired = trimesh.Trimesh(
                vertices=vertices.copy(), faces=all_faces, process=False
            )
            repaired.remove_unreferenced_vertices()
            if not repaired.is_watertight or not repaired.is_winding_consistent:
                return None
            # La TOPOLOGIE doit etre preservee a l'identique. Sans ce
            # garde-fou, mesure directe sur turkish.stl : le retrait finit
            # par DETACHER un morceau (euler -14 -> -13, 1 corps -> 2), et
            # TetGen echoue alors sur `recoversubfaces` SANS plus signaler
            # la moindre auto-intersection - un echec devenu indechiffrable.
            # Une reparation qui change le nombre de corps ou la
            # caracteristique d'Euler ne repare plus : elle mutile.
            if require_euler_preserved and repaired.euler_number != mesh.euler_number:
                return None
            if len(repaired.split(only_watertight=False)) != _body_count_before:
                return None
            if _count_self_intersecting_faces(repaired) >= len(initial_bad):
                return None
            return repaired

        # Croissance CIBLEE : uniquement les faces bordant les boucles en
        # echec. Un trou a 3 sommets n'a qu'un seul patch possible - la face
        # qu'on vient justement de retirer - donc son bord DOIT s'eloigner
        # du croisement pour devenir rebouchable.
        grew = False
        for loop in failed:
            size = len(loop)
            for index in range(size):
                key = tuple(sorted((loop[index], loop[(index + 1) % size])))
                for face_index in edge_to_faces.get(key, ()):
                    if face_index not in removed:
                        removed.add(face_index)
                        grew = True
        if not grew:
            return None
    return None


def _clean_and_repair(mesh: trimesh.Trimesh) -> tuple[trimesh.Trimesh, IngestReport]:
    was_watertight_initially = bool(mesh.is_watertight)

    vertex_count_before = mesh.vertices.shape[0]
    face_count_before = mesh.faces.shape[0]

    mesh.merge_vertices()
    # API trimesh >= 5.0 : ces methodes retournent un masque booleen a
    # appliquer via update_faces, plutot que de modifier le mesh en place.
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.update_faces(mesh.unique_faces())

    dominant_shell_fragments_dropped = 0
    # Tente d'ISOLER LA COQUE DOMINANTE (voir `_extract_dominant_shell`
    # juste en dessous) DES CE STADE, meme si `mesh.is_watertight` rapporte
    # encore True - le critere de declenchement n'est pas "le mesh est
    # non-watertight", c'est "une coque domine tres largement le reste"
    # (seuil de volume interne a la fonction, garde-fou explicite contre
    # toute perte de geometrie legitime multi-corps). Mesure sur
    # res/test_models/turkish.stl : 16 fragments minuscules (76 faces sur
    # 10180, 0.31% du volume, normales inversees - probablement des
    # artefacts d'export CSG) accoles au corps principal par des jonctions
    # non-manifold rendaient le mesh COMBINE non-watertight - mais un cas
    # structurellement similaire (fragments qui se touchent/penetrent sans
    # creer d'arete partagee non-manifold detectable par le simple
    # comptage d'aretes de trimesh) peut rester `is_watertight=True` tout
    # en beneficiant du meme nettoyage. Appele avant toute reparation
    # standard : `_extract_dominant_shell` opere sur les coques deja
    # INDIVIDUELLEMENT watertight (`mesh.split(only_watertight=False)`),
    # jamais besoin de reparer quoi que ce soit pour les detecter.
    dominant, dropped_count = _extract_dominant_shell(mesh)
    if dominant is not None:
        mesh = dominant
        dominant_shell_fragments_dropped = dropped_count

    if not mesh.is_watertight:
        trimesh.repair.fill_holes(mesh)
        trimesh.repair.fix_normals(mesh)
        trimesh.repair.fix_winding(mesh)

    if not mesh.is_watertight:
        # Dernier recours : `_repair_with_pymeshfix` (MeshFix) - corrige
        # generiquement les auto-intersections ET les trous non-manifold
        # qu'aucun des deux mecanismes ci-dessus n'a resolus. Meme
        # discipline opt-in que le reste (dependance optionnelle, jamais
        # d'exception si absente ou en echec, mesh conserve tel quel dans
        # ce cas) - mais accepte son resultat SEULEMENT s'il ne s'effondre
        # pas radicalement par rapport a l'entree (garde-fou explicite,
        # motive par le meme incident turkish.stl ci-dessus : un resultat
        # watertight n'est pas forcement un resultat FIDELE).
        pymeshfix_repaired = _repair_with_pymeshfix(mesh)
        if (
            pymeshfix_repaired is not None
            and pymeshfix_repaired.is_watertight
            and len(pymeshfix_repaired.faces) >= 0.5 * len(mesh.faces)
        ):
            mesh = pymeshfix_repaired

    # Les auto-intersections ne sont PAS reparees ici : leur seul juge fiable
    # est TetGen lui-meme (voir _tetrahedralize). Un compteur geometrique est
    # renseigne a titre informatif pour le rapport, jamais comme declencheur.
    self_intersecting_before = _count_self_intersecting_faces(mesh)

    vertices_merged = max(0, vertex_count_before - mesh.vertices.shape[0])
    faces_removed = max(0, face_count_before - mesh.faces.shape[0])

    report = IngestReport(
        was_watertight_initially=was_watertight_initially,
        is_watertight_after_repair=bool(mesh.is_watertight),
        faces_removed_degenerate=faces_removed,
        vertices_merged=vertices_merged,
        self_intersecting_faces_before=self_intersecting_before,
        tetrahedra_generated=0,  # rempli plus tard, apres tetraedrisation
        source_face_count=mesh.faces.shape[0],
        source_vertex_count=mesh.vertices.shape[0],
        dominant_shell_fragments_dropped=dominant_shell_fragments_dropped,
    )
    return mesh, report


def _extract_dominant_shell(mesh: trimesh.Trimesh) -> tuple[Optional[trimesh.Trimesh], int]:
    """
    Isole la plus grosse coque INDIVIDUELLEMENT watertight d'un mesh
    combine non-watertight, si elle domine tres largement les autres
    composantes (seuil `_DOMINANT_SHELL_FACE_RATIO`/`_DOMINANT_SHELL_VOLUME_RATIO`).

    Motive par un cas reel mesure sur `res/test_models/turkish.stl` : le
    mesh combine (10180 faces) n'est pas watertight a cause de 16
    fragments minuscules (76 faces, 0.31% du volume, normales inversees -
    probablement des artefacts d'export CSG) accoles au corps principal
    par des jonctions non-manifold - mais `mesh.split(only_watertight=False)`
    revele que le corps principal (10104 faces) EST DEJA individuellement
    watertight, avec les bonnes dimensions. Les fragments ne sont pas des
    "trous a combler", ce sont des morceaux distincts et negligeables.

    Retourne (None, 0) si aucune composante ne domine assez clairement
    (garde-fou : jamais applique a un mesh legitimement multi-corps ou a
    des composantes de taille comparable - dans ce cas, l'appelant retombe
    sur le repli suivant de la cascade, jamais un silence).
    """
    try:
        components = mesh.split(only_watertight=False)
    except Exception:
        return None, 0

    if len(components) <= 1:
        return None, 0

    watertight_components = [c for c in components if c.is_watertight]
    if not watertight_components:
        return None, 0

    dominant = max(watertight_components, key=lambda c: abs(c.volume))
    total_volume = sum(abs(c.volume) for c in components)

    volume_ratio = abs(dominant.volume) / max(total_volume, 1e-9)

    if volume_ratio < _DOMINANT_SHELL_VOLUME_RATIO:
        return None, 0

    # R-013 : une vraie piece posee a cote de la grande (volume POSITIF non negligeable) ne doit jamais etre
    # jetee en silence. Les fragments de turkish.stl sont soit a volume signe negatif (normales inversees),
    # soit de l'ordre de 1e-6 du volume total ; une piece imprimable n'est ni l'un ni l'autre.
    for c in components:
        if c is not dominant and float(c.volume) > 0 and abs(c.volume) / max(total_volume, 1e-9) >= _KEPT_PIECE_VOLUME_RATIO:
            return None, 0

    dropped_count = len(components) - 1
    return dominant, dropped_count


# Seuil de dominance pour `_extract_dominant_shell` : la coque retenue doit
# representer au moins 95% du VOLUME total pour que les fragments restants
# soient consideres negligeables plutot que de la geometrie legitime a
# preserver. Le volume est le signal pertinent, PAS le nombre de faces
# (une petite coque peut avoir proportionnellement autant de triangles
# qu'une grande - mesure directement sur un cas synthetique pendant la
# conception de ce fix : un fragment 1000x plus petit en volume peut porter
# le MEME nombre de faces qu'un cube primitif de reference, rendant un
# ratio de faces non pertinent comme critere de dominance). Mesure sur
# turkish.stl : la coque dominante y represente 99.69% du volume total,
# tres au-dela de ce seuil (marge large, pas un seuil calibre au plus
# juste).
_DOMINANT_SHELL_VOLUME_RATIO = 0.95


# Une composante a volume positif d'au moins cette part du volume total est une piece, pas un fragment.
_KEPT_PIECE_VOLUME_RATIO = 1e-3
