"""
pipeline/step1_ingest.py

Etape 1 du pipeline : ingestion du mesh source (STL/3MF) et preparation de
toutes les structures dont les etapes suivantes ont besoin.

Cette etape est la seule du pipeline qui touche a l'I/O disque cote geometrie.
Elle ne depend d'aucune autre etape - c'est le point d'entree du pipeline et
le fichier le plus simple a tester en isolation (aucune donnee amont requise,
juste un chemin de fichier).

Sortie : (MeshData, TetMesh, IngestReport)
"""

from __future__ import annotations

import os
import pickle
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Optional

import math

import numpy as np
import scipy.sparse as sp
import scipy.spatial
import trimesh
import trimesh.transformations as trimesh_tf
from trimesh import geometry as trimesh_geometry
from trimesh import grouping as trimesh_grouping
from trimesh.creation import triangulate_polygon
from trimesh.path import polygons as trimesh_polygons

from core import _cpp
from core.mesh_types import JunctionZone, MeshData, TetMesh
from core.spatial import KDTree, RayIntersector
from pipeline._self_intersection_resolve import resolve_self_intersections


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
    # Fragments negligeables (coques individuellement watertight mais
    # minuscules par rapport a la coque dominante) ecartes pour rendre le
    # mesh watertight - donnee EXPLICITE, jamais une perte silencieuse de
    # geometrie. 0 dans l'immense majorite des cas (mesh deja propre) -
    # mesure sur res/test_models/turkish.stl : 16 fragments (76 faces sur
    # 10180, 0.31% du volume) accoles au corps principal par des jonctions
    # non-manifold, faisant echouer la reparation standard ET pymeshfix
    # (qui reduisait le mesh entier a 10 faces - voir CLAUDE.md).
    dominant_shell_fragments_dropped: int = 0


class MeshIngestError(Exception):
    """Levee quand le mesh source ne peut pas etre charge ou exploite."""


# Tolerance du test de separation par plan (_triangles_cross) : en dessous,
# un sommet est considere SUR le plan de l'autre triangle, pas d'un cote.
# Deux faces adjacentes coplanaires se toucheraient sinon a tort.
_COPLANAR_EPS = 1e-9


# ---------------------------------------------------------------------------
# Etape 1 - Chargement du fichier brut
# ---------------------------------------------------------------------------

def _load_raw_mesh(mesh_path: str, allow_multi_mesh: bool) -> trimesh.Trimesh:
    loaded = trimesh.load(mesh_path, force="mesh" if not allow_multi_mesh else None)

    if isinstance(loaded, trimesh.Scene):
        if not allow_multi_mesh:
            raise MeshIngestError(
                f"'{mesh_path}' contient plusieurs objets (Scene) mais "
                "allow_multi_mesh=False. Fusionnez le fichier ou activez ce parametre."
            )
        geometries = list(loaded.geometry.values())
        if len(geometries) == 0:
            raise MeshIngestError(f"'{mesh_path}' ne contient aucune geometrie.")
        loaded = trimesh.util.concatenate(geometries)

    if not isinstance(loaded, trimesh.Trimesh):
        raise MeshIngestError(
            f"'{mesh_path}' n'a pas produit un Trimesh exploitable "
            f"(type obtenu: {type(loaded)})."
        )

    if loaded.faces.shape[0] == 0:
        raise MeshIngestError(f"'{mesh_path}' est charge mais ne contient aucune face.")

    return loaded


# ---------------------------------------------------------------------------
# Etape 2bis - Auto-intersections
#
# TetGen REFUSE ("The input surface mesh contain self-intersections") tout
# maillage dont des triangles se traversent, meme parfaitement watertight,
# a winding coherent, d'un seul tenant et sans face degeneree - c'est
# exactement le cas de res/test_models/death.stl (mesure : 5024 faces
# fautives sur 12662, soit 39.7%, reparties sur tout le modele).
#
# Aucune des reparations deja en place (merge_vertices, nondegenerate_faces,
# unique_faces, fill_holes, fix_normals, fix_winding) ne traite ce defaut :
# elles corrigent la TOPOLOGIE, pas la geometrie qui s'interpenetre. Une
# union booleenne manifold3d ne le traite pas non plus (mesure : rend le
# maillage a l'identique, meme nombre de faces, meme volume).
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Reparation LOCALE des auto-intersections (retrait minimal + rebouchage 3D)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Etape 2 - Nettoyage et reparation
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Etape 4 - Structures d'acceleration (BVH + KDTree)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Etape 5 - Assemblage du MeshData
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Etape 6 - Tetraedrisation volumique (tetgen)
# ---------------------------------------------------------------------------

_INTERIOR_POINT_SEED = 0  # deterministe - reproductibilite des seeds de region
_INTERIOR_POINT_INITIAL_SAMPLES = 64
_INTERIOR_POINT_MAX_SAMPLES = 4096  # plafond - au-dela, echec explicite


def _find_guaranteed_interior_point(
    shell: trimesh.Trimesh, fallback_point: np.ndarray, context_label: str
) -> np.ndarray:
    """
    Retourne un point GARANTI a l'interieur de `shell` (verifie via
    `shell.contains()`), pour toute coque watertight arbitraire (convexe ou
    non). Essaie d'abord `fallback_point` (typiquement un centroide deja
    calcule par l'appelant, chemin rapide - reste correct sur la tres grande
    majorite des coques, notamment tout ce qui est convexe ou proche de
    convexe). Si `fallback_point` echoue le test d'inclusion (mesure
    concrete : forme non convexe ou coque a plusieurs composantes
    disjointes, cf. CLAUDE.md limite #42 - un centroide de coque "rest" hors
    du solide reel sur un modele en T, qui casse SILENCIEUSEMENT
    `TetGen.add_region` : la region dont le seed est hors-coque n'est plus
    localisable dans la triangulation, et son volume cible est ignore au
    profit de l'autre region, appliquee alors a tout le maillage), retombe
    sur un tirage par rejet (`trimesh.sample.volume_mesh`, seed FIXE pour la
    reproductibilite - meme discipline que le determinisme deja exige du tri
    DFS spatial de step5), avec une escalade deterministe du nombre
    d'essais si le tirage initial ne produit aucun point accepte (mesure :
    `volume_mesh` peut renvoyer 0 point sur une coque fine par rapport a sa
    bounding box).

    Retourne TOUJOURS un seul point INDIVIDUELLEMENT valide, jamais une
    moyenne de plusieurs points acceptes - piege mesure concretement : la
    moyenne de plusieurs points valides peut elle-meme retomber hors d'une
    coque non convexe, meme si chacun des points moyennes est
    individuellement a l'interieur.

    Args:
        shell: coque watertight (suppose deja etanche cote appelant, meme
            convention que le reste de ce fichier a ce stade du pipeline).
        fallback_point: point candidat a essayer en premier (ex: centroide).
        context_label: texte court identifiant l'appelant, uniquement pour
            un message d'erreur lisible (ex: "coque 'rest'").

    Returns:
        Point (3,) verifie a l'interieur de `shell`.

    Raises:
        MeshIngestError: coque vide/degeneree (aucune face, bounding box
            indefinie - sans ce garde-fou, `trimesh.sample.volume_mesh`
            plante avec un `TypeError` cryptique, mesure concretement) avant
            meme d'essayer l'echantillonnage, ou aucun point trouve apres
            escalade jusqu'a `_INTERIOR_POINT_MAX_SAMPLES` essais (coque
            quasi degeneree, ratio volume/bounding-box extremement faible).
    """
    if shell.faces.shape[0] == 0 or shell.bounds is None:
        raise MeshIngestError(
            f"{context_label} : coque vide ou degeneree, impossible d'y "
            "chercher un point interieur."
        )

    if shell.contains([fallback_point])[0]:
        return np.asarray(fallback_point, dtype=np.float64)

    sample_count = _INTERIOR_POINT_INITIAL_SAMPLES
    while sample_count <= _INTERIOR_POINT_MAX_SAMPLES:
        points = trimesh.sample.volume_mesh(shell, count=sample_count, seed=_INTERIOR_POINT_SEED)
        if len(points) > 0:
            return np.asarray(points[0], dtype=np.float64)
        sample_count *= 8

    raise MeshIngestError(
        f"{context_label} : impossible de trouver un point garanti a "
        f"l'interieur apres {_INTERIOR_POINT_MAX_SAMPLES} essais "
        "d'echantillonnage volumique (coque quasi degeneree - ratio "
        "volume/bounding-box extremement faible)."
    )


def _find_internal_cavities(mesh: trimesh.Trimesh) -> list[np.ndarray]:
    """
    Detecte les cavites internes du mesh (ex: boite creuse) et retourne un
    point representatif a l'interieur de chacune, a passer a
    tetgen.TetGen.add_hole(). Sans ça, tetgen ne sait pas qu'un volume
    interne doit rester vide et peut generer des tetraedres qui traversent
    la cavite (mesure : jusqu'a ~6% des tetraedres sur une boite creuse
    simple, centroides verifies hors du solide reel).

    Methode : separe le mesh en coques connexes etanches
    (`mesh.split(only_watertight=True)`). Le signe du volume d'une coque
    (propriete trimesh, deja orientee) distingue frontiere solide (positif)
    de cavite (negatif - ses normales, heritees d'une operation booleenne
    de difference, pointent vers l'interieur du vide). Le centroide de
    chaque coque-cavite sert de point representatif ; si ce centroide tombe
    hors de sa propre coque (cavite non convexe, ex: forme en tore), repli
    via `_find_guaranteed_interior_point` (rejection sampling deterministe)
    plutot qu'un echec immediat.
    """
    shells = mesh.split(only_watertight=True)
    if len(shells) <= 1:
        return []

    hole_points = []
    for shell in shells:
        if shell.volume >= 0:
            continue  # coque solide (frontiere exterieure ou paroi normale), pas une cavite

        centroid = shell.vertices.mean(axis=0)
        hole_points.append(
            _find_guaranteed_interior_point(shell, centroid, "cavite interne (volume signe negatif)")
        )

    return hole_points


def _shared_cut_cap(
    open_mesh: trimesh.Trimesh, plane_origin: np.ndarray, plane_normal: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Triangule UNE SEULE FOIS le contour de bord cree par une coupe
    `slice_plane(..., cap=False)`, pour etre reutilisee identique des deux
    cotes de la coupe - voir `_split_for_platform_refinement` pour pourquoi
    deux triangulations independantes (une par `slice_plane(cap=True)`)
    divergent sur un contour non convexe.

    Reproduit la logique interne du parametre `cap=True` de
    `trimesh.intersections.slice_mesh_plane` (dedup des sommets, extraction
    des aretes de bord sur le plan, `triangulate_polygon`), mais appliquee
    UNE FOIS a la place de deux fois independamment.

    Retourne (cap_vertices (K,3), cap_faces (T,3)) dans le repere du monde,
    triangulation orientee de facon a ce que la normale des triangles suive
    `plane_normal` (meme convention que `slice_plane(cap=True)` : le cap
    ferme le mesh du cote NEGATIF de `plane_normal`, cote qui a ete tranche).
    """
    vertices = open_mesh.vertices
    faces = open_mesh.faces

    # slice_plane(cap=False) laisse des sommets numeriquement dupliques a
    # chaque intersection arete/plan, sans les fusionner - meme etape de
    # dedup que le code interne de slice_mesh_plane avant de calculer un cap.
    unique, inverse = trimesh_grouping.unique_rows(vertices)
    vertices = vertices[unique]
    faces = inverse[faces]
    faces = faces[(faces[:, :1] != faces[:, 1:]).all(axis=1)]

    dots = np.dot(vertices - plane_origin, plane_normal)
    on_plane = np.abs(dots) < 1e-8

    edges = trimesh_geometry.faces_to_edges(faces)
    edges.sort(axis=1)
    edges = edges[on_plane[edges].all(axis=1)]
    edges = edges[edges[:, 0] != edges[:, 1]]

    # arete de bord du contour de coupe = apparait exactement une fois
    # (une arete interieure a une face qui serait par coincidence entierement
    # sur le plan apparaitrait deux fois, une par face adjacente).
    unique_edge_idx = trimesh_grouping.group_rows(edges, require_count=1)
    edges_boundary = edges[unique_edge_idx]

    to_2D = trimesh_geometry.plane_transform(origin=plane_origin, normal=-plane_normal)
    to_3D = np.linalg.inv(to_2D)
    vertices_2D = trimesh_tf.transform_points(vertices, to_2D)

    all_v2 = []
    all_f = []
    offset = 0
    for poly in trimesh_polygons.edges_to_polygons(edges_boundary, vertices_2D[:, :2]):
        v2, f2 = triangulate_polygon(poly, force_vertices=True)
        all_v2.append(v2)
        all_f.append(f2 + offset)
        offset += len(v2)

    if not all_v2:
        return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64)

    v2_all = np.vstack(all_v2)
    f_all = np.vstack(all_f)
    v3_all = trimesh_tf.transform_points(trimesh.util.stack_3D(v2_all), to_3D)

    return np.asarray(v3_all, dtype=np.float64), np.asarray(f_all, dtype=np.int64)


def _split_for_platform_refinement(
    mesh: trimesh.Trimesh,
    near_platform_margin: float,
    platform_z: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Decoupe `mesh` en un PLC (piecewise linear complex) combine, deux
    regions de part et d'autre de `z = platform_z + near_platform_margin`
    ("slab" pres du plateau, "rest" au-dela), pour alimenter
    `TetGen.add_region` avec un point-seed par region - mecanisme retenu
    pour le raffinement local (voir CLAUDE.md pour le mecanisme bgmesh/metric
    ecarte : silencieusement ignore par le binding tetgen installe des qu'un
    `switches=` non vide est passe, et provoque un segfault natif
    reproductible si on force `-m` malgre tout).

    Le contour de coupe est triangule UNE SEULE FOIS (`_shared_cut_cap`) et
    n'apparait qu'UNE SEULE fois dans le PLC combine, comme facette INTERNE
    partagee entre les deux regions - PAS deux fois via deux appels
    independants a `slice_plane(cap=True)`. Correction d'un bug reel mesure
    sur `res/test_models/death.stl` (benchy avec disque de support rond) :
    deux appels independants a `slice_plane(cap=True)` (un pour "slab", un
    pour "rest") trianguleraient le MEME contour (memes sommets XY, ecart
    mesure = 0.0) mais de facon INDEPENDANTE - `triangulate_polygon` choisit
    ses diagonales differemment a chaque appel des que le contour est non
    convexe. Mesure directe : `l_shape_overhang.stl` (8 sommets au contour)
    n'a que 2 triangles divergents sur 6, sous le seuil de tolerance de
    TetGen ; `death.stl` (47 sommets, disque de support + jambes du benchy)
    en a 12 sur 47, ce que TetGen refuse net ("The input surface mesh
    contain self-intersections"). Le repli en cascade de `_tetrahedralize`
    (limite #43) masquait l'echec en desactivant purement et simplement le
    raffinement plateau - reintroduisant l'erreur de discretisation FEM de
    la limite #34 exactement dans la zone qui devrait rester plate (le
    disque de support), d'ou les couches fragmentees en « escalier »
    observees visuellement. Corrige : triangulation unique + une seule copie
    dans le PLC final, verifie sans aucun warning TetGen sur les deux
    modeles (`death.stl` : 112698 tetraedres avec raffinement plateau reel,
    contre 107037 en mode degrade avant ce fix ; `l_shape_overhang.stl` :
    1677 tetraedres, identique a la mesure deja documentee - non-regression
    confirmee).

    Retourne (combined_vertices (N,3), combined_faces (M,3), slab_seed (3,),
    rest_seed (3,)). Le seed de chaque region part du centroide de sa coque
    ouverte (coordonnee z explicitement repoussee vers l'interieur de sa
    region pour ne jamais tomber exactement sur la frontiere de coupe),
    VERIFIE via `shell.contains()` (`_find_guaranteed_interior_point`) sur
    une coque watertight temporaire (coque ouverte + cap oriente vers
    l'interieur) construite UNIQUEMENT pour ce test d'inclusion - jamais
    pour le PLC final, qui lui ne contient qu'une seule copie du cap. Le
    test d'inclusion echoue et retombe sur un tirage par rejet (deterministe)
    sur forme non convexe (CLAUDE.md limite #42 : le centroide de la coque
    "rest" tombe hors du solide sur le modele en T reel,
    `l-triple-axis-overhanged-shape.stl`, et un seed hors-coque casse
    silencieusement `TetGen.add_region`).
    """
    cut_z = platform_z + near_platform_margin
    plane_origin = np.array([0.0, 0.0, cut_z])
    plane_normal_up = np.array([0.0, 0.0, 1.0])

    slab_open = mesh.slice_plane(plane_origin=plane_origin, plane_normal=-plane_normal_up, cap=False)
    rest_open = mesh.slice_plane(plane_origin=plane_origin, plane_normal=plane_normal_up, cap=False)

    cap_vertices, cap_faces = _shared_cut_cap(rest_open, plane_origin, plane_normal_up)

    n_slab = len(slab_open.vertices)
    n_rest = len(rest_open.vertices)

    combined_vertices_dup = np.vstack([slab_open.vertices, rest_open.vertices, cap_vertices])
    combined_faces_dup = np.vstack([
        slab_open.faces,
        rest_open.faces + n_slab,
        cap_faces + n_slab + n_rest,  # une seule copie - facette interne partagee
    ])

    # slab_open/rest_open portent chacun leur PROPRE copie numerique des
    # sommets du contour de coupe (memes coordonnees que cap_vertices,
    # indices differents) - sans cette fusion globale, le cap resterait
    # topologiquement deconnecte des bords ouverts malgre des coordonnees
    # identiques.
    unique, inverse = trimesh_grouping.unique_rows(combined_vertices_dup)
    combined_vertices = combined_vertices_dup[unique]
    combined_faces = inverse[combined_faces_dup]
    combined_faces = combined_faces[(combined_faces[:, :1] != combined_faces[:, 1:]).all(axis=1)]

    # Coques watertight temporaires (coque ouverte + cap oriente vers
    # l'interieur de CETTE coque), uniquement pour le test d'inclusion des
    # seeds - jamais utilisees pour le PLC final envoye a tetgen.
    slab_for_seed = trimesh.Trimesh(
        vertices=np.vstack([slab_open.vertices, cap_vertices]),
        faces=np.vstack([slab_open.faces, cap_faces[:, ::-1] + n_slab]),
        process=True,
    )
    rest_for_seed = trimesh.Trimesh(
        vertices=np.vstack([rest_open.vertices, cap_vertices]),
        faces=np.vstack([rest_open.faces, cap_faces + n_rest]),
        process=True,
    )

    slab_seed_candidate = slab_open.vertices.mean(axis=0)
    slab_seed_candidate[2] = min(slab_seed_candidate[2], platform_z + near_platform_margin * 0.5)
    rest_seed_candidate = rest_open.vertices.mean(axis=0)
    rest_seed_candidate[2] = max(rest_seed_candidate[2], platform_z + near_platform_margin * 1.5)

    slab_seed = _find_guaranteed_interior_point(
        slab_for_seed, slab_seed_candidate, "coque 'slab' (pres du plateau)"
    )
    rest_seed = _find_guaranteed_interior_point(
        rest_for_seed, rest_seed_candidate, "coque 'rest' (reste du modele)"
    )

    return (
        np.asarray(combined_vertices, dtype=np.float64),
        np.asarray(combined_faces, dtype=np.int64),
        np.asarray(slab_seed, dtype=np.float64),
        np.asarray(rest_seed, dtype=np.float64),
    )


def _merged_junction_seed(mesh: trimesh.Trimesh, junction_zones: list[JunctionZone]) -> np.ndarray:
    """
    Fusionne toutes les `JunctionZone` detectees en UN SEUL seed interne
    (centroide pondere par rayon), pour alimenter une UNIQUE region
    `TetGen.add_region` supplementaire plutot que N regions separees.

    Necessaire car `TetGen.add_region` (binding 0.8.4) ne se generalise PAS
    correctement a N regions actives simultanement - mesure directe
    (scripts/measure_add_region_n_zones.py, session de la limite #48) : sur
    4 regions (1 "rest" + 3 zones fines), AUCUNE des 3 zones fines
    n'atteint le facteur de densite attendu face a "rest" (l'une d'elles,
    junction_y, se retrouve meme MOINS dense que la region censee etre
    grossiere) - une generalisation du bug deja documente a 2 regions
    (limite #42), pas juste une variante mineure. Le fallback a 2 regions
    (1 "rest" + 1 zone fine fusionnee) fonctionne correctement, mesure sur
    le meme modele (scripts/measure_add_region_merged_zone.py) : facteur
    ~48x de densite locale pour la zone fusionnee contre la zone "rest" et
    un point de controle hors de toute zone fine, aucune fuite du volume
    fin hors de son rayon d'influence.

    Consequence acceptee : deux jonctions multi-region tres eloignees
    l'une de l'autre partagent le meme volume cible fin, meme si l'une
    seule en aurait besoin - imprecision mineure comparee au fait de payer
    le cout d'un maillage fin partout (`resolution` globalement raffine),
    et strictement meilleure que le comportement actuel (aucun raffinement
    du tout aux jonctions, limite #47).

    Le point retourne n'est PAS garanti a l'interieur du mesh (une moyenne
    de centroides individuellement valides peut retomber hors d'un solide
    non convexe, piege deja documente pour `_find_guaranteed_interior_point`)
    - l'appelant doit le passer par cette fonction avant de l'utiliser comme
    seed `add_region`.
    """
    centroids = np.array([zone.centroid for zone in junction_zones])
    radii = np.array([zone.radius for zone in junction_zones])
    weights = radii / radii.sum()
    return np.average(centroids, axis=0, weights=weights)


def _find_rest_seed_outside_junction_zones(
    mesh: trimesh.Trimesh, junction_zones: list[JunctionZone]
) -> np.ndarray:
    """
    Seed "rest" garanti a la fois A L'INTERIEUR du mesh ET LE PLUS ELOIGNE
    POSSIBLE de toute `JunctionZone` - bug reel mesure en ecrivant les tests
    de cette fonction (pas suppose) : le centroide naif du mesh (`mesh.
    vertices.mean(axis=0)`) peut tomber A L'INTERIEUR du rayon d'une zone de
    jonction (mesure sur le tripode synthetique : centroide a 3.42mm du
    centre d'une zone de rayon 6.0mm), auquel cas TetGen applique le volume
    FIN de la region jonction a la region "rest" tout entiere - explosion
    de densite mesuree sur tout le maillage (905 -> 44985 tetraedres sur le
    tripode, loin de toute jonction : 34 -> 2302 tetraedres locaux), pas
    seulement pres de la jonction ciblee.

    Essaie d'abord le centroide naif (chemin rapide, correct des que le
    mesh est assez grand par rapport aux zones de jonction). Sinon, cherche
    parmi TOUS les sommets du mesh (deterministe, pas un tirage aleatoire)
    celui qui MAXIMISE la distance minimale a toutes les zones - garantit de
    trouver le meilleur candidat possible plutot que d'esperer un tirage
    aleatoire chanceux, cas mesure directement sur res/test_models/death.stl
    (motif festonne tres etendu, une zone fusionnee de rayon 47mm couvre une
    portion importante du volume) : un rejection-sampling par tirage
    aleatoire y echouait apres 4096 essais malgre l'existence de sommets
    valides (les coins/extremites du mesh, hors de la zone dense mais rares
    relativement au volume total echantillonne par `trimesh.sample.
    volume_mesh`).
    """
    candidate = mesh.vertices.mean(axis=0)

    def _min_dist_to_zones(point: np.ndarray) -> float:
        return min(
            (np.linalg.norm(point - zone.centroid) - zone.radius for zone in junction_zones),
            default=np.inf,
        )

    if mesh.contains([candidate])[0] and _min_dist_to_zones(candidate) > 0.0:
        return np.asarray(candidate, dtype=np.float64)

    # Les sommets du mesh sont sur la SURFACE, jamais utilisables tels
    # quels comme seed volumique TetGen - on cherche le sommet le plus
    # eloigne de toutes les zones, puis on le POUSSE vers le centroide
    # global du mesh (garanti interieur pour un solide simplement connexe,
    # meme discipline que le "candidate" naif ci-dessus) jusqu'a ce que le
    # point pousse soit lui-meme a l'interieur du mesh.
    vertices = mesh.vertices
    min_dists = np.array([
        min((np.linalg.norm(v - zone.centroid) - zone.radius for zone in junction_zones), default=np.inf)
        for v in vertices
    ])
    best_idx = int(np.argmax(min_dists))

    if min_dists[best_idx] <= 0.0:
        raise MeshIngestError(
            "Impossible de trouver un seed 'rest' hors de toutes les zones de "
            "jonction : meme le sommet du mesh le plus eloigne "
            f"(distance={min_dists[best_idx]:.2f}mm au bord de la zone la plus "
            "proche) reste a l'interieur d'une zone - les zones de jonction "
            "couvrent la quasi-totalite du volume sur ce modele."
        )

    surface_point = vertices[best_idx]
    for fraction in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5):
        pushed = surface_point + fraction * (candidate - surface_point)
        if mesh.contains([pushed])[0] and _min_dist_to_zones(pushed) > 0.0:
            return np.asarray(pushed, dtype=np.float64)

    raise MeshIngestError(
        "Impossible de pousser le seed 'rest' (sommet le plus eloigne des "
        "zones de jonction) a l'interieur du mesh sans re-tomber dans une "
        "zone de jonction - geometrie non convexe non couverte par cette "
        "heuristique."
    )


# Sous-chaines de messages RuntimeError leves par TetGen (relayes via
# `_tetrahedralize_once_isolated`) que la cascade de reparation reactive
# (limite #43) doit traiter comme RECUPERABLES - une reparation
# (pymeshfix) ou un repli (sans raffinement plateau) a une vraie chance de
# les resoudre, contrairement a une erreur structurelle non liee au mesh.
# "self-intersection" couvrait le refus explicite initial de TetGen ; le
# crash `recoversubfaces` est une VARIANTE differente mais tout aussi
# fatale du meme probleme sous-jacent (deja documentee dans CLAUDE.md
# limite #44 comme telle, jamais reliee au code jusqu'a ce fix) - mesure
# concretement sur res/test_models/turkish.stl : la cascade ne se
# declenchait PAS sur ce message, remontant l'erreur brute a l'appelant
# sans jamais essayer la reparation pymeshfix, alors qu'elle etait deja
# tentee (et fonctionnelle) pour "self-intersection".
_RECOVERABLE_TETGEN_ERROR_SUBSTRINGS = (
    "self-intersection",
    "recoversubfaces",
    # Meme classe d'echec interne TetGen que "recoversubfaces" (limite
    # #44), variante distincte rencontree en traitant turkish.stl avec
    # raffinement plateau actif (`near_platform_resolution` non None) : le
    # PLC combine plateau/reste (`_split_for_platform_refinement`) peut
    # faire echouer une AUTRE phase interne de recouvrement de frontiere
    # que la reconstruction de sous-facettes. Meme repli deja en place
    # (cascade de reparation, puis "sans raffinement plateau" en dernier
    # recours) s'applique sans changement.
    "split_segment",
)


def _is_recoverable_tetgen_error(exc: RuntimeError) -> bool:
    message = str(exc).lower()
    return any(substring in message for substring in _RECOVERABLE_TETGEN_ERROR_SUBSTRINGS)


_TETGEN_BRACKETED_INDICES = re.compile(r"\[([\d,\s]+)\]")


def _parse_tetgen_conflict_vertices(tetgen_output: str) -> set[int]:
    """
    Extrait les indices de sommets nommes explicitement par TetGen dans un
    message d'avertissement/erreur ("Two segments exactly intersect. 1st
    seg [869,836]...", "A vertex lies on a facet. vertex : [4340]...",
    "Two facets exactly intersect. 1st facet triangle: [...]...").

    Approche volontairement generique (pas un parseur par type de message) :
    TOUT groupe `[a,b,c]` d'entiers separes par virgule dans la sortie est
    un indice de sommet dans TOUS les formats de message observes (segment,
    facette, sommet isole) - `tag(-1)` voisin utilise des parentheses, pas
    des crochets, donc aucune confusion mesuree. Un parseur specifique par
    motif de phrase serait plus fragile a la moindre variation de
    formatage entre versions de TetGen.
    """
    indices: set[int] = set()
    for match in _TETGEN_BRACKETED_INDICES.finditer(tetgen_output):
        for token in match.group(1).split(","):
            token = token.strip()
            if token.isdigit():
                indices.add(int(token))
    return indices


_GEOMETRIC_RESOLUTION_MAX_TETGEN_PROBES = 5


def _probe_tetgen_conflicts(mesh: trimesh.Trimesh) -> str:
    """
    Sonde TetGen en mode DETECTION (`switches='d'`, leger - pas de
    maillage de qualite) pour recuperer sa sortie textuelle complete,
    porteuse des indices de sommets/segments/facettes en conflit.

    Mesure directe sur turkish.stl : la sortie de la VRAIE tetraedrisation
    (`_tetrahedralize_once_isolated`, mode qualite `pq1.2a...`) s'arrete
    net juste apres "Recovering boundaries...", AVANT les avertissements
    eux-memes, meme redirigee vers de vrais fichiers - tampon de la libc du
    binding jamais vide sur ce chemin de code precis. Le mode DETECTION,
    plus leger, a ete mesure comme flushant correctement dans tous les cas
    testes - d'ou ce sondage SEPARE, dans son propre sous-processus isole
    (`pipeline/_tetgen_detect_worker.py`), uniquement pour obtenir le
    texte, jamais pour produire un maillage tetraedrique reel.

    Best-effort, ne leve jamais : retourne une chaine vide si le
    sous-processus crashe ou si les fichiers de sortie sont illisibles -
    l'appelant traite deja ce cas comme "aucun sommet nommable".
    """
    try:
        with tempfile.TemporaryDirectory(prefix="apex_tetgen_probe_") as tmp_dir:
            input_path = os.path.join(tmp_dir, "input.pkl")
            stdout_path = os.path.join(tmp_dir, "stdout.txt")
            stderr_path = os.path.join(tmp_dir, "stderr.txt")

            with open(input_path, "wb") as f:
                pickle.dump(
                    {
                        "vertices": np.asarray(mesh.vertices, dtype=np.float64),
                        "faces": np.asarray(mesh.faces, dtype=np.int64),
                    },
                    f,
                )

            with open(stdout_path, "wb") as stdout_file, open(stderr_path, "wb") as stderr_file:
                subprocess.run(
                    [sys.executable, "-m", "pipeline._tetgen_detect_worker", input_path],
                    cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    stdout=stdout_file,
                    stderr=stderr_file,
                )

            with open(stdout_path, "rb") as f:
                text = f.read().decode("utf-8", errors="replace")
            with open(stderr_path, "rb") as f:
                text += f.read().decode("utf-8", errors="replace")
            return text
    except OSError:
        return ""


def _finish_geometric_resolution(
    resolved: trimesh.Trimesh,
    resolution: float,
) -> trimesh.Trimesh:
    """
    Affine iterativement le resultat de `resolve_self_intersections` en
    utilisant TetGen LUI-MEME comme oracle (principe deja acte, limite
    #43 : "TetGen est le seul juge fiable de ce qui le bloque").

    A chaque iteration : tente une VRAIE tetraedrisation (subprocess isole,
    `_tetrahedralize_once_isolated`) ; si TetGen la refuse pour une raison
    RECUPERABLE et nomme des sommets en conflit dans sa sortie textuelle
    (via `_probe_tetgen_conflicts`), cible une reparation locale
    (`_repair_local_self_intersections`, `extra_seed_faces` = les faces
    portant ces sommets) sur EXACTEMENT ce que TetGen vient de refuser, au
    lieu d'un repli generique. Mesure sur turkish.stl : 8 auto-
    intersections locales (issues du soudage de trous residuels de
    `resolve_self_intersections`) rendaient TetGen incapable d'aboutir
    malgre un SAT=0 global (contact EXACT entre deux segments de faces
    partageant un sommet, donc invisible au test SAT qui exclut ces
    paires) ; ce mecanisme reduit le residu jusqu'a un dernier segment
    nomme, cible et repare en 1 iteration supplementaire - TetGen aboutit
    alors sur la geometrie COMPLETE (330306 tetraedres, z=[30.26,90.78],
    contre l'echec total precedent), pour un ecart de volume negligeable.

    Deliberement SANS raffinement plateau/jonction ici (jamais
    `near_platform_resolution`/`junction_refinement_resolution`), meme si
    l'appelant final en utilise un : ce raffinement construit un PLC
    COMBINE avec de nouveaux sommets injectes au plan de coupe
    (`_split_for_platform_refinement`), dont la numerotation ne correspond
    PAS a celle de `current` - un conflit nomme par TetGen sur CE PLC
    combine ne peut donc pas etre traduit fidelement en `extra_seed_faces`
    sur le mesh de surface qu'on repare ici. Bug mesure en integrant ce
    mecanisme : avec raffinement plateau actif, le sondage ne trouvait
    jamais aucun sommet exploitable (mauvais referentiel), la fonction
    rendait la main sans le moindre progres, et la tentative reelle de
    l'appelant (elle, avec raffinement) echouait alors identiquement a
    chaque fois. La question posee ici est volontairement plus etroite :
    « ce mesh de SURFACE est-il tetraedrisable EN SOI » - la compatibilite
    avec le raffinement plateau reste de la responsabilite de la cascade
    normale de `_tetrahedralize`, qui sait deja s'en degrader proprement
    (repli "sans raffinement plateau", limite #44).

    Best-effort, ne leve jamais : rend le meilleur mesh obtenu (le dernier
    `resolved` valide) des qu'un probe reussit, qu'une reparation echoue,
    qu'aucun sommet n'est nommable, ou que le budget d'iterations est
    epuise - l'appelant (`_tetrahedralize`) refera de toute facon une VRAIE
    tentative sur ce resultat dans sa cascade normale, avec son cortege
    habituel de replis (pymeshfix, sans raffinement plateau).
    """
    current = resolved
    for _probe in range(_GEOMETRIC_RESOLUTION_MAX_TETGEN_PROBES):
        try:
            _tetrahedralize_once_isolated(
                current, resolution, None, 5.0, 0.0, None, None,
            )
            return current  # deja verifie tetraedrisable tel quel (sans raffinement)
        except _TetgenCrashError:
            return current  # rien a cibler depuis un crash natif
        except RuntimeError as exc:
            if not _is_recoverable_tetgen_error(exc):
                return current
            # `exc.tetgen_output` (capture depuis la VRAIE tentative) s'est
            # mesuree peu fiable - tronquee avant les avertissements
            # nommes eux-memes sur ce chemin de code precis (voir
            # `_probe_tetgen_conflicts`). Un sondage DEDIE en mode
            # detection, plus leger, recupere le texte complet de facon
            # fiable.
            tetgen_output = _probe_tetgen_conflicts(current)
            conflict_vertices = _parse_tetgen_conflict_vertices(tetgen_output)
            if not conflict_vertices:
                return current
            faces = np.asarray(current.faces)
            seed_faces = [
                i for i, face in enumerate(faces)
                if len(conflict_vertices.intersection(int(v) for v in face)) >= 2
            ]
            if not seed_faces:
                return current
            repaired = _repair_local_self_intersections(
                current, extra_seed_faces=seed_faces, require_euler_preserved=False,
            )
            if repaired is None:
                return current
            current = repaired
    return current


def _tetrahedralize(
    mesh: trimesh.Trimesh,
    resolution: float,
    near_platform_resolution: Optional[float] = None,
    near_platform_margin: float = 5.0,
    platform_z: float = 0.0,
    junction_zones: Optional[list[JunctionZone]] = None,
    junction_refinement_resolution: Optional[float] = None,
) -> TetMesh:
    """
    Genere un maillage tetraedrique volumique respectant la frontiere du mesh
    de surface, via TetGen. `resolution` est converti en volume maximal
    autorise par tetraedre (parametre -a de tetgen), qui reste le plafond
    global.

    `near_platform_resolution` (opt-in, None = desactive par defaut - meme
    discipline que `direction_blend_radius=0.0` en step4) raffine LOCALEMENT
    la densite pres du plateau (z < platform_z + near_platform_margin), via
    deux regions TetGen (`_split_for_platform_refinement` + `add_region`,
    une pour la tranche pres du plateau a `near_platform_resolution`, une
    pour le reste a `resolution`) - PAS via un maillage de fond
    (`bgmesh`/`metric`), mecanisme teste et ecarte : silencieusement ignore
    par le binding tetgen installe (0.8.4) des qu'un `switches=` non vide
    est passe, et provoque un segfault natif reproductible si on force `-m`
    malgre tout (teste sous plusieurs formes, voir CLAUDE.md). Corrige une
    erreur de discretisation FEM mesuree pres du plateau (CLAUDE.md) : un
    tetraedre P1 assez grand pour s'etendre du contact plateau jusqu'a une
    zone ou la direction cible a deja beaucoup tourne ne peut representer
    qu'un seul gradient constant sur toute son etendue, produisant un
    artefact d'« escalier » visible dans l'inclinaison des couches - affine
    localement plutot que globalement pour ne pas payer le cout d'un
    maillage fin partout ou ce n'est pas necessaire (mesure sur
    l_shape_overhang.stl : amelioration angulaire quasi identique a un
    raffinement global complet, pour ~5.4x moins de tetraedres).

    Bug corrige (CLAUDE.md limite #42) : `add_region` ignorait
    silencieusement le volume de la region "rest" (donc `resolution`
    perdait tout effet hors de la tranche pres du plateau) quand le seed de
    cette region tombait hors de sa propre coque sur une geometrie non
    convexe - `_split_for_platform_refinement` garantit desormais un seed
    valide pour les deux regions (`_find_guaranteed_interior_point`).

    `junction_zones`/`junction_refinement_resolution` (opt-in, None =
    desactive par defaut - meme discipline) raffinent LOCALEMENT les
    jonctions multi-region a forte divergence directionnelle (CLAUDE.md
    limite #47/#48), generalisant le mecanisme ci-dessus au-dela du seul
    cas plateau. Contrairement au plateau, une jonction n'a pas de plan de
    coupe unique : pas de decoupe de surface, juste un seed interne
    supplementaire (voir `_merged_junction_seed`) - PLUSIEURS zones
    detectees sont FUSIONNEES en un seul seed avant d'etre passees a
    `add_region`, jamais une par zone : mesure directe
    (scripts/measure_add_region_n_zones.py) que le bug limite #42 se
    generalise a N regions actives simultanement (aucune des zones fines
    n'atteignait le facteur de densite attendu, une meme MOINS dense que
    la region censee etre grossiere) - le fallback fusionne fonctionne
    correctement (scripts/measure_add_region_merged_zone.py, facteur ~48x).
    """
    if not mesh.is_watertight:
        raise MeshIngestError(
            "Impossible de tetraedriser un mesh non etanche apres reparation. "
            "Verifiez la geometrie source."
        )

    # Strategie de repli, du moins invasif au plus invasif. TetGen est le
    # seul juge fiable des auto-intersections qui le bloquent : ce sont
    # souvent des contacts EXACTS ("Two segments exactly intersect", "A
    # vertex lies on a facet") que ni trimesh (watertight/winding
    # coherents), ni manifold3d (rend le maillage a l'identique), ni un test
    # geometrique de traversee ne signalent - mesure sur
    # res/test_models/death.stl : 3 triangles fautifs sur 12662, zero
    # detecte par le test SAT de _count_self_intersecting_faces. D'ou une
    # reparation REACTIVE (on repare parce que TetGen a refuse), jamais
    # preventive : un maillage deja accepte n'est jamais modifie.
    attempts: list[tuple[str, trimesh.Trimesh, Optional[float]]] = [
        ("tel quel", mesh, near_platform_resolution),
    ]

    # Resolution GEOMETRIQUE des auto-intersections (decoupe le long de la
    # courbe d'intersection + extraction de l'enveloppe exterieure par
    # nombre d'enroulement) - tentee EN PREMIER, avant pymeshfix :
    # contrairement a pymeshfix (qui reecrit tout le maillage), elle
    # preserve la topologie et l'aire de la geometrie non affectee, donc
    # strictement moins invasive. Best-effort : ne leve jamais, retourne
    # None si le defaut n'est pas de cette nature ou si le resultat ne
    # satisfait pas ses propres post-conditions (voir docstring du module).
    # Meme garde-fou anti-effondrement que pour pymeshfix ci-dessous (limite
    # #51/#52) : un resultat radicalement plus petit que l'entree est rejete
    # plutot qu'accepte silencieusement. Mesure sur turkish.stl (session
    # ayant introduit ce mecanisme) : reduit les auto-intersections de 50 a
    # 0 faces SAT et de 51 a 2 segments nommes par TetGen, mais ne suffit
    # pas encore a lui seul a faire accepter ce modele precis (recouvrements
    # coplanaires non traites) - reste neanmoins utile en cascade, avant le
    # repli plus destructeur pymeshfix.
    resolved = resolve_self_intersections(mesh)
    resolved_collapsed = False
    if resolved is not None and len(resolved.faces) < 0.5 * len(mesh.faces):
        resolved = None
        resolved_collapsed = True
    if resolved is not None:
        # `resolve_self_intersections` peut rendre un mesh etanche et
        # connexe mais avec quelques auto-intersections RESIDUELLES,
        # introduites par son propre soudage de trous minuscules
        # (`_close_small_boundary_loops` - fusionner des sommets proches
        # peut faire chevaucher deux triangles fins voisins). Ce n'est plus
        # le defaut d'origine (recouvrement CSG transversal), c'est un
        # artefact local du soudage - nettoye ici avec le meme mecanisme
        # que ci-dessous (retrait + retriangulation Klincsek), sans la
        # garde de caracteristique d'Euler stricte (inadaptee : le soudage
        # de trous a deja legitimement change l'euler avant cet appel).
        # Mesure sur turkish.stl : 8 auto-intersections residuelles -> 0,
        # corps unique preserve, ecart de volume 2.9e-6%.
        residual_sat = _find_self_intersecting_faces(resolved)
        if residual_sat:
            cleaned_resolved = _repair_local_self_intersections(
                resolved, require_euler_preserved=False,
            )
            if cleaned_resolved is not None:
                resolved = cleaned_resolved
        # Le test SAT ci-dessus ignore volontairement toute paire de faces
        # partageant un sommet (contact par construction) - or un contact
        # EXACT entre deux faces qui partagent DEJA un sommet ailleurs
        # (ex: un eventail de triangles fins autour d'un sommet issu du
        # soudage) lui est invisible, alors que TetGen le refuse.
        # `_finish_geometric_resolution` boucle sur TetGen lui-meme comme
        # oracle final pour ce residu : mesure sur turkish.stl, un dernier
        # segment nomme (2 faces partageant le sommet 890) restait apres
        # le nettoyage SAT ci-dessus ; cible et repare en 1 iteration
        # supplementaire, TetGen aboutit alors sur la geometrie COMPLETE.
        resolved = _finish_geometric_resolution(resolved, resolution)
        attempts.append(("apres resolution geometrique", resolved, near_platform_resolution))

    # Meme garde-fou anti-effondrement que `_clean_and_repair` (limite #51) :
    # MeshFix peut rendre un mesh watertight en supprimant une portion
    # entiere de la geometrie plutot qu'en la reparant. Bug reel mesure sur
    # `turkish.stl` : ce deuxieme appel a `_repair_with_pymeshfix`, propre a
    # la cascade de tetraedrisation et jusqu'ici SANS garde-fou (seul
    # `_clean_and_repair` en avait un), acceptait un resultat passant de
    # 10104 a 1488 faces - watertight, donc silencieusement accepte par la
    # cascade, mais amputant tout le haut du modele (z tronque de 90.78 a
    # 63.26mm). TetGen tetraedrise alors ce mesh mutile SANS AUCUNE erreur
    # ni warning (il est parfaitement valide, juste faux) : le pipeline
    # entier (step2 a step8) traite ensuite fidelement un modele auquel il
    # manque les deux tiers de sa hauteur, sans jamais lever d'exception -
    # silence total jusqu'a l'inspection visuelle du GLB final.
    repaired = _repair_with_pymeshfix(mesh)
    pymeshfix_collapsed = False
    if repaired is not None and len(repaired.faces) < 0.5 * len(mesh.faces):
        # Effondrement detecte : ce resultat n'est utilisable NULLE PART
        # dans la cascade, ni comme tentative dediee ni comme source de
        # repli pour "sans raffinement plateau" ci-dessous (meme mesh
        # mutile, memes deux tiers de hauteur disparus).
        repaired = None
        pymeshfix_collapsed = True
    if repaired is not None:
        attempts.append(("apres reparation pymeshfix", repaired, near_platform_resolution))

    # Dernier recours : abandonner le raffinement local pres du plateau.
    # `_split_for_platform_refinement` recolle deux coques obtenues par
    # slice_plane(cap=True) ; leurs faces de bouchage sont coplanaires et
    # superposees le long de la frontiere de coupe, ce que TetGen peut
    # refuser sur une geometrie dense qui touche ce plan en de nombreux
    # points - mesure sur death.stl repare : 107037 tetraedres SANS
    # raffinement, refus net AVEC. Mieux vaut un maillage un peu plus
    # grossier pres du plateau qu'aucun maillage du tout.
    #
    # Ordre de preference : `resolved` (resolution geometrique, la moins
    # invasive) > `repaired` (pymeshfix) > `mesh` (brut). Bug reel mesure
    # en integrant ce repli avec `resolved` : la source ne considerait
    # QUE `repaired`/`mesh`, jamais `resolved` - sur turkish.stl, avec
    # `near_platform_resolution` actif, l'attempt "apres resolution
    # geometrique" echouait pour une raison PROPRE au raffinement plateau
    # ("split_segment", nouvelle variante limite #44 sur le PLC combine
    # avec ce mesh precis), mais ce repli retombait alors sur `mesh` NON
    # REPARE (51 auto-intersections brutes) au lieu du mesh deja nettoye
    # par la resolution geometrique - echec garanti pour une raison sans
    # aucun rapport avec celle qui vient d'etre reglee.
    if near_platform_resolution is not None:
        source = resolved if resolved is not None else (
            repaired if repaired is not None else mesh
        )
        attempts.append(("sans raffinement plateau", source, None))

    last_error: Optional[Exception] = None
    for _label, candidate, npr in attempts:
        try:
            return _tetrahedralize_once_isolated(
                candidate, resolution, npr, near_platform_margin, platform_z,
                junction_zones, junction_refinement_resolution,
            )
        except _TetgenCrashError as exc:
            # Crash NATIF (segfault/access violation) dans le sous-processus
            # isole, jamais une RuntimeError Python normale - traite comme
            # le meme signal d'echec que "self-intersection" ci-dessous : le
            # repli suivant de la cascade (reparation pymeshfix, puis sans
            # raffinement plateau) a la meme chance de reussir, un crash sur
            # auto-intersection etant la cause la plus probable mesuree
            # (voir CLAUDE.md limite #43 et le worker isole plus bas).
            # DOIT rester AVANT `except RuntimeError` ci-dessous :
            # `_TetgenCrashError` en est une sous-classe, et Python retient
            # le PREMIER except applicable - inverser l'ordre ferait
            # tomber tout crash natif dans la branche "self-intersection"
            # (dont le message ne contient jamais ce texte), le relevant
            # immediatement au lieu de tenter le repli suivant. Bug reel
            # mesure sur `death.stl` modifie par l'utilisateur (2340
            # triangles) : le crash etait bien detecte (code de sortie
            # 3221225477 = 0xC0000005, access violation Windows) mais
            # remontait quand meme jusqu'a l'appelant au lieu de declencher
            # la reparation pymeshfix, a cause de cet ordre incorrect.
            last_error = exc
        except RuntimeError as exc:
            if not _is_recoverable_tetgen_error(exc):
                raise
            last_error = exc

    if repaired is not None:
        hint = ""
    elif pymeshfix_collapsed:
        hint = (
            " La reparation `pymeshfix` a produit un resultat watertight "
            "mais radicalement plus petit que le mesh source (effondrement "
            "geometrique, rejete par garde-fou) - le defaut source n'a pas "
            "de reparation automatique fiable sur ce mesh."
        )
    else:
        hint = (
            " La reparation automatique n'a pas pu s'appliquer (dependance "
            "optionnelle `pymeshfix` absente ou en echec) - l'installer peut suffire."
        )
    if resolved_collapsed:
        hint += (
            " La resolution geometrique des auto-intersections a produit un "
            "resultat radicalement plus petit que le mesh source (rejete par "
            "garde-fou)."
        )
    raise MeshIngestError(
        f"Le maillage contient des auto-intersections (ou un autre defaut) que "
        f"TetGen refuse ou sur lequel il plante, et que la reparation "
        f"automatique n'a pas suffi a corriger : {last_error}.{hint}"
    )


class _TetgenCrashError(RuntimeError):
    """
    Leve par `_tetrahedralize_once_isolated` quand le sous-processus TetGen
    se termine par un crash natif (code de sortie non nul, typiquement
    SIGSEGV) plutot qu'une exception Python normale - distinct d'une
    `RuntimeError` "self-intersection" habituelle UNIQUEMENT pour le
    message d'erreur ; traitee de facon identique par la cascade de repli
    de `_tetrahedralize` (meme `except` logique, capture separee pour
    documenter explicitement la distinction).
    """


class _TetgenNamedConflictError(RuntimeError):
    """
    `RuntimeError` normale (echec Python, pas un crash) MAIS accompagnee du
    texte complet imprime par TetGen (`.tetgen_output`), qui peut nommer
    precisement les sommets/segments/facettes en conflit ("Two segments
    exactly intersect", "A vertex lies on a facet", etc. - meme motif que
    la limite #43). Ces avertissements ne survivent normalement pas au
    passage par le sous-processus isole (`str(exc)` seul remonte, message
    generique) - `_tetrahedralize_worker.py` les capture desormais via
    `redirect_stdout`/`redirect_stderr` autour de l'appel reel. Utilisee
    par `_tetrahedralize` pour cibler une reparation locale
    (`_repair_local_self_intersections(..., extra_seed_faces=...)`) sur
    exactement les faces nommees, au lieu d'un repli generique.
    """

    def __init__(self, message: str, tetgen_output: str = ""):
        super().__init__(message)
        self.tetgen_output = tetgen_output


def _tetrahedralize_once_isolated(
    mesh: trimesh.Trimesh,
    resolution: float,
    near_platform_resolution: Optional[float],
    near_platform_margin: float,
    platform_z: float,
    junction_zones: Optional[list[JunctionZone]],
    junction_refinement_resolution: Optional[float],
) -> TetMesh:
    """
    Fait tourner UNE tentative de `_tetrahedralize_once` dans un
    sous-processus isole (`pipeline/_tetrahedralize_worker.py`), pour
    survivre a un crash NATIF de TetGen (segfault) que ni `except
    RuntimeError` ni aucun autre mecanisme Python ne peut intercepter
    in-process - un tel crash tue tout le process qui l'a appele.

    Mesure concretement sur une version modifiee de `res/test_models/death.stl`
    (2340 triangles) : le process se terminait avec le code de sortie 139
    (SIGSEGV) des la 1ere tentative de la cascade de reparation existante,
    juste apres les memes avertissements TetGen deja documentes a la limite
    #43 ("Two segments exactly intersect" etc.) - AVANT que la reparation
    pymeshfix n'ait la moindre chance de s'executer, puisque le process qui
    aurait du l'appeler etait deja mort. L'isolation par sous-processus
    s'applique donc a CHAQUE tentative de la cascade (y compris la 1ere,
    "tel quel"), jamais seulement en dernier recours - c'est la seule facon
    de garantir qu'un crash sur le mesh original n'empeche pas d'essayer la
    reparation ensuite.

    Cout accepte : overhead de lancement d'un sous-processus Python par
    tentative (import de tetgen/scipy/trimesh a chaque fois) - largement
    inferieur au temps de la tetraedrisation elle-meme sur tout modele non
    trivial, et strictement necessaire pour la garantie de robustesse
    demandee (voir CLAUDE.md limite #46/#49, section "point ouvert
    death.stl").

    Serialisation deliberement minimale (vertices/faces bruts, pas un objet
    Trimesh complet) : plus robuste au pickling inter-version que de
    transporter des objets trimesh a travers la frontiere de process.

    Raises:
        RuntimeError: relayee telle quelle si le worker a lui-meme leve une
            RuntimeError Python normale (ex: message "self-intersection" de
            TetGen) - meme comportement qu'un appel in-process, la cascade
            de `_tetrahedralize` la traite exactement comme avant.
        _TetgenCrashError: si le sous-processus s'est termine par un crash
            natif (code de sortie non nul, ou fichier de sortie absent/
            illisible) plutot qu'une exception Python.
        MeshIngestError: relayee telle quelle si le worker a leve cette
            exception (entree structurellement invalide, pas un defaut de
            maillage a reparer).
    """
    with tempfile.TemporaryDirectory(prefix="apex_tetgen_") as tmp_dir:
        input_path = os.path.join(tmp_dir, "input.pkl")
        output_path = os.path.join(tmp_dir, "output.pkl")
        stdout_path = os.path.join(tmp_dir, "stdout.txt")
        stderr_path = os.path.join(tmp_dir, "stderr.txt")

        with open(input_path, "wb") as f:
            pickle.dump(
                {
                    "vertices": np.asarray(mesh.vertices, dtype=np.float64),
                    "faces": np.asarray(mesh.faces, dtype=np.int64),
                    "resolution": resolution,
                    "near_platform_resolution": near_platform_resolution,
                    "near_platform_margin": near_platform_margin,
                    "platform_z": platform_z,
                    "junction_zones": junction_zones,
                    "junction_refinement_resolution": junction_refinement_resolution,
                },
                f,
            )

        # Les avertissements nommes de TetGen ("Two segments exactly
        # intersect [869,836]...") sortent de l'extension C du binding par
        # ses descripteurs OS bruts, en contournant l'objet Python
        # `sys.stdout` - une redirection Python (`contextlib.redirect_*`,
        # essayee et rejetee) ne les voit jamais. `subprocess.run(...,
        # capture_output=True)` (PIPE) ne suffit PAS non plus : mesure
        # directe sur turkish.stl, le texte s'arretait net juste apres
        # "Recovering boundaries..." (avant les avertissements eux-memes),
        # signature d'un tampon de la libc du binding jamais vide vers un
        # PIPE avant la fin du sous-processus. Rediriger vers de VRAIS
        # FICHIERS resout le probleme (verifie : texte complet recupere,
        # y compris les indices de sommets) - la nature du descripteur
        # (fichier vs pipe) change le comportement de bufferisation cote
        # OS, meme si Python ne voit aucune difference de son point de vue.
        with open(stdout_path, "wb") as stdout_file, open(stderr_path, "wb") as stderr_file:
            completed = subprocess.run(
                [sys.executable, "-m", "pipeline._tetrahedralize_worker", input_path, output_path],
                cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                stdout=stdout_file,
                stderr=stderr_file,
            )

        if completed.returncode != 0 or not os.path.isfile(output_path):
            with open(stderr_path, "rb") as f:
                stderr_tail = f.read().decode("utf-8", errors="replace")[-500:]
            raise _TetgenCrashError(
                f"Le sous-processus TetGen s'est termine anormalement "
                f"(code de sortie {completed.returncode}, probablement un "
                f"crash natif) - stderr : {stderr_tail}"
            )

        with open(output_path, "rb") as f:
            result = pickle.load(f)

        if not result["ok"]:
            error_type = result.get("error_type", "RuntimeError")
            message = result["error"]
            if error_type == "MeshIngestError":
                raise MeshIngestError(message)
            with open(stdout_path, "rb") as f:
                tetgen_output = f.read().decode("utf-8", errors="replace")
            with open(stderr_path, "rb") as f:
                tetgen_output += f.read().decode("utf-8", errors="replace")
            raise _TetgenNamedConflictError(message, tetgen_output)

    tet_vertices = result["tet_vertices"]
    tetrahedra = result["tetrahedra"]
    tet_adjacency = _build_tet_adjacency(tetrahedra, len(tet_vertices))
    return TetMesh(
        tet_vertices=np.asarray(tet_vertices, dtype=np.float64),
        tetrahedra=np.asarray(tetrahedra, dtype=np.int64),
        tet_adjacency=tet_adjacency,
    )


def _tetrahedralize_once(
    mesh: trimesh.Trimesh,
    resolution: float,
    near_platform_resolution: Optional[float],
    near_platform_margin: float,
    platform_z: float,
    junction_zones: Optional[list[JunctionZone]] = None,
    junction_refinement_resolution: Optional[float] = None,
) -> TetMesh:
    """
    Une passe de tetraedrisation, sans aucune tentative de reparation :
    laisse remonter telle quelle l'erreur de TetGen. Separee de
    `_tetrahedralize` pour que celle-ci puisse la rejouer sur un maillage
    repare sans dupliquer la construction des regions/trous.

    Jusqu'a 3 regions TetGen actives simultanement : "slab" (pres du
    plateau), "junction" (jonctions multi-region fusionnees en UN seul
    seed, cf. `_merged_junction_seed`), "rest" (le reste, a `resolution`).
    Contrairement au plateau, une jonction n'a pas de plan de coupe : le
    PLC combine (`combined_vertices`/`combined_faces`) reste celui du
    plateau seul (ou le mesh brut si le plateau n'est pas actif) - le seed
    de jonction est juste un point interne supplementaire, jamais une
    decoupe de surface.
    """
    import tetgen

    max_tet_volume = resolution ** 3
    hole_points = _find_internal_cavities(mesh)

    merged_junction_seed = None
    if junction_zones and junction_refinement_resolution is not None:
        merged_junction_seed = _merged_junction_seed(mesh, junction_zones)

    if near_platform_resolution is not None:
        combined_vertices, combined_faces, slab_seed, rest_seed = _split_for_platform_refinement(
            mesh, near_platform_margin, platform_z
        )
        if merged_junction_seed is not None:
            # Garde-fou EXPLICITE (pas un silence) : si slab_seed ou
            # rest_seed tombent a l'interieur du rayon d'une zone de
            # jonction, TetGen appliquerait le volume FIN de la jonction a
            # toute la region plateau/rest - meme bug mesure et corrige
            # pour le cas jonction seule (_find_rest_seed_outside_junction_zones).
            # Non re-derive automatiquement ici (deplacer slab_seed/rest_seed
            # casserait la garantie deja etablie par _split_for_platform_refinement
            # d'etre chacun a l'interieur de SA PROPRE coque) - signale plutot
            # qu'un resultat silencieusement corrompu.
            for seed_name, seed in [("slab_seed", slab_seed), ("rest_seed", rest_seed)]:
                for zone in junction_zones:
                    if np.linalg.norm(seed - zone.centroid) <= zone.radius:
                        raise MeshIngestError(
                            f"{seed_name} ({seed}) tombe a l'interieur d'une zone de "
                            f"jonction (centre={zone.centroid}, rayon={zone.radius}) - "
                            "raffinement plateau et raffinement jonction se chevauchent "
                            "sur ce modele, non supporte (near_platform_margin ou les "
                            "zones detectees devraient etre ajustes)."
                        )
        tet = tetgen.TetGen(combined_vertices, combined_faces)
        for hole_point in hole_points:
            tet.add_hole(hole_point)
        if merged_junction_seed is not None:
            # ORDRE CRITIQUE quand une 3e region (jonction) est active,
            # mesure pas suppose (scripts/measure_add_region_*.py, session
            # limite #48) : la region la plus GROSSIERE ("rest") doit etre
            # add_region(1, ...), sinon la densite explose sur TOUT le
            # maillage (mesure sur le tripode synthetique : far<3mm passe de
            # 34 a 1150 tetraedres selon l'ordre, memes seeds/volumes). Ce
            # meme defaut existe deja avec seulement 2 regions
            # (slab=1,rest=2 vs rest=1,slab=2 : far<3mm 847 vs 29 sur un
            # modele synthetique) mais N'EST PAS corrige ici pour le cas
            # SANS jonction (slab_seed/rest_seed) : `l_shape_overhang.stl`
            # (modele de reference du projet, deja valide) produit
            # exactement 1677 tetraedres avec l'ordre actuel slab=1/rest=2 -
            # identique au chiffre deja documente CLAUDE.md limite #36,
            # aucune regression mesuree sur ce cas. Inverser l'ordre pour
            # TOUS les appels aurait ete un changement de comportement non
            # mesure sur le cas deja valide - reserve au nouveau cas (3
            # regions) ou le defaut est confirme.
            tet.add_region(1, rest_seed, max_tet_volume)
            tet.add_region(2, slab_seed, near_platform_resolution ** 3)
            tet.add_region(3, merged_junction_seed, junction_refinement_resolution ** 3)
        else:
            tet.add_region(1, slab_seed, near_platform_resolution ** 3)
            tet.add_region(2, rest_seed, max_tet_volume)
        # -A active les regions (chaque region utilise SON PROPRE max_vol
        # declare via add_region) ; -a sans valeur numerique explicite ne
        # fixe PAS de plafond global unique - verifie fonctionnel, aucun
        # crash, resultat mesure (CLAUDE.md).
        tet.tetrahedralize(switches="pq1.2Aa")
    elif merged_junction_seed is not None:
        # Pas de raffinement plateau actif : 2 regions seulement ("junction"
        # + "rest"), sur le mesh brut (pas de PLC combine necessaire, aucune
        # decoupe de surface pour une jonction - voir _merged_junction_seed).
        # Seed "rest" garanti HORS de toute zone de jonction (pas seulement
        # a l'interieur du mesh) - voir docstring de
        # _find_rest_seed_outside_junction_zones pour le bug mesure sinon.
        rest_seed = _find_rest_seed_outside_junction_zones(mesh, junction_zones)
        tet = tetgen.TetGen(mesh.vertices, mesh.faces)
        for hole_point in hole_points:
            tet.add_hole(hole_point)
        # ORDRE CRITIQUE, mesure pas suppose (scripts/measure_add_region_*.py) :
        # "rest" (grossier) DOIT etre la region 1, "junction" (fin) la region 2
        # - inverser l'ordre (junction=1, rest=2) fait EXPLOSER la densite sur
        # tout le maillage (mesure : 3080 tetraedres corrects avec rest=1, vs
        # 46000+ avec junction=1, memes seeds, meme mesh, meme volumes cibles).
        # Comportement TetGen non documente, mais reproductible a l'identique
        # sur ce binding (0.8.4) - meme discipline que la limite #42 (TetGen
        # est le seul juge fiable de son propre comportement d'add_region).
        tet.add_region(1, rest_seed, max_tet_volume)
        tet.add_region(2, merged_junction_seed, junction_refinement_resolution ** 3)
        tet.tetrahedralize(switches="pq1.2Aa")
    else:
        tet = tetgen.TetGen(mesh.vertices, mesh.faces)
        for hole_point in hole_points:
            tet.add_hole(hole_point)
        tet.tetrahedralize(switches=f"pq1.2a{max_tet_volume:.10f}")

    # Lire depuis les attributs de l'objet plutot que la valeur de retour :
    # le nombre de valeurs retournees par tetrahedralize() varie selon la
    # version de tetgen installee, tandis que .node/.elem sont stables.
    tet_vertices = tet.node
    tetrahedra = tet.elem

    tet_adjacency = _build_tet_adjacency(tetrahedra, len(tet_vertices))

    return TetMesh(
        tet_vertices=np.asarray(tet_vertices, dtype=np.float64),
        tetrahedra=np.asarray(tetrahedra, dtype=np.int64),
        tet_adjacency=tet_adjacency,
    )


def _build_tet_adjacency(tetrahedra: np.ndarray, node_count: int) -> sp.csr_matrix:
    """
    Construit la matrice d'adjacence creuse noeud-noeud du maillage
    tetraedrique (1 si deux noeuds partagent une arete d'un tetraedre).
    Reutilisee par l'assemblage du Laplacien discret a l'etape 4 du pipeline.
    """
    if _cpp.enabled():
        # L'ordre des aretes n'a aucune importance ici : `csr_matrix` trie ses
        # indices a la construction, donc la matrice finale est identique quel
        # que soit l'ordre fourni.
        edge_array = _cpp.module().tet_unique_edges(np.asarray(tetrahedra, dtype=np.int64))
        if len(edge_array) == 0:
            return sp.csr_matrix((node_count, node_count))
        rows = edge_array[:, 0]
        cols = edge_array[:, 1]
    else:
        edges = set()
        for tet in tetrahedra:
            for i in range(4):
                for j in range(i + 1, 4):
                    a, b = int(tet[i]), int(tet[j])
                    edges.add((min(a, b), max(a, b)))

        if len(edges) == 0:
            return sp.csr_matrix((node_count, node_count))

        rows, cols = zip(*edges)
        rows = np.array(rows)
        cols = np.array(cols)

    data = np.ones(len(rows) * 2)
    all_rows = np.concatenate([rows, cols])
    all_cols = np.concatenate([cols, rows])

    adjacency = sp.csr_matrix(
        (data, (all_rows, all_cols)), shape=(node_count, node_count)
    )
    return adjacency


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    mesh_path: str,
    resolution: float,
    allow_multi_mesh: bool = False,
    near_platform_resolution: Optional[float] = None,
    near_platform_margin: float = 5.0,
    platform_z: Optional[float] = None,
    junction_zones: Optional[list[JunctionZone]] = None,
    junction_refinement_resolution: Optional[float] = None,
) -> tuple[MeshData, TetMesh, IngestReport]:
    """
    Charge, nettoie et prepare le mesh source pour le reste du pipeline.

    Args:
        mesh_path: chemin vers le fichier STL/3MF source.
        resolution: taille cible de tetraedre (FieldSolverParams.resolution),
            utilisee pour controler la densite de la tetraedrisation.
        allow_multi_mesh: si True, fusionne les geometries multiples d'une
            Scene plutot que de lever une erreur.
        near_platform_resolution: taille de tetraedre cible pres du plateau,
            plus fine que `resolution` (opt-in, None = desactive par defaut,
            comportement historique de resolution globale unique inchange).
            Corrige une erreur de discretisation FEM mesuree pres du plateau
            (CLAUDE.md) sans payer le cout d'un maillage fin partout.
        near_platform_margin: hauteur (mm) au-dessus de `platform_z` sur
            laquelle le raffinement s'applique, avec une transition lisse
            vers `resolution` au-dela. Ignore si `near_platform_resolution`
            est None.
        platform_z: niveau du plateau. `None` (defaut) = DETECTION
            AUTOMATIQUE depuis le mesh nettoye (`min(z)` des sommets) -
            meme discipline de derivation geometrique que
            `_find_guaranteed_interior_point`/`_find_internal_cavities`
            plutot que d'exiger un parametre externe suppose correct. Une
            valeur explicite reste un override total. La valeur REELLEMENT
            utilisee (fournie ou deduite) est toujours remontee via
            `IngestReport.platform_z_used` - jamais un silence sur une
            valeur derivee, meme discipline que
            `self_intersecting_faces_before`. Corrige un bug reel mesure
            sur `res/test_models/death.stl` (min(z)=1.245, pas 0.0) : avec
            l'ancien defaut `platform_z=0.0` fige, `_compute_platform_contact`
            (step3) ne detectait AUCUNE face de contact plateau sur ce
            modele, cassant en cascade l'ancrage plateau du solveur FEM
            (step4) - voir CLAUDE.md.
        junction_zones: zones de jonction multi-region a raffiner localement
            (opt-in, None = desactive par defaut), typiquement produites par
            `pipeline.step1b_junction_detection.detect_junction_zones` sur
            une premiere passe legere (mesh + step2 + step3) orchestree par
            `main.run_pipeline` - `step1_ingest.run()` lui-meme ne les
            calcule jamais, il les recoit deja prealablement calculees.
            Corrige la non-monotonie du champ scalaire aux jonctions
            (CLAUDE.md limite #47/#48), portee laissee explicitement ouverte
            par le mecanisme `near_platform_resolution` ci-dessus (limite
            #36 : "traiter le cas jonction... non fait").
        junction_refinement_resolution: taille de tetraedre cible aux zones
            de jonction, plus fine que `resolution` (opt-in, None =
            desactive par defaut). Ignore si `junction_zones` est vide ou
            None.

    Returns:
        Tuple (MeshData, TetMesh, IngestReport).

    Raises:
        MeshIngestError: si le fichier est illisible, vide, ou ne peut pas
            etre reduit a un volume tetraedrisable meme apres reparation.
    """
    raw_mesh = _load_raw_mesh(mesh_path, allow_multi_mesh)
    cleaned_mesh, report = _clean_and_repair(raw_mesh)

    resolved_platform_z = (
        float(cleaned_mesh.vertices[:, 2].min()) if platform_z is None else platform_z
    )
    report.platform_z_used = resolved_platform_z

    mesh_data = _build_mesh_data(cleaned_mesh)
    tet_mesh = _tetrahedralize(
        cleaned_mesh, resolution, near_platform_resolution, near_platform_margin, resolved_platform_z,
        junction_zones, junction_refinement_resolution,
    )

    report.tetrahedra_generated = tet_mesh.tet_count

    return mesh_data, tet_mesh, report