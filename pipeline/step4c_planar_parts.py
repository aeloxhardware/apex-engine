"""
Etape 4c - PARTIES A COUCHES PLANES (2026-10-19, opt-in : `FieldSolverParams.planar_parts`).

Alternative sans champ scalaire resolu : le volume est decoupe en parties (region de step3 x
composante connexe, comme le multibranche), et chaque partie est tranchee en PLANS
perpendiculaires a SA direction de croissance, orientee en S'ELOIGNANT de son attache (le
plateau pour une racine, sinon les parties de profondeur moindre) puis projetee sur le demi-plan
radial atteignable.

DECOUPAGE DES PARTIES PENDANTES. Une partie tranchee en plans ne doit contenir aucun
tetraedre qui commence en l'air : un tetraedre est ATTEIGNABLE s'il est ancre (au plateau ou
contre une partie de profondeur moindre) ou s'il touche un tetraedre atteignable de la meme
partie dont le centre est plus bas selon la direction de la partie. Les tetraedres non atteignables
(le crochet sous un bras horizontal) forment des parties filles, orientees en s'eloignant de leur
porteur, jusqu'a ce que tout soit atteignable (`MAX_SPLIT_ROUNDS`).

Realisation : step5 est reutilise tel quel avec un maillage dont les noeuds sont DUPLIQUES par
partie et un champ lineaire par partie, phi_b(x) = d_b . (x - x_min_b) + debut_b, ou debut_b
cumule les plages des parties precedentes dans l'ordre d'impression (profondeur, indice) : deux
parties ne partagent jamais une isovaleur, et l'ordre par isovaleur EST l'ordre des parties.

LIMITES : coutures planes entre parties (pas de raccord lisse) ; les `Constraint[]` ne sont pas
honorees ; le maillage rendu REMPLACE celui de step1 pour les etapes suivantes ; step5b (ordre
contraint par la tete) brouille l'ordre des parties : a desactiver avec cette voie.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.csgraph as spgraph

from core.mesh_types import MeshData, RegionDependencyGraph, ScalarField, TetMesh
from pipeline import step4_field_solver as s4
from pipeline.step4b_multibranch import _branch_labels

MAX_SPLIT_ROUNDS = 4
# Une naissance en l'air (composante sans ancre qui apparait en montant selon la direction de la partie)
# n'est consideree comme PENDANTE que si elle ne rejoint le reste qu'apres au moins cette hauteur : en
# dessous, c'est la rugosite du maillage d'une paroi plane.
PLANAR_MIN_PERSISTENCE_MM = 1.0
# Morceau de tranche de Reeb fondu dans son voisin en dessous de cette taille.
REEB_MIN_COMPONENT_TETS = 15
# Une branche est coupee a un coude quand sa direction de croissance tourne de plus de cet angle.
REEB_BEND_DEG = 35.0
# La tranche fait au moins ce multiple de la distance mediane entre centres de tetraedres voisins.
REEB_BAND_TET_FACTOR = 3.0


@dataclass
class PlanarPartsDiagnostic:
    initial_parts: int = 0
    final_parts: int = 0
    split_rounds: int = 0
    hanging_tets: list[int] = field(default_factory=list)


LAST_DIAGNOSTIC = PlanarPartsDiagnostic()


def _analyze(label, tets, verts, adjacency, is_platform):
    """Graphe des parties, racines, profondeurs depuis le plateau."""
    B = int(label.max()) + 1
    coo = sp.triu(adjacency, k=1).tocoo()
    la, lb = label[coo.row], label[coo.col]
    diff = la != lb
    bg = sp.csr_matrix((np.ones(int(diff.sum())), (la[diff], lb[diff])), shape=(B, B))
    bg = ((bg + bg.T) > 0).astype(float)
    branch_tets = [np.where(label == b)[0] for b in range(B)]
    branch_nodes = [np.unique(tets[ids]) for ids in branch_tets]
    roots = [b for b in range(B) if int(is_platform[branch_nodes[b]].sum()) >= 3]
    depth = (spgraph.dijkstra(bg, directed=False, indices=roots, min_only=True, unweighted=True)
             if roots else np.full(B, np.inf))
    while np.isinf(depth).any():
        cand = np.where(np.isinf(depth))[0]
        b0 = int(min(cand, key=lambda b: verts[branch_nodes[b], 2].min()))
        extra = spgraph.dijkstra(bg, directed=False, indices=[b0], min_only=True, unweighted=True)
        depth = np.where(np.isinf(depth), extra, depth)
        depth[b0] = 0.0
    return B, bg, branch_tets, branch_nodes, roots, depth


def _part_directions(B, bg, branch_nodes, roots, depth, verts, is_platform, base_dir, up, radial_origin):
    """Une direction par partie : celle de sa region, orientee en s'eloignant de l'attache."""
    part_dir = np.zeros((B, 3))
    for b in range(B):
        d = np.asarray(base_dir[b], dtype=float)
        n = np.linalg.norm(d)
        d = d / n if n > 1e-12 else up.copy()
        nb = branch_nodes[b]
        if b in roots:
            att = nb[is_platform[nb]]
        else:
            parents = [branch_nodes[q] for q in range(B) if depth[q] < depth[b] and bg[b, q] > 0]
            att = np.intersect1d(nb, np.concatenate(parents)) if parents else np.empty(0, dtype=int)
        anchor = verts[att].mean(axis=0) if len(att) else verts[nb][np.argmin(verts[nb] @ up)]
        away = (verts[nb].mean(axis=0) - anchor) @ d
        if away < 0.0 or (abs(away) < 1e-9 and d @ up < 0.0):
            d = -d
        part_dir[b] = d
    if radial_origin is not None:
        ctr = np.array([verts[branch_nodes[b]].mean(axis=0) for b in range(B)])
        proxy = TetMesh(tet_vertices=ctr, tetrahedra=np.tile(np.arange(B)[:, None], (1, 4)))
        part_dir, _ = s4._project_to_reachable_radial_tilt(proxy, part_dir, radial_origin, up)
    return part_dir


def _anchored_tets(label, depth, parent, tets, adjacency, is_platform):
    """Ancres : plateau (racine), contact avec une partie de profondeur moindre (partie d'origine) ou avec
    la partie PARENT (partie fille)."""
    T = len(label)
    coo = sp.triu(adjacency, k=1).tocoo()
    a, b = coo.row, coo.col
    la, lb = label[a], label[b]
    anchored = np.zeros(T, dtype=bool)
    for x, lx, ly in ((a, la, lb), (b, lb, la)):
        child_anchor = (parent[lx] >= 0) & (ly == parent[lx])
        orig_anchor = (parent[lx] < 0) & (parent[ly] < 0) & (depth[ly] < depth[lx])
        anchored[x[child_anchor | orig_anchor]] = True
    anchored |= is_platform[tets].any(axis=1) & (depth[label] == 0) & (parent[label] < 0)
    return anchored


def _hanging_tets(label, branch_tets, part_dir, centers, adjacency, anchored, min_tets, persistence):
    """Tetraedres PENDANTS : balayage des tetraedres de chaque partie du plus bas au plus haut selon
    sa direction, en regroupant par connexite (union-find). Un groupe ne vient d'une ancre que s'il en
    contient une ; un groupe ne en l'air qui ne rejoint un groupe plus ancien qu'apres `persistence` mm
    est pendant (le crochet sous un bras) : ses tetraedres sont rendus avec un identifiant de groupe.
    Rend hang_id (T,) : -1 = non pendant."""
    T = len(label)
    indptr, indices = adjacency.indptr, adjacency.indices
    hang = np.full(T, -1, dtype=np.int64)
    next_id = 0
    for b, ids in enumerate(branch_tets):
        if len(ids) == 0:
            continue
        phi = centers[ids] @ part_dir[b]
        order = ids[np.argsort(phi, kind="stable")]
        phi_all = {}
        parent_uf = {}
        birth = {}
        anch = {}
        members = {}

        def find(x):
            while parent_uf[x] != x:
                parent_uf[x] = parent_uf[parent_uf[x]]
                x = parent_uf[x]
            return x
        any_anchor = bool(anchored[ids].any())
        for t in order:
            t = int(t)
            pt = float(centers[t] @ part_dir[b])
            roots = set()
            for n in indices[indptr[t]:indptr[t + 1]]:
                n = int(n)
                if n in parent_uf and label[n] == b:
                    roots.add(find(n))
            parent_uf[t] = t
            birth[t] = pt
            anch[t] = bool(anchored[t]) or (not any_anchor and len(parent_uf) == 1)
            members[t] = [t]
            if not roots:
                continue
            cand = list(roots) + [t]
            # doyen : un groupe ancre, sinon le plus ancien (naissance la plus basse)
            elder = min(cand, key=lambda r: (0 if anch[r] else 1, birth[r], r))
            for r in cand:
                if r == elder:
                    continue
                if (not anch[r]) and (pt - birth[r] >= persistence) and len(members[r]) >= min_tets:
                    hang[members[r]] = next_id
                    next_id += 1
                else:
                    members[elder].extend(members[r])
                parent_uf[r] = elder
                members[r] = []
            anch[elder] = any(anch[r] for r in cand)
            if t != elder and parent_uf[t] != elder:
                parent_uf[t] = elder
    return hang


def _reeb_labels(tet_mesh, adjacency, platform_nodes, band, centers):
    """Squelette topologique : branches du graphe de Reeb de la distance geodesique au plateau.

    Les tetraedres sont ranges par tranches de distance (largeur `band`) ; dans chaque tranche, les
    composantes connexes sont les noeuds du graphe ; deux composantes de tranches voisines qui se
    touchent sont reliees. Une BRANCHE est une chaine maximale de composantes qui n'ont qu'un
    successeur et dont le successeur n'a qu'un predecesseur : elle s'arrete a chaque separation
    (un bras qui part d'un pilier, un crochet qui descend d'un bras) et a chaque jonction.
    Rend (label par tetraedre, direction par tetraedre = gradient geodesique moyen de sa branche)."""
    tets = tet_mesh.tetrahedra
    T = len(tets)
    d = s4._geodesic_distance_to_platform(tet_mesh, platform_nodes)
    finite = np.isfinite(d)
    if not finite.all():
        d = np.where(finite, d, (d[finite].max() if finite.any() else 0.0) + band)
    f = d[tets].mean(axis=1)
    coo = sp.triu(adjacency, k=1).tocoo()
    a, b = coo.row, coo.col
    band = max(band, REEB_BAND_TET_FACTOR * float(np.median(np.linalg.norm(centers[a] - centers[b], axis=1))))
    k = np.floor(f / band).astype(np.int64)
    # composantes de chaque tranche ; les petits morceaux (decoupe irreguliere des tranches dans le
    # maillage) sont fondus dans leur voisin, sinon ils cassent les chaines en faux embranchements
    comp = _branch_labels(tet_mesh, lambda x, y: k[x] == k[y], adjacency, REEB_MIN_COMPONENT_TETS)
    C = int(comp.max()) + 1
    fc = np.bincount(comp, weights=f, minlength=C) / np.maximum(np.bincount(comp, minlength=C), 1)
    ca, cb = comp[a], comp[b]
    diff = ca != cb
    ca, cb = ca[diff], cb[diff]
    lo = np.where(fc[ca] < fc[cb], ca, cb)
    hi = np.where(fc[ca] < fc[cb], cb, ca)
    pairs = np.unique(np.stack([lo, hi], axis=1), axis=0) if len(lo) else np.zeros((0, 2), dtype=np.int64)
    n_succ = np.bincount(pairs[:, 0], minlength=C)
    n_pred = np.bincount(pairs[:, 1], minlength=C)
    chain = (n_succ[pairs[:, 0]] == 1) & (n_pred[pairs[:, 1]] == 1)
    succ = np.full(C, -1, dtype=np.int64)
    succ[pairs[chain, 0]] = pairs[chain, 1]
    has_pred = np.zeros(C, dtype=bool)
    has_pred[pairs[chain, 1]] = True
    vol_c = np.bincount(comp, minlength=C).astype(float)
    ctr = np.stack([np.bincount(comp, weights=centers[:, i], minlength=C) for i in range(3)], axis=1) / vol_c[:, None]
    # parcours de chaque chaine depuis sa tete ; nouvelle branche a chaque coude
    branch_of_comp = np.full(C, -1, dtype=np.int64)
    seg_dir = []
    cos_bend = np.cos(np.radians(REEB_BEND_DEG))
    nb = 0
    for head in np.where(~has_pred)[0]:
        c = int(head)
        branch_of_comp[c] = nb
        start = ctr[c]
        seg_dir.append(np.zeros(3))
        while succ[c] >= 0:
            n = int(succ[c])
            step = ctr[n] - ctr[c]
            sn = np.linalg.norm(step)
            cur = ctr[c] - start
            cn = np.linalg.norm(cur)
            if sn > 1e-9 and cn > 1e-9 and (step @ cur) / (sn * cn) < cos_bend:
                seg_dir[nb] = cur / cn
                nb += 1
                seg_dir.append(np.zeros(3))
                start = ctr[c]
            branch_of_comp[n] = nb
            c = n
        cur = ctr[c] - start
        cn = np.linalg.norm(cur)
        seg_dir[nb] = cur / cn if cn > 1e-9 else np.zeros(3)
        nb += 1
    for c in np.where(branch_of_comp < 0)[0]:       # composantes en cycle (rare) : branche propre
        branch_of_comp[c] = nb
        seg_dir.append(np.zeros(3))
        nb += 1
    label = branch_of_comp[comp].astype(np.int64)
    # direction : de la tete a la queue de la branche ; a defaut (branche d'une seule tranche),
    # gradient geodesique moyen (pondere par le volume)
    shape_gradients, volumes = s4._compute_shape_gradients(tet_mesh)
    grad = s4._compute_per_tet_gradient(tet_mesh, shape_gradients, d)
    gn = np.linalg.norm(grad, axis=1)
    unit = np.where(gn[:, None] > 1e-12, grad / np.maximum(gn, 1e-12)[:, None], 0.0)
    acc = np.zeros((nb, 3))
    np.add.at(acc, label, unit * volumes[:, None])
    sd = np.array(seg_dir)
    use = np.linalg.norm(sd, axis=1) > 1e-9
    acc[use] = sd[use]
    return label, acc


def build_planar_parts(tet_mesh: TetMesh, mesh_data: MeshData, rdg: RegionDependencyGraph, params,
                       overhang_map=None) -> tuple[TetMesh, ScalarField]:
    global LAST_DIAGNOSTIC
    diag = PlanarPartsDiagnostic()
    T = tet_mesh.tet_count
    tets = tet_mesh.tetrahedra
    verts = tet_mesh.tet_vertices
    centers = verts[tets].mean(axis=1)
    N = tet_mesh.node_count

    _, nearest_face = mesh_data.kdtree.query(centers)
    nearest_face = np.asarray(nearest_face, dtype=np.int64)
    face_region = np.full(mesh_data.face_count, -1, dtype=np.int64)
    face_dir = np.zeros((mesh_data.face_count, 3))
    face_dir[:, 2] = 1.0
    for rid, region in rdg.regions.items():
        face_region[region.face_ids] = rid
        face_dir[region.face_ids] = region.target_direction
    tet_region = face_region[nearest_face]
    tet_dir = face_dir[nearest_face]

    adjacency = s4._build_tet_tet_adjacency(tet_mesh)
    min_tets = int(params.multibranch_min_branch_tets)
    platform_nodes = s4._identify_platform_nodes(tet_mesh, mesh_data, rdg)
    if getattr(params, "planar_parts_partition", "region") == "reeb":
        reeb, reeb_dir = _reeb_labels(tet_mesh, adjacency, platform_nodes, float(params.planar_reeb_band_mm), centers)
        label = _branch_labels(tet_mesh, lambda a, b: reeb[a] == reeb[b], adjacency, min_tets)
        # direction de la partie : celle de la branche de Reeb majoritaire (en volume) qu'elle contient
        tet_dir = reeb_dir[reeb]
    else:
        label = _branch_labels(tet_mesh, lambda a, b: tet_region[a] == tet_region[b], adjacency, min_tets)
    is_platform = np.zeros(N, dtype=bool)
    is_platform[platform_nodes] = True
    up = np.asarray(params.up_vector if params.up_vector is not None else [0.0, 0.0, 1.0], dtype=float)
    diag.initial_parts = int(label.max()) + 1

    B, bg, branch_tets, branch_nodes, roots, depth = _analyze(label, tets, verts, adjacency, is_platform)
    base_dir = [tet_dir[ids].mean(axis=0) if len(ids) else up for ids in branch_tets]
    part_dir = _part_directions(B, bg, branch_nodes, roots, depth, verts, is_platform, base_dir, up,
                                params.radial_tilt_origin)
    depth = np.asarray(depth, dtype=np.int64)
    parent = np.full(B, -1, dtype=np.int64)
    for round_ in range(MAX_SPLIT_ROUNDS):
        anchored = _anchored_tets(label, depth, parent, tets, adjacency, is_platform)
        hang = _hanging_tets(label, branch_tets, part_dir, centers, adjacency, anchored, min_tets,
                             PLANAR_MIN_PERSISTENCE_MM)
        ids = np.where(hang >= 0)[0]
        if len(ids) == 0:
            break
        comp = hang
        # une partie fille par (partie d'origine, composante) : parent explicite, direction opposee
        pairs = np.unique(np.stack([label[ids], comp[ids]], axis=1), axis=0)
        new_label = label.copy()
        new_dirs, new_depth, new_parent = [], [], []
        for k, (par, c) in enumerate(pairs):
            sel = ids[(label[ids] == par) & (comp[ids] == c)]
            new_label[sel] = B + k
            d = -part_dir[par]
            new_dirs.append(d)
            new_depth.append(int(depth[par]) + 1)
            new_parent.append(int(par))
        label = new_label
        part_dir = np.vstack([part_dir, np.array(new_dirs)])
        depth = np.concatenate([depth, new_depth])
        parent = np.concatenate([parent, new_parent])
        B += len(pairs)
        branch_tets = [np.where(label == b)[0] for b in range(B)]
        diag.split_rounds += 1
        diag.hanging_tets.append(int(len(ids)))
    branch_nodes = [np.unique(tets[ids]) for ids in branch_tets]
    if params.radial_tilt_origin is not None and B > 0:
        ctr = np.array([verts[branch_nodes[b]].mean(axis=0) if len(branch_nodes[b]) else np.zeros(3)
                        for b in range(B)])
        proxy = TetMesh(tet_vertices=ctr, tetrahedra=np.tile(np.arange(B)[:, None], (1, 4)))
        part_dir, _ = s4._project_to_reachable_radial_tilt(proxy, part_dir, params.radial_tilt_origin, up)

    diag.final_parts = B
    LAST_DIAGNOSTIC = diag

    new_verts, new_tets, values, ranges = [], np.zeros_like(tets), [], []
    directions = np.zeros((T, 3))
    base = 0
    for b in range(B):
        ids = branch_tets[b]
        nodes = branch_nodes[b]
        lut = np.full(N, -1, dtype=np.int64)
        lut[nodes] = base + np.arange(len(nodes))
        new_tets[ids] = lut[tets[ids]]
        pts = verts[nodes]
        d = part_dir[b]
        proj = pts @ d
        values.append(proj - proj.min())
        ranges.append(float(proj.max() - proj.min()))
        new_verts.append(pts)
        directions[ids] = d
        base += len(nodes)
    order = sorted(range(B), key=lambda b: (depth[b], b))
    start = np.zeros(B)
    acc = 0.0
    for b in order:
        start[b] = acc
        acc += ranges[b] + 1.0
    for b in range(B):
        values[b] = values[b] + start[b]
    new_mesh = TetMesh(tet_vertices=np.vstack(new_verts), tetrahedra=new_tets)
    field_ = ScalarField(values=np.concatenate(values), gradient=directions)
    field_.branch_labels = label
    field_.branch_depths = {int(b): int(depth[b]) for b in range(B)}
    return new_mesh, field_
