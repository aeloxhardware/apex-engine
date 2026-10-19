"""
Etape 4b - CHAMP SCALAIRE MULTIBRANCHE (2026-10-19, opt-in : `FieldSolverParams.multibranch`).

Le defaut de fond du champ unique (limites #51 a #54) : un seul phi doit croitre depuis
le plateau ET croitre vers le bas sous un plafond suspendu ; le terme de lissage arbitre
en ecrasant le second, et `enforce_upward_growth` ne fait que masquer le conflit.

Ici le volume est decoupe en BRANCHES, resolues l'une apres l'autre dans l'ordre de leur
profondeur depuis le plateau :

  1. chaque tetraedre herite de la region de step3 de la face de surface la plus proche ;
  2. une branche = une composante connexe (par faces) d'une meme region ; les branches de
     moins de `multibranch_min_branch_tets` tetraedres sont fondues dans leur voisine ;
  3. le graphe des branches est parcouru en largeur depuis les branches au contact du
     plateau (profondeur = nombre de sauts) ;
  4. chaque branche resout SON systeme (lissage + « le gradient suit la direction de
     croissance »), avec pour Dirichlet DUR les noeuds deja valorises par les branches
     precedentes (interface avec le parent) ; une racine est ancree au plateau.

phi est donc continu a travers les interfaces, mais chaque branche croit en
S'ELOIGNANT de son attache : une coque suspendue croit vers le bas depuis le pilier,
sans que le plateau lui impose de monter.

Direction de croissance : `multibranch_direction_mode` = "region" (cible de step3, retournee
si elle pointe vers l'attache) ou "geodesic" (gradient de la distance geodesique depuis
l'attache, a l'interieur de la branche).

LIMITES (signalees, pas masquees) : les `Constraint[]` utilisateur ne sont pas honorees ;
`enforce_upward_growth` n'est pas applique ; la projection radiale atteignable
(`radial_tilt_origin`) l'est.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.csgraph as spgraph

from core.mesh_types import MeshData, RegionDependencyGraph, ScalarField, TetMesh
from pipeline import step4_field_solver as s4


@dataclass
class MultibranchDiagnostic:
    branch_count: int = 0
    max_depth: int = 0
    root_count: int = 0
    free_root_count: int = 0
    flipped_branch_count: int = 0
    geodesic_fallback_tets: int = 0
    branch_tet_counts: list[int] = field(default_factory=list)
    branch_depths: list[int] = field(default_factory=list)


LAST_DIAGNOSTIC: MultibranchDiagnostic = MultibranchDiagnostic()


def _branch_labels(tet_mesh: TetMesh, same_fn, adjacency: sp.csr_matrix,
                   min_tets: int) -> np.ndarray:
    """Composantes connexes selon `same_fn(a, b)` (aretes gardees), puis fusion des petites
    dans leur voisine."""
    T = tet_mesh.tet_count
    coo = sp.triu(adjacency, k=1).tocoo()
    a, b = coo.row, coo.col
    same = same_fn(a, b)
    g = sp.csr_matrix((np.ones(int(same.sum())), (a[same], b[same])), shape=(T, T))
    _n, label = spgraph.connected_components(g, directed=False)
    label = label.astype(np.int64)

    for _ in range(50):
        sizes = np.bincount(label)
        small = sizes < min_tets
        if not small.any():
            break
        la, lb = label[a], label[b]
        diff = la != lb
        pa = np.concatenate([la[diff], lb[diff]])
        pb = np.concatenate([lb[diff], la[diff]])
        if len(pa) == 0:
            break
        key = pa.astype(np.int64) * (label.max() + 1) + pb
        uniq, counts = np.unique(key, return_counts=True)
        ua, ub = uniq // (label.max() + 1), uniq % (label.max() + 1)
        # un voisin n'est retenu que s'il est au moins aussi gros (evite les echanges A<->B)
        ok = small[ua] & ((sizes[ub] > sizes[ua]) | ((sizes[ub] == sizes[ua]) & (ub < ua)))
        if not ok.any():
            break
        ua, ub, counts = ua[ok], ub[ok], counts[ok]
        order = np.lexsort((-counts, ua))
        ua, ub = ua[order], ub[order]
        first = np.concatenate([[True], ua[1:] != ua[:-1]])
        target = np.arange(label.max() + 1)
        target[ua[first]] = ub[first]
        label = target[label]
    _u, label = np.unique(label, return_inverse=True)
    return label.astype(np.int64)


def _local_graph(tet_mesh: TetMesh, tet_ids: np.ndarray) -> sp.csr_matrix:
    tets = tet_mesh.tetrahedra[tet_ids]
    pairs = np.vstack([tets[:, [0, 1]], tets[:, [0, 2]], tets[:, [0, 3]],
                       tets[:, [1, 2]], tets[:, [1, 3]], tets[:, [2, 3]]])
    lengths = np.linalg.norm(tet_mesh.tet_vertices[pairs[:, 0]] - tet_mesh.tet_vertices[pairs[:, 1]], axis=1)
    n = tet_mesh.node_count
    return sp.coo_matrix((np.concatenate([lengths, lengths]),
                          (np.concatenate([pairs[:, 0], pairs[:, 1]]),
                           np.concatenate([pairs[:, 1], pairs[:, 0]]))), shape=(n, n)).tocsr()


def _gradient_weights(tet_mesh, mesh_data, overhang_map, rdg, params, nearest_face):
    """Poids du terme de direction PAR TETRAEDRE, comme `build_static_system` : poids de base
    partout, majore autour des surplombs (germes propages sur `overhang_influence_depth`)."""
    effective_risk = np.asarray(overhang_map.risk_score, dtype=np.float64).copy()
    if rdg.platform_contact_faces is not None:
        effective_risk[np.asarray(rdg.platform_contact_faces, dtype=bool)] = 0.0
    tet_risk = effective_risk[nearest_face]
    bf_nodes, bf_tets, nearest_surface, face_pairs = s4._map_boundary_faces_to_surface(tet_mesh, mesh_data)
    risky = np.zeros(tet_mesh.tet_count, dtype=bool)
    if len(bf_tets) > 0:
        seeds = np.unique(bf_tets[effective_risk[nearest_surface] > 0.0])
        if len(seeds) > 0:
            risky = s4._spread_risk_influence(s4._build_tet_tet_adjacency(tet_mesh, face_pairs), seeds,
                                              params.overhang_influence_depth)
    w = np.full(tet_mesh.tet_count, params.gradient_baseline_weight, dtype=np.float64)
    w[risky] += params.gradient_overhang_weight * tet_risk[risky]
    return w


def solve_multibranch_field(tet_mesh: TetMesh, mesh_data: MeshData,
                            region_dependency_graph: RegionDependencyGraph,
                            params, overhang_map=None) -> ScalarField:
    global LAST_DIAGNOSTIC
    diag = MultibranchDiagnostic()
    T, N = tet_mesh.tet_count, tet_mesh.node_count
    tets = tet_mesh.tetrahedra
    verts = tet_mesh.tet_vertices
    shape_gradients, volumes = s4._compute_shape_gradients(tet_mesh)
    centers = verts[tets].mean(axis=1)

    # 1. region par tetraedre
    _, nearest_face = mesh_data.kdtree.query(centers)
    nearest_face = np.asarray(nearest_face, dtype=np.int64)
    face_region = np.full(mesh_data.face_count, -1, dtype=np.int64)
    face_dir = np.zeros((mesh_data.face_count, 3))
    face_dir[:, 2] = 1.0
    for rid, region in region_dependency_graph.regions.items():
        face_region[region.face_ids] = rid
        face_dir[region.face_ids] = region.target_direction
    tet_region = face_region[nearest_face]
    region_dirs = face_dir[nearest_face]
    if params.multibranch_blend_directions:
        up0 = np.asarray(params.up_vector if params.up_vector is not None else [0.0, 0.0, 1.0], dtype=float)
        region_dirs = s4._blend_tet_directions(tet_mesh, mesh_data, region_dependency_graph, up0, params)

    # 2. branches
    adjacency = s4._build_tet_tet_adjacency(tet_mesh)
    platform_nodes = s4._identify_platform_nodes(tet_mesh, mesh_data, region_dependency_graph)
    flow = None
    if params.multibranch_partition == "flow":
        dist = s4._heat_method_distance(tet_mesh, shape_gradients, volumes, platform_nodes, params)
        dist = np.where(np.isfinite(dist), dist, np.nanmax(dist[np.isfinite(dist)]) if np.isfinite(dist).any() else 0.0)
        g = s4._compute_per_tet_gradient(tet_mesh, shape_gradients, dist)
        gn = np.linalg.norm(g, axis=1)
        flow = np.zeros_like(g)
        okf = gn > 1e-9
        flow[okf] = g[okf] / gn[okf, None]
        # Quantification du flux sur la sphere (cases d'ouverture ~ angle) : une courbe
        # progressive (coude de 90 deg sur quelques tetraedres) est ainsi decoupee en
        # tronçons, ce qu'un seuil sur l'angle entre VOISINS ne ferait jamais (derive).
        from core.geometry_utils import fibonacci_sphere_directions
        ang = np.radians(max(params.multibranch_flow_angle_deg, 5.0))
        centers_dir = fibonacci_sphere_directions(max(int(4.0 * np.pi / (ang * ang)), 6))
        flow_bin = np.argmax(flow @ np.asarray(centers_dir).T, axis=1)
        flow_bin[~okf] = -1

        def same_fn(a, b):
            return okf[a] & okf[b] & (flow_bin[a] == flow_bin[b])
        region_dirs = np.where(okf[:, None], flow, region_dirs)
    else:
        def same_fn(a, b):
            return tet_region[a] == tet_region[b]
    label = _branch_labels(tet_mesh, same_fn, adjacency, int(params.multibranch_min_branch_tets))
    B = int(label.max()) + 1
    diag.branch_count = B

    # 3. graphe des branches, racines, profondeurs
    coo = sp.triu(adjacency, k=1).tocoo()
    la, lb = label[coo.row], label[coo.col]
    diff = la != lb
    bg = sp.csr_matrix((np.ones(int(diff.sum())), (la[diff], lb[diff])), shape=(B, B))
    bg = ((bg + bg.T) > 0).astype(float)

    is_platform = np.zeros(N, dtype=bool)
    is_platform[platform_nodes] = True

    branch_tets = [np.where(label == b)[0] for b in range(B)]
    branch_nodes = [np.unique(tets[ids]) for ids in branch_tets]
    roots = [b for b in range(B) if int(is_platform[branch_nodes[b]].sum()) >= 3]
    diag.root_count = len(roots)
    if roots:
        depth = spgraph.dijkstra(bg, directed=False, indices=roots, min_only=True, unweighted=True)
    else:
        depth = np.full(B, np.inf)
    free_roots = []
    # composantes deconnectees de tout plateau : on ancre leur point le plus bas
    while np.isinf(depth).any():
        cand = np.where(np.isinf(depth))[0]
        b0 = int(min(cand, key=lambda b: verts[branch_nodes[b], 2].min()))
        free_roots.append(b0)
        extra = spgraph.dijkstra(bg, directed=False, indices=[b0], min_only=True, unweighted=True)
        depth = np.where(np.isinf(depth), extra, depth)
        depth[b0] = 0.0
    diag.free_root_count = len(free_roots)
    all_roots = set(roots) | set(free_roots)
    order = sorted(range(B), key=lambda b: (depth[b], b))
    diag.max_depth = int(depth.max()) if B else 0
    diag.branch_tet_counts = [len(x) for x in branch_tets]
    diag.branch_depths = [int(depth[b]) for b in range(B)]

    # 4. noeuds fixes / inconnus par branche (topologiques : ne dependent pas des valeurs)
    assigned = np.zeros(N, dtype=bool)
    fixed_nodes: dict[int, np.ndarray] = {}
    anchor_values: dict[int, np.ndarray] = {}
    for b in order:
        nb = branch_nodes[b]
        if b in roots:
            fx = nb[assigned[nb] | is_platform[nb]]
        elif b in free_roots:
            low = nb[verts[nb, 2] <= verts[nb, 2].min() + 1e-9 * (np.ptp(verts[:, 2]) or 1.0)]
            fx = np.union1d(low, nb[assigned[nb]])
        else:
            fx = nb[assigned[nb]]
        fixed_nodes[b] = fx
        assigned[nb] = True

    # 5. directions de croissance par tetraedre
    risk_per_tet = None
    if overhang_map is not None:
        eff = np.asarray(overhang_map.risk_score, dtype=np.float64).copy()
        if region_dependency_graph.platform_contact_faces is not None:
            eff[np.asarray(region_dependency_graph.platform_contact_faces, dtype=bool)] = 0.0
        risk_per_tet = eff[nearest_face]
    mode = params.multibranch_direction_mode
    directions = region_dirs.copy()
    fallback = 0
    flipped = 0
    for b in order:
        ids = branch_tets[b]
        fx = fixed_nodes[b]
        if len(fx) == 0:
            continue
        anchor = verts[fx].mean(axis=0)
        d_b = region_dirs[ids]
        if mode == "geodesic":
            graph = _local_graph(tet_mesh, ids)
            dist = spgraph.dijkstra(graph, directed=False, indices=fx, min_only=True)
            nb = branch_nodes[b]
            ok_nodes = np.isfinite(dist)
            filled = np.where(ok_nodes, dist, np.nanmax(np.where(ok_nodes, dist, np.nan)) if ok_nodes.any() else 0.0)
            grad = np.einsum("tij,tj->ti", shape_gradients[ids], filled[tets[ids]])
            norm = np.linalg.norm(grad, axis=1)
            usable = (norm > 1e-9) & ok_nodes[tets[ids]].all(axis=1)
            d_b = d_b.copy()
            d_b[usable] = grad[usable] / norm[usable, None]
            fallback += int((~usable).sum())
            # le repli (region) est lui aussi oriente ci-dessous
        if params.multibranch_vertical_if_safe is not None and risk_per_tet is not None:
            if float((risk_per_tet[ids] > 0.0).mean()) < params.multibranch_vertical_if_safe:
                d_b = np.tile(np.asarray(params.up_vector if params.up_vector is not None
                                         else [0.0, 0.0, 1.0], dtype=float), (len(ids), 1))
                directions[ids] = d_b
                continue
        if mode == "up":      # signe comme `enforce_upward_growth` : isosurfaces identiques, sens montant
            up = np.asarray(params.up_vector if params.up_vector is not None else [0.0, 0.0, 1.0], dtype=float)
            flip = d_b @ up < 0.0
            d_b = d_b.copy()
            d_b[flip] *= -1.0
            directions[ids] = d_b
            continue
        if mode == "flow_mean":
            m = (d_b * volumes[ids, None]).sum(axis=0)
            nm = np.linalg.norm(m)
            if nm > 1e-9:
                d_b = np.tile(m / nm, (len(ids), 1))
            directions[ids] = d_b
            continue
        if mode == "flow":
            directions[ids] = d_b
            continue
        away = np.einsum("ij,ij->i", d_b, centers[ids] - anchor)
        if mode == "geodesic":
            bad = ~usable
            d_b[bad & (away < 0.0)] *= -1.0       # repli : oriente tetraedre par tetraedre
        elif float((away * volumes[ids]).sum()) < 0.0:
            d_b = -d_b
            flipped += 1
        directions[ids] = d_b
    diag.flipped_branch_count = flipped
    diag.geodesic_fallback_tets = fallback

    if params.multibranch_smoothing_iterations > 0:
        directions = s4._iterative_vector_smoothing(tet_mesh, directions, params.multibranch_smoothing_iterations,
                                                    params.geodesic_smoothing_step)
    if params.radial_tilt_origin is not None:
        up = params.up_vector if params.up_vector is not None else np.array([0.0, 0.0, 1.0])
        directions, _ = s4._project_to_reachable_radial_tilt(tet_mesh, directions, params.radial_tilt_origin, up)

    # 6. resolution branche par branche
    phi = np.full(N, np.nan)
    ws = float(params.laplacian_smoothness_weight)
    if overhang_map is not None and params.multibranch_use_risk_weights:
        tet_w = _gradient_weights(tet_mesh, mesh_data, overhang_map, region_dependency_graph, params, nearest_face)
    else:
        tet_w = np.full(T, float(params.gradient_baseline_weight))
    for b in order:
        ids = branch_tets[b]
        G = shape_gradients[ids]
        V = volumes[ids]
        d = directions[ids]
        nrm = np.linalg.norm(d, axis=1, keepdims=True)
        nrm[nrm == 0.0] = 1.0
        d = d / nrm
        uniq, inv = np.unique(tets[ids], return_inverse=True)
        loc = inv.reshape(-1, 4)
        fx = fixed_nodes[b]
        is_fixed = np.isin(uniq, fx)
        values_fixed = np.zeros(len(uniq))
        got = ~np.isnan(phi[uniq])
        values_fixed[got] = phi[uniq][got]
        if b in all_roots:
            plat = is_platform[uniq] & ~got
            values_fixed[plat] = params.platform_anchor_value
        # Une racine libre ancre ses noeuds les plus bas a 0
        if not is_fixed.any():
            continue
        # noeuds fixes sans valeur (ne devrait pas arriver) : 0
        wg = tet_w[ids]
        local = (V * (ws + wg))[:, None, None] * np.einsum("tki,tkj->tij", G, G)
        rhs_loc = (V * wg)[:, None] * np.einsum("tki,tk->ti", G, d)
        rows = np.repeat(loc, 4, axis=1).ravel()
        cols = np.tile(loc, (1, 4)).ravel()
        n = len(uniq)
        A = sp.coo_matrix((local.reshape(-1), (rows, cols)), shape=(n, n)).tocsr()
        rhs = np.zeros(n)
        np.add.at(rhs, loc.ravel(), rhs_loc.ravel())
        free = ~is_fixed
        if free.any():
            Auu = A[free][:, free]
            Auf = A[free][:, is_fixed]
            bu = rhs[free] - Auf @ values_fixed[is_fixed]
            xu = s4._solve_linear_system(Auu.tocsr(), bu, params)
            vals = values_fixed.copy()
            vals[free] = xu
        else:
            vals = values_fixed
        newly = np.isnan(phi[uniq])
        phi[uniq[newly]] = vals[newly]
        # les noeuds fixes issus du plateau gardent leur valeur d'ancrage
    phi = np.where(np.isnan(phi), 0.0, phi)
    LAST_DIAGNOSTIC = diag
    return ScalarField(values=phi, gradient=s4._compute_per_tet_gradient(tet_mesh, shape_gradients, phi),
                       branch_labels=label if params.multibranch_order_by_branch else None,
                       branch_depths={int(b): int(depth[b]) for b in range(B)} if params.multibranch_order_by_branch else None)
