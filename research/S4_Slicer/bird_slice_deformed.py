"""
Tranchage des couches Bird DANS L'ESPACE DEFORME, puis retour de la
trajectoire dans la piece reelle (plan approuve le 2026-09-26).

Pourquoi : l'etape 6 aplatit chaque couche sur un plan pour y calculer
parois et remplissage. Les couches Bird sont tres courbes dans la piece
reelle (189/284 a plus de 3 mm de leur plan sur le Benchy, 65 qui se
recouvrent une fois aplaties) : l'approximation casse. Dans l'espace
deforme, en revanche, chaque couche est EXACTEMENT un plan z = constante :
les etapes 5 et 6 y travaillent sans approximation. C'est l'astuce de Bird
(Cura tranche des plans), sans Cura.

Retour : un point de la trajectoire est localise dans son tetraedre deforme
(coordonnees barycentriques) et replace au meme endroit relatif dans le
tetraedre d'origine. La deformation etant lineaire par tetraedre, c'est
exact, sauf dans les zones de plissement ou deux tetraedres deformes se
recouvrent (le point est alors ambigu, et compte comme tel).

Ponts DESACTIVES dans cette version : un pont dans l'espace deforme n'a pas
de sens physique. A reevaluer dans la piece reelle.

A lancer avec le venv d'APEX :
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe bird_slice_deformed.py <modele> <passe> [variante, ex. _pin]
      [--skin] [--skin-thickness=<mm>] [--first-layer=<motif>] [--skin-exclude=<json de P5>]
  --skin : corps CREUSE sous les zones de lissage non-planaire retenues par
           surface_skin.py (loop_*_skin_pass<k>.npz), de l'epaisseur des
           couches de lissage (--skin-thickness, 0,2 mm par defaut) ;
  --first-layer : premiere couche pleine avec ce motif (surface_patterns,
           SANS valeur par defaut : sans l'option, l'etape 6 la traite comme
           les autres). Sortie suffixee _skin avec --skin.
Sortie : output_models/loop_<modele>_deformed_slice_pass<k>.pkl, lisible par
render_toolpaths.py (toolpath_2d dans le plan de la couche, toolpath_3d et
contour_3d dans la piece reelle, repere du STL).
"""
import os
import pickle
import sys
import time

import numpy as np
import trimesh
from scipy.spatial import cKDTree
from shapely.geometry import Polygon

from bird_to_apex import barycentric, frame_offset, tet_gradients

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)

from core.mesh_types import OverhangMap, Region, RegionDependencyGraph, ScalarField, TetMesh  # noqa: E402
from main import PipelineParams  # noqa: E402
from pipeline.step1_ingest import _build_mesh_data  # noqa: E402
from pipeline.step5_layers import run as layers_run  # noqa: E402
from pipeline.step6_toolpath2d import _classify_contours, _lift_to_3d, _project_to_local_2d, run as toolpath2d_run  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
_OPPOSITE = [((1, 2, 3), 0), ((0, 2, 3), 1), ((0, 1, 3), 2), ((0, 1, 2), 3)]
# Longueur maximale d'un segment de trajectoire avant son retour dans la piece
# reelle (voir densify) : bien plus court qu'une arete de tetraedre.
MAX_SEGMENT_MM = 0.5


def boundary_faces(cells, reference_points):
    """Faces de tetraedres qui n'apparaissent qu'une fois = la surface.
    Orientees vers l'exterieur a l'aide du sommet oppose, evalue sur la piece
    d'ORIGINE (toujours valide, contrairement a l'espace deforme ou quelques
    tetraedres se retournent)."""
    tri = np.concatenate([cells[:, list(f)] for f, _ in _OPPOSITE])
    opp = np.concatenate([cells[:, o] for _, o in _OPPOSITE])
    keys = np.sort(tri, axis=1)
    _, first, counts = np.unique(keys, axis=0, return_index=True, return_counts=True)
    keep = first[counts == 1]
    tri, opp = tri[keep], opp[keep]
    p = reference_points
    normal = np.cross(p[tri[:, 1]] - p[tri[:, 0]], p[tri[:, 2]] - p[tri[:, 0]])
    inward = np.einsum("ij,ij->i", normal, p[opp] - p[tri[:, 0]]) > 0
    tri[inward] = tri[inward][:, [0, 2, 1]]
    return tri


class DeformedToOriginal:
    """Ramene des points de l'espace deforme vers la piece d'origine."""

    def __init__(self, cells, deformed, original, k=24):
        self.cells, self.deformed, self.original, self.k = cells, deformed, original, k
        self.tets = deformed[cells]
        self.tree = cKDTree(self.tets.mean(axis=1))

    def __call__(self, points):
        _, cand = self.tree.query(points, k=self.k)
        best_score = np.full(len(points), -np.inf)
        best_bary = np.zeros((len(points), 4))
        best_tet = np.zeros(len(points), dtype=np.int64)
        containing = np.zeros(len(points), dtype=np.int64)
        for j in range(self.k):
            t = cand[:, j]
            with np.errstate(all="ignore"):
                bary = barycentric(points, self.tets[t])
            score = np.nan_to_num(bary.min(axis=1), nan=-np.inf)
            # STRICTEMENT a l'interieur : un point sur la face commune de deux
            # tetraedres voisins est "dans" les deux, ce n'est pas une ambiguite.
            # Seul un recouvrement (plissement) place un point au coeur de deux.
            containing += score > 1e-6
            better = score > best_score
            best_score[better], best_bary[better], best_tet[better] = score[better], bary[better], t[better]
        mapped = np.einsum("ij,ijk->ik", best_bary, self.original[self.cells[best_tet]])
        return mapped, best_tet, best_score, containing

    def map_polyline(self, points, outside_margin=0.1):
        """
        Comme __call__, mais pour une POLYLIGNE (points consecutifs proches) :
        le tetraedre de chaque point est choisi pour que le chemin ramene soit
        le plus CONTINU possible dans la piece reelle.

        Dans un plissement, l'espace deforme se recouvre : un point y est au
        coeur de deux tetraedres venus de deux epaisseurs differentes de la
        piece. Choisir point par point (le mieux contenu) fait alors sauter la
        trajectoire d'une epaisseur a l'autre (sauts de 2 a 15 mm mesures).
        Ici, candidats = tetraedres qui contiennent le point (ou, pour un
        point hors de tout tetraedre, ceux a `outside_margin` pres du moins
        exterieur) ; on retient la suite de candidats qui minimise la somme
        des distances entre points ramenes consecutifs (Viterbi, exact).
        """
        n = len(points)
        _, cand = self.tree.query(points, k=self.k)
        with np.errstate(all="ignore"):
            bary = np.stack([barycentric(points, self.tets[cand[:, j]]) for j in range(self.k)], axis=1)
        score = np.nan_to_num(bary.min(axis=2), nan=-np.inf)                     # (n,k)
        mapped = np.einsum("njv,njvk->njk", bary, self.original[self.cells[cand]])  # (n,k,3)
        valid = score > -1e-6
        none = ~valid.any(axis=1)
        if none.any():
            best = score[none].max(axis=1, keepdims=True)
            valid[none] = score[none] >= best - outside_margin
        containing = (score > 1e-6).sum(axis=1)

        cost = np.where(valid[0], 0.0, np.inf)
        back = np.zeros((n, self.k), dtype=np.int64)
        for i in range(1, n):
            step = np.linalg.norm(mapped[i - 1][:, None, :] - mapped[i][None, :, :], axis=2)  # (k,k)
            total = cost[:, None] + step
            back[i] = total.argmin(axis=0)
            cost = np.where(valid[i], total.min(axis=0), np.inf)
        choice = np.zeros(n, dtype=np.int64)
        choice[-1] = int(cost.argmin())
        for i in range(n - 1, 0, -1):
            choice[i - 1] = back[i, choice[i]]
        rows = np.arange(n)
        return mapped[rows, choice], cand[rows, choice], score[rows, choice], containing


def material_polygons(loops_2d):
    """Polygones de matiere d'une couche (parite des boucles imbriquees)."""
    polys = []
    for outer, holes in _classify_contours(loops_2d):
        g = Polygon(outer, holes=holes).buffer(0)
        polys.extend(g.geoms if g.geom_type == "MultiPolygon" else ([g] if not g.is_empty else []))
    return polys


def main():
    flags = {a.split("=")[0]: (a.split("=", 1)[1] if "=" in a else True)
             for a in sys.argv[1:] if a.startswith("--")}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    model_name = args[0] if args else "benchy upsidedown tilted"
    pass_index = int(args[1]) if len(args) > 1 else 4
    # variante de la boucle : "" (d'origine), "_pin" (base figee), "_a0.5_pin"...
    variant = args[2] if len(args) > 2 else ""
    use_skin = "--skin" in flags
    skin_thickness = float(flags.get("--skin-thickness", 0.2))
    first_layer_pattern = flags.get("--first-layer")
    stl = os.path.join(HERE, "input_models", f"{model_name}.stl")
    if not os.path.exists(stl):                  # comme bird_loop.model_path
        stl = os.path.join(ROOT, "res", "test_models", f"{model_name}.stl")
    prefix = os.path.join(HERE, "output_models", "loop_" + model_name.replace(" ", "_") + variant)
    params = PipelineParams()
    w, density, walls_n = params.toolpath2d.nozzle_width, params.toolpath2d.infill_density, params.toolpath2d.wall_count

    cells = np.load(f"{prefix}_cells.npy")
    p0 = np.load(f"{prefix}_pass0.npy")
    pk = np.load(f"{prefix}_pass{pass_index}.npy")
    offset = frame_offset(stl, prefix)      # repere Bird (recentre) -> repere du STL

    # --- maillage deforme, surface, champ -------------------------------------
    faces = boundary_faces(cells, p0)
    surface = trimesh.Trimesh(pk, faces, process=False)
    surface.remove_unreferenced_vertices()
    mesh_data = _build_mesh_data(surface)
    tet_mesh = TetMesh(tet_vertices=pk, tetrahedra=cells)

    z = pk[:, 2].copy()
    with np.errstate(all="ignore"):
        grad_check = tet_gradients(pk, cells, z)
    ok = np.isfinite(grad_check).all(axis=1)
    dev = np.abs(grad_check[ok] - [0.0, 0.0, 1.0]).max()
    print(f"controle gradient : {ok.sum()}/{len(cells)} tetraedres calculables, "
          f"ecart max a (0,0,1) = {dev:.1e}")
    # Champ lineaire en z deforme : son gradient vaut exactement (0,0,1)
    # partout, y compris sur les tetraedres ecrases ou le calcul deborde.
    gradient = np.tile([0.0, 0.0, 1.0], (len(cells), 1))
    scalar_field = ScalarField(values=z, gradient=gradient)

    graph = RegionDependencyGraph(
        regions={0: Region(id=0, face_ids=np.arange(mesh_data.face_count),
                           target_direction=np.array([0.0, 0.0, 1.0]))},
        dependencies=[],
    )
    m = mesh_data.face_count
    no_bridges = OverhangMap(risk_score=np.zeros(m), is_overhang=np.zeros(m, bool),
                             is_bridge=np.zeros(m, bool), bridge_span_direction=np.full((m, 3), np.nan))

    t0 = time.time()
    layers = layers_run(tet_mesh, mesh_data, scalar_field, graph, params.layers)
    if use_skin:
        import skin_carve
        skin = np.load(f"{prefix}_skin_pass{pass_index}.npz")
        if "--skin-exclude" in flags:
            # zones retirees par P5 (collision avec la vraie tete) : leur
            # creusement est annule -- la surface y reste en couches normales
            import json
            excluded = json.load(open(flags["--skin-exclude"]))["excluded_zones"]
            skin = {key: skin[key] for key in skin.files}
            skin["zone"] = np.where(np.isin(skin["zone"], excluded), -1, skin["zone"])
            print(f"zones retirees (non creusees) : {len(excluded)}")
        psi = skin_carve.skin_psi(p0, skin, skin_thickness)
        isovalues = [L.isovalue for L in layers]
        if len(set(isovalues)) != len(isovalues):
            raise SystemExit("plusieurs couches par isovaleur : creusement non gere dans ce cas")
        todo = [i for i, L in enumerate(layers) if skin_carve.layer_is_touched(cells, pk, psi, L.isovalue)]
        touched = len(todo)
        workers = int(flags.get("--workers", max(1, min(10, (os.cpu_count() or 2) - 2))))
        if workers <= 1 or len(todo) < 2:
            exact = skin_carve.ExactDistance(skin, cells, pk, p0, skin_thickness)
            for i in todo:
                region = skin_carve.carved_region(cells, pk, psi, layers[i].isovalue, exact=exact,
                                                  step=skin_thickness)
                layers[i].contour_3d = skin_carve.region_loops(region, layers[i].isovalue)
        else:
            # Couches independantes une fois psi calcule (refonte P3 du
            # 2026-09-30) : un processus par coeur, chacun construit sa propre
            # ExactDistance une fois ; resultat identique a la boucle serie.
            from concurrent.futures import ProcessPoolExecutor
            skin_d = {key: np.asarray(skin[key]) for key in (skin.files if hasattr(skin, "files") else skin)}
            with ProcessPoolExecutor(workers, initializer=skin_carve.carve_worker_init,
                                     initargs=(skin_d, cells, pk, p0, psi, skin_thickness)) as pool:
                iso = [layers[i].isovalue for i in todo]
                for i, loops in zip(todo, pool.map(skin_carve.carve_worker_layer, iso,
                                                   chunksize=max(1, len(todo) // (4 * workers)))):
                    layers[i].contour_3d = loops
            print(f"creusement parallele : {workers} processus")
        print(f"creusement sous les zones de surface (epaisseur {skin_thickness} mm) : "
              f"{touched}/{len(layers)} couches recalculees ({time.time() - t0:.0f} s)")
    layers = toolpath2d_run(layers, mesh_data, no_bridges, params.toolpath2d)
    if first_layer_pattern:
        import skin_carve
        first = min(layers, key=lambda L: L.isovalue)
        loops = [_project_to_local_2d(np.asarray(c), first.local_frame) for c in first.contour_3d if len(c) >= 4]
        paths = skin_carve.first_layer_toolpath(
            loops, first_layer_pattern, params.toolpath2d.nozzle_width,
            params.toolpath2d.wall_count, params.toolpath2d.degenerate_area_ratio)
        first.toolpath_2d = [np.column_stack([p, np.zeros(len(p))]) for p in paths]
        print(f"premiere couche : motif {first_layer_pattern}, {len(paths)} passes")
    print(f"etapes 5+6 dans l'espace deforme : {len(layers)} couches ({time.time() - t0:.0f} s)")

    # --- verification 1 : planeite, recouvrement ------------------------------
    spread = np.array([np.ptp(np.vstack(L.contour_3d)[:, 2]) for L in layers if L.contour_3d])
    from collections import Counter
    from shapely.validation import explain_validity
    reasons, invalid = Counter(), 0
    for L in layers:
        bad = [explain_validity(Polygon(c[:, :2])) for c in L.contour_3d if len(c) >= 3]
        bad = [r.split("[")[0].strip() for r in bad if r != "Valid Geometry"]
        invalid += bool(bad)
        reasons.update(bad)
    print(f"planeite : ecart en z max par couche = {spread.max():.1e} mm ; "
          f"couches avec un contour invalide : {invalid}/{len(layers)} (avant : 65/284) ; "
          f"raisons : {dict(reasons)}")

    # --- verification 2 : remplissage mesure / attendu ------------------------
    ratios, split_pieces, empty_pieces, holes_total, holes_bordered = [], 0, 0, 0, 0
    for L in layers:
        # dans le REPERE LOCAL de la couche, celui de toolpath_2d (origine et
        # axes propres a chaque couche, meme si le plan est horizontal)
        loops = [_project_to_local_2d(np.asarray(c), L.local_frame) for c in L.contour_3d if len(c) >= 3]
        polys = material_polygons(loops) if loops else []
        area = sum(p.area for p in polys)
        tp = L.toolpath_2d or []
        lengths = [np.linalg.norm(np.diff(t[:, :2], axis=0), axis=1).sum() for t in tp]
        # remplissage = polylignes ouvertes qui ne sont pas un contour brut
        # (les parois offsettees sont des anneaux fermes ; depuis le zigzag,
        # le remplissage n'est plus reconnaissable a ses 2 points)
        raw = cKDTree(np.vstack(loops)) if loops else None
        is_infill = [
            len(t) >= 2 and not np.allclose(t[0, :2], t[-1, :2])
            and (raw is None or raw.query(t[:, :2])[0].max() > 1e-6)
            for t in tp
        ]
        infill = sum(l for f, l in zip(is_infill, lengths) if f)
        walls = sum(lengths) - infill
        expected = max(0.0, area - walls * w) * density / w
        if expected > 20.0:          # couches ou du remplissage est attendu
            ratios.append(infill / expected)
        # Couverture (defauts corriges le 2026-09-26 dans _offset_walls) :
        #  - chaque morceau d'une zone coupee en deux par la 1re paroi doit
        #    contenir de la trajectoire ;
        #  - chaque trou doit etre borde par son contour.
        from shapely.geometry import MultiPoint
        tp_pts = MultiPoint(np.vstack([t[:, :2] for t in tp])) if tp else None
        for p in polys:
            g = p.buffer(-w, join_style="mitre")
            if g.geom_type == "MultiPolygon":
                for q in g.geoms:
                    if q.area > 0.5:
                        split_pieces += 1
                        empty_pieces += tp_pts is None or not q.buffer(1e-6).intersects(tp_pts)
        # trous : sur les boucles BRUTES, telles que l'etape 6 les recoit (la
        # reparation buffer(0) des boucles qui se croisent cree de faux trous
        # et deplace les sommets)
        if tp and loops:
            tree = cKDTree(np.vstack([t[:, :2] for t in tp]))
            for _outer, holes in _classify_contours(loops):
                for hole in holes:
                    holes_total += 1
                    holes_bordered += (tree.query(hole[:, :2])[0] < 1e-6).mean() >= 0.95
    ratios = np.array(ratios)
    if len(ratios):
        print(f"remplissage mesure / attendu, sur {len(ratios)} couches qui en attendent : "
              f"mediane {np.median(ratios):.2f}, 10e-90e centiles {np.percentile(ratios, 10):.2f}-"
              f"{np.percentile(ratios, 90):.2f} (couche 185 avant : 0,27)")
    else:
        print("remplissage mesure / attendu : aucune couche assez grande pour en attendre")
    print(f"morceaux (zone coupee par la 1re paroi, > 0,5 mm2) : {split_pieces}, "
          f"dont SANS trajectoire : {empty_pieces}")
    print(f"trous : {holes_total}, bordes par leur contour : {holes_bordered}")

    # --- retour dans la piece reelle ------------------------------------------
    # Une droite de l'espace deforme devient une COURBE dans la piece reelle
    # (la deformation n'est lineaire que par tetraedre). Ramener seulement les
    # extremites d'une ligne de remplissage (2 points) donnerait une corde qui
    # traverse la piece au lieu de suivre la couche : on subdivise d'abord.
    def densify(points, step=MAX_SEGMENT_MM):
        out = [points[:1]]
        for a, b in zip(points[:-1], points[1:]):
            n = max(1, int(np.ceil(np.linalg.norm(b - a) / step)))
            out.append(a + (b - a) * (np.arange(1, n + 1)[:, None] / n))
        return np.vstack(out)

    back = DeformedToOriginal(cells, pk, p0)

    def map_refined(points, rounds=6):
        """Ramene une polyligne deformee, puis recoupe en deux tout segment
        qui depasse 2 x MAX_SEGMENT_MM une fois ramene (zone que la
        deformation a comprimee : un petit pas deforme y couvre plusieurs mm
        de piece reelle) et recommence. Un segment qui ne raccourcit pas en
        le recoupant est un vrai saut (pli), laisse tel quel."""
        result = back.map_polyline(points)
        for _ in range(rounds):
            lengths = np.linalg.norm(np.diff(result[0], axis=0), axis=1)
            long = np.flatnonzero(lengths > 2 * MAX_SEGMENT_MM)
            if not len(long):
                break
            midpoints = (points[long] + points[long + 1]) / 2
            points = np.insert(points, long + 1, midpoints, axis=0)
            result = back.map_polyline(points)
        return (points, *result)

    def outside_distance(points, tets, score):
        """Distance (mm, espace deforme) de chaque point a son tetraedre :
        0 dedans ; sinon ecart au point le plus proche obtenu en ramenant ses
        coordonnees barycentriques dans le tetraedre (approximation par
        exces de la vraie distance, suffisante pour trier bord / dehors)."""
        dist = np.zeros(len(points))
        m = score < -1e-6
        if m.any():
            corners = back.tets[tets[m]]
            with np.errstate(all="ignore"):
                b = np.clip(barycentric(points[m], corners), 0.0, None)
            b /= b.sum(axis=1, keepdims=True)
            dist[m] = np.linalg.norm(points[m] - np.einsum("nv,nvk->nk", b, corners), axis=1)
        return dist

    out_layers, n_pts, n_out, n_amb, n_cut, n_jump, all_contour = [], 0, 0, 0, 0, 0, []
    for L in layers:
        tp3, tps = [], []
        for t in (L.toolpath_2d or []):
            refined, mapped, tet, score, containing = map_refined(
                densify(_lift_to_3d(np.asarray(t), L.local_frame)))
            score = outside_distance(refined, tet, score)
            n_pts += len(mapped)
            n_out += int((score > 0).sum())
            n_amb += int((containing >= 2).sum())
            # Un point a plus d'une demi-largeur de cordon HORS du volume
            # deforme n'est pas de la matiere : la couche l'a recu parce que son
            # contour se croise lui-meme dans un pli, et la reparation du
            # polygone y a cree une fausse zone pleine. Extrapole, il sortait
            # de la piece (jusqu'a 4 mm, couche 65 du Benchy). On coupe la
            # trajectoire a ces points au lieu de les deposer.
            keep = score <= w / 2
            n_cut += int((~keep).sum())
            edges = np.flatnonzero(np.diff(np.concatenate([[0], keep.astype(np.int8), [0]])))
            for start, stop in zip(edges[::2], edges[1::2]):
                run, run_score = mapped[start:stop], score[start:stop]
                # un segment que map_refined n'a pas su raccourcir est un vrai
                # saut (changement d'epaisseur dans un pli) : l'extruder
                # tracerait une droite a travers la piece -> trajet a vide
                jumps = np.flatnonzero(
                    np.linalg.norm(np.diff(run, axis=0), axis=1) > 2 * MAX_SEGMENT_MM) + 1
                n_jump += len(jumps)
                for piece, piece_score in zip(np.split(run, jumps), np.split(run_score, jumps)):
                    if len(piece) >= 2:
                        tp3.append(piece + offset)
                        tps.append(piece_score)
        contours = [back.map_polyline(np.asarray(c))[0] + offset for c in L.contour_3d]
        all_contour.extend(contours)
        out_layers.append({"contour_3d": contours, "isovalue": L.isovalue,
                           "toolpath_2d": [np.asarray(t) for t in (L.toolpath_2d or [])],
                           "toolpath_3d": tp3, "toolpath_score": tps})
    print(f"retour : {n_pts} points, hors de tout tetraedre {n_out} ({n_out / n_pts:.2%}), "
          f"ambigus (plissement) {n_amb} ({n_amb / n_pts:.2%}), coupes (a plus de w/2 hors du volume) {n_cut} ({n_cut / n_pts:.2%}), sauts transformes en trajets a vide {n_jump}")
    seg = np.concatenate([np.linalg.norm(np.diff(t, axis=0), axis=1)
                          for L in out_layers for t in L["toolpath_3d"] if len(t) > 1])
    print(f"segments ramenes : median {np.median(seg):.2f} mm, 99e centile "
          f"{np.percentile(seg, 99):.2f} mm, max {seg.max():.2f} mm (avant subdivision : max 39,2)")

    original = trimesh.load(stl)
    pts = np.vstack(all_contour)
    pts = pts[np.random.default_rng(0).choice(len(pts), min(20000, len(pts)), replace=False)]
    _, dist, _ = trimesh.proximity.closest_point(original, pts)
    print(f"contours ramenes -> surface du STL d'origine : mediane {np.median(dist):.3f} mm, "
          f"95e centile {np.percentile(dist, 95):.3f} mm, max {dist.max():.3f} mm")

    out = f"{prefix}_deformed_slice{'_skin' if use_skin else ''}_pass{pass_index}.pkl"
    with open(out, "wb") as f:
        pickle.dump(out_layers, f)
    print(f"sortie : {out}")


if __name__ == "__main__":
    main()
