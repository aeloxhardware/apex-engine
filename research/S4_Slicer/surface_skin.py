"""
Lissage non-planaire -- P2 (version 2, "P2-bis") : quelles faces de la surface
d'ORIGINE recoivent une couche de surface, en quels morceaux, a quel moment,
avec quelle orientation de tete.

Pour chaque face exterieure (hors contact plateau) :
  - MOMENT PORTEUR : rang de la couche qui porte le sommet le plus haut de la
    face (hauteur deformee phi) -- des cette couche finie, toute la matiere
    sous la face est deposee ;
  - ORIENTATION : la buse pointe vers -n (n = normale sortante). Angles
    d'axes par l'IK de l'etape 8 d'Apex, une fois par NORMALE DISTINCTE (les
    faces subdivisees d'une meme face d'origine partagent la leur) ;
  - DERNIER MOMENT ACCESSIBLE : la matiere imprimee ne fait que croitre
    (piece inter phi <= phi(t)), donc une face bloquee par l'objet a un
    moment l'est a tous les suivants. Un seul calcul par face : la tete
    placee au centre de la face, le plus petit phi de ses points tombant dans
    la piece donne la derniere couche avant blocage. Le plateau (z < 0) ne
    depend pas du moment. Tete placee pour des milliers de faces a la fois
    (matrice de rotation par face, meme convention que l'etape 9 : axe le
    plus interne d'abord, translation ensuite).
Une face est RETENUE si son dernier moment accessible >= son moment porteur.

Morceaux (decisions de l'utilisateur, 2026-09-27) :
  1. coupe aux ARETES VIVES : angle ENTRE LES FACES < EDGE_FACE_ANGLE_DEG
     (150 deg = 30 deg entre normales). Si le lissage des facettes est active
     (--facet-smoothing=<deg entre normales>), une arete sous ce seuil ne coupe
     plus ;
  2. BANDES DE MOMENT en nombre MINIMAL : chaque face a un intervalle
     [moment porteur, dernier moment accessible] ; le plus petit ensemble de
     moments qui touche tous les intervalles (glouton exact, par extremite
     droite croissante) ;
  3. frontieres de bandes LISSES : la surface est coupee le long des lignes
     de niveau phi = isovaleur des moments choisis (phi est lineaire sur
     chaque face : la coupe est exacte), au lieu de suivre les facettes ;
  4. chaque bande s'imprime AU PLUS TOT : au plus grand moment porteur de ses
     faces (toujours <= au dernier moment accessible de chacune).
Un morceau plus petit que MIN_ZONE_AREA est abandonne (et ne sera pas creuse).

Oracle : surface_skin_v1.py (version 1, coupe recursive en deux moities,
collision retestee a chaque decoupe ; 11 min 30 sur tilted benchy).

Repere : Bird (loop_*_pass*.npy), plateau = plan z = 0 (repere machine d'Apex).

Usage (venv d'Apex) :
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe -u surface_skin.py <modele> <passe> [variante]
        [--edge-angle=<deg entre faces>] [--facet-smoothing=<deg entre normales>] [--micro-inradius=<mm>]
        [--head-tolerance=<deg, defaut 45>]
Sortie : output_models/loop_<modele><variante>_skin_pass<k>.npz (maillage de
surface COUPE aux frontieres de bandes) + _skin_map_pass<k>.ply (statut) +
_skin_zones_pass<k>.ply (un morceau par couleur).
"""
import safe_io  # noqa: F401 -- ecritures robustes aux verrous transitoires
import os
import pickle
import sys
import time

import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import main as apex_main  # noqa: E402
from bird_slice_deformed import boundary_faces  # noqa: E402
from bird_to_apex import barycentric  # noqa: E402
from skin_carve import check_memory  # noqa: E402
import surface_paths as sp  # noqa: E402
from pipeline.step8_machine_toolpath import (  # noqa: E402
    Step8Params, _RawPoint, _resolve_all_axis_angles_coupled, _resolve_axis_angles,
)
from pipeline.step9_collision_sim import _filter_nozzle_contact_zone, _surface_sample_points  # noqa: E402

NOZZLE_WIDTH = 0.4
TIP_EXCLUSION = 2.0 * NOZZLE_WIDTH        # la buse touche la surface
HEAD_SAMPLE_SPACING = 5.0                 # comme l'etape 9 (limite #26)
MIN_ZONE_AREA = 4.0 * NOZZLE_WIDTH ** 2   # plus petit : pas de couche de surface
HEAD_TOLERANCE_DEG = 45.0                 # ecart maximal tete / normale (provisoire, decision utilisateur)
MICRO_INRADIUS = 1.5 * NOZZLE_WIDTH       # micro-zone : plus grande distance au bord en dessous
                                          # (trop etroite pour un cordon de bord + un motif)
MAX_EDGE = 5.0 * NOZZLE_WIDTH             # subdivision de la surface analysee
PLATE_CONTACT_TOL = 1e-3                  # faces posees sur le plateau (cf. bird_core)
EDGE_FACE_ANGLE_DEG = 150.0               # arete vive : angle ENTRE LES FACES en dessous
FACE_BATCH = 256                          # faces par paquet (memoire bornee)

STATUS = {0: "retenue", 1: "inatteignable", 2: "collision plateau", 3: "collision objet",
          4: "contact plateau", 5: "zone trop petite", 6: "micro-zone non lissee"}
STATUS_COLORS = {0: (40, 170, 70), 1: (150, 0, 150), 2: (230, 140, 0), 3: (220, 30, 30),
                 4: (120, 120, 120), 5: (90, 140, 230), 6: (0, 190, 200)}


class PrintedSolid:
    """La piece d'ORIGINE avec, en chaque point, sa hauteur deformee phi : la
    matiere imprimee au moment t est {phi <= phi(t)}. Localisation d'un point
    dans les tetraedres d'origine (k plus proches par centre)."""

    def __init__(self, cells, original, deformed, k=16):
        self.tets = original[cells]
        self.phi = deformed[cells][:, :, 2]
        self.tree = cKDTree(self.tets.mean(axis=1))
        self.k = k
        self.reach = float(np.linalg.norm(self.tets - self.tets.mean(axis=1, keepdims=True), axis=2).max())
        self.lo, self.hi = original.min(axis=0), original.max(axis=0)
        self.tolerance = -1e-9

    def phi_at(self, points):
        """phi aux points (NaN hors de la piece)."""
        out = np.full(len(points), np.nan)
        box = np.all((points >= self.lo) & (points <= self.hi), axis=1)
        idx_box = np.flatnonzero(box)
        if not len(idx_box):
            return out
        dist, cand = self.tree.query(points[idx_box], k=self.k)
        pending = dist[:, 0] <= self.reach
        for j in range(self.k):
            idx = np.flatnonzero(pending)
            if not len(idx):
                break
            t = cand[idx, j]
            with np.errstate(all="ignore"):
                bary = barycentric(points[idx_box[idx]], self.tets[t])
            inside = bary.min(axis=1) >= self.tolerance
            hit = idx[inside]
            out[idx_box[hit]] = np.einsum("nv,nv->n", bary[inside], self.phi[t[inside]])
            pending[hit] = False
        return out


def layer_rank(values, isovalues):
    """Rang de la couche qui porte la hauteur deformee `values`."""
    return np.clip(np.searchsorted(isovalues, values - 1e-9, side="left"), 0, len(isovalues) - 1)


def face_adjacency(faces):
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(faces)), 3)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    e, o = edges[order], owner[order]
    same = np.all(e[1:] == e[:-1], axis=1)
    a, b = o[:-1][same], o[1:][same]
    n = len(faces)
    return coo_matrix((np.ones(len(a) * 2), (np.r_[a, b], np.r_[b, a])), shape=(n, n)).tocsr()


def rotation_matrices(angles_deg, axes):
    """(F,3,3) : meme composition que step9._positioned_toolhead_vertices
    (rotations en ordre INVERSE de la liste d'axes, la plus interne d'abord)."""
    f = len(angles_deg)
    total = np.tile(np.eye(3), (f, 1, 1))
    for i, spec in reversed(list(enumerate(axes))):
        k = np.asarray(spec.axis, float)
        k = k / np.linalg.norm(k)
        kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        a = np.radians(angles_deg[:, i])[:, None, None]
        r = np.eye(3) + np.sin(a) * kx + (1 - np.cos(a)) * (kx @ kx)
        total = r @ total                     # appliquee APRES les plus internes
    return total


def last_accessible_moment(head, centers, angles, axes, solid, isovalues, face_ids, plate_margin=0.0):
    """(plateau, dernier_moment) par face de `face_ids`. dernier_moment = n-1
    si aucun point de tete ne tombe dans la piece ; -1 si bloquee des le
    debut."""
    n = len(isovalues)
    plate = np.zeros(len(face_ids), bool)
    last = np.full(len(face_ids), n - 1)
    for start in range(0, len(face_ids), FACE_BATCH):
        ids = face_ids[start:start + FACE_BATCH]
        check_memory(f"collision, paquet {start // FACE_BATCH}")
        rot = rotation_matrices(angles[ids], axes)
        pts = np.einsum("fij,pj->fpi", rot, head) + centers[ids][:, None, :]
        plate[start:start + len(ids)] = (pts[:, :, 2] < plate_margin).any(axis=1)
        phi = solid.phi_at(pts.reshape(-1, 3)).reshape(len(ids), -1)
        with np.errstate(all="ignore"):
            low = np.nanmin(np.where(np.isnan(phi), np.inf, phi), axis=1)
        inside = np.isfinite(low)
        # bloquee au moment m ssi min phi <= iso[m] + 1e-9 (meme test que la v1)
        last[start:start + len(ids)][inside] = np.searchsorted(isovalues, low[inside] - 1e-9, side="left") - 1
    return plate, last


_ACCESS = {}


def _access_init(head, centers, angles, axes, cells, original, deformed, isovalues, plate_margin):
    _ACCESS.update(head=head, centers=centers, angles=angles, axes=axes, isovalues=isovalues,
                   plate_margin=plate_margin, solid=PrintedSolid(cells, original, deformed))


def _access_chunk(face_ids):
    a = _ACCESS
    return last_accessible_moment(a["head"], a["centers"], a["angles"], a["axes"], a["solid"], a["isovalues"],
                                  face_ids, a["plate_margin"])


def last_accessible_moment_parallel(head, centers, angles, axes, cells, original, deformed, isovalues, face_ids,
                                    plate_margin=0.0, workers=None):
    """Meme resultat que last_accessible_moment, faces reparties sur plusieurs
    processus (faces independantes)."""
    workers = workers or max(1, min(10, (os.cpu_count() or 2) - 2))
    if workers <= 1 or len(face_ids) < 4 * FACE_BATCH:
        solid = PrintedSolid(cells, original, deformed)
        return last_accessible_moment(head, centers, angles, axes, solid, isovalues, face_ids, plate_margin)
    from concurrent.futures import ProcessPoolExecutor
    chunks = np.array_split(face_ids, workers * 4)
    with ProcessPoolExecutor(workers, initializer=_access_init,
                             initargs=(head, centers, angles, axes, cells, original, deformed, isovalues,
                                       plate_margin)) as pool:
        parts = list(pool.map(_access_chunk, chunks))
    return np.concatenate([p for p, _l in parts]), np.concatenate([l for _p, l in parts])


def minimal_stabbing(lo, hi):
    """Plus petit ensemble de moments touchant tous les intervalles [lo, hi]
    (glouton exact : extremite droite croissante)."""
    points = []
    for i in np.argsort(hi, kind="stable"):
        if not points or points[-1] < lo[i]:
            points.append(int(hi[i]))
    return np.array(points)


def split_along_level(vertices, faces, phi, face_mask, level):
    """Coupe les faces de `face_mask` traversees par phi = level (phi lineaire
    sur chaque face : coupe exacte). Sommets crees partages entre faces
    voisines (cle = arete). Rend (vertices, faces, phi, parent) ; parent =
    index de la face coupee d'origine pour chaque nouvelle face."""
    s = phi[faces] - level
    pos = s >= 0
    cross = face_mask & (pos.any(axis=1) & ~pos.all(axis=1))
    parent = np.arange(len(faces))
    if not cross.any():
        return vertices, faces, phi, parent
    ci = np.flatnonzero(cross)
    f, sp = faces[ci], pos[ci]
    # sommet isole en premier (celui dont le signe differe des deux autres), orientation conservee
    alone = np.where(sp.sum(axis=1) == 1, np.argmax(sp, axis=1), np.argmin(sp, axis=1))
    roll = np.stack([alone, (alone + 1) % 3, (alone + 2) % 3], axis=1)
    f = np.take_along_axis(f, roll, axis=1)
    a, b, c = f[:, 0], f[:, 1], f[:, 2]
    new_v, new_phi, key_index = [], [], {}

    def edge_point(p, q):
        out = np.empty(len(p), int)
        for j, (u, w) in enumerate(zip(p, q)):
            key = (min(u, w), max(u, w))
            if key not in key_index:
                su, sw = phi[u] - level, phi[w] - level
                t = su / (su - sw)
                if t <= 1e-9:
                    key_index[key] = u
                elif t >= 1 - 1e-9:
                    key_index[key] = w
                else:
                    key_index[key] = len(vertices) + len(new_v)
                    new_v.append(vertices[u] + t * (vertices[w] - vertices[u]))
                    new_phi.append(level)
            out[j] = key_index[key]
        return out

    pab, pac = edge_point(a, b), edge_point(a, c)
    tri = np.concatenate([np.stack([a, pab, pac], 1), np.stack([pab, b, c], 1), np.stack([pab, c, pac], 1)])
    par = np.tile(ci, 3)
    good = (tri[:, 0] != tri[:, 1]) & (tri[:, 1] != tri[:, 2]) & (tri[:, 0] != tri[:, 2])
    keep = np.flatnonzero(~cross)
    vertices = np.vstack([vertices, np.array(new_v).reshape(-1, 3)])
    phi = np.concatenate([phi, np.array(new_phi)])
    return vertices, np.vstack([faces[keep], tri[good]]), phi, np.concatenate([keep, par[good]])


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    model_name = args[0] if args else "benchy upsidedown tilted"
    k = int(args[1]) if len(args) > 1 else 4
    variant = args[2] if len(args) > 2 else ""
    edge_face_angle = float(flags.get("--edge-angle", EDGE_FACE_ANGLE_DEG))
    facet_smoothing = float(flags["--facet-smoothing"]) if "--facet-smoothing" in flags else None
    micro_inradius = float(flags.get("--micro-inradius", MICRO_INRADIUS))
    head_tolerance = float(flags.get("--head-tolerance", HEAD_TOLERANCE_DEG))
    trim_file = flags.get("--skin-trim")
    trim_radius = float(flags.get("--trim-radius", 2.0))
    prefix = os.path.join(HERE, "output_models", "loop_" + model_name.replace(" ", "_") + variant)
    timings = {}
    t0 = t_start = time.time()

    cells = np.load(f"{prefix}_cells.npy")
    p0 = np.load(f"{prefix}_pass0.npy")
    pk = np.load(f"{prefix}_pass{k}.npy")
    layers = pickle.load(open(f"{prefix}_deformed_slice_pass{k}.pkl", "rb"))
    isovalues = np.sort(np.unique([layer["isovalue"] for layer in layers]))

    surf_vertices, faces = trimesh.remesh.subdivide_to_size(p0, boundary_faces(cells, p0), MAX_EDGE)
    solid = PrintedSolid(cells, p0, pk)
    solid.tolerance = -1e-6                 # sommets de surface : sur la frontiere
    surf_phi = solid.phi_at(surf_vertices)
    solid.tolerance = -1e-9
    lost = np.isnan(surf_phi)
    if lost.any():                          # sommet non localise : plus proche sommet d'origine
        surf_phi[lost] = pk[cKDTree(p0).query(surf_vertices[lost])[1], 2]
    tri = surf_vertices[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    areas = 0.5 * np.linalg.norm(normals, axis=1)
    normals /= np.maximum(2 * areas, 1e-12)[:, None]
    centers = tri.mean(axis=1)
    n_faces = len(faces)
    status = np.zeros(n_faces, dtype=np.int8)
    status[tri[:, :, 2].max(axis=1) < p0[:, 2].min() + PLATE_CONTACT_TOL] = 4
    moment = layer_rank(surf_phi[faces].max(axis=1), isovalues)
    timings["surface et phi"] = time.time() - t0

    # --- orientation : IK de l'etape 8 --------------------------------------
    # Machine actuelle (limite #59 d'Apex) : C est IMPOSE par la position de la
    # face, B seul oriente la tete, dans le plan radial -> une IK par face (la
    # position compte). Machine a C libre : une IK par normale distincte.
    # Tolerance tete / normale : --head-tolerance (defaut provisoire 45 deg,
    # decision utilisateur du 2026-09-28, a ajuster).
    t0 = time.time()
    machine = apex_main._default_machine_profile()
    axes = machine.rotation_axes
    step8 = Step8Params()
    todo = np.flatnonzero(status == 0)
    angles = np.zeros((n_faces, len(axes)))
    residual = np.zeros(n_faces)
    if machine.position_coupled_axes:
        raw = [_RawPoint(position=centers[i], desired_direction=normals[i], layer_id=0, segment_index=0,
                         color_id=None, feedrate=1.0, extrusion=0.0, is_travel=False) for i in todo]
        a, err = _resolve_all_axis_angles_coupled(raw, axes, np.array([0.0, 0.0, 1.0]), machine, step8)
        angles[todo], residual[todo] = np.asarray(a), np.asarray(err)
        # Correctif du 2026-09-30 : C depend de la POSITION, donc l'ecart tete /
        # normale varie sur la face. Juge au seul centre, une face acceptee
        # pouvait depasser la tolerance a ses sommets (points > 45 deg vus en
        # P5). L'ecart retenu est le PIRE des trois sommets et du centre.
        corners = tri[todo].reshape(-1, 3)
        raw_c = [_RawPoint(position=q, desired_direction=normals[todo[j // 3]], layer_id=0, segment_index=0,
                           color_id=None, feedrate=1.0, extrusion=0.0, is_travel=False)
                 for j, q in enumerate(corners)]
        _a, err_c = _resolve_all_axis_angles_coupled(raw_c, axes, np.array([0.0, 0.0, 1.0]), machine, step8)
        residual[todo] = np.maximum(residual[todo], np.asarray(err_c).reshape(-1, 3).max(axis=1))
        ik_label = f"IK couplee a la position ({len(todo)} faces, centre + 3 sommets)"
    else:
        reference = np.array([0.0, 0.0, -1.0])
        uniq, inverse = np.unique(np.round(normals[todo], 9), axis=0, return_inverse=True)
        u_angles = np.zeros((len(uniq), len(axes)))
        u_resid = np.zeros(len(uniq))
        previous = None
        for i, n in enumerate(uniq):
            a, err = _resolve_axis_angles(-n, axes, reference, step8, previous)
            u_angles[i], u_resid[i], previous = a, err, a
        angles[todo], residual[todo] = u_angles[inverse.ravel()], u_resid[inverse.ravel()]
        ik_label = f"IK ({len(uniq)} normales distinctes)"
    status[(status == 0) & (np.degrees(residual) > head_tolerance)] = 1
    timings[ik_label] = time.time() - t0

    # --- dernier moment accessible (plateau, objet) ---------------------------
    t0 = time.time()
    head = _surface_sample_points(machine.nozzle_model, HEAD_SAMPLE_SPACING)
    # Memes regles que l'etape 9 et P5 (skin_print.py) : zone de contact de la
    # buse exclue (nozzle_contact_clearance, 3 mm) et MARGE de securite
    # (safety_margin, 1 mm) -- la tete est gonflee de la marge le long des
    # normales de son maillage (des deux cotes : maillage non etanche, sens des
    # normales non garanti) et le plateau est un plan releve de la marge. Sans
    # cela, P5 retirait 30 zones sur 48 que ce pre-calcul avait acceptees.
    sim = apex_main.PipelineParams().collision
    head = _filter_nozzle_contact_zone(head, sim.nozzle_contact_clearance)
    _c, _d, tri_id = trimesh.proximity.closest_point(machine.nozzle_model, head)
    hn = machine.nozzle_model.face_normals[tri_id]
    head = np.vstack([head, head + sim.safety_margin * hn, head - sim.safety_margin * hn])
    plate_margin = sim.safety_margin
    last = np.full(n_faces, -1)
    todo = np.flatnonzero(status == 0)
    # Ce calcul ne depend PAS des rognages de skin_loop : identique a chaque
    # tour (291 s sur Squirtle). Cache cle par une empreinte de TOUTES ses
    # entrees (2026-09-30) ; parallelise par paquets de faces sinon.
    import hashlib
    key = hashlib.sha256()
    for arr in (todo, centers[todo], angles[todo], head, isovalues, cells, p0, pk):
        key.update(np.ascontiguousarray(arr).tobytes())
    key.update(repr(plate_margin).encode())
    key = key.hexdigest()
    cache = f"{prefix}_skin_access_cache_pass{k}.npz"
    cached = np.load(cache) if os.path.isfile(cache) else None
    if cached is not None and str(cached["key"]) == key:
        plate, last[todo] = cached["plate"], cached["last"]
        timings_note = "cache"
    else:
        plate, last[todo] = last_accessible_moment_parallel(head, centers, angles, axes, cells, p0, pk, isovalues,
                                                            todo, plate_margin)
        np.savez(cache, key=key, plate=plate, last=last[todo])
        timings_note = "calcule"
    status[todo[plate]] = 2
    status[(status == 0) & (last < moment)] = 3
    # Rognage decide par P5 (decision utilisateur du 2026-09-28) : les faces a
    # moins de trim_radius d'un point ou la VRAIE tete touche (verification le
    # long des trajectoires, skin_print.py) passent en collision objet. Le
    # fichier cumule les points de toutes les iterations de skin_loop.py.
    if trim_file:
        import json
        trim_pts = np.asarray(json.load(open(trim_file))["points"], dtype=float).reshape(-1, 3)
        if len(trim_pts):
            near = cKDTree(trim_pts).query(tri.reshape(-1, 3))[0].reshape(-1, 3).min(axis=1)
            near = np.minimum(near, cKDTree(trim_pts).query(centers)[0])
            trimmed = (status == 0) & (near <= trim_radius)
            status[trimmed] = 3
            print(f"rognage P5 : {len(trim_pts)} points de collision, {trimmed.sum()} faces retirees "
                  f"(rayon {trim_radius} mm, {areas[trimmed].sum():.1f} mm2)", flush=True)
    timings[f"collision ({len(todo)} faces, {len(head)} points de tete, {timings_note})"] = time.time() - t0

    # --- morceaux : aretes vives, puis bandes minimales a frontieres lisses ---
    t0 = time.time()
    kept = np.flatnonzero(status == 0)
    adj = face_adjacency(faces[kept]).tocoo()
    cos_n = np.einsum("ij,ij->i", normals[kept][adj.row], normals[kept][adj.col])
    deviation = np.degrees(np.arccos(np.clip(cos_n, -1, 1)))      # angle entre normales
    cut_above = 180.0 - edge_face_angle
    if facet_smoothing is not None:
        cut_above = max(cut_above, facet_smoothing)
    link = deviation <= cut_above
    graph = coo_matrix((np.ones(link.sum()), (adj.row[link], adj.col[link])), shape=(len(kept),) * 2)
    _n, piece_of_kept = connected_components(graph, directed=False)
    piece = np.full(n_faces, -1)
    piece[kept] = piece_of_kept
    n_edge_pieces = piece_of_kept.max() + 1 if len(kept) else 0

    band_moments = {}                       # morceau -> moments choisis (croissants)
    levels_needed = {}                      # isovaleur -> faces (d'origine) a couper
    for p in range(n_edge_pieces):
        ids = np.flatnonzero(piece == p)
        points = minimal_stabbing(moment[ids], last[ids])
        band_moments[p] = points
        for m in points[:-1]:
            levels_needed.setdefault(float(isovalues[m]), []).append(ids)

    # coupe exacte le long des lignes de niveau (frontieres lisses)
    v, f, phi = surf_vertices, faces, surf_phi.copy()
    origin = np.arange(n_faces)             # face d'origine de chaque face courante
    for level in sorted(levels_needed):
        mask_orig = np.zeros(n_faces, bool)
        mask_orig[np.concatenate(levels_needed[level])] = True
        v, f, phi, par = split_along_level(v, f, phi, mask_orig[origin], level)
        origin = origin[par]
    timings[f"coupe le long de {len(levels_needed)} lignes de niveau"] = time.time() - t0

    # attributs herites de la face d'origine ; moment porteur recalcule
    t0 = time.time()
    status2, piece2, last2 = status[origin], piece[origin], last[origin]
    moment2 = layer_rank(phi[f].max(axis=1), isovalues)
    band = np.full(len(f), -1)
    for p in range(n_edge_pieces):
        sel = np.flatnonzero(piece2 == p)
        pts = band_moments[p]
        # plus petit moment choisi >= moment porteur ET <= dernier moment accessible
        j = np.searchsorted(pts, moment2[sel], side="left")
        ok = j < len(pts)
        ok[ok] = pts[j[ok]] <= last2[sel][ok]
        # repli (rare) : premier moment choisi valide pour la face
        for q in np.flatnonzero(~ok):
            valid = np.flatnonzero((pts >= moment2[sel[q]]) & (pts <= last2[sel[q]]))
            j[q] = valid[0] if len(valid) else -1
        band[sel] = np.where(j >= 0, j, -1)
    # zones = composantes connexes de (morceau, bande)
    tri2 = v[f]
    areas2 = 0.5 * np.linalg.norm(np.cross(tri2[:, 1] - tri2[:, 0], tri2[:, 2] - tri2[:, 0]), axis=1)
    zone = np.full(len(f), -1)
    zone_moment = []
    kept2 = np.flatnonzero((status2 == 0) & (band >= 0))
    adj2 = face_adjacency(f[kept2]).tocoo()
    same = (piece2[kept2][adj2.row] == piece2[kept2][adj2.col]) & (band[kept2][adj2.row] == band[kept2][adj2.col])
    g2 = coo_matrix((np.ones(same.sum()), (adj2.row[same], adj2.col[same])), shape=(len(kept2),) * 2)
    _n, comp = connected_components(g2, directed=False)
    for c in range(comp.max() + 1 if len(kept2) else 0):
        ids = kept2[comp == c]
        if areas2[ids].sum() < MIN_ZONE_AREA:
            status2[ids] = 5
            continue
        zone[ids] = len(zone_moment)
        zone_moment.append(int(moment2[ids].max()))          # au plus tot
    zone_moment = np.array(zone_moment, dtype=int)
    # micro-zones : rattachees a une voisine si la fusion reste accessible,
    # sinon NON lissees (decision utilisateur du 2026-09-28) -- statut 6, hors zone,
    # donc ni creusees par P3 ni tracees par P4
    t_micro = time.time()
    zone, merged_moment, zone_merges, micro_left = sp.merge_micro_zones(
        v, f, areas2, zone, zone_moment, last2, micro_inradius)
    for z in micro_left:
        status2[zone == z] = 6
        zone[zone == z] = -1
    live = sorted(set(zone[zone >= 0].tolist()))
    remap = np.full(len(zone_moment), -1)
    remap[live] = np.arange(len(live))
    zone = np.where(zone >= 0, remap[np.maximum(zone, 0)], -1)
    zone_moment = np.array([merged_moment[z] for z in live], dtype=int)
    zone_merges = np.array([(remap[q] if remap[q] >= 0 else -1, zz, m) for zz, q, _a, _b, m, _l in zone_merges],
                           dtype=int).reshape(-1, 3)
    print(f"micro-zones (distance au bord < {micro_inradius:.2f} mm) : {len(zone_merges)} rattachees, "
          f"{sum(1 for r in micro_left.values() if r == 'isolee')} isolees et "
          f"{sum(1 for r in micro_left.values() if r != 'isolee')} a voisines inaccessibles non lissees "
          f"({areas2[status2 == 6].sum():.1f} mm2) | {time.time() - t_micro:.1f} s", flush=True)
    ok_moment = all((last2[zone == z] >= zone_moment[z]).all() for z in range(len(zone_moment)))
    timings["bandes et zones"] = time.time() - t0

    # --- sorties -------------------------------------------------------------
    normals2, angles2 = normals[origin], angles[origin]
    np.savez(f"{prefix}_skin_pass{k}.npz", vertices=v, faces=f, status=status2, moment=moment2, zone=zone,
             zone_moment=zone_moment, angles=angles2, normals=normals2, areas=areas2, isovalues=isovalues,
             nozzle_width=NOZZLE_WIDTH, last_access=last2, origin=origin, edge_piece=piece2,
             surface_vertices=surf_vertices, surface_faces=faces, surface_status=status,
             edge_face_angle=edge_face_angle, zone_merges=zone_merges, micro_inradius=micro_inradius)
    colors = np.array([STATUS_COLORS[s] + (255,) for s in status2], dtype=np.uint8)
    trimesh.Trimesh(v, f, face_colors=colors, process=False).export(f"{prefix}_skin_map_pass{k}.ply")
    rng = np.random.default_rng(1)
    palette = rng.integers(60, 230, size=(max(len(zone_moment), 1), 3))
    zc = np.where(zone[:, None] >= 0, palette[np.maximum(zone, 0)], 225)
    zc = np.column_stack([zc, np.full(len(f), 255)]).astype(np.uint8)
    trimesh.Trimesh(v, f, face_colors=zc, process=False).export(f"{prefix}_skin_zones_pass{k}.ply")

    total = areas2.sum()
    print(f"{model_name}{variant} passe {k} : {n_faces} faces ({len(f)} apres coupe), {len(isovalues)} couches")
    for code, label in STATUS.items():
        m = status2 == code
        print(f"  {label:<18} {m.sum():>6} faces  {areas2[m].sum() / total:>6.1%} de la surface")
    n_bands = sum(len(b) for b in band_moments.values())
    print(f"  morceaux apres coupe aux aretes (< {edge_face_angle:.0f} deg entre faces"
          f"{'' if facet_smoothing is None else f', lissage {facet_smoothing:.0f} deg'}) : {n_edge_pieces} ; "
          f"bandes (minimum) : {n_bands} ; zones : {len(zone_moment)}"
          f" (moments {zone_moment.min() if len(zone_moment) else '-'} a {zone_moment.max() if len(zone_moment) else '-'})")
    print(f"  chaque zone accessible a son moment : {ok_moment}")
    print(f"  sommets de surface non localises dans un tetraedre : {int(lost.sum())}")
    print("  temps : " + " ; ".join(f"{name} {dt:.1f} s" for name, dt in timings.items())
          + f" ; total {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
