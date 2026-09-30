"""
Lissage non-planaire -- P2 : quelles faces de la surface d'ORIGINE recoivent
une couche de surface, a quel moment, avec quelle orientation de tete.

Pour chaque face exterieure (hors contact plateau) :
  - MOMENT au plus tot : rang de la couche qui porte le sommet le plus haut de
    la face (hauteur deformee phi) -- des cette couche finie, toute la matiere
    sous la face est deposee ;
  - ORIENTATION : la buse pointe vers -n (n = normale sortante). Angles d'axes
    par l'IK de l'etape 8 d'Apex (_resolve_axis_angles, repli par grille
    compris), sur la cinematique reelle (main._default_machine_profile) ;
  - PRE-CALCUL de collision, avant l'etape 9 (qui ne verra que les survivantes) :
      plateau : un point de la tete sous le plan du plateau (z < 0) ;
      objet   : un point de la tete DANS la matiere deja imprimee a ce moment,
                c'est-a-dire dans la piece ET a phi <= phi du moment. C'est
                ce qui accepte une surface fermee plus tard (interieur du
                toit du Benchy) : seule la matiere deja la compte.
    Les points de la tete a moins de TIP_EXCLUSION de la pointe sont ignores :
    la buse touche la surface par definition (meme exclusion que le cordon en
    cours a l'etape 9).
Zones = faces retenues voisines. Une zone s'imprime au moment ou sa derniere
couche porteuse est finie, si toutes ses faces y sont encore accessibles ;
sinon elle est coupee en deux bandes de moment, recursivement. Une zone plus
petite que MIN_ZONE_AREA est abandonnee (et ne sera pas creusee).

Aucune constante propre a un modele : les seuils sont relatifs a la largeur
de buse ; le plateau est le plan z = 0 du repere Bird, qui est le repere
machine d'Apex (platform_origin = 0, axe C = Z).

Repere : Bird (loop_*_pass*.npy), PAS celui du STL.

Usage (venv d'Apex) :
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe surface_skin.py <modele> <passe> [variante]
Entrees : output_models/loop_<modele><variante>_{cells,pass0,pass<k>}.npy et
          _deformed_slice_pass<k>.pkl (isovaleurs des couches, bird_slice_deformed.py).
Sortie  : output_models/loop_<modele><variante>_skin_pass<k>.npz
          + _skin_map_pass<k>.ply (faces colorees par statut, pour validation).
"""
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
from pipeline.step8_machine_toolpath import Step8Params, _resolve_axis_angles  # noqa: E402
from pipeline.step9_collision_sim import _positioned_toolhead_vertices, _surface_sample_points  # noqa: E402

NOZZLE_WIDTH = 0.4
TIP_EXCLUSION = 2.0 * NOZZLE_WIDTH        # la buse touche la surface
HEAD_SAMPLE_SPACING = 5.0                 # comme l'etape 9 (limite #26)
MIN_ZONE_AREA = 4.0 * NOZZLE_WIDTH ** 2   # plus petit : pas de couche de surface
# Arete maximale de la surface analysee. Un STL grossier (l_shape : 44 faces
# de plusieurs cm) donnait UN point de test et UN moment par grande face ;
# la surface est subdivisee jusqu'a cette taille (les faces deja petites ne
# bougent pas, la geometrie est inchangee).
MAX_EDGE = 5.0 * NOZZLE_WIDTH
PLATE_CONTACT_TOL = 1e-3                  # faces posees sur le plateau (cf. bird_core)

STATUS = {0: "retenue", 1: "inatteignable", 2: "collision plateau", 3: "collision objet",
          4: "contact plateau", 5: "zone trop petite"}
STATUS_COLORS = {0: (40, 170, 70), 1: (150, 0, 150), 2: (230, 140, 0), 3: (220, 30, 30),
                 4: (120, 120, 120), 5: (90, 140, 230)}


class PrintedSolid:
    """La matiere deja imprimee au moment t : les points de la piece
    d'ORIGINE dont la hauteur deformee phi est <= phi(t). Localisation d'un
    point dans les tetraedres d'origine (k plus proches par centre)."""

    def __init__(self, cells, original, deformed, k=16):
        self.tets = original[cells]
        self.phi = deformed[cells][:, :, 2]
        self.tree = cKDTree(self.tets.mean(axis=1))
        self.k = k
        self.reach = float(np.linalg.norm(self.tets - self.tets.mean(axis=1, keepdims=True), axis=2).max())
        self.tolerance = -1e-9

    def phi_at(self, points):
        """phi aux points (NaN hors de la piece)."""
        dist, cand = self.tree.query(points, k=self.k)
        out = np.full(len(points), np.nan)
        tolerance = self.tolerance
        pending = dist[:, 0] <= self.reach
        for j in range(self.k):
            idx = np.flatnonzero(pending)
            if not len(idx):
                break
            t = cand[idx, j]
            with np.errstate(all="ignore"):
                bary = barycentric(points[idx], self.tets[t])
            inside = bary.min(axis=1) >= tolerance
            hit = idx[inside]
            out[hit] = np.einsum("nv,nv->n", bary[inside], self.phi[t[inside]])
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


def main():
    model_name = sys.argv[1] if len(sys.argv) > 1 else "benchy upsidedown tilted"
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    variant = sys.argv[3] if len(sys.argv) > 3 else ""
    prefix = os.path.join(HERE, "output_models", "loop_" + model_name.replace(" ", "_") + variant)
    t_start = time.time()

    cells = np.load(f"{prefix}_cells.npy")
    p0 = np.load(f"{prefix}_pass0.npy")
    pk = np.load(f"{prefix}_pass{k}.npy")
    layers = pickle.load(open(f"{prefix}_deformed_slice_pass{k}.pkl", "rb"))
    isovalues = np.sort(np.unique([layer["isovalue"] for layer in layers]))

    # surface d'origine = faces frontiere des tetraedres, subdivisee (MAX_EDGE)
    surf_vertices, faces = trimesh.remesh.subdivide_to_size(p0, boundary_faces(cells, p0), MAX_EDGE)
    solid = PrintedSolid(cells, p0, pk)
    solid.tolerance = -1e-6                 # sommets de surface : sur la frontiere
    surf_phi = solid.phi_at(surf_vertices)
    solid.tolerance = -1e-9
    tri = surf_vertices[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    areas = 0.5 * np.linalg.norm(normals, axis=1)
    normals /= np.maximum(2 * areas, 1e-12)[:, None]
    centers = tri.mean(axis=1)
    n_faces = len(faces)
    status = np.zeros(n_faces, dtype=np.int8)

    on_plate = tri[:, :, 2].max(axis=1) < p0[:, 2].min() + PLATE_CONTACT_TOL
    status[on_plate] = 4
    lost = np.isnan(surf_phi)
    if lost.any():                          # sommet non localise : plus proche sommet d'origine
        surf_phi[lost] = pk[cKDTree(p0).query(surf_vertices[lost])[1], 2]
    moment = layer_rank(surf_phi[faces].max(axis=1), isovalues)

    # --- orientation : IK de l'etape 8, buse vers -n --------------------------
    machine = apex_main._default_machine_profile()
    axes = machine.rotation_axes
    reference = np.array([0.0, 0.0, -1.0])
    step8 = Step8Params()
    angles = np.zeros((n_faces, len(axes)))
    residual = np.zeros(n_faces)
    previous = None
    for i in np.flatnonzero(status == 0):
        a, err = _resolve_axis_angles(-normals[i], axes, reference, step8, previous)
        angles[i], residual[i], previous = a, err, a
    unreachable = (status == 0) & (np.degrees(residual) > step8.max_unreachable_angular_error_deg)
    status[unreachable] = 1

    # --- pre-calcul de collision ---------------------------------------------
    head = _surface_sample_points(machine.nozzle_model, HEAD_SAMPLE_SPACING)
    head = head[np.linalg.norm(head, axis=1) > TIP_EXCLUSION]

    def collides(face_ids, at_moment):
        """Code de statut (0, 2 ou 3) de chaque face, tete placee au centre de
        la face, contre la matiere imprimee au moment `at_moment` (par face)."""
        out = np.zeros(len(face_ids), dtype=np.int8)
        for j, f in enumerate(face_ids):
            pts = _positioned_toolhead_vertices(head, centers[f], list(angles[f]), axes)
            if (pts[:, 2] < 0.0).any():
                out[j] = 2
                continue
            phi = solid.phi_at(pts)
            if np.any(phi <= isovalues[at_moment[j]] + 1e-9):
                out[j] = 3
        return out

    candidates = np.flatnonzero(status == 0)
    status[candidates] = collides(candidates, moment[candidates])
    t_pre = time.time() - t_start

    # --- zones et moment d'impression ----------------------------------------
    adjacency = face_adjacency(faces)
    zone_of = np.full(n_faces, -1)
    zone_moment = []

    def assign(face_ids):
        """Place un ensemble connexe de faces retenues : au moment de sa
        derniere couche porteuse si toutes y sont accessibles, sinon coupe en
        deux bandes de moment."""
        if areas[face_ids].sum() < MIN_ZONE_AREA:
            status[face_ids] = 5
            return
        when = int(moment[face_ids].max())
        blocked = collides(face_ids, np.full(len(face_ids), when)) != 0
        if not blocked.any() or moment[face_ids].min() == when:
            if blocked.any():                      # bande reduite a un seul moment
                status[face_ids[blocked]] = 3
                face_ids = face_ids[~blocked]
                if not len(face_ids) or areas[face_ids].sum() < MIN_ZONE_AREA:
                    status[face_ids] = 5
                    return
            zone_of[face_ids] = len(zone_moment)
            zone_moment.append(when)
            return
        # milieu entre le premier et le dernier moment (qui different ici) :
        # les deux bandes sont non vides -- la mediane ne le garantit pas
        # quand la plupart des faces partagent le moment maximal
        cut = (moment[face_ids].min() + moment[face_ids].max()) // 2
        for band in (face_ids[moment[face_ids] <= cut], face_ids[moment[face_ids] > cut]):
            split(band)

    def split(face_ids):
        if not len(face_ids):
            return
        sub = adjacency[face_ids][:, face_ids]
        n_comp, labels = connected_components(sub, directed=False)
        for c in range(n_comp):
            assign(face_ids[labels == c])

    split(np.flatnonzero(status == 0))
    zone_moment = np.array(zone_moment, dtype=int)

    # --- sorties -------------------------------------------------------------
    np.savez(f"{prefix}_skin_pass{k}.npz", vertices=surf_vertices, faces=faces, status=status, moment=moment, zone=zone_of,
             zone_moment=zone_moment, angles=angles, normals=normals, areas=areas,
             isovalues=isovalues, nozzle_width=NOZZLE_WIDTH)
    colors = np.array([STATUS_COLORS[s] + (255,) for s in status], dtype=np.uint8)
    trimesh.Trimesh(surf_vertices, faces, face_colors=colors, process=False).export(f"{prefix}_skin_map_pass{k}.ply")

    total = areas.sum()
    print(f"{model_name}{variant} passe {k} : {n_faces} faces, {len(isovalues)} couches")
    for code, label in STATUS.items():
        m = status == code
        print(f"  {label:<18} {m.sum():>6} faces  {areas[m].sum() / total:>6.1%} de la surface")
    print(f"  sommets de surface non localises dans un tetraedre : {int(lost.sum())}")
    print(f"  zones : {len(zone_moment)} (moments {zone_moment.min() if len(zone_moment) else '-'}"
          f" a {zone_moment.max() if len(zone_moment) else '-'})")
    print(f"  temps : pre-calcul {t_pre:.0f} s, total {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
