"""
P4, mesure prealable : de combien le DEPLIAGE a plat deforme-t-il les zones de
lissage ? (decide la compensation des motifs deplies : Hilbert, spirales,
concentriques carre / triangulaire / octagramme.)

Pour chaque zone de surface_skin.py (loop_*_skin_pass<k>.npz), d'abord coupee
aux ARETES VIVES (angle entre faces < seuil, 150 deg par defaut = 30 deg
entre normales, decision de l'utilisateur), chaque morceau est deplie par
LSCM (applications conformes aux moindres carres, scipy.sparse), puis remis a
l'echelle pour conserver son aire. Par triangle : valeurs singulieres s1 >= s2
du passage 3D -> plan. Un motif trace a l'espacement e dans le plan aura sur
la piece un espacement e / s (s = echelle locale) : c'est cet ecart qu'on
mesure, pondere par l'aire.

Usage (venv d'Apex) :
  ../../.venv/Scripts/python.exe measure_unfold.py <modele> <passe> <variante> [angle_entre_faces_deg]
"""
import os
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import lsqr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from surface_skin import face_adjacency  # noqa: E402


def split_at_sharp_edges(faces, normals, face_ids, min_face_angle_deg):
    """Composantes connexes de la zone, sans traverser une arete vive (angle
    ENTRE LES FACES < min_face_angle_deg, soit un angle entre normales >
    180 - min_face_angle_deg)."""
    sub = face_adjacency(faces[face_ids]).tocoo()
    cos_limit = np.cos(np.radians(180.0 - min_face_angle_deg))
    keep = np.einsum("ij,ij->i", normals[face_ids][sub.row], normals[face_ids][sub.col]) >= cos_limit
    graph = coo_matrix((np.ones(keep.sum()), (sub.row[keep], sub.col[keep])), shape=(len(face_ids),) * 2)
    _n, labels = connected_components(graph, directed=False)
    return [face_ids[labels == c] for c in range(labels.max() + 1)]


def lscm(vertices, faces):
    """Depliage conforme (Levy et al. 2002). Rend uv (V,2) sur les sommets
    utilises (index locaux) et les sommets locaux."""
    used, local = np.unique(faces, return_inverse=True)
    f = local.reshape(-1, 3)
    p = vertices[used]
    a, b, c = p[f[:, 0]], p[f[:, 1]], p[f[:, 2]]
    e1 = b - a
    n = np.cross(e1, c - a)
    area2 = np.linalg.norm(n, axis=1)
    x_axis = e1 / np.linalg.norm(e1, axis=1)[:, None]
    y_axis = np.cross(n / area2[:, None], x_axis)
    # coordonnees locales 2D des trois sommets (a = origine)
    q = np.zeros((len(f), 3, 2))
    q[:, 1, 0] = np.linalg.norm(e1, axis=1)
    q[:, 2, 0] = np.einsum("ij,ij->i", c - a, x_axis)
    q[:, 2, 1] = np.einsum("ij,ij->i", c - a, y_axis)
    # W_j = (x_k - x_i) + i (y_k - y_i) pour (i, j, k) cyclique, / sqrt(aire)
    w = np.empty((len(f), 3), complex)
    for j, (i_, k_) in enumerate([(1, 2), (2, 0), (0, 1)]):
        d = q[:, k_] - q[:, i_]
        w[:, j] = (d[:, 0] + 1j * d[:, 1])
    w /= np.sqrt(area2 / 2)[:, None]
    nv = len(used)
    rows = np.repeat(np.arange(len(f)), 3)
    m = coo_matrix((w.ravel(), (rows, f.ravel())), shape=(len(f), nv)).tocsr()
    # deux sommets epingles : les plus eloignes (approximativement)
    i0 = int(np.argmax(np.linalg.norm(p - p.mean(0), axis=1)))
    i1 = int(np.argmax(np.linalg.norm(p - p[i0], axis=1)))
    pinned = np.array([i0, i1])
    pin_uv = np.array([[0.0, 0.0], [np.linalg.norm(p[i1] - p[i0]), 0.0]])
    free = np.setdiff1d(np.arange(nv), pinned)
    mr, mi = m.real.tocsc(), m.imag.tocsc()
    # u + i v : Re et Im de M (u + i v) = 0
    from scipy.sparse import bmat
    a_free = bmat([[mr[:, free], -mi[:, free]], [mi[:, free], mr[:, free]]]).tocsr()
    a_pin = bmat([[mr[:, pinned], -mi[:, pinned]], [mi[:, pinned], mr[:, pinned]]]).tocsr()
    rhs = -a_pin @ np.r_[pin_uv[:, 0], pin_uv[:, 1]]
    sol = lsqr(a_free, rhs, atol=1e-12, btol=1e-12, iter_lim=20000)[0]
    uv = np.zeros((nv, 2))
    uv[pinned] = pin_uv
    uv[free, 0], uv[free, 1] = sol[:len(free)], sol[len(free):]
    return uv, f, q, area2 / 2


def stretch(uv, f, q, area):
    """Valeurs singulieres (s1 >= s2) du passage triangle 3D -> plan, apres
    remise a l'echelle de l'aire totale ; et aire du triangle deplie signee."""
    t = uv[f]
    du = np.stack([t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]], axis=2)       # (T,2,2) colonnes = aretes planes
    dq = np.stack([q[:, 1] - q[:, 0], q[:, 2] - q[:, 0]], axis=2)       # aretes locales 3D
    jac = du @ np.linalg.inv(dq)
    s = np.linalg.svd(jac, compute_uv=False)
    area_uv = 0.5 * np.abs(np.linalg.det(du))
    scale = np.sqrt(area.sum() / max(area_uv.sum(), 1e-30))
    flipped = np.linalg.det(du) * np.sign(np.linalg.det(du).sum()) < 0
    return s * scale, flipped


def weighted_percentile(values, weights, pct):
    order = np.argsort(values)
    cum = np.cumsum(weights[order]) / weights.sum()
    return np.interp(np.asarray(pct) / 100.0, cum, values[order])


def main():
    model, k, variant = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    min_face_angle = float(sys.argv[4]) if len(sys.argv) > 4 else 150.0
    prefix = os.path.join(HERE, "output_models", "loop_" + model.replace(" ", "_") + variant)
    skin = np.load(f"{prefix}_skin_pass{k}.npz")
    vertices, faces, normals, zone = skin["vertices"], skin["faces"], skin["normals"], skin["zone"]
    t0 = time.time()
    pieces = []
    for z in range(zone.max() + 1):
        pieces += split_at_sharp_edges(faces, normals, np.flatnonzero(zone == z), min_face_angle)
    all_s, all_w, rows = [], [], []
    flips = 0
    for ids in pieces:
        if len(ids) < 2:
            continue
        uv, f, q, area = lscm(vertices, faces[ids])
        s, flipped = stretch(uv, f, q, area)
        flips += int(flipped.sum())
        # echelle locale (moyenne geometrique) et anisotropie (conformite)
        local = np.sqrt(s[:, 0] * s[:, 1])
        dev = np.abs(1.0 / local - 1.0)                  # ecart relatif d'espacement du motif
        all_s.append(dev)
        all_w.append(area)
        rows.append((area.sum(), weighted_percentile(dev, area, 95), (s[:, 0] / s[:, 1]).max()))
    dev, w = np.concatenate(all_s), np.concatenate(all_w)
    rows = np.array(rows)
    print(f"{model}{variant} passe {k} : {zone.max() + 1} zones -> {len(pieces)} morceaux apres coupe aux aretes "
          f"(< {min_face_angle:.0f} deg entre faces) ; depliage {time.time() - t0:.1f} s ; triangles retournes {flips}")
    p50, p90, p95, p99 = weighted_percentile(dev, w, [50, 90, 95, 99])
    print(f"ecart d'espacement d'un motif deplie (pondere par l'aire) : median {p50:.1%}, 90e {p90:.1%}, 95e {p95:.1%}, 99e {p99:.1%}")
    for lim in (0.05, 0.10, 0.20):
        print(f"  surface ou l'ecart depasse {lim:.0%} : {w[dev > lim].sum() / w.sum():.1%}")
    big = rows[:, 0] >= np.percentile(rows[:, 0], 90)
    print(f"morceaux : aire mediane {np.median(rows[:, 0]):.1f} mm2, max {rows[:, 0].max():.1f} ; "
          f"les 10 % plus grands : ecart 95e centile median {np.median(rows[big, 1]):.1%}, max {rows[big, 1].max():.1%}")
    print(f"morceaux dont le 95e centile depasse 10 % : {(rows[:, 1] > 0.10).sum()} sur {len(rows)}, "
          f"soit {rows[rows[:, 1] > 0.10, 0].sum() / rows[:, 0].sum():.1%} de la surface")


if __name__ == "__main__":
    main()
