# Porte de pipeline/bird/skin/surface_paths.py (branche bird_old, 2026-10-28) pour la voie chaleur : code de calcul copie tel quel (il ne depend d'aucun champ).
# Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer).
# Porte de research/S4_Slicer/surface_paths.py (2026-10-07, plan bird-into-apex-cpp) : code de calcul
# copie TEL QUEL, seuls les imports changent (egalite au bit avec la recherche).
# Fait partie de la voie Bird d'Apex : derive de S4 Slicer (Joshua Bird, GPL-3.0).
"""
Lissage non-planaire -- P4 : trajectoires des couches de lissage, tracees
DIRECTEMENT sur la surface (decision de l'utilisateur, 2026-09-27).

Principe : 9 motifs sur 10 sont des LIGNES DE NIVEAU d'une distance mesuree
sur la surface depuis une graine (bord, centre, forme englobante, ligne de
coupe), espacees de la largeur de buse. Une distance geodesique a un gradient
de norme 1 sur la surface : l'espacement est exact par construction, meme sur
une forme bombee (seule la forme du motif se deforme, jamais l'espacement).

Etape 4a (ce fichier, pour l'instant) : outils de surface.
  - ZoneMesh : sous-maillage d'une zone de P2-bis, raffine (aretes <= REFINE)
    -- mesure : a 0,38 mm d'arete, espacement 0,401-0,425 mm pour 0,4 vise ;
    a 0,19 mm, 0,401-0,408 (demi-sphere R = 10) ;
  - distances geodesiques par la methode de la CHALEUR SIGNEE (potpourri3d,
    MeshSignedHeatSolver) depuis des courbes (bord, graines) ou des sommets.
    Mesure de la POSITION des lignes (ecart a l'exact median / 95e / max) :
      demi-sphere : signee 0,002/0,016/0,031 ; fast marching 0,007/0,035/0,066 ;
                    chaleur simple 0,045/0,148/0,213 (lignes toutes GLISSEES :
                    espacement juste, premier anneau et recouvrement faux) ;
      disque en triangles-aiguilles : signee 0,059/0,228/0,357 ; chaleur
                    0,152/0,435/0,720 ; fast marching 0,348/1,100/1,278 (le fast
                    marching se degrade sur les triangles obtus) ;
    raccourcir la diffusion de la chaleur simple EMPIRE (t = 0,1 : x6) ;
    courbe OUVERTE traversant la zone de bord a bord (graine du rectiligne) :
    distance SIGNEE puis valeur absolue, 0,001 mm (demi-cylindre) -- en mode
    non signe, la chaleur signee donne n'importe quoi (3,5 mm d'ecart median) ;
    depuis un POINT, la chaleur signee est mediocre (0,13-0,17 mm) : la
    distance au centre est le RAYON des coordonnees polaires (0,000-0,011) ;
  - coordonnees polaires geodesiques (log map, methode de la chaleur
    vectorielle) ;
  - lignes de niveau (marching triangles) rendues en 3D avec leur face et
    leurs coordonnees barycentriques (pour la normale, en P4 plus loin).

Dependance : potpourri3d (MIT ; numpy, scipy), venv d'Apex.
"""
import os
import sys

import numpy as np
import potpourri3d as pp3d
import trimesh

from element_ardor.skin.surface_mesh import (  # noqa: F401  (réexport : imports existants)
    GOOD_MESH_QUALITY,
    NOZZLE_WIDTH,
    PatternResult,
    Polyline,
    REFINE,
    SEED_JITTER,
    ZoneMesh,
    boundary_edges,
    boundary_loops,
    split_bowtie_vertices,
)
from element_ardor.skin.surface_patterns import (  # noqa: F401  (réexport : imports existants)
    _connector,
    _nearest_value,
    _nested,
    aligned_rectilinear,
    archimedean,
    clip_keep_above,
    concentric_circles,
    concentric_conform,
    concentric_shape,
    concentric_square,
    concentric_triangle,
    densify,
    distance_from_center,
    first_level,
    gaps_between,
    line_direction,
    merge_polylines,
    monotonic_line,
    monotonic_strips,
    octagram_spiral,
    order_by_level,
    parallel_lines,
    path_neighbor_gaps,
    rectilinear,
    shape_field,
    signed_line_distance,
    spiral,
    zigzag_strips,
    zone_center,
)
from element_ardor.skin.surface_hilbert import (  # noqa: F401  (réexport : imports existants)
    HILBERT_STRETCH_TOLERANCE,
    _hilbert_d2xy,
    arap_stretch,
    arap_unfold,
    geodesic_grid_coords,
    hilbert,
    hilbert_geodesic,
    hilbert_stretch,
    locate_uv,
)


# --- cordon de bord (decision utilisateur 2026-09-28) --------------------------
EDGE_OVERLAP = 0.25        # recouvrement motif / cordon de bord, en largeurs (usage des slicers plans ;
                           # mesure : 0,15 -> -0,5 a -2,6 pts de couverture ; a fixer apres essai)


def with_edge_ring(zm, pattern, width=NOZZLE_WIDTH, overlap=0.2, edge_rings=1, edge_overlap=EDGE_OVERLAP,
                   d_edge=None, **kw):
    """Motif precede de `edge_rings` cordons qui suivent le bord de la zone :
    lignes de niveau de la distance au bord MESUREE SUR LA SURFACE (le cordon
    epouse la courbe de la piece), le premier a w/2 - recouvrement du bord,
    les suivants tous les w. Le motif commence un cordon plus loin, moins
    `edge_overlap` largeur (le conforme continue exactement au pas w, ses
    anneaux etant paralleles au cordon). Le motif recoit un recouvrement
    EFFECTIF tel que son premier niveau tombe la : ses raccords (zigzag)
    longent sa propre limite, jamais le cordon deja depose."""
    d_edge = zm.distance_to_boundary() if d_edge is None else d_edge
    fn = pattern if callable(pattern) else PATTERNS[pattern]
    if edge_rings <= 0:
        return fn(zm, width, overlap, d_edge=d_edge, **kw)
    base = first_level(width, overlap)
    rings = [(q, pl.points) for q, pl in zm.lines_multi(d_edge, base + width * np.arange(edge_rings))]
    rings.sort(key=lambda kp: kp[0])
    inner = base + width * (edge_rings - (0.0 if pattern == "concentric_conform" else edge_overlap))
    res = fn(zm, width, width / 2.0 - inner, d_edge=d_edge, **kw)
    ring_paths = order_by_level(rings)
    # dernier cordon ferme : il commence et finit au point le plus proche du
    # debut du motif (sinon un trajet a vide traversait la zone)
    if ring_paths and res.paths:
        last = ring_paths[-1]
        if len(last) > 2 and np.linalg.norm(last[0] - last[-1]) < 1e-6:
            i = int(np.argmin(np.linalg.norm(last - res.paths[0][0], axis=1)))
            ring_paths[-1] = np.vstack([last[i:], last[1:i + 1]])
    return PatternResult(ring_paths + list(res.paths), res.level_lines, d_edge)


HILBERT_WARN_FRACTION = 0.10         # avertissement si plus de 10 % de la surface depasse

HILBERT_OPTIONS = {
    "A": ("Garder tel quel",
          "Le motif reste parfaitement regulier, mais sur les parties bombees les lignes "
          "s'ecartent ou se serrent un peu : de petits jours ou de legeres surepaisseurs "
          "peuvent apparaitre a ces endroits."),
    "B": ("Decouper la zone",
          "La zone est partagee en morceaux assez plats pour que les lignes restent bien "
          "espacees. Le motif reste regulier dans chaque morceau, mais on verra des raccords "
          "entre les morceaux."),
    "C": ("Resserrer la grille",
          "L'ecart entre les lignes redevient juste en moyenne, mais les carres du motif se "
          "deforment en losanges. Sur une forme tres bombee, les pointes des losanges se "
          "touchent presque : les lignes s'y serrent et de petits manques peuvent apparaitre."),
}
# mesure (test du 2026-09-28, part des lignes a w +- 5 %, reference plane 94 %) :
#   demi-sphere A 27 %  B 90 %  C 74 % ; coque Benchy A 88 %  B 93 %  C 88 % ;
#   dino A 52 %  B 92 %  C 74 %. B coute 7 a 18 s par grande zone (ARAP par morceau).


def split_zone(zm, min_area=25 * NOZZLE_WIDTH ** 2):
    """Remede B : coupe recursive en deux (plan par le barycentre, normal a
    l'axe principal de la zone) jusqu'a ce que chaque morceau se deplie avec
    moins de HILBERT_WARN_FRACTION de surface etiree, ou devienne trop petit.
    Rend [(ZoneMesh du morceau, uv ARAP, part etiree)]. Coupe au niveau des
    triangles du maillage de calcul (0,2 mm) : limite en escalier fin, lisse
    ensuite par le cordon de bord."""
    uv, frac = hilbert_stretch(zm)
    tri = zm.V[zm.F]
    area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    if frac <= HILBERT_WARN_FRACTION or area.sum() < 2 * min_area:
        return [(zm, uv, frac)]
    c = tri.mean(axis=1)
    g = (c * area[:, None]).sum(0) / area.sum()
    _u, _s, vt = np.linalg.svd((c - g) * np.sqrt(area)[:, None], full_matrices=False)
    side = (c - g) @ vt[0] >= 0
    out = []
    for keep in (side, ~side):
        piece = ZoneMesh(zm.V, zm.F[keep], refine=np.inf)
        piece.parent = zm.parent[keep][piece.parent]
        out.extend(split_zone(piece, min_area))
    return out


def hilbert_split(zm, width=NOZZLE_WIDTH, overlap=0.2, angle_deg=0.0, d_edge=None, piece_rings=0,
                  piece_overlap=None, **_):
    """Hilbert, remede B : un Hilbert par morceau, chacun avec sa propre
    distance au bord et, s'il est demande, son propre cordon de bord : le
    long d'une coupe, les deux cordons voisins sont a w l'un de l'autre."""
    piece_overlap = EDGE_OVERLAP if piece_overlap is None else piece_overlap
    paths = []
    for piece, uv, _frac in split_zone(zm):
        paths.extend(with_edge_ring(piece, hilbert, width, overlap, piece_rings, piece_overlap,
                                    angle_deg=angle_deg, uv=uv).paths)
    return PatternResult(paths, [], zm.distance_to_boundary() if d_edge is None else d_edge)


# --- micro-zones (decision utilisateur 2026-09-28) -----------------------------
def merge_micro_zones(vertices, faces, areas, zone, zone_moment, last_access, micro_inradius):
    """Rattache chaque micro-zone (plus grande distance au bord, mesuree sur
    la surface, < micro_inradius : trop etroite pour un cordon de bord et un
    motif) a une zone VOISINE (arete commune) ou la fusion reste imprimable :
    la zone fusionnee s'imprime au plus tard des deux moments, qui doit rester
    <= au dernier moment accessible de TOUTES ses faces. Parmi les voisines
    admissibles : la plus longue frontiere commune. Repete jusqu'a stabilite.
    Rend (zone, moments {z: m}, journal [(z, cible, m_z, m_cible, m, frontiere)],
    restantes {z: "isolee" | "voisines inaccessibles"}) ; les restantes ne
    sont PAS lissees (decision utilisateur) : c'est a l'appelant de les retirer."""
    zone = zone.copy()
    moment = {z: int(m) for z, m in enumerate(zone_moment)}
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    fid = np.tile(np.arange(len(faces)), 3)
    key = np.sort(e, axis=1)
    o = np.lexsort((key[:, 1], key[:, 0]))
    key, fid = key[o], fid[o]
    same = np.all(key[1:] == key[:-1], axis=1)
    pa, pb = fid[:-1][same], fid[1:][same]                     # paires de faces voisines
    elen = np.linalg.norm(vertices[key[:-1][same, 0]] - vertices[key[:-1][same, 1]], axis=1)
    radius = {}

    def inradius(z):
        if z not in radius:
            try:
                radius[z] = float(ZoneMesh(vertices, faces[zone == z]).distance_to_boundary().max())
            except Exception:                                  # zone degeneree : traitee comme micro
                radius[z] = 0.0
        return radius[z]

    log = []
    changed = True
    while changed:
        changed = False
        live = [z for z in moment if (zone == z).any()]
        for z in sorted(live, key=lambda q: areas[zone == q].sum()):
            if not (zone == z).any() or inradius(z) >= micro_inradius:
                continue
            za, zb = zone[pa], zone[pb]
            touch = (za == z) ^ (zb == z)
            border = {}
            for q, w in zip(np.where(za == z, zb, za)[touch], elen[touch]):
                if q >= 0:
                    border[int(q)] = border.get(int(q), 0.0) + float(w)
            best = None
            for q, length in border.items():
                m = max(moment[z], moment[q])
                if m <= last_access[(zone == z) | (zone == q)].min() and (best is None or length > best[1]):
                    best = (q, length, m)
            if best is None:
                continue
            q, length, m = best
            zone[zone == z] = q
            log.append((z, q, moment[z], moment[q], m, length))
            moment[q] = m
            radius.pop(q, None)                                # la cible a grandi : rayon a recalculer
            changed = True
    remaining = {}
    for z in [z for z in moment if (zone == z).any()]:
        if inradius(z) < micro_inradius:
            za, zb = zone[pa], zone[pb]
            touch = (za == z) ^ (zb == z)
            has_neighbour = (np.where(za == z, zb, za)[touch] >= 0).any()
            remaining[z] = "voisines inaccessibles" if has_neighbour else "isolee"
    return zone, moment, log, remaining


# --- registre des motifs (memes noms que surface_patterns.PATTERNS) ----------
PATTERNS = {
    "concentric_conform": concentric_conform,
    "concentric": concentric_circles,
    "rectilinear": rectilinear,
    "aligned_rectilinear": aligned_rectilinear,
    "monotonic_line": monotonic_line,
    "concentric_square": concentric_square,
    "concentric_triangle": concentric_triangle,
    "archimedean": archimedean,
    "octagram_spiral": octagram_spiral,
    "hilbert": hilbert,
}
