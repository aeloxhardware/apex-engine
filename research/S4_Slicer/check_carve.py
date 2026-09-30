"""
Verification de P3 (corps creuse sous les zones de lissage non-planaire).

Chaque point des contours ramenes dans la piece (sortie _deformed_slice_skin)
est classe :
  surface  : sur la surface d'origine non lissee (< TOL_SURF) ;
  sous-peau: a l'epaisseur de lissage (+- 25 %) de la face RETENUE la plus
             proche -- la ou le corps doit s'arreter ;
  autre    : ni l'un ni l'autre -> anomalie, localisee par couche.

Usage (venv d'Apex) :
  ../../.venv/Scripts/python.exe check_carve.py <modele> <passe> <variante> [epaisseur]
"""
import os
import pickle
import sys

import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bird_to_apex import frame_offset  # noqa: E402

TOL_SURF = 0.05


def main():
    model, k, variant = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    thickness = float(sys.argv[4]) if len(sys.argv) > 4 else 0.2
    prefix = os.path.join(HERE, "output_models", "loop_" + model.replace(" ", "_") + variant)
    stl = os.path.join(HERE, "input_models", f"{model}.stl")
    offset = frame_offset(stl, prefix)
    layers = pickle.load(open(f"{prefix}_deformed_slice_skin_pass{k}.pkl", "rb"))
    skin = np.load(f"{prefix}_skin_pass{k}.npz")
    keep = (skin["status"] == 0) & (skin["zone"] >= 0)
    surface = trimesh.Trimesh(skin["vertices"], skin["faces"], process=False)
    retained = trimesh.Trimesh(skin["vertices"], skin["faces"][keep], process=False)

    # INVARIANT : creuser n'ajoute jamais de matiere. Region creusee (vue de
    # dessus, parite des contours) incluse dans la region de reference, a la
    # tolerance pres. Un trou bouche ne laisse aucun point de contour : seul ce
    # test le voit.
    from shapely.geometry import Polygon
    ref = {L["isovalue"]: L for L in pickle.load(open(f"{prefix}_deformed_slice_pass{k}.pkl", "rb"))}

    def region(loops):
        out = None
        for c in loops:
            c = np.asarray(c)
            if len(c) < 4:
                continue
            p = Polygon(c[:, :2]).buffer(0)
            out = p if out is None else out.symmetric_difference(p)
        return out

    added = []
    for i, L in enumerate(layers):
        R = ref.get(L["isovalue"])
        if R is None or not L["contour_3d"]:
            continue
        a, b = region(L["contour_3d"]), region(R["contour_3d"])
        if a is None or b is None:
            continue
        extra = a.difference(b.buffer(0.05)).area      # 0,05 mm : projection d'une couche courbe
        added.append((extra, i))
    added.sort(reverse=True)
    bad = [(e, i) for e, i in added if e > 0.05]
    print(f"matiere AJOUTEE par le creusement (> 0,05 mm2) : {len(bad)} couches sur {len(added)}"
          + ("" if not bad else " ; pires : " + ", ".join(f"{i + 1} ({e:.2f} mm2)" for e, i in bad[:8])))

    pts, layer_of = [], []
    for i, L in enumerate(layers):
        for c in L["contour_3d"]:
            pts.append(np.asarray(c)[:, :3] - offset)       # repere Bird
            layer_of.append(np.full(len(c), i))
    pts, layer_of = np.vstack(pts), np.concatenate(layer_of)
    rng = np.random.default_rng(0)
    if len(pts) > 60000:
        pick = rng.choice(len(pts), 60000, replace=False)
        pts, layer_of = pts[pick], layer_of[pick]
    _c, d_surf, _t = trimesh.proximity.closest_point(surface, pts)
    _c, d_skin, _t = trimesh.proximity.closest_point(retained, pts)

    on_surface = d_surf < TOL_SURF
    under_skin = ~on_surface & (np.abs(d_skin - thickness) <= 0.25 * thickness)
    other = ~on_surface & ~under_skin
    print(f"{len(pts)} points de contour : surface {on_surface.mean():.1%}, sous-peau {under_skin.mean():.1%}, "
          f"autre {other.mean():.1%}")
    if under_skin.any():
        print(f"  sous-peau : distance a la face retenue mediane {np.median(d_skin[under_skin]):.3f} mm "
              f"(visee {thickness} mm)")
    if other.any():
        print(f"  autre : distance a la surface mediane {np.median(d_surf[other]):.3f} mm, 95e centile "
              f"{np.percentile(d_surf[other], 95):.3f} mm ; distance a la peau mediane {np.median(d_skin[other]):.3f} mm")
        worst = np.bincount(layer_of[other], minlength=len(layers))
        top = np.argsort(-worst)[:8]
        print("  couches les plus touchees : " + ", ".join(f"{t + 1} ({worst[t]})" for t in top if worst[t]))
        bands = [(0, 0.1), (0.1, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 1e9)]
        print("  repartition (distance a la surface) : " + ", ".join(
            f"{a}-{b if b < 1e9 else '+'} mm: {((d_surf[other] >= a) & (d_surf[other] < b)).sum()}" for a, b in bands))


if __name__ == "__main__":
    main()
