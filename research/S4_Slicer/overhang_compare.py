"""Part de la surface d'origine en surplomb (> 45 deg vers le bas par rapport a
la direction d'empilement = gradient de phi du tetraedre porteur de la face),
faces posees sur le plateau exclues. usage : overhang_compare.py <modele> <passe> <variante>..."""
import sys, numpy as np
sys.path.insert(0, ".")
from bird_slice_deformed import boundary_faces
from bird_to_apex import tet_gradients
model, k = sys.argv[1], sys.argv[2]
for v in ["plan"] + sys.argv[3:]:
    pre = "output_models/loop_" + model.replace(" ", "_") + (sys.argv[3] if v == "plan" else v)
    cells, p0 = np.load(pre + "_cells.npy"), np.load(pre + "_pass0.npy")
    faces = boundary_faces(cells, p0)
    tri = p0[faces]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); area = np.linalg.norm(n, axis=1) / 2
    n /= np.maximum(2 * area, 1e-12)[:, None]
    # tetraedre porteur de chaque face : celui qui contient ses 3 sommets
    key = {tuple(sorted(f)): i for i, f in enumerate(faces)}
    owner = np.zeros(len(faces), int)
    for t, c in enumerate(cells):
        for f in ((c[0], c[1], c[2]), (c[0], c[1], c[3]), (c[0], c[2], c[3]), (c[1], c[2], c[3])):
            i = key.get(tuple(sorted(f)))
            if i is not None: owner[i] = t
    if v == "plan":
        d = np.tile([0, 0, 1.0], (len(faces), 1))
    else:
        g = tet_gradients(p0, cells, np.load(f"{pre}_pass{k}.npy")[:, 2]); d = g[owner]
        d /= np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-12)
    plate = tri[:, :, 2].max(1) < p0[:, 2].min() + 0.25
    over = ((n * -d).sum(1) > np.cos(np.radians(45))) & ~plate
    print(f"{v:12s}: surplomb {area[over].sum() / area[~plate].sum():6.2%} de la surface ({area[over].sum():.0f} mm2)")
