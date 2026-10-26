"""Export VTK (legacy ASCII ou binaire big-endian) d'un maillage tétraédrique, lisible par ParaView.

    from vtk_export import write_vtk
    write_vtk("out.vtk", nodes, tets, point_data={"phi": phi}, cell_data={"quality": q})
"""
import numpy as np


def write_vtk(path, nodes, tets, point_data=None, cell_data=None, title="TetraNet"):
    nodes = np.asarray(nodes, dtype=">f8").reshape(-1, 3)
    tets = np.asarray(tets, dtype=">i4").reshape(-1, 4)
    T = len(tets)
    with open(path, "wb") as f:
        f.write(f"# vtk DataFile Version 3.0\n{title}\nBINARY\nDATASET UNSTRUCTURED_GRID\n".encode())
        f.write(f"POINTS {len(nodes)} double\n".encode())
        f.write(nodes.tobytes())
        f.write(f"\nCELLS {T} {5 * T}\n".encode())
        f.write(np.hstack([np.full((T, 1), 4, dtype=">i4"), tets]).tobytes())
        f.write(f"\nCELL_TYPES {T}\n".encode())
        f.write(np.full(T, 10, dtype=">i4").tobytes())
        for header, data, n in (("POINT_DATA", point_data, len(nodes)), ("CELL_DATA", cell_data, T)):
            if not data:
                continue
            f.write(f"\n{header} {n}\n".encode())
            for name, arr in data.items():
                arr = np.asarray(arr, dtype=">f8").reshape(n)
                f.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n".encode())
                f.write(arr.tobytes())
                f.write(b"\n")


def tet_quality(nodes, tets):
    """Volume signé, rapport d'arêtes (max/min) et angle dièdre minimal (degrés) par tétraèdre."""
    nodes = np.asarray(nodes, float)
    tets = np.asarray(tets)
    p = nodes[tets]
    vol = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]) / 6.0
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    L = np.stack([np.linalg.norm(p[:, a] - p[:, b], axis=1) for a, b in pairs], axis=1)
    # angles dièdres : pour chaque arête (a,b), faces opposées (c) et (d)
    dih = []
    for a, b in pairs:
        c, d = [i for i in range(4) if i not in (a, b)]
        e = p[:, b] - p[:, a]
        e = e / np.maximum(np.linalg.norm(e, axis=1, keepdims=True), 1e-300)
        u = p[:, c] - p[:, a]; u = u - np.einsum("ij,ij->i", u, e)[:, None] * e
        v = p[:, d] - p[:, a]; v = v - np.einsum("ij,ij->i", v, e)[:, None] * e
        cs = np.einsum("ij,ij->i", u, v) / np.maximum(np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1), 1e-300)
        dih.append(np.degrees(np.arccos(np.clip(cs, -1, 1))))
    dih = np.stack(dih, axis=1)
    return vol, L.max(axis=1) / np.maximum(L.min(axis=1), 1e-300), dih.min(axis=1)
