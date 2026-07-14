"""Utilitaires de traitement de maillages : simplification de géométrie 3D."""

import math

try:
    import trimesh
    _HAS_TRIMESH = True
except ImportError:
    _HAS_TRIMESH = False


def simplify_mesh(triangles, ratio):
    """Simplifie un maillage de triangles selon un ratio donné en réduisant la résolution.

    L'objectif est de réduire le nombre de faces tout en conservant la géométrie globale,
    et non de simplement supprimer des triangles au hasard.
    """
    if ratio >= 0.99:
        return triangles

    total = len(triangles)
    if total < 16:
        return triangles

    target_count = int(total * ratio)
    target_count = max(12, target_count)
    if target_count >= total:
        return triangles

    if _HAS_TRIMESH:
        try:
            import trimesh
            vertices = []
            faces = []
            v_map = {}
            for tri in triangles:
                face_idx = []
                for pt in tri:
                    key = (round(pt[0], 6), round(pt[1], 6), round(pt[2], 6))
                    if key not in v_map:
                        v_map[key] = len(vertices)
                        vertices.append(pt)
                    face_idx.append(v_map[key])
                faces.append(face_idx)

            mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
            simplified = mesh.simplify_quadric_decimation(target_count)

            new_triangles = []
            for face in simplified.faces:
                a, b, c = face
                if a == b or b == c or a == c:
                    continue
                new_triangles.append([
                    simplified.vertices[a].tolist(),
                    simplified.vertices[b].tolist(),
                    simplified.vertices[c].tolist()
                ])
            if len(new_triangles) >= max(12, int(total * 0.25)):
                return new_triangles
        except Exception as e:
            print(f"[Simplify] Erreur avec Trimesh: {e}. Utilisation du repli de secours.")

    # Réduction native par fusion de sommets en grille et élimination des faces dégénérées.
    all_points = [pt for tri in triangles for pt in tri]
    xs = [pt[0] for pt in all_points]
    ys = [pt[1] for pt in all_points]
    zs = [pt[2] for pt in all_points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)
    dx = max(max_x - min_x, 1e-6)
    dy = max(max_y - min_y, 1e-6)
    dz = max(max_z - min_z, 1e-6)

    # Taille de cellule approximative pour obtenir une densité raisonnable.
    volume = dx * dy * dz
    density = max(1.0, (target_count / total) * 8.0)
    cell_size = max(1e-4, math.pow(volume / max(target_count, 1), 1.0 / 3.0) * 1.2 / density)

    cluster_map = {}
    cluster_vertices = {}
    vertex_map = {}
    next_id = 0

    for pt in all_points:
        key = (
            int((pt[0] - min_x) / cell_size),
            int((pt[1] - min_y) / cell_size),
            int((pt[2] - min_z) / cell_size)
        )
        if key not in cluster_vertices:
            cluster_vertices[key] = [0.0, 0.0, 0.0, 0]
        cluster_vertices[key][0] += pt[0]
        cluster_vertices[key][1] += pt[1]
        cluster_vertices[key][2] += pt[2]
        cluster_vertices[key][3] += 1

    for key, sums in cluster_vertices.items():
        count = sums[3]
        cluster_map[key] = (sums[0] / count, sums[1] / count, sums[2] / count)
        vertex_map[key] = next_id
        next_id += 1

    def quantize(pt):
        return (
            int((pt[0] - min_x) / cell_size),
            int((pt[1] - min_y) / cell_size),
            int((pt[2] - min_z) / cell_size)
        )

    simplified = []
    for tri in triangles:
        ids = [vertex_map.get(quantize(pt)) for pt in tri]
        if None in ids:
            continue
        if ids[0] == ids[1] or ids[1] == ids[2] or ids[0] == ids[2]:
            continue
        simplified.append([
            cluster_map[quantize(tri[0])],
            cluster_map[quantize(tri[1])],
            cluster_map[quantize(tri[2])]
        ])

    if len(simplified) < 12 or len(simplified) >= total:
        return triangles

    # Si on n'a pas suffisamment réduit, retirer quelques triangles en fin de liste.
    if len(simplified) > target_count:
        return simplified[:target_count]

    return simplified
