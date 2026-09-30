"""
Mesure du taux de surplomb REEL (pas seulement l'angle de normale) sur un
mesh deforme : pour chaque face dont la normale depasse le seuil d'angle,
teste si elle a bien de la matiere EN DESSOUS via un ray-cast vertical
descendant -- meme principe que pipeline/step2_overhang.py du projet Apex
(limite documentee : un angle de normale seul ne dit rien du contexte
spatial, une face en "surplomb" au sens de l'angle peut tres bien reposer
sur la matiere juste en dessous d'elle).

Usage: python measure_real_overhang.py <fichier.stl> [seuil_deg] [distance_max_support]
"""
import sys
import numpy as np
import trimesh

path = sys.argv[1]
threshold_apex_deg = float(sys.argv[2]) if len(sys.argv) > 2 else 45.0
# Distance max a laquelle un "hit" en dessous compte comme un vrai support --
# sans ce plafond, un impact a l'autre bout du modele (ex: 29mm plus bas sur
# un modele de 51mm de haut) est compte comme "supporte" alors que la face
# est en realite suspendue au sens d'une impression reelle. Meme principe
# que bridge_max_gap_distance dans step2_overhang.py cote Apex.
max_support_distance = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0

mesh = trimesh.load(path)
up = np.array([0.0, 0.0, 1.0])

normals = mesh.face_normals
areas = mesh.area_faces
face_centers = mesh.triangles_center

cos_angle = normals @ up
angle_from_vertical_deg = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
# convention Apex : angle depuis l'HORIZONTALE, donc seuil_S4 = 90 + seuil_Apex
threshold_from_vertical = 90 + threshold_apex_deg

is_overhang_by_angle = angle_from_vertical_deg > threshold_from_vertical
n_overhang_by_angle = is_overhang_by_angle.sum()
area_overhang_by_angle = areas[is_overhang_by_angle].sum()

print(f"Modele : {path}")
print(f"Seuil Apex : {threshold_apex_deg} deg (angle_from_vertical > {threshold_from_vertical} deg)")
print(f"Faces totales : {len(mesh.faces)}, surface totale : {areas.sum():.2f}")
print()
print(f"[TEST ANGLE SEUL] faces en surplomb : {n_overhang_by_angle} ({n_overhang_by_angle/len(mesh.faces)*100:.2f}%)")
print(f"[TEST ANGLE SEUL] surface en surplomb : {area_overhang_by_angle:.2f} ({area_overhang_by_angle/areas.sum()*100:.2f}%)")
print()

if n_overhang_by_angle == 0:
    print("Aucune face en surplomb selon l'angle -- rien a tester par ray-cast.")
    sys.exit(0)

# Ray-cast vertical descendant depuis chaque face en surplomb, meme logique
# que step2_overhang.py : si le rayon ne touche rien avant une distance
# significative, la face est REELLEMENT suspendue dans le vide.
#
# Le decalage doit se faire selon la direction du RAYON (-Z), jamais selon
# la normale de la face -- un decalage selon la normale peut pousser
# l'origine A L'INTERIEUR du solide sur une face inclinee, d'ou l'auto-impact
# immediat detecte comme "support" (bug mesure et corrige : sur un cas de
# controle a vrai porte-a-faux absolu, l'ancienne version rapportait a tort
# 100% de faux positifs).
overhang_indices = np.where(is_overhang_by_angle)[0]
DECAL = 1e-3
origins = face_centers[overhang_indices] + np.array([0, 0, -DECAL])
directions = np.tile(-up, (len(overhang_indices), 1))

locations, index_ray, index_tri = mesh.ray.intersects_location(
    ray_origins=origins, ray_directions=directions, multiple_hits=False
)

# Garde-fou 1 : un impact a distance quasi nulle de l'origine (meme apres
# decalage) est encore un artefact numerique, pas un vrai support.
# Garde-fou 2 (celui qui manquait) : un impact trouve tres LOIN en dessous
# n'est pas non plus un vrai support pour l'impression -- il n'y a que du
# VIDE entre la face et cette matiere, donc la face est bel et bien
# suspendue a l'echelle qui compte. Sans ce plafond, "touche N'IMPORTE QUOI,
# n'importe ou en dessous" est traite a tort comme "est supporte".
MIN_HIT_DISTANCE = 1e-6
hit_ray_indices = set()
hit_distances = {}
for loc, ridx in zip(locations, index_ray):
    dist = np.linalg.norm(loc - origins[ridx])
    if MIN_HIT_DISTANCE < dist <= max_support_distance:
        hit_ray_indices.add(ridx)
        hit_distances[ridx] = min(hit_distances.get(ridx, dist), dist)

truly_suspended = np.array([i not in hit_ray_indices for i in range(len(overhang_indices))])

n_truly_suspended = truly_suspended.sum()
area_truly_suspended = areas[overhang_indices][truly_suspended].sum()
area_false_positive = areas[overhang_indices][~truly_suspended].sum()

print(f"[TEST SUPPORT REEL] parmi les faces en surplomb par angle :")
print(f"  reellement suspendues (rien en dessous) : {n_truly_suspended} ({n_truly_suspended/n_overhang_by_angle*100:.2f}%)")
print(f"  fausses alertes (matiere juste en dessous) : {n_overhang_by_angle - n_truly_suspended} ({(n_overhang_by_angle-n_truly_suspended)/n_overhang_by_angle*100:.2f}%)")
print()
print(f"  surface REELLEMENT en surplomb : {area_truly_suspended:.2f} ({area_truly_suspended/areas.sum()*100:.2f}% du total)")
print(f"  surface faux-positif (supportee) : {area_false_positive:.2f} ({area_false_positive/areas.sum()*100:.2f}% du total)")
