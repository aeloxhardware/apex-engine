# Profils de machine (R-030 point 4)

Un profil = un fichier `<id>.json`, au format de `core.mesh_types.MachineProfile` (voir `config/params.schema.json`).
Ajouter une machine ou une tête = ajouter un fichier ici (et `res/3d/toolheads/<toolhead_id>/<toolhead_id>.glb`),
sans toucher au code. `index.json` nomme le profil utilisé quand l'appelant n'en désigne aucun.

## aetheris-4 (machine actuelle : R C Z B, pas d'axe Y)
- `rotation_axes` : C (plateau, axe Z, tour complet) puis B (tête, axe X, 0° à 185°).
- `position_coupled_axes = [0]` : C est imposé par la position (limite #59). `radial_direction` = +R à angles
  nuls : B > 0 penche la normale de couche vers -Y, c'est-à-dire vers +R. Après l'upgrade (axe linéaire de plus) :
  `position_coupled_axes = []`.
- `max_angular_velocity` = 1666,7 °/s et `max_angular_acceleration` = 1300 °/s², en rad : `config.g` actif
  (`M203 B100000` /min, `M201 B1300`, confirmé par l'utilisateur le 2026-10-06). Un seul couple pour tous les axes de
  rotation : il ne s'applique qu'à B (C, couplé à la position, n'est jamais lissé). Y (600 °/s) et Z (6,7 mm/s) ne
  sont pas modélisés.
- `platform_radius` = 150 mm : plateau de 300 mm de diamètre (utilisateur, 2026-09-29).
- `toolhead_id` = `aetheris-4` : BIQU H2 V2S Lite, buse longue (`res/3d/toolheads/aetheris-4/`).
