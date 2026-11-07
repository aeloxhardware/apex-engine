"""
Lissage non-planaire sur la voie CHALEUR (chantier ouvert le 2026-10-28).

Reprise, adaptee au champ φ de la methode de la chaleur, de la chaine de lissage de la voie Bird
(`pipeline/bird/skin/`, branche `bird_old`) :
  - `zones`  : P2, quelles faces de la surface d'origine peuvent recevoir une couche de surface ;
  - `carve`  : P3, creusement du corps sous les zones retenues (decoupe des isosurfaces de φ le long de ψ = 0).
Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer).
"""
