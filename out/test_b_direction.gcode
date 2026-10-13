; Apex -- essai A SEC du sens de B (sans extrusion, sans chauffe)
; pointe immobile en R=60 mm, Z=40 mm ; pivot a 46 mm de la pointe
; b_sign = 1 : B firmware = 1 x B d'Apex
; A VERIFIER : 1) la pointe ne bouge pas ; 2) B positif penche le corps de la tete vers +R (buse vers le centre)
G21
G90
G00 X60.0000 Y0.0000 Z40.0000 B0 F1500 ; B Apex 0
G00 X75.7329 Y0.0000 Z37.2259 B20 F1500 ; B Apex 20
G00 X89.5682 Y0.0000 Z29.2380 B40 F1500 ; B Apex 40
G00 X99.8372 Y0.0000 Z17.0000 B60 F1500 ; B Apex 60
G00 X89.5682 Y0.0000 Z29.2380 B40 F1500 ; B Apex 40
G00 X75.7329 Y0.0000 Z37.2259 B20 F1500 ; B Apex 20
G00 X60.0000 Y0.0000 Z40.0000 B0 F1500 ; B Apex 0
; fin
