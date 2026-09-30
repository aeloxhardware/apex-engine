"""Generateur G-code AUTONOME pour la machine Apex 4 axes.

But : sortir tout de suite un G-code imprimable en forme de L, sans passer
par le pipeline geometrique (step1..step10). N'importe RIEN du projet.

Geometrie : un L en forme de barreau CIRCULAIRE (section ronde, pas carree)
-- choix delibere pour exploiter la cinematique polaire de la machine :
un cercle CENTRE SUR L'AXE DE ROTATION (R=0) se trace en maintenant R
CONSTANT pendant que C balaie, le mouvement le plus natif qui soit pour
cette machine (aucune approximation/subdivision necessaire pour un cercle
centre, contrairement a un carre -- cf. tilted_frame_point). Contour seul
(perimetre exterieur, milieu creux), sauf la toute premiere couche pleine
(remplissage en spirale concentrique, pour l'adherence au plateau).

- Branche 1 (pilier) : monte en Z, tete vers le bas (B=0). Couches 1-4
  discretes (adherence), le reste en spirale continue ("vase mode").
- Coude : la section se REORIENTE progressivement (le CENTRE de la section
  avance dans la direction de depot courante d(beta), pas seulement la
  tete qui s'incline sur place) -- couches a epaisseur variable en Z.
- Branche 2 (surplomb) : B constant (OVERHANG_ANGLE_DEG), section inclinee
  du meme angle, spine avancant purement en R (le barreau reste
  geometriquement horizontal).

Sortie firmware : G00 Y<C> X<R> Z<Z> B<B> F<feedrate>  (+ E pour l'extrusion)
  X (firmware "X") = rayon plateau, mm
  Y (firmware "Y") = rotation plateau, degres
  Z              = hauteur, mm
  B              = rotation de la tete, deg (0 = vers le bas, 90 = horizontale)

Usage : python scripts/gen_l_overhang_gcode.py [chemin_sortie]
Defaut : out/l_overhang.gcode
"""

from __future__ import annotations

import math
import os
import sys

# --------------------------------------------------------------------------
# Parametres
# --------------------------------------------------------------------------

CIRCLE_DIAMETER = 25.0      # diametre de la section circulaire, mm
PILLAR_HEIGHT = 30.0        # hauteur de la branche verticale avant le surplomb, mm
OVERHANG_LENGTH = 60.0      # longueur de la branche en surplomb, mm
OVERHANG_ANGLE_DEG = 45.0   # angle de tete pendant l'impression du surplomb, deg
                            # (0=vers le bas, 90=horizontal) ; 45 = seuil
                            # classique d'auto-support, la tige reste
                            # geometriquement horizontale (pilotee en R),
                            # seule l'inclinaison de depot change
LAYER_HEIGHT = 0.2          # hauteur de couche, mm
EXTRUSION_WIDTH = 0.4       # largeur de cordon, mm
FILAMENT_DIAMETER = 1.75    # mm
NOZZLE_TEMP = 195           # deg C
# Vitesses imposees par l'utilisateur (erreur de config assumee) :
#   couches lentes a 100 mm/s, autres couches 300 mm/s. F est en mm/min.
PRINT_FEEDRATE = 300 * 60        # 18000 mm/min = 300 mm/s
TRAVEL_FEEDRATE = 300 * 60       # 18000 mm/min (deplacement a vide)
FIRST_LAYER_FEEDRATE = 100 * 60  # 6000 mm/min = 100 mm/s (couches lentes)
SOLID_BASE_LAYER_COUNT = 4       # nombre de couches LENTES et DISCRETES (pas de
                                  # spirale continue) en debut de branche verticale,
                                  # imprimees a FIRST_LAYER_FEEDRATE. SEULE la toute
                                  # premiere a un remplissage (spirale concentrique,
                                  # adherence au plateau) ; les suivantes ne sont que
                                  # le contour. cf. emit_solid_base_layers.

# Longueur (le long du spine, pas Z pur) sur laquelle B passe de 0 a
# OVERHANG_ANGLE_DEG. DOIT scaler avec CIRCLE_DIAMETER, pas une constante
# figee : le spine du coude ne s'eloigne de l'axe R=0 que de
# ELBOW_Z_SPAN * moyenne(sin(beta), 0..OVERHANG_ANGLE_DEG) ; pour que la
# section (rayon CIRCLE_RADIUS) cesse d'enjamber l'axe a la fin du coude
# (sinon elle continue de tracer un tour complet a chaque couche, meme
# inclinee -- c'est le bug "continue de tourner dans le vide"), il faut
# ELBOW_Z_SPAN * moy(sin(beta)) >= CIRCLE_RADIUS*cos(OVERHANG_ANGLE_DEG).
# Mesure (avec l'ancienne section carree, meme ordre de grandeur) :
# ELBOW_Z_SPAN=CIRCLE_DIAMETER donne une marge de ~5% au-dessus du seuil
# exact ; x1.1 donne une marge confortable (~18%) sans trop rallonger la
# transition. Un CERCLE n'a pas de coin qui deborde plus loin qu'un autre
# point (contrairement a la diagonale d'un carre) : cette marge est donc
# au moins aussi bonne, jamais pire, qu'avec la section carree d'origine.
ELBOW_Z_SPAN = CIRCLE_DIAMETER * 1.1

# Section CIRCULAIRE centree sur l'axe de rotation du plateau (R=0) :
# exploite pleinement la cinematique polaire -- un cercle centre se trace
# a RAYON MACHINE CONSTANT pendant que C balaie, le mouvement le plus
# natif possible pour cette machine (voir generate_circle_loop).
CIRCLE_RADIUS = CIRCLE_DIAMETER / 2.0

# --------------------------------------------------------------------------
# Extrusion
# --------------------------------------------------------------------------

_FILAMENT_AREA = math.pi * (FILAMENT_DIAMETER / 2.0) ** 2
# section deposee par mm de trajet = hauteur de couche * largeur de cordon
_BEAD_AREA = LAYER_HEIGHT * EXTRUSION_WIDTH


def e_for_length(path_len_mm: float) -> float:
    """Longueur de filament (mm) a extruder pour path_len_mm de trajet."""
    return path_len_mm * _BEAD_AREA / _FILAMENT_AREA


# --------------------------------------------------------------------------
# Helpers cinematique
# --------------------------------------------------------------------------

def xy_to_rc(x: float, y: float, prev_c_deg: float) -> tuple[float, float]:
    """Cartesien (x, y) machine -> (R, C degres).

    R = rayon depuis le centre du plateau. C = angle du plateau.
    Si (x, y) tombe pile au centre (R ~ 0), C est indetermine : on herite
    le C du point precedent (cf. limite #19 du projet).
    """
    r = math.hypot(x, y)
    if r < 1e-6:
        return 0.0, prev_c_deg
    c = math.degrees(math.atan2(y, x))
    return r, c


def polyline_length(points: list[tuple[float, float]]) -> float:
    total = 0.0
    for i in range(1, len(points)):
        total += math.hypot(points[i][0] - points[i - 1][0],
                            points[i][1] - points[i - 1][1])
    return total


# Distance max (mm) entre deux points consecutifs d'une trajectoire, une
# fois convertis en (R, C). Necessaire des que la trajectoire n'est PAS un
# cercle a rayon machine constant (donc dans le coude et le surplomb, ou le
# centre de la section s'est deplace hors de R=0) : le firmware interpole
# R et C lineairement entre deux points commandes, ce qui ne suit une
# courbe quelconque que si les points sont assez rapproches. Un cercle
# centre sur l'axe (le pilier, R=0) n'a PAS besoin de cette subdivision
# pour etre exact (voir generate_circle_loop) -- on la garde quand meme
# partout par simplicite et pour un maillage de points homogene.
MAX_CHORD_MM = 1.0


def subdivide_loop_xy(loop_xy: list[tuple[float, float]],
                       max_chord: float = MAX_CHORD_MM,
                       ) -> list[tuple[float, float]]:
    """Subdivise chaque segment d'une boucle fermee en corde <= max_chord.

    Necessaire des qu'on va convertir les points en (R, C) : le firmware
    interpole R et C lineairement entre deux points commandes, ce qui ne
    correspond a une ligne droite CARTESIENNE que si les points sont assez
    rapproches. loop_xy doit deja etre ferme (dernier point == premier).
    """
    out = [loop_xy[0]]
    for i in range(1, len(loop_xy)):
        x0, y0 = loop_xy[i - 1]
        x1, y1 = loop_xy[i]
        seg_len = math.hypot(x1 - x0, y1 - y0)
        n = max(1, math.ceil(seg_len / max_chord))
        for k in range(1, n + 1):
            t = k / n
            out.append((x0 + (x1 - x0) * t, y0 + (y1 - y0) * t))
    return out


def tilted_frame_point(n: float, t: float, beta_deg: float,
                        r_center: float, z_center: float,
                        ) -> tuple[float, float, float]:
    """Point (n,t) d'une section, dans le repere INCLINE de beta.

    Generique : (n,t) peuvent parametrer un cercle (generate_circle_loop)
    ou n'importe quelle autre forme locale 2D -- cette fonction ne
    connait pas la forme de la section, seulement son inclinaison.

    n = position le long de la direction PERPENDICULAIRE a la direction de
        depot courante, dans le plan (R,Z) -- c'est la direction dans
        laquelle la SECTION s'etend "en profondeur" (a beta=0, n devient un
        simple decalage radial ; a beta=90, un decalage en -Z).
    t = position tangentielle (perpendiculaire a n ET a la direction de
        depot), traitee comme un decalage cartesien Y -- approximation
        petit angle deja utilisee ailleurs dans ce script (le surplomb
        n'a pas de vrai axe lineaire tangentiel sur cette machine).
    (r_center, z_center) = centre de la section (le long du "spine"
        d'impression).

    A beta=0, ceci degenere EXACTEMENT vers la section plate du pilier
    (X=r_center+n, Y=t, Z=z_center constant) -- pas une approximation a
    cet endroit precis, la jonction pilier/coude est donc exacte.
    """
    beta = math.radians(beta_deg)
    x = r_center + n * math.cos(beta)
    y = t
    z = z_center - n * math.sin(beta)
    return x, y, z


def generate_circle_loop(radius: float, max_chord: float = MAX_CHORD_MM,
                          ) -> list[tuple[float, float]]:
    """Boucle fermee (n,t) d'un CERCLE de rayon radius, centre a l'origine.

    Remplace la section carree d'origine : exploite pleinement la
    cinematique polaire. Un cercle centre sur l'axe de rotation (R=0, le
    cas du pilier a beta=0) se trace a RAYON MACHINE CONSTANT pendant que
    C balaie -- AUCUNE approximation de corde n'est meme necessaire dans
    ce cas precis (contrairement a un carre, qui degenerait en cercle
    sans subdivision -- ici on VEUT un cercle). La subdivision fine est
    quand meme appliquee pour un maillage de points homogene et parce que
    le coude/surplomb (centre hors R=0) en ont, eux, reellement besoin.
    """
    n_points = max(8, math.ceil(2.0 * math.pi * radius / max_chord))
    return [(radius * math.cos(2.0 * math.pi * i / n_points),
             radius * math.sin(2.0 * math.pi * i / n_points))
            for i in range(n_points + 1)]  # +1 : referme la boucle


def generate_spiral_infill(radius: float, extrusion_width: float,
                            max_chord: float = MAX_CHORD_MM,
                            ) -> list[tuple[float, float]]:
    """Remplissage plein en SPIRALE CONCENTRIQUE (Archimede) d'un disque.

    r(theta) = theta * extrusion_width / (2*pi) : le pas radial par tour
    complet vaut exactement extrusion_width (jointif, remplissage 100%,
    sans le vide/recouvrement qu'un pas mal choisi laisserait). La spirale
    part du CENTRE (r=0, la singularite R=0/C indetermine) -- c'est
    INEVITABLE pour un remplissage sans trou au centre d'un disque, et
    xy_to_rc gere deja ce cas (r<1e-6 -> herite le C precedent).

    Bien plus naturel que le zigzag pour cette cinematique : R croit en
    continu pendant que C balaie en continu, jamais d'aller-retour
    cartesien. Pas angulaire ADAPTATIF (fin pres du centre ou r est petit,
    plus large loin du centre) pour garder une corde <= max_chord partout,
    sans exploser le nombre de points pres du centre.
    """
    pitch = extrusion_width / (2.0 * math.pi)  # dr/dtheta
    theta_max = radius / pitch

    points: list[tuple[float, float]] = [(0.0, 0.0)]
    theta = 0.0
    # rayon plancher pour le calcul du pas angulaire : evite une division
    # par ~0 tout pres du centre (le premier pas resterait borne et sense).
    min_r_for_step = extrusion_width / 2.0
    while theta < theta_max:
        r = theta * pitch
        d_theta = max_chord / max(r, min_r_for_step)
        theta = min(theta + d_theta, theta_max)
        r = theta * pitch
        points.append((r * math.cos(theta), r * math.sin(theta)))
    return points


# --------------------------------------------------------------------------
# Generation du G-code
# --------------------------------------------------------------------------

class GCodeWriter:
    def __init__(self) -> None:
        self.lines: list[str] = []
        self.e = 0.0            # position absolue de l'extrudeur
        self.last_c = 0.0       # dernier C emis (pour la singularite)
        self.last_pos = None    # (R, C, Z, B)

    def raw(self, line: str = "") -> None:
        self.lines.append(line)

    def comment(self, text: str) -> None:
        self.lines.append(f"; {text}")

    def move(self, *, r: float, c: float, z: float, b: float,
             feed: float, extrude_len: float | None = None) -> None:
        """Emet un G00 Y<C> X<R> Z<Z> B<B> F<feed> [E<e>].

        extrude_len = longueur de trajet a extruder (mm) ; None = deplacement
        a vide (pas de E).
        """
        parts = ["G00"]
        parts.append(f"Y{c:.4f}")
        parts.append(f"X{r:.4f}")
        parts.append(f"Z{z:.4f}")
        parts.append(f"B{b:.4f}")
        if extrude_len is not None and extrude_len > 0.0:
            self.e += e_for_length(extrude_len)
            parts.append(f"E{self.e:.5f}")
        parts.append(f"F{feed:.0f}")
        self.lines.append(" ".join(parts))
        self.last_c = c
        self.last_pos = (r, c, z, b)

    def travel_to_xy(self, x: float, y: float, z: float, b: float) -> None:
        r, c = xy_to_rc(x, y, self.last_c)
        # Continuite angulaire, comme pour les points d'extrusion : sans ce
        # recalage, ce travel ramenerait C dans [-180,180] brut, c'est-a-dire
        # ferait faire au plateau un tour complet EN ARRIERE a chaque
        # changement de couche pour "revenir" a un point ou il se trouve
        # deja (a un multiple de 360 pres). En le recalant au plus proche de
        # la position courante, le plateau continue de tourner dans le meme
        # sens sur toute la branche : le travel entre deux couches devient
        # un simple lift en Z, sans aucune rotation inutile. C'est le point
        # meme de la cinematique polaire par rapport au cartesien.
        c = _closest_angle(c, self.last_c)
        self.move(r=r, c=c, z=z, b=b, feed=TRAVEL_FEEDRATE)

    def dump(self) -> str:
        return "\n".join(self.lines) + "\n"


def emit_header(g: GCodeWriter) -> None:
    g.comment("=== L overhang - genere par scripts/gen_l_overhang_gcode.py ===")
    g.comment("Machine Apex 4 axes : G00 Y<C> X<R> Z<Z> B<B> F")
    g.comment(f"section circulaire diametre {CIRCLE_DIAMETER} mm (contour "
              f"seul + spirale couche 1), pilier {PILLAR_HEIGHT} mm, "
              f"surplomb {OVERHANG_LENGTH} mm a {OVERHANG_ANGLE_DEG} deg, "
              f"hauteur totale {PILLAR_HEIGHT + CIRCLE_DIAMETER} mm")
    g.comment(f"layer_height {LAYER_HEIGHT} mm, extrusion_width "
              f"{EXTRUSION_WIDTH} mm")
    g.raw()
    g.raw("G21")            # unites mm
    g.raw("G90")            # positionnement absolu
    g.raw("M83")            # ... mais extrudeur en relatif ? non : on gere E absolu
    g.raw("G92 E0")         # reset extrudeur
    g.raw(f"M109 S{NOZZLE_TEMP}")   # chauffe buse ET attend (pas de M104)
    g.raw("M82")            # extrudeur absolu (coherent avec self.e)
    g.raw("G92 E0")
    g.raw("G00 B0.0000 F3000")   # tete vers le bas
    g.raw("G00 Y0.0000 X0.0000 Z5.0000 F3000")  # au centre, au-dessus du plateau
    g.raw()


def emit_footer(g: GCodeWriter) -> None:
    g.raw()
    g.comment("fin")
    g.raw("G00 Z90.0000 F3000")   # remonte
    g.raw("M84")                  # moteurs off


def _emit_xy_polyline(g: GCodeWriter, points_xy: list[tuple[float, float]],
                       z: float, b: float, feed: float) -> None:
    """Emet une polyligne (x,y) subdivisee, a Z et B fixes (couche discrete).

    points_xy ne doit PAS inclure le point de depart courant : le premier
    point de points_xy est rejoint par un travel_to_xy (a vide), puis
    chaque segment suivant est extrude.
    """
    x0, y0 = points_xy[0]
    g.travel_to_xy(x0, y0, z, b=b)
    for k in range(1, len(points_xy)):
        x, y = points_xy[k]
        xp, yp = points_xy[k - 1]
        seg_len = math.hypot(x - xp, y - yp)
        r, c = xy_to_rc(x, y, g.last_c)
        c = _closest_angle(c, g.last_c)
        g.move(r=r, c=c, z=z, b=b, feed=feed, extrude_len=seg_len)


def emit_solid_base_layers(g: GCodeWriter) -> None:
    """Les SOLID_BASE_LAYER_COUNT premieres couches : discretes, B=0.

    Par souci d'adherence au plateau, ces couches ne sont pas la spirale
    "vase mode" du reste du pilier : ce sont de vraies couches a Z fixe.
    SEULE la toute premiere a un remplissage -- une spirale concentrique
    (generate_spiral_infill), le mode de remplissage le plus naturel pour
    cette cinematique (R croit en continu pendant que C balaie). Les
    couches suivantes (2..SOLID_BASE_LAYER_COUNT) ne sont que le contour :
    la premiere couche a deja fait adherer toute la surface au plateau,
    inutile de remplir les couches du dessus.

    Le contour (cercle) est trace en dernier a chaque couche, pour que la
    jonction avec la spirale qui suit (emit_vertical_branch) parte du
    meme point.
    """
    g.comment(f"--- {SOLID_BASE_LAYER_COUNT} couches de base (adherence), B=0 ---")
    loop_xy = generate_circle_loop(CIRCLE_RADIUS)
    infill_xy = generate_spiral_infill(CIRCLE_RADIUS, EXTRUSION_WIDTH)

    for i in range(1, SOLID_BASE_LAYER_COUNT + 1):
        z = i * LAYER_HEIGHT
        g.comment(f"couche de base {i}/{SOLID_BASE_LAYER_COUNT}  Z={z:.2f}")
        if i == 1:
            _emit_xy_polyline(g, infill_xy, z, b=0.0, feed=FIRST_LAYER_FEEDRATE)
        _emit_xy_polyline(g, loop_xy, z, b=0.0, feed=FIRST_LAYER_FEEDRATE)


def emit_vertical_branch(g: GCodeWriter) -> None:
    """Reste du pilier, VASE MODE : une seule spirale continue, B=0.

    Part de Z=SOLID_BASE_LAYER_COUNT*LAYER_HEIGHT (fin des couches de base,
    emit_solid_base_layers) jusqu'a Z=PILLAR_HEIGHT. Pas de "couches" au
    sens palier : Z monte en continu, proportionnel a la fraction de tour
    parcourue le long du contour, jamais par saut. Le point de depart
    est deja la ou emit_solid_base_layers s'est arrete (meme point de
    depart du cercle) : pas de travel_to_xy ici, la continuite est
    geometrique.

    Le cercle est centre sur l'axe (R=0) : chaque point est a
    RAYON MACHINE CONSTANT (=CIRCLE_RADIUS), seul C balaie -- la spirale
    entiere est donc un simple HELICE (R constant, C et Z montant tous
    les deux en continu), le mouvement le plus natif pour cette machine.
    """
    z_start = SOLID_BASE_LAYER_COUNT * LAYER_HEIGHT
    z_span = PILLAR_HEIGHT - z_start
    n_layers = int(round(z_span / LAYER_HEIGHT))
    g.comment("--- reste du pilier (helice continue), B=0 ---")
    loop_xy = generate_circle_loop(CIRCLE_RADIUS)
    # le dernier point de loop_xy == le premier (boucle fermee) : ne pas le
    # recompter comme un pas de tour distinct, sinon Z avancerait deux fois
    # au meme endroit a la jonction entre deux tours.
    steps_per_turn = len(loop_xy) - 1
    total_steps = steps_per_turn * n_layers

    step_index = 0
    for turn in range(n_layers):
        for k in range(1, len(loop_xy)):
            step_index += 1
            x, y = loop_xy[k]
            xp, yp = loop_xy[k - 1]
            seg_len = math.hypot(x - xp, y - yp)
            r, c = xy_to_rc(x, y, g.last_c)
            # continuite angulaire : on ne veut pas de saut de +/-360
            c = _closest_angle(c, g.last_c)
            # Z proportionnel a la fraction de spirale deja parcourue, PAS
            # a "quelle couche" : c'est ce qui fait une rampe continue au
            # lieu d'un palier.
            z = z_start + (step_index / total_steps) * z_span
            g.move(r=r, c=c, z=z, b=0.0, feed=PRINT_FEEDRATE, extrude_len=seg_len)


def _closest_angle(target_deg: float, ref_deg: float) -> float:
    """Ramene target a moins de 180 deg de ref (evite un saut C de 360)."""
    while target_deg - ref_deg > 180.0:
        target_deg -= 360.0
    while target_deg - ref_deg < -180.0:
        target_deg += 360.0
    return target_deg


def emit_elbow(g: GCodeWriter) -> tuple[float, float]:
    """Coude : la section se REORIENTE progressivement, pas juste B.

    Contrairement a une version anterieure qui ne faisait QUE monter B en
    retracant le meme cercle plat centre sur l'axe (donc tournant en pleine
    rotation, "dans le vide", meme une fois la tete inclinee) : ici le
    CENTRE de la section (le "spine" d'impression) avance a chaque couche
    dans la direction de depot courante d(beta) = (sin(beta), cos(beta))
    en (R,Z). A beta=0 le spine monte pur +Z (identique au pilier) ; a
    mesure que beta croit vers OVERHANG_ANGLE_DEG, le spine part de plus
    en plus en R et de moins en moins en Z -- EXACTEMENT les couches a
    epaisseur variable en Z demandees (dZ_par_couche = LAYER_HEIGHT*cos(beta),
    strictement decroissant).

    Consequence organique (pas un cas particulier code en dur) : des que
    le spine s'eloigne de l'axe R=0 au-dela du rayon de la section, celle-ci
    ne "boucle" plus autour de l'axe -- le balayage C se resserre de
    lui-meme en oscillation bornee au lieu d'un tour complet. C'est la
    reponse directe a "continue de tourner dans le vide".

    Retourne (r_end, z_end) : la position REELLE ou le coude s'arrete,
    pour que la branche en surplomb reparte exactement de la, au lieu
    d'une position arbitraire recalculee independamment.
    """
    g.comment(f"--- coude (spine oriente), B 0 -> {OVERHANG_ANGLE_DEG} deg ---")
    n_layers = int(round(ELBOW_Z_SPAN / LAYER_HEIGHT))
    loop_nt = generate_circle_loop(CIRCLE_RADIUS)

    r_center = 0.0
    z_center = PILLAR_HEIGHT
    # point de depart = exactement la ou le pilier s'est arrete (beta=0,
    # r_center=0 degenere vers le cercle plat, cf. docstring de
    # tilted_frame_point) : pas un recalage arbitraire, jonction exacte.
    prev_x, prev_y, prev_z = tilted_frame_point(
        loop_nt[0][0], loop_nt[0][1], 0.0, r_center, z_center)

    for turn in range(n_layers):
        beta_start = OVERHANG_ANGLE_DEG * turn / n_layers
        beta_end = OVERHANG_ANGLE_DEG * (turn + 1) / n_layers
        beta_mid = (beta_start + beta_end) / 2.0

        # avance du spine dans la direction de depot COURANTE (au milieu
        # du pas, pour rester centre sur l'intervalle) : c'est ce qui fait
        # une couche a epaisseur variable en Z plutot qu'un simple B
        # cosmetique.
        beta_mid_rad = math.radians(beta_mid)
        r_center += LAYER_HEIGHT * math.sin(beta_mid_rad)
        z_center += LAYER_HEIGHT * math.cos(beta_mid_rad)

        # la section de CETTE couche est inclinee a beta_end (orientation
        # de tete a la fin du pas, coherent avec le sens de B qui monte).
        # k=0 == dernier point du tour precedent (boucle fermee) : jamais
        # reemis, meme convention que emit_vertical_branch.
        for k in range(1, len(loop_nt)):
            n_k, t_k = loop_nt[k]
            x, y, z = tilted_frame_point(n_k, t_k, beta_end, r_center, z_center)
            seg_len = math.sqrt((x - prev_x) ** 2 + (y - prev_y) ** 2
                                 + (z - prev_z) ** 2)
            r, c = xy_to_rc(x, y, g.last_c)
            c = _closest_angle(c, g.last_c)
            g.move(r=r, c=c, z=z, b=beta_end, feed=PRINT_FEEDRATE,
                   extrude_len=seg_len)
            prev_x, prev_y, prev_z = x, y, z

    return r_center, z_center


def emit_horizontal_branch(g: GCodeWriter, r_start: float, z_center: float) -> None:
    """Branche en surplomb : B=OVERHANG_ANGLE_DEG CONSTANT, section INCLINEE.

    (r_start, z_center) viennent du point d'arrivee REEL du coude
    (emit_elbow) : la jonction est exacte, pas une position recalculee
    independamment.

    Beta est desormais CONSTANT (le coude a deja fini de tourner la tete) :
    la tige reste GEOMETRIQUEMENT horizontale, donc le spine avance
    PUREMENT en R, a Z CONSTANT -- contrairement au coude ou le spine
    suivait la direction de depot d(beta) (qui a une composante Z tant que
    beta varie). La SECTION, elle, reste inclinee a beta=OVERHANG_ANGLE_DEG
    (meme mapping tilted_frame_point que le coude, continuite garantie).

    Pas radial entre deux couches : LAYER_HEIGHT / sin(beta), PAS
    LAYER_HEIGHT tout court. Justification : l'epaisseur REELLE de couche
    (mesuree perpendiculairement au plan incline) doit rester LAYER_HEIGHT ;
    or avancer le spine de dR en restant a Z constant ne deplace le plan
    incline (dont la normale a une composante sin(beta) selon R) que de
    dR*sin(beta) perpendiculairement a lui-meme. Poser dR*sin(beta) =
    LAYER_HEIGHT donne dR = LAYER_HEIGHT/sin(beta).
    """
    beta_rad = math.radians(OVERHANG_ANGLE_DEG)
    dr_per_layer = LAYER_HEIGHT / math.sin(beta_rad)
    n_layers = int(round(OVERHANG_LENGTH / dr_per_layer))

    g.comment(f"--- surplomb (section inclinee {OVERHANG_ANGLE_DEG} deg, "
              f"B constant), R avance de {dr_per_layer:.4f} mm/couche ---")
    loop_nt = generate_circle_loop(CIRCLE_RADIUS)

    r_center = r_start
    prev_x, prev_y, prev_z = tilted_frame_point(
        loop_nt[0][0], loop_nt[0][1], OVERHANG_ANGLE_DEG, r_center, z_center)

    for turn in range(n_layers):
        r_center += dr_per_layer  # spine PUR en R, Z_center inchange
        # k=0 == dernier point du tour precedent (boucle fermee) : jamais
        # reemis, meme convention que le coude/pilier.
        for k in range(1, len(loop_nt)):
            n_k, t_k = loop_nt[k]
            x, y, z = tilted_frame_point(n_k, t_k, OVERHANG_ANGLE_DEG,
                                          r_center, z_center)
            seg_len = math.sqrt((x - prev_x) ** 2 + (y - prev_y) ** 2
                                 + (z - prev_z) ** 2)
            r, c = xy_to_rc(x, y, g.last_c)
            c = _closest_angle(c, g.last_c)
            g.move(r=r, c=c, z=z, b=OVERHANG_ANGLE_DEG, feed=PRINT_FEEDRATE,
                   extrude_len=seg_len)
            prev_x, prev_y, prev_z = x, y, z


def main() -> None:
    out_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join("out", "l_overhang.gcode")
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    g = GCodeWriter()
    emit_header(g)
    emit_solid_base_layers(g)
    emit_vertical_branch(g)
    r_end, z_end = emit_elbow(g)
    emit_horizontal_branch(g, r_start=r_end, z_center=z_end)
    emit_footer(g)

    text = g.dump()
    with open(out_path, "w", encoding="ascii") as f:
        f.write(text)

    n_moves = sum(1 for ln in g.lines if ln.startswith("G00") and "E" in ln)
    print(f"ecrit : {out_path}")
    print(f"  {len(g.lines)} lignes, {n_moves} segments d'extrusion")
    print(f"  filament total : {g.e:.1f} mm")
    print(f"  hauteur piece  : {PILLAR_HEIGHT + CIRCLE_DIAMETER:.1f} mm")


if __name__ == "__main__":
    main()
