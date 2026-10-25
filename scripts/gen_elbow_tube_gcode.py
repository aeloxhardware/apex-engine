"""G-code AUTONOME : tube coude en surplomb, machine Apex 4 axes.

Ne depend de RIEN du projet (ni pipeline, ni core) : un seul fichier,
stdlib seule. Remplace scripts/gen_l_overhang_gcode.py, dont l'impression
du 1er septembre a rate -- voir "LE DEFAUT CORRIGE" ci-dessous.

PIECE : un tube ROND de diametre TUBE_DIAMETER.
  - branche verticale  : PILLAR_HEIGHT mm (axe du tube vertical)
  - coude              : l'axe du tube bascule de 90 deg (vertical ->
                         horizontal), la tete s'incline de 0 a 45 deg
  - branche horizontale: OVERHANG_LENGTH mm, en surplomb pur

SORTIE FIRMWARE : G00 Y<C> X<R> Z<Z> B<B> [E<e>] F<f>
  X = R, rayon plateau (mm) ; Y = C, rotation plateau (deg) ;
  Z = hauteur (mm) ; B = inclinaison tete (deg, 0 = buse vers le bas).

=============================================================================
LE DEFAUT CORRIGE (impression ratee du 1er septembre)
=============================================================================
Symptome rapporte : "une mini-bouteille dont le goulot se tordait et se
resserrait", plateau parti en butee. Trace dans le G-code de septembre :
derniere ligne Y111600, soit 310 TOURS cumules.

Cause : l'ancien script tracait la section du tube comme un cercle CENTRE
SUR L'AXE DE ROTATION (R=0). Tracer un cercle centre, c'est maintenir R
constant et faire balayer C sur 360 deg -- tres elegant tant que la piece
est un pilier vertical, et c'est CE QUI CASSE des que la tete s'incline :
la section, elle, cesse d'etre centree sur l'axe (elle part a l'horizontale
avec la branche), mais le script continuait de lui faire faire un tour
complet par couche en reduisant son rayon machine. D'ou le col qui se
resserre en vrille.

Correction : UN SEUL REGIME DE MOUVEMENT SUR TOUTE LA PIECE -- l'oscillation.
Jamais de rotation continue, nulle part, y compris sur le pilier et sur la
premiere couche ou elle serait pourtant tentante (et "naturelle" pour un
plateau tournant). Deux raisons, dans cet ordre :

  1. Un CHANGEMENT de regime en cours de piece est precisement ce qui a
     casse en septembre : la rotation continue tenait tant que la section
     restait centree, puis la geometrie cessait de le permettre sans que le
     mouvement change -- d'ou le goulot en vrille.
  2. La cinematique interdit l'oscillation sur une piece centree. La tete
     n'a que deux libertes de position (R et Z) : elle balaie un DEMI-PLAN
     vertical fixe, et R ne franchit jamais le centre. Un point (x,y) n'est
     atteignable que si le plateau l'amene dans ce plan, donc C = atan2(y,x)
     est IMPOSE par le point. Contour qui entoure l'axe => C parcourt
     360 deg, point. Aucune oscillation de R ou de B n'atteint un point
     situe derriere l'axe.

D'ou PART_OFFSET_R : on decale la piece pour que son contour n'entoure
jamais l'axe. C'est ce decalage, et lui seul, qui rend l'oscillation
possible.

=============================================================================
CONVENTION DE TETE (verifiee avec l'utilisateur, pas supposee)
=============================================================================
B > 0 penche la buse VERS LE CENTRE du plateau (vers -R). La matiere pousse
donc dans la direction opposee : le surplomb s'etend vers +R, vers
l'exterieur du plateau. Toute la geometrie ci-dessous en decoule.

=============================================================================
CE QUE CE SCRIPT NE PEUT PAS GARANTIR (a verifier avant de lancer)
=============================================================================
1. PIVOT DE LA TETE. Si l'axe de rotation B ne passe pas exactement par la
   pointe de la buse, incliner la tete DEPLACE la pointe, et il faut
   compenser. Regler PIVOT_ABOVE_TIP_MM (0 = le firmware compense deja, ou
   le pivot est a la pointe). Se mesure en 30 s : jog B de 0 a 90 deg et
   regarder de combien la pointe se deplace.
2. SEMANTIQUE DE F. Le firmware voit Y (degres) comme un axe lineaire : la
   vitesse reelle de la buse depend donc du rayon. Ce script calcule F par
   segment pour tenir une vitesse MATIERE constante (voir _feed_for).
   Si votre firmware ignore les axes rotatifs dans son calcul de vitesse,
   mettre CONSTANT_FEED = True pour retrouver le comportement de septembre.

Usage : python scripts/gen_elbow_tube_gcode.py [sortie.gcode] [--speed N] [--pivot L] [--full-infill] [--planar-cut [--cut-offset mm]]
  --speed N : multiplie TOUTES les vitesses par N (defaut 1)
  --full-infill : remplissage 100 % sur toutes les couches (anneaux concentriques)
  --planar-cut : plans horizontaux coupes a 45 deg puis couches inclinees (numpy, trimesh, shapely)
  --pivot L : distance pivot B -> pointe de buse, mm (defaut 0)
"""

from __future__ import annotations

import math
import os
import sys

# ---------------------------------------------------------------------------
# Parametres piece
# ---------------------------------------------------------------------------

TUBE_DIAMETER = 10.0        # diametre EXTERIEUR du tube, mm
PILLAR_HEIGHT = 15.0        # branche verticale (demande utilisateur), mm
OVERHANG_LENGTH = 30.0      # branche horizontale en surplomb (idem), mm
OVERHANG_ANGLE_DEG = 45.0   # inclinaison de tete atteinte progressivement

# Rayon de courbure de l'axe du tube dans le coude, mm. Fixe la douceur du
# raccord : trop petit, le tube se pince a l'interieur du virage ; trop
# grand, le coude mange toute la hauteur disponible. TUBE_DIAMETER donne un
# coude franc mais sain (rayon interieur = TUBE_DIAMETER/2, jamais nul).
BEND_RADIUS = TUBE_DIAMETER

# Decalage de TOUTE la piece par rapport a l'axe de rotation du plateau, mm.
#
# CE PARAMETRE N'EST PAS COSMETIQUE : il conditionne le regime de mouvement.
# La tete n'a que deux libertes de position (R et Z), donc elle ne balaie
# qu'un DEMI-PLAN vertical fixe ; R ne franchit jamais le centre. Un point
# de la piece a (x,y) n'est atteignable que si le plateau l'amene dans ce
# plan, c'est-a-dire a C = atan2(y, x) : C est IMPOSE par le point, jamais
# libre. Si le contour entoure l'axe, atan2 parcourt 360 deg et le plateau
# DOIT faire un tour complet par couche -- aucune oscillation de R ou de B
# ne peut aller chercher un point situe derriere l'axe.
#
# Mesure (voir l'historique de mise au point) : contour du tube D10, le
# decalage minimal qui supprime tout tour complet est ~9 mm.
#   decalage  0 mm -> C balaie 353 deg (tour complet, piece sur l'axe)
#   decalage  5 mm -> C balaie 177 deg (le contour frole encore l'axe)
#   decalage 10 mm -> C balaie  60 deg -> OSCILLATION sur toute la piece
#   decalage 15 mm -> C balaie  39 deg
PART_OFFSET_R = 15.0

# ---------------------------------------------------------------------------
# Parametres procede
# ---------------------------------------------------------------------------

LAYER_HEIGHT = 0.2
EXTRUSION_WIDTH = 0.4
FILAMENT_DIAMETER = 1.75
NOZZLE_TEMP = 210
# None = AUCUNE commande de plateau dans le fichier (ni M140, ni M190, ni
# extinction en fin) : le plateau n'est pas chauffe par ce G-code.
BED_TEMP = None

# False = la SEULE commande de temperature du fichier est M109 (chauffer
# la buse et attendre). Consequence : la buse RESTE CHAUDE apres
# l'impression, a eteindre a la main. True = ajoute M104 S0 en fin.
NOZZLE_OFF_AT_END = False

PRINT_SPEED = 25.0          # vitesse MATIERE visee, mm/s (vraie vitesse buse)
FIRST_LAYER_SPEED = 12.0    # idem, couches d'adherence
TRAVEL_SPEED = 60.0

SOLID_BASE_LAYER_COUNT = 3  # couches lentes a Z fixe (adherence plateau)

# Compensation du pivot de tete : distance (mm) entre l'axe de rotation B
# et la pointe de la buse. 0 = pivot a la pointe, ou firmware qui compense
# deja (cas le plus courant). Voir l'avertissement en tete de fichier.
PIVOT_ABOVE_TIP_MM = 0.0

# Repli si le firmware ne compte PAS les axes rotatifs dans sa vitesse.
CONSTANT_FEED = False
CONSTANT_FEED_VALUE = 1800.0   # mm/min, utilise seulement si CONSTANT_FEED

MAX_CHORD_MM = 0.6          # subdivision max entre deux points commandes

# Multiplicateur global de TOUTES les vitesses (impression, premiere couche,
# trajets a vide). Surcharge en ligne de commande : --speed 2
SPEED_FACTOR = 1.0
PLANAR_CUT = False          # --planar-cut : plans horizontaux coupes a 45 deg, puis couches inclinees
CUT_OFFSET = 0.0            # --cut-offset : decalage du plan de coupe le long de sa normale, mm
FULL_INFILL = False         # --full-infill : remplissage 100 % (anneaux concentriques)

_TUBE_RADIUS = TUBE_DIAMETER / 2.0
_FILAMENT_AREA = math.pi * (FILAMENT_DIAMETER / 2.0) ** 2


# ---------------------------------------------------------------------------
# Geometrie de la piece
# ---------------------------------------------------------------------------
#
# Repere piece : X = direction du surplomb (vers +R), Z = hauteur, Y
# lateral. L'axe de rotation du plateau est la droite X=Y=0.
#
# L'axe du tube ("spine") est parametre par une abscisse curviligne s :
#   - s in [0, PILLAR_HEIGHT]      : montee verticale sur l'axe machine
#   - coude                        : arc de cercle de rayon BEND_RADIUS,
#                                    centre a (BEND_RADIUS, 0, PILLAR_HEIGHT)
#   - puis                         : horizontale vers +X
#
# La tete, elle, s'incline de 0 a OVERHANG_ANGLE_DEG : l'axe du tube tourne
# de 90 deg, la tete seulement de 45. Les deux ne sont PAS lies -- c'est
# exactement ce qui fait du surplomb un surplomb.

_BEND_CENTER = (BEND_RADIUS, PILLAR_HEIGHT)       # (x, z), l'axe Y passe par la
_BEND_ARC = 0.5 * math.pi * BEND_RADIUS
_HORIZONTAL_Z = PILLAR_HEIGHT + BEND_RADIUS       # hauteur de l'axe horizontal
# L'axe est prolonge de _TUBE_RADIUS au-dela de la longueur demandee : les
# plans de coupe sont inclines a 45 deg, donc le DERNIER plan qui coupe
# encore le tube laisse une calotte de matiere au-dela de lui. Sans cette
# rallonge, la branche est tronquee de _TUBE_RADIUS (mesure : portee 40 mm
# au lieu des 45 attendus) -- la piece finit plus courte que commandee.
_SPINE_TOTAL = PILLAR_HEIGHT + _BEND_ARC + OVERHANG_LENGTH + _TUBE_RADIUS


def spine_at(s: float) -> tuple[float, float, float, float]:
    """Point de l'axe du tube et sa tangente, a l'abscisse curviligne s.

    Retourne (x, z, tx, tz) : position dans le plan (X,Z) et tangente
    unitaire. Y vaut toujours 0 (l'axe du tube reste dans le plan Y=0 ;
    c'est la SECTION qui s'etend en Y).
    """
    if s <= PILLAR_HEIGHT:
        return 0.0, s, 0.0, 1.0
    if s <= PILLAR_HEIGHT + _BEND_ARC:
        a = (s - PILLAR_HEIGHT) / BEND_RADIUS          # angle parcouru, rad
        cx, cz = _BEND_CENTER
        return (cx - BEND_RADIUS * math.cos(a),
                cz + BEND_RADIUS * math.sin(a),
                math.sin(a), math.cos(a))
    d = s - PILLAR_HEIGHT - _BEND_ARC
    return BEND_RADIUS + d, _HORIZONTAL_Z, 1.0, 0.0


def tube_point(s: float, phi: float,
               radius: float | None = None) -> tuple[float, float, float]:
    """Point de la SURFACE du tube : section circulaire a l'abscisse s.

    phi = 0 pointe vers l'exterieur du virage (dos du tube), phi = pi vers
    l'interieur. La section est un vrai cercle de rayon _TUBE_RADIUS dans le
    plan PERPENDICULAIRE a la tangente : c'est ce qui fait un tube ROND,
    et non le cylindre oblique ecrase que produisait l'ancien script.
    """
    rad = _TUBE_RADIUS if radius is None else radius
    x, z, tx, tz = spine_at(s)
    # normale dans le plan (X,Z), perpendiculaire a la tangente
    nx, nz = tz, -tx
    return (x + rad * math.cos(phi) * nx,
            rad * math.sin(phi),
            z + rad * math.cos(phi) * nz)


# ---------------------------------------------------------------------------
# Champ de decoupe : ce qui definit les COUCHES
# ---------------------------------------------------------------------------
#
# Une couche = une surface de niveau du champ height_field ci-dessous. Le
# champ encode l'inclinaison progressive de la tete :
#
#   - sous le coude   : plans horizontaux      (B = 0)
#   - dans le coude   : plans pivotant autour d'une droite fixe, la tete
#                       passant de 0 a 45 deg  (eventail de plans)
#   - au-dela         : plans PARALLELES a 45 deg, qui tranchent la branche
#                       horizontale en ellipses (B = 45 constant)
#
# L'eventail pivote autour de la droite Y a (x=BEND_RADIUS, z=PILLAR_HEIGHT),
# c'est-a-dire l'axe meme du coude : c'est la seule position ou les plans ne
# se croisent jamais a l'interieur de la matiere (sans quoi deux couches se
# chevaucheraient).

_FAN_PIVOT = _BEND_CENTER                      # (x, z)
_BETA_MAX = math.radians(OVERHANG_ANGLE_DEG)


def _fan_angle(x: float, z: float) -> tuple[float, float]:
    """Angle du point autour du pivot de l'eventail, et distance au pivot.

    Reference : 0 rad = rayon pivot -> DEBUT du coude (direction -X depuis
    le pivot), croissant vers +Z puis +X. Sur l'axe du tube, l'angle vaut
    donc 0 a l'entree du coude et +90 deg a sa sortie ; il est NEGATIF dans
    le pilier. Se tromper de reference ici renvoie B=0 en plein coude (et
    le tube repart en vrille au lieu de basculer) -- c'est exactement le
    bug corrige apres mesure sur la premiere version de ce script.
    """
    px, pz = _FAN_PIVOT
    dx, dz = x - px, z - pz
    if dx == 0.0 and dz == 0.0:
        return 0.0, 0.0
    return math.atan2(dz, -dx), math.hypot(dx, dz)


def height_field(x: float, z: float) -> float:
    """Hauteur de depot en un point : la couche k est le niveau k*LAYER_HEIGHT.

    Construit pour etre CONTINU et monotone le long du tube, condition pour
    que les couches s'empilent sans se croiser.
    """
    _px, pz = _FAN_PIVOT
    ang, rho = _fan_angle(x, z)

    if ang <= 0.0:
        return z                                   # plans horizontaux
    if ang >= _BETA_MAX:
        # plans paralleles inclines de _BETA_MAX, mesures depuis le plan
        # limite de l'eventail (continuite exacte en ang == _BETA_MAX)
        return pz + _BETA_MAX * BEND_RADIUS + rho * math.sin(ang - _BETA_MAX)
    # eventail : hauteur = arc parcouru, normalise au rayon du coude pour
    # que l'epaisseur de couche vaille LAYER_HEIGHT sur l'AXE du tube
    return pz + ang * BEND_RADIUS


def head_angle_deg(x: float, z: float) -> float:
    """Inclinaison B de la tete en un point : 0 sous le coude, 45 au-dela."""
    ang, _rho = _fan_angle(x, z)
    return math.degrees(min(max(ang, 0.0), _BETA_MAX))


def _bracket_s(phi: float, target: float,
               radius: float | None = None) -> float | None:
    """Trouve s tel que height_field(surface(s, phi)) == target (bissection).

    Retourne None si le niveau ne coupe pas cette generatrice. Le champ est
    monotone le long de chaque generatrice (verifie par construction : la
    normale des couches garde une composante positive sur la tangente du
    tube partout, tant que l'inclinaison reste <= 45 deg), donc une simple
    bissection suffit -- pas besoin de solveur.
    """
    lo, hi = 0.0, _SPINE_TOTAL

    def f(s: float) -> float:
        x, _y, z = tube_point(s, phi, radius)
        return height_field(x, z) - target

    flo, fhi = f(lo), f(hi)
    if flo > 0.0 or fhi < 0.0:
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Cinematique machine
# ---------------------------------------------------------------------------

def to_machine(x: float, y: float, z: float, b_deg: float,
               prev_c: float) -> tuple[float, float, float]:
    """(x, y, z) piece + inclinaison -> (R, C, Z) machine.

    C est ramene au plus proche de prev_c : sans ce recalage, un contour qui
    traverse +-180 deg ferait faire au plateau un tour complet en arriere.
    C'est aussi ce qui rend l'OSCILLATION lisible dans le fichier : C reste
    au voisinage de sa valeur precedente au lieu de sauter.
    """
    # Decalage applique ICI et nulle part ailleurs : toute la geometrie de
    # la piece est decrite dans un repere centre sur l'axe du tube, et ne
    # connait pas le plateau. Un seul endroit ou le decalage existe = aucun
    # risque qu'une partie du code l'applique et une autre l'oublie.
    x = x + PART_OFFSET_R

    r = math.hypot(x, y)
    c = prev_c if r < 1e-9 else math.degrees(math.atan2(y, x))
    c = _nearest_angle(c, prev_c)

    if PIVOT_ABOVE_TIP_MM:
        # Le firmware positionne le PIVOT, pas la pointe. Calibre a B=0, la
        # commande P coincide avec la pointe ; incline de B (vers -R), la
        # pointe arrive en :
        #     pointe = P + (-L sin B, +L (1 - cos B))      [en (R, Z)]
        # c'est-a-dire vers le centre ET plus haut -- symptome observe en
        # impression : chaque couche s'ecarte du tube a l'oppose du bras.
        # Pour poser la pointe en T, il faut donc commander
        #     P = T + (+L sin B, -L (1 - cos B)).
        # (Premiere version : signe de Z inverse, qui doublait la montee au
        # lieu de l'annuler. Verifie par simulation directe, voir historique.)
        b = math.radians(b_deg)
        r += PIVOT_ABOVE_TIP_MM * math.sin(b)
        z -= PIVOT_ABOVE_TIP_MM * (1.0 - math.cos(b))
    return r, c, z


def _nearest_angle(target: float, ref: float) -> float:
    while target - ref > 180.0:
        target -= 360.0
    while target - ref < -180.0:
        target += 360.0
    return target


def _feed_for(true_len: float, dr: float, dc: float, dz: float,
              speed_mm_s: float) -> float:
    """F (mm/min) pour que la BUSE avance a speed_mm_s dans la matiere.

    Le firmware voit Y en degres comme un axe lineaire : sa notion de
    distance est sqrt(dR^2 + dC^2 + dZ^2), qui n'a rien a voir avec la
    distance reelle parcourue par la buse des que le plateau tourne. On
    corrige donc F par le rapport des deux. Sans ca, la vitesse reelle
    depend du rayon -- c'est ce qui rendait les 300 mm/s de l'ancien script
    imprimables par accident au rayon ou il travaillait, et seulement la.
    """
    # SPEED_FACTOR applique ICI et nulle part ailleurs : c'est le seul point
    # par lequel passe chaque F du fichier (impression, premiere couche,
    # trajets a vide, mode CONSTANT_FEED). Un facteur applique ailleurs
    # laisserait forcement une vitesse de cote.
    if CONSTANT_FEED:
        return CONSTANT_FEED_VALUE * SPEED_FACTOR
    speed_mm_s *= SPEED_FACTOR
    firmware_len = math.sqrt(dr * dr + dc * dc + dz * dz)
    if true_len < 1e-9:
        return speed_mm_s * 60.0
    return speed_mm_s * 60.0 * (firmware_len / true_len)


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------

class GCode:
    def __init__(self) -> None:
        self.lines: list[str] = []
        self.e = 0.0
        self.c = 0.0
        self.pos: tuple[float, float, float] | None = None   # piece (x,y,z)
        self.c_min = 0.0
        self.c_max = 0.0
        self.trace = None      # liste optionnelle (x, y, z, extrude, couche), pour les apercus
        self.layer = 0

    def raw(self, s: str = "") -> None:
        self.lines.append(s)

    def comment(self, s: str) -> None:
        self.lines.append(f"; {s}")

    def go(self, x: float, y: float, z: float, b: float, speed: float,
           extrude: bool, layer_thickness: float = LAYER_HEIGHT) -> None:
        """Un point de trajectoire, exprime dans le repere PIECE."""
        r, c, mz = to_machine(x, y, z, b, self.c)
        true_len = 0.0
        if self.pos is not None:
            px, py, pz = self.pos
            true_len = math.dist((x, y, z), (px, py, pz))

        dr = dc = dz = 0.0
        if self.lines and self.pos is not None:
            dr, dc, dz = abs(r - self._r), abs(c - self.c), abs(mz - self._z)

        parts = ["G00", f"Y{c:.4f}", f"X{r:.4f}", f"Z{mz:.4f}", f"B{b:.4f}"]
        if extrude and true_len > 0.0:
            # Section deposee = epaisseur LOCALE de couche x largeur de
            # cordon. L'epaisseur varie dans le coude (les plans sont en
            # eventail : plus epais au dos du tube qu'au ventre), ce que
            # l'ancien script ignorait -> sous-extrusion au dos.
            self.e += true_len * layer_thickness * EXTRUSION_WIDTH / _FILAMENT_AREA
            parts.append(f"E{self.e:.5f}")
        parts.append(f"F{_feed_for(true_len, dr, dc, dz, speed):.0f}")

        self.lines.append(" ".join(parts))
        self.c, self._r, self._z = c, r, mz
        self.pos = (x, y, z)
        self.c_min, self.c_max = min(self.c_min, c), max(self.c_max, c)
        if self.trace is not None:
            self.trace.append((x, y, z, extrude and true_len > 0.0, self.layer))

    def dump(self) -> str:
        return "\n".join(self.lines) + "\n"


def _local_thickness(x: float, z: float) -> float:
    """Epaisseur reelle de couche en un point (eventail du coude).

    Dans l'eventail, deux plans consecutifs s'ecartent proportionnellement a
    la distance au pivot : LAYER_HEIGHT sur l'axe du tube, moins au ventre
    du virage, plus au dos. Hors eventail, epaisseur constante.
    """
    ang, rho = _fan_angle(x, z)
    if ang <= 0.0 or ang >= _BETA_MAX:
        return LAYER_HEIGHT
    return LAYER_HEIGHT * rho / BEND_RADIUS


def _section_points(level: float, radius: float | None = None
                    ) -> list[tuple[float, float, float]]:
    """Contour ferme d'une couche : intersection du niveau avec le tube.

    Pour chaque generatrice phi du tube, on cherche l'abscisse s ou le champ
    atteint le niveau demande. C'est ce qui donne de VRAIES ellipses dans la
    branche horizontale (plan a 45 deg coupant un tube rond), la ou l'ancien
    script dessinait un cercle dans le plan incline -- donc un tube ecrase.
    """
    n = max(48, int(math.ceil(2.0 * math.pi * _TUBE_RADIUS / MAX_CHORD_MM)))
    pts: list[tuple[float, float, float]] = []
    for i in range(n + 1):                      # +1 : referme le contour
        phi = 2.0 * math.pi * (i % n) / n
        s = _bracket_s(phi, level, radius)
        if s is None:
            return []
        pts.append(tube_point(s, phi, radius))
    return pts


def _emit_section(g: GCode, pts: list[tuple[float, float, float]],
                  speed: float, first_point_travel: bool) -> None:
    """Trace un contour ferme, en subdivisant pour rester fidele en (R,C)."""
    if first_point_travel:
        x, y, z = pts[0]
        b = head_angle_deg(x, z)
        g.go(x, y, z + 1.0, b, TRAVEL_SPEED, extrude=False)
        g.go(x, y, z, b, TRAVEL_SPEED, extrude=False)

    for i in range(1, len(pts)):
        x0, y0, z0 = pts[i - 1]
        x1, y1, z1 = pts[i]
        seg = math.dist((x0, y0, z0), (x1, y1, z1))
        steps = max(1, int(math.ceil(seg / MAX_CHORD_MM)))
        for k in range(1, steps + 1):
            t = k / steps
            x, y, z = (x0 + (x1 - x0) * t, y0 + (y1 - y0) * t,
                       z0 + (z1 - z0) * t)
            g.go(x, y, z, head_angle_deg(x, z), speed, extrude=True,
                 layer_thickness=_local_thickness(x, z))


def _emit_base_disc(g: GCode, level: float) -> None:
    """Premiere couche : disque plein, en VA-ET-VIENT (jamais une spirale).

    Une spirale d'Archimede centree serait le mouvement le plus naturel
    d'un plateau tournant -- et c'est precisement ce qu'on s'interdit ici :
    elle impose une rotation CONTINUE (12 tours mesures pour ce seul
    disque), alors que toute la piece doit rester en oscillation. Un
    balayage en va-et-vient parallele a Y garde C dans la meme plage etroite
    que les couches suivantes.

    ORDRE : contour D'ABORD, remplissage ensuite -- exception assumee a la
    regle "interieur vers exterieur" de CLAUDE.md, limitee a cette couche
    (la seule de la piece qui ait un remplissage). Pour une premiere couche :
      - le contour, pose sur un plateau propre, fixe le diametre du tube
        dont toute la paroi repose ensuite sur lui ;
      - chaque ligne de remplissage vient buter contre un cordon deja
        colle, au lieu d'avoir des extremites libres qui se relevent.
    """
    g.layer += 1
    g.comment("couche 1 : contour puis disque plein (va-et-vient), B=0")
    pts = [(_TUBE_RADIUS * math.cos(2.0 * math.pi * i / 96),
            _TUBE_RADIUS * math.sin(2.0 * math.pi * i / 96), level)
           for i in range(97)]
    _emit_section(g, pts, FIRST_LAYER_SPEED, first_point_travel=True)

    n = max(2, int(math.ceil(2.0 * _TUBE_RADIUS / EXTRUSION_WIDTH)))
    first = True
    for i in range(n + 1):
        x = -_TUBE_RADIUS + 2.0 * _TUBE_RADIUS * i / n
        half = math.sqrt(max(0.0, _TUBE_RADIUS ** 2 - x * x))
        if half < 1e-6:
            continue
        y0, y1 = (-half, half) if i % 2 == 0 else (half, -half)
        # rejoint le debut de la ligne a vide, puis la trace en extrudant
        g.go(x, y0, level, 0.0, TRAVEL_SPEED, extrude=False)
        _emit_span(g, x, y0, x, y1, level, FIRST_LAYER_SPEED, first)
        first = False


def _emit_span(g: GCode, x0: float, y0: float, x1: float, y1: float,
               z: float, speed: float, _first: bool) -> None:
    """Un segment droit, subdivise pour rester fidele une fois en (R, C)."""
    seg = math.hypot(x1 - x0, y1 - y0)
    steps = max(1, int(math.ceil(seg / MAX_CHORD_MM)))
    for k in range(1, steps + 1):
        t = k / steps
        g.go(x0 + (x1 - x0) * t, y0 + (y1 - y0) * t, z, 0.0, speed,
             extrude=True)



# ---------------------------------------------------------------------------
# Mode --planar-cut : couches planes COUPEES a 45 deg, puis couches inclinees
# ---------------------------------------------------------------------------
#
# Au lieu de l'eventail continu, les plans horizontaux sont tronques
# progressivement par un plan a OVERHANG_ANGLE_DEG (chaque couche plane est
# plus courte que la precedente de LAYER_HEIGHT / tan(angle) : un escalier a
# 45 deg, imprimable sans support). Les couches inclinees, PARALLELES a ce plan,
# commencent ensuite posees dessus. Le plan de coupe passe par le pivot de
# l'ancien eventail (centre du coude) : il coupe alors le tube selon sa section
# circulaire a mi-coude, ou les deux parois sont a 45 deg de la verticale.
#
# Contrairement au mode historique (contours generes generatrice par
# generatrice), ce mode tranche un VRAI maillage du tube : necessite numpy,
# trimesh et shapely, importes seulement ici.

def _planar_cut_run(g: "GCode") -> None:
    import numpy as np
    import trimesh
    from shapely.geometry import Polygon, box
    from shapely.ops import polygonize, unary_union
    from shapely.geometry import LineString

    beta = _BETA_MAX
    nrm = np.array([math.sin(beta), 0.0, math.cos(beta)])
    u_inc = np.array([math.cos(beta), 0.0, -math.sin(beta)])
    w = np.array([0.0, 1.0, 0.0])
    c0 = _FAN_PIVOT[0] * nrm[0] + _FAN_PIVOT[1] * nrm[2] + CUT_OFFSET

    # maillage ferme du tube
    ds = 0.2
    ns = int(math.ceil(_SPINE_TOTAL / ds))
    nphi = 128
    verts = []
    for si in range(ns + 1):
        sv = _SPINE_TOTAL * si / ns
        for pi in range(nphi):
            verts.append(tube_point(sv, 2.0 * math.pi * pi / nphi))
    faces = []
    for si in range(ns):
        for pi in range(nphi):
            a = si * nphi + pi
            b = si * nphi + (pi + 1) % nphi
            c = (si + 1) * nphi + pi
            d = (si + 1) * nphi + (pi + 1) % nphi
            faces += [(a, b, c), (b, d, c)]
    for si, sv in ((0, 0.0), (ns, _SPINE_TOTAL)):
        x, z, _tx, _tz = spine_at(sv)
        verts.append((x, 0.0, z))
        ci = len(verts) - 1
        for pi in range(nphi):
            faces.append((ci, si * nphi + pi, si * nphi + (pi + 1) % nphi))
    mesh = trimesh.Trimesh(np.array(verts), np.array(faces), process=False)
    v_all = mesh.vertices @ nrm
    v_max = float(v_all.max())

    def slice_polys(normal, origin, u):
        seg = trimesh.intersections.mesh_plane(mesh, normal, origin)
        if len(seg) == 0:
            return None
        lines = []
        for a, b in seg:
            pa = (round(float((a - origin) @ u), 6), round(float(a @ w), 6))
            pb = (round(float((b - origin) @ u), 6), round(float(b @ w), 6))
            if pa != pb:
                lines.append(LineString([pa, pb]))
        polys = list(polygonize(unary_union(lines)))
        if not polys:
            return None
        return unary_union(polys)

    def rings_of(geom):
        """Anneaux concentriques (pas = largeur de cordon), du centre au bord."""
        out = []
        pieces = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
        for poly in pieces:
            if poly.is_empty or poly.area < EXTRUSION_WIDTH ** 2:
                continue
            levels = []
            k = 0
            while True:
                cur = poly.buffer(-k * EXTRUSION_WIDTH) if k else poly
                if cur.is_empty:
                    break
                levels.append((k * EXTRUSION_WIDTH, cur))
                k += 1
            if not FULL_INFILL:
                levels = levels[:1]
            else:
                kd = k * EXTRUSION_WIDTH - EXTRUSION_WIDTH * 0.5
                core = poly.buffer(-kd)
                if not core.is_empty:
                    levels.append((kd, core))
            for _d, cur in reversed(sorted(levels, key=lambda t: t[0])):
                for gg in (cur.geoms if hasattr(cur, "geoms") else [cur]):
                    for ring in [gg.exterior, *gg.interiors]:
                        coords = list(ring.simplify(0.01).coords)
                        if len(coords) >= 4:
                            out.append(coords)
        return out

    def to3(coords, origin, u):
        return [tuple(origin + a * u + bb * w) for a, bb in coords]

    def emit_layer(coords_list, origin, u, b_deg, speed, lift):
        g.layer += 1
        first = True
        for coords in coords_list:
            pts = to3(coords, origin, u)
            if g.pos is not None and not (first and lift):
                # depart le plus proche de la position courante
                px, py, pz = g.pos
                i0 = min(range(len(pts) - 1),
                         key=lambda i: math.dist(pts[i], (px, py, pz)))
                pts = pts[i0:-1] + pts[:i0] + [pts[i0]]
            x, y, z = pts[0]
            if first and lift:
                g.go(x, y, z + 2.0, b_deg, TRAVEL_SPEED, extrude=False)
                g.go(x, y, z, b_deg, TRAVEL_SPEED, extrude=False)
            else:
                g.go(x, y, z, b_deg, TRAVEL_SPEED, extrude=False)
            first = False
            for i in range(1, len(pts)):
                x0, y0, z0 = pts[i - 1]
                x1, y1, z1 = pts[i]
                seg = math.dist(pts[i - 1], pts[i])
                steps = max(1, int(math.ceil(seg / MAX_CHORD_MM)))
                for kk in range(1, steps + 1):
                    t = kk / steps
                    g.go(x0 + (x1 - x0) * t, y0 + (y1 - y0) * t,
                         z0 + (z1 - z0) * t, b_deg, speed, extrude=True)

    e_z = np.array([0.0, 0.0, 1.0])
    u_pl = np.array([1.0, 0.0, 0.0])

    g.raw()
    g.comment("--- couches planes (B=0) tronquees par le plan a "
              f"{OVERHANG_ANGLE_DEG:g} deg ---")
    k = 2
    n_planar = 0
    while True:
        z = k * LAYER_HEIGHT
        geom = slice_polys(e_z, np.array([0.0, 0.0, z]), u_pl)
        if geom is None:
            break
        xc = (c0 - z * nrm[2]) / nrm[0]
        geom = geom.intersection(box(-1e3, -1e3, xc, 1e3))
        if geom.is_empty or geom.area < EXTRUSION_WIDTH ** 2:
            break
        speed = FIRST_LAYER_SPEED if k <= SOLID_BASE_LAYER_COUNT else PRINT_SPEED
        emit_layer(rings_of(geom), np.array([0.0, 0.0, z]), u_pl, 0.0, speed,
                   lift=(k == 2))
        n_planar += 1
        k += 1

    g.raw()
    g.comment(f"--- couches inclinees (B={OVERHANG_ANGLE_DEG:g}), paralleles "
              "au plan de coupe ---")
    j = 1
    n_incl = 0
    while True:
        c = c0 + j * LAYER_HEIGHT
        if c > v_max - 1e-6:
            break
        origin = nrm * c
        geom = slice_polys(nrm, origin, u_inc)
        if geom is not None and not geom.is_empty and geom.area >= EXTRUSION_WIDTH ** 2:
            emit_layer(rings_of(geom), origin, u_inc, OVERHANG_ANGLE_DEG,
                       PRINT_SPEED, lift=(n_incl == 0))
            n_incl += 1
        j += 1
    g.comment(f"planar-cut : {n_planar + 1} couches planes, {n_incl} inclinees, "
              f"plan de coupe c0={c0:.3f}")


def build() -> GCode:
    g = GCode()
    total_layers = int(math.floor(_max_field_level() / LAYER_HEIGHT))

    g.comment("=== tube coude en surplomb - gen_elbow_tube_gcode.py ===")
    g.comment("G00 Y<C> X<R> Z<Z> B<B> F  |  B>0 = buse penchee vers -R")
    g.comment(f"tube rond D{TUBE_DIAMETER} mm, pilier {PILLAR_HEIGHT} mm, "
              f"surplomb {OVERHANG_LENGTH} mm a {OVERHANG_ANGLE_DEG} deg")
    g.comment(f"coude R{BEND_RADIUS} mm, couche {LAYER_HEIGHT} mm, "
              f"cordon {EXTRUSION_WIDTH} mm")
    g.comment(f"pivot B -> pointe : {PIVOT_ABOVE_TIP_MM:g} mm"
              + ("  (AUCUNE compensation)" if not PIVOT_ABOVE_TIP_MM else ""))
    g.comment(f"vitesse x{SPEED_FACTOR:g} : impression "
              f"{PRINT_SPEED * SPEED_FACTOR:g} mm/s, couche 1 "
              f"{FIRST_LAYER_SPEED * SPEED_FACTOR:g} mm/s, vide "
              f"{TRAVEL_SPEED * SPEED_FACTOR:g} mm/s")
    g.raw()
    g.raw("G21")
    g.raw("G90")
    g.raw("M82")
    g.raw("G92 E0")
    if BED_TEMP is not None:
        g.raw(f"M140 S{BED_TEMP}")
    g.raw(f"M109 S{NOZZLE_TEMP}")
    if BED_TEMP is not None:
        g.raw(f"M190 S{BED_TEMP}")
    g.raw("G00 B0.0000 F1200")
    g.raw("G00 Y0.0000 X0.0000 Z5.0000 F1200")
    g.raw()

    phase = None
    if PLANAR_CUT:
        _emit_base_disc(g, LAYER_HEIGHT)
        _planar_cut_run(g)
    for k in (range(1, total_layers + 1) if not PLANAR_CUT else ()):
        level = k * LAYER_HEIGHT
        if k == 1:
            _emit_base_disc(g, level)
            continue

        pts = _section_points(level)
        if not pts:
            break

        b_here = max(head_angle_deg(x, z) for x, _y, z in pts)
        new_phase = ("pilier" if b_here <= 1e-6 else
                     "coude" if b_here < OVERHANG_ANGLE_DEG - 1e-6 else
                     "surplomb")
        if new_phase != phase:
            phase = new_phase
            g.raw()
            if phase == "pilier":
                g.comment("--- pilier : B=0, R et C oscillent ---")
            elif phase == "coude":
                g.comment("--- coude : B monte vers 45 deg, R et C oscillent ---")
            else:
                g.comment("--- surplomb : B=45 constant, R et C oscillent ---")

        speed = FIRST_LAYER_SPEED if k <= SOLID_BASE_LAYER_COUNT else PRINT_SPEED
        if FULL_INFILL:
            # remplissage 100 % : anneaux concentriques, du centre vers la
            # paroi (regle "interieur vers exterieur"), pas = largeur de cordon.
            rings = []
            r = _TUBE_RADIUS - EXTRUSION_WIDTH
            while r >= EXTRUSION_WIDTH / 2.0 - 1e-9:
                inner = _section_points(level, r)
                if inner:
                    rings.append(inner)
                r -= EXTRUSION_WIDTH
            rings.reverse()
            first = True
            for ring in rings:
                x, y, z = ring[0]
                if first:
                    first = False
                    _emit_section(g, ring, speed, first_point_travel=(k == 2))
                else:
                    g.go(x, y, z, head_angle_deg(x, z), TRAVEL_SPEED,
                         extrude=False)
                    _emit_section(g, ring, speed, first_point_travel=False)
            if rings:
                x, y, z = pts[0]
                g.go(x, y, z, head_angle_deg(x, z), TRAVEL_SPEED, extrude=False)
            _emit_section(g, pts, speed,
                          first_point_travel=(k == 2 and not rings))
        else:
            _emit_section(g, pts, speed, first_point_travel=(k == 2))

    g.raw()
    g.comment("fin")
    # Degagement : remonter D'ABORD, redresser la tete ENSUITE. Avec un pivot
    # au-dessus de la pointe, redresser B en place fait decrire a la pointe
    # un arc de L*sin(B) vers l'exterieur et L*(1-cos B) vers le bas (32 mm
    # et 13 mm pour L=45), au ras de la piece juste imprimee. Monter tout
    # droit en gardant B eloigne la pointe du cote non encore imprime (la
    # derniere couche est la plus haute du champ), sans rien croiser.
    g.raw(f"G00 Z{_HORIZONTAL_Z + 30.0:.4f} F1200")
    g.raw("G00 B0.0000 F1200")
    if NOZZLE_OFF_AT_END:
        g.raw("M104 S0")
    if BED_TEMP is not None:
        g.raw("M140 S0")
    g.raw("M84")
    return g


def _max_field_level() -> float:
    """Niveau de champ le plus eleve atteint par la piece.

    C'est le NIVEAU (au sens de height_field), pas la hauteur Z : au bout de
    la branche horizontale le tube n'est qu'a z=25 alors que le champ y vaut
    ~51, les plans etant inclines a 45 deg. Confondre les deux tronque la
    piece avant meme le coude -- mesure sur la premiere version : 125
    couches au lieu de 256, aucun surplomb genere.

    Balaye la surface entiere (pas seulement l'axe) : le maximum est atteint
    sur une generatrice, jamais au centre.
    """
    best = 0.0
    for i in range(400):
        s = _SPINE_TOTAL * i / 399.0
        for j in range(72):
            x, _y, z = tube_point(s, 2.0 * math.pi * j / 72.0)
            best = max(best, height_field(x, z))
    return best


def main() -> None:
    global SPEED_FACTOR, FULL_INFILL, PLANAR_CUT, CUT_OFFSET
    args = sys.argv[1:]
    if "--planar-cut" in args:
        PLANAR_CUT = True
        args.remove("--planar-cut")
    if "--cut-offset" in args:
        i = args.index("--cut-offset")
        CUT_OFFSET = float(args[i + 1])
        del args[i:i + 2]
    if "--full-infill" in args:
        FULL_INFILL = True
        args.remove("--full-infill")
    if "--speed" in args:
        i = args.index("--speed")
        SPEED_FACTOR = float(args[i + 1])
        if SPEED_FACTOR <= 0.0:
            raise SystemExit("--speed doit etre > 0")
        del args[i:i + 2]
    if "--pivot" in args:
        global PIVOT_ABOVE_TIP_MM
        i = args.index("--pivot")
        PIVOT_ABOVE_TIP_MM = float(args[i + 1])
        del args[i:i + 2]
    out = args[0] if args else os.path.join("out", "elbow_tube.gcode")
    if os.path.dirname(out):
        os.makedirs(os.path.dirname(out), exist_ok=True)

    g = build()
    with open(out, "w", encoding="ascii") as f:
        f.write(g.dump())

    n_extrude = sum(1 for ln in g.lines if ln.startswith("G00") and " E" in ln)
    print(f"ecrit : {out}")
    print(f"  {len(g.lines)} lignes, {n_extrude} segments d'extrusion")
    print(f"  filament : {g.e:.1f} mm")
    print(f"  hauteur piece : {_HORIZONTAL_Z + _TUBE_RADIUS:.1f} mm, "
          f"portee {BEND_RADIUS + OVERHANG_LENGTH + _TUBE_RADIUS:.1f} mm")
    print(f"  plage de C : {g.c_min:.1f} a {g.c_max:.1f} deg "
          f"(amplitude {g.c_max - g.c_min:.1f})")


if __name__ == "__main__":
    main()
