"""Texte d'erreur de tétraédralisation (R-016, apex/ROADMAP.md), construit à partir du diagnostic
chiffré de TetraNet. Le mailleur renvoie des NOMBRES ; le texte français vit ici, pour pouvoir
être traduit et reformulé sans toucher au C++.

Règle de rédaction de l'utilisateur : chaque défaut cité s'accorde en nombre (« un trou » /
« des trous », « une face retournée » / « des faces retournées »…) ; « point incohérent » est
au singulier pour 0 et 1, au pluriel à partir de 2. Seuls les défauts réellement constatés
sont cités, du plus grave au moins grave.

Le même texte existe côté interface (apex-ui) : tests/test_surface_message.py le verrouille.
"""
from __future__ import annotations


def plural(n: int) -> bool:
    """Français : 0 et 1 prennent le singulier."""
    return n >= 2


def _phrase(n: int, singular: str, plural_form: str) -> str:
    return plural_form if plural(n) else singular


def defects(d: dict, self_intersections: int = 0) -> list[str]:
    out: list[str] = []
    holes = int(d.get("surface_holes", 0))
    if holes > 0:
        out.append(_phrase(holes, "un trou", "des trous"))
    flipped = int(d.get("surface_flipped_edges", 0))
    if flipped > 0:  # trois arêtes incohérentes = une seule face isolée à l'envers
        out.append("une face retournée" if flipped <= 3 else "des faces retournées")
    if self_intersections > 0:
        out.append(_phrase(self_intersections, "une auto-intersection", "des auto-intersections"))
    non_manifold = int(d.get("surface_nonmanifold_edges", 0))
    if non_manifold > 0:
        out.append(_phrase(non_manifold, "une arête partagée par plus de deux faces",
                           "des arêtes partagées par plus de deux faces"))
    duplicates = int(d.get("surface_duplicate_faces", 0))
    if duplicates > 0:
        out.append(_phrase(duplicates, "une face en double", "des faces en double"))
    return out


def _join(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " et " + items[-1]


def surface_error_message(d: dict, volume_pct: float, self_intersections: int = 0) -> str:
    """`d` : diagnostic de TetraNet (`ambiguous_points` = points dont le nombre d'enroulement ne
    tranche pas entre intérieur et extérieur), `volume_pct` : volume du maillage obtenu en % du
    volume calculé sur la surface."""
    found = defects(d, self_intersections) or ["un défaut de surface"]
    n = int(d.get("ambiguous_points", d.get("sign_conflicts", 0)))
    points = f"{n:,}".replace(",", " ") + (" points incohérents" if plural(n) else " point incohérent")
    return ("Erreur de tétraédralisation : la surface du modèle présente " + _join(found) +
            ", au-delà du seuil de tolérance, pour un total de " + points +
            ", et un maillage volumique d'environ " + f"{volume_pct:.0f}" +
            " % du volume calculé sur la surface. Lancez les procédures de réparation automatique "
            "du maillage triangulé et assurez-vous que le fichier est complet et fermé.")
