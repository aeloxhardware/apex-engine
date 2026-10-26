"""Texte d'erreur de tétraédralisation, construit à partir du diagnostic chiffré de TetraNet.

Le mailleur renvoie des nombres ; le texte français vit ici (puis dans le moteur / l'interface),
pour qu'il puisse être traduit et reformulé sans toucher au C++. Chaque défaut s'accorde en
nombre : « un trou » / « des trous », « une face retournée » / « des faces retournées »…
"""

def plural(n):
    """Français : 0 et 1 prennent le singulier."""
    return n >= 2

def phrase(n, un, des):
    return un if not plural(n) else des

def defects(d, self_intersections=0):
    """Défauts constatés, du plus grave au moins grave, chacun accordé selon son nombre."""
    out = []
    holes = d.get("surface_holes", 0)
    if holes > 0: out.append(phrase(holes, "un trou", "des trous"))
    flipped = d.get("surface_flipped_edges", 0)
    if flipped > 0:  # trois arêtes incohérentes = une seule face isolée à l'envers
        out.append("une face retournée" if flipped <= 3 else "des faces retournées")
    if self_intersections > 0: out.append(phrase(self_intersections, "une auto-intersection", "des auto-intersections"))
    nm = d.get("surface_nonmanifold_edges", 0)
    if nm > 0: out.append(phrase(nm, "une arête partagée par plus de deux faces", "des arêtes partagées par plus de deux faces"))
    dup = d.get("surface_duplicate_faces", 0)
    if dup > 0: out.append(phrase(dup, "une face en double", "des faces en double"))
    return out

def join_fr(items):
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " et " + items[-1]

def message(d, volume_pct, self_intersections=0):
    found = defects(d, self_intersections) or ["un défaut de surface"]
    n = d["sign_conflicts"]
    pts = f"{n:,}".replace(",", " ") + (" points incohérents" if plural(n) else " point incohérent")
    return ("Erreur de tétraédralisation : la surface du modèle présente " + join_fr(found) +
            ", au-delà du seuil de tolérance, pour un total de " + pts + ", et un maillage volumique d'environ " +
            f"{volume_pct:.0f}" + " % du volume calculé sur la surface. Lancez les procédures de réparation automatique "
            "du maillage triangulé et assurez-vous que le fichier est complet et fermé.")

if __name__ == "__main__":
    import json, subprocess, sys
    cli, stl = sys.argv[1:3]
    r = json.loads(subprocess.run([cli, stl, "--cell", "1.0", "--alpha", "0.5"], capture_output=True, text=True).stdout.strip().splitlines()[-1])
    print({k: r[k] for k in ("status", "conflict_ratio", "warning", "surface_holes", "surface_open_edges", "surface_flipped_edges", "surface_nonmanifold_edges", "surface_duplicate_faces")})
    if r["status"] == 4: print(message(r, 100 * r["volume_ratio"]))
