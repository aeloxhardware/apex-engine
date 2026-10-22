"""Graphe 3D interactif (plotly) de la distance de Dijkstra au plateau dans le maillage tetraedrique.

Montre : la piece (translucide), les centres de tetraedres colores par distance geodesique au plateau,
et les plus courts chemins (arbre de Dijkstra) depuis les noeuds les plus eloignes.
Usage : python scripts/plot_dijkstra_interactive.py <modele.stl> <sortie.html> [resolution] [nb_chemins]
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra
import plotly.graph_objects as go
from main import PipelineParams
from pipeline import step4_field_solver as s4
from scripts.export_slicing_preview import run_slicing_only


def build_graph(tm):
    t = tm.tetrahedra
    pairs = np.vstack([t[:, [0, 1]], t[:, [0, 2]], t[:, [0, 3]], t[:, [1, 2]], t[:, [1, 3]], t[:, [2, 3]]])
    w = np.linalg.norm(tm.tet_vertices[pairs[:, 0]] - tm.tet_vertices[pairs[:, 1]], axis=1)
    n = tm.node_count
    return sp.coo_matrix((np.r_[w, w], (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])), shape=(n, n)).tocsr()


def main(stl, out, resolution=1.5, n_paths=6):
    r = run_slicing_only(stl, PipelineParams(resolution=resolution))
    tm, md, rg = r.tet_mesh, r.mesh_data, r.region_dependency_graph
    plat = s4._identify_platform_nodes(tm, md, rg)
    V = tm.tet_vertices
    dist, pred, src = dijkstra(build_graph(tm), directed=False, indices=plat, min_only=True, return_predecessors=True)
    finite = np.isfinite(dist)
    # chemins : noeuds les plus eloignes, espaces d'au moins 15 % de la diagonale
    diag = np.linalg.norm(V.max(0) - V.min(0))
    chosen = []
    for i in np.argsort(-np.where(finite, dist, -1)):
        if len(chosen) >= n_paths:
            break
        if all(np.linalg.norm(V[i] - V[j]) > 0.15 * diag for j in chosen):
            chosen.append(i)
    fig = go.Figure()
    fig.add_trace(go.Mesh3d(x=md.vertices[:, 0], y=md.vertices[:, 1], z=md.vertices[:, 2],
                            i=md.faces[:, 0], j=md.faces[:, 1], k=md.faces[:, 2],
                            color="lightgrey", opacity=0.12, name="piece", hoverinfo="skip"))
    c = V[tm.tetrahedra].mean(axis=1)
    dt = dist[tm.tetrahedra].mean(axis=1)
    keep = np.isfinite(dt)
    rng = np.random.default_rng(0)
    idx = np.flatnonzero(keep)
    if len(idx) > 40000:
        idx = rng.choice(idx, 40000, replace=False)
    fig.add_trace(go.Scatter3d(x=c[idx, 0], y=c[idx, 1], z=c[idx, 2], mode="markers",
                               marker=dict(size=2, color=dt[idx], colorscale="Viridis", opacity=0.55,
                                           colorbar=dict(title="distance<br>Dijkstra (mm)")),
                               name="tetraedres", text=[f"{d:.1f} mm" for d in dt[idx]], hoverinfo="text+x+y+z"))
    fig.add_trace(go.Scatter3d(x=V[plat, 0], y=V[plat, 1], z=V[plat, 2], mode="markers",
                               marker=dict(size=2, color="red"), name="noeuds plateau (sources)", visible="legendonly"))
    for k, i in enumerate(chosen):
        path = [i]
        while pred[path[-1]] >= 0:
            path.append(pred[path[-1]])
        P = V[path]
        fig.add_trace(go.Scatter3d(x=P[:, 0], y=P[:, 1], z=P[:, 2], mode="lines+markers",
                                   line=dict(color="white", width=6), marker=dict(size=3, color="white"),
                                   name=f"chemin {k+1} ({dist[i]:.1f} mm)", legendgroup=f"c{k}"))
        # fleches dans le sens d'impression : du plateau vers le noeud eloigne (inverse du parcours de Dijkstra)
        if len(path) > 1:
            P = V[path]
            d = np.empty_like(P)
            d[:-1] = P[:-1] - P[1:]   # sens plateau -> noeud eloigne = inverse du parcours de Dijkstra
            d[-1] = d[-2]             # dernier sommet (plateau) : meme direction que le segment voisin
            d /= np.linalg.norm(d, axis=1, keepdims=True) + 1e-12
            fig.add_trace(go.Cone(x=P[:, 0], y=P[:, 1], z=P[:, 2], u=d[:, 0], v=d[:, 1], w=d[:, 2],
                                  sizemode="absolute", sizeref=0.024 * diag, anchor="cm", showscale=False,
                                  colorscale=[[0, "red"], [1, "red"]], name=f"sens d'impression {k+1}",
                                  legendgroup=f"c{k}", showlegend=False, hoverinfo="skip"))
    fig.update_layout(template="plotly_dark", scene=dict(aspectmode="data"), title=os.path.basename(stl) +
                      " : distance de Dijkstra au plateau et plus courts chemins", legend=dict(itemsizing="constant"))
    fig.write_html(out, include_plotlyjs=True)
    print("ecrit", out, "tetraedres", len(tm.tetrahedra), "chemins", len(chosen))


if __name__ == "__main__":
    a = sys.argv
    main(a[1], a[2], float(a[3]) if len(a) > 3 else 1.5, int(a[4]) if len(a) > 4 else 6)
