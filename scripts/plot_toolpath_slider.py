"""Trajectoires par couche (ordre d'impression) avec curseur sur la derniere couche affichee.
Usage : python scripts/plot_toolpath_slider.py <modele.stl> <sortie.html> <variante> [resolution]
variante : ref | heat (champ = distance au plateau par la methode de la chaleur)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
from main import PipelineParams
from scripts.export_slicing_preview import run_slicing_only
from pipeline.step6_toolpath2d import _lift_to_3d


def main(stl, out, variant, resolution):
    p = PipelineParams(resolution=resolution)
    from scripts.measure_matrix import apply_variant
    apply_variant(p, variant)
    r = run_slicing_only(stl, p)
    layers = [L for L in r.colored_layers if L.toolpath_2d]
    n = len(layers)
    fig = go.Figure()
    md = r.mesh_data
    fig.add_trace(go.Mesh3d(x=md.vertices[:, 0], y=md.vertices[:, 1], z=md.vertices[:, 2], i=md.faces[:, 0], j=md.faces[:, 1],
                            k=md.faces[:, 2], color="lightgrey", opacity=0.07, hoverinfo="skip", name="piece"))
    cols = sample_colorscale("Turbo", list(np.linspace(0, 1, max(n, 2))))
    for k, L in enumerate(layers):
        xs, ys, zs = [], [], []
        for pl in L.toolpath_2d:
            P = _lift_to_3d(np.asarray(pl), L.local_frame)
            xs += list(P[:, 0]) + [None]; ys += list(P[:, 1]) + [None]; zs += list(P[:, 2]) + [None]
        fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", line=dict(color=cols[k], width=2),
                                   name=f"couche {k+1}", showlegend=False, hoverinfo="name", visible=True))
    steps = []
    for m in range(n + 1):
        vis = [True] + [kk < m for kk in range(n)]
        steps.append(dict(method="update", label=str(m), args=[{"visible": vis}, {"title": f"{os.path.basename(stl)} [{variant}] : couches 1 a {m} sur {n} (ordre d'impression)"}]))
    fig.update_layout(template="plotly_dark", scene=dict(aspectmode="data"),
                      sliders=[dict(active=n, currentvalue=dict(prefix="derniere couche affichee : "), pad=dict(t=30), steps=steps)],
                      title=f"{os.path.basename(stl)} [{variant}] : couches 1 a {n} sur {n} (ordre d'impression)")
    fig.write_html(out, include_plotlyjs=True)
    print("ecrit", out, "couches", n)


if __name__ == "__main__":
    a = sys.argv
    main(a[1], a[2], a[3], float(a[4]) if len(a) > 4 else 0.8)
