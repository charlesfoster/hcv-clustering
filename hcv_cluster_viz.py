#!/usr/bin/env python3
"""Cluster-network visualization, shared by the CLI (--plot-network) and the Streamlit GUI."""

from __future__ import annotations

import networkx as nx
import plotly.colors
import plotly.graph_objects as go


def filter_singletons(
    node_rows: list[dict[str, str]], edge_rows: list[dict[str, str]]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    kept_ids = {row["sample_id"] for row in node_rows if int(row["cluster_size"]) > 1}
    filtered_nodes = [row for row in node_rows if row["sample_id"] in kept_ids]
    filtered_edges = [row for row in edge_rows if row["source"] in kept_ids and row["target"] in kept_ids]
    return filtered_nodes, filtered_edges


def build_network_figure(node_rows: list[dict[str, str]], edge_rows: list[dict[str, str]]) -> go.Figure:
    graph = nx.Graph()
    for row in node_rows:
        graph.add_node(row["sample_id"], cluster_id=row["cluster_id"], cluster_size=int(row["cluster_size"]))
    for row in edge_rows:
        if row["source"] in graph and row["target"] in graph:
            graph.add_edge(row["source"], row["target"], distance=float(row["distance"]))

    layout = nx.spring_layout(graph, seed=42, k=1 / max(len(graph), 1) ** 0.5)

    edge_x: list[float | None] = []
    edge_y: list[float | None] = []
    for source, target in graph.edges():
        x0, y0 = layout[source]
        x1, y1 = layout[target]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, mode="lines", line=dict(width=1, color="#999999"), hoverinfo="none"
    )

    multi_cluster_ids = sorted(
        {data["cluster_id"] for _, data in graph.nodes(data=True) if data["cluster_size"] > 1}
    )
    palette = plotly.colors.qualitative.Plotly
    color_by_cluster = {cid: palette[i % len(palette)] for i, cid in enumerate(multi_cluster_ids)}

    node_x, node_y, node_color, node_size, node_text = [], [], [], [], []
    for sample_id, data in graph.nodes(data=True):
        x, y = layout[sample_id]
        node_x.append(x)
        node_y.append(y)
        singleton = data["cluster_size"] == 1
        node_color.append("#c7c7c7" if singleton else color_by_cluster[data["cluster_id"]])
        node_size.append(8 if singleton else 13)
        node_text.append(f"{sample_id}<br>cluster: {data['cluster_id']}<br>cluster size: {data['cluster_size']}")

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers",
        marker=dict(color=node_color, size=node_size, line=dict(width=1, color="#333333")),
        text=node_text,
        hoverinfo="text",
    )

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        showlegend=False,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=10, r=10, t=10, b=10),
        height=600,
        plot_bgcolor="rgba(0,0,0,0)",
    )
    return fig
