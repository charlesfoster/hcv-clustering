#!/usr/bin/env python3
"""Network visualization shared by the CLI and Streamlit GUI.

The clustering tables remain the source of graph topology. Optional metadata only
changes presentation and hover text.
"""

from __future__ import annotations

import hashlib
import html
import math
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping, Sequence

import networkx as nx
import plotly.colors
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from hcv_cluster_metadata import AGE_RANGE_ORDER, MISSING_VALUE

PositionMap = dict[str, tuple[float, float]]
LayoutMode = Literal["spring", "components", "component_packed"]
NodeSpacing = Literal["compact", "normal", "expanded"]

_SYMBOLS = (
    "circle", "square", "diamond", "triangle-up", "triangle-down", "pentagon",
    "hexagon", "star", "cross", "x", "hourglass", "bowtie",
)
_PALETTE = tuple(plotly.colors.qualitative.Dark24)
_OUTLINE_PALETTE = _PALETTE[8:] + _PALETTE[:8]
_CENTER_SYMBOLS = ("circle", "x", "cross", "diamond", "star")
_SPACING_CONFIG = {
    "compact": (0.10, 0.45),
    "normal": (0.18, 0.75),
    "expanded": (0.28, 1.20),
}


def filter_singletons(
    node_rows: list[dict[str, Any]], edge_rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    kept_ids = {row["sample_id"] for row in node_rows if int(row["cluster_size"]) > 1}
    filtered_nodes = [row for row in node_rows if row["sample_id"] in kept_ids]
    filtered_edges = [row for row in edge_rows if row["source"] in kept_ids and row["target"] in kept_ids]
    return filtered_nodes, filtered_edges


def _graph_from_rows(
    node_rows: Sequence[Mapping[str, Any]], edge_rows: Sequence[Mapping[str, Any]]
) -> nx.Graph:
    graph = nx.Graph()
    # NetworkX arrays follow insertion order, so sort first for reproducibility.
    for row in sorted(node_rows, key=lambda item: str(item["sample_id"])):
        graph.add_node(str(row["sample_id"]), **dict(row))
    valid_ids = set(graph)
    valid_edges = [
        row for row in edge_rows if str(row["source"]) in valid_ids and str(row["target"]) in valid_ids
    ]
    for row in sorted(
        valid_edges,
        key=lambda item: tuple(sorted((str(item["source"]), str(item["target"])))),
    ):
        graph.add_edge(str(row["source"]), str(row["target"]), distance=float(row["distance"]))
    return graph


def _stable_seed(value: str, seed: int) -> int:
    return int.from_bytes(hashlib.sha256(f"{seed}:{value}".encode()).digest()[:4], "big")


def _separate_overlapping_nodes(
    positions: Mapping[str, Sequence[float]], min_distance: float
) -> PositionMap:
    """Deterministically separate node centres while preserving the layout centroid."""
    ordered = sorted(positions)
    if len(ordered) < 2:
        return {
            node: (float(positions[node][0]), float(positions[node][1]))
            for node in ordered
        }
    original_center = (
        sum(float(positions[node][0]) for node in ordered) / len(ordered),
        sum(float(positions[node][1]) for node in ordered) / len(ordered),
    )
    current = {
        node: [float(positions[node][0]), float(positions[node][1])]
        for node in ordered
    }
    for _iteration in range(200):
        deltas = {node: [0.0, 0.0] for node in ordered}
        found_overlap = False
        for index, left in enumerate(ordered):
            for right in ordered[index + 1:]:
                dx = current[right][0] - current[left][0]
                dy = current[right][1] - current[left][1]
                distance = math.hypot(dx, dy)
                if distance >= min_distance:
                    continue
                found_overlap = True
                if distance < 1e-12:
                    angle = (
                        _stable_seed(f"{left}|{right}", 42) / (2**32 - 1)
                    ) * 2 * math.pi
                    unit_x, unit_y = math.cos(angle), math.sin(angle)
                else:
                    unit_x, unit_y = dx / distance, dy / distance
                push = (min_distance - distance) / 2 + 1e-6
                deltas[left][0] -= unit_x * push
                deltas[left][1] -= unit_y * push
                deltas[right][0] += unit_x * push
                deltas[right][1] += unit_y * push
        if not found_overlap:
            break
        for node in ordered:
            current[node][0] += deltas[node][0]
            current[node][1] += deltas[node][1]

    current_center = (
        sum(current[node][0] for node in ordered) / len(ordered),
        sum(current[node][1] for node in ordered) / len(ordered),
    )
    shift_x = original_center[0] - current_center[0]
    shift_y = original_center[1] - current_center[1]
    return {
        node: (current[node][0] + shift_x, current[node][1] + shift_y)
        for node in ordered
    }


def _normalized_component_layout(
    graph: nx.Graph,
    nodes: Sequence[str],
    seed: int,
    min_distance: float,
) -> PositionMap:
    ordered = sorted(nodes)
    if len(ordered) == 1:
        return {ordered[0]: (0.0, 0.0)}
    raw = nx.spring_layout(
        graph.subgraph(ordered),
        seed=_stable_seed("|".join(ordered), seed),
        k=1 / math.sqrt(len(ordered)),
    )
    extent = max(max(abs(float(point[0])), abs(float(point[1]))) for point in raw.values()) or 1.0
    radius = min(1.35, 0.65 + 0.12 * math.sqrt(len(ordered)))
    normalized = {
        node: (float(raw[node][0]) / extent * radius, float(raw[node][1]) / extent * radius)
        for node in ordered
    }
    separated = _separate_overlapping_nodes(normalized, min_distance)
    min_x = min(point[0] for point in separated.values())
    max_x = max(point[0] for point in separated.values())
    min_y = min(point[1] for point in separated.values())
    max_y = max(point[1] for point in separated.values())
    center_x, center_y = (min_x + max_x) / 2, (min_y + max_y) / 2
    return {node: (x - center_x, y - center_y) for node, (x, y) in separated.items()}


def _display_value(value: Any) -> str:
    if value is None or (isinstance(value, str) and not value.strip()):
        return MISSING_VALUE
    return str(value)


def _field_label(field: str) -> str:
    return field.replace("_", " ")


def _category_sort_key(value: str) -> tuple[int, int | str]:
    if value in AGE_RANGE_ORDER:
        return (0, AGE_RANGE_ORDER.index(value))
    if value == MISSING_VALUE:
        return (2, value)
    return (1, value.casefold())


def compute_network_layout(
    node_rows: Sequence[Mapping[str, Any]],
    edge_rows: Sequence[Mapping[str, Any]],
    *,
    mode: LayoutMode = "spring",
    group_by: str | None = None,
    seed: int = 42,
    node_spacing: NodeSpacing = "normal",
) -> PositionMap:
    """Return deterministic positions that callers can cache and reuse.

    ``spring`` preserves the original whole-graph behaviour. ``components`` (also
    ``component_packed``) lays out connected components independently. Supplying a
    field such as ``group_by='genotype'`` places component grids in group regions.
    """
    if node_spacing not in _SPACING_CONFIG:
        raise ValueError(f"Unknown node spacing: {node_spacing!r}")
    min_distance, component_gap = _SPACING_CONFIG[node_spacing]
    graph = _graph_from_rows(node_rows, edge_rows)
    if not graph:
        return {}
    if mode == "spring":
        raw = nx.spring_layout(graph, seed=seed, k=1 / math.sqrt(len(graph)))
        return _separate_overlapping_nodes(raw, min_distance)
    if mode not in ("components", "component_packed"):
        raise ValueError(f"Unknown layout mode: {mode!r}")
    if group_by is not None and any(group_by not in data for _, data in graph.nodes(data=True)):
        raise ValueError(f"Cannot group layout by missing field {group_by!r}")

    components = [tuple(sorted(component)) for component in nx.connected_components(graph)]
    components.sort(key=lambda component: (-len(component), component))
    grouped: dict[str, list[tuple[str, ...]]] = {}
    for component in components:
        if group_by is None:
            label = ""
        else:
            values = sorted({_display_value(graph.nodes[node].get(group_by)) for node in component})
            label = " / ".join(values)
        grouped.setdefault(label, []).append(component)

    positions: PositionMap = {}
    group_offset_x = 0.0
    for label in sorted(grouped, key=_category_sort_key):
        group_components = grouped[label]
        columns = max(1, math.ceil(math.sqrt(len(group_components))))
        rows = math.ceil(len(group_components) / columns)
        local_layouts = [
            _normalized_component_layout(graph, component, seed, min_distance)
            for component in group_components
        ]
        dimensions = []
        for local in local_layouts:
            x_values = [point[0] for point in local.values()]
            y_values = [point[1] for point in local.values()]
            dimensions.append(
                (
                    max(max(x_values) - min(x_values), min_distance),
                    max(max(y_values) - min(y_values), min_distance),
                )
            )
        column_widths = [
            max(
                (dimensions[index][0] for index in range(column, len(dimensions), columns)),
                default=min_distance,
            )
            for column in range(columns)
        ]
        row_heights = [
            max(
                dimensions[index][1]
                for index in range(row * columns, min((row + 1) * columns, len(dimensions)))
            )
            for row in range(rows)
        ]
        x_centers: list[float] = []
        cursor_x = group_offset_x
        for width in column_widths:
            x_centers.append(cursor_x + width / 2)
            cursor_x += width + component_gap
        total_height = sum(row_heights) + component_gap * max(rows - 1, 0)
        y_centers: list[float] = []
        cursor_y = total_height / 2
        for height in row_heights:
            y_centers.append(cursor_y - height / 2)
            cursor_y -= height + component_gap
        for index, component in enumerate(group_components):
            center_x = x_centers[index % columns]
            center_y = y_centers[index // columns]
            for node, (x, y) in local_layouts[index].items():
                positions[node] = (x + center_x, y + center_y)
        group_offset_x = cursor_x + component_gap * 2
    return positions


def ordered_categories(values: Iterable[Any]) -> tuple[str, ...]:
    """Return unique display values with age bins and missing values ordered sensibly."""
    return tuple(sorted({_display_value(value) for value in values}, key=_category_sort_key))


def build_category_mapping(
    values: Iterable[Any],
    *,
    channel: Literal["color", "symbol", "outline", "center"] = "color",
) -> dict[str, str]:
    """Build an input-order-independent mapping suitable for reuse across plots."""
    categories = ordered_categories(values)
    if channel == "symbol":
        return {value: _SYMBOLS[index % len(_SYMBOLS)] for index, value in enumerate(categories)}
    if channel == "center":
        return {
            value: _CENTER_SYMBOLS[index % len(_CENTER_SYMBOLS)]
            for index, value in enumerate(categories)
        }
    palette = _OUTLINE_PALETTE if channel == "outline" else _PALETTE
    # Sorted sequential assignment avoids hash collisions within a plot. To keep
    # colours identical across filtered subsets, build once from the complete node
    # table and pass the result through ``encoding_maps``.
    nonmissing = [value for value in categories if value != MISSING_VALUE]
    mapping = {value: palette[index % len(palette)] for index, value in enumerate(nonmissing)}
    if MISSING_VALUE in categories:
        mapping[MISSING_VALUE] = "#bdbdbd"
    return mapping


def _all_numeric(values: Sequence[str]) -> bool:
    try:
        return all(math.isfinite(float(value)) for value in values)
    except ValueError:
        return False


def get_encoding_warnings(
    node_rows: Sequence[Mapping[str, Any]],
    *,
    color_by: str | None = None,
    symbol_by: str | None = None,
    size_by: str | None = None,
    outline_by: str | None = None,
    center_by: str | None = None,
    size_order: Sequence[Any] | None = None,
) -> list[str]:
    """Describe encodings that are likely to be difficult to read."""
    result: list[str] = []
    for field, label, recommended_max in (
        (color_by, "colour", 20),
        (symbol_by, "shape", 8),
        (outline_by, "outline", 8),
        (center_by, "centre mark", 4),
    ):
        if field:
            count = len(ordered_categories(row.get(field) for row in node_rows))
            if count > recommended_max:
                result.append(
                    f"{field!r} has {count} categories; {label} is clearest with "
                    f"{recommended_max} or fewer"
                )
    if symbol_by:
        count = len(ordered_categories(row.get(symbol_by) for row in node_rows))
        if count > len(_SYMBOLS):
            result.append(
                f"{symbol_by!r} exceeds the {len(_SYMBOLS)} distinct marker shapes; "
                "filter categories or choose another shape field"
            )
    if center_by:
        count = len(ordered_categories(row.get(center_by) for row in node_rows))
        if count > len(_CENTER_SYMBOLS):
            result.append(
                f"{center_by!r} exceeds the {len(_CENTER_SYMBOLS)} distinct centre marks; "
                "filter categories or choose another centre field"
            )
    if size_by:
        values = [_display_value(row.get(size_by)) for row in node_rows]
        nonmissing = [value for value in values if value != MISSING_VALUE]
        is_age_range = all(value in AGE_RANGE_ORDER for value in nonmissing)
        if nonmissing and not _all_numeric(nonmissing) and not is_age_range and size_order is None:
            result.append(
                f"{size_by!r} is categorical and has no inherent size order; provide size_order "
                "or use it for colour/shape instead"
            )
    return result


def _require_fields(node_rows: Sequence[Mapping[str, Any]], fields: Iterable[str | None]) -> None:
    available = {field for row in node_rows for field in row}
    missing = sorted({field for field in fields if field and field not in available})
    if missing:
        raise ValueError("Unknown node metadata field(s): " + ", ".join(missing))


def _size_encoding(values: Sequence[str]) -> tuple[list[float], list[tuple[str, float]]]:
    nonmissing = [value for value in values if value != MISSING_VALUE]
    if not nonmissing:
        return [8.0] * len(values), [(MISSING_VALUE, 8.0)]
    if _all_numeric(nonmissing):
        numeric = [float(value) for value in nonmissing]
        low, high = min(numeric), max(numeric)

        def numeric_size(value: str) -> float:
            if value == MISSING_VALUE:
                return 8.0
            if low == high:
                return 13.0
            return 9.0 + (float(value) - low) / (high - low) * 10.0

        sizes = [numeric_size(value) for value in values]
        unique = sorted(set(numeric))
        selected = unique if len(unique) <= 4 else [unique[0], unique[len(unique) // 2], unique[-1]]
        legend = [(f"{value:g}", numeric_size(str(value))) for value in selected]
        if MISSING_VALUE in values:
            legend.append((MISSING_VALUE, 8.0))
        return sizes, legend
    categories = ordered_categories(values)
    if all(category in AGE_RANGE_ORDER or category == MISSING_VALUE for category in categories):
        mapping = {AGE_RANGE_ORDER[0]: 9.0, AGE_RANGE_ORDER[1]: 14.0, AGE_RANGE_ORDER[2]: 19.0}
        mapping[MISSING_VALUE] = 8.0
    else:
        raise ValueError(
            "Categorical node size requires an explicit order; use age_range, a numeric field, "
            "or pass size_order"
        )
    return [mapping[value] for value in values], [(value, mapping[value]) for value in categories]


def _ordered_size_encoding(
    values: Sequence[str], size_order: Sequence[Any]
) -> tuple[list[float], list[tuple[str, float]]]:
    order = [_display_value(value) for value in size_order]
    if len(order) != len(set(order)):
        raise ValueError("size_order contains duplicate values")
    present = set(values).difference({MISSING_VALUE})
    missing = sorted(present.difference(order), key=_category_sort_key)
    if missing:
        raise ValueError("size_order does not include: " + ", ".join(missing))
    ordered_present = [value for value in order if value in present]
    mapping = {
        value: 13.0 if len(ordered_present) == 1 else 9.0 + index / (len(ordered_present) - 1) * 10.0
        for index, value in enumerate(ordered_present)
    }
    if MISSING_VALUE in values:
        mapping[MISSING_VALUE] = 8.0
    legend_order = ordered_present + ([MISSING_VALUE] if MISSING_VALUE in values else [])
    return [mapping[value] for value in values], [(value, mapping[value]) for value in legend_order]


def _supplied_encoding_map(
    encoding_maps: Mapping[str, Mapping[str, str]], channel: str, field: str
) -> Mapping[str, str] | None:
    """Resolve a channel-qualified map, with field-only keys kept for compatibility."""
    return encoding_maps.get(f"{channel}:{field}") or encoding_maps.get(field)


def _legend_trace(
    *, name: str, group: str, title: str, color: str = "#777777", symbol: str = "circle",
    size: float = 11, line_color: str = "#333333", line_width: float = 1,
) -> go.Scatter:
    return go.Scatter(
        x=[None], y=[None], mode="markers",
        marker=dict(color=color, symbol=symbol, size=size, line=dict(color=line_color, width=line_width)),
        name=name, legendgroup=group, legendgrouptitle_text=title,
        hoverinfo="skip", showlegend=True,
    )


def _cluster_label_positions(
    ordered_nodes: Sequence[str],
    node_data: Sequence[Mapping[str, Any]],
    layout: Mapping[str, Sequence[float]],
) -> tuple[list[float], list[float], list[str]]:
    """Place one label just outside each cluster's node bounding box.

    Candidate positions on all four sides are scored against nodes from other
    clusters and labels already placed. This is intentionally deterministic so
    metadata restyling does not make labels jump around.
    """
    grouped_nodes: dict[str, list[str]] = {}
    for node, data in zip(ordered_nodes, node_data, strict=True):
        if int(data["cluster_size"]) <= 1:
            continue
        grouped_nodes.setdefault(_display_value(data.get("cluster_id")), []).append(node)
    if not ordered_nodes:
        return [], [], []

    all_x = [float(layout[node][0]) for node in ordered_nodes]
    all_y = [float(layout[node][1]) for node in ordered_nodes]
    extent = max(max(all_x) - min(all_x), max(all_y) - min(all_y), 1.0)
    offset = max(0.16, extent * 0.025)
    placed: list[tuple[float, float]] = []
    label_x: list[float] = []
    label_y: list[float] = []
    labels: list[str] = []

    for cluster_id in sorted(grouped_nodes, key=_category_sort_key):
        cluster_nodes = grouped_nodes[cluster_id]
        xs = [float(layout[node][0]) for node in cluster_nodes]
        ys = [float(layout[node][1]) for node in cluster_nodes]
        center_x = (min(xs) + max(xs)) / 2
        center_y = (min(ys) + max(ys)) / 2
        candidates = (
            (center_x, max(ys) + offset),
            (max(xs) + offset, center_y),
            (center_x, min(ys) - offset),
            (min(xs) - offset, center_y),
        )
        cluster_node_set = set(cluster_nodes)
        other_nodes = [node for node in ordered_nodes if node not in cluster_node_set]

        def clearance(candidate: tuple[float, float]) -> float:
            obstacles = [
                (float(layout[node][0]), float(layout[node][1])) for node in other_nodes
            ] + placed
            if not obstacles:
                return float("inf")
            return min(math.hypot(candidate[0] - x, candidate[1] - y) for x, y in obstacles)

        # ``max`` keeps the first (above-cluster) position when scores tie.
        chosen = max(candidates, key=clearance)
        label_x.append(chosen[0])
        label_y.append(chosen[1])
        labels.append(cluster_id)
        placed.append(chosen)
    return label_x, label_y, labels


def build_network_figure(
    node_rows: list[dict[str, Any]],
    edge_rows: list[dict[str, Any]],
    *,
    positions: Mapping[str, Sequence[float]] | None = None,
    layout_mode: LayoutMode = "spring",
    group_by: str | None = None,
    node_spacing: NodeSpacing = "normal",
    color_by: str | None = None,
    symbol_by: str | None = None,
    size_by: str | None = None,
    outline_by: str | None = None,
    center_by: str | None = None,
    hover_fields: Sequence[str] | None = None,
    encoding_maps: Mapping[str, Mapping[str, str]] | None = None,
    size_order: Sequence[Any] | None = None,
    show_cluster_labels: bool = False,
    height: int = 600,
) -> go.Figure:
    """Build an interactive, vector-safe network with independent metadata channels.

    Defaults retain the original cluster-coloured look. Callers can cache the result
    of :func:`compute_network_layout` and supply it as ``positions`` so restyling
    metadata never moves nodes. IDs are hover-only because node mode is ``markers``.
    Reusable maps should use keys such as ``color:location``, ``symbol:status`` and
    ``outline:status``; unqualified field keys remain accepted for compatibility.
    """
    requested_hover = list(hover_fields) if hover_fields is not None else [
        field for field in dict.fromkeys(field for row in node_rows for field in row)
        if field not in {"sample_id", "cluster_id", "cluster_size"}
    ]
    _require_fields(
        node_rows,
        (color_by, symbol_by, size_by, outline_by, center_by, group_by, *requested_hover),
    )
    graph = _graph_from_rows(node_rows, edge_rows)
    if positions is None:
        layout = compute_network_layout(
            node_rows,
            edge_rows,
            mode=layout_mode,
            group_by=group_by,
            node_spacing=node_spacing,
        )
    else:
        missing_positions = sorted(set(graph).difference(positions))
        if missing_positions:
            raise ValueError("Positions missing node(s): " + ", ".join(missing_positions))
        layout = {node: (float(positions[node][0]), float(positions[node][1])) for node in graph}

    edge_x: list[float | None] = []
    edge_y: list[float | None] = []
    for source, target in graph.edges():
        x0, y0 = layout[source]
        x1, y1 = layout[target]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, mode="lines", line=dict(width=1, color="#999999"),
        hoverinfo="none", name="Network edges", showlegend=False,
    )

    ordered_nodes = list(graph.nodes())
    node_data = [graph.nodes[node] for node in ordered_nodes]
    multi_cluster_ids = sorted({str(data["cluster_id"]) for data in node_data if int(data["cluster_size"]) > 1})
    legacy_palette = plotly.colors.qualitative.Plotly
    legacy_colors = {value: legacy_palette[index % len(legacy_palette)] for index, value in enumerate(multi_cluster_ids)}
    supplied_maps = encoding_maps or {}
    legends: list[go.Scatter] = []

    if color_by:
        values = [_display_value(data.get(color_by)) for data in node_data]
        mapping = dict(
            _supplied_encoding_map(supplied_maps, "color", color_by)
            or build_category_mapping(values, channel="color")
        )
        unknown = sorted(set(values).difference(mapping))
        if unknown:
            raise ValueError(f"Colour map for {color_by!r} lacks: {', '.join(unknown)}")
        node_colors = [mapping[value] for value in values]
        legends += [
            _legend_trace(
                name=value,
                group=f"color:{color_by}",
                title=f"Colour: {_field_label(color_by)}",
                color=mapping[value],
            )
            for value in ordered_categories(values)
        ]
    else:
        node_colors = [
            "#c7c7c7" if int(data["cluster_size"]) == 1 else legacy_colors[str(data["cluster_id"])]
            for data in node_data
        ]

    if symbol_by:
        values = [_display_value(data.get(symbol_by)) for data in node_data]
        if len(set(values)) > len(_SYMBOLS):
            raise ValueError(f"{symbol_by!r} has more than {len(_SYMBOLS)} categories and cannot be represented with distinct marker shapes")
        mapping = dict(
            _supplied_encoding_map(supplied_maps, "symbol", symbol_by)
            or build_category_mapping(values, channel="symbol")
        )
        unknown = sorted(set(values).difference(mapping))
        if unknown:
            raise ValueError(f"Symbol map for {symbol_by!r} lacks: {', '.join(unknown)}")
        node_symbols = [mapping[value] for value in values]
        legends += [
            _legend_trace(
                name=value,
                group=f"symbol:{symbol_by}",
                title=f"Shape: {_field_label(symbol_by)}",
                symbol=mapping[value],
            )
            for value in ordered_categories(values)
        ]
    else:
        node_symbols = ["circle"] * len(node_data)

    if size_by:
        values = [_display_value(data.get(size_by)) for data in node_data]
        node_sizes, size_legend = (
            _ordered_size_encoding(values, size_order)
            if size_order is not None
            else _size_encoding(values)
        )
        legends += [
            _legend_trace(
                name=label,
                group=f"size:{size_by}",
                title=f"Size: {_field_label(size_by)}",
                size=size,
            )
            for label, size in size_legend
        ]
    else:
        node_sizes = [8.0 if int(data["cluster_size"]) == 1 else 13.0 for data in node_data]

    if outline_by:
        values = [_display_value(data.get(outline_by)) for data in node_data]
        mapping = dict(
            _supplied_encoding_map(supplied_maps, "outline", outline_by)
            or build_category_mapping(values, channel="outline")
        )
        unknown = sorted(set(values).difference(mapping))
        if unknown:
            raise ValueError(f"Outline map for {outline_by!r} lacks: {', '.join(unknown)}")
        outline_colors = [mapping[value] for value in values]
        outline_widths = [2.5] * len(node_data)
        legends += [
            _legend_trace(
                name=value,
                group=f"outline:{outline_by}",
                title=f"Outline: {_field_label(outline_by)}",
                color="#ffffff",
                line_color=mapping[value],
                line_width=3,
            )
            for value in ordered_categories(values)
        ]
    else:
        outline_colors = ["#333333"] * len(node_data)
        outline_widths = [1.0] * len(node_data)

    center_trace: go.Scatter | None = None
    if center_by:
        values = [_display_value(data.get(center_by)) for data in node_data]
        if len(set(values)) > len(_CENTER_SYMBOLS):
            raise ValueError(
                f"{center_by!r} has more than {len(_CENTER_SYMBOLS)} categories and "
                "cannot be represented with distinct centre marks"
            )
        mapping = dict(
            _supplied_encoding_map(supplied_maps, "center", center_by)
            or build_category_mapping(values, channel="center")
        )
        unknown = sorted(set(values).difference(mapping))
        if unknown:
            raise ValueError(f"Centre map for {center_by!r} lacks: {', '.join(unknown)}")
        legends += [
            _legend_trace(
                name=value,
                group=f"center:{center_by}",
                title=f"Centre: {_field_label(center_by)}",
                color="#222222",
                symbol=mapping[value],
                size=8,
                line_width=0,
            )
            for value in ordered_categories(values)
        ]

    hover_text: list[str] = []
    for sample_id, data in zip(ordered_nodes, node_data, strict=True):
        lines = [f"<b>{html.escape(sample_id)}</b>"]
        lines.append(f"cluster: {html.escape(_display_value(data.get('cluster_id')))}")
        lines.append(f"cluster size: {html.escape(_display_value(data.get('cluster_size')))}")
        for field in requested_hover:
            lines.append(f"{html.escape(field)}: {html.escape(_display_value(data.get(field)))}")
        hover_text.append("<br>".join(lines))

    node_trace = go.Scatter(
        x=[layout[node][0] for node in ordered_nodes],
        y=[layout[node][1] for node in ordered_nodes],
        mode="markers",
        marker=dict(color=node_colors, size=node_sizes, symbol=node_symbols, line=dict(width=outline_widths, color=outline_colors)),
        text=hover_text, hovertemplate="%{text}<extra></extra>", name="Samples", showlegend=False,
    )
    if center_by:
        center_trace = go.Scatter(
            x=[layout[node][0] for node in ordered_nodes],
            y=[layout[node][1] for node in ordered_nodes],
            mode="markers",
            marker=dict(
                color="#222222",
                size=6,
                symbol=[mapping[value] for value in values],
                line=dict(width=0),
            ),
            text=hover_text,
            hovertemplate="%{text}<extra></extra>",
            name="Centre marks",
            showlegend=False,
        )
    warnings = get_encoding_warnings(
        node_rows,
        color_by=color_by,
        symbol_by=symbol_by,
        size_by=size_by,
        outline_by=outline_by,
        center_by=center_by,
        size_order=size_order,
    )
    annotations: list[dict[str, Any]] = []
    if group_by:
        grouped_nodes: dict[str, list[str]] = {}
        for node, data in zip(ordered_nodes, node_data, strict=True):
            grouped_nodes.setdefault(_display_value(data.get(group_by)), []).append(node)
        for group_label in sorted(grouped_nodes, key=_category_sort_key):
            group_nodes = grouped_nodes[group_label]
            group_x = [layout[node][0] for node in group_nodes]
            group_y = [layout[node][1] for node in group_nodes]
            annotations.append(
                {
                    "x": (min(group_x) + max(group_x)) / 2,
                    "y": max(group_y) + 0.55,
                    "text": f"{html.escape(_field_label(group_by))}: {html.escape(group_label)}",
                    "showarrow": False,
                    "font": {"size": 14},
                }
            )
    data_traces = [edge_trace, node_trace]
    if center_trace is not None:
        data_traces.append(center_trace)
    if show_cluster_labels:
        label_x, label_y, cluster_labels = _cluster_label_positions(
            ordered_nodes, node_data, layout
        )
        data_traces.append(
            go.Scatter(
                x=label_x,
                y=label_y,
                mode="text",
                text=cluster_labels,
                textposition="middle center",
                textfont=dict(size=12, color="#222222"),
                hoverinfo="skip",
                name="Cluster IDs",
                showlegend=False,
                cliponaxis=False,
            )
        )
    fig = go.Figure(data=[*data_traces, *legends])
    fig.update_layout(
        showlegend=bool(legends), legend=dict(tracegroupgap=12),
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        margin=dict(l=10, r=10, t=10, b=10), height=height,
        plot_bgcolor="rgba(0,0,0,0)",
        annotations=annotations,
        meta={"encoding_warnings": warnings, "positions": {node: list(layout[node]) for node in ordered_nodes}},
    )
    return fig


def build_small_multiples_figure(
    node_rows: list[dict[str, Any]],
    edge_rows: list[dict[str, Any]],
    fields: Sequence[str],
    *,
    positions: Mapping[str, Sequence[float]] | None = None,
    layout_mode: LayoutMode = "spring",
    group_by: str | None = None,
    node_spacing: NodeSpacing = "normal",
    hover_fields: Sequence[str] | None = None,
    encoding_maps: Mapping[str, Mapping[str, str]] | None = None,
    show_cluster_labels: bool = False,
    columns: int = 2,
) -> go.Figure:
    """Render up to four metadata views with identical node positions.

    Each panel uses colour—the most readily compared visual channel—for one field.
    The resulting Plotly figure remains vector-safe for SVG export.
    """
    selected_fields = list(dict.fromkeys(fields))
    if not selected_fields:
        raise ValueError("Select at least one metadata field for small multiples")
    if len(selected_fields) > 4:
        raise ValueError("Small multiples support at most four metadata fields")
    _require_fields(node_rows, selected_fields)
    if columns < 1:
        raise ValueError("Small-multiple column count must be positive")
    if positions is None:
        layout = compute_network_layout(
            node_rows,
            edge_rows,
            mode=layout_mode,
            group_by=group_by,
            node_spacing=node_spacing,
        )
    else:
        layout = {
            str(node): (float(point[0]), float(point[1]))
            for node, point in positions.items()
        }

    panel_columns = min(columns, len(selected_fields))
    panel_rows = math.ceil(len(selected_fields) / panel_columns)
    figure = make_subplots(
        rows=panel_rows,
        cols=panel_columns,
        subplot_titles=[_field_label(field) for field in selected_fields],
        horizontal_spacing=0.12,
        vertical_spacing=0.12,
    )
    for index, field in enumerate(selected_fields):
        row = index // panel_columns + 1
        column = index % panel_columns + 1
        panel = build_network_figure(
            node_rows,
            edge_rows,
            positions=layout,
            color_by=field,
            hover_fields=hover_fields,
            encoding_maps=encoding_maps,
            show_cluster_labels=show_cluster_labels,
            height=420,
        )
        for trace in panel.data:
            figure.add_trace(trace, row=row, col=column)

    figure.update_xaxes(visible=False)
    figure.update_yaxes(visible=False)
    figure.update_layout(
        showlegend=True,
        legend=dict(tracegroupgap=12, x=1.02, xanchor="left", y=1, yanchor="top"),
        margin=dict(l=30, r=300, t=60, b=30),
        width=max(1000, 520 * panel_columns + 300),
        height=max(520, 460 * panel_rows),
        plot_bgcolor="rgba(0,0,0,0)",
        meta={
            "small_multiple_fields": selected_fields,
            "positions": {node: list(point) for node, point in layout.items()},
        },
    )
    return figure


def export_figure_svg(figure: go.Figure, path: str | Path) -> Path:
    """Export editable vector SVG, rejecting any embedded raster fallback."""
    output_path = Path(path)
    svg = figure.to_image(format="svg")
    svg_bytes = svg.encode("utf-8") if isinstance(svg, str) else bytes(svg)
    lowered = svg_bytes.lower()
    if b"<svg" not in lowered:
        raise ValueError("Plotly did not produce an SVG document")
    if b"<image" in lowered or b"data:image/" in lowered:
        raise ValueError("SVG export contains an embedded raster image")
    output_path.write_bytes(svg_bytes)
    return output_path
