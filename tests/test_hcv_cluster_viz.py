import hcv_cluster_viz
import pytest


def _metadata_network():
    nodes = [
        {
            "sample_id": "A", "cluster_id": "1a_C0001", "cluster_size": "2",
            "genotype": "1a", "location": "P1", "status": "yes", "age_range": "0–30",
            "indigenous": "yes",
        },
        {
            "sample_id": "B", "cluster_id": "1a_C0001", "cluster_size": "2",
            "genotype": "1a", "location": "P2", "status": "no", "age_range": "31–60",
            "indigenous": "no",
        },
        {
            "sample_id": "C", "cluster_id": "3a_C0001", "cluster_size": "2",
            "genotype": "3a", "location": "P1", "status": "yes", "age_range": "61+",
            "indigenous": "no",
        },
        {
            "sample_id": "D", "cluster_id": "3a_C0001", "cluster_size": "2",
            "genotype": "3a", "location": "(missing)", "status": "no", "age_range": "(missing)",
            "indigenous": "(missing)",
        },
    ]
    edges = [
        {"source": "A", "target": "B", "distance": "0.001"},
        {"source": "C", "target": "D", "distance": "0.002"},
    ]
    return nodes, edges


def test_filter_singletons_drops_size_one_clusters() -> None:
    nodes = [
        {"sample_id": "A", "cluster_id": "C0001", "cluster_size": "2"},
        {"sample_id": "B", "cluster_id": "C0001", "cluster_size": "2"},
        {"sample_id": "C", "cluster_id": "C0002", "cluster_size": "1"},
    ]
    edges = [
        {"source": "A", "target": "B", "distance": "0.001"},
    ]

    filtered_nodes, filtered_edges = hcv_cluster_viz.filter_singletons(nodes, edges)

    assert {row["sample_id"] for row in filtered_nodes} == {"A", "B"}
    assert filtered_edges == edges


def test_filter_singletons_drops_edges_touching_removed_nodes() -> None:
    nodes = [
        {"sample_id": "A", "cluster_id": "C0001", "cluster_size": "1"},
        {"sample_id": "B", "cluster_id": "C0002", "cluster_size": "1"},
    ]
    edges = [
        {"source": "A", "target": "B", "distance": "0.02"},
    ]

    filtered_nodes, filtered_edges = hcv_cluster_viz.filter_singletons(nodes, edges)

    assert filtered_nodes == []
    assert filtered_edges == []


def test_build_network_figure_has_edge_and_node_trace() -> None:
    nodes = [
        {"sample_id": "A", "cluster_id": "C0001", "cluster_size": "2"},
        {"sample_id": "B", "cluster_id": "C0001", "cluster_size": "2"},
    ]
    edges = [
        {"source": "A", "target": "B", "distance": "0.001"},
    ]

    figure = hcv_cluster_viz.build_network_figure(nodes, edges)

    assert len(figure.data) == 2
    assert len(figure.data[1].x) == 2
    assert figure.data[1].mode == "markers"
    assert figure.layout.showlegend is False


def test_metadata_channels_have_independent_legends_and_hover_only_ids() -> None:
    nodes, edges = _metadata_network()
    figure = hcv_cluster_viz.build_network_figure(
        nodes,
        edges,
        color_by="location",
        symbol_by="status",
        size_by="age_range",
        outline_by="indigenous",
        hover_fields=["genotype", "location"],
    )
    node_trace = figure.data[1]
    assert node_trace.mode == "markers"
    assert "text" not in node_trace.mode
    assert node_trace.text[0].startswith("<b>A</b><br>cluster:")
    assert "genotype: 1a" in node_trace.text[0]
    assert "location: P1" in node_trace.text[0]
    assert figure.layout.showlegend is True
    titles = {
        trace.legendgrouptitle.text
        for trace in figure.data[2:]
        if trace.legendgrouptitle and trace.legendgrouptitle.text
    }
    assert titles == {
        "Colour: location", "Shape: status", "Size: age_range", "Outline: indigenous",
    }
    assert list(node_trace.marker.size[:3]) == [9.0, 14.0, 19.0]


def test_missing_values_are_explicit_in_hover_and_legend() -> None:
    nodes, edges = _metadata_network()
    figure = hcv_cluster_viz.build_network_figure(
        nodes, edges, color_by="location", hover_fields=["age_range"],
    )
    assert any(trace.name == "(missing)" for trace in figure.data[2:])
    assert "age_range: (missing)" in figure.data[1].text[3]


def test_category_order_and_mapping_are_stable() -> None:
    values = ["Prison B", "Prison A", "(missing)"]
    assert hcv_cluster_viz.build_category_mapping(values) == hcv_cluster_viz.build_category_mapping(reversed(values))
    assert hcv_cluster_viz.ordered_categories(["61+", "0–30", "31–60", "(missing)"]) == (
        "0–30", "31–60", "61+", "(missing)",
    )


def test_positions_can_be_reused_when_encodings_change() -> None:
    nodes, edges = _metadata_network()
    positions = hcv_cluster_viz.compute_network_layout(nodes, edges, mode="components")
    first = hcv_cluster_viz.build_network_figure(nodes, edges, positions=positions, color_by="location")
    second = hcv_cluster_viz.build_network_figure(nodes, edges, positions=positions, color_by="genotype")
    assert list(first.data[1].x) == list(second.data[1].x)
    assert list(first.data[1].y) == list(second.data[1].y)


def test_component_layout_is_deterministic_and_handles_combined_genotypes() -> None:
    nodes, edges = _metadata_network()
    positions = hcv_cluster_viz.compute_network_layout(
        nodes, edges, mode="component_packed", group_by="genotype",
    )
    reversed_positions = hcv_cluster_viz.compute_network_layout(
        list(reversed(nodes)), list(reversed(edges)), mode="component_packed", group_by="genotype",
    )
    assert positions == reversed_positions
    assert max(positions[node][0] for node in ("A", "B")) < min(
        positions[node][0] for node in ("C", "D")
    )


def test_grouped_combined_figure_labels_each_genotype_region() -> None:
    nodes, edges = _metadata_network()

    figure = hcv_cluster_viz.build_network_figure(
        nodes,
        edges,
        layout_mode="components",
        group_by="genotype",
    )

    assert {annotation.text for annotation in figure.layout.annotations} == {
        "genotype: 1a",
        "genotype: 3a",
    }


def test_shape_cardinality_is_rejected_before_ambiguous_symbol_reuse() -> None:
    nodes = [
        {"sample_id": str(index), "cluster_id": str(index), "cluster_size": "1", "kind": str(index)}
        for index in range(13)
    ]
    warnings = hcv_cluster_viz.get_encoding_warnings(nodes, symbol_by="kind")
    assert any("distinct marker shapes" in warning for warning in warnings)
    with pytest.raises(ValueError, match="distinct marker shapes"):
        hcv_cluster_viz.build_network_figure(nodes, [], symbol_by="kind")


def test_unknown_encoding_field_and_incomplete_positions_are_rejected() -> None:
    nodes, edges = _metadata_network()
    with pytest.raises(ValueError, match="Unknown node metadata"):
        hcv_cluster_viz.build_network_figure(nodes, edges, color_by="not_a_field")
    with pytest.raises(ValueError, match="Positions missing"):
        hcv_cluster_viz.build_network_figure(nodes, edges, positions={"A": (0, 0)})


def test_channel_qualified_maps_support_same_field_on_multiple_channels() -> None:
    nodes, edges = _metadata_network()
    maps = {
        "color:status": {"yes": "#ff0000", "no": "#0000ff"},
        "symbol:status": {"yes": "square", "no": "diamond"},
        "outline:status": {"yes": "#111111", "no": "#eeeeee"},
    }
    figure = hcv_cluster_viz.build_network_figure(
        nodes, edges, color_by="status", symbol_by="status", outline_by="status",
        encoding_maps=maps,
    )
    markers = figure.data[1].marker
    assert list(markers.color) == ["#ff0000", "#0000ff", "#ff0000", "#0000ff"]
    assert list(markers.symbol) == ["square", "diamond", "square", "diamond"]
    assert list(markers.line.color) == ["#111111", "#eeeeee", "#111111", "#eeeeee"]


def test_unordered_categorical_size_requires_explicit_order() -> None:
    nodes, edges = _metadata_network()
    assert any("no inherent size order" in warning for warning in hcv_cluster_viz.get_encoding_warnings(nodes, size_by="status"))
    with pytest.raises(ValueError, match="explicit order"):
        hcv_cluster_viz.build_network_figure(nodes, edges, size_by="status")
    figure = hcv_cluster_viz.build_network_figure(
        nodes, edges, size_by="status", size_order=["no", "yes"],
    )
    assert list(figure.data[1].marker.size) == [19.0, 9.0, 19.0, 9.0]


def test_svg_helper_rejects_raster_and_writes_vector(monkeypatch, tmp_path) -> None:
    figure = hcv_cluster_viz.build_network_figure(*_metadata_network())
    monkeypatch.setattr(figure, "to_image", lambda **_: b'<svg><g class="trace"><path/><text>Legend</text></g></svg>')
    output = hcv_cluster_viz.export_figure_svg(figure, tmp_path / "network.svg")
    assert "<path" in output.read_text(encoding="utf-8")

    monkeypatch.setattr(figure, "to_image", lambda **_: b'<svg><image href="data:image/png;base64,x"/></svg>')
    with pytest.raises(ValueError, match="embedded raster"):
        hcv_cluster_viz.export_figure_svg(figure, tmp_path / "bad.svg")
