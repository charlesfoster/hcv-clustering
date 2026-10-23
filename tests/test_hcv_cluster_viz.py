from hcv_cluster import viz
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


@pytest.mark.parametrize(
    ("nodes", "edges", "expected_ids", "expected_edge_count"),
    [
        pytest.param(
            [
                {"sample_id": "A", "cluster_id": "C0001", "cluster_size": "2"},
                {"sample_id": "B", "cluster_id": "C0001", "cluster_size": "2"},
                {"sample_id": "C", "cluster_id": "C0002", "cluster_size": "1"},
            ],
            [{"source": "A", "target": "B", "distance": "0.001"}],
            {"A", "B"},
            1,
            id="drops-size-one-clusters",
        ),
        pytest.param(
            [
                {"sample_id": "A", "cluster_id": "C0001", "cluster_size": "1"},
                {"sample_id": "B", "cluster_id": "C0002", "cluster_size": "1"},
            ],
            [{"source": "A", "target": "B", "distance": "0.02"}],
            set(),
            0,
            id="drops-edges-touching-removed-nodes",
        ),
    ],
)
def test_filter_singletons(nodes, edges, expected_ids, expected_edge_count) -> None:
    filtered_nodes, filtered_edges = viz.filter_singletons(nodes, edges)

    assert {row["sample_id"] for row in filtered_nodes} == expected_ids
    assert filtered_edges == edges[:expected_edge_count]


def test_metadata_channels_have_independent_legends_and_hover_only_ids() -> None:
    nodes, edges = _metadata_network()
    figure = viz.build_network_figure(
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
        "Colour: location", "Shape: status", "Size: age range", "Outline: indigenous",
    }
    assert list(node_trace.marker.size[:3]) == [9.0, 14.0, 19.0]


def test_center_mark_adds_vector_overlay_and_independent_legend() -> None:
    nodes, edges = _metadata_network()

    figure = viz.build_network_figure(
        nodes,
        edges,
        color_by="location",
        symbol_by="age_range",
        outline_by="indigenous",
        center_by="status",
    )

    center_trace = figure.data[2]
    assert center_trace.name == "Centre marks"
    assert center_trace.mode == "markers"
    assert list(center_trace.x) == list(figure.data[1].x)
    assert center_trace.text[0].startswith("<b>A</b>")
    titles = {
        trace.legendgrouptitle.text
        for trace in figure.data[3:]
        if trace.legendgrouptitle and trace.legendgrouptitle.text
    }
    assert "Centre: status" in titles


def test_missing_values_are_explicit_in_hover_and_legend() -> None:
    nodes, edges = _metadata_network()
    figure = viz.build_network_figure(
        nodes, edges, color_by="location", hover_fields=["age_range"],
    )
    assert any(trace.name == "(missing)" for trace in figure.data[2:])
    assert "age_range: (missing)" in figure.data[1].text[3]


def test_category_order_and_mapping_are_stable() -> None:
    values = ["Prison B", "Prison A", "(missing)"]
    assert viz.build_category_mapping(values) == viz.build_category_mapping(reversed(values))
    assert viz.ordered_categories(["61+", "0–30", "31–60", "(missing)"]) == (
        "0–30", "31–60", "61+", "(missing)",
    )


def test_positions_can_be_reused_when_encodings_change() -> None:
    nodes, edges = _metadata_network()
    positions = viz.compute_network_layout(nodes, edges, mode="components")
    first = viz.build_network_figure(nodes, edges, positions=positions, color_by="location")
    second = viz.build_network_figure(nodes, edges, positions=positions, color_by="genotype")
    assert list(first.data[1].x) == list(second.data[1].x)
    assert list(first.data[1].y) == list(second.data[1].y)


def test_component_layout_is_deterministic_and_handles_combined_genotypes() -> None:
    nodes, edges = _metadata_network()
    positions = viz.compute_network_layout(
        nodes, edges, mode="component_packed", group_by="genotype",
    )
    reversed_positions = viz.compute_network_layout(
        list(reversed(nodes)), list(reversed(edges)), mode="component_packed", group_by="genotype",
    )
    assert positions == reversed_positions
    assert max(positions[node][0] for node in ("A", "B")) < min(
        positions[node][0] for node in ("C", "D")
    )


def test_collision_aware_layout_enforces_expanded_minimum_spacing() -> None:
    nodes = [
        {"sample_id": str(index), "cluster_id": "C1", "cluster_size": "12"}
        for index in range(12)
    ]
    positions = viz.compute_network_layout(
        nodes, [], mode="spring", node_spacing="expanded"
    )

    distances = [
        ((positions[left][0] - positions[right][0]) ** 2
         + (positions[left][1] - positions[right][1]) ** 2) ** 0.5
        for left in positions
        for right in positions
        if left < right
    ]
    assert min(distances) >= 0.279


def test_grouped_combined_figure_labels_each_genotype_region() -> None:
    nodes, edges = _metadata_network()

    figure = viz.build_network_figure(
        nodes,
        edges,
        layout_mode="components",
        group_by="genotype",
    )

    assert {annotation.text for annotation in figure.layout.annotations} == {
        "genotype: 1a",
        "genotype: 3a",
    }


def test_cluster_id_labels_are_one_vector_text_item_per_non_singleton_cluster() -> None:
    nodes, edges = _metadata_network()
    nodes.append({"sample_id": "E", "cluster_id": "1a_C0002", "cluster_size": "1"})

    figure = viz.build_network_figure(
        nodes, edges, show_cluster_labels=True
    )

    label_trace = next(trace for trace in figure.data if trace.name == "Cluster IDs")
    assert label_trace.mode == "text"
    # The singleton cluster (1a_C0002) gets no label.
    assert list(label_trace.text) == ["1a_C0001", "3a_C0001"]
    assert len(label_trace.x) == 2
    assert label_trace.hoverinfo == "skip"


def test_shape_cardinality_is_rejected_before_ambiguous_symbol_reuse() -> None:
    nodes = [
        {"sample_id": str(index), "cluster_id": str(index), "cluster_size": "1", "kind": str(index)}
        for index in range(13)
    ]
    warnings = viz.get_encoding_warnings(nodes, symbol_by="kind")
    assert any("distinct marker shapes" in warning for warning in warnings)
    with pytest.raises(ValueError, match="distinct marker shapes"):
        viz.build_network_figure(nodes, [], symbol_by="kind")


def test_unknown_encoding_field_and_incomplete_positions_are_rejected() -> None:
    nodes, edges = _metadata_network()
    with pytest.raises(ValueError, match="Unknown node metadata"):
        viz.build_network_figure(nodes, edges, color_by="not_a_field")
    with pytest.raises(ValueError, match="Positions missing"):
        viz.build_network_figure(nodes, edges, positions={"A": (0, 0)})


def test_channel_qualified_maps_support_same_field_on_multiple_channels() -> None:
    nodes, edges = _metadata_network()
    maps = {
        "color:status": {"yes": "#ff0000", "no": "#0000ff"},
        "symbol:status": {"yes": "square", "no": "diamond"},
        "outline:status": {"yes": "#111111", "no": "#eeeeee"},
    }
    figure = viz.build_network_figure(
        nodes, edges, color_by="status", symbol_by="status", outline_by="status",
        encoding_maps=maps,
    )
    markers = figure.data[1].marker
    assert list(markers.color) == ["#ff0000", "#0000ff", "#ff0000", "#0000ff"]
    assert list(markers.symbol) == ["square", "diamond", "square", "diamond"]
    assert list(markers.line.color) == ["#111111", "#eeeeee", "#111111", "#eeeeee"]


def test_unordered_categorical_size_requires_explicit_order() -> None:
    nodes, edges = _metadata_network()
    assert any("no inherent size order" in warning for warning in viz.get_encoding_warnings(nodes, size_by="status"))
    with pytest.raises(ValueError, match="explicit order"):
        viz.build_network_figure(nodes, edges, size_by="status")
    figure = viz.build_network_figure(
        nodes, edges, size_by="status", size_order=["no", "yes"],
    )
    assert list(figure.data[1].marker.size) == [19.0, 9.0, 19.0, 9.0]


def test_small_multiples_reuse_positions_and_keep_hover_only_labels() -> None:
    nodes, edges = _metadata_network()
    positions = viz.compute_network_layout(nodes, edges, mode="components")

    figure = viz.build_small_multiples_figure(
        nodes,
        edges,
        ["location", "age_range"],
        positions=positions,
        hover_fields=["genotype"],
        show_cluster_labels=True,
    )

    sample_traces = [trace for trace in figure.data if trace.name == "Samples"]
    assert len(sample_traces) == 2
    assert list(sample_traces[0].x) == list(sample_traces[1].x)
    assert all(trace.mode == "markers" for trace in sample_traces)
    assert figure.layout.meta["small_multiple_fields"] == ["location", "age_range"]
    assert {annotation.text for annotation in figure.layout.annotations} == {
        "location", "age range",
    }
    # One cluster-label trace is repeated per panel.
    label_traces = [trace for trace in figure.data if trace.name == "Cluster IDs"]
    assert len(label_traces) == 2
    assert all(list(trace.text) == ["1a_C0001", "3a_C0001"] for trace in label_traces)


def test_svg_helper_rejects_raster_and_writes_vector(monkeypatch, tmp_path) -> None:
    figure = viz.build_network_figure(*_metadata_network())
    monkeypatch.setattr(figure, "to_image", lambda **_: b'<svg><g class="trace"><path/><text>Legend</text></g></svg>')
    output = viz.export_figure_svg(figure, tmp_path / "network.svg")
    assert "<path" in output.read_text(encoding="utf-8")

    monkeypatch.setattr(figure, "to_image", lambda **_: b'<svg><image href="data:image/png;base64,x"/></svg>')
    with pytest.raises(ValueError, match="embedded raster"):
        viz.export_figure_svg(figure, tmp_path / "bad.svg")
