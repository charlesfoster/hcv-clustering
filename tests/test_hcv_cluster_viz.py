import hcv_cluster_viz


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
