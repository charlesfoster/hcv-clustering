from pathlib import Path

import hcv_workflow


def test_parse_distance_table_with_header(tmp_path: Path) -> None:
    distances = tmp_path / "distances.csv"
    distances.write_text("ID1,ID2,Distance\nA,B,0.01\nA,C,0.02\n", encoding="utf-8")

    rows = hcv_workflow.parse_distance_table(distances)

    assert rows == [
        hcv_workflow.DistanceRow("A", "B", 0.01),
        hcv_workflow.DistanceRow("A", "C", 0.02),
    ]


def test_threshold_distances_excludes_self_links_by_default() -> None:
    rows = [
        hcv_workflow.DistanceRow("A", "A", 0.0),
        hcv_workflow.DistanceRow("A", "B", 0.014),
        hcv_workflow.DistanceRow("A", "C", 0.016),
    ]

    links = hcv_workflow.threshold_distances(rows, threshold=0.015)

    assert links == [hcv_workflow.LinkRow("A", "B", 0.014)]


def test_connected_components_include_isolates() -> None:
    links = [
        hcv_workflow.LinkRow("A", "B", 0.01),
        hcv_workflow.LinkRow("B", "C", 0.01),
    ]

    assignments = hcv_workflow.connected_components(links, nodes=["A", "B", "C", "D"])

    assert assignments["A"] == assignments["B"] == assignments["C"]
    assert assignments["D"] != assignments["A"]


def test_read_nodes_file_accepts_plain_lists(tmp_path: Path) -> None:
    nodes = tmp_path / "nodes.txt"
    nodes.write_text("A\nB\nA\n", encoding="utf-8")

    assert hcv_workflow.read_nodes_file(nodes) == ["A", "B"]


def test_write_clusters_reports_component_size(tmp_path: Path) -> None:
    output = tmp_path / "clusters.csv"
    assignments = {"A": "C0001", "B": "C0001", "D": "C0002"}

    hcv_workflow.write_clusters(output, assignments)

    assert output.read_text(encoding="utf-8").splitlines() == [
        "sample_id,cluster_id,cluster_size",
        "A,C0001,2",
        "B,C0001,2",
        "D,C0002,1",
    ]


def test_detected_pass_genotypes_and_merged_cluster_ids(tmp_path: Path) -> None:
    rows = [
        {"assignment_status": "pass", "assigned_genotype": "1a"},
        {"assignment_status": "qc_fail", "assigned_genotype": ""},
        {"assignment_status": "pass", "assigned_genotype": "3a"},
    ]
    assert hcv_workflow.detected_pass_genotypes(rows) == ["1a", "3a"]

    cluster_1a = tmp_path / "1a_clusters.csv"
    cluster_3a = tmp_path / "3a_clusters.csv"
    cluster_1a.write_text("sample_id,cluster_id,cluster_size\nA,C0001,1\n", encoding="utf-8")
    cluster_3a.write_text("sample_id,cluster_id,cluster_size\nB,C0001,1\n", encoding="utf-8")

    output = tmp_path / "merged.csv"
    hcv_workflow.merge_cluster_tables([("1a", cluster_1a), ("3a", cluster_3a)], output)

    assert output.read_text(encoding="utf-8").splitlines() == [
        "sample_id,genotype,cluster_id,source_cluster_id,cluster_size",
        "A,1a,1a_C0001,C0001,1",
        "B,3a,3a_C0001,C0001,1",
    ]
