from __future__ import annotations

import argparse
import csv
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from hcv_cluster import cli, plotting


def _write_rows(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_network_tables(
    clusters: Path,
    links: Path,
    *,
    combined: bool = False,
) -> None:
    fields = ("sample_id", "cluster_id", "cluster_size")
    rows: list[dict[str, object]] = [
        {"sample_id": "a", "cluster_id": "cluster_001", "cluster_size": 2},
        {"sample_id": "b", "cluster_id": "cluster_001", "cluster_size": 2},
    ]
    if combined:
        fields = ("sample_id", "genotype", "cluster_id", "source_cluster_id", "cluster_size")
        rows = [
            {
                **row,
                "genotype": "1a",
                "cluster_id": f"1a_{row['cluster_id']}",
                "source_cluster_id": row["cluster_id"],
            }
            for row in rows
        ]
    _write_rows(clusters, fields, rows)
    _write_rows(
        links,
        ("source", "target", "distance"),
        [{"source": "a", "target": "b", "distance": 0.01}],
    )


def _plot_args(**overrides: object) -> argparse.Namespace:
    values: dict[str, object] = {
        "plot_network": "html",
        "plot_scope": "per-genotype",
        "plot_combined_layout": "packed",
        "plot_node_spacing": "normal",
        "plot_hide_singletons": False,
        "plot_color_by": None,
        "plot_symbol_by": None,
        "plot_size_by": None,
        "plot_size_order": None,
        "plot_outline_by": None,
        "plot_center_by": None,
        "plot_small_multiple_field": None,
        "plot_hover_field": None,
        "distance": "tn93",
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_small_multiples_reject_composite_encodings() -> None:
    args = cli.build_parser().parse_args(
        [
            "run",
            "--input",
            "samples.fasta",
            "--plot-small-multiple-field",
            "location",
            "--plot-color-by",
            "location",
        ]
    )

    with pytest.raises(ValueError, match="cannot be combined"):
        plotting.validate_plot_fields(
            args, [{"sample_id": "A", "location": "P1"}]
        )


def test_unknown_plot_field_is_rejected_before_genotyping(tmp_path: Path, monkeypatch) -> None:
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("sample_id,location\na,P1\n", encoding="utf-8")
    called = False

    def fake_genotyping(_argv):
        nonlocal called
        called = True
        return 0

    monkeypatch.setattr(cli.genotyping, "main", fake_genotyping)
    args = cli.build_parser().parse_args(
        [
            "run", "--input", "samples.fasta", "--metadata", str(metadata),
            "--plot-color-by", "prison_name", "--outdir", str(tmp_path / "results"),
        ]
    )

    with pytest.raises(ValueError, match="Unknown plot metadata field.*prison_name"):
        cli.command_run(args)
    assert called is False


def test_dry_run_validates_metadata_but_does_not_persist_it(tmp_path: Path) -> None:
    metadata = tmp_path / "input_metadata.csv"
    metadata.write_text("sample_id,age\na,30\n", encoding="utf-8")
    outdir = tmp_path / "results"
    args = cli.build_parser().parse_args(
        [
            "run", "--input", "samples.fasta", "--metadata", str(metadata),
            "--outdir", str(outdir), "--dry-run",
        ]
    )

    assert cli.command_run(args) == 0
    assert not outdir.exists()


def test_metadata_is_persisted_normalized_before_pipeline_work(tmp_path: Path, monkeypatch) -> None:
    metadata = tmp_path / "input_metadata.csv"
    metadata.write_text("sample_id,age,location\na,31,P1\n", encoding="utf-8")
    outdir = tmp_path / "results"
    fasta = tmp_path / "samples.fasta"
    fasta.write_text(">a\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(cli.genotyping, "main", lambda _argv: 7)
    args = cli.build_parser().parse_args(
        [
            "run", "--input", str(fasta), "--metadata", str(metadata),
            "--outdir", str(outdir),
        ]
    )

    assert cli.command_run(args) == 7
    persisted = (outdir / "metadata.csv").read_text(encoding="utf-8")
    assert persisted.splitlines()[0] == "sample_id,age,location,age_range"
    assert "a,31,P1,31–60" in persisted


@pytest.mark.parametrize(
    ("plot_format", "png_calls", "html_calls", "svg_calls"),
    [
        ("png", 1, 0, 0),
        ("html", 0, 1, 0),
        ("svg", 0, 0, 1),
        ("both", 1, 1, 0),
        ("all", 1, 1, 1),
    ],
)
def test_render_cluster_plot_dispatches_formats(
    tmp_path: Path,
    monkeypatch,
    plot_format: str,
    png_calls: int,
    html_calls: int,
    svg_calls: int,
) -> None:
    clusters, links = tmp_path / "clusters.csv", tmp_path / "links.csv"
    _write_network_tables(clusters, links)
    figure = SimpleNamespace(write_image=Mock(), write_html=Mock())
    build = Mock(return_value=figure)
    export_svg = Mock()
    monkeypatch.setattr(plotting.viz, "build_network_figure", build)
    monkeypatch.setattr(plotting.viz, "export_figure_svg", export_svg)

    plotting.render_cluster_plot(
        plot_format,
        False,
        tmp_path,
        clusters,
        links,
        "snp_",
        genotype="1a",
        metadata_rows=[{"sample_id": "a", "location": "P1"}],
        color_by="location",
        hover_fields=("sample_id", "location"),
    )

    assert figure.write_image.call_count == png_calls
    assert figure.write_html.call_count == html_calls
    assert export_svg.call_count == svg_calls
    node_rows = build.call_args.args[0]
    assert node_rows[0]["genotype"] == "1a"
    assert node_rows[0]["location"] == "P1"
    assert node_rows[1]["location"] == "(missing)"
    assert build.call_args.kwargs["hover_fields"] == ("location",)
    if svg_calls:
        assert export_svg.call_args.args[1] == tmp_path / "snp_network.svg"


def test_combined_render_uses_merged_tables_and_by_genotype_layout(
    tmp_path: Path, monkeypatch
) -> None:
    _write_network_tables(tmp_path / "clusters.csv", tmp_path / "links.csv", combined=True)
    _write_network_tables(
        tmp_path / "clusters.snp.csv", tmp_path / "links.snp.csv", combined=True
    )
    per_clusters = tmp_path / "by_genotype" / "1a" / "clusters.csv"
    per_links = tmp_path / "by_genotype" / "1a" / "links.csv"
    _write_network_tables(per_clusters, per_links)
    snp_per_clusters = tmp_path / "by_genotype" / "1a" / "snp_clusters.csv"
    snp_per_links = tmp_path / "by_genotype" / "1a" / "snp_links.csv"
    _write_network_tables(snp_per_clusters, snp_per_links)
    render = Mock()
    monkeypatch.setattr(plotting, "render_cluster_plot", render)

    plotting.render_requested_plots(
        _plot_args(
            distance="both", plot_scope="combined", plot_combined_layout="by-genotype"
        ),
        tmp_path,
        [],
        [("1a", per_clusters)],
        [("1a", per_links)],
        [("1a", snp_per_clusters)],
        [("1a", snp_per_links)],
    )

    assert render.call_count == 2
    assert render.call_args_list[0].args[2:6] == (
        tmp_path,
        tmp_path / "clusters.csv",
        tmp_path / "links.csv",
        "",
    )
    assert render.call_args_list[1].args[2:6] == (
        tmp_path,
        tmp_path / "clusters.snp.csv",
        tmp_path / "links.snp.csv",
        "snp_",
    )
    assert all(
        call.kwargs["combined_layout"] == "by-genotype" for call in render.call_args_list
    )


def test_plot_scope_both_dispatches_per_genotype_and_combined(tmp_path: Path, monkeypatch) -> None:
    _write_network_tables(tmp_path / "clusters.csv", tmp_path / "links.csv", combined=True)
    per_clusters = tmp_path / "by_genotype" / "1a" / "clusters.csv"
    per_links = tmp_path / "by_genotype" / "1a" / "links.csv"
    _write_network_tables(per_clusters, per_links)
    render = Mock()
    monkeypatch.setattr(plotting, "render_cluster_plot", render)

    plotting.render_requested_plots(
        _plot_args(plot_scope="both"),
        tmp_path,
        [],
        [("1a", per_clusters)],
        [("1a", per_links)],
        [],
        [],
    )

    assert render.call_count == 2
    assert render.call_args_list[0].kwargs["genotype"] == "1a"
    assert render.call_args_list[1].kwargs["combined_layout"] == "packed"


def test_snp_only_combined_plot_uses_primary_root_paths(tmp_path: Path, monkeypatch) -> None:
    _write_network_tables(tmp_path / "clusters.csv", tmp_path / "links.csv", combined=True)
    render = Mock()
    monkeypatch.setattr(plotting, "render_cluster_plot", render)

    plotting.render_requested_plots(
        _plot_args(distance="snp", plot_scope="combined"),
        tmp_path,
        [],
        [],
        [],
        [],
        [],
    )

    assert render.call_count == 1
    assert render.call_args.args[2:6] == (
        tmp_path,
        tmp_path / "clusters.csv",
        tmp_path / "links.csv",
        "",
    )


@pytest.mark.parametrize(
    ("builder", "render_kwargs", "expected_args", "expected_kwargs"),
    [
        pytest.param(
            "build_network_figure",
            {"combined_layout": "by-genotype"},
            {},
            {
                "layout_mode": "component_packed",
                "group_by": "genotype",
                "hover_fields": ("genotype",),
            },
            id="combined-layout-component-packing",
        ),
        pytest.param(
            "build_small_multiples_figure",
            {
                "combined_layout": "packed",
                "small_multiple_fields": ("genotype",),
                "node_spacing": "expanded",
            },
            {2: ("genotype",)},
            {"layout_mode": "component_packed", "node_spacing": "expanded"},
            id="small-multiples-shared-renderer",
        ),
    ],
)
def test_render_cluster_plot_passes_layout_options_to_renderer(
    tmp_path: Path,
    monkeypatch,
    builder: str,
    render_kwargs: dict[str, object],
    expected_args: dict[int, object],
    expected_kwargs: dict[str, object],
) -> None:
    clusters, links = tmp_path / "clusters.csv", tmp_path / "links.csv"
    _write_network_tables(clusters, links, combined=True)
    figure = SimpleNamespace(write_image=Mock(), write_html=Mock())
    build = Mock(return_value=figure)
    monkeypatch.setattr(plotting.viz, builder, build)

    plotting.render_cluster_plot(
        "html", False, tmp_path, clusters, links, "", **render_kwargs
    )

    for index, value in expected_args.items():
        assert build.call_args.args[index] == value
    for name, value in expected_kwargs.items():
        assert build.call_args.kwargs[name] == value
