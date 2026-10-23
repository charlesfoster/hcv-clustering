#!/usr/bin/env python3
"""Network plots for hcv-cluster: render per-genotype and combined cluster networks.

Shared by ``hcv-cluster run`` (flags on the run command) and the ``plot`` stage used by the
Nextflow workflow, so both draw the same figures from the same tables.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import metadata, stages, variants, viz, workflow


def add_plot_arguments(parser: argparse.ArgumentParser, adv) -> None:
    """The --plot-* options. ``adv`` hides an option's help unless advanced help is requested."""
    parser.add_argument(
        "--plot-network",
        metavar="CHOICE",
        choices=("none", "png", "html", "svg", "both", "all"),
        default="none",
        help=adv(
            "Render cluster network plots. png: static image. html: interactive. svg: "
            "editable vector image. both: PNG+HTML (backward compatible). all: PNG+HTML+SVG. "
            "Default: none."
        ),
    )
    parser.add_argument(
        "--plot-scope",
        metavar="CHOICE",
        choices=("per-genotype", "combined", "both"),
        default="per-genotype",
        help=adv(
            "Plot each genotype separately (default), all genotypes in one presentation-only "
            "network, or both. Combined plots never add cross-genotype links."
        ),
    )
    parser.add_argument(
        "--plot-combined-layout",
        metavar="CHOICE",
        choices=("packed", "by-genotype"),
        default="packed",
        help=adv(
            "Combined-network layout: compact connected-component packing (default), or "
            "spatially group packed components by genotype."
        ),
    )
    parser.add_argument(
        "--plot-node-spacing",
        metavar="CHOICE",
        choices=("compact", "normal", "expanded"),
        default="normal",
        help=adv(
            "Minimum collision-aware node separation: compact, normal (default), or expanded"
        ),
    )
    parser.add_argument(
        "--plot-color-by",
        metavar="FIELD",
        default=None,
        help=adv("Metadata field encoded as node colour"),
    )
    parser.add_argument(
        "--plot-symbol-by",
        metavar="FIELD",
        default=None,
        help=adv("Low-cardinality metadata field encoded as node shape"),
    )
    parser.add_argument(
        "--plot-size-by",
        metavar="FIELD",
        default=None,
        help=adv(
            "Numeric or intrinsically ordered metadata field encoded as node size; use "
            "--plot-size-order for other categorical fields"
        ),
    )
    parser.add_argument(
        "--plot-size-order",
        action="append",
        metavar="VALUE",
        default=None,
        help=adv(
            "Category order from smallest to largest for --plot-size-by; repeat once per value"
        ),
    )
    parser.add_argument(
        "--plot-outline-by",
        metavar="FIELD",
        default=None,
        help=adv("Low-cardinality metadata field encoded as node outline colour"),
    )
    parser.add_argument(
        "--plot-center-by",
        metavar="FIELD",
        default=None,
        help=adv("Low-cardinality metadata field encoded as a small centre mark"),
    )
    parser.add_argument(
        "--plot-small-multiple-field",
        action="append",
        metavar="FIELD",
        default=None,
        help=adv(
            "Render the same network in a colour panel for this metadata field; repeat for "
            "up to four fields. Cannot be combined with composite node encodings."
        ),
    )
    parser.add_argument(
        "--plot-hover-field",
        action="append",
        metavar="FIELD",
        default=None,
        help=adv(
            "Metadata field included in interactive hover text; repeat for multiple fields. "
            "Defaults to all supplied metadata; sample_id is always included."
        ),
    )
    parser.add_argument(
        "--plot-hide-singletons",
        action="store_true",
        help=adv("Omit singleton (unclustered) sequences from --plot-network output"),
    )


_PLOT_BUILTIN_FIELDS = frozenset(
    {
        "sample_id",
        "genotype",
        "cluster_id",
        "source_cluster_id",
        "cluster_size",
        variants.BASE_SAMPLE_FIELD,
        variants.MIXED_INFECTION_FIELD,
    }
)


def requested_plot_fields(args: argparse.Namespace) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(
            field
            for field in (
                args.plot_color_by,
                args.plot_symbol_by,
                args.plot_size_by,
                args.plot_outline_by,
                args.plot_center_by,
                *(args.plot_small_multiple_field or ()),
                *(args.plot_hover_field or ()),
            )
            if field
        )
    )


def validate_plot_fields(
    args: argparse.Namespace,
    metadata_rows: Sequence[Mapping[str, Any]],
) -> None:
    requested = set(requested_plot_fields(args))
    available = _PLOT_BUILTIN_FIELDS | set(metadata.metadata_fields(metadata_rows))
    unknown = sorted(requested.difference(available))
    if unknown:
        available_text = ", ".join(sorted(available))
        metadata_hint = " Supply the field in --metadata." if not metadata_rows else ""
        raise ValueError(
            f"Unknown plot metadata field(s): {', '.join(unknown)}. "
            f"Available fields: {available_text}.{metadata_hint}"
        )
    if args.plot_size_order and not args.plot_size_by:
        raise ValueError("--plot-size-order requires --plot-size-by")
    if args.plot_size_order and len(args.plot_size_order) != len(set(args.plot_size_order)):
        raise ValueError("--plot-size-order contains duplicate values")
    if args.plot_small_multiple_field:
        if len(dict.fromkeys(args.plot_small_multiple_field)) > 4:
            raise ValueError("--plot-small-multiple-field may select at most four fields")
        composite_fields = (
            args.plot_color_by,
            args.plot_symbol_by,
            args.plot_size_by,
            args.plot_outline_by,
            args.plot_center_by,
        )
        if any(composite_fields):
            raise ValueError(
                "Small-multiple fields cannot be combined with composite node encodings"
            )


def prepare_plot_rows(
    clusters_csv: Path,
    links_csv: Path,
    *,
    genotype: str | None,
    hide_singletons: bool,
    metadata_rows: Sequence[Mapping[str, Any]],
    cluster_labels: Mapping[str, str] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    node_rows: list[dict[str, Any]] = viz.apply_persistent_labels(
        stages.read_csv_rows(clusters_csv), cluster_labels
    )
    edge_rows: list[dict[str, Any]] = stages.read_csv_rows(links_csv)
    for row in node_rows:
        if genotype is not None:
            row.setdefault("genotype", genotype)
            row.setdefault("source_cluster_id", row.get("cluster_id", ""))
    if hide_singletons:
        node_rows, edge_rows = viz.filter_singletons(node_rows, edge_rows)
    if metadata_rows:
        node_rows, _report = metadata.join_metadata(node_rows, metadata_rows)
    return node_rows, edge_rows


def render_cluster_plot(
    plot_network: str,
    hide_singletons: bool,
    output_dir: Path,
    clusters_csv: Path,
    links_csv: Path,
    file_prefix: str,
    *,
    genotype: str | None = None,
    metadata_rows: Sequence[Mapping[str, Any]] = (),
    color_by: str | None = None,
    symbol_by: str | None = None,
    size_by: str | None = None,
    size_order: Sequence[str] | None = None,
    outline_by: str | None = None,
    center_by: str | None = None,
    small_multiple_fields: Sequence[str] = (),
    hover_fields: Sequence[str] | None = None,
    combined_layout: str | None = None,
    node_spacing: str = "normal",
    encoding_maps: Mapping[str, Mapping[str, str]] | None = None,
    cluster_labels: Mapping[str, str] | None = None,
) -> None:
    if plot_network == "none":
        return
    node_rows, edge_rows = prepare_plot_rows(
        clusters_csv,
        links_csv,
        genotype=genotype,
        hide_singletons=hide_singletons,
        metadata_rows=metadata_rows,
        cluster_labels=cluster_labels,
    )
    if not node_rows:
        return
    layout_mode = "component_packed" if combined_layout is not None else "spring"
    group_by = "genotype" if combined_layout == "by-genotype" else None
    requested_hover_list = [
        field for field in (hover_fields or ()) if field != "sample_id"
    ]
    if combined_layout is not None and "genotype" not in requested_hover_list:
        requested_hover_list.insert(0, "genotype")
    requested_hover = tuple(requested_hover_list)
    if small_multiple_fields:
        fig = viz.build_small_multiples_figure(
            node_rows,
            edge_rows,
            small_multiple_fields,
            layout_mode=layout_mode,
            group_by=group_by,
            node_spacing=node_spacing,
            hover_fields=requested_hover,
            encoding_maps=encoding_maps,
        )
    else:
        fig = viz.build_network_figure(
            node_rows,
            edge_rows,
            layout_mode=layout_mode,
            group_by=group_by,
            node_spacing=node_spacing,
            color_by=color_by,
            symbol_by=symbol_by,
            size_by=size_by,
            size_order=size_order,
            outline_by=outline_by,
            center_by=center_by,
            hover_fields=requested_hover,
            encoding_maps=encoding_maps,
        )
    if plot_network in ("png", "both", "all"):
        try:
            fig.write_image(output_dir / f"{file_prefix}network.png", scale=2)
        except Exception as exc:
            print(f"WARNING: could not export PNG network plot for {output_dir.name}: {exc}")
    if plot_network in ("html", "both", "all"):
        try:
            fig.write_html(output_dir / f"{file_prefix}network.html", include_plotlyjs=True)
        except Exception as exc:
            print(f"WARNING: could not export HTML network plot for {output_dir.name}: {exc}")
    if plot_network in ("svg", "all"):
        try:
            viz.export_figure_svg(fig, output_dir / f"{file_prefix}network.svg")
        except Exception as exc:
            print(f"WARNING: could not export SVG network plot for {output_dir.name}: {exc}")


def summarize_ids(sample_ids: Sequence[str], limit: int = 5) -> str:
    ordered = sorted(set(sample_ids))
    displayed = ", ".join(ordered[:limit])
    if len(ordered) > limit:
        displayed += f", ... (+{len(ordered) - limit} more)"
    return displayed


def warn_plot_metadata(
    node_rows: Sequence[Mapping[str, Any]],
    metadata_rows: Sequence[Mapping[str, Any]],
    args: argparse.Namespace,
) -> None:
    if metadata_rows:
        node_ids = {metadata.metadata_key(row) for row in node_rows}
        metadata_ids = {str(row["sample_id"]) for row in metadata_rows}
        missing = sorted(node_ids.difference(metadata_ids))
        unknown = sorted(metadata_ids.difference(node_ids))
        if missing:
            print(
                f"WARNING: metadata is missing for {len(missing)} plotted sample(s): "
                f"{summarize_ids(missing)}",
                file=sys.stderr,
            )
        if unknown:
            print(
                f"WARNING: metadata contains {len(unknown)} sample(s) not represented in "
                f"the plotted networks: {summarize_ids(unknown)}",
                file=sys.stderr,
            )
    warnings = viz.get_encoding_warnings(
        node_rows,
        color_by=args.plot_color_by,
        symbol_by=args.plot_symbol_by,
        size_by=args.plot_size_by,
        outline_by=args.plot_outline_by,
        center_by=args.plot_center_by,
        size_order=args.plot_size_order,
    )
    for field in args.plot_small_multiple_field or ():
        warnings.extend(
            viz.get_encoding_warnings(node_rows, color_by=field)
        )
    for warning in dict.fromkeys(warnings):
        print(f"WARNING: {warning}", file=sys.stderr)


def validate_plot_encodings(
    node_rows: Sequence[Mapping[str, Any]], args: argparse.Namespace
) -> None:
    warnings = viz.get_encoding_warnings(
        node_rows,
        color_by=args.plot_color_by,
        symbol_by=args.plot_symbol_by,
        size_by=args.plot_size_by,
        outline_by=args.plot_outline_by,
        center_by=args.plot_center_by,
        size_order=args.plot_size_order,
    )
    unusable = [
        warning
        for warning in warnings
        if "distinct marker shapes" in warning
        or "distinct centre marks" in warning
        or "has no inherent size order" in warning
    ]
    if unusable:
        raise ValueError("Invalid plot encoding: " + "; ".join(unusable))


def render_requested_plots(
    args: argparse.Namespace,
    outdir: Path,
    metadata_rows: Sequence[Mapping[str, Any]],
    cluster_tables: list[tuple[str, Path]],
    link_tables: list[tuple[str, Path]],
    snp_cluster_tables: list[tuple[str, Path]],
    snp_link_tables: list[tuple[str, Path]],
) -> None:
    """Render requested per-genotype and/or combined plots from completed tables."""
    if args.plot_network == "none":
        return

    metric_tables: list[
        tuple[str, list[tuple[str, Path]], list[tuple[str, Path]], Path, Path]
    ] = [("", cluster_tables, link_tables, outdir / "clusters.csv", outdir / "links.csv")]
    if args.distance == "both":
        metric_tables.append(
            (
                "snp_",
                snp_cluster_tables,
                snp_link_tables,
                outdir / "clusters.snp.csv",
                outdir / "links.snp.csv",
            )
        )

    # Persistent labels are tracked for the primary metric only; its per-genotype tables
    # carry no label column, so they borrow the labels from the merged table.
    primary_labels = viz.persistent_label_map(stages.read_csv_rows(outdir / "clusters.csv"))

    # Build mappings from the union of every plotted metric so categories keep the
    # same colour/shape/outline as users switch plots.
    union_by_sample: dict[str, dict[str, Any]] = {}
    for prefix, _clusters, _links, merged_clusters, merged_links in metric_tables:
        metric_nodes, _metric_edges = prepare_plot_rows(
            merged_clusters,
            merged_links,
            genotype=None,
            hide_singletons=args.plot_hide_singletons,
            metadata_rows=metadata_rows,
            cluster_labels=primary_labels if not prefix else None,
        )
        for row in metric_nodes:
            union_by_sample.setdefault(str(row["sample_id"]), row)
    union_nodes = list(union_by_sample.values())
    validate_plot_encodings(union_nodes, args)
    warn_plot_metadata(union_nodes, metadata_rows, args)

    encoding_maps: dict[str, Mapping[str, str]] = {}
    for channel, field in (
        ("color", args.plot_color_by),
        ("symbol", args.plot_symbol_by),
        ("outline", args.plot_outline_by),
        ("center", args.plot_center_by),
    ):
        if field:
            encoding_maps[f"{channel}:{field}"] = viz.build_category_mapping(
                (row.get(field) for row in union_nodes), channel=channel
            )
    for field in args.plot_small_multiple_field or ():
        encoding_maps[f"color:{field}"] = viz.build_category_mapping(
            (row.get(field) for row in union_nodes), channel="color"
        )
    effective_hover_fields = (
        tuple(args.plot_hover_field)
        if args.plot_hover_field is not None
        else metadata.metadata_fields(metadata_rows)
    )
    render_options = {
        "metadata_rows": metadata_rows,
        "color_by": args.plot_color_by,
        "symbol_by": args.plot_symbol_by,
        "size_by": args.plot_size_by,
        "size_order": args.plot_size_order,
        "outline_by": args.plot_outline_by,
        "center_by": args.plot_center_by,
        "small_multiple_fields": tuple(args.plot_small_multiple_field or ()),
        "hover_fields": effective_hover_fields,
        "node_spacing": args.plot_node_spacing,
        "encoding_maps": encoding_maps,
    }

    for file_prefix, per_gt_clusters, per_gt_links, merged_clusters, merged_links in metric_tables:
        if args.plot_scope in ("per-genotype", "both"):
            links_by_genotype = dict(per_gt_links)
            for genotype, clusters_path in per_gt_clusters:
                genotype_dir = outdir / "by_genotype" / workflow.safe_genotype_name(genotype)
                render_cluster_plot(
                    args.plot_network,
                    args.plot_hide_singletons,
                    genotype_dir,
                    clusters_path,
                    links_by_genotype[genotype],
                    file_prefix,
                    genotype=genotype,
                    cluster_labels=primary_labels if not file_prefix else None,
                    **render_options,
                )
        if args.plot_scope in ("combined", "both"):
            render_cluster_plot(
                args.plot_network,
                args.plot_hide_singletons,
                outdir,
                merged_clusters,
                merged_links,
                file_prefix,
                combined_layout=args.plot_combined_layout,
                cluster_labels=primary_labels if not file_prefix else None,
                **render_options,
            )
