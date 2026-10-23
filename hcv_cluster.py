#!/usr/bin/env python3
"""User-facing CLI for HCV genomic clustering with genotype-aware thresholds."""

from __future__ import annotations

import argparse
import contextlib
import csv
import datetime
import logging
import os
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import assign_hcv_genotypes_from_fasta
import hcv_cluster_audit
import hcv_cluster_labels
import hcv_cluster_metadata
import hcv_cluster_prep
import hcv_cluster_stages
import hcv_cluster_variants
import hcv_cluster_viz
import hcv_workflow


from hcv_cluster_methods import (  # noqa: F401  (re-exported for existing callers)
    REGION_THRESHOLDS_TN93,
    REGION_THRESHOLDS_SNP,
    EVIDENCE_ANCHORED_REGIONS_TN93,
    EVIDENCE_ANCHORED_REGIONS_SNP,
    REGION_THRESHOLDS,
    FALLBACK_THRESHOLD,
    FALLBACK_THRESHOLD_SNP,
    canonical_threshold_region,
    resolve_threshold,
    region_threshold_is_evidence_based,
    region_threshold_is_explicit,
    compute_snp_distances_detailed,
    compute_snp_distances,
    write_snp_csv,
    write_snp_links_csv,
    _snp_count_summary,
)


@contextlib.contextmanager
def _quiet_logging():
    logging.disable(logging.INFO)
    try:
        yield
    finally:
        logging.disable(logging.NOTSET)


def _read_qc_stats(qc_path: Path) -> dict[str, int]:
    with qc_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {
        "total": len(rows),
        "passed": sum(1 for r in rows if r.get("passed_qc") == "true"),
    }


def _read_cluster_stats(clusters_path: Path) -> dict[str, int]:
    with clusters_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {"total": 0, "singletons": 0, "n_multi": 0, "largest": 0}
    singletons = sum(1 for r in rows if int(r["cluster_size"]) == 1)
    multi_ids = {r["cluster_id"] for r in rows if int(r["cluster_size"]) > 1}
    largest = max(int(r["cluster_size"]) for r in rows)
    return {"total": len(rows), "singletons": singletons, "n_multi": len(multi_ids), "largest": largest}


def _fasta_has_records(path: Path) -> bool:
    with path.open(encoding="utf-8") as handle:
        return any(line.startswith(">") for line in handle)


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


_PLOT_BUILTIN_FIELDS = frozenset(
    {
        "sample_id",
        "genotype",
        "cluster_id",
        "source_cluster_id",
        "cluster_size",
        hcv_cluster_variants.BASE_SAMPLE_FIELD,
        hcv_cluster_variants.MIXED_INFECTION_FIELD,
    }
)


def _requested_plot_fields(args: argparse.Namespace) -> tuple[str, ...]:
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


def _validate_plot_fields(
    args: argparse.Namespace,
    metadata_rows: Sequence[Mapping[str, Any]],
) -> None:
    requested = set(_requested_plot_fields(args))
    available = _PLOT_BUILTIN_FIELDS | set(hcv_cluster_metadata.metadata_fields(metadata_rows))
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


def _prepare_plot_rows(
    clusters_csv: Path,
    links_csv: Path,
    *,
    genotype: str | None,
    hide_singletons: bool,
    metadata_rows: Sequence[Mapping[str, Any]],
    cluster_labels: Mapping[str, str] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    node_rows: list[dict[str, Any]] = hcv_cluster_viz.apply_persistent_labels(
        _read_csv_rows(clusters_csv), cluster_labels
    )
    edge_rows: list[dict[str, Any]] = _read_csv_rows(links_csv)
    for row in node_rows:
        if genotype is not None:
            row.setdefault("genotype", genotype)
            row.setdefault("source_cluster_id", row.get("cluster_id", ""))
    if hide_singletons:
        node_rows, edge_rows = hcv_cluster_viz.filter_singletons(node_rows, edge_rows)
    if metadata_rows:
        node_rows, _report = hcv_cluster_metadata.join_metadata(node_rows, metadata_rows)
    return node_rows, edge_rows


def _render_cluster_plot(
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
    node_rows, edge_rows = _prepare_plot_rows(
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
        fig = hcv_cluster_viz.build_small_multiples_figure(
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
        fig = hcv_cluster_viz.build_network_figure(
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
            hcv_cluster_viz.export_figure_svg(fig, output_dir / f"{file_prefix}network.svg")
        except Exception as exc:
            print(f"WARNING: could not export SVG network plot for {output_dir.name}: {exc}")


def _summarize_ids(sample_ids: Sequence[str], limit: int = 5) -> str:
    ordered = sorted(set(sample_ids))
    displayed = ", ".join(ordered[:limit])
    if len(ordered) > limit:
        displayed += f", ... (+{len(ordered) - limit} more)"
    return displayed


def _warn_plot_metadata(
    node_rows: Sequence[Mapping[str, Any]],
    metadata_rows: Sequence[Mapping[str, Any]],
    args: argparse.Namespace,
) -> None:
    if metadata_rows:
        node_ids = {hcv_cluster_metadata.metadata_key(row) for row in node_rows}
        metadata_ids = {str(row["sample_id"]) for row in metadata_rows}
        missing = sorted(node_ids.difference(metadata_ids))
        unknown = sorted(metadata_ids.difference(node_ids))
        if missing:
            print(
                f"WARNING: metadata is missing for {len(missing)} plotted sample(s): "
                f"{_summarize_ids(missing)}",
                file=sys.stderr,
            )
        if unknown:
            print(
                f"WARNING: metadata contains {len(unknown)} sample(s) not represented in "
                f"the plotted networks: {_summarize_ids(unknown)}",
                file=sys.stderr,
            )
    warnings = hcv_cluster_viz.get_encoding_warnings(
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
            hcv_cluster_viz.get_encoding_warnings(node_rows, color_by=field)
        )
    for warning in dict.fromkeys(warnings):
        print(f"WARNING: {warning}", file=sys.stderr)


def _validate_plot_encodings(
    node_rows: Sequence[Mapping[str, Any]], args: argparse.Namespace
) -> None:
    warnings = hcv_cluster_viz.get_encoding_warnings(
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


def _render_requested_plots(
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
    primary_labels = hcv_cluster_viz.persistent_label_map(_read_csv_rows(outdir / "clusters.csv"))

    # Build mappings from the union of every plotted metric so categories keep the
    # same colour/shape/outline as users switch plots.
    union_by_sample: dict[str, dict[str, Any]] = {}
    for prefix, _clusters, _links, merged_clusters, merged_links in metric_tables:
        metric_nodes, _metric_edges = _prepare_plot_rows(
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
    _validate_plot_encodings(union_nodes, args)
    _warn_plot_metadata(union_nodes, metadata_rows, args)

    encoding_maps: dict[str, Mapping[str, str]] = {}
    for channel, field in (
        ("color", args.plot_color_by),
        ("symbol", args.plot_symbol_by),
        ("outline", args.plot_outline_by),
        ("center", args.plot_center_by),
    ):
        if field:
            encoding_maps[f"{channel}:{field}"] = hcv_cluster_viz.build_category_mapping(
                (row.get(field) for row in union_nodes), channel=channel
            )
    for field in args.plot_small_multiple_field or ():
        encoding_maps[f"color:{field}"] = hcv_cluster_viz.build_category_mapping(
            (row.get(field) for row in union_nodes), channel="color"
        )
    effective_hover_fields = (
        tuple(args.plot_hover_field)
        if args.plot_hover_field is not None
        else hcv_cluster_metadata.metadata_fields(metadata_rows)
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
                genotype_dir = outdir / "by_genotype" / hcv_workflow.safe_genotype_name(genotype)
                _render_cluster_plot(
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
            _render_cluster_plot(
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


def _label_track(args: argparse.Namespace) -> str:
    """The primary metric drives the cluster tables, so it is the tracked clustering."""
    return "snp" if args.distance == "snp" else "tn93"


def _label_definition(args: argparse.Namespace) -> dict[str, object]:
    return hcv_cluster_stages.label_definition(
        _label_track(args),
        args.region,
        args.threshold,
        args.region_strategy,
        args.ambiguities,
        args.ambiguity_fraction,
        args.min_overlap,
        args.extra_tn93_args,
        args.snp_count_threshold,
    )


def _check_label_options(args: argparse.Namespace) -> None:
    """Fail before any pipeline work if the registry cannot be used for this run."""
    if args.update_registry and args.registry is None:
        raise ValueError("--update-registry requires --registry")
    if args.label_date is not None:
        try:
            args.label_date = datetime.date.fromisoformat(args.label_date).isoformat()
        except ValueError:
            raise ValueError(f"--label-date must be YYYY-MM-DD, got {args.label_date!r}") from None
    if args.registry is not None and args.registry.is_dir():
        hcv_cluster_labels.load_registry(args.registry).check_definition(
            _label_track(args), _label_definition(args)
        )


def _apply_cluster_labels(args: argparse.Namespace, outdir: Path, clusters_csv: Path) -> None:
    """Add persistent ``cluster_label`` to the merged cluster table and write label reports."""
    hcv_cluster_stages.label_clusters(
        clusters_csv,
        registry=args.registry,
        update_registry=args.update_registry,
        track=_label_track(args),
        definition=_label_definition(args),
        prefix=args.label_prefix,
        label_date=args.label_date,
        out_clusters=clusters_csv,
        out_labels=outdir / "cluster_labels.csv",
        out_changes=outdir / "label_changes.csv",
        outdir_label=str(outdir),
    )


def _write_sample_outputs(
    args: argparse.Namespace,
    outdir: Path,
    input_path: Path,
    metadata_rows: Sequence[Mapping[str, Any]],
    *,
    genotypes: Sequence[str] = (),
    clusters_csv: Path | None = None,
) -> None:
    """Run the audit stage over the files this run wrote (same code the workflow uses)."""
    by_genotype = outdir / "by_genotype"
    paths = {
        gt: by_genotype / hcv_workflow.safe_genotype_name(gt) for gt in genotypes
    }
    variants_csv = outdir / hcv_cluster_variants.VARIANTS_FILENAME
    resolved = outdir / "input.variants_resolved.fasta"
    hcv_cluster_stages.write_sample_outputs(
        outdir,
        input_fasta=input_path,
        genotypes_csv=outdir / "genotypes.csv",
        variants_csv=variants_csv if variants_csv.exists() else None,
        resolved_fasta=resolved if resolved.exists() else None,
        qc={gt: d / "prep.qc.csv" for gt, d in paths.items()},
        validation={gt: d / "prep.genotype_validation.csv" for gt, d in paths.items()},
        decisions={gt: d / "variant_decisions.csv" for gt, d in paths.items() if (d / "variant_decisions.csv").exists()},
        clusters_csv=clusters_csv,
        metadata_rows=metadata_rows,
        min_query_coverage=args.min_query_coverage,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        genotype_filter=args.genotype,
        reference_map=args.reference_map,
    )


def command_run(args: argparse.Namespace) -> int:
    outdir: Path = args.outdir
    input_path: Path = Path(args.input)
    genotype_fasta_dir = outdir / "genotype_fastas"
    genotype_csv = outdir / "genotypes.csv"
    metadata_rows = (
        hcv_cluster_metadata.load_metadata_csv(args.metadata) if args.metadata is not None else []
    )
    _validate_plot_fields(args, metadata_rows)
    _check_label_options(args)

    if args.dry_run:
        genotype_argv_display = [
            "assign_hcv_genotypes_from_fasta",
            "--input", str(input_path),
            "--output-csv", str(genotype_csv),
            "--outdir", str(genotype_fasta_dir),
            "--write-genotype-fastas",
        ]
        print(hcv_workflow.shell_join(genotype_argv_display))
        genotype_hint = args.genotype or ["<each detected pass genotype>"]
        for genotype in genotype_hint:
            genotype_dir = outdir / "by_genotype" / hcv_workflow.safe_genotype_name(genotype)
            prefix = genotype_dir / "prep"
            genotype_fasta = genotype_fasta_dir / f"hcv_{genotype}.fasta"
            clustering_fasta = Path(f"{prefix}.clustering.fasta")
            threshold = resolve_threshold(args.region, args.threshold)
            prep_cmd = [
                "hcv_cluster_prep", "prep-align",
                "--input", str(genotype_fasta),
                "--genotype", genotype,
                "--out-prefix", str(prefix),
            ]
            if args.reuse_alignments:
                cached_prefix = (
                    args.reuse_alignments
                    / "by_genotype"
                    / hcv_workflow.safe_genotype_name(genotype)
                    / "prep"
                )
                cached_alignment = Path(f"{cached_prefix}.aligned.fasta")
                # Only print a flag that would actually work; the real run applies the
                # same existence check, and prep-align errors on a missing cache file.
                if cached_alignment.exists():
                    prep_cmd.extend(["--cached-alignment", str(cached_alignment)])
            print(hcv_workflow.shell_join(prep_cmd))
            if args.distance in ("tn93", "both"):
                tn93_csv = genotype_dir / "tn93.csv"
                tn93_cmd = [
                    args.tn93,
                    *hcv_workflow.shell_split(args.extra_tn93_args),
                    "-t", str(threshold),
                    "-a", args.ambiguities,
                ]
                if args.ambiguity_fraction is not None:
                    tn93_cmd.extend(["-g", str(args.ambiguity_fraction)])
                if args.min_overlap is not None:
                    tn93_cmd.extend(["-l", str(args.min_overlap)])
                tn93_cmd.extend(["-o", str(tn93_csv), str(clustering_fasta)])
                print(hcv_workflow.shell_join(tn93_cmd))
            if args.distance in ("snp", "both"):
                snp_csv = genotype_dir / "snp.csv"
                print(hcv_workflow.shell_join(["write_snp_csv", str(snp_csv), str(clustering_fasta)]))
        return 0

    if args.reuse_alignments and not args.reuse_alignments.is_dir():
        raise RuntimeError(
            f"Alignment results directory does not exist: {args.reuse_alignments}"
        )

    outdir.mkdir(parents=True, exist_ok=True)
    if args.metadata is not None:
        hcv_cluster_metadata.write_metadata_csv(metadata_rows, outdir / "metadata.csv")
    suppress = _quiet_logging if not args.verbose else contextlib.nullcontext

    # --- Repeated FASTA IDs (e.g. mixed infections) ---
    _input_records, variant_plan = hcv_cluster_stages.resolve_ids(input_path, None, None)
    genotype_input = input_path
    if variant_plan.has_duplicates:
        genotype_input = outdir / "input.variants_resolved.fasta"
        hcv_cluster_prep.write_fasta_records(variant_plan.records, genotype_input)
        hcv_cluster_variants.write_variants_csv(variant_plan, outdir / hcv_cluster_variants.VARIANTS_FILENAME)
        repeats = sum(1 for v in variant_plan.variants if v.status == "dropped_identical_repeat")
        mixed = variant_plan.mixed_sample_ids()
        print(
            f"NOTE: repeated FASTA IDs found. {len(mixed)} sample(s) have distinct sequences "
            f"(renamed ID{hcv_cluster_variants.VARIANT_SEPARATOR}1, "
            f"ID{hcv_cluster_variants.VARIANT_SEPARATOR}2, ...); {repeats} identical repeat(s) dropped. "
            f"Same-genotype variants are reduced to one per sample by {args.duplicate_selection}. "
            f"Details: {hcv_cluster_variants.VARIANTS_FILENAME}"
        )

    # --- Genotyping ---
    print(f"Checking subtypes for: {input_path}")

    genotype_argv = [
        "--input", str(genotype_input),
        "--output-csv", str(genotype_csv),
        "--outdir", str(genotype_fasta_dir),
        "--panel-fasta", str(args.panel_fasta),
        "--minimap2", args.minimap2,
        "--preset", args.minimap2_preset,
        "--max-secondary", str(args.max_secondary),
        "--min-query-coverage", str(args.min_query_coverage),
        "--min-identity", str(args.min_identity),
        "--close-hit-fraction", str(args.close_hit_fraction),
        "--write-genotype-fastas",
    ]
    if args.extra_minimap2_args:
        genotype_argv.extend(["--extra-minimap2-args", args.extra_minimap2_args])
    if args.keep_paf:
        genotype_argv.append("--keep-paf")
    rc = assign_hcv_genotypes_from_fasta.main(genotype_argv)
    if rc:
        return rc

    genotype_rows = hcv_workflow.read_genotype_assignments(genotype_csv)
    genotypes = hcv_workflow.detected_pass_genotypes(genotype_rows)
    if args.genotype:
        requested = set(args.genotype)
        genotypes = [genotype for genotype in genotypes if genotype in requested]
    if not genotypes:
        _write_sample_outputs(args, outdir, input_path, metadata_rows)
        raise RuntimeError("No passing genotype assignments were available to process")

    # A detected subtype with no configured reference must not abort the whole run and
    # discard every other subtype's results; skip it with a warning, as hcv_workflow does.
    reference_catalog = hcv_cluster_prep.load_reference_catalog(args.reference_map)
    unsupported = [genotype for genotype in genotypes if genotype not in reference_catalog]
    if unsupported:
        unsupported_counts = {
            gt: sum(1 for r in genotype_rows if r.get("assigned_genotype") == gt)
            for gt in unsupported
        }
        detail = ", ".join(f"{gt} ({unsupported_counts[gt]})" for gt in sorted(unsupported))
        print(
            f"WARNING: no reference genome configured for subtype(s) {detail}; "
            "these sequences are skipped. Supply one with --reference-map to include them.",
            file=sys.stderr,
        )
        genotypes = [genotype for genotype in genotypes if genotype in reference_catalog]
    if not genotypes:
        _write_sample_outputs(args, outdir, input_path, metadata_rows)
        raise RuntimeError("No passing genotype assignments have a configured reference genome")

    gt_counts = {
        gt: sum(1 for r in genotype_rows if r.get("assigned_genotype") == gt)
        for gt in genotypes
    }
    gt_summary = ", ".join(f"{gt} ({gt_counts[gt]})" for gt in genotypes)
    print(f"--> {len(genotypes)} subtype(s) found: {gt_summary}")
    if args.snp_count_threshold is not None:
        print(
            f"WARNING: --snp-count-threshold {args.snp_count_threshold} links on raw SNP "
            "count. Pairs are compared over differing numbers of sites (N/gap positions "
            "are skipped), so the same count means different divergence for different "
            "pairs, and no HCV clustering threshold in the literature is defined this "
            "way. Prefer the p-distance default for anything reportable."
        )
    if args.threshold is None and not region_threshold_is_evidence_based(args.region):
        if region_threshold_is_explicit(args.region):
            print(
                f"NOTE: the {resolve_threshold(args.region, None):.3g} default for region "
                f"'{args.region}' is a conservative operational convention, not a "
                "region-validated cutoff. Sensitivity analysis is recommended."
            )
        else:
            print(
                f"WARNING: no HCV-specific clustering threshold evidence for region '{args.region}'; "
                f"using {FALLBACK_THRESHOLD:.3g} as a provisional fallback. "
                "See docs/threshold_rationale.md and consider passing --threshold explicitly."
            )
    print()

    cluster_tables: list[tuple[str, Path]] = []
    snp_cluster_tables: list[tuple[str, Path]] = []
    link_tables: list[tuple[str, Path]] = []
    snp_link_tables: list[tuple[str, Path]] = []
    for genotype in genotypes:
        print(f"Analysing subtype: {genotype}")

        genotype_dir = outdir / "by_genotype" / hcv_workflow.safe_genotype_name(genotype)
        genotype_dir.mkdir(parents=True, exist_ok=True)
        prefix = genotype_dir / "prep"
        genotype_fasta = genotype_fasta_dir / f"hcv_{genotype}.fasta"

        prep_argv = [
            "prep-align",
            "--input", str(genotype_fasta),
            "--genotype", genotype,
            "--out-prefix", str(prefix),
            "--refs-dir", str(args.refs_dir),
            "--min-coverage", str(args.min_coverage),
            "--region-strategy", args.region_strategy,
            "--mafft", args.mafft,
            "--threads", str(args.threads),
            "--genotype-validation", args.genotype_validation,
            "--genotype-validation-kmer-size", str(args.genotype_validation_kmer_size),
        ]
        if args.reference_map:
            prep_argv.extend(["--reference-map", str(args.reference_map)])
        if args.ncbi_email:
            prep_argv.extend(["--ncbi-email", args.ncbi_email])
        if args.force_download:
            prep_argv.append("--force-download")
        if args.verbose:
            prep_argv.append("--verbose")
        if args.region:
            prep_argv.extend(["--region", args.region])
        if args.keep_temp:
            prep_argv.append("--keep-temp")
        if args.reuse_alignments:
            cached_prefix = (
                args.reuse_alignments
                / "by_genotype"
                / hcv_workflow.safe_genotype_name(genotype)
                / "prep"
            )
            cached_alignment = Path(f"{cached_prefix}.aligned.fasta")
            if cached_alignment.exists():
                prep_argv.extend(["--cached-alignment", str(cached_alignment)])
                reusable, total = hcv_cluster_prep.count_reusable_sequences(
                    genotype_fasta, cached_alignment
                )
                print(
                    f"--> Reusing {reusable}/{total} cached alignments from "
                    f"{cached_alignment}; {total - reusable} to align"
                )
            else:
                print(f"--> No cached alignment for subtype {genotype}; aligning all sequences")
        with suppress():
            rc = hcv_cluster_prep.main(prep_argv)
        if rc:
            return rc

        qc_path = Path(f"{prefix}.qc.csv")
        validation_path = Path(f"{prefix}.genotype_validation.csv")
        if validation_path.exists():
            print(
                f"--> {len(_read_csv_rows(validation_path))} sequence(s) excluded: their k-mer "
                f"profile conflicts with subtype {genotype} (see {validation_path.name})"
            )
        if qc_path.exists():
            qc_st = _read_qc_stats(qc_path)
            pct = f"{qc_st['passed'] / qc_st['total'] * 100:.0f}%" if qc_st["total"] else "n/a"
            print(f"--> {qc_st['total']} sequences, {qc_st['passed']} passed QC ({pct})")

        clustering_fasta = Path(f"{prefix}.clustering.fasta")
        if variant_plan.has_duplicates:
            dropped = hcv_cluster_stages.select_variants(
                variant_plan, genotype_rows, genotype, qc_path, args.duplicate_selection,
                clustering_fasta, out_decisions=genotype_dir / "variant_decisions.csv",
            )
            if dropped:
                print(
                    f"--> {len(dropped)} same-sample variant(s) dropped in favour of a "
                    f"better sequence ({args.duplicate_selection})"
                )

        if not _fasta_has_records(clustering_fasta):
            print(f"--> No sequences left to cluster for subtype {genotype}; skipping\n")
            continue

        if args.distance in ("tn93", "both"):
            threshold = resolve_threshold(args.region, args.threshold, distance="tn93")
            tn93_csv = genotype_dir / "tn93.csv"
            links_csv = genotype_dir / "links.csv"
            clusters_csv = genotype_dir / "clusters.csv"
            hcv_cluster_stages.cluster_sequences(
                "tn93",
                clustering_fasta,
                out_distances=tn93_csv,
                out_links=links_csv,
                out_clusters=clusters_csv,
                threshold=threshold,
                tn93=args.tn93,
                ambiguities=args.ambiguities,
                ambiguity_fraction=args.ambiguity_fraction,
                min_overlap=args.min_overlap,
                extra_tn93_args=args.extra_tn93_args,
                quiet=not args.verbose,
            )
            cluster_tables.append((genotype, clusters_csv))
            link_tables.append((genotype, links_csv))
            cl_st = _read_cluster_stats(clusters_csv)
            metric = "TN93" if args.distance == "both" else "TN93"
            print(
                f"--> {metric} (threshold {threshold:.4g}): "
                f"{cl_st['n_multi']} cluster(s), {cl_st['singletons']} singleton(s), "
                f"largest: {cl_st['largest']}"
            )
        if args.distance in ("snp", "both"):
            snp_threshold = resolve_threshold(args.region, args.threshold, distance="snp")
            snp_csv = genotype_dir / "snp.csv"
            if args.distance == "snp":
                snp_links_csv = genotype_dir / "links.csv"
                snp_clusters_csv = genotype_dir / "clusters.csv"
            else:
                snp_links_csv = genotype_dir / "snp_links.csv"
                snp_clusters_csv = genotype_dir / "snp_clusters.csv"
            snp_links = hcv_cluster_stages.cluster_sequences(
                "snp",
                clustering_fasta,
                out_distances=snp_csv,
                out_links=snp_links_csv,
                out_clusters=snp_clusters_csv,
                threshold=snp_threshold,
                snp_count_threshold=args.snp_count_threshold,
            )
            if args.distance == "snp":
                cluster_tables.append((genotype, snp_clusters_csv))
                link_tables.append((genotype, snp_links_csv))
            else:
                snp_cluster_tables.append((genotype, snp_clusters_csv))
                snp_link_tables.append((genotype, snp_links_csv))
            snp_st = _read_cluster_stats(snp_clusters_csv)
            if args.snp_count_threshold is not None:
                snp_threshold_label = f"<={args.snp_count_threshold} SNPs"
            else:
                snp_threshold_label = f"threshold {snp_threshold:.4g}"
            print(
                f"--> SNP   ({snp_threshold_label}): "
                f"{snp_st['n_multi']} cluster(s), {snp_st['singletons']} singleton(s), "
                f"largest: {snp_st['largest']}"
            )
            print(f"          {_snp_count_summary(snp_links)}")
        print()

    primary_clusters_csv = outdir / "clusters.csv"
    primary_links_csv = outdir / "links.csv"
    variants_csv = outdir / hcv_cluster_variants.VARIANTS_FILENAME
    hcv_cluster_stages.merge_tables(
        dict(cluster_tables), dict(link_tables), primary_clusters_csv, primary_links_csv, variants_csv
    )
    if args.distance == "both":
        hcv_cluster_stages.merge_tables(
            dict(snp_cluster_tables), dict(snp_link_tables),
            outdir / "clusters.snp.csv", outdir / "links.snp.csv", variants_csv,
        )
    for clusters_path in [path for _, path in cluster_tables + snp_cluster_tables]:
        hcv_cluster_stages.annotate_clusters(clusters_path, variants_csv)

    if args.registry is not None:
        _apply_cluster_labels(args, outdir, primary_clusters_csv)

    _write_sample_outputs(
        args, outdir, input_path, metadata_rows,
        genotypes=genotypes, clusters_csv=primary_clusters_csv,
    )

    _render_requested_plots(
        args,
        outdir,
        metadata_rows,
        cluster_tables,
        link_tables,
        snp_cluster_tables,
        snp_link_tables,
    )

    print(f"Done. Results written to: {outdir}")
    return 0


def command_refs(args: argparse.Namespace) -> int:
    suppress = _quiet_logging if not args.verbose else contextlib.nullcontext
    label = f"genotype {args.genotype}" if args.genotype != "all" else "all genotypes"
    print(f"Caching reference sequences for: {label}")
    prep_refs_argv = [
        "prep-refs",
        "--genotype", args.genotype,
        "--refs-dir", str(args.refs_dir),
    ]
    if args.reference_map:
        prep_refs_argv.extend(["--reference-map", str(args.reference_map)])
    if args.ncbi_email:
        prep_refs_argv.extend(["--ncbi-email", args.ncbi_email])
    if args.force_download:
        prep_refs_argv.append("--force-download")
    if args.verbose:
        prep_refs_argv.append("-v")
    with suppress():
        rc = hcv_cluster_prep.main(prep_refs_argv)
    if rc == 0:
        print(f"--> Done (stored in: {args.refs_dir})")
    return rc


def build_parser(show_advanced: bool = False) -> argparse.ArgumentParser:
    adv = lambda h: h if show_advanced else argparse.SUPPRESS  # noqa: E731

    parser = argparse.ArgumentParser(
        prog="hcv-cluster",
        description="HCV genomic clustering with genotype-aware thresholds.",
        epilog="Use --help-advanced to show all tuning options.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command")

    run_parser = subparsers.add_parser(
        "run",
        help="Run the full HCV clustering pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    run_parser.add_argument(
        "-i", "--input",
        metavar="PATH",
        required=True,
        type=Path,
        help="Input multi-FASTA of HCV consensus sequences",
    )
    run_parser.add_argument(
        "-o", "--outdir",
        metavar="PATH",
        default=Path("results"),
        type=Path,
        help="Output directory (default: results)",
    )
    run_parser.add_argument(
        "--metadata",
        metavar="CSV",
        type=Path,
        default=None,
        help=(
            "Optional sample metadata CSV keyed by sample_id. Normalized metadata is "
            "saved as metadata.csv in the results directory."
        ),
    )
    run_parser.add_argument(
        "-t", "--threshold",
        metavar="FLOAT",
        type=float,
        default=None,
        help=(
            "Maximum distance for cluster linking. Default: region-dependent "
            "(0.03 for e1-e2, the default region, and e1-e2-nohvr1; 0.045 for "
            "core-e2 with HVR1 included; 0.015 for ns5b; 0.03 for any other region — "
            "see docs/threshold_rationale.md)."
        ),
    )
    run_parser.add_argument(
        "-d", "--distance",
        metavar="CHOICE",
        choices=("tn93", "snp", "both"),
        default="tn93",
        help="Distance metric(s). tn93: TN93 evolutionary distance (default). snp: normalized SNP fraction. both: compute both; TN93 drives primary outputs.",
    )
    run_parser.add_argument(
        "--ambiguities",
        metavar="STR",
        default="average",
        help=(
            "TN93 ambiguity handling: average (default), resolve, skip, gapmm. "
            "N is always masked to a gap before clustering regardless of this "
            "setting (neutralizes low-depth-masked positions); average then "
            "gives real IUPAC ambiguity codes (R, Y, etc.) proportional "
            "treatment. See docs/threshold_rationale.md."
        ),
    )
    run_parser.add_argument(
        "-r", "--region",
        metavar="EXPR",
        default="e1-e2",
        help=(
            "Reference-anchored region expression (default: e1-e2, i.e. full "
            "E1-through-E2 with HVR1 included — see docs/threshold_rationale.md). "
            "Individual regions: core e1 e2 p7 ns2 ns3 ns4a ns4b ns5a ns5b. "
            "Append -nohvr1 to e1-e2 or core-e2 to mask the first 81 nt of E2. "
            "Named presets: structural (=core-e2-nohvr1), envelope (=e1-e2), "
            "nonstructural (=ns2-ns5b), cds (whole coding sequence). "
            "Range syntax: first-last selects all regions from first through last "
            "(e.g. e1-e2, ns3-ns5b). "
            "Union syntax: a+b includes both (e.g. core-e2-nohvr1+ns3, e1-e2+ns5a)."
        ),
    )
    run_parser.add_argument(
        "-T", "--threads",
        metavar="INT",
        type=int,
        default=1,
        help="MAFFT alignment threads (default: 1)",
    )
    run_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose/debug logging",
    )
    run_parser.add_argument(
        "--help-advanced",
        action="store_true",
        help="Show all advanced tuning options and exit",
    )
    run_parser.add_argument(
        "--genotype",
        action="append",
        metavar="GENOTYPE",
        help=adv("Restrict processing to this genotype after auto-detection. May be repeated."),
    )
    run_parser.add_argument(
        "--registry",
        metavar="DIR",
        type=Path,
        default=None,
        help=(
            "Persistent cluster-label registry directory. Clusters of two or more sequences "
            "keep their label (e.g. HCV-2608-001) across runs; singletons are never labelled. "
            "Read-only unless --update-registry is given."
        ),
    )
    run_parser.add_argument(
        "--update-registry",
        action="store_true",
        help=(
            "Save this run's labels to --registry: name new clusters, record merges and "
            "splits. Leave off for exploratory runs (different thresholds, test data)."
        ),
    )
    run_parser.add_argument(
        "--label-prefix",
        metavar="STR",
        default=hcv_cluster_labels.DEFAULT_PREFIX,
        help=adv("Prefix for newly created cluster labels (default: HCV)"),
    )
    run_parser.add_argument(
        "--label-date",
        metavar="YYYY-MM-DD",
        default=None,
        help=adv("Naming date for new cluster labels (default: today); sets the YYMM part"),
    )
    run_parser.add_argument(
        "--duplicate-selection",
        choices=hcv_cluster_variants.SELECTION_MODES,
        default="region-coverage",
        help=(
            "When one FASTA ID has several distinct sequences with the same genotype, keep "
            "only the one with the best clustering-region coverage (ties: most complete), "
            "or only the most complete. Variants of different genotypes are all kept. "
            "(default: region-coverage)"
        ),
    )
    run_parser.add_argument(
        "--min-coverage",
        metavar="FLOAT",
        type=float,
        default=0.7,
        help=adv("Minimum non-gap non-N selected-region coverage fraction (default: 0.7)"),
    )
    run_parser.add_argument(
        "--min-identity",
        metavar="FLOAT",
        type=float,
        default=0.75,
        help=adv(
            "Minimum minimap2 gap-compressed identity for genotype pass "
            "(default: 0.75)"
        ),
    )
    run_parser.add_argument(
        "--min-query-coverage",
        metavar="FLOAT",
        type=float,
        default=0.50,
        help=adv(
            "Minimum fraction of callable (non-N) query bases aligned for genotype pass "
            "after combining compatible split alignments (default: 0.50)"
        ),
    )
    run_parser.add_argument(
        "--panel-fasta",
        metavar="PATH",
        type=Path,
        default=assign_hcv_genotypes_from_fasta.DEFAULT_PANEL,
        help=adv("Reference panel FASTA for genotyping"),
    )
    run_parser.add_argument(
        "--minimap2",
        metavar="STR",
        default="minimap2",
        help=adv("minimap2 executable path/name"),
    )
    run_parser.add_argument(
        "--minimap2-preset",
        metavar="STR",
        default="asm20",
        help=adv("minimap2 preset (default: asm20)"),
    )
    run_parser.add_argument(
        "--extra-minimap2-args",
        metavar="STR",
        default="",
        help=adv("Extra minimap2 arguments as a shell string"),
    )
    run_parser.add_argument(
        "--max-secondary",
        metavar="INT",
        type=int,
        default=50,
        help=adv("Maximum secondary hits for minimap2 (default: 50)"),
    )
    run_parser.add_argument(
        "--close-hit-fraction",
        metavar="FLOAT",
        type=float,
        default=0.98,
        help=adv("Score fraction threshold for reporting close genotype hits (default: 0.98)"),
    )
    run_parser.add_argument(
        "--keep-paf",
        action="store_true",
        help=adv("Keep raw minimap2 PAF file"),
    )
    run_parser.add_argument(
        "--refs-dir",
        metavar="PATH",
        type=Path,
        default=hcv_cluster_prep.DEFAULT_REFS_DIR,
        help=adv("Directory of reference GenBank/FASTA files (default: the set bundled with the repository)"),
    )
    run_parser.add_argument(
        "--reference-map",
        metavar="PATH",
        type=Path,
        default=None,
        help=adv("JSON file extending or overriding genotype-to-accession mappings"),
    )
    run_parser.add_argument(
        "--ncbi-email",
        metavar="STR",
        default=os.environ.get("NCBI_EMAIL"),
        help=adv("Email for NCBI E-utilities"),
    )
    run_parser.add_argument(
        "--force-download",
        action="store_true",
        help=adv("Re-download cached NCBI reference records"),
    )
    run_parser.add_argument(
        "--region-strategy",
        metavar="CHOICE",
        choices=("fixed", "max-usable"),
        default="fixed",
        help=adv("Region selection strategy (default: fixed)"),
    )
    run_parser.add_argument(
        "--mafft",
        metavar="STR",
        default="mafft",
        help=adv("MAFFT executable path/name"),
    )
    run_parser.add_argument(
        "--genotype-validation",
        metavar="CHOICE",
        choices=("auto", "headers", "kmer", "none"),
        default="auto",
        help=adv("Genotype validation mode (default: auto)"),
    )
    run_parser.add_argument(
        "--genotype-validation-kmer-size",
        metavar="INT",
        type=int,
        default=15,
        help=adv("k-mer size for reference-anchor genotype validation (default: 15)"),
    )
    run_parser.add_argument(
        "--keep-temp",
        action="store_true",
        help=adv("Keep temporary MAFFT files"),
    )
    run_parser.add_argument(
        "--snp-count-threshold",
        metavar="N",
        type=int,
        default=None,
        help=(
            "Link SNP pairs on absolute differing-site count (e.g. 30) instead of "
            "p-distance. The count is reported either way; see docs/threshold_rationale.md "
            "for why an absolute count is not comparable across pairs"
        ),
    )
    run_parser.add_argument(
        "--reuse-alignments",
        metavar="RESULTS_DIR",
        type=Path,
        default=None,
        help=adv(
            "Reuse unchanged sequences from per-genotype prep.aligned.fasta files "
            "under a previous results directory; align only new or changed sequences"
        ),
    )
    run_parser.add_argument(
        "--tn93",
        metavar="STR",
        default="tn93",
        help=adv("tn93 executable path/name"),
    )
    run_parser.add_argument(
        "--ambiguity-fraction",
        metavar="FLOAT",
        type=float,
        default=None,
        help=adv("TN93 maximum resolvable ambiguity fraction (-g flag)"),
    )
    run_parser.add_argument(
        "--min-overlap",
        metavar="INT",
        type=int,
        default=None,
        help=adv("Minimum pairwise overlap passed to TN93 (-l flag)"),
    )
    run_parser.add_argument(
        "--extra-tn93-args",
        metavar="STR",
        default="",
        help=adv("Extra TN93 arguments as a shell string"),
    )
    run_parser.add_argument(
        "--dry-run",
        action="store_true",
        help=adv("Print commands without executing them"),
    )
    run_parser.add_argument(
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
    run_parser.add_argument(
        "--plot-scope",
        metavar="CHOICE",
        choices=("per-genotype", "combined", "both"),
        default="per-genotype",
        help=adv(
            "Plot each genotype separately (default), all genotypes in one presentation-only "
            "network, or both. Combined plots never add cross-genotype links."
        ),
    )
    run_parser.add_argument(
        "--plot-combined-layout",
        metavar="CHOICE",
        choices=("packed", "by-genotype"),
        default="packed",
        help=adv(
            "Combined-network layout: compact connected-component packing (default), or "
            "spatially group packed components by genotype."
        ),
    )
    run_parser.add_argument(
        "--plot-node-spacing",
        metavar="CHOICE",
        choices=("compact", "normal", "expanded"),
        default="normal",
        help=adv(
            "Minimum collision-aware node separation: compact, normal (default), or expanded"
        ),
    )
    run_parser.add_argument(
        "--plot-color-by",
        metavar="FIELD",
        default=None,
        help=adv("Metadata field encoded as node colour"),
    )
    run_parser.add_argument(
        "--plot-symbol-by",
        metavar="FIELD",
        default=None,
        help=adv("Low-cardinality metadata field encoded as node shape"),
    )
    run_parser.add_argument(
        "--plot-size-by",
        metavar="FIELD",
        default=None,
        help=adv(
            "Numeric or intrinsically ordered metadata field encoded as node size; use "
            "--plot-size-order for other categorical fields"
        ),
    )
    run_parser.add_argument(
        "--plot-size-order",
        action="append",
        metavar="VALUE",
        default=None,
        help=adv(
            "Category order from smallest to largest for --plot-size-by; repeat once per value"
        ),
    )
    run_parser.add_argument(
        "--plot-outline-by",
        metavar="FIELD",
        default=None,
        help=adv("Low-cardinality metadata field encoded as node outline colour"),
    )
    run_parser.add_argument(
        "--plot-center-by",
        metavar="FIELD",
        default=None,
        help=adv("Low-cardinality metadata field encoded as a small centre mark"),
    )
    run_parser.add_argument(
        "--plot-small-multiple-field",
        action="append",
        metavar="FIELD",
        default=None,
        help=adv(
            "Render the same network in a colour panel for this metadata field; repeat for "
            "up to four fields. Cannot be combined with composite node encodings."
        ),
    )
    run_parser.add_argument(
        "--plot-hover-field",
        action="append",
        metavar="FIELD",
        default=None,
        help=adv(
            "Metadata field included in interactive hover text; repeat for multiple fields. "
            "Defaults to all supplied metadata; sample_id is always included."
        ),
    )
    run_parser.add_argument(
        "--plot-hide-singletons",
        action="store_true",
        help=adv("Omit singleton (unclustered) sequences from --plot-network output"),
    )
    run_parser.set_defaults(func=command_run)

    refs_parser = subparsers.add_parser(
        "refs",
        help="Download and cache reference records",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    refs_parser.add_argument(
        "--genotype",
        default="all",
        help="Genotype to prepare, or 'all' (default: all)",
    )
    refs_parser.add_argument(
        "--refs-dir",
        metavar="PATH",
        type=Path,
        default=hcv_cluster_prep.DEFAULT_REFS_DIR,
        help="Reference directory (default: the set bundled with the repository)",
    )
    refs_parser.add_argument(
        "--reference-map",
        metavar="PATH",
        type=Path,
        default=None,
    )
    refs_parser.add_argument(
        "--ncbi-email",
        metavar="STR",
        default=os.environ.get("NCBI_EMAIL"),
    )
    refs_parser.add_argument(
        "--force-download",
        action="store_true",
    )
    refs_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
    )
    refs_parser.set_defaults(func=command_refs)

    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    show_advanced = "--help-advanced" in argv

    parser = build_parser(show_advanced=show_advanced)

    if not argv or argv[0] in {"-h", "--help", "--help-advanced"}:
        parser.print_help()
        return 0

    if show_advanced:
        # Replace --help-advanced with --help so argparse handles it natively
        # for the active subcommand without tripping on required args.
        argv = [("--help" if a == "--help-advanced" else a) for a in argv]
        try:
            parser.parse_args(argv)
        except SystemExit as exc:
            return int(exc.code) if exc.code is not None else 0
        return 0

    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.error("a subcommand is required")
    try:
        return args.func(args)
    except (RuntimeError, ValueError, hcv_cluster_prep.HcvPrepError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
