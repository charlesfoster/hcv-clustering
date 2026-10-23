#!/usr/bin/env python3
"""File-in/file-out pipeline stages.

Every step that is not an external tool (minimap2, MAFFT, TN93) lives here as a plain
function over files, and as a subcommand of ``hcv-stage``. ``hcv_cluster.py run`` calls
these same functions in one process, and the Nextflow workflow calls them as separate
processes, so the two cannot drift apart. No scientific logic belongs in the workflow
definition: thresholds, region handling and sample bookkeeping are all decided here.

Stage order (``genotype`` and ``align`` are the existing assign_hcv_genotypes_from_fasta.py
and ``hcv_cluster_prep.py prep-align`` commands)::

    resolve-ids -> genotype -> align (per subtype) -> select-variants (per subtype)
        -> cluster (per subtype and method) -> merge -> label -> audit
"""

from __future__ import annotations

import argparse
import csv
import datetime
import shutil
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import hcv_cluster_audit
import hcv_cluster_labels
import hcv_cluster_metadata
import hcv_cluster_methods
import hcv_cluster_prep
import hcv_cluster_variants
import hcv_workflow


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_labelled_paths(values: Sequence[str] | None) -> dict[str, Path]:
    """Parse repeated ``GENOTYPE=PATH`` arguments."""
    parsed: dict[str, Path] = {}
    for value in values or ():
        genotype, separator, path = value.partition("=")
        if not separator or not genotype or not path:
            raise ValueError(f"Expected GENOTYPE=PATH, got {value!r}")
        parsed[genotype] = Path(path)
    return parsed


# --- resolve-ids ------------------------------------------------------------------------


def resolve_ids(
    input_path: Path, out_fasta: Path | None, out_variants: Path | None
) -> tuple[list[hcv_cluster_prep.FastaRecord], hcv_cluster_variants.VariantPlan]:
    """Drop exact repeats and rename distinct same-ID sequences (mixed infections)."""
    input_records = hcv_cluster_variants.read_input_records(input_path)
    plan = hcv_cluster_variants.resolve_duplicate_ids(input_records)
    if out_fasta is not None:
        hcv_cluster_prep.write_fasta_records(plan.records, out_fasta)
    if out_variants is not None:
        hcv_cluster_variants.write_variants_csv(plan, out_variants)
    return input_records, plan


# --- plan-genotypes ---------------------------------------------------------------------


def plan_genotypes(
    genotype_rows: Sequence[Mapping[str, str]],
    genotype_filter: Sequence[str] | None,
    catalog: Mapping[str, Any],
) -> tuple[list[str], dict[str, int]]:
    """Subtypes to analyse, and the sequence count of passing subtypes with no reference.

    A subtype without a configured reference is skipped, never fatal: one rare subtype
    must not discard every other subtype's results.
    """
    genotypes = hcv_workflow.detected_pass_genotypes(genotype_rows)
    if genotype_filter:
        genotypes = [g for g in genotypes if g in set(genotype_filter)]
    unsupported = {
        g: sum(1 for row in genotype_rows if row.get("assigned_genotype") == g)
        for g in genotypes
        if g not in catalog
    }
    return [g for g in genotypes if g not in unsupported], unsupported


def warn_unsupported(unsupported: Mapping[str, int]) -> None:
    if unsupported:
        detail = ", ".join(f"{gt} ({n})" for gt, n in sorted(unsupported.items()))
        print(
            f"WARNING: no reference genome configured for subtype(s) {detail}; "
            "these sequences are skipped. Supply one with --reference-map to include them.",
            file=sys.stderr,
        )


def filter_genotype_fastas(
    runnable: Sequence[str], fasta_dir: Path, skipped_dir: Path
) -> None:
    """Move per-subtype FASTAs that will not be analysed out of the way."""
    skipped_dir.mkdir(parents=True, exist_ok=True)
    for fasta in sorted(fasta_dir.glob("hcv_*.fasta")):
        if fasta.stem.removeprefix("hcv_") not in runnable:
            shutil.move(str(fasta), skipped_dir / fasta.name)


# --- select-variants --------------------------------------------------------------------


def select_variants(
    plan: hcv_cluster_variants.VariantPlan,
    genotype_rows: Sequence[Mapping[str, str]],
    genotype: str,
    qc_path: Path,
    mode: str,
    clustering_fasta: Path,
    out_fasta: Path | None = None,
    out_decisions: Path | None = None,
) -> list[hcv_cluster_variants.VariantRecord]:
    """Keep one usable variant per sample within ``genotype``; return those dropped."""
    hcv_cluster_variants.record_genotypes(plan, genotype_rows)
    dropped = hcv_cluster_variants.select_same_genotype(plan, genotype, qc_path, mode)
    target = out_fasta or clustering_fasta
    if target != clustering_fasta:
        shutil.copyfile(clustering_fasta, target)
    if dropped:
        hcv_cluster_variants.remove_from_fasta(target, {variant.sequence_id for variant in dropped})
    if out_decisions is not None:
        hcv_cluster_variants.write_decisions(plan, genotype, out_decisions)
    return dropped


# --- cluster ----------------------------------------------------------------------------


def cluster_sequences(
    method: str,
    fasta: Path,
    *,
    out_distances: Path,
    out_links: Path,
    out_clusters: Path,
    threshold: float,
    snp_count_threshold: int | None = None,
    tn93: str = "tn93",
    ambiguities: str = "average",
    ambiguity_fraction: float | None = None,
    min_overlap: int | None = None,
    extra_tn93_args: str = "",
    quiet: bool = True,
) -> list[dict[str, str | int | float]]:
    """Distances, linking and connected components for one subtype with one method.

    Returns the SNP links (empty for TN93) so callers can describe them in whole SNPs.
    """
    links: list[dict[str, str | int | float]] = []
    if not any(line.startswith(">") for line in fasta.read_text().splitlines()):
        # Nothing survived QC for this subtype: leave empty tables rather than failing.
        for path, header in (
            (out_distances, "source,target,distance\n"),
            (out_links, "source,target,distance\n"),
            (out_clusters, "sample_id,cluster_id,cluster_size\n"),
        ):
            path.write_text(header)
        return links
    if method == "tn93":
        rc = hcv_workflow.run_tn93(
            argparse.Namespace(
                input=fasta,
                output=out_distances,
                tn93=tn93,
                tn93_threshold=threshold,
                ambiguities=ambiguities,
                ambiguity_fraction=ambiguity_fraction,
                min_overlap=min_overlap,
                quiet=quiet,
                extra_tn93_args=extra_tn93_args,
                dry_run=False,
            )
        )
        if rc:
            raise RuntimeError(f"tn93 failed with exit code {rc}")
        hcv_workflow.command_link(
            argparse.Namespace(
                distances=out_distances, output=out_links, threshold=threshold, include_self=False
            )
        )
    elif method == "snp":
        detailed = hcv_cluster_methods.compute_snp_distances_detailed(fasta)
        hcv_cluster_methods.write_snp_csv(out_distances, fasta, detailed)
        links = hcv_cluster_methods.write_snp_links_csv(
            out_links, detailed, threshold=threshold, snp_count_threshold=snp_count_threshold
        )
    else:
        raise ValueError(f"Unknown clustering method: {method}")
    hcv_workflow.command_cluster(
        argparse.Namespace(links=out_links, output=out_clusters, nodes_fasta=fasta, nodes_file=None)
    )
    return links


# --- merge ------------------------------------------------------------------------------


def merge_tables(
    cluster_tables: Mapping[str, Path],
    link_tables: Mapping[str, Path],
    out_clusters: Path,
    out_links: Path,
    variants_csv: Path | None = None,
) -> None:
    """Merge per-subtype tables into genotype-qualified, run-wide ones."""
    hcv_workflow.merge_cluster_tables(sorted(cluster_tables.items()), out_clusters)
    hcv_workflow.merge_link_tables(sorted(link_tables.items()), out_links)
    annotate_clusters(out_clusters, variants_csv)


def annotate_clusters(clusters_csv: Path, variants_csv: Path | None) -> None:
    """Add base_sample_id / mixed_infection (identity values when nothing was renamed)."""
    plan = (
        hcv_cluster_variants.read_variants_csv(variants_csv)
        if variants_csv is not None
        else hcv_cluster_variants.VariantPlan(records=[])
    )
    hcv_cluster_variants.annotate_cluster_table(clusters_csv, plan)


# --- label ------------------------------------------------------------------------------


def label_definition(
    metric: str,
    region: str,
    threshold: float | None,
    region_strategy: str = "fixed",
    ambiguities: str = "average",
    ambiguity_fraction: float | None = None,
    min_overlap: int | None = None,
    extra_tn93_args: str = "",
    snp_count_threshold: int | None = None,
) -> dict[str, object]:
    """What a persistent label is valid for: changing any of this changes who clusters."""
    definition: dict[str, object] = {
        "metric": metric,
        "region": hcv_cluster_methods.canonical_threshold_region(region),
        "threshold": hcv_cluster_methods.resolve_threshold(region, threshold, distance=metric),
        "region_strategy": region_strategy,
    }
    if metric == "snp" and snp_count_threshold is not None:
        definition["snp_count_threshold"] = snp_count_threshold
    if metric == "tn93":
        definition.update(
            ambiguities=ambiguities,
            ambiguity_fraction=ambiguity_fraction,
            min_overlap=min_overlap,
            extra_tn93_args=extra_tn93_args,
        )
    return definition


def label_clusters(
    clusters_csv: Path,
    *,
    registry: Path,
    update_registry: bool,
    track: str,
    definition: Mapping[str, object],
    prefix: str,
    label_date: str | None,
    out_clusters: Path,
    out_labels: Path,
    out_changes: Path,
    registry_out: Path | None = None,
    outdir_label: str = "",
) -> hcv_cluster_labels.LabelResult:
    """Add the persistent ``cluster_label`` column and write label reports."""
    rows = read_csv_rows(clusters_csv)
    with clusters_csv.open(newline="", encoding="utf-8") as handle:
        fieldnames = list(csv.DictReader(handle).fieldnames or ())
    base = hcv_cluster_variants.BASE_SAMPLE_FIELD
    partition = {(row["genotype"], row.get(base) or row["sample_id"]): row["cluster_id"] for row in rows}
    now = datetime.datetime.now()
    result, table, _registry = hcv_cluster_labels.label_run(
        registry,
        track,
        definition,
        partition,
        run_id=now.strftime("%Y%m%dT%H%M%S"),
        run_date=label_date or now.strftime("%Y-%m-%d"),
        outdir=outdir_label,
        commit=update_registry,
        prefix=prefix,
        output_directory=registry_out,
    )
    if "cluster_label" not in fieldnames:
        fieldnames.insert(fieldnames.index("cluster_id") + 1, "cluster_label")
    for row in rows:
        row["cluster_label"] = result.labels[(row["genotype"], row.get(base) or row["sample_id"])]
    with out_clusters.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    label_fields = ["label", "track", "genotype", "created_date", "split_from",
                    "absorbed_labels", "n_members_this_run", "n_members_total"]
    with out_labels.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=label_fields, restval="")
        writer.writeheader()
        writer.writerows(table)
    event_fields = ["run_id", "date", "track", "event", "label", "related_label", "detail"]
    with out_changes.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=event_fields)
        writer.writeheader()
        writer.writerows(result.events)

    counts = {kind: sum(1 for e in result.events if e["event"] == kind) for kind in ("created", "split", "merged")}
    labelled = [member for member, label in result.labels.items() if label]
    summary = (
        f"Cluster labels: {len({partition[member] for member in labelled})} labelled cluster(s) "
        f"covering {len(labelled)} sequence(s); {counts['created'] + counts['split']} new, "
        f"{counts['merged']} merged"
    )
    if update_registry:
        print(f"{summary}. Registry updated: {registry_out or registry}")
    else:
        print(
            f"{summary}. Registry NOT updated (read-only); new clusters are shown as "
            f"{hcv_cluster_labels.UNASSIGNED}. Re-run with --update-registry to name them."
        )
    return result


# --- audit ------------------------------------------------------------------------------


def derive_excluded_subtypes(
    genotype_rows: Sequence[Mapping[str, str]],
    processed: Sequence[str],
    genotype_filter: Sequence[str] | None,
    catalog: Mapping[str, Any],
) -> dict[str, str]:
    """Subtypes that passed genotyping but were not analysed, with the reason."""
    excluded: dict[str, str] = {}
    requested = set(genotype_filter or ())
    for genotype in hcv_workflow.detected_pass_genotypes(genotype_rows):
        if genotype in processed:
            continue
        if requested and genotype not in requested:
            excluded[genotype] = f"Subtype {genotype} excluded by --genotype ({', '.join(genotype_filter or ())})"
        elif genotype not in catalog:
            excluded[genotype] = (
                f"No reference genome configured for subtype {genotype}; supply one with --reference-map"
            )
        else:
            excluded[genotype] = f"Subtype {genotype} was not analysed"
    return excluded


def write_sample_outputs(
    outdir: Path,
    *,
    input_fasta: Path,
    genotypes_csv: Path,
    variants_csv: Path | None = None,
    resolved_fasta: Path | None = None,
    qc: Mapping[str, Path] | None = None,
    validation: Mapping[str, Path] | None = None,
    decisions: Mapping[str, Path] | None = None,
    clusters_csv: Path | None = None,
    metadata_rows: Sequence[Mapping[str, Any]] = (),
    min_query_coverage: float,
    min_identity: float,
    min_coverage: float,
    genotype_filter: Sequence[str] | None = None,
    reference_map: Path | None = None,
) -> list[hcv_cluster_audit.SampleAudit]:
    """Write the sample audit, annotated metadata tables and labelled survivor FASTA.

    Reads everything from files, so it can run as its own pipeline step. Also used when a
    run stops early (nothing passed genotyping), when the explanation matters most.
    """
    qc = qc or {}
    input_records = hcv_cluster_variants.read_input_records(input_fasta)
    resolved = (
        hcv_cluster_prep.read_fasta_records(resolved_fasta)
        if resolved_fasta is not None and resolved_fasta.exists()
        else input_records
    )
    plan = (
        hcv_cluster_variants.read_variants_csv(variants_csv, resolved)
        if variants_csv is not None
        else hcv_cluster_variants.VariantPlan(records=list(resolved))
    )
    genotype_rows = hcv_workflow.read_genotype_assignments(genotypes_csv)
    hcv_cluster_variants.record_genotypes(plan, genotype_rows)
    for path in (decisions or {}).values():
        hcv_cluster_variants.apply_decisions(plan, path)

    catalog = hcv_cluster_prep.load_reference_catalog(reference_map)
    excluded = derive_excluded_subtypes(genotype_rows, list(qc), genotype_filter, catalog)
    validation_failures = {
        genotype: {row["sample_id"]: row["reason"] for row in read_csv_rows(path)}
        for genotype, path in (validation or {}).items()
        if path.exists()
    }
    audits = hcv_cluster_audit.build_audit(
        input_records,
        plan,
        genotype_rows,
        min_query_coverage=min_query_coverage,
        min_identity=min_identity,
        min_coverage=min_coverage,
        excluded_subtypes=excluded,
        qc_by_genotype={genotype: hcv_cluster_audit.read_qc_rows(path) for genotype, path in qc.items()},
        validation_failures=validation_failures,
        cluster_rows=read_csv_rows(clusters_csv) if clusters_csv is not None and clusters_csv.exists() else (),
    )
    outdir.mkdir(parents=True, exist_ok=True)
    hcv_cluster_audit.write_audit_csv(audits, outdir / hcv_cluster_audit.AUDIT_FILENAME)
    hcv_cluster_audit.write_metadata_tables(
        audits,
        metadata_rows,
        outdir / hcv_cluster_audit.METADATA_ALL_FILENAME,
        outdir / hcv_cluster_audit.METADATA_SURVIVING_FILENAME,
    )
    hcv_cluster_audit.write_surviving_fasta(
        audits, plan.records, outdir / hcv_cluster_audit.SURVIVING_FASTA_FILENAME
    )
    if plan.has_duplicates:
        hcv_cluster_variants.write_variants_csv(plan, outdir / hcv_cluster_variants.VARIANTS_FILENAME)
    print(f"Sample audit: {hcv_cluster_audit.summarize(audits)} ({hcv_cluster_audit.AUDIT_FILENAME})")
    absent = sorted(
        {str(row["sample_id"]) for row in metadata_rows}.difference(a.base_sample_id for a in audits)
    )
    if absent:
        shown = ", ".join(absent[:5]) + (f", ... (+{len(absent) - 5} more)" if len(absent) > 5 else "")
        print(
            f"NOTE: metadata has {len(absent)} sample(s) not present in the input FASTA: {shown}",
            file=sys.stderr,
        )
    return audits


# --- command line -----------------------------------------------------------------------


def _load_metadata(path: Path | None) -> list[dict[str, Any]]:
    return hcv_cluster_metadata.load_metadata_csv(path) if path is not None else []


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hcv-stage",
        description="Pipeline stages used by the Nextflow workflow (and by hcv-cluster run).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("resolve-ids", help="Drop exact repeats; rename distinct same-ID sequences")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--out-fasta", type=Path, required=True)
    p.add_argument("--out-variants", type=Path, required=True)
    p.set_defaults(func=_cmd_resolve_ids)

    p = sub.add_parser("plan-genotypes", help="Set aside per-subtype FASTAs that cannot be analysed")
    p.add_argument("--genotypes", type=Path, required=True)
    p.add_argument("--fasta-dir", type=Path, required=True)
    p.add_argument("--skipped-dir", type=Path, default=Path("genotype_fastas_skipped"))
    p.add_argument("--genotype", action="append", default=None)
    p.add_argument("--reference-map", type=Path, default=None)
    p.set_defaults(func=_cmd_plan_genotypes)

    p = sub.add_parser("select-variants", help="Keep one usable same-genotype variant per sample")
    p.add_argument("--variants", type=Path, required=True)
    p.add_argument("--genotypes", type=Path, required=True)
    p.add_argument("--genotype", required=True)
    p.add_argument("--qc", type=Path, required=True)
    p.add_argument("--clustering-fasta", type=Path, required=True)
    p.add_argument("--mode", choices=hcv_cluster_variants.SELECTION_MODES, default="region-coverage")
    p.add_argument("--out-fasta", type=Path, required=True)
    p.add_argument("--out-decisions", type=Path, required=True)
    p.set_defaults(func=_cmd_select_variants)

    p = sub.add_parser("cluster", help="Distances, links and clusters for one subtype and method")
    p.add_argument("--method", choices=("tn93", "snp"), required=True)
    p.add_argument("--input", type=Path, required=True, help="Clustering FASTA")
    p.add_argument("--out-distances", type=Path, required=True)
    p.add_argument("--out-links", type=Path, required=True)
    p.add_argument("--out-clusters", type=Path, required=True)
    p.add_argument("--region", default="e1-e2")
    p.add_argument("--threshold", type=float, default=None)
    p.add_argument("--snp-count-threshold", type=int, default=None)
    p.add_argument("--tn93", default="tn93")
    p.add_argument("--ambiguities", default="average")
    p.add_argument("--ambiguity-fraction", type=float, default=None)
    p.add_argument("--min-overlap", type=int, default=None)
    p.add_argument("--extra-tn93-args", default="")
    p.add_argument("--variants", type=Path, default=None,
                   help="sequence_variants.csv, to add base_sample_id / mixed_infection columns")
    p.add_argument("--verbose", action="store_true")
    p.set_defaults(func=_cmd_cluster)

    p = sub.add_parser("merge", help="Merge per-subtype cluster and link tables")
    p.add_argument("--clusters", nargs="+", default=[], metavar="GENOTYPE=PATH")
    p.add_argument("--links", nargs="+", default=[], metavar="GENOTYPE=PATH")
    p.add_argument("--tables-dir", type=Path, default=None,
                   help="Instead of --clusters/--links: a directory of <genotype>.<method>.clusters.csv "
                        "and <genotype>.<method>.links.csv files")
    p.add_argument("--method", default="tn93", help="Method name used in --tables-dir file names")
    p.add_argument("--variants", type=Path, default=None)
    p.add_argument("--out-clusters", type=Path, required=True)
    p.add_argument("--out-links", type=Path, required=True)
    p.set_defaults(func=_cmd_merge)

    p = sub.add_parser("label", help="Assign persistent cluster labels from a registry")
    p.add_argument("--clusters", type=Path, required=True)
    p.add_argument("--registry", type=Path, required=True)
    p.add_argument("--update-registry", action="store_true")
    p.add_argument("--registry-out", type=Path, default=None,
                   help="Write the updated registry here instead of in place (what a workflow "
                        "process needs, since its input directory is read-only)")
    p.add_argument("--metric", choices=("tn93", "snp"), default="tn93")
    p.add_argument("--region", default="e1-e2")
    p.add_argument("--threshold", type=float, default=None)
    p.add_argument("--region-strategy", default="fixed")
    p.add_argument("--ambiguities", default="average")
    p.add_argument("--ambiguity-fraction", type=float, default=None)
    p.add_argument("--min-overlap", type=int, default=None)
    p.add_argument("--extra-tn93-args", default="")
    p.add_argument("--snp-count-threshold", type=int, default=None)
    p.add_argument("--label-prefix", default=hcv_cluster_labels.DEFAULT_PREFIX)
    p.add_argument("--label-date", default=None)
    p.add_argument("--out-clusters", type=Path, required=True)
    p.add_argument("--out-labels", type=Path, default=Path("cluster_labels.csv"))
    p.add_argument("--out-changes", type=Path, default=Path("label_changes.csv"))
    p.set_defaults(func=_cmd_label)

    p = sub.add_parser("audit", help="Sample audit, metadata tables and labelled survivor FASTA")
    p.add_argument("--input", type=Path, required=True, help="Original input FASTA")
    p.add_argument("--resolved-fasta", type=Path, default=None)
    p.add_argument("--variants", type=Path, default=None)
    p.add_argument("--genotypes", type=Path, required=True)
    p.add_argument("--qc", nargs="*", default=[], metavar="GENOTYPE=PATH")
    p.add_argument("--validation", nargs="*", default=[], metavar="GENOTYPE=PATH")
    p.add_argument("--decisions", nargs="*", default=[], metavar="GENOTYPE=PATH")
    p.add_argument("--clusters", type=Path, default=None)
    p.add_argument("--metadata", type=Path, default=None)
    p.add_argument("--min-query-coverage", type=float, default=0.5)
    p.add_argument("--min-identity", type=float, default=0.75)
    p.add_argument("--min-coverage", type=float, default=0.7)
    p.add_argument("--genotype", action="append", default=None)
    p.add_argument("--reference-map", type=Path, default=None)
    p.add_argument("--outdir", type=Path, default=Path("."))
    p.set_defaults(func=_cmd_audit)
    return parser


def _cmd_resolve_ids(args: argparse.Namespace) -> int:
    _records, plan = resolve_ids(args.input, args.out_fasta, args.out_variants)
    if plan.has_duplicates:
        repeats = sum(1 for v in plan.variants if v.status == "dropped_identical_repeat")
        print(
            f"NOTE: repeated FASTA IDs found. {len(plan.mixed_sample_ids())} sample(s) have distinct "
            f"sequences; {repeats} identical repeat(s) dropped."
        )
    return 0


def _cmd_plan_genotypes(args: argparse.Namespace) -> int:
    runnable, unsupported = plan_genotypes(
        hcv_workflow.read_genotype_assignments(args.genotypes),
        args.genotype,
        hcv_cluster_prep.load_reference_catalog(args.reference_map),
    )
    warn_unsupported(unsupported)
    filter_genotype_fastas(runnable, args.fasta_dir, args.skipped_dir)
    print(f"{len(runnable)} subtype(s) to analyse: {', '.join(runnable) or 'none'}")
    return 0


def _cmd_select_variants(args: argparse.Namespace) -> int:
    plan = hcv_cluster_variants.read_variants_csv(args.variants)
    dropped = select_variants(
        plan,
        hcv_workflow.read_genotype_assignments(args.genotypes),
        args.genotype,
        args.qc,
        args.mode,
        args.clustering_fasta,
        args.out_fasta,
        args.out_decisions,
    )
    if dropped:
        print(f"{len(dropped)} same-sample variant(s) dropped in favour of a better sequence ({args.mode})")
    return 0


def _cmd_cluster(args: argparse.Namespace) -> int:
    threshold = hcv_cluster_methods.resolve_threshold(args.region, args.threshold, distance=args.method)
    cluster_sequences(
        args.method,
        args.input,
        out_distances=args.out_distances,
        out_links=args.out_links,
        out_clusters=args.out_clusters,
        threshold=threshold,
        snp_count_threshold=args.snp_count_threshold,
        tn93=args.tn93,
        ambiguities=args.ambiguities,
        ambiguity_fraction=args.ambiguity_fraction,
        min_overlap=args.min_overlap,
        extra_tn93_args=args.extra_tn93_args,
        quiet=not args.verbose,
    )
    annotate_clusters(args.out_clusters, args.variants)
    print(f"{args.method} threshold {threshold:.4g}")
    return 0


def _tables_in_dir(directory: Path, method: str, kind: str) -> dict[str, Path]:
    suffix = f".{method}.{kind}.csv"
    return {path.name.removesuffix(suffix): path for path in sorted(directory.glob(f"*{suffix}"))}


def _cmd_merge(args: argparse.Namespace) -> int:
    if args.tables_dir is not None:
        clusters = _tables_in_dir(args.tables_dir, args.method, "clusters")
        links = _tables_in_dir(args.tables_dir, args.method, "links")
    else:
        clusters, links = parse_labelled_paths(args.clusters), parse_labelled_paths(args.links)
    merge_tables(
        clusters,
        links,
        args.out_clusters,
        args.out_links,
        args.variants,
    )
    return 0


def _cmd_label(args: argparse.Namespace) -> int:
    label_clusters(
        args.clusters,
        registry=args.registry,
        update_registry=args.update_registry,
        track=args.metric,
        definition=label_definition(
            args.metric,
            args.region,
            args.threshold,
            args.region_strategy,
            args.ambiguities,
            args.ambiguity_fraction,
            args.min_overlap,
            args.extra_tn93_args,
            args.snp_count_threshold,
        ),
        prefix=args.label_prefix,
        label_date=args.label_date,
        out_clusters=args.out_clusters,
        out_labels=args.out_labels,
        out_changes=args.out_changes,
        registry_out=args.registry_out,
    )
    return 0


def _cmd_audit(args: argparse.Namespace) -> int:
    write_sample_outputs(
        args.outdir,
        input_fasta=args.input,
        genotypes_csv=args.genotypes,
        variants_csv=args.variants,
        resolved_fasta=args.resolved_fasta,
        qc=parse_labelled_paths(args.qc),
        validation=parse_labelled_paths(args.validation),
        decisions=parse_labelled_paths(args.decisions),
        clusters_csv=args.clusters,
        metadata_rows=_load_metadata(args.metadata),
        min_query_coverage=args.min_query_coverage,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        genotype_filter=args.genotype,
        reference_map=args.reference_map,
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (RuntimeError, ValueError, hcv_cluster_prep.HcvPrepError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
