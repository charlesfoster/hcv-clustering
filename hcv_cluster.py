#!/usr/bin/env python3
"""User-facing CLI for HCV genomic clustering with genotype-aware thresholds."""

from __future__ import annotations

import argparse
import contextlib
import csv
import logging
import os
import sys
from pathlib import Path

import assign_hcv_genotypes_from_fasta
import hcv_cluster_prep
import hcv_cluster_viz
import hcv_workflow


# Evidence-based per-region defaults; see docs/threshold_rationale.md for citations
# and the reasoning behind these specific values. Deliberately NOT genotype-split:
# the HCV Core-E2 literature applies a single threshold across genotypes, and the
# only Australian genotype-specific number available (different region and metric)
# does not support a 1a/1b-vs-rest split, so a genotype split isn't defensible as
# evidence-based.
#
# TN93: "core-e2-nohvr1" (the default region) matches Bartlett et al. 2017's
# Core-early-E2 TN93 pairwise/connected-components network directly (same metric,
# same clustering algorithm as this pipeline) — the strongest available precedent.
# "core-e2" (HVR1 included) and "ns5b" have no HVR1-inclusive/NS5B TN93 precedent,
# so they reuse Lamoury et al. 2015's uncorrected-p-distance values as a (slightly
# conservative, since TN93 >= p-distance for the same pair) approximation.
REGION_THRESHOLDS_TN93: dict[str, float] = {
    "core-e2-nohvr1": 0.03,
    "core-e2": 0.045,
    "ns5b": 0.015,
}

# SNP: this pipeline's SNP distance is uncorrected p-distance (fraction of
# ACGT-comparable sites that differ) -- exactly the metric Lamoury et al. 2015 used
# (MEGA v6 p-distance, partial deletion). Their region thresholds are therefore a
# direct match, not an approximation, for all three regions below.
REGION_THRESHOLDS_SNP: dict[str, float] = {
    "core-e2-nohvr1": 0.03,
    "core-e2": 0.045,
    "ns5b": 0.015,
}

# Backward-compat aliases: TN93 was the only metric these names covered originally.
REGION_THRESHOLDS: dict[str, float] = REGION_THRESHOLDS_TN93
FALLBACK_THRESHOLD: float = REGION_THRESHOLDS_TN93["core-e2-nohvr1"]
FALLBACK_THRESHOLD_SNP: float = REGION_THRESHOLDS_SNP["core-e2-nohvr1"]


def resolve_threshold(region: str, user_threshold: float | None, distance: str = "tn93") -> float:
    if user_threshold is not None:
        return user_threshold
    canonical_region = hcv_cluster_prep.expand_region_expression(region)
    if distance == "snp":
        return REGION_THRESHOLDS_SNP.get(canonical_region, FALLBACK_THRESHOLD_SNP)
    return REGION_THRESHOLDS_TN93.get(canonical_region, FALLBACK_THRESHOLD)


def region_threshold_is_evidence_based(region: str, distance: str = "tn93") -> bool:
    canonical_region = hcv_cluster_prep.expand_region_expression(region)
    table = REGION_THRESHOLDS_SNP if distance == "snp" else REGION_THRESHOLDS_TN93
    return canonical_region in table


def compute_snp_distances_detailed(
    fasta_path: Path,
) -> list[tuple[str, str, int, int, float]]:
    records = hcv_cluster_prep.read_fasta_records(fasta_path)
    results: list[tuple[str, str, int, int, float]] = []
    acgt = frozenset("ACGT")
    for i in range(len(records)):
        id_i = records[i].header.split()[0]
        seq_i = records[i].sequence.upper()
        for j in range(i + 1, len(records)):
            id_j = records[j].header.split()[0]
            seq_j = records[j].sequence.upper()
            snp_count = 0
            comparable_sites = 0
            for char_i, char_j in zip(seq_i, seq_j):
                if char_i in acgt and char_j in acgt:
                    comparable_sites += 1
                    if char_i != char_j:
                        snp_count += 1
            normalized = snp_count / comparable_sites if comparable_sites > 0 else 1.0
            results.append((id_i, id_j, snp_count, comparable_sites, normalized))
    return results


def compute_snp_distances(fasta_path: Path) -> list[hcv_workflow.DistanceRow]:
    detailed = compute_snp_distances_detailed(fasta_path)
    return [
        hcv_workflow.DistanceRow(id_i, id_j, normalized)
        for id_i, id_j, _snp_count, _comparable_sites, normalized in detailed
    ]


def write_snp_csv(path: Path, fasta_path: Path) -> list[hcv_workflow.DistanceRow]:
    detailed = compute_snp_distances_detailed(fasta_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("source", "target", "distance", "snp_count", "comparable_sites"),
        )
        writer.writeheader()
        for id_i, id_j, snp_count, comparable_sites, normalized in detailed:
            writer.writerow(
                {
                    "source": id_i,
                    "target": id_j,
                    "distance": f"{normalized:.10g}",
                    "snp_count": snp_count,
                    "comparable_sites": comparable_sites,
                }
            )
    return [
        hcv_workflow.DistanceRow(id_i, id_j, normalized)
        for id_i, id_j, _snp_count, _comparable_sites, normalized in detailed
    ]


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


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _render_cluster_plot(
    plot_network: str,
    hide_singletons: bool,
    genotype_dir: Path,
    clusters_csv: Path,
    links_csv: Path,
    file_prefix: str,
) -> None:
    if plot_network == "none":
        return
    node_rows = _read_csv_rows(clusters_csv)
    edge_rows = _read_csv_rows(links_csv)
    if hide_singletons:
        node_rows, edge_rows = hcv_cluster_viz.filter_singletons(node_rows, edge_rows)
    if not node_rows:
        return
    try:
        fig = hcv_cluster_viz.build_network_figure(node_rows, edge_rows)
        if plot_network in ("png", "both"):
            fig.write_image(genotype_dir / f"{file_prefix}network.png", scale=2)
        if plot_network in ("html", "both"):
            fig.write_html(genotype_dir / f"{file_prefix}network.html", include_plotlyjs=True)
    except Exception as exc:
        print(f"WARNING: could not render network plot for {genotype_dir.name}: {exc}")


def command_run(args: argparse.Namespace) -> int:
    outdir: Path = args.outdir
    input_path: Path = Path(args.input)
    genotype_fasta_dir = outdir / "genotype_fastas"
    genotype_csv = outdir / "genotypes.csv"

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
                prep_cmd.extend(
                    ["--cached-alignment", str(Path(f"{cached_prefix}.aligned.fasta"))]
                )
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
    suppress = _quiet_logging if not args.verbose else contextlib.nullcontext

    # --- Genotyping ---
    print(f"Checking subtypes for: {input_path}")

    genotype_argv = [
        "--input", str(input_path),
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
        raise RuntimeError("No passing genotype assignments were available to process")

    gt_counts = {
        gt: sum(1 for r in genotype_rows if r.get("assigned_genotype") == gt)
        for gt in genotypes
    }
    gt_summary = ", ".join(f"{gt} ({gt_counts[gt]})" for gt in genotypes)
    print(f"--> {len(genotypes)} subtype(s) found: {gt_summary}")
    if args.threshold is None and not region_threshold_is_evidence_based(args.region):
        print(
            f"WARNING: no HCV-specific clustering threshold evidence for region '{args.region}'; "
            f"using the core-E2 default ({FALLBACK_THRESHOLD:.3g}) as a starting point. "
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
                print(f"--> Reusing unchanged sequences from: {cached_alignment}")
            else:
                print(f"--> No cached alignment for subtype {genotype}; aligning all sequences")
        with suppress():
            rc = hcv_cluster_prep.main(prep_argv)
        if rc:
            return rc

        qc_path = Path(f"{prefix}.qc.csv")
        if qc_path.exists():
            qc_st = _read_qc_stats(qc_path)
            pct = f"{qc_st['passed'] / qc_st['total'] * 100:.0f}%" if qc_st["total"] else "n/a"
            print(f"--> {qc_st['total']} sequences, {qc_st['passed']} passed QC ({pct})")

        clustering_fasta = Path(f"{prefix}.clustering.fasta")

        if args.distance in ("tn93", "both"):
            threshold = resolve_threshold(args.region, args.threshold, distance="tn93")
            tn93_csv = genotype_dir / "tn93.csv"
            links_csv = genotype_dir / "links.csv"
            clusters_csv = genotype_dir / "clusters.csv"
            tn93_args = argparse.Namespace(
                input=clustering_fasta,
                output=tn93_csv,
                tn93=args.tn93,
                tn93_threshold=threshold,
                ambiguities=args.ambiguities,
                ambiguity_fraction=args.ambiguity_fraction,
                min_overlap=args.min_overlap,
                quiet=not args.verbose,
                extra_tn93_args=args.extra_tn93_args,
                dry_run=False,
            )
            rc = hcv_workflow.run_tn93(tn93_args)
            if rc:
                return rc
            link_args = argparse.Namespace(
                distances=tn93_csv,
                output=links_csv,
                threshold=threshold,
                include_self=False,
            )
            hcv_workflow.command_link(link_args)
            cluster_args = argparse.Namespace(
                links=links_csv,
                output=clusters_csv,
                nodes_fasta=clustering_fasta,
                nodes_file=None,
            )
            hcv_workflow.command_cluster(cluster_args)
            cluster_tables.append((genotype, clusters_csv))
            link_tables.append((genotype, links_csv))
            cl_st = _read_cluster_stats(clusters_csv)
            metric = "TN93" if args.distance == "both" else "TN93"
            print(
                f"--> {metric} (threshold {threshold:.4g}): "
                f"{cl_st['n_multi']} cluster(s), {cl_st['singletons']} singleton(s), "
                f"largest: {cl_st['largest']}"
            )
            _render_cluster_plot(
                args.plot_network, args.plot_hide_singletons, genotype_dir, clusters_csv, links_csv, ""
            )

        if args.distance in ("snp", "both"):
            snp_threshold = resolve_threshold(args.region, args.threshold, distance="snp")
            snp_csv = genotype_dir / "snp.csv"
            write_snp_csv(snp_csv, clustering_fasta)
            if args.distance == "snp":
                snp_links_csv = genotype_dir / "links.csv"
                snp_clusters_csv = genotype_dir / "clusters.csv"
            else:
                snp_links_csv = genotype_dir / "snp_links.csv"
                snp_clusters_csv = genotype_dir / "snp_clusters.csv"
            link_args = argparse.Namespace(
                distances=snp_csv,
                output=snp_links_csv,
                threshold=snp_threshold,
                include_self=False,
            )
            hcv_workflow.command_link(link_args)
            cluster_args = argparse.Namespace(
                links=snp_links_csv,
                output=snp_clusters_csv,
                nodes_fasta=clustering_fasta,
                nodes_file=None,
            )
            hcv_workflow.command_cluster(cluster_args)
            if args.distance == "snp":
                cluster_tables.append((genotype, snp_clusters_csv))
                link_tables.append((genotype, snp_links_csv))
            else:
                snp_cluster_tables.append((genotype, snp_clusters_csv))
                snp_link_tables.append((genotype, snp_links_csv))
            snp_st = _read_cluster_stats(snp_clusters_csv)
            print(
                f"--> SNP   (threshold {snp_threshold:.4g}): "
                f"{snp_st['n_multi']} cluster(s), {snp_st['singletons']} singleton(s), "
                f"largest: {snp_st['largest']}"
            )
            snp_file_prefix = "" if args.distance == "snp" else "snp_"
            _render_cluster_plot(
                args.plot_network,
                args.plot_hide_singletons,
                genotype_dir,
                snp_clusters_csv,
                snp_links_csv,
                snp_file_prefix,
            )

        print()

    hcv_workflow.merge_cluster_tables(cluster_tables, outdir / "clusters.csv")
    hcv_workflow.merge_link_tables(link_tables, outdir / "links.csv")
    if args.distance == "both":
        hcv_workflow.merge_cluster_tables(snp_cluster_tables, outdir / "clusters.snp.csv")
        hcv_workflow.merge_link_tables(snp_link_tables, outdir / "links.snp.csv")

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
        "-t", "--threshold",
        metavar="FLOAT",
        type=float,
        default=None,
        help=(
            "Maximum distance for cluster linking. Default: region-dependent "
            "(0.03 for core-e2-nohvr1, the default region; 0.045 for core-e2 "
            "with HVR1 included; 0.015 for ns5b; 0.03 for any other region — "
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
        default="core-e2-nohvr1",
        help=(
            "Reference-anchored region expression (default: core-e2-nohvr1, i.e. "
            "Core-through-E2 with HVR1 masked — see docs/threshold_rationale.md). "
            "Individual regions: core e1 e2 p7 ns2 ns3 ns4a ns4b ns5a ns5b. "
            "core-e2 is the same span with HVR1 included. "
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
        "--min-coverage",
        metavar="FLOAT",
        type=float,
        default=0.8,
        help=adv("Minimum non-gap non-N selected-region coverage fraction (default: 0.8)"),
    )
    run_parser.add_argument(
        "--min-identity",
        metavar="FLOAT",
        type=float,
        default=0.75,
        help=adv("Minimum minimap2 alignment identity for genotype pass (default: 0.75)"),
    )
    run_parser.add_argument(
        "--min-query-coverage",
        metavar="FLOAT",
        type=float,
        default=0.50,
        help=adv("Minimum minimap2 query coverage fraction for genotype pass (default: 0.50)"),
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
        default="asm10",
        help=adv("minimap2 preset (default: asm10)"),
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
        default=Path("refs"),
        help=adv("Directory for cached reference GenBank/FASTA files (default: refs)"),
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
        choices=("none", "png", "html", "both"),
        default="none",
        help=adv(
            "Render a cluster network plot per genotype (and metric, if --distance both) "
            "into by_genotype/<genotype>/. png: static image. html: interactive, opens in a "
            "browser. both: write both. Default: none."
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
        default=Path("refs"),
        help="Reference cache directory (default: refs)",
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
