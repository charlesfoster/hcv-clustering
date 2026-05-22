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
import hcv_workflow


GENOTYPE_THRESHOLDS: dict[str, float] = {
    "1a": 0.012,
    "1b": 0.012,
    "2a": 0.015,
    "2b": 0.015,
    "3a": 0.015,
}
FALLBACK_THRESHOLD: float = 0.015


def resolve_threshold(genotype: str, user_threshold: float | None) -> float:
    if user_threshold is not None:
        return user_threshold
    return GENOTYPE_THRESHOLDS.get(genotype.lower(), FALLBACK_THRESHOLD)


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
            threshold = resolve_threshold(genotype, args.threshold)
            prep_cmd = [
                "hcv_cluster_prep", "prep-align",
                "--input", str(genotype_fasta),
                "--genotype", genotype,
                "--out-prefix", str(prefix),
            ]
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
    print()

    cluster_tables: list[tuple[str, Path]] = []
    snp_cluster_tables: list[tuple[str, Path]] = []

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
        with suppress():
            rc = hcv_cluster_prep.main(prep_argv)
        if rc:
            return rc

        qc_path = Path(f"{prefix}.qc.csv")
        if qc_path.exists():
            qc_st = _read_qc_stats(qc_path)
            pct = f"{qc_st['passed'] / qc_st['total'] * 100:.0f}%" if qc_st["total"] else "n/a"
            print(f"--> {qc_st['total']} sequences, {qc_st['passed']} passed QC ({pct})")

        threshold = resolve_threshold(genotype, args.threshold)
        clustering_fasta = Path(f"{prefix}.clustering.fasta")

        if args.distance in ("tn93", "both"):
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
            cl_st = _read_cluster_stats(clusters_csv)
            metric = "TN93" if args.distance == "both" else "TN93"
            print(
                f"--> {metric} (threshold {threshold:.4g}): "
                f"{cl_st['n_multi']} cluster(s), {cl_st['singletons']} singleton(s), "
                f"largest: {cl_st['largest']}"
            )

        if args.distance in ("snp", "both"):
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
                threshold=threshold,
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
            else:
                snp_cluster_tables.append((genotype, snp_clusters_csv))
            snp_st = _read_cluster_stats(snp_clusters_csv)
            print(
                f"--> SNP   (threshold {threshold:.4g}): "
                f"{snp_st['n_multi']} cluster(s), {snp_st['singletons']} singleton(s), "
                f"largest: {snp_st['largest']}"
            )

        print()

    hcv_workflow.merge_cluster_tables(cluster_tables, outdir / "clusters.csv")
    if args.distance == "both":
        hcv_workflow.merge_cluster_tables(snp_cluster_tables, outdir / "clusters.snp.csv")

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
        help="Maximum distance for cluster linking. Default: 0.012 for 1a/1b, 0.015 for all others.",
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
        default="resolve",
        help="TN93 ambiguity handling: resolve (default), average, skip, gapmm.",
    )
    run_parser.add_argument(
        "-r", "--region",
        metavar="EXPR",
        default="core-e2",
        help=(
            "Reference-anchored region expression (default: core-e2). "
            "Individual regions: core e1 e2 p7 ns2 ns3 ns4a ns4b ns5a ns5b. "
            "Named presets: structural (=core-e2), envelope (=e1-e2), "
            "nonstructural (=ns2-ns5b), cds (whole coding sequence). "
            "Range syntax: first-last selects all regions from first through last "
            "(e.g. e1-e2, ns3-ns5b). "
            "Union syntax: a+b includes both (e.g. core-e2+ns3, e1-e2+ns5a)."
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
