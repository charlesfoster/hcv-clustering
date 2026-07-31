#!/usr/bin/env python3
"""Unified HCV preprocessing, genotype, TN93, linking, and clustering CLI."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import assign_hcv_genotypes_from_fasta
import hcv_cluster_prep


@dataclass(frozen=True)
class DistanceRow:
    source: str
    target: str
    distance: float


@dataclass(frozen=True)
class LinkRow:
    source: str
    target: str
    distance: float | None = None


def read_fasta_ids(path: Path) -> list[str]:
    records = hcv_cluster_prep.read_fasta_records(path)
    hcv_cluster_prep.reject_duplicate_ids(records)
    return [record.header.split()[0] for record in records]


def read_nodes_file(path: Path) -> list[str]:
    nodes: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if not sample.strip():
            return nodes
        dialect = sniff_dialect(sample)
        reader = csv.DictReader(handle, dialect=dialect)
        if reader.fieldnames:
            field = first_present(reader.fieldnames, ("node", "sample_id", "sequence_id", "id"))
            if field is not None:
                for row in reader:
                    value = (row.get(field) or "").strip()
                    if value:
                        nodes.append(value)
                return unique_in_order(nodes)

    with path.open(encoding="utf-8") as handle:
        return unique_in_order(line.strip().split()[0] for line in handle if line.strip())


def unique_in_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def first_present(fieldnames: Iterable[str], candidates: Iterable[str]) -> str | None:
    by_lower = {field.lower(): field for field in fieldnames}
    for candidate in candidates:
        if candidate.lower() in by_lower:
            return by_lower[candidate.lower()]
    return None


def sniff_dialect(sample: str) -> csv.Dialect:
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        dialect = csv.excel_tab if "\t" in sample else csv.excel
        return dialect


def parse_distance_table(path: Path) -> list[DistanceRow]:
    with path.open(newline="", encoding="utf-8") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if not sample.strip():
            return []
        dialect = sniff_dialect(sample)
        reader = csv.reader(handle, dialect=dialect)
        rows = [row for row in reader if row and any(cell.strip() for cell in row)]

    if not rows:
        return []

    header = [cell.strip() for cell in rows[0]]
    source_i, target_i, distance_i, data_start = distance_column_indices(header)
    parsed: list[DistanceRow] = []
    for line_number, row in enumerate(rows[data_start:], start=data_start + 1):
        if max(source_i, target_i, distance_i) >= len(row):
            raise ValueError(f"Distance row {line_number} has too few columns in {path}")
        source = row[source_i].strip()
        target = row[target_i].strip()
        raw_distance = row[distance_i].strip()
        if not source or not target:
            continue
        try:
            distance = float(raw_distance)
        except ValueError as exc:
            raise ValueError(f"Invalid distance '{raw_distance}' on row {line_number} in {path}") from exc
        parsed.append(DistanceRow(source, target, distance))
    return parsed


def distance_column_indices(header: list[str]) -> tuple[int, int, int, int]:
    lowered = [field.lower() for field in header]
    source_names = ("source", "query", "id1", "seq1", "sequence_1", "name1", "sample_id_1")
    target_names = ("target", "reference", "id2", "seq2", "sequence_2", "name2", "sample_id_2")
    distance_names = ("distance", "tn93", "tn93_distance", "dist")

    source_i = next((lowered.index(name) for name in source_names if name in lowered), None)
    target_i = next((lowered.index(name) for name in target_names if name in lowered), None)
    distance_i = next((lowered.index(name) for name in distance_names if name in lowered), None)
    if source_i is not None and target_i is not None and distance_i is not None:
        return source_i, target_i, distance_i, 1

    try:
        float(header[2])
    except (IndexError, ValueError):
        pass
    else:
        return 0, 1, 2, 0

    raise ValueError(
        "Could not identify distance table columns. Expected a header containing "
        "source/id1/seq1, target/id2/seq2, and distance/tn93 columns, or a three-column table."
    )


def threshold_distances(rows: Iterable[DistanceRow], threshold: float, include_self: bool = False) -> list[LinkRow]:
    if threshold < 0:
        raise ValueError("--threshold must be non-negative")
    links: list[LinkRow] = []
    for row in rows:
        if not include_self and row.source == row.target:
            continue
        if row.distance <= threshold:
            links.append(LinkRow(row.source, row.target, row.distance))
    return links


def write_links(path: Path, links: Iterable[LinkRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("source", "target", "distance"))
        writer.writeheader()
        for link in links:
            writer.writerow(
                {
                    "source": link.source,
                    "target": link.target,
                    "distance": "" if link.distance is None else f"{link.distance:.10g}",
                }
            )


def parse_links(path: Path) -> list[LinkRow]:
    with path.open(newline="", encoding="utf-8") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if not sample.strip():
            return []
        dialect = sniff_dialect(sample)
        reader = csv.DictReader(handle, dialect=dialect)
        if reader.fieldnames:
            source_field = first_present(reader.fieldnames, ("source", "id1", "seq1", "sample_id_1"))
            target_field = first_present(reader.fieldnames, ("target", "id2", "seq2", "sample_id_2"))
            distance_field = first_present(reader.fieldnames, ("distance", "tn93", "tn93_distance", "dist"))
            if source_field and target_field:
                links: list[LinkRow] = []
                for row in reader:
                    source = (row.get(source_field) or "").strip()
                    target = (row.get(target_field) or "").strip()
                    distance = None
                    if distance_field and (row.get(distance_field) or "").strip():
                        distance = float(row[distance_field])
                    if source and target:
                        links.append(LinkRow(source, target, distance))
                return links

    with path.open(encoding="utf-8") as handle:
        links = []
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").replace(",", "\t").split("\t")
            if len(parts) < 2:
                raise ValueError(f"Link row {line_number} has fewer than two columns in {path}")
            links.append(LinkRow(parts[0].strip(), parts[1].strip(), None))
        return links


def connected_components(links: Iterable[LinkRow], nodes: Iterable[str] = ()) -> dict[str, str]:
    graph: dict[str, set[str]] = defaultdict(set)
    for node in nodes:
        graph[node]
    for link in links:
        graph[link.source].add(link.target)
        graph[link.target].add(link.source)

    assignments: dict[str, str] = {}
    cluster_index = 1
    for start in sorted(graph):
        if start in assignments:
            continue
        cluster_id = f"C{cluster_index:04d}"
        cluster_index += 1
        queue: deque[str] = deque([start])
        assignments[start] = cluster_id
        while queue:
            node = queue.popleft()
            for neighbor in sorted(graph[node]):
                if neighbor not in assignments:
                    assignments[neighbor] = cluster_id
                    queue.append(neighbor)
    return assignments


def write_clusters(path: Path, assignments: dict[str, str]) -> None:
    sizes: dict[str, int] = defaultdict(int)
    for cluster_id in assignments.values():
        sizes[cluster_id] += 1

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("sample_id", "cluster_id", "cluster_size"))
        writer.writeheader()
        for sample_id in sorted(assignments):
            cluster_id = assignments[sample_id]
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "cluster_id": cluster_id,
                    "cluster_size": sizes[cluster_id],
                }
            )


def read_genotype_assignments(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def detected_pass_genotypes(rows: Iterable[dict[str, str]]) -> list[str]:
    return sorted(
        {
            (row.get("assigned_genotype") or "").strip()
            for row in rows
            if (row.get("assignment_status") or "").strip() == "pass"
            and (row.get("assigned_genotype") or "").strip()
        }
    )


def safe_genotype_name(genotype: str) -> str:
    return "".join(char if char.isalnum() or char in {"-", "_", "."} else "_" for char in genotype)


def merge_cluster_tables(inputs: Iterable[tuple[str, Path]], output: Path) -> None:
    rows: list[dict[str, str | int]] = []
    for genotype, path in inputs:
        safe_genotype = safe_genotype_name(genotype)
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                source_cluster_id = row["cluster_id"]
                rows.append(
                    {
                        "sample_id": row["sample_id"],
                        "genotype": genotype,
                        "cluster_id": f"{safe_genotype}_{source_cluster_id}",
                        "source_cluster_id": source_cluster_id,
                        "cluster_size": int(row["cluster_size"]),
                    }
                )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sample_id", "genotype", "cluster_id", "source_cluster_id", "cluster_size"),
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (str(row["genotype"]), str(row["sample_id"]))))


def merge_link_tables(inputs: Iterable[tuple[str, Path]], output: Path) -> None:
    rows: list[dict[str, str]] = []
    for _genotype, path in inputs:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rows.append({"source": row["source"], "target": row["target"], "distance": row["distance"]})

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("source", "target", "distance"))
        writer.writeheader()
        writer.writerows(rows)


def run_tn93(args: argparse.Namespace) -> int:
    executable = shutil.which(args.tn93) if not Path(args.tn93).exists() else args.tn93
    if executable is None:
        raise RuntimeError(f"Could not find TN93 executable '{args.tn93}'. Install it with `pixi install`.")

    command = [str(executable), *shell_split(args.extra_tn93_args)]
    if args.tn93_threshold is not None:
        command.extend(["-t", str(args.tn93_threshold)])
    if args.ambiguities:
        command.extend(["-a", args.ambiguities])
    if args.ambiguity_fraction is not None:
        command.extend(["-g", str(args.ambiguity_fraction)])
    if args.min_overlap is not None:
        command.extend(["-l", str(args.min_overlap)])
    if args.quiet:
        command.append("-q")
    command.extend(["-o", str(args.output), str(args.input)])

    if args.dry_run:
        print(" ".join(command))
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, stderr=subprocess.PIPE, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"TN93 failed with exit code {completed.returncode}: "
            f"{completed.stderr.strip() or 'no stderr captured'}"
        )
    if not args.output.exists():
        raise RuntimeError(f"TN93 completed but did not create expected output: {args.output}")
    return 0


def shell_split(value: str) -> list[str]:
    if not value:
        return []
    return shlex.split(value)


def shell_join(command: Iterable[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def command_link(args: argparse.Namespace) -> int:
    rows = parse_distance_table(args.distances)
    links = threshold_distances(rows, args.threshold, args.include_self)
    write_links(args.output, links)
    return 0


def command_cluster(args: argparse.Namespace) -> int:
    nodes = []
    if args.nodes_fasta:
        nodes.extend(read_fasta_ids(args.nodes_fasta))
    if args.nodes_file:
        nodes.extend(read_nodes_file(args.nodes_file))
    links = parse_links(args.links)
    assignments = connected_components(links, unique_in_order(nodes))
    write_clusters(args.output, assignments)
    return 0


def command_all(args: argparse.Namespace) -> int:
    args.outdir.mkdir(parents=True, exist_ok=True)
    genotype_fasta_dir = args.outdir / "genotype_fastas"
    genotype_csv = args.outdir / "genotypes.csv"
    clusters_csv = args.outdir / "clusters.csv"

    if args.dry_run:
        print(
            shell_join(
                [
                    "hcv-workflow",
                    "genotype",
                    "--input",
                    str(args.input),
                    "--output-csv",
                    str(genotype_csv),
                    "--outdir",
                    str(genotype_fasta_dir),
                    "--write-genotype-fastas",
                ]
            )
        )
        genotype_hint = args.genotype or ["<each detected pass genotype>"]
        for genotype in genotype_hint:
            safe_genotype = safe_genotype_name(genotype)
            genotype_dir = args.outdir / "by_genotype" / safe_genotype
            prefix = genotype_dir / args.prefix
            genotype_fasta = genotype_fasta_dir / f"hcv_{genotype}.fasta"
            clustering_fasta = Path(f"{prefix}.clustering.fasta")
            tn93_csv = genotype_dir / "tn93.csv"
            links_csv = genotype_dir / "links.csv"
            genotype_clusters_csv = genotype_dir / "clusters.csv"
            print(
                shell_join(
                    [
                        "hcv-workflow",
                        "prep",
                        "--input",
                        str(genotype_fasta),
                        "--genotype",
                        genotype,
                        "--out-prefix",
                        str(prefix),
                    ]
                )
            )
            tn93_cmd = [
                args.tn93,
                *shell_split(args.extra_tn93_args),
                "-t",
                str(args.tn93_threshold if args.tn93_threshold is not None else args.threshold),
                "-a",
                args.ambiguities,
            ]
            if args.ambiguity_fraction is not None:
                tn93_cmd.extend(["-g", str(args.ambiguity_fraction)])
            if args.min_overlap is not None:
                tn93_cmd.extend(["-l", str(args.min_overlap)])
            if args.quiet:
                tn93_cmd.append("-q")
            tn93_cmd.extend(["-o", str(tn93_csv), str(clustering_fasta)])
            print(shell_join(tn93_cmd))
            print(
                shell_join(
                    [
                        "hcv-workflow",
                        "link",
                        "--distances",
                        str(tn93_csv),
                        "--output",
                        str(links_csv),
                        "--threshold",
                        str(args.threshold),
                    ]
                )
            )
            print(
                shell_join(
                    [
                        "hcv-workflow",
                        "cluster",
                        "--links",
                        str(links_csv),
                        "--nodes-fasta",
                        str(clustering_fasta),
                        "--output",
                        str(genotype_clusters_csv),
                    ]
                )
            )
        print(f"# merge per-genotype cluster tables -> {clusters_csv}")
        return 0

    genotype_argv = [
        "--input",
        str(args.input),
        "--output-csv",
        str(genotype_csv),
        "--outdir",
        str(genotype_fasta_dir),
        "--panel-fasta",
        str(args.panel_fasta),
        "--minimap2",
        args.minimap2,
        "--preset",
        args.minimap2_preset,
        "--max-secondary",
        str(args.max_secondary),
        "--min-query-coverage",
        str(args.min_query_coverage),
        "--min-identity",
        str(args.min_identity),
        "--close-hit-fraction",
        str(args.close_hit_fraction),
        "--write-genotype-fastas",
    ]
    if args.extra_minimap2_args:
        genotype_argv.extend(["--extra-minimap2-args", args.extra_minimap2_args])
    if args.keep_paf:
        genotype_argv.append("--keep-paf")
    rc = assign_hcv_genotypes_from_fasta.main(genotype_argv)
    if rc:
        return rc

    genotype_rows = read_genotype_assignments(genotype_csv)
    genotypes = detected_pass_genotypes(genotype_rows)
    if args.genotype:
        requested = set(args.genotype)
        genotypes = [genotype for genotype in genotypes if genotype in requested]
    if not genotypes:
        raise RuntimeError("No passing genotype assignments were available to process")

    cluster_tables: list[tuple[str, Path]] = []
    for genotype in genotypes:
        safe_genotype = safe_genotype_name(genotype)
        genotype_dir = args.outdir / "by_genotype" / safe_genotype
        genotype_dir.mkdir(parents=True, exist_ok=True)
        prefix = genotype_dir / args.prefix
        genotype_fasta = genotype_fasta_dir / f"hcv_{genotype}.fasta"
        if not genotype_fasta.exists():
            raise RuntimeError(f"Expected per-genotype FASTA was not created: {genotype_fasta}")

        prep_argv = [
            "prep-align",
            "--input",
            str(genotype_fasta),
            "--genotype",
            genotype,
            "--out-prefix",
            str(prefix),
            "--refs-dir",
            str(args.refs_dir),
            "--min-coverage",
            str(args.min_coverage),
            "--region-strategy",
            args.region_strategy,
            "--mafft",
            args.mafft,
            "--threads",
            str(args.threads),
            "--genotype-validation",
            args.genotype_validation,
            "--genotype-validation-kmer-size",
            str(args.genotype_validation_kmer_size),
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
        rc = hcv_cluster_prep.main(prep_argv)
        if rc:
            return rc

        clustering_fasta = Path(f"{prefix}.clustering.fasta")
        tn93_csv = genotype_dir / "tn93.csv"
        links_csv = genotype_dir / "links.csv"
        genotype_clusters_csv = genotype_dir / "clusters.csv"
        tn93_args = argparse.Namespace(
            input=clustering_fasta,
            output=tn93_csv,
            tn93=args.tn93,
            tn93_threshold=args.tn93_threshold if args.tn93_threshold is not None else args.threshold,
            ambiguities=args.ambiguities,
            ambiguity_fraction=args.ambiguity_fraction,
            min_overlap=args.min_overlap,
            quiet=args.quiet,
            extra_tn93_args=args.extra_tn93_args,
            dry_run=False,
        )
        rc = run_tn93(tn93_args)
        if rc:
            return rc

        link_args = argparse.Namespace(
            distances=tn93_csv,
            output=links_csv,
            threshold=args.threshold,
            include_self=False,
        )
        command_link(link_args)
        cluster_args = argparse.Namespace(
            links=links_csv,
            output=genotype_clusters_csv,
            nodes_fasta=clustering_fasta,
            nodes_file=None,
        )
        command_cluster(cluster_args)
        cluster_tables.append((genotype, genotype_clusters_csv))

    merge_cluster_tables(cluster_tables, clusters_csv)
    return 0


def add_link_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("link", help="Threshold a TN93 distance table into linked pairs")
    parser.add_argument("--distances", required=True, type=Path, help="TN93 pairwise distance CSV/TSV")
    parser.add_argument("--output", required=True, type=Path, help="Output linked-pairs CSV")
    parser.add_argument("--threshold", required=True, type=float, help="Maximum TN93 distance to link")
    parser.add_argument("--include-self", action="store_true", help="Keep self links if present in the distance table")
    parser.set_defaults(func=command_link)


def add_cluster_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("cluster", help="Build connected-component clusters from linked pairs")
    parser.add_argument("--links", required=True, type=Path, help="Linked-pairs CSV/TSV edge list")
    parser.add_argument("--output", required=True, type=Path, help="Output cluster membership CSV")
    parser.add_argument("--nodes-fasta", type=Path, help="FASTA whose sequence IDs should be included as nodes")
    parser.add_argument("--nodes-file", type=Path, help="Optional node table/list for including isolates")
    parser.set_defaults(func=command_cluster)


def add_tn93_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("tn93", help="Run TN93 on an aligned FASTA")
    parser.add_argument("--input", required=True, type=Path, help="Aligned FASTA input")
    parser.add_argument("--output", required=True, type=Path, help="Output TN93 distance table")
    parser.add_argument("--tn93", default="tn93", help="TN93 executable path/name")
    parser.add_argument(
        "--tn93-threshold",
        type=float,
        default=1.0,
        help="TN93 reporting threshold passed with -t. Use a value at least as high as the link threshold.",
    )
    parser.add_argument(
        "--ambiguities",
        default="average",
        help="TN93 ambiguity handling passed with -a: average, resolve, skip, gapmm, or an ambiguity list",
    )
    parser.add_argument("--ambiguity-fraction", type=float, help="TN93 maximum resolvable ambiguity fraction passed with -g")
    parser.add_argument("--min-overlap", type=int, help="Minimum pairwise overlap passed to TN93 with -l")
    parser.add_argument("--quiet", action="store_true", help="Pass -q to TN93")
    parser.add_argument("--extra-tn93-args", default="", help="Extra TN93 arguments as one shell-like string")
    parser.add_argument("--dry-run", action="store_true", help="Print the TN93 command without running it")
    parser.set_defaults(func=run_tn93)


def add_all_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("all", help="Run genotype, prep, TN93, linking, and clustering")
    parser.add_argument("--input", required=True, type=Path, help="Raw input FASTA")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--prefix", default="prep", help="Prefix for prep outputs inside --outdir")
    parser.add_argument(
        "--genotype",
        action="append",
        help="Optional genotype/subtype to process after detection. May be repeated; defaults to all detected pass genotypes.",
    )
    parser.add_argument("--threshold", required=True, type=float, help="Maximum TN93 distance to link")
    parser.add_argument("--panel-fasta", type=Path, default=assign_hcv_genotypes_from_fasta.DEFAULT_PANEL)
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--minimap2-preset", default="asm10")
    parser.add_argument("--extra-minimap2-args", default="")
    parser.add_argument("--max-secondary", type=int, default=50)
    parser.add_argument("--min-query-coverage", type=float, default=0.50)
    parser.add_argument("--min-identity", type=float, default=0.75)
    parser.add_argument("--close-hit-fraction", type=float, default=0.98)
    parser.add_argument("--keep-paf", action="store_true")
    parser.add_argument("--refs-dir", type=Path, default=Path("refs"))
    parser.add_argument("--reference-map", type=Path)
    parser.add_argument("--ncbi-email")
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--region")
    parser.add_argument("--region-strategy", choices=("fixed", "max-usable"), default="fixed")
    parser.add_argument("--min-coverage", type=float, default=0.8)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--genotype-validation", choices=("auto", "headers", "kmer", "none"), default="auto")
    parser.add_argument("--genotype-validation-kmer-size", type=int, default=15)
    parser.add_argument("--keep-temp", action="store_true")
    parser.add_argument("--tn93", default="tn93")
    parser.add_argument("--tn93-threshold", type=float, help="TN93 reporting threshold; defaults to --threshold in all mode")
    parser.add_argument("--ambiguities", default="average")
    parser.add_argument("--ambiguity-fraction", type=float)
    parser.add_argument("--min-overlap", type=int)
    parser.add_argument("--quiet", action="store_true", help="Pass -q to TN93")
    parser.add_argument("--extra-tn93-args", default="")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.set_defaults(func=command_all)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hcv-workflow",
        description="Unified HCV genotype, preprocessing, TN93 distance, linking, and clustering workflow.",
    )
    subparsers = parser.add_subparsers(dest="command")
    subparsers.add_parser("genotype", help="Assign HCV genotypes using the existing minimap2 workflow")
    subparsers.add_parser("prep", help="Run the existing genotype-specific prep-align workflow")
    subparsers.add_parser("prep-refs", help="Download/cache reference records for preprocessing")
    add_tn93_parser(subparsers)
    add_link_parser(subparsers)
    add_cluster_parser(subparsers)
    add_all_parser(subparsers)
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help"}:
        parser = build_parser()
        parser.print_help()
        return 0
    if argv[0] == "genotype":
        return assign_hcv_genotypes_from_fasta.main(argv[1:])
    if argv[0] == "prep":
        return hcv_cluster_prep.main(["prep-align", *argv[1:]])
    if argv[0] == "prep-refs":
        return hcv_cluster_prep.main(["prep-refs", *argv[1:]])

    parser = build_parser()
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
