#!/usr/bin/env python3
import argparse
import csv
import gzip
import logging
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict
from pathlib import Path

import hcv_cluster_prep


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PANEL = SCRIPT_DIR / "reference_data" / "hcv_references.fasta"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Assign HCV genotype/subtype to each sequence in a FASTA by competitive "
            "minimap2 alignment against an HCV reference panel."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA/FASTA.GZ containing sequences to assign")
    parser.add_argument("-o", "--output-csv", required=True, help="Output CSV with one genotype assignment per sequence")
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory used for optional per-genotype FASTA files and temporary output context",
    )
    parser.add_argument(
        "--panel-fasta",
        default=str(DEFAULT_PANEL),
        help=f"Reference-panel FASTA. Default: {DEFAULT_PANEL}",
    )
    parser.add_argument("--minimap2", default="minimap2", help="minimap2 executable path/name")
    parser.add_argument(
        "--preset",
        default="asm10",
        help="minimap2 preset for query FASTA versus reference panel. Default: asm10",
    )
    parser.add_argument(
        "--extra-minimap2-args",
        default="",
        help="Additional minimap2 arguments as a single shell-like string, e.g. '--secondary=no'",
    )
    parser.add_argument(
        "--max-secondary",
        type=int,
        default=50,
        help="Maximum secondary/reference-panel hits minimap2 may report per sequence. Default: 50",
    )
    parser.add_argument(
        "--min-query-coverage",
        type=float,
        default=0.50,
        help="Minimum best-hit aligned query fraction required for a pass assignment. Default: 0.50",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.75,
        help="Minimum best-hit identity required for a pass assignment. Default: 0.75",
    )
    parser.add_argument(
        "--close-hit-fraction",
        type=float,
        default=0.98,
        help="Report non-winning references with alignment score at least this fraction of the winner. Default: 0.98",
    )
    parser.add_argument(
        "--write-genotype-fastas",
        action="store_true",
        help="Write assigned sequences into per-genotype FASTA files such as hcv_1a.fasta",
    )
    parser.add_argument(
        "--keep-paf",
        action="store_true",
        help="Keep the raw minimap2 PAF beside the output CSV for inspection",
    )
    return parser.parse_args(argv)


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_fasta(path):
    records = OrderedDict()
    current_id = None
    current_description = ""
    chunks = []

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = {
                        "description": current_description,
                        "sequence": "".join(chunks),
                    }
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Blank FASTA header at line {line_number}: {path}")
                current_id = header.split()[0]
                if current_id in records:
                    raise ValueError(f"Duplicate FASTA sequence ID '{current_id}' in {path}")
                current_description = header
                chunks = []
            elif current_id is None:
                raise ValueError(f"Sequence data found before first FASTA header at line {line_number}: {path}")
            else:
                chunks.append(line.strip())

    if current_id is not None:
        records[current_id] = {"description": current_description, "sequence": "".join(chunks)}

    return records


def normalize_records(records):
    """Strip alignment-gap characters and non-IUPAC bases before minimap2 sees the
    sequences, so pre-aligned/exported FASTA input doesn't inflate query_length and
    tank query_coverage. Mirrors hcv_cluster_prep.normalize_input_record."""
    logger = logging.getLogger(__name__)
    normalized = OrderedDict()
    for query_id, record in records.items():
        fasta_record = hcv_cluster_prep.FastaRecord(header=query_id, sequence=record["sequence"])
        try:
            sequence = hcv_cluster_prep.normalize_input_record(fasta_record, logger).sequence
        except hcv_cluster_prep.HcvPrepError:
            # Empty after stripping gaps/whitespace; leave blank so it falls through
            # to the existing "no_alignment" qc_fail path instead of crashing the batch.
            sequence = ""
        normalized[query_id] = {"description": record["description"], "sequence": sequence}
    return normalized


def parse_reference_ids(path):
    reference_ids = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                reference_ids.append(line[1:].strip().split()[0])
    return reference_ids


def extract_genotype(reference_id):
    parts = reference_id.split("_")
    if not parts or not parts[0]:
        raise ValueError(f"Unable to parse genotype/subtype from reference ID: {reference_id}")
    return parts[0]


def genotype_group(genotype):
    return "".join(char for char in genotype if char.isdigit())


def shell_split(value):
    if not value:
        return []
    import shlex

    return shlex.split(value)


def run_minimap2(args, query_fasta_path, paf_path):
    executable = shutil.which(args.minimap2) if not Path(args.minimap2).exists() else args.minimap2
    if executable is None:
        raise RuntimeError(
            f"Could not find minimap2 executable '{args.minimap2}'. Install minimap2 or pass --minimap2 /path/to/minimap2."
        )

    command = [
        str(executable),
        "-x",
        args.preset,
        "-c",
        "-N",
        str(args.max_secondary),
        *shell_split(args.extra_minimap2_args),
        str(args.panel_fasta),
        str(query_fasta_path),
    ]
    with open(paf_path, "w", encoding="utf-8") as paf_handle:
        completed = subprocess.run(command, stdout=paf_handle, stderr=subprocess.PIPE, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "minimap2 failed with exit code "
            f"{completed.returncode}: {completed.stderr.strip() or 'no stderr captured'}"
        )


def parse_optional_tags(fields):
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_paf(path):
    hits_by_query = defaultdict(list)
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF line {line_number} has fewer than 12 fields: {path}")

            tags = parse_optional_tags(fields[12:])
            query_id = fields[0]
            query_length = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            target_id = fields[5]
            target_length = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            matching_bases = int(fields[9])
            alignment_block_length = int(fields[10])
            mapq = int(fields[11])
            alignment_score = int(tags.get("AS", matching_bases))
            primary = tags.get("tp", "P") != "S"

            query_aligned_bases = query_end - query_start
            target_aligned_bases = target_end - target_start
            identity = matching_bases / alignment_block_length if alignment_block_length else 0.0
            query_coverage = query_aligned_bases / query_length if query_length else 0.0
            target_coverage = target_aligned_bases / target_length if target_length else 0.0

            hits_by_query[query_id].append(
                {
                    "query_id": query_id,
                    "query_length": query_length,
                    "target_id": target_id,
                    "target_length": target_length,
                    "query_aligned_bases": query_aligned_bases,
                    "target_aligned_bases": target_aligned_bases,
                    "matching_bases": matching_bases,
                    "alignment_block_length": alignment_block_length,
                    "alignment_score": alignment_score,
                    "mapq": mapq,
                    "identity": identity,
                    "query_coverage": query_coverage,
                    "target_coverage": target_coverage,
                    "primary": primary,
                }
            )
    return hits_by_query


def sort_hits(hits):
    return sorted(
        hits,
        key=lambda hit: (
            not hit["primary"],
            -hit["alignment_score"],
            -hit["matching_bases"],
            -hit["query_aligned_bases"],
            -hit["identity"],
            -hit["mapq"],
            hit["target_id"],
        ),
    )


def fmt_float(value):
    return f"{value:.6f}"


def build_assignment_rows(records, hits_by_query, args):
    rows = []
    for query_id, record in records.items():
        sequence_length = len(record["sequence"])
        hits = sort_hits(hits_by_query.get(query_id, []))
        if not hits:
            rows.append(
                {
                    "sequence_id": query_id,
                    "assigned_genotype": "",
                    "genotype_group": "",
                    "assignment_status": "qc_fail",
                    "qc_fail_reason": "no_alignment",
                    "best_ref": "",
                    "best_ref_genotype": "",
                    "query_length": sequence_length,
                    "best_alignment_score": 0,
                    "second_alignment_score": "",
                    "score_margin": "",
                    "score_margin_fraction": "",
                    "query_coverage": "",
                    "target_coverage": "",
                    "identity": "",
                    "mapq": "",
                    "matching_bases": 0,
                    "alignment_block_length": 0,
                    "close_hits": "",
                    "other_potential_genotypes": "",
                }
            )
            continue

        best = hits[0]
        best_genotype = extract_genotype(best["target_id"])
        second_score = hits[1]["alignment_score"] if len(hits) > 1 else None
        score_margin = best["alignment_score"] - second_score if second_score is not None else None
        score_margin_fraction = score_margin / best["alignment_score"] if score_margin is not None and best["alignment_score"] else None
        close_score = best["alignment_score"] * args.close_hit_fraction
        close_hits = [hit for hit in hits[1:] if hit["alignment_score"] >= close_score]
        other_genotypes = [
            extract_genotype(hit["target_id"])
            for hit in close_hits
            if extract_genotype(hit["target_id"]) != best_genotype
        ]

        fail_reasons = []
        if best["query_coverage"] < args.min_query_coverage:
            fail_reasons.append("low_query_coverage")
        if best["identity"] < args.min_identity:
            fail_reasons.append("low_identity")

        rows.append(
            {
                "sequence_id": query_id,
                "assigned_genotype": "" if fail_reasons else best_genotype,
                "genotype_group": "" if fail_reasons else genotype_group(best_genotype),
                "assignment_status": "qc_fail" if fail_reasons else "pass",
                "qc_fail_reason": ";".join(fail_reasons),
                "best_ref": best["target_id"],
                "best_ref_genotype": best_genotype,
                "query_length": best["query_length"],
                "best_alignment_score": best["alignment_score"],
                "second_alignment_score": "" if second_score is None else second_score,
                "score_margin": "" if score_margin is None else score_margin,
                "score_margin_fraction": "" if score_margin_fraction is None else fmt_float(score_margin_fraction),
                "query_coverage": fmt_float(best["query_coverage"]),
                "target_coverage": fmt_float(best["target_coverage"]),
                "identity": fmt_float(best["identity"]),
                "mapq": best["mapq"],
                "matching_bases": best["matching_bases"],
                "alignment_block_length": best["alignment_block_length"],
                "close_hits": ";".join(hit["target_id"] for hit in close_hits),
                "other_potential_genotypes": ";".join(dict.fromkeys(other_genotypes)),
            }
        )
    return rows


def write_csv(path, rows):
    columns = [
        "sequence_id",
        "assigned_genotype",
        "genotype_group",
        "assignment_status",
        "qc_fail_reason",
        "best_ref",
        "best_ref_genotype",
        "query_length",
        "best_alignment_score",
        "second_alignment_score",
        "score_margin",
        "score_margin_fraction",
        "query_coverage",
        "target_coverage",
        "identity",
        "mapq",
        "matching_bases",
        "alignment_block_length",
        "close_hits",
        "other_potential_genotypes",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def fasta_wrap(sequence, width=80):
    for index in range(0, len(sequence), width):
        yield sequence[index : index + width]


def write_records_fasta(path, records):
    with open(path, "w", encoding="utf-8") as handle:
        for record in records.values():
            handle.write(f">{record['description']}\n")
            for line in fasta_wrap(record["sequence"]):
                handle.write(f"{line}\n")


def write_genotype_fastas(outdir, records, rows):
    by_genotype = defaultdict(list)
    for row in rows:
        genotype = row["assigned_genotype"]
        if genotype:
            by_genotype[genotype].append(row["sequence_id"])

    for genotype, sequence_ids in sorted(by_genotype.items()):
        output_path = outdir / f"hcv_{genotype}.fasta"
        with open(output_path, "w", encoding="utf-8") as handle:
            for sequence_id in sequence_ids:
                record = records[sequence_id]
                handle.write(f">{record['description']}\n")
                for line in fasta_wrap(record["sequence"]):
                    handle.write(f"{line}\n")


def main(argv=None):
    args = parse_args(argv)
    args.input = Path(args.input)
    args.output_csv = Path(args.output_csv)
    args.outdir = Path(args.outdir)
    args.panel_fasta = Path(args.panel_fasta)

    if not args.input.exists():
        raise ValueError(f"Input FASTA does not exist: {args.input}")
    if not args.panel_fasta.exists():
        raise ValueError(f"Reference-panel FASTA does not exist: {args.panel_fasta}")
    if not 0.0 <= args.close_hit_fraction <= 1.0:
        raise ValueError("--close-hit-fraction must be between 0 and 1")

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)

    records = parse_fasta(args.input)
    if not records:
        raise ValueError(f"Input FASTA contains no sequences: {args.input}")
    records = normalize_records(records)
    if not parse_reference_ids(args.panel_fasta):
        raise ValueError(f"Reference-panel FASTA contains no sequences: {args.panel_fasta}")

    keep_paf_path = args.output_csv.with_suffix(".paf")
    with tempfile.TemporaryDirectory(prefix="hcv_genotype_", dir=args.outdir) as tmpdir:
        query_fasta_path = Path(tmpdir) / "normalized_query.fasta"
        write_records_fasta(query_fasta_path, records)
        paf_path = keep_paf_path if args.keep_paf else Path(tmpdir) / "competitive.paf"
        run_minimap2(args, query_fasta_path, paf_path)
        hits_by_query = parse_paf(paf_path)

    rows = build_assignment_rows(records, hits_by_query, args)
    write_csv(args.output_csv, rows)
    if args.write_genotype_fastas:
        write_genotype_fastas(args.outdir, records, rows)


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
