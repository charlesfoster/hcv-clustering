#!/usr/bin/env python3
"""Update the master FASTA from the output of an h2seq run.

Standard library only. Takes the h2seq results folder and the samplesheet that was given to that run
(its ``sample`` column lists the samples), and updates ``h2seq_sequences_master.fasta``:

* each sample's IUPAC consensus is added with the sample ID as the FASTA header (``>SAMPLE``). A mixed infection has several consensus files (main, alt1, ...); each is added under the same ID.
* a sample in the samplesheet that h2seq reported but produced no consensus for (no usable reads, no
  mapped reads, ...) gets an all-N record, so the master FASTA has a record for every sample on which
  sequencing was attempted. The clustering pipeline then lists it as excluded, with the reason.
* a sample that is already in the master FASTA is replaced only if the new consensus has MORE called
  (non-N) bases than the existing one; the old records move to the archive FASTA. Otherwise the
  existing records are kept, so a failed re-run never overwrites a good sequence.

Before anything is written, the master FASTA is copied to a timestamped ``.bak`` file, and every
decision is appended to ``<master>.update_log.tsv``. Existing records are never reordered or
rewritten, except the ones replaced as described above.

Example:
    python3 update_master_fasta.py --results ~/runs/2026-11-09 --samplesheet samplesheet.csv \\
        --master h2seq_sequences_master.fasta
    python3 update_master_fasta.py --results ... --samplesheet ... --master ... --dry-run
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

CONSENSUS_DIR = "all_consensus_genomes"
SUMMARY_FILE = "combined_results_summary.csv"
DEFAULT_MASTER = "h2seq_sequences_master.fasta"
DEFAULT_N_LENGTH = 9600  # used for an all-N record only when the run does not give a reference length
LINE_WIDTH = 60
LOG_FIELDS = ("time", "results", "sample", "action", "detail")


@dataclass
class Record:
    sample: str
    component: str  # main, alt1, ...
    read_type: str  # long, short, or n/a
    sequence: str
    placeholder: bool = False

    @property
    def called(self) -> int:
        return sum(1 for base in self.sequence if base not in "Nn-")


@dataclass
class MasterBlock:
    sample: str
    text: str  # the record exactly as it is in the file


def read_samplesheet(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        names = {(name or "").strip().lower(): name for name in reader.fieldnames or []}
        column = names.get("sample") or names.get("sample_id")
        if column is None:
            raise SystemExit(f"error: {path} has no 'sample' column (found: {', '.join(names) or 'none'})")
        samples: list[str] = []
        for row in reader:
            sample = (row[column] or "").strip()
            if sample and sample not in samples:
                samples.append(sample)
    if not samples:
        raise SystemExit(f"error: {path} lists no samples")
    return samples


def read_summary(results: Path) -> dict[str, list[dict[str, str]]]:
    path = results / SUMMARY_FILE
    if not path.is_file():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows: dict[str, list[dict[str, str]]] = {}
        for row in csv.DictReader(handle):
            rows.setdefault((row.get("sample_id") or "").strip(), []).append(row)
    return rows


def read_first_sequence(path: Path) -> str:
    sequence: list[str] = []
    seen_header = False
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if seen_header:
                break
            seen_header = True
        elif seen_header:
            sequence.append(line.strip())
    return "".join(sequence).upper()


def find_consensus(results: Path, sample: str) -> dict[str, list[Record]]:
    """The IUPAC consensus files for one sample, by read type. ``*.simple.consensus_*`` (majority-allele) files are skipped."""
    pattern = re.compile(rf"(?:^|\.){re.escape(sample)}\.(?:.*\.)?consensus_(main|alt\d+)\.fa$")
    found: dict[str, list[Record]] = {}
    for path in sorted((results / CONSENSUS_DIR).glob("*.fa")):
        match = pattern.search(path.name)
        if not match or ".simple." in path.name:
            continue
        read_type = "long" if "long" in path.name else "short" if "short" in path.name else "n/a"
        found.setdefault(read_type, []).append(Record(sample, match.group(1), read_type, read_first_sequence(path)))
    for records in found.values():
        records.sort(key=lambda r: (r.component != "main", r.component))
    return found


def reference_length(rows: list[dict[str, str]], default: int) -> int:
    for row in rows:
        value = (row.get("reference_length") or "").strip()
        if value.isdigit():
            return int(value)
    return default


def records_for(results: Path, sample: str, summary: dict[str, list[dict[str, str]]], default_length: int) -> tuple[list[Record], str]:
    """The records to add for one sample, and a note. Empty when the sample is not in this run at all."""
    by_type = find_consensus(results, sample)
    if by_type:
        note = ""
        if len(by_type) > 1:  # long and short reads both made a consensus: keep the more complete read type
            best = max(by_type, key=lambda t: by_type[t][0].called)
            note = f"consensus from both read types; kept {best}, the more complete"
            by_type = {best: by_type[best]}
        return next(iter(by_type.values())), note
    rows = summary.get(sample)
    if not rows:
        return [], "not in this run's results"
    reasons = ";".join(sorted({(r.get("qc_fail_reason") or r.get("qc_status") or "").strip() for r in rows} - {""}))
    length = reference_length(rows, default_length)
    return [Record(sample, "main", "n/a", "N" * length, placeholder=True)], f"no consensus ({reasons or 'no reason given'}); all-N of {length}"


def read_master(path: Path) -> list[MasterBlock]:
    if not path.is_file():
        return []
    blocks: list[MasterBlock] = []
    sample, lines = None, []
    for line in path.read_text().splitlines(keepends=True):
        if line.startswith(">"):
            if sample is not None:
                blocks.append(MasterBlock(sample, "".join(lines)))
            sample, lines = line[1:].split()[0] if line[1:].split() else "", [line]
        elif sample is not None:
            lines.append(line)
    if sample is not None:
        blocks.append(MasterBlock(sample, "".join(lines)))
    return blocks


def block_called(block: MasterBlock) -> int:
    return sum(1 for line in block.text.splitlines()[1:] for base in line.strip() if base not in "Nn-")


def format_record(record: Record) -> str:
    """The header is the sample ID alone. Where each record came from is in the update log."""
    wrapped = "\n".join(record.sequence[i : i + LINE_WIDTH] for i in range(0, len(record.sequence), LINE_WIDTH))
    return f">{record.sample}\n{wrapped}\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Update the master FASTA from an h2seq run.")
    parser.add_argument("--results", type=Path, required=True, help="h2seq output folder (has all_consensus_genomes/)")
    parser.add_argument("--samplesheet", type=Path, required=True, help="the samplesheet given to that h2seq run")
    parser.add_argument("--master", type=Path, default=Path(DEFAULT_MASTER), help=f"master FASTA (default: {DEFAULT_MASTER})")
    parser.add_argument("--archive", type=Path, default=None, help="FASTA that replaced records move to (default: <master>_archive.fasta)")
    parser.add_argument("--default-length", type=int, default=DEFAULT_N_LENGTH,
                        help=f"length of an all-N record when the run gives no reference length (default: {DEFAULT_N_LENGTH})")
    parser.add_argument("--dry-run", action="store_true", help="say what would change; write nothing")
    args = parser.parse_args(argv)

    if not (args.results / CONSENSUS_DIR).is_dir():
        print(f"error: {args.results} has no {CONSENSUS_DIR}/ folder: is it an h2seq output folder?", file=sys.stderr)
        return 2
    if not args.samplesheet.is_file():
        print(f"error: samplesheet not found: {args.samplesheet}", file=sys.stderr)
        return 2
    samples = read_samplesheet(args.samplesheet)
    summary = read_summary(args.results)
    master = read_master(args.master)

    plan: list[tuple[str, str, str, list[Record]]] = []  # sample, action, detail, records to add
    for sample in samples:
        records, note = records_for(args.results, sample, summary, args.default_length)
        if not records:
            plan.append((sample, "skipped", note, []))
            continue
        existing = [b for b in master if b.sample == sample]
        if not existing:
            plan.append((sample, "added_placeholder" if records[0].placeholder else "added", note, records))
            continue
        old_best, new_best = max(block_called(b) for b in existing), max(r.called for r in records)
        if new_best > old_best:
            plan.append((sample, "replaced", f"{note}; called bases {old_best} -> {new_best}".strip("; "), records))
        else:
            plan.append((sample, "kept_existing", f"existing has {old_best} called bases, new has {new_best}", []))

    counts: dict[str, int] = {}
    for sample, action, detail, _ in plan:
        counts[action] = counts.get(action, 0) + 1
        print(f"{action:18s} {sample}  {detail}".rstrip())
    if all(action == "skipped" for _, action, _, _ in plan):
        print("error: none of the samplesheet's samples is in this run's results: wrong results folder? Nothing changed.", file=sys.stderr)
        return 2
    print("Summary: " + ", ".join(f"{n} {a}" for a, n in sorted(counts.items())))
    changes = [item for item in plan if item[1] in ("added", "added_placeholder", "replaced")]
    if args.dry_run or not changes:
        print("Dry run: nothing written." if args.dry_run else "Nothing to change.")
        return 0

    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    backup = args.master.with_name(f"{args.master.name}.{stamp}.bak")
    had_master = args.master.is_file()
    if had_master:
        shutil.copy2(args.master, backup)
    replaced = {sample for sample, action, _, _ in plan if action == "replaced"}
    archive = args.archive or args.master.with_name(f"{args.master.stem}_archive.fasta")
    if replaced:
        with archive.open("a") as handle:
            handle.writelines(b.text for b in master if b.sample in replaced)
    kept = "".join(b.text if b.text.endswith("\n") else b.text + "\n" for b in master if b.sample not in replaced)
    added = "".join(format_record(r) for _, _, _, records in changes for r in records)
    args.master.parent.mkdir(parents=True, exist_ok=True)
    args.master.write_text(kept + added)

    log = args.master.with_name(f"{args.master.name}.update_log.tsv")
    new_log = not log.is_file()
    with log.open("a", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if new_log:
            writer.writerow(LOG_FIELDS)
        for sample, action, detail, _ in plan:
            writer.writerow([stamp, str(args.results), sample, action, detail])
    print(f"Updated {args.master} ({'backup ' + backup.name if had_master else 'new file'}; log {log.name}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
