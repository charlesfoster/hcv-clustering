#!/usr/bin/env python3
"""Handle repeated FASTA identifiers (e.g. mixed infections) for hcv-cluster runs.

A sample may appear several times in the input FASTA under one identifier. Exact
repeats are dropped. Distinct sequences are renamed ``SAMPLE__v1``, ``SAMPLE__v2``,
... so every downstream step keeps a unique key, and each renamed sequence remembers
its original sample identifier for the metadata join.

Variants assigned to different genotypes are all kept (they land in different
networks). Among variants sharing a genotype only one is kept, chosen by either:

* ``region-coverage``: best coverage of the selected clustering region, then the
  most complete sequence when coverage ties;
* ``completeness``: the most complete sequence (most unambiguous A/C/G/T bases).

Remaining ties fall to input order.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping

import assign_hcv_genotypes_from_fasta
import hcv_cluster_prep

VARIANT_SEPARATOR = "__v"
SELECTION_MODES = ("region-coverage", "completeness")
BASE_SAMPLE_FIELD = "base_sample_id"
MIXED_INFECTION_FIELD = "mixed_infection"
VARIANTS_FILENAME = "sequence_variants.csv"
VARIANT_FIELDS = (
    "sequence_id",
    BASE_SAMPLE_FIELD,
    "variant",
    "n_variants",
    "input_record",
    "completeness",
    "genotype",
    "region_coverage",
    "status",
    "detail",
)


@dataclass
class VariantRecord:
    sequence_id: str
    sample_id: str
    variant: str
    n_variants: int
    input_record: int
    completeness: int
    genotype: str = ""
    region_coverage: float | None = None
    status: str = "kept"
    detail: str = ""

    def as_row(self) -> dict[str, str | int]:
        return {
            "sequence_id": self.sequence_id,
            BASE_SAMPLE_FIELD: self.sample_id,
            "variant": self.variant,
            "n_variants": self.n_variants,
            "input_record": self.input_record,
            "completeness": self.completeness,
            "genotype": self.genotype,
            "region_coverage": "" if self.region_coverage is None else f"{self.region_coverage:.6f}",
            "status": self.status,
            "detail": self.detail,
        }


@dataclass
class VariantPlan:
    """Result of resolving repeated identifiers in one input FASTA."""

    records: list[hcv_cluster_prep.FastaRecord]
    variants: list[VariantRecord] = field(default_factory=list)

    @property
    def has_duplicates(self) -> bool:
        return bool(self.variants)

    def base_sample_ids(self) -> dict[str, str]:
        return {variant.sequence_id: variant.sample_id for variant in self.variants}

    def mixed_sample_ids(self) -> set[str]:
        return {variant.sample_id for variant in self.variants if variant.n_variants > 1}


def _normalized_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper().replace("U", "T")


def completeness(sequence: str) -> int:
    """Count unambiguous nucleotides; gaps, Ns and IUPAC ambiguity codes do not count."""
    return sum(1 for base in _normalized_sequence(sequence) if base in "ACGT")


def read_input_records(path: Path) -> list[hcv_cluster_prep.FastaRecord]:
    """Read FASTA records (plain or .gz) without rejecting repeated identifiers."""
    records: list[hcv_cluster_prep.FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []
    with assign_hcv_genotypes_from_fasta.open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(hcv_cluster_prep.FastaRecord(header, "".join(chunks)))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Blank FASTA header at line {line_number}: {path}")
                chunks = []
            elif header is None:
                raise ValueError(f"Sequence data found before first FASTA header at line {line_number}: {path}")
            else:
                chunks.append(line)
    if header is not None:
        records.append(hcv_cluster_prep.FastaRecord(header, "".join(chunks)))
    return records


def resolve_duplicate_ids(records: list[hcv_cluster_prep.FastaRecord]) -> VariantPlan:
    """Drop exact repeats and give distinct same-ID sequences unique variant names."""
    by_id: dict[str, list[tuple[int, hcv_cluster_prep.FastaRecord]]] = {}
    for index, record in enumerate(records, start=1):
        by_id.setdefault(record.header.split()[0], []).append((index, record))
    if all(len(group) == 1 for group in by_id.values()):
        return VariantPlan(records=list(records))

    existing_ids = set(by_id)
    output: list[hcv_cluster_prep.FastaRecord] = []
    variants: list[VariantRecord] = []
    for sample_id, group in by_id.items():
        if len(group) == 1:
            continue
        distinct: dict[str, int] = {}
        for index, record in group:
            distinct.setdefault(_normalized_sequence(record.sequence), index)
        n_variants = len(distinct)
        variant_by_index = {index: number for number, index in enumerate(distinct.values(), start=1)}
        for index, record in group:
            key = _normalized_sequence(record.sequence)
            first_index = distinct[key]
            if index != first_index:
                variants.append(
                    VariantRecord(
                        sequence_id="",
                        sample_id=sample_id,
                        variant="",
                        n_variants=n_variants,
                        input_record=index,
                        completeness=completeness(record.sequence),
                        status="dropped_identical_repeat",
                        detail=f"identical to input record {first_index}",
                    )
                )
                continue
            if n_variants == 1:
                sequence_id, variant = sample_id, ""
            else:
                variant = f"v{variant_by_index[index]}"
                sequence_id = f"{sample_id}{VARIANT_SEPARATOR}{variant_by_index[index]}"
                if sequence_id in existing_ids:
                    raise ValueError(
                        f"Cannot rename repeated FASTA ID '{sample_id}' to '{sequence_id}': "
                        "that ID already exists in the input"
                    )
            variants.append(
                VariantRecord(
                    sequence_id=sequence_id,
                    sample_id=sample_id,
                    variant=variant,
                    n_variants=n_variants,
                    input_record=index,
                    completeness=completeness(record.sequence),
                )
            )

    kept_names = {
        variant.input_record: variant.sequence_id for variant in variants if variant.status == "kept"
    }
    dropped = {variant.input_record for variant in variants if variant.status != "kept"}
    for index, record in enumerate(records, start=1):
        if index in dropped:
            continue
        name = kept_names.get(index)
        output.append(hcv_cluster_prep.FastaRecord(name, record.sequence) if name else record)
    return VariantPlan(records=output, variants=variants)


def record_genotypes(plan: VariantPlan, genotype_rows: Iterable[Mapping[str, str]]) -> None:
    by_sequence = {(row.get("sequence_id") or "").strip(): row for row in genotype_rows}
    for variant in plan.variants:
        if variant.status != "kept":
            continue
        row = by_sequence.get(variant.sequence_id, {})
        variant.genotype = (row.get("assigned_genotype") or "").strip()
        if (row.get("assignment_status") or "").strip() != "pass":
            variant.status = "genotype_not_assigned"
            variant.detail = (row.get("qc_fail_reason") or "").strip()


def _read_region_coverage(qc_path: Path) -> dict[str, float]:
    coverage: dict[str, float] = {}
    if not qc_path.exists():
        return coverage
    with qc_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                coverage[row["sample_id"]] = float(row.get("coverage_fraction") or 0)
            except ValueError:
                coverage[row["sample_id"]] = 0.0
    return coverage


def select_same_genotype(
    plan: VariantPlan,
    genotype: str,
    qc_path: Path,
    mode: str,
) -> list[VariantRecord]:
    """Keep one variant per sample within ``genotype``; return the variants dropped."""
    if mode not in SELECTION_MODES:
        raise ValueError(f"Unknown duplicate selection mode: {mode}")
    coverage = _read_region_coverage(qc_path)
    groups: dict[str, list[VariantRecord]] = {}
    for variant in plan.variants:
        if variant.status == "kept" and variant.genotype == genotype:
            variant.region_coverage = coverage.get(variant.sequence_id, 0.0)
            groups.setdefault(variant.sample_id, []).append(variant)

    dropped: list[VariantRecord] = []
    for group in groups.values():
        if len(group) < 2:
            continue
        if mode == "region-coverage":
            # Coverage is compared at the precision qc.csv reports, so equal coverage ties.
            ranked = sorted(
                group,
                key=lambda v: (-round(v.region_coverage or 0.0, 6), -v.completeness, v.input_record),
            )
        else:
            ranked = sorted(group, key=lambda v: (-v.completeness, v.input_record))
        winner = ranked[0]
        for loser in ranked[1:]:
            loser.status = "dropped_same_genotype"
            loser.detail = f"{mode}: {winner.sequence_id} kept for genotype {genotype}"
            dropped.append(loser)
    return dropped


def remove_from_fasta(path: Path, sequence_ids: set[str]) -> None:
    if not sequence_ids or not path.exists():
        return
    records = hcv_cluster_prep.read_fasta_records(path)
    kept = [record for record in records if record.header.split()[0] not in sequence_ids]
    hcv_cluster_prep.write_fasta_records(kept, path)


def annotate_cluster_table(path: Path, plan: VariantPlan) -> None:
    """Add base_sample_id and mixed_infection columns to a clusters CSV in place."""
    if not path.exists():
        return
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or ())
        rows = list(reader)
    base_ids = plan.base_sample_ids()
    mixed = plan.mixed_sample_ids()
    for extra in (BASE_SAMPLE_FIELD, MIXED_INFECTION_FIELD):
        if extra not in fieldnames:
            insert_at = fieldnames.index("sample_id") + 1 if extra == BASE_SAMPLE_FIELD else len(fieldnames)
            fieldnames.insert(insert_at, extra)
    for row in rows:
        base = base_ids.get(row["sample_id"], row["sample_id"])
        row[BASE_SAMPLE_FIELD] = base
        row[MIXED_INFECTION_FIELD] = "true" if base in mixed else "false"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_variants_csv(plan: VariantPlan, path: Path) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=VARIANT_FIELDS)
        writer.writeheader()
        for variant in sorted(plan.variants, key=lambda v: v.input_record):
            writer.writerow(variant.as_row())
    return path
