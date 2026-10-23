#!/usr/bin/env python3
"""Sample-level audit trail for hcv-cluster runs.

Every record in the input FASTA gets one row in ``sample_audit.csv`` recording the
outcome of each pipeline stage and, for samples that were excluded, which stage
removed them and why. The same classification drives the metadata tables and the
labelled FASTA of surviving samples, so the three can never disagree.

Stages, in the order a sequence meets them:

* ``input``: exact repeats of an earlier record with the same FASTA ID are dropped;
* ``genotyping``: the minimap2 panel assignment must pass its coverage/identity QC;
* ``subtype``: the assigned subtype must be requested (``--genotype``) and have a
  configured reference genome;
* ``genotype_validation``: a k-mer check must not find the sequence closer to a different
  subtype than the one assigned;
* ``region_qc``: the clustering region must be sufficiently covered after alignment;
* ``variant_selection``: among same-genotype variants of one sample, only one is kept;
* ``clustering``: the survivor must appear in the merged cluster table.

A failed stage marks every later stage ``not_run``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from . import metadata, prep, variants

AUDIT_FILENAME = "sample_audit.csv"
METADATA_ALL_FILENAME = "metadata.all_samples.csv"
METADATA_SURVIVING_FILENAME = "metadata.surviving_samples.csv"
SURVIVING_FASTA_FILENAME = "surviving_samples.fasta"
LABEL_SEPARATOR = "|"

PASS = "pass"
FAIL = "fail"
NOT_RUN = "not_run"

STAGES = (
    "input",
    "genotyping",
    "subtype",
    "genotype_validation",
    "region_qc",
    "variant_selection",
    "clustering",
)

AUDIT_FIELDS = (
    "sample_id",
    variants.BASE_SAMPLE_FIELD,
    variants.MIXED_INFECTION_FIELD,
    "input_record",
    "overall_status",
    "failed_stage",
    "overall_reason",
    "genotype",
    "cluster_id",
    "cluster_label",
    "cluster_size",
    *(f"{stage}_{suffix}" for stage in STAGES for suffix in ("status", "reason")),
)

# Pipeline-derived columns prepended to the user's metadata in the metadata tables.
METADATA_PIPELINE_FIELDS = (
    "sample_id",
    variants.BASE_SAMPLE_FIELD,
    variants.MIXED_INFECTION_FIELD,
    "genotype",
    "overall_status",
    "failed_stage",
    "cluster_id",
    "cluster_label",
    "cluster_size",
)

_GENOTYPING_FAIL_TEXT = {
    "no_alignment": "no alignment to any sequence in the genotype reference panel",
}


@dataclass
class SampleAudit:
    sample_id: str
    base_sample_id: str
    mixed_infection: bool
    input_record: int
    genotype: str = ""
    cluster_id: str = ""
    # None when the run has no label registry; "" for a singleton, which is never labelled.
    cluster_label: str | None = None
    cluster_size: str = ""
    is_repeat: bool = False

    def __post_init__(self) -> None:
        self.stages: dict[str, tuple[str, str]] = {stage: (NOT_RUN, "") for stage in STAGES}

    @property
    def failed_stage(self) -> str:
        return next((stage for stage in STAGES if self.stages[stage][0] == FAIL), "")

    @property
    def suffix_label(self) -> str:
        """Cluster name for FASTA headers: persistent label when tracked, else run-local ID."""
        if self.cluster_label is None:
            return self.cluster_id
        return self.cluster_label or "singleton"

    @property
    def survived(self) -> bool:
        return not self.failed_stage

    def record(self, stage: str, passed: bool, reason: str) -> bool:
        """Record a stage outcome; return whether the sample continues to the next stage."""
        self.stages[stage] = (PASS if passed else FAIL, reason)
        return passed

    def as_row(self) -> dict[str, str | int]:
        failed = self.failed_stage
        row: dict[str, str | int] = {
            "sample_id": self.sample_id,
            variants.BASE_SAMPLE_FIELD: self.base_sample_id,
            variants.MIXED_INFECTION_FIELD: "true" if self.mixed_infection else "false",
            "input_record": self.input_record,
            "overall_status": FAIL if failed else PASS,
            "failed_stage": failed,
            "overall_reason": self.stages[failed][1] if failed else "Passed all stages",
            "genotype": self.genotype,
            "cluster_id": self.cluster_id,
            "cluster_label": self.cluster_label or "",
            "cluster_size": self.cluster_size,
        }
        for stage in STAGES:
            status, reason = self.stages[stage]
            row[f"{stage}_status"] = status
            row[f"{stage}_reason"] = reason
        return row


def _genotyping_outcome(
    row: Mapping[str, str] | None, min_query_coverage: float, min_identity: float
) -> tuple[bool, str]:
    if row is None:
        return False, "no genotype assignment was produced for this sequence"
    genotype = (row.get("assigned_genotype") or "").strip()
    if (row.get("assignment_status") or "").strip() == "pass" and genotype:
        return True, (
            f"Assigned genotype {genotype}; callable query coverage "
            f"{row.get('query_coverage', '')} >= {min_query_coverage:g}; "
            f"identity {row.get('identity', '')} >= {min_identity:g}"
        )
    parts: list[str] = []
    for code in (row.get("qc_fail_reason") or "").split(";"):
        code = code.strip()
        if code == "low_query_coverage":
            parts.append(
                f"callable query coverage {row.get('query_coverage', '')} < {min_query_coverage:g}"
            )
        elif code == "low_identity":
            parts.append(f"identity {row.get('identity', '')} < {min_identity:g}")
        elif code:
            parts.append(_GENOTYPING_FAIL_TEXT.get(code, code))
    return False, "; ".join(parts) or "genotype assignment failed"


def _region_qc_outcome(row: Mapping[str, str], min_coverage: float) -> tuple[bool, str]:
    passed = (row.get("passed_qc") or "").strip() == "true"
    coverage = float(row.get("coverage_fraction") or 0)
    comparator = ">=" if passed else "<"
    return passed, (
        f"{'Passed' if passed else 'Failed'}: {row.get('selected_region', 'region')} coverage "
        f"{coverage:.3f} {comparator} {min_coverage:g}; "
        f"{row.get('non_gap_bases', '?')}/{row.get('region_length', '?')} callable bases; "
        f"{row.get('n_bases', '?')} N; {row.get('gap_bases', '?')} gaps"
    )


def read_qc_rows(qc_path: Path) -> dict[str, dict[str, str]]:
    if not qc_path.exists():
        return {}
    with qc_path.open(newline="", encoding="utf-8") as handle:
        return {row["sample_id"]: row for row in csv.DictReader(handle)}


def build_audit(
    input_records: Sequence[prep.FastaRecord],
    variant_plan: variants.VariantPlan,
    genotype_rows: Iterable[Mapping[str, str]],
    *,
    min_query_coverage: float,
    min_identity: float,
    min_coverage: float,
    excluded_subtypes: Mapping[str, str] | None = None,
    qc_by_genotype: Mapping[str, Mapping[str, Mapping[str, str]]] | None = None,
    validation_failures: Mapping[str, Mapping[str, str]] | None = None,
    cluster_rows: Iterable[Mapping[str, str]] = (),
) -> list[SampleAudit]:
    """Classify every input record against each pipeline stage.

    ``excluded_subtypes`` maps a subtype that passed genotyping but was not analysed to
    the reason. ``qc_by_genotype`` maps a subtype to its parsed ``prep.qc.csv`` rows, and
    ``validation_failures`` maps a subtype to ``{sequence_id: reason}`` from genotype validation.
    """
    excluded_subtypes = excluded_subtypes or {}
    qc_by_genotype = qc_by_genotype or {}
    validation_failures = validation_failures or {}
    genotype_by_sequence = {(row.get("sequence_id") or "").strip(): row for row in genotype_rows}
    clusters = {row["sample_id"]: row for row in cluster_rows}
    mixed = variant_plan.mixed_sample_ids()
    variants_by_record = {variant.input_record: variant for variant in variant_plan.variants}

    audits: list[SampleAudit] = []
    for index, record in enumerate(input_records, start=1):
        original_id = record.header.split()[0]
        variant = variants_by_record.get(index)
        sequence_id = variant.sequence_id if variant is not None else original_id
        audit = SampleAudit(
            sample_id=sequence_id or original_id,
            base_sample_id=original_id,
            mixed_infection=original_id in mixed,
            input_record=index,
        )
        audits.append(audit)

        if variant is not None and variant.status == "dropped_identical_repeat":
            audit.is_repeat = True
            audit.record("input", False, f"Redundant copy: {variant.detail}; the first copy was analysed")
            continue
        audit.record(
            "input",
            True,
            f"Distinct sequence kept as {sequence_id} ({variant.n_variants} sequences share this ID)"
            if variant is not None and variant.variant
            else "Unique FASTA ID",
        )

        row = genotype_by_sequence.get(sequence_id)
        passed, reason = _genotyping_outcome(row, min_query_coverage, min_identity)
        audit.genotype = (row.get("assigned_genotype") or "").strip() if row else ""
        if not audit.record("genotyping", passed, reason):
            continue

        genotype = audit.genotype
        if genotype in excluded_subtypes:
            audit.record("subtype", False, excluded_subtypes[genotype])
            continue
        audit.record("subtype", True, f"Subtype {genotype} analysed")

        conflict = validation_failures.get(genotype, {}).get(sequence_id)
        if conflict is not None:
            audit.record("genotype_validation", False, conflict)
            continue
        audit.record("genotype_validation", True, "No conflicting subtype found")

        qc_row = qc_by_genotype.get(genotype, {}).get(sequence_id)
        if qc_row is None:
            audit.record("region_qc", False, "no region QC result was produced for this sequence")
            continue
        passed, reason = _region_qc_outcome(qc_row, min_coverage)
        if not audit.record("region_qc", passed, reason):
            continue

        if variant is not None and variant.status == "dropped_same_genotype":
            audit.record("variant_selection", False, variant.detail)
            continue
        audit.record(
            "variant_selection",
            True,
            "Only sequence kept for this sample and genotype"
            if variant is not None and variant.n_variants > 1
            else "No competing variants",
        )

        cluster = clusters.get(sequence_id)
        if cluster is None:
            # Only reachable before clustering has run (e.g. an early exit).
            audit.record("clustering", False, "not present in the cluster table")
            continue
        audit.cluster_id = cluster.get("cluster_id", "")
        audit.cluster_label = cluster.get("cluster_label")
        audit.cluster_size = cluster.get("cluster_size", "")
        size = int(audit.cluster_size or 0)
        sequences = int(cluster.get("cluster_sequences") or size)
        name = audit.cluster_label or audit.cluster_id
        if size > 1:
            reason = (
                f"Member of {name} ({size} patients, {sequences} sequences)"
                if sequences != size
                else f"Member of {name} ({size} sequences)"
            )
        elif sequences > 1:
            reason = f"Not a cluster: all {sequences} linked sequences are from one patient"
        else:
            reason = "Singleton (no sequence within threshold)"
        audit.record("clustering", True, reason)
    return audits


def write_audit_csv(audits: Sequence[SampleAudit], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDS)
        writer.writeheader()
        for audit in audits:
            writer.writerow(audit.as_row())
    return path


def summarize(audits: Sequence[SampleAudit]) -> str:
    """One-line description of how many records survived and where the rest were lost."""
    failures: dict[str, int] = {}
    for audit in audits:
        if audit.failed_stage:
            failures[audit.failed_stage] = failures.get(audit.failed_stage, 0) + 1
    survived = sum(1 for audit in audits if audit.survived)
    text = f"{survived}/{len(audits)} input records survived"
    if failures:
        text += "; excluded at " + ", ".join(
            f"{stage} ({failures[stage]})" for stage in STAGES if stage in failures
        )
    return text


def write_metadata_tables(
    audits: Sequence[SampleAudit],
    metadata_rows: Sequence[Mapping[str, Any]],
    all_path: Path,
    surviving_path: Path,
) -> tuple[Path, Path]:
    """Write metadata for every analysed sequence, and for survivors only.

    Rows are keyed by pipeline sequence ID (renamed mixed-infection variants keep their
    original ID in ``base_sample_id``, which is what the user's metadata is joined on).
    Redundant exact repeats are not distinct sequences and appear only in the audit.
    """
    nodes = [
        {
            "sample_id": audit.sample_id,
            variants.BASE_SAMPLE_FIELD: audit.base_sample_id,
            variants.MIXED_INFECTION_FIELD: "true" if audit.mixed_infection else "false",
            "genotype": audit.genotype,
            "overall_status": PASS if audit.survived else FAIL,
            "failed_stage": audit.failed_stage,
            "cluster_id": audit.cluster_id,
            "cluster_label": audit.cluster_label or "",
            "cluster_size": audit.cluster_size,
        }
        for audit in audits
        if not audit.is_repeat
    ]
    joined = (
        metadata.join_metadata(nodes, metadata_rows)[0] if metadata_rows else nodes
    )
    fields = (*METADATA_PIPELINE_FIELDS, *metadata.metadata_fields(metadata_rows))
    for path, rows in (
        (all_path, joined),
        (surviving_path, [row for row in joined if row["overall_status"] == PASS]),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    return all_path, surviving_path


def labelled_header(header: str, genotype: str, cluster_id: str) -> str:
    """Append ``|genotype|cluster`` to the ID token, keeping any trailing description."""
    identifier, _, description = header.partition(" ")
    labelled = LABEL_SEPARATOR.join((identifier, genotype, cluster_id))
    return f"{labelled} {description}".rstrip()


def write_surviving_fasta(
    audits: Sequence[SampleAudit],
    records: Sequence[prep.FastaRecord],
    path: Path,
) -> Path:
    """Write the original (unaligned) sequences of survivors with genotype+cluster labels.

    ``records`` is the variant-resolved record list, whose IDs match the audit's
    ``sample_id``.
    """
    by_id = {audit.sample_id: audit for audit in audits if audit.survived and not audit.is_repeat}
    labelled = [
        prep.FastaRecord(
            labelled_header(
                record.header,
                by_id[record.header.split()[0]].genotype,
                by_id[record.header.split()[0]].suffix_label,
            ),
            record.sequence,
        )
        for record in records
        if record.header.split()[0] in by_id
    ]
    prep.write_fasta_records(labelled, path)
    return path
