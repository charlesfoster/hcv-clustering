#!/usr/bin/env python3
"""Clustering methods for hcv-cluster: region thresholds and the in-process SNP distance.

Kept apart from the CLI so the stage commands (used by the Nextflow workflow) and
``hcv_cluster run`` share exactly one implementation.
"""

from __future__ import annotations

import csv
from pathlib import Path

from . import prep, workflow


# Explicit per-region defaults; some are provisional operational conventions and
# others are literature-anchored. See docs/threshold_rationale.md. Deliberately NOT genotype-split:
# the HCV Core-E2 literature applies a single threshold across genotypes, and the
# only Australian genotype-specific number available (different region and metric)
# does not support a 1a/1b-vs-rest split, so a genotype split isn't defensible as
# evidence-based.
#
# TN93: "core-e2-nohvr1" is informed by Bartlett et al. 2017's
# Core-early-E2 TN93 pairwise/connected-components network (same metric and
# graph algorithm, but a shorter fragment) — the closest available precedent.
# "core-e2" (HVR1 included) and "ns5b" have no HVR1-inclusive/NS5B TN93 precedent,
# so they reuse Lamoury et al. 2015's uncorrected-p-distance values as a (slightly
# conservative, since TN93 >= p-distance for the same pair) approximation.
REGION_THRESHOLDS_TN93: dict[str, float] = {
    "e1-e2": 0.03,
    "e1-e2-nohvr1": 0.03,
    "core-e2-nohvr1": 0.03,
    "core-e2": 0.045,
    "ns5b": 0.015,
}

# SNP: this pipeline's SNP distance is uncorrected p-distance (fraction of
# ACGT-comparable sites that differ), the metric Lamoury et al. 2015 used (with a
# different missing-data rule). The metric is directly comparable, but the study's
# shorter fragments and ClusterPicker method mean the region thresholds remain
# literature anchors rather than exact validation of this workflow.
REGION_THRESHOLDS_SNP: dict[str, float] = {
    "e1-e2": 0.03,
    "e1-e2-nohvr1": 0.03,
    "core-e2-nohvr1": 0.03,
    "core-e2": 0.045,
    "ns5b": 0.015,
}

# An explicit operational default is not the same as a region-specific validated
# cutoff. E1-E2 uses 0.03 for continuity with collaborator practice and must keep
# the provisional warning until it has been calibrated for this exact region.
EVIDENCE_ANCHORED_REGIONS_TN93 = frozenset({"core-e2-nohvr1", "core-e2", "ns5b"})
EVIDENCE_ANCHORED_REGIONS_SNP = frozenset({"core-e2-nohvr1", "core-e2", "ns5b"})

# Backward-compat aliases: TN93 was the only metric these names covered originally.
REGION_THRESHOLDS: dict[str, float] = REGION_THRESHOLDS_TN93
FALLBACK_THRESHOLD: float = REGION_THRESHOLDS_TN93["e1-e2"]
FALLBACK_THRESHOLD_SNP: float = REGION_THRESHOLDS_SNP["e1-e2"]


def canonical_threshold_region(region: str) -> str:
    """Normalize equivalent region spellings before choosing a default cutoff."""
    canonical = prep.expand_region_expression(region)
    equivalent_expressions = {
        "e1+e2": "e1-e2",
        "e1,e2": "e1-e2",
        "core-e1+e2": "core-e2",
        "core-e1,e2": "core-e2",
        "core+e1+e2": "core-e2",
        "core,e1,e2": "core-e2",
    }
    return equivalent_expressions.get(canonical, canonical)


def resolve_threshold(region: str, user_threshold: float | None, distance: str = "tn93") -> float:
    if user_threshold is not None:
        return user_threshold
    canonical_region = canonical_threshold_region(region)
    if distance == "snp":
        return REGION_THRESHOLDS_SNP.get(canonical_region, FALLBACK_THRESHOLD_SNP)
    return REGION_THRESHOLDS_TN93.get(canonical_region, FALLBACK_THRESHOLD)


def region_threshold_is_evidence_based(region: str, distance: str = "tn93") -> bool:
    canonical_region = canonical_threshold_region(region)
    anchored = EVIDENCE_ANCHORED_REGIONS_SNP if distance == "snp" else EVIDENCE_ANCHORED_REGIONS_TN93
    return canonical_region in anchored


def region_threshold_is_explicit(region: str, distance: str = "tn93") -> bool:
    """Whether the region has a named default rather than the generic fallback."""
    canonical_region = canonical_threshold_region(region)
    table = REGION_THRESHOLDS_SNP if distance == "snp" else REGION_THRESHOLDS_TN93
    return canonical_region in table


def compute_snp_distances_detailed(
    fasta_path: Path,
) -> list[tuple[str, str, int, int, float]]:
    records = prep.read_fasta_records(fasta_path)
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


def compute_snp_distances(fasta_path: Path) -> list[workflow.DistanceRow]:
    detailed = compute_snp_distances_detailed(fasta_path)
    return [
        workflow.DistanceRow(id_i, id_j, normalized)
        for id_i, id_j, _snp_count, _comparable_sites, normalized in detailed
    ]


def write_snp_csv(
    path: Path,
    fasta_path: Path,
    detailed: list[tuple[str, str, int, int, float]] | None = None,
) -> list[workflow.DistanceRow]:
    if detailed is None:
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
        workflow.DistanceRow(id_i, id_j, normalized)
        for id_i, id_j, _snp_count, _comparable_sites, normalized in detailed
    ]


def write_snp_links_csv(
    path: Path,
    detailed: list[tuple[str, str, int, int, float]],
    threshold: float,
    snp_count_threshold: int | None = None,
) -> list[dict[str, str | int | float]]:
    """Write SNP links, carrying the absolute SNP count alongside the p-distance.

    With snp_count_threshold set, pairs are linked on the raw number of differing
    sites ("within 30 SNPs") instead of on p-distance. The count is always written
    either way, so a link can be described in SNPs even when p-distance defined it.
    """
    if snp_count_threshold is not None and snp_count_threshold < 0:
        raise ValueError("--snp-count-threshold must be non-negative")

    links: list[dict[str, str | int | float]] = []
    for source, target, snp_count, comparable_sites, distance in detailed:
        if source == target:
            continue
        if snp_count_threshold is not None:
            linked = snp_count <= snp_count_threshold
        else:
            linked = distance <= threshold
        if linked:
            links.append(
                {
                    "source": source,
                    "target": target,
                    "distance": distance,
                    "snp_count": snp_count,
                    "comparable_sites": comparable_sites,
                }
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        csv_writer = csv.DictWriter(
            handle, fieldnames=("source", "target", "distance", "snp_count", "comparable_sites")
        )
        csv_writer.writeheader()
        for link in links:
            csv_writer.writerow({**link, "distance": f"{float(link['distance']):.10g}"})
    return links


def _snp_count_summary(links: list[dict[str, str | int | float]]) -> str:
    """Describe linked pairs in whole SNPs, which is how people talk about them."""
    if not links:
        return "no linked pairs"
    counts = sorted(int(link["snp_count"]) for link in links)
    sites = sorted(int(link["comparable_sites"]) for link in links)
    middle = counts[len(counts) // 2]
    return (
        f"linked pairs differ by {counts[0]}-{counts[-1]} SNPs (median {middle}) "
        f"over {sites[0]}-{sites[-1]} comparable sites"
    )
