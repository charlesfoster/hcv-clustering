#!/usr/bin/env python3
"""Prepare HCV consensus sequences for genotype-specific MicrobeTrace clustering."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner, substitution_matrices
    from Bio.Seq import Seq as BioSeq
except ImportError:  # pragma: no cover - handled at runtime with a clearer message
    SeqIO = None
    PairwiseAligner = None
    substitution_matrices = None
    BioSeq = None


LOGGER_NAME = "hcv_cluster_prep"
EUTILS_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
IUPAC_DNA = set("ACGTRYSWKMBDHVN")


@dataclass(frozen=True)
class ReferenceSpec:
    genotype: str
    accession: str


@dataclass(frozen=True)
class ReferenceBundle:
    spec: ReferenceSpec
    genbank_path: Path
    fasta_path: Path
    coords_path: Path
    regions_path: Path


@dataclass(frozen=True)
class CoreE2Coords:
    start: int  # 1-based inclusive reference coordinate
    end: int  # 1-based inclusive reference coordinate
    source: str
    detail: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class RegionSegment:
    name: str
    start: int  # 1-based inclusive reference coordinate
    end: int  # 1-based inclusive reference coordinate
    source: str
    detail: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class RegionSelection:
    expression: str
    strategy: str
    segments: tuple[RegionSegment, ...]
    detail: str = ""

    @property
    def length(self) -> int:
        return sum(segment.length for segment in self.segments)

    @property
    def start(self) -> int:
        return min(segment.start for segment in self.segments)

    @property
    def end(self) -> int:
        return max(segment.end for segment in self.segments)

    @property
    def segment_text(self) -> str:
        return ";".join(f"{segment.name}:{segment.start}-{segment.end}" for segment in self.segments)


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


@dataclass(frozen=True)
class CoverageMetrics:
    extracted_length: int
    non_gap_bases: int
    n_bases: int
    gap_bases: int
    coverage_fraction: float
    n_fraction: float
    gap_fraction: float


@dataclass(frozen=True)
class WindowCandidate:
    segment: RegionSegment
    support_count: int
    mean_coverage: float


@dataclass(frozen=True)
class CdsInfo:
    feature: Any
    translation: str


class HcvPrepError(RuntimeError):
    """Raised for user-facing preprocessing failures."""


DEFAULT_REFERENCE_CATALOG = {
    "1a": ReferenceSpec("1a", "EF407457.1"),
    "1b": ReferenceSpec("1b", "EU781828.1"),
    "2a": ReferenceSpec("2a", "D00944.1"),
    "2b": ReferenceSpec("2b", "D10988.1"),
    "3a": ReferenceSpec("3a", "D17763.1"),
}

CANONICAL_REGION_ORDER = (
    "core",
    "e1",
    "e2",
    "p7",
    "ns2",
    "ns3",
    "ns4a",
    "ns4b",
    "ns5a",
    "ns5b",
)
CANONICAL_REGION_INDEX = {name: index for index, name in enumerate(CANONICAL_REGION_ORDER)}
REGION_ALIASES = {
    "c": "core",
    "core": "core",
    "capsid": "core",
    "e1": "e1",
    "e2": "e2",
    "e2/ns1": "e2",
    "ns1": "e2",
    "p7": "p7",
    "ns2": "ns2",
    "ns3": "ns3",
    "ns4a": "ns4a",
    "ns4-a": "ns4a",
    "ns4b": "ns4b",
    "ns4-b": "ns4b",
    "ns5a": "ns5a",
    "ns5-a": "ns5a",
    "ns5b": "ns5b",
    "ns5-b": "ns5b",
}
REGION_PRESETS = {
    "structural": "core-e2-nohvr1",
    "envelope": "e1-e2",
    "nonstructural": "ns2-ns5b",
    "all-cds": "cds",
    "whole-coding": "cds",
    "coding": "cds",
    "polyprotein": "cds",
}

# HVR1 is the N-terminal 27 aa (81 nt) of E2, immediately following the E1/E2
# cleavage site (H77 polyprotein residues 384-410; E2 itself starts at residue
# 384). See docs/threshold_rationale.md for citations. "core-e2-nohvr1" excises
# exactly this span from the contiguous Core-through-E2 region; "core-e2" keeps
# the full, HVR1-inclusive span unchanged.
HVR1_LENGTH_NT = 81

DIRECT_BOUNDARY_FEATURE_TYPES = {
    "mat_peptide",
    "CDS",
    "gene",
    "misc_feature",
    "region",
    "sig_peptide",
}


def setup_logging(log_path: Path | None = None, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger(LOGGER_NAME)
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def require_biopython() -> None:
    if SeqIO is None or BioSeq is None:
        raise HcvPrepError(
            "Biopython is required. Install it with `pixi install`, "
            "`conda install -c conda-forge biopython`, or `pip install biopython`."
        )


def normalize_genotype(genotype: str) -> str:
    return genotype.strip().lower()


def load_reference_catalog(reference_map_path: Path | None) -> dict[str, ReferenceSpec]:
    catalog = dict(DEFAULT_REFERENCE_CATALOG)
    if reference_map_path is None:
        return catalog

    with reference_map_path.open() as handle:
        raw = json.load(handle)

    if not isinstance(raw, dict):
        raise HcvPrepError("--reference-map must be a JSON object keyed by genotype")

    for genotype, value in raw.items():
        normalized = normalize_genotype(genotype)
        if isinstance(value, str):
            accession = value
        elif isinstance(value, dict) and isinstance(value.get("accession"), str):
            accession = value["accession"]
        else:
            raise HcvPrepError(
                "Each --reference-map entry must be either an accession string "
                "or an object with an 'accession' field"
            )
        catalog[normalized] = ReferenceSpec(normalized, accession.strip())

    return catalog


def get_reference_spec(genotype: str, catalog: dict[str, ReferenceSpec]) -> ReferenceSpec:
    normalized = normalize_genotype(genotype)
    if normalized not in catalog:
        allowed = ", ".join(sorted(catalog))
        raise HcvPrepError(f"Unsupported genotype '{genotype}'. Allowed values: {allowed}")
    return catalog[normalized]


def accession_stem(accession: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", accession)


def reference_paths(spec: ReferenceSpec, refs_dir: Path) -> ReferenceBundle:
    stem = accession_stem(spec.accession)
    return ReferenceBundle(
        spec=spec,
        genbank_path=refs_dir / f"{stem}.gb",
        fasta_path=refs_dir / f"{stem}.fasta",
        coords_path=refs_dir / f"{stem}.core_e2.json",
        regions_path=refs_dir / f"{stem}.regions.json",
    )


def fetch_ncbi_record(accession: str, rettype: str, email: str | None) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": rettype,
        "retmode": "text",
        "tool": "hcv_cluster_prep",
    }
    if email:
        params["email"] = email

    url = f"{EUTILS_EFETCH}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "hcv_cluster_prep/0.1 (NCBI E-utilities)"},
    )

    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if attempt < 3:
                time.sleep(2 ** (attempt - 1))

    raise HcvPrepError(f"Failed to download {accession} from NCBI: {last_error}")


def normalize_reference_sequence(sequence: str) -> str:
    return re.sub(r"\s+", "", sequence).upper().replace("U", "T")


def download_reference(
    spec: ReferenceSpec,
    refs_dir: Path,
    force: bool,
    email: str | None,
    logger: logging.Logger,
) -> ReferenceBundle:
    require_biopython()
    refs_dir.mkdir(parents=True, exist_ok=True)
    bundle = reference_paths(spec, refs_dir)

    if force or not bundle.genbank_path.exists():
        logger.info("Downloading GenBank record %s for genotype %s", spec.accession, spec.genotype)
        text = fetch_ncbi_record(spec.accession, "gbwithparts", email)
        if not text.lstrip().startswith("LOCUS"):
            text = fetch_ncbi_record(spec.accession, "gb", email)
        if not text.lstrip().startswith("LOCUS"):
            raise HcvPrepError(f"NCBI did not return a GenBank record for {spec.accession}")
        bundle.genbank_path.write_text(text)
    else:
        logger.info("Using cached GenBank record %s", bundle.genbank_path)

    try:
        record = SeqIO.read(str(bundle.genbank_path), "genbank")
    except Exception as exc:  # noqa: BLE001
        raise HcvPrepError(f"Could not parse cached GenBank record {bundle.genbank_path}: {exc}") from exc

    if force or not bundle.fasta_path.exists():
        write_fasta_records(
            [FastaRecord(spec.accession, normalize_reference_sequence(str(record.seq)))],
            bundle.fasta_path,
        )
        logger.info("Wrote reference FASTA %s", bundle.fasta_path)
    else:
        logger.info("Using cached reference FASTA %s", bundle.fasta_path)

    return bundle


def write_reference_region_cache(
    bundle: ReferenceBundle,
    regions: dict[str, RegionSegment],
    logger: logging.Logger,
) -> None:
    core_e2 = resolve_region_selection("core-e2", regions, strategy="fixed")
    core_e2_segment = core_e2.segments[0]
    payload = {
        "genotype": bundle.spec.genotype,
        "reference_accession": bundle.spec.accession,
        "core_e2_start": core_e2_segment.start,
        "core_e2_end": core_e2_segment.end,
        "core_e2_length": core_e2_segment.length,
        "coordinate_source": core_e2_segment.source,
        "coordinate_detail": core_e2_segment.detail,
    }
    bundle.coords_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    regions_payload = {
        "genotype": bundle.spec.genotype,
        "reference_accession": bundle.spec.accession,
        "regions": {
            name: {
                "start": segment.start,
                "end": segment.end,
                "length": segment.length,
                "source": segment.source,
                "detail": segment.detail,
            }
            for name, segment in sorted(regions.items(), key=lambda item: (item[1].start, item[0]))
        },
    }
    bundle.regions_path.write_text(json.dumps(regions_payload, indent=2, sort_keys=True) + "\n")
    logger.info(
        "Core-E2 coordinates for %s: %s-%s (%s nt; %s)",
        bundle.spec.accession,
        core_e2_segment.start,
        core_e2_segment.end,
        core_e2_segment.length,
        core_e2_segment.source,
    )
    logger.info("Wrote region metadata: %s", bundle.regions_path)


def parse_reference_region_map(
    genbank_path: Path,
    refs_dir: Path,
    catalog: dict[str, ReferenceSpec],
    force_download: bool,
    email: str | None,
    logger: logging.Logger,
) -> dict[str, RegionSegment]:
    """Find HCV protein region bounds from features, then transfer AA bounds if needed.

    The preferred path uses explicit mature peptide/gene/product annotations. If those
    are absent, the fallback aligns the target polyprotein translation to a cached
    reference with direct annotations and maps amino-acid boundaries back through the
    target CDS location. No mature-region nucleotide coordinates are hard-coded.
    """
    require_biopython()
    try:
        record = SeqIO.read(str(genbank_path), "genbank")
    except Exception as exc:  # noqa: BLE001
        raise HcvPrepError(f"Could not parse GenBank record {genbank_path}: {exc}") from exc

    regions: dict[str, RegionSegment] = {}
    cds = cds_segment_from_record(record)
    if cds is not None:
        regions["cds"] = cds

    direct_regions = regions_from_feature_annotations(record)
    regions.update(direct_regions)

    missing_regions = [name for name in CANONICAL_REGION_ORDER if name not in regions]
    if missing_regions:
        if direct_regions:
            logger.info(
                "Direct mature-region annotations in %s are missing: %s; trying amino-acid boundary transfer",
                genbank_path,
                ", ".join(missing_regions),
            )
        else:
            logger.warning(
                "No direct mature-region feature boundaries found in %s; trying amino-acid boundary transfer",
                genbank_path,
            )
        transferred = regions_from_aa_template(
            record=record,
            current_genbank_path=genbank_path,
            refs_dir=refs_dir,
            catalog=catalog,
            force_download=force_download,
            email=email,
            logger=logger,
        )
        for name, segment in transferred.items():
            regions.setdefault(name, segment)

    if "core" not in regions or "e2" not in regions:
        raise HcvPrepError(
            f"Could not determine core-E2 coordinates for {genbank_path}. "
            "The GenBank record lacks direct core/E2 annotations and amino-acid "
            "boundary transfer did not find a usable annotated template."
        )
    return regions


def parse_reference_coords(
    genbank_path: Path,
    refs_dir: Path,
    catalog: dict[str, ReferenceSpec],
    force_download: bool,
    email: str | None,
    logger: logging.Logger,
) -> CoreE2Coords:
    regions = parse_reference_region_map(
        genbank_path=genbank_path,
        refs_dir=refs_dir,
        catalog=catalog,
        force_download=force_download,
        email=email,
        logger=logger,
    )
    selection = resolve_region_selection("core-e2", regions, strategy="fixed")
    segment = selection.segments[0]
    return CoreE2Coords(
        start=segment.start,
        end=segment.end,
        source=segment.source,
        detail=segment.detail,
    )


def feature_text(feature: Any) -> str:
    values = [str(feature.type)]
    for key in ("product", "gene", "note", "standard_name", "label", "region_name", "function"):
        for value in feature.qualifiers.get(key, []):
            values.append(str(value))
    return " ".join(values).lower()


def feature_bounds_0_based(feature: Any) -> tuple[int, int]:
    location = feature.location
    parts = getattr(location, "parts", None) or [location]
    starts = [int(part.start) for part in parts]
    ends = [int(part.end) for part in parts]
    return min(starts), max(ends)


def feature_length(feature: Any) -> int:
    start, end = feature_bounds_0_based(feature)
    return end - start


def is_broad_polyprotein_feature(feature: Any) -> bool:
    text = feature_text(feature)
    return "polyprotein" in text and feature_length(feature) > 2000


def qualifier_values(feature: Any, key: str) -> set[str]:
    return {str(value).strip().lower() for value in feature.qualifiers.get(key, [])}


def qualifier_values_for_keys(feature: Any, keys: Iterable[str]) -> set[str]:
    values: set[str] = set()
    for key in keys:
        values.update(qualifier_values(feature, key))
    return values


def is_direct_boundary_candidate(feature: Any, max_length: int) -> bool:
    return (
        feature.type in DIRECT_BOUNDARY_FEATURE_TYPES
        and feature_length(feature) <= max_length
        and not is_broad_polyprotein_feature(feature)
    )


def is_core_feature(feature: Any) -> bool:
    if not is_direct_boundary_candidate(feature, 1000):
        return False
    text = feature_text(feature)
    if re.search(r"(^|[^a-z0-9])(core|capsid)([^a-z0-9]|$)", text):
        return True
    boundary_names = qualifier_values_for_keys(feature, ("gene", "product", "label", "standard_name"))
    return "c" in boundary_names


def is_e2_feature(feature: Any) -> bool:
    if not is_direct_boundary_candidate(feature, 2500):
        return False
    text = feature_text(feature)
    if re.search(r"(^|[^a-z0-9])e2([^a-z0-9]|$)", text):
        return True
    if "envelope glycoprotein 2" in text or "envelope protein 2" in text:
        return True
    boundary_names = qualifier_values_for_keys(feature, ("gene", "product", "label", "standard_name"))
    return "e2" in boundary_names


def normalize_region_name(name: str) -> str:
    return re.sub(r"\s+", "", name.strip().lower().replace("_", "-"))


def canonical_region_name(name: str) -> str | None:
    normalized = normalize_region_name(name)
    if normalized in REGION_ALIASES:
        return REGION_ALIASES[normalized]
    compact = normalized.replace("-", "")
    return REGION_ALIASES.get(compact)


def canonical_region_from_feature(feature: Any) -> str | None:
    if not is_direct_boundary_candidate(feature, 2500):
        return None

    names = qualifier_values_for_keys(feature, ("gene", "product", "label", "standard_name", "region_name"))
    for raw_name in sorted(names, key=len):
        exact = canonical_region_name(raw_name)
        if exact is not None:
            return exact

    text = feature_text(feature)
    if re.search(r"(^|[^a-z0-9])(core|capsid)([^a-z0-9]|$)", text):
        return "core"
    for name in ("ns5b", "ns5a", "ns4b", "ns4a", "ns3", "ns2", "p7", "e2", "e1"):
        if re.search(rf"(^|[^a-z0-9]){re.escape(name)}([^a-z0-9]|$)", text):
            return name
    if "envelope glycoprotein 2" in text or "envelope protein 2" in text:
        return "e2"
    if "envelope glycoprotein 1" in text or "envelope protein 1" in text:
        return "e1"
    return None


def feature_summary(feature: Any) -> str:
    start, end = feature_bounds_0_based(feature)
    product = feature.qualifiers.get("product", [""])[0]
    gene = feature.qualifiers.get("gene", [""])[0]
    label = product or gene or feature.type
    return f"{feature.type}:{label}:{start + 1}-{end}"


def regions_from_feature_annotations(record: Any) -> dict[str, RegionSegment]:
    regions: dict[str, RegionSegment] = {}
    for feature in record.features:
        name = canonical_region_from_feature(feature)
        if name is None:
            continue
        start, end = feature_bounds_0_based(feature)
        segment = RegionSegment(
            name=name,
            start=start + 1,
            end=end,
            source="genbank_feature_annotation",
            detail=feature_summary(feature),
        )
        existing = regions.get(name)
        if existing is None or segment.length > existing.length:
            regions[name] = segment
    return regions


def coords_from_feature_annotations(record: Any) -> CoreE2Coords | None:
    regions = regions_from_feature_annotations(record)
    if "core" not in regions or "e2" not in regions:
        return None
    core = regions["core"]
    e2 = regions["e2"]
    if e2.end <= core.start:
        return None

    return CoreE2Coords(
        start=core.start,
        end=e2.end,
        source="genbank_feature_annotation",
        detail=f"core={core.detail}; e2={e2.detail}",
    )


def cds_segment_from_record(record: Any) -> RegionSegment | None:
    cds = find_polyprotein_cds(record)
    if cds is None:
        return None
    start, end = feature_bounds_0_based(cds.feature)
    return RegionSegment(
        name="cds",
        start=start + 1,
        end=end,
        source="polyprotein_cds_feature",
        detail=feature_summary(cds.feature),
    )


def find_polyprotein_cds(record: Any) -> CdsInfo | None:
    candidates: list[tuple[int, Any, str]] = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        translation = get_cds_translation(record, feature)
        if not translation:
            continue
        text = feature_text(feature)
        score = len(translation)
        if "polyprotein" in text:
            score += 100000
        candidates.append((score, feature, translation))

    if not candidates:
        return None
    _, feature, translation = max(candidates, key=lambda item: item[0])
    return CdsInfo(feature=feature, translation=translation)


def get_cds_translation(record: Any, feature: Any) -> str:
    translations = feature.qualifiers.get("translation", [])
    if translations:
        return re.sub(r"\s+", "", str(translations[0])).replace("*", "")

    try:
        cds_seq = str(feature.extract(record.seq)).upper().replace("U", "T")
    except Exception:  # noqa: BLE001
        return ""

    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1
    if codon_start > 0:
        cds_seq = cds_seq[codon_start:]
    if len(cds_seq) < 3 or BioSeq is None:
        return ""
    return str(BioSeq(cds_seq).translate(to_stop=True)).replace("*", "")


def cds_genome_positions(feature: Any) -> list[int]:
    location = feature.location
    parts = list(getattr(location, "parts", None) or [location])
    strand = location.strand
    if strand is None and parts:
        strand = parts[0].strand
    strand = strand or 1

    positions: list[int] = []
    if strand == -1:
        for part in sorted(parts, key=lambda item: int(item.start), reverse=True):
            positions.extend(range(int(part.end) - 1, int(part.start) - 1, -1))
    else:
        for part in sorted(parts, key=lambda item: int(item.start)):
            positions.extend(range(int(part.start), int(part.end)))

    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1
    return positions[codon_start:]


def coords_to_aa_bounds(coords: CoreE2Coords, cds_feature: Any) -> tuple[int, int] | None:
    positions = cds_genome_positions(cds_feature)
    offset_by_position = {position: offset for offset, position in enumerate(positions)}
    start_offset = offset_by_position.get(coords.start - 1)
    end_offset = offset_by_position.get(coords.end - 1)
    if start_offset is None or end_offset is None:
        return None
    return start_offset // 3, (end_offset // 3) + 1


def aa_bounds_to_reference_coords(
    aa_start: int,
    aa_end: int,
    cds_feature: Any,
) -> tuple[int, int] | None:
    positions = cds_genome_positions(cds_feature)
    start_offset = aa_start * 3
    end_offset = (aa_end * 3) - 1
    if start_offset < 0 or end_offset >= len(positions):
        return None
    start_pos = positions[start_offset]
    end_pos = positions[end_offset]
    return min(start_pos, end_pos) + 1, max(start_pos, end_pos) + 1


def regions_from_aa_template(
    record: Any,
    current_genbank_path: Path,
    refs_dir: Path,
    catalog: dict[str, ReferenceSpec],
    force_download: bool,
    email: str | None,
    logger: logging.Logger,
) -> dict[str, RegionSegment]:
    target_cds = find_polyprotein_cds(record)
    if target_cds is None:
        logger.warning("No usable polyprotein CDS found for amino-acid boundary transfer")
        return {}

    ensure_template_genbank_files(
        current_genbank_path=current_genbank_path,
        refs_dir=refs_dir,
        catalog=catalog,
        force_download=force_download,
        email=email,
        logger=logger,
    )

    candidates: list[tuple[int, Path, Any, CdsInfo, dict[str, RegionSegment]]] = []
    for template_path in sorted(refs_dir.glob("*.gb")):
        if template_path.resolve() == current_genbank_path.resolve():
            continue
        try:
            template_record = SeqIO.read(str(template_path), "genbank")
        except Exception:  # noqa: BLE001
            continue

        template_regions = regions_from_feature_annotations(template_record)
        if not template_regions:
            continue

        template_cds = find_polyprotein_cds(template_record)
        if template_cds is None:
            continue

        candidates.append((len(template_regions), template_path, template_record, template_cds, template_regions))

    for _, template_path, _template_record, template_cds, template_regions in sorted(
        candidates,
        key=lambda item: (-item[0], str(item[1])),
    ):
        transferred: dict[str, RegionSegment] = {}
        for name, template_segment in template_regions.items():
            template_aa_bounds = coords_to_aa_bounds(template_segment, template_cds.feature)
            if template_aa_bounds is None:
                continue

            mapped = transfer_aa_bounds(
                template_aa=template_cds.translation,
                target_aa=target_cds.translation,
                template_aa_start=template_aa_bounds[0],
                template_aa_end=template_aa_bounds[1],
            )
            if mapped is None:
                continue

            mapped_coords = aa_bounds_to_reference_coords(mapped[0], mapped[1], target_cds.feature)
            if mapped_coords is None:
                continue

            transferred[name] = RegionSegment(
                name=name,
                start=mapped_coords[0],
                end=mapped_coords[1],
                source="polyprotein_aa_boundary_transfer",
                detail=(
                    f"template={template_path.name}; "
                    f"template_region={template_segment.name}:{template_segment.start}-{template_segment.end}; "
                    f"target_aa={mapped[0] + 1}-{mapped[1]}"
                ),
            )

        if transferred:
            return transferred

    return {}


def coords_from_aa_template(
    record: Any,
    current_genbank_path: Path,
    refs_dir: Path,
    catalog: dict[str, ReferenceSpec],
    force_download: bool,
    email: str | None,
    logger: logging.Logger,
) -> CoreE2Coords | None:
    regions = regions_from_aa_template(
        record=record,
        current_genbank_path=current_genbank_path,
        refs_dir=refs_dir,
        catalog=catalog,
        force_download=force_download,
        email=email,
        logger=logger,
    )
    if "core" not in regions or "e2" not in regions:
        return None
    selection = resolve_region_selection("core-e2", regions, strategy="fixed")
    segment = selection.segments[0]
    return CoreE2Coords(segment.start, segment.end, segment.source, segment.detail)


def ensure_template_genbank_files(
    current_genbank_path: Path,
    refs_dir: Path,
    catalog: dict[str, ReferenceSpec],
    force_download: bool,
    email: str | None,
    logger: logging.Logger,
) -> None:
    if has_annotated_template(current_genbank_path, refs_dir):
        return

    logger.info("Downloading known references as annotation templates for boundary transfer")
    for spec in sorted(catalog.values(), key=lambda item: item.genotype):
        bundle = reference_paths(spec, refs_dir)
        if bundle.genbank_path.resolve() == current_genbank_path.resolve():
            continue
        download_reference(spec, refs_dir, force_download, email, logger)


def has_annotated_template(current_genbank_path: Path, refs_dir: Path) -> bool:
    for path in refs_dir.glob("*.gb"):
        if path.resolve() == current_genbank_path.resolve():
            continue
        try:
            record = SeqIO.read(str(path), "genbank")
        except Exception:  # noqa: BLE001
            continue
        if regions_from_feature_annotations(record):
            return True
    return False


def transfer_aa_bounds(
    template_aa: str,
    target_aa: str,
    template_aa_start: int,
    template_aa_end: int,
) -> tuple[int, int] | None:
    if PairwiseAligner is None:
        return None

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    try:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:  # noqa: BLE001
        aligner.match_score = 1
        aligner.mismatch_score = -1

    alignments = aligner.align(target_aa, template_aa)
    if len(alignments) == 0:
        return None
    alignment = alignments[0]

    target_start = map_template_position_at_or_after(alignment, template_aa_start)
    target_end_last = map_template_position_at_or_before(alignment, template_aa_end - 1)
    if target_start is None or target_end_last is None or target_end_last < target_start:
        return None
    return target_start, target_end_last + 1


def map_template_position_at_or_after(alignment: Any, template_pos: int) -> int | None:
    target_blocks, template_blocks = alignment.aligned
    for target_block, template_block in zip(target_blocks, template_blocks):
        target_start, _ = int(target_block[0]), int(target_block[1])
        template_start, template_end = int(template_block[0]), int(template_block[1])
        if template_start <= template_pos < template_end:
            return target_start + (template_pos - template_start)
        if template_pos < template_start:
            return target_start
    return None


def map_template_position_at_or_before(alignment: Any, template_pos: int) -> int | None:
    target_blocks, template_blocks = alignment.aligned
    last_target_position: int | None = None
    for target_block, template_block in zip(target_blocks, template_blocks):
        target_start, target_end = int(target_block[0]), int(target_block[1])
        template_start, template_end = int(template_block[0]), int(template_block[1])
        if template_start <= template_pos < template_end:
            return target_start + (template_pos - template_start)
        if template_end <= template_pos:
            last_target_position = target_end - 1
        if template_start > template_pos:
            break
    return last_target_position


def resolve_region_selection(
    expression: str,
    regions: dict[str, RegionSegment],
    strategy: str,
) -> RegionSelection:
    expanded = expand_region_expression(expression)
    selected: list[RegionSegment] = []
    for part in split_region_expression(expanded):
        selected.extend(resolve_region_part(part, regions))

    if not selected:
        raise HcvPrepError(f"Region expression '{expression}' did not select any reference bases")

    merged = merge_region_segments(selected)
    return RegionSelection(expression=expanded, strategy=strategy, segments=tuple(merged))


def expand_region_expression(expression: str) -> str:
    normalized = normalize_region_name(expression)
    return REGION_PRESETS.get(normalized, normalized)


def split_region_expression(expression: str) -> list[str]:
    return [part for part in re.split(r"[+,]", expression) if part]


def validate_region_syntax(expression: str) -> str | None:
    """Check that `expression` is a well-formed region expression, without needing
    a loaded reference. Returns None if valid, else a human-readable error message."""
    expanded = expand_region_expression(expression)
    parts = split_region_expression(expanded)
    if not parts:
        return f"Region expression '{expression}' is empty"

    known = ", ".join(CANONICAL_REGION_ORDER)
    presets = ", ".join(sorted(REGION_PRESETS))
    for part in parts:
        if part in {"cds", "polyprotein", "core-e2-nohvr1"} or part in REGION_PRESETS:
            continue
        if "-" in part:
            left, right = part.split("-", 1)
            left_name = canonical_region_name(left)
            right_name = canonical_region_name(right)
            if left_name is not None and right_name is not None:
                if CANONICAL_REGION_INDEX[left_name] > CANONICAL_REGION_INDEX[right_name]:
                    return f"Region range '{part}' is reversed ({left_name} comes after {right_name})"
                continue
        if canonical_region_name(part) is not None:
            continue
        return f"Unknown region '{part}'. Known regions: {known}. Presets: {presets}, cds."
    return None


def resolve_region_part(part: str, regions: dict[str, RegionSegment]) -> list[RegionSegment]:
    preset = REGION_PRESETS.get(part)
    if preset is not None:
        return resolve_region_part(preset, regions)

    if part == "core-e2-nohvr1":
        return resolve_core_e2_nohvr1(regions)

    if part in {"cds", "polyprotein"}:
        return [require_region(part, regions)]

    if "-" in part:
        left, right = part.split("-", 1)
        left_name = canonical_region_name(left)
        right_name = canonical_region_name(right)
        if left_name is not None and right_name is not None:
            return [resolve_region_range(left_name, right_name, regions)]

    name = canonical_region_name(part)
    if name is None:
        raise HcvPrepError(
            f"Unknown region '{part}'. Known names include: {known_region_help(regions)}"
        )
    return [require_region(name, regions)]


def require_region(name: str, regions: dict[str, RegionSegment]) -> RegionSegment:
    if name in regions:
        return regions[name]
    raise HcvPrepError(
        f"Reference annotations do not provide region '{name}'. "
        f"Available regions: {known_region_help(regions)}"
    )


def resolve_region_range(
    left_name: str,
    right_name: str,
    regions: dict[str, RegionSegment],
) -> RegionSegment:
    if CANONICAL_REGION_INDEX[left_name] > CANONICAL_REGION_INDEX[right_name]:
        raise HcvPrepError(f"Region range '{left_name}-{right_name}' is reversed")

    left = require_region(left_name, regions)
    right = require_region(right_name, regions)
    start = min(left.start, right.start)
    end = max(left.end, right.end)
    return RegionSegment(
        name=f"{left_name}-{right_name}",
        start=start,
        end=end,
        source=combine_region_sources((left, right)),
        detail=f"{left_name}={left.start}-{left.end}; {right_name}={right.start}-{right.end}",
    )


def mask_hvr1_from_e2(e2: RegionSegment) -> RegionSegment:
    """E2 with its N-terminal HVR1 (first HVR1_LENGTH_NT nt, immediately following
    the E1/E2 cleavage site) excised, matching Bartlett et al. 2017's Core-early-E2
    definition. See docs/threshold_rationale.md."""
    trimmed_start = e2.start + HVR1_LENGTH_NT
    if trimmed_start > e2.end:
        raise HcvPrepError(
            f"E2 region ({e2.start}-{e2.end}, {e2.length} nt) is shorter than the "
            f"{HVR1_LENGTH_NT} nt HVR1 mask; cannot construct 'core-e2-nohvr1' for this reference"
        )
    return RegionSegment(
        name="e2-nohvr1",
        start=trimmed_start,
        end=e2.end,
        source=e2.source,
        detail=f"{e2.detail}; HVR1 (first {HVR1_LENGTH_NT} nt of E2) masked",
    )


def resolve_core_e2_nohvr1(regions: dict[str, RegionSegment]) -> list[RegionSegment]:
    core_through_e1 = resolve_region_range("core", "e1", regions)
    e2 = require_region("e2", regions)
    return [core_through_e1, mask_hvr1_from_e2(e2)]


def merge_region_segments(segments: Iterable[RegionSegment]) -> list[RegionSegment]:
    sorted_segments = sorted(segments, key=lambda segment: (segment.start, segment.end, segment.name))
    merged: list[RegionSegment] = []
    for segment in sorted_segments:
        if not merged or segment.start > merged[-1].end + 1:
            merged.append(segment)
            continue

        previous = merged[-1]
        merged[-1] = RegionSegment(
            name=f"{previous.name}+{segment.name}",
            start=previous.start,
            end=max(previous.end, segment.end),
            source=combine_region_sources((previous, segment)),
            detail="; ".join(detail for detail in (previous.detail, segment.detail) if detail),
        )
    return merged


def combine_region_sources(segments: Iterable[RegionSegment]) -> str:
    sources = sorted({segment.source for segment in segments if segment.source})
    return "+".join(sources) if sources else "unknown"


def known_region_help(regions: dict[str, RegionSegment]) -> str:
    names = [name for name in CANONICAL_REGION_ORDER if name in regions]
    if "cds" in regions:
        names.append("cds")
    if "core" in regions and "e1" in regions and "e2" in regions:
        names.append("core-e2-nohvr1")
    names.extend(sorted(REGION_PRESETS))
    return ", ".join(dict.fromkeys(names))


def coverage_reason_slug(selection: RegionSelection) -> str:
    if selection.expression == "core-e2" and selection.strategy == "fixed":
        return "low_core_e2_coverage"
    slug = re.sub(r"[^a-z0-9]+", "_", selection.expression.lower()).strip("_")
    if selection.strategy == "max-usable":
        slug = "max_usable" if not slug else f"max_usable_{slug}"
    return f"low_{slug}_coverage" if slug else "low_selected_region_coverage"


def parse_fasta_lines(lines: Iterable[str], source: str) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(FastaRecord(header=header, sequence="".join(chunks)))
            header = line[1:]
            if header == "":
                raise HcvPrepError(f"Empty FASTA header in {source} at line {line_number}")
            chunks = []
        else:
            if header is None:
                raise HcvPrepError(f"Sequence data before first FASTA header in {source} at line {line_number}")
            chunks.append(re.sub(r"\s+", "", line))

    if header is not None:
        records.append(FastaRecord(header=header, sequence="".join(chunks)))
    if not records:
        raise HcvPrepError(f"No FASTA records found in {source}")
    return records


def read_fasta_records(path: Path) -> list[FastaRecord]:
    with path.open() as handle:
        return parse_fasta_lines(handle, str(path))


def parse_fasta_text(text: str, source: str) -> list[FastaRecord]:
    return parse_fasta_lines(text.splitlines(), source)


def write_fasta_records(records: Iterable[FastaRecord], path: Path, line_width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            sequence = record.sequence
            for start in range(0, len(sequence), line_width):
                handle.write(sequence[start : start + line_width] + "\n")


def reject_duplicate_ids(records: list[FastaRecord]) -> None:
    seen: dict[str, int] = {}
    for index, record in enumerate(records, start=1):
        if record.header in seen:
            raise HcvPrepError(
                f"Duplicate FASTA ID '{record.header}' at record {index}; "
                f"previous occurrence was record {seen[record.header]}"
            )
        seen[record.header] = index


def normalize_input_record(record: FastaRecord, logger: logging.Logger) -> FastaRecord:
    raw = re.sub(r"\s+", "", record.sequence).upper().replace("U", "T")
    normalized: list[str] = []
    removed_gaps = 0
    invalid = 0
    for char in raw:
        if char in {"-", "."}:
            removed_gaps += 1
            continue
        if char in IUPAC_DNA:
            normalized.append(char)
        else:
            normalized.append("N")
            invalid += 1

    if not normalized:
        raise HcvPrepError(f"Sequence '{record.header}' is empty after normalization")
    if removed_gaps:
        logger.debug("Removed %s pre-existing gap characters from '%s'", removed_gaps, record.header)
    if invalid:
        logger.warning("Converted %s unsupported characters to N in '%s'", invalid, record.header)
    return FastaRecord(header=record.header, sequence="".join(normalized))


def normalize_input_records(records: list[FastaRecord], logger: logging.Logger) -> list[FastaRecord]:
    return [normalize_input_record(record, logger) for record in records]


HEADER_GENOTYPE_PATTERNS = (
    re.compile(r"(?:^|[|;\s,])(?:genotype|subtype|gt|hcv_gt|hcv-genotype)[:= ]+([0-9][a-z]?)", re.I),
    re.compile(r"(?:^|[|;\s,])GT([0-9][a-z]?)(?:$|[|;\s,])", re.I),
)


def genotype_from_header(header: str) -> str | None:
    for pattern in HEADER_GENOTYPE_PATTERNS:
        match = pattern.search(header)
        if match:
            return normalize_genotype(match.group(1))
    return None


def validate_input_genotypes(
    records: list[FastaRecord],
    requested_genotype: str,
    mode: str,
    kmer_size: int,
    refs_dir: Path,
    catalog: dict[str, ReferenceSpec],
    force_download: bool,
    email: str | None,
    logger: logging.Logger,
) -> None:
    if mode == "none":
        logger.warning("Skipping sequence-level genotype validation")
        return

    header_values = {record.header: genotype_from_header(record.header) for record in records}
    present = {header: value for header, value in header_values.items() if value is not None}
    if present:
        missing = [header for header, value in header_values.items() if value is None]
        mismatches = [
            (header, value)
            for header, value in present.items()
            if value != requested_genotype
        ]
        if missing:
            raise HcvPrepError(
                "Some FASTA headers include genotype metadata but others do not. "
                f"Missing genotype tags: {', '.join(missing[:5])}"
            )
        if mismatches:
            preview = ", ".join(f"{header}={value}" for header, value in mismatches[:5])
            raise HcvPrepError(
                f"FASTA header genotype metadata does not match --genotype {requested_genotype}: {preview}"
            )
        logger.info("Validated genotype %s from FASTA header metadata", requested_genotype)
        return

    if mode == "headers":
        raise HcvPrepError(
            "No genotype metadata found in FASTA headers. Use headers like "
            "'sample genotype=1a' or choose --genotype-validation kmer/auto/none."
        )

    validate_genotypes_by_kmer(
        records=records,
        requested_genotype=requested_genotype,
        kmer_size=kmer_size,
        refs_dir=refs_dir,
        catalog=catalog,
        force_download=force_download,
        email=email,
        logger=logger,
    )


def validate_genotypes_by_kmer(
    records: list[FastaRecord],
    requested_genotype: str,
    kmer_size: int,
    refs_dir: Path,
    catalog: dict[str, ReferenceSpec],
    force_download: bool,
    email: str | None,
    logger: logging.Logger,
) -> None:
    if kmer_size < 8:
        raise HcvPrepError("--genotype-validation-kmer-size must be at least 8")

    reference_kmers: dict[str, set[str]] = {}
    for spec in sorted(catalog.values(), key=lambda item: item.genotype):
        bundle = download_reference(spec, refs_dir, force_download, email, logger)
        ref_records = read_fasta_records(bundle.fasta_path)
        reference_kmers[spec.genotype] = dna_kmers(ref_records[0].sequence, kmer_size)

    failures: list[str] = []
    for record in records:
        query_kmers = dna_kmers(record.sequence, kmer_size)
        if len(query_kmers) < 20:
            failures.append(f"{record.header}: too few unambiguous {kmer_size}-mers")
            continue

        scores = {
            genotype: len(query_kmers & ref_kmers) / len(query_kmers)
            for genotype, ref_kmers in reference_kmers.items()
        }
        best_genotype, best_score = max(scores.items(), key=lambda item: item[1])
        requested_score = scores.get(requested_genotype, 0.0)
        margin = best_score - requested_score
        if best_genotype != requested_genotype and margin >= 0.005 and best_score >= 0.02:
            failures.append(
                f"{record.header}: closest={best_genotype} "
                f"({best_score:.3f}), requested={requested_genotype} ({requested_score:.3f})"
            )

    if failures:
        preview = "; ".join(failures[:10])
        raise HcvPrepError(
            "Genotype validation failed. Do not align mixed genotypes in one run. "
            f"{preview}"
        )

    logger.info("Validated genotype %s by %s-mer comparison to reference anchors", requested_genotype, kmer_size)


def dna_kmers(sequence: str, k: int) -> set[str]:
    seq = normalize_reference_sequence(sequence)
    kmers: set[str] = set()
    for index in range(0, len(seq) - k + 1):
        kmer = seq[index : index + k]
        if set(kmer) <= {"A", "C", "G", "T"}:
            kmers.add(kmer)
    return kmers


def run_mafft_addfragments(
    records: list[FastaRecord],
    reference_fasta: Path,
    reference_accession: str,
    mafft_command: str,
    threads: int,
    keep_temp: bool,
    logger: logging.Logger,
) -> tuple[str, dict[str, str]]:
    mafft_exe = shutil.which(mafft_command)
    if mafft_exe is None:
        raise HcvPrepError(
            f"MAFFT executable '{mafft_command}' was not found on PATH. "
            "Install MAFFT or pass --mafft /path/to/mafft."
        )
    if threads < 1:
        raise HcvPrepError("--threads must be at least 1")

    surrogate_to_header = {
        f"seq_{index:06d}": record.header
        for index, record in enumerate(records, start=1)
    }
    surrogate_records = [
        FastaRecord(surrogate, records[index].sequence)
        for index, surrogate in enumerate(surrogate_to_header)
    ]

    if keep_temp:
        temp_dir = Path(tempfile.mkdtemp(prefix="hcv_cluster_prep_mafft_"))
        logger.info("Keeping temporary MAFFT files in %s", temp_dir)
        return run_mafft_in_tempdir(
            temp_dir,
            surrogate_records,
            surrogate_to_header,
            reference_fasta,
            reference_accession,
            mafft_exe,
            threads,
            logger,
        )

    with tempfile.TemporaryDirectory(prefix="hcv_cluster_prep_mafft_") as temp_dir_name:
        return run_mafft_in_tempdir(
            Path(temp_dir_name),
            surrogate_records,
            surrogate_to_header,
            reference_fasta,
            reference_accession,
            mafft_exe,
            threads,
            logger,
        )


def run_mafft_in_tempdir(
    temp_dir: Path,
    surrogate_records: list[FastaRecord],
    surrogate_to_header: dict[str, str],
    reference_fasta: Path,
    reference_accession: str,
    mafft_exe: str,
    threads: int,
    logger: logging.Logger,
) -> tuple[str, dict[str, str]]:
    fragments_path = temp_dir / "fragments.fasta"
    mafft_stdout_path = temp_dir / "mafft.aligned.raw.fasta"
    write_fasta_records(surrogate_records, fragments_path)

    base_cmd = [
        mafft_exe,
        "--quiet",
        "--thread",
        str(threads),
        "--keeplength",
        "--addfragments",
        str(fragments_path),
        str(reference_fasta),
    ]
    cmd = [base_cmd[0], "--anysymbol", *base_cmd[1:]]
    logger.info("Running MAFFT --addfragments for %s input sequences", len(surrogate_records))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if proc.returncode != 0 and "anysymbol" in proc.stderr.lower():
        logger.warning("MAFFT rejected --anysymbol; retrying without it")
        proc = subprocess.run(base_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if proc.returncode != 0:
        raise HcvPrepError(f"MAFFT failed with exit code {proc.returncode}: {proc.stderr.strip()}")
    if not proc.stdout.strip():
        raise HcvPrepError("MAFFT completed but produced no alignment output")

    mafft_stdout_path.write_text(proc.stdout)
    aligned_records = parse_fasta_text(proc.stdout, "MAFFT stdout")
    aligned_by_header = {record.header: record.sequence for record in aligned_records}

    if reference_accession in aligned_by_header:
        reference_aligned = aligned_by_header[reference_accession]
    else:
        reference_aligned = aligned_records[0].sequence
        logger.warning(
            "Reference accession '%s' was not found in MAFFT output; using first aligned sequence as reference",
            reference_accession,
        )

    query_alignment: dict[str, str] = {}
    missing: list[str] = []
    for surrogate, original_header in surrogate_to_header.items():
        aligned = aligned_by_header.get(surrogate)
        if aligned is None:
            missing.append(original_header)
        else:
            query_alignment[original_header] = aligned

    if missing:
        raise HcvPrepError(f"MAFFT output is missing aligned sequences for: {', '.join(missing[:10])}")

    return reference_aligned, query_alignment


def reuse_cached_alignment(
    records: list[FastaRecord],
    cache_path: Path,
    reference_fasta: Path,
    reference_accession: str,
    logger: logging.Logger,
) -> tuple[str, dict[str, str], list[FastaRecord]]:
    """Reuse unchanged records from a previous reference-anchored alignment."""
    if not cache_path.exists():
        raise HcvPrepError(f"Cached alignment does not exist: {cache_path}")
    cached_records = read_fasta_records(cache_path)
    reject_duplicate_ids(cached_records)
    cached_by_header = {record.header: record.sequence for record in cached_records}

    reference_aligned = cached_by_header.pop(reference_accession, None)
    if reference_aligned is None:
        raise HcvPrepError(
            f"Cached alignment {cache_path} does not contain reference {reference_accession}"
        )

    reference_records = read_fasta_records(reference_fasta)
    expected_reference = normalize_reference_sequence(reference_records[0].sequence)
    if reference_aligned.replace("-", "").upper() != expected_reference:
        raise HcvPrepError(
            f"Cached alignment {cache_path} was not built with reference {reference_accession}"
        )

    alignment_length = len(reference_aligned)
    if any(len(sequence) != alignment_length for sequence in cached_by_header.values()):
        raise HcvPrepError(f"Cached alignment {cache_path} contains inconsistent sequence lengths")

    reused: dict[str, str] = {}
    to_align: list[FastaRecord] = []
    for record in records:
        cached = cached_by_header.get(record.header)
        if cached is not None and cached.replace("-", "").upper() == record.sequence:
            reused[record.header] = cached
        else:
            to_align.append(record)

    logger.info(
        "Reusing %s/%s unchanged sequences from %s; %s require alignment",
        len(reused),
        len(records),
        cache_path,
        len(to_align),
    )
    return reference_aligned, reused, to_align


def reference_position_columns(reference_aligned: str) -> dict[int, int]:
    reference_position = 0
    columns: dict[int, int] = {}
    for column, char in enumerate(reference_aligned):
        if char == "-":
            continue
        reference_position += 1
        columns[reference_position] = column
    return columns


def reference_region_columns(reference_aligned: str, segment: RegionSegment | CoreE2Coords) -> tuple[int, int]:
    columns = reference_position_columns(reference_aligned)
    start_column = columns.get(segment.start)
    end_reference_column = columns.get(segment.end)
    if start_column is None or end_reference_column is None:
        raise HcvPrepError(
            f"Could not map reference coordinates {segment.start}-{segment.end} "
            "onto the MAFFT alignment"
        )
    return start_column, end_reference_column + 1


def extract_region_selection(
    reference_aligned: str,
    query_alignment: dict[str, str],
    selection: RegionSelection,
) -> dict[str, str]:
    segment_columns = [
        reference_region_columns(reference_aligned, segment)
        for segment in selection.segments
    ]
    extracted: dict[str, str] = {}
    for sample_id, aligned_sequence in query_alignment.items():
        if len(aligned_sequence) != len(reference_aligned):
            raise HcvPrepError(f"Aligned sequence length mismatch for '{sample_id}'")
        extracted[sample_id] = "".join(
            aligned_sequence[start_column:end_column]
            for start_column, end_column in segment_columns
        )
    return extracted


def extract_core_e2(
    reference_aligned: str,
    query_alignment: dict[str, str],
    coords: CoreE2Coords,
) -> dict[str, str]:
    segment = RegionSegment("core-e2", coords.start, coords.end, coords.source, coords.detail)
    selection = RegionSelection("core-e2", "fixed", (segment,))
    return extract_region_selection(reference_aligned, query_alignment, selection)


def compute_coverage(trimmed_sequence: str, expected_length: int) -> CoverageMetrics:
    gap_bases = sum(1 for char in trimmed_sequence if char == "-")
    n_bases = sum(1 for char in trimmed_sequence if char.upper() == "N")
    non_gap_bases = sum(1 for char in trimmed_sequence if char != "-" and char.upper() != "N")
    coverage = non_gap_bases / expected_length if expected_length else 0.0
    return CoverageMetrics(
        extracted_length=len(trimmed_sequence),
        non_gap_bases=non_gap_bases,
        n_bases=n_bases,
        gap_bases=gap_bases,
        coverage_fraction=coverage,
        n_fraction=n_bases / expected_length if expected_length else 0.0,
        gap_fraction=gap_bases / expected_length if expected_length else 0.0,
    )


def mask_n_as_gap(sequence: str) -> str:
    """Replace N (case-insensitive) with '-' for the clustering FASTA only. tn93
    excludes alignment gaps from every -a mode uniformly, but N is handled very
    differently per mode (e.g. `average` can manufacture large distance purely
    from a low-depth-masked block with zero real information — see
    docs/threshold_rationale.md). Masking N as a gap here neutralizes it the
    same way regardless of --ambiguities, while leaving real IUPAC ambiguity
    codes (R, Y, etc.) untouched so their proportional treatment is preserved."""
    return sequence.replace("N", "-").replace("n", "-")


def select_max_usable_region(
    records: list[FastaRecord],
    reference_aligned: str,
    query_alignment: dict[str, str],
    search_selection: RegionSelection,
    min_coverage: float,
    window_size: int,
    window_step: int,
    min_window_length: int,
    min_taxa_fraction: float,
    window_coverage: float,
    logger: logging.Logger,
) -> RegionSelection:
    if window_size < 1:
        raise HcvPrepError("--max-usable-window-size must be at least 1")
    if window_step < 1:
        raise HcvPrepError("--max-usable-window-step must be at least 1")
    if min_window_length < 1:
        raise HcvPrepError("--max-usable-min-window-length must be at least 1")
    if not 0 <= min_taxa_fraction <= 1:
        raise HcvPrepError("--max-usable-min-taxa-fraction must be between 0 and 1")
    if not 0 <= window_coverage <= 1:
        raise HcvPrepError("--max-usable-window-coverage must be between 0 and 1")

    sample_ids = [record.header for record in records]
    min_taxa = max(1, math.ceil(len(sample_ids) * min_taxa_fraction))
    candidates = score_window_candidates(
        records=records,
        reference_aligned=reference_aligned,
        query_alignment=query_alignment,
        search_selection=search_selection,
        window_size=window_size,
        window_step=window_step,
        min_window_length=min_window_length,
        window_coverage=window_coverage,
    )
    if not candidates:
        raise HcvPrepError(
            "The max-usable selector found no candidate windows. "
            "Try a smaller --max-usable-min-window-length or a broader --region search space."
        )

    best_selection: RegionSelection | None = None
    best_score: tuple[int, int, int, int, int] | None = None
    best_retained = 0
    best_threshold = 0
    for taxa_threshold in range(min_taxa, len(sample_ids) + 1):
        window_segments = [
            candidate.segment
            for candidate in candidates
            if candidate.support_count >= taxa_threshold
        ]
        if not window_segments:
            continue

        merged = merge_region_segments(window_segments)
        candidate_selection = RegionSelection(
            expression=search_selection.expression,
            strategy="max-usable",
            segments=tuple(merged),
            detail=f"window_taxa_threshold={taxa_threshold}",
        )
        extracted = extract_region_selection(reference_aligned, query_alignment, candidate_selection)
        retained = sum(
            1
            for sample_id in sample_ids
            if compute_coverage(extracted[sample_id], candidate_selection.length).coverage_fraction >= min_coverage
        )
        if retained < min_taxa:
            continue

        score = (
            retained * candidate_selection.length,
            retained,
            candidate_selection.length,
            taxa_threshold,
            -len(candidate_selection.segments),
        )
        if best_score is None or score > best_score:
            best_score = score
            best_selection = candidate_selection
            best_retained = retained
            best_threshold = taxa_threshold

    if best_selection is None:
        raise HcvPrepError(
            "The max-usable selector could not find reproducible windows that retain "
            f"at least {min_taxa}/{len(sample_ids)} taxa at the requested coverage settings."
        )

    relabeled = tuple(
        RegionSegment(
            name=f"max{index:02d}",
            start=segment.start,
            end=segment.end,
            source="input_alignment_max_usable",
            detail=(
                f"search={search_selection.expression}; window_size={window_size}; "
                f"window_step={window_step}; min_window_length={min_window_length}; "
                f"window_coverage={window_coverage}; taxa_threshold={best_threshold}"
            ),
        )
        for index, segment in enumerate(best_selection.segments, start=1)
    )
    selected = RegionSelection(
        expression=search_selection.expression,
        strategy="max-usable",
        segments=relabeled,
        detail=best_selection.detail,
    )
    logger.info(
        "Max-usable selector chose %s nt across %s segment(s), retaining %s/%s taxa at %.3f coverage",
        selected.length,
        len(selected.segments),
        best_retained,
        len(sample_ids),
        min_coverage,
    )
    logger.info("Max-usable selected segments: %s", selected.segment_text)
    return selected


def score_window_candidates(
    records: list[FastaRecord],
    reference_aligned: str,
    query_alignment: dict[str, str],
    search_selection: RegionSelection,
    window_size: int,
    window_step: int,
    min_window_length: int,
    window_coverage: float,
) -> list[WindowCandidate]:
    columns = reference_position_columns(reference_aligned)
    candidates: list[WindowCandidate] = []
    for segment in generate_candidate_windows(search_selection, window_size, window_step, min_window_length):
        sample_coverages = [
            segment_coverage_fraction(record.header, query_alignment, columns, segment)
            for record in records
        ]
        support_count = sum(1 for coverage in sample_coverages if coverage >= window_coverage)
        mean_coverage = sum(sample_coverages) / len(sample_coverages) if sample_coverages else 0.0
        candidates.append(WindowCandidate(segment, support_count, mean_coverage))
    return candidates


def generate_candidate_windows(
    search_selection: RegionSelection,
    window_size: int,
    window_step: int,
    min_window_length: int,
) -> Iterable[RegionSegment]:
    for search_segment in search_selection.segments:
        window_start = search_segment.start
        window_index = 1
        while window_start <= search_segment.end:
            window_end = min(window_start + window_size - 1, search_segment.end)
            if window_end - window_start + 1 >= min_window_length:
                yield RegionSegment(
                    name=f"{search_segment.name}_win{window_index:03d}",
                    start=window_start,
                    end=window_end,
                    source="input_alignment_candidate_window",
                    detail=f"search={search_segment.name}:{search_segment.start}-{search_segment.end}",
                )
            if window_end == search_segment.end:
                break
            window_start += window_step
            window_index += 1


def segment_coverage_fraction(
    sample_id: str,
    query_alignment: dict[str, str],
    columns: dict[int, int],
    segment: RegionSegment,
) -> float:
    aligned = query_alignment[sample_id]
    usable = 0
    for reference_position in range(segment.start, segment.end + 1):
        column = columns.get(reference_position)
        if column is None:
            raise HcvPrepError(
                f"Could not map reference coordinate {reference_position} onto the MAFFT alignment"
            )
        base = aligned[column].upper()
        if base != "-" and base != "N":
            usable += 1
    return usable / segment.length if segment.length else 0.0


def write_outputs(
    records: list[FastaRecord],
    query_alignment: dict[str, str],
    reference_aligned: str,
    extracted: dict[str, str],
    selection: RegionSelection,
    spec: ReferenceSpec,
    out_prefix: Path,
    min_coverage: float,
    logger: logging.Logger,
) -> None:
    clustering_path = Path(f"{out_prefix}.clustering.fasta")
    qc_path = Path(f"{out_prefix}.qc.csv")
    aligned_path = Path(f"{out_prefix}.aligned.fasta")

    retained: list[FastaRecord] = []
    qc_rows: list[dict[str, str | int | float | bool]] = []
    removal_reason = coverage_reason_slug(selection)
    is_core_e2_family = selection.expression in {"core-e2", "core-e2-nohvr1"}
    core_e2_start = selection.start if is_core_e2_family else ""
    core_e2_end = selection.end if is_core_e2_family else ""
    core_e2_length = selection.length if is_core_e2_family else ""
    for record in records:
        sample_id = record.header
        trimmed = extracted[sample_id]
        metrics = compute_coverage(trimmed, selection.length)
        passed = metrics.coverage_fraction >= min_coverage
        reason = "" if passed else removal_reason
        if passed:
            retained.append(FastaRecord(sample_id, mask_n_as_gap(trimmed)))
        qc_rows.append(
            {
                "sample_id": sample_id,
                "genotype": spec.genotype,
                "reference_accession": spec.accession,
                "selected_region": selection.expression,
                "region_strategy": selection.strategy,
                "region_segments": selection.segment_text,
                "region_start": selection.start,
                "region_end": selection.end,
                "region_length": selection.length,
                "core_e2_start": core_e2_start,
                "core_e2_end": core_e2_end,
                "core_e2_length": core_e2_length,
                "extracted_length": metrics.extracted_length,
                "non_gap_bases": metrics.non_gap_bases,
                "n_bases": metrics.n_bases,
                "gap_bases": metrics.gap_bases,
                "coverage_fraction": f"{metrics.coverage_fraction:.6f}",
                "n_fraction": f"{metrics.n_fraction:.6f}",
                "gap_fraction": f"{metrics.gap_fraction:.6f}",
                "passed_qc": "true" if passed else "false",
                "removal_reason": reason,
            }
        )

    write_fasta_records(retained, clustering_path)
    write_qc_csv(qc_rows, qc_path)

    aligned_records = [FastaRecord(spec.accession, reference_aligned)]
    aligned_records.extend(FastaRecord(record.header, query_alignment[record.header]) for record in records)
    write_fasta_records(aligned_records, aligned_path)

    logger.info("Wrote clustering FASTA: %s (%s/%s retained)", clustering_path, len(retained), len(records))
    logger.info("Wrote QC CSV: %s", qc_path)
    logger.info("Wrote full reference-anchored alignment: %s", aligned_path)
    if not retained:
        logger.warning("No sequences passed the %.3f selected-region coverage threshold", min_coverage)


def write_qc_csv(rows: list[dict[str, str | int | float | bool]], path: Path) -> None:
    fieldnames = [
        "sample_id",
        "genotype",
        "reference_accession",
        "selected_region",
        "region_strategy",
        "region_segments",
        "region_start",
        "region_end",
        "region_length",
        "core_e2_start",
        "core_e2_end",
        "core_e2_length",
        "extracted_length",
        "non_gap_bases",
        "n_bases",
        "gap_bases",
        "coverage_fraction",
        "n_fraction",
        "gap_fraction",
        "passed_qc",
        "removal_reason",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def prepare_one_reference(
    spec: ReferenceSpec,
    refs_dir: Path,
    catalog: dict[str, ReferenceSpec],
    force_download: bool,
    email: str | None,
    logger: logging.Logger,
) -> tuple[ReferenceBundle, dict[str, RegionSegment]]:
    bundle = download_reference(spec, refs_dir, force_download, email, logger)
    regions = parse_reference_region_map(
        genbank_path=bundle.genbank_path,
        refs_dir=refs_dir,
        catalog=catalog,
        force_download=force_download,
        email=email,
        logger=logger,
    )
    write_reference_region_cache(bundle, regions, logger)
    return bundle, regions


def command_prep_refs(args: argparse.Namespace) -> int:
    logger = setup_logging(Path(args.log) if args.log else None, args.verbose)
    catalog = load_reference_catalog(Path(args.reference_map) if args.reference_map else None)
    refs_dir = Path(args.refs_dir)

    if normalize_genotype(args.genotype) == "all":
        specs = [catalog[key] for key in sorted(catalog)]
    else:
        specs = [get_reference_spec(args.genotype, catalog)]

    for spec in specs:
        prepare_one_reference(
            spec=spec,
            refs_dir=refs_dir,
            catalog=catalog,
            force_download=args.force_download,
            email=args.ncbi_email,
            logger=logger,
        )
    return 0


def command_prep_align(args: argparse.Namespace) -> int:
    out_prefix = Path(args.out_prefix)
    logger = setup_logging(Path(f"{out_prefix}.log"), args.verbose)
    if not 0 <= args.min_coverage <= 1:
        raise HcvPrepError("--min-coverage must be between 0 and 1")

    catalog = load_reference_catalog(Path(args.reference_map) if args.reference_map else None)
    spec = get_reference_spec(args.genotype, catalog)
    refs_dir = Path(args.refs_dir)

    bundle, regions = prepare_one_reference(
        spec=spec,
        refs_dir=refs_dir,
        catalog=catalog,
        force_download=args.force_download,
        email=args.ncbi_email,
        logger=logger,
    )

    input_path = Path(args.input)
    records = read_fasta_records(input_path)
    reject_duplicate_ids(records)
    normalized_records = normalize_input_records(records, logger)
    logger.info("Read %s unique input sequences from %s", len(normalized_records), input_path)

    validate_input_genotypes(
        records=normalized_records,
        requested_genotype=spec.genotype,
        mode=args.genotype_validation,
        kmer_size=args.genotype_validation_kmer_size,
        refs_dir=refs_dir,
        catalog=catalog,
        force_download=args.force_download,
        email=args.ncbi_email,
        logger=logger,
    )

    cached_reference: str | None = None
    query_alignment: dict[str, str] = {}
    records_to_align = normalized_records
    if args.cached_alignment:
        cache_path = Path(args.cached_alignment)
        cached_reference, query_alignment, records_to_align = reuse_cached_alignment(
            records=normalized_records,
            cache_path=cache_path,
            reference_fasta=bundle.fasta_path,
            reference_accession=spec.accession,
            logger=logger,
        )

    if records_to_align:
        reference_aligned, new_alignment = run_mafft_addfragments(
            records=records_to_align,
            reference_fasta=bundle.fasta_path,
            reference_accession=spec.accession,
            mafft_command=args.mafft,
            threads=args.threads,
            keep_temp=args.keep_temp,
            logger=logger,
        )
        if cached_reference is not None and reference_aligned != cached_reference:
            raise HcvPrepError("New and cached alignments use different reference coordinates")
        query_alignment.update(new_alignment)
    elif cached_reference is not None:
        reference_aligned = cached_reference
    else:  # pragma: no cover - read_fasta_records rejects an empty input earlier
        raise HcvPrepError("No sequences were available to align")

    query_alignment = {
        record.header: query_alignment[record.header]
        for record in normalized_records
    }
    region_expression = args.region
    if region_expression is None:
        region_expression = "cds" if args.region_strategy == "max-usable" else "core-e2-nohvr1"

    search_selection = resolve_region_selection(
        region_expression,
        regions,
        strategy=args.region_strategy,
    )
    if args.region_strategy == "max-usable":
        selection = select_max_usable_region(
            records=normalized_records,
            reference_aligned=reference_aligned,
            query_alignment=query_alignment,
            search_selection=search_selection,
            min_coverage=args.min_coverage,
            window_size=args.max_usable_window_size,
            window_step=args.max_usable_window_step,
            min_window_length=args.max_usable_min_window_length,
            min_taxa_fraction=args.max_usable_min_taxa_fraction,
            window_coverage=(
                args.max_usable_window_coverage
                if args.max_usable_window_coverage is not None
                else args.min_coverage
            ),
            logger=logger,
        )
    else:
        selection = search_selection
        logger.info(
            "Selected fixed region %s: %s (%s nt)",
            selection.expression,
            selection.segment_text,
            selection.length,
        )

    extracted = extract_region_selection(reference_aligned, query_alignment, selection)
    write_outputs(
        records=normalized_records,
        query_alignment=query_alignment,
        reference_aligned=reference_aligned,
        extracted=extracted,
        selection=selection,
        spec=spec,
        out_prefix=out_prefix,
        min_coverage=args.min_coverage,
        logger=logger,
    )
    logger.info(
        "Reference is used only as a genotype-specific coordinate anchor; "
        "this preprocessing does not infer transmission."
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--refs-dir", default="refs", help="Directory for cached reference GenBank/FASTA files")
    common.add_argument("--reference-map", help="JSON file extending or overriding genotype-to-accession mappings")
    common.add_argument("--ncbi-email", default=os.environ.get("NCBI_EMAIL"), help="Email passed to NCBI E-utilities")
    common.add_argument("--force-download", action="store_true", help="Re-download cached NCBI reference records")
    common.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging")

    parser = argparse.ArgumentParser(
        prog="hcv_cluster_prep",
        description="Genotype-specific HCV region preprocessing for MicrobeTrace clustering.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    refs_parser = subparsers.add_parser(
        "prep-refs",
        parents=[common],
        help="Download/cache reference records and determine mature-region coordinates",
    )
    refs_parser.add_argument("--genotype", default="all", help="Genotype to prepare, or 'all'")
    refs_parser.add_argument("--log", help="Optional log file for reference preparation")
    refs_parser.set_defaults(func=command_prep_refs)

    align_parser = subparsers.add_parser(
        "prep-align",
        parents=[common],
        help="Align input sequences to a genotype reference and export region-trimmed QC outputs",
    )
    align_parser.add_argument("--input", required=True, help="Input multi-FASTA of unique HCV consensus IDs")
    align_parser.add_argument("--genotype", required=True, help="Run genotype; built-ins: 1a, 1b, 2a, 2b, 3a")
    align_parser.add_argument("--out-prefix", required=True, help="Prefix for .clustering.fasta, .qc.csv, .aligned.fasta, .log")
    align_parser.add_argument(
        "--region",
        help=(
            "Reference-anchored region expression. Defaults to core-e2-nohvr1 (Core-through-E2, "
            "HVR1 masked) for fixed mode and cds for max-usable mode. Examples: core, e1-e2, "
            "core-e2-nohvr1, core-e2 (HVR1 included), ns3, ns5a-ns5b, core-e2-nohvr1+ns3, cds."
        ),
    )
    align_parser.add_argument(
        "--region-strategy",
        choices=("fixed", "max-usable"),
        default="fixed",
        help="Use a named fixed region or select reproducible high-coverage windows from --region",
    )
    align_parser.add_argument(
        "--max-usable-region",
        action="store_const",
        const="max-usable",
        dest="region_strategy",
        help="Alias for --region-strategy max-usable",
    )
    align_parser.add_argument("--min-coverage", type=float, default=0.8, help="Minimum non-gap non-N selected-region coverage")
    align_parser.add_argument(
        "--max-usable-window-size",
        type=int,
        default=1000,
        help="Window size, in reference bases, for max-usable region selection",
    )
    align_parser.add_argument(
        "--max-usable-window-step",
        type=int,
        default=1000,
        help="Window step, in reference bases, for max-usable region selection",
    )
    align_parser.add_argument(
        "--max-usable-min-window-length",
        type=int,
        default=200,
        help="Minimum final partial window length for max-usable region selection",
    )
    align_parser.add_argument(
        "--max-usable-min-taxa-fraction",
        type=float,
        default=0.8,
        help="Minimum fraction of taxa the max-usable selector must retain",
    )
    align_parser.add_argument(
        "--max-usable-window-coverage",
        type=float,
        default=None,
        help="Per-window coverage required for a taxon to support a max-usable candidate window; defaults to --min-coverage",
    )
    align_parser.add_argument("--mafft", default="mafft", help="MAFFT executable path/name")
    align_parser.add_argument("--threads", type=int, default=1, help="MAFFT thread count")
    align_parser.add_argument(
        "--cached-alignment",
        help=(
            "Previous .aligned.fasta to reuse. Unchanged sample IDs are taken from the "
            "cache; only new or changed sequences are sent to MAFFT."
        ),
    )
    align_parser.add_argument(
        "--genotype-validation",
        choices=("auto", "headers", "kmer", "none"),
        default="auto",
        help="Validate requested genotype by FASTA header tags, k-mer reference comparison, or skip",
    )
    align_parser.add_argument(
        "--genotype-validation-kmer-size",
        type=int,
        default=15,
        help="K-mer size for reference-anchor genotype validation",
    )
    align_parser.add_argument("--keep-temp", action="store_true", help="Keep temporary MAFFT input/output files")
    align_parser.set_defaults(func=command_prep_align)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except HcvPrepError as exc:
        logger = logging.getLogger(LOGGER_NAME)
        if logger.handlers:
            logger.error("%s", exc)
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
