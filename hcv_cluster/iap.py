#!/usr/bin/env python3
"""The IAP model: pairwise SNP counts, collection-date differences and a relatedness index.

For every pair of sequences within one subtype that both have a collection date, write
the raw number of differing nucleotides, the time between the two collections in years,
their ratio and the relatedness index, with the subtype and a note where a value is empty. Each subtype gets its own table
and ``merge_pairs`` concatenates them. Pairs are never made across subtypes: each subtype
is aligned separately, so a raw difference between them is not defined (and they are far
too divergent for it to mean anything).

The relatedness index is the probability of seeing more differences than were observed
between two sequences from one transmission chain collected that far apart, under a
negative-binomial model (``relatedness_index``). The model is calibrated for the full
E1-E2 region (HVR1 included), so ``check_region`` refuses any other region.

``raw_dist`` counts the sites where the two sequences differ. A site where either has a gap or N
is not compared, so indels add nothing. A site where either carries an IUPAC ambiguity code is,
by default, counted as the chance that the two differ if each takes one of its possible bases
at random (R against A is 0.5): the best stand-in for the difference between majority-base
calls, which is what the model was fitted on. ``--iap-ambiguity skip`` ignores those sites instead.
Pairs are flagged where the count may mislead: ``low_coverage`` (the pair compares under 95% of
the sites of the most complete sequence, so the count is low) and ``noisy`` (a sequence has many
ambiguity codes, as a mixed or poor sample does).
"""

from __future__ import annotations

import calendar
import csv
import datetime
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.stats import nbinom

from . import metadata, methods, prep

IAP_FILENAME = "iap_pairs.csv"
DEFAULT_DATE_COLUMN = "sample_collection_date"
DEFAULT_ALPHA = 0.05  # a link is IAP-supported when its relatedness index is at least this
DATE_FORMATS = "d/m/yyyy, yyyy-mm-dd or d-Mon-yy (e.g. 5-Mar-24)"
CALIBRATED_REGION = "e1-e2"
IAP_FIELDS = ("seqid_1", "seqid_2", "raw_dist", "time_diff_years", "snp_rate", "subtype", "relatedness_index", "note",
              "compared_sites", "ambiguous_sites", "flags")
AMBIGUITY_MODES = ("average", "skip")
DEFAULT_AMBIGUITY = "average"
LOW_OVERLAP = 0.95  # a pair compares less than this share of the most complete sequence's sites: low_coverage
NOISY_AMBIGUOUS_SITES = 10  # a sequence with at least this many ambiguity codes in the region: noisy
FLAG_LOW_COVERAGE = "low_coverage"
FLAG_NOISY = "noisy"
NOTE_TOO_FAR_APART = "too_far_apart"  # beyond the time range of the model: no relatedness_index
NOTE_SAME_DAY = "same_day"  # no time between the samples: no snp_rate

# Fitted parameters of the model: the expected number of differences between two related
# sequences is A * years + B, with a constant variance of C ** 2.
MODEL_A = 4.417171969933791
MODEL_B = 9.634037423658565
MODEL_C = 10.382634613596712

_DAY_FIRST = re.compile(r"(\d{1,2})/(\d{1,2})/(\d{4})")
_ISO = re.compile(r"(\d{4})-(\d{2})-(\d{2})")
_DAY_MONTH_NAME = re.compile(r"(\d{1,2})[- ]([A-Za-z]{3,9})[- ](\d{2}|\d{4})")
_MONTH_NAMES = (
    "january", "february", "march", "april", "may", "june",
    "july", "august", "september", "october", "november", "december",
)


# IUPAC codes and the bases they stand for. N, gaps and anything else carry no information.
_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "K": "GT", "M": "AC", "S": "CG",
          "W": "AT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}


def _base_probabilities(sequences: Sequence[str]) -> np.ndarray:
    """(sequences, sites, 4): the chance of each of A, C, G, T at each site; all zero where there is no information."""
    table = np.zeros((256, 4), dtype=np.float32)
    for code, bases in _IUPAC.items():
        for base in bases:
            table[ord(code), "ACGT".index(base)] = 1 / len(bases)
    width = max(len(sequence) for sequence in sequences)
    codes = np.full((len(sequences), width), ord("-"), dtype=np.uint8)
    for row, sequence in enumerate(sequences):
        codes[row, : len(sequence)] = np.frombuffer(sequence.upper().encode("ascii", "replace"), dtype=np.uint8)
    return table[codes]


@dataclass
class PairCount:
    id_1: str
    id_2: str
    distance: float  # differing sites; fractional where ambiguity codes are counted by chance
    compared: int  # sites both sequences have a base (or a code) at
    ambiguous: int  # of those, sites where either carries an ambiguity code
    flags: str = ""


def count_pairs(fasta: Path, ambiguity: str = DEFAULT_AMBIGUITY) -> list[PairCount]:
    """Differences between every pair of sequences in a clustering FASTA, with the flags described above."""
    if ambiguity not in AMBIGUITY_MODES:
        raise ValueError(f"--iap-ambiguity must be one of {', '.join(AMBIGUITY_MODES)}, got {ambiguity!r}")
    records = prep.read_fasta_records(fasta)
    if len(records) < 2:
        return []
    ids = [record.header.split()[0] for record in records]
    prob = _base_probabilities([record.sequence for record in records])
    called = prob.sum(axis=2) > 0
    ambiguous = called & (prob.max(axis=2) < 0.999)
    most_complete = max(int(called.sum(axis=1).max()), 1)
    codes_per_sequence = ambiguous.sum(axis=1)
    pairs: list[PairCount] = []
    for i in range(len(records) - 1):
        j = slice(i + 1, None)
        both = called[i] & called[j]
        agree = np.einsum("lk,mlk->ml", prob[i], prob[j])  # the chance the two sites hold the same base
        either = ambiguous[i] | ambiguous[j]
        if ambiguity == "average":
            comparable = both
            distance = ((1 - agree) * comparable).sum(axis=1)
        else:
            comparable = both & ~either
            distance = ((agree < 0.5) & comparable).sum(axis=1)
        for offset, (dist, comp) in enumerate(zip(distance, comparable.sum(axis=1))):
            k = i + 1 + offset
            flags = []
            if both[offset].sum() / most_complete < LOW_OVERLAP:
                flags.append(FLAG_LOW_COVERAGE)
            if max(codes_per_sequence[i], codes_per_sequence[k]) >= NOISY_AMBIGUOUS_SITES:
                flags.append(FLAG_NOISY)
            pairs.append(PairCount(ids[i], ids[k], float(dist), int(comp), int((both[offset] & either[offset]).sum()), ";".join(flags)))
    return pairs


def _count_text(value: float) -> str:
    return str(int(round(value))) if abs(value - round(value)) < 1e-9 else f"{value:.3f}"


@dataclass
class PairSummary:
    pairs: int = 0
    skipped_pairs: int = 0  # pairs involving a sequence with no collection date
    undated_sequences: int = 0
    unscored_pairs: int = 0  # too far apart in time for the model to give an index
    low_coverage_pairs: int = 0
    noisy_pairs: int = 0


def check_region(region: str, region_strategy: str = "fixed") -> None:
    """Refuse to produce IAP input for any region but the one the model was calibrated on."""
    if methods.canonical_threshold_region(region) != CALIBRATED_REGION or region_strategy != "fixed":
        chosen = region if region_strategy == "fixed" else f"{region} ({region_strategy} strategy)"
        raise ValueError(
            f"The IAP model is calibrated for the full {CALIBRATED_REGION} region (HVR1 included, "
            f"fixed region strategy), not '{chosen}'. Use --region {CALIBRATED_REGION}, "
            "or run without --iap-model"
        )


def relatedness_index(raw_dist: Sequence[float], time_years: Sequence[float]) -> np.ndarray:
    """Chance that two sequences from one transmission chain differ at more sites than observed.

    The differences follow a negative binomial with mean ``MODEL_A * years + MODEL_B`` and
    variance ``MODEL_C ** 2``. That needs the mean to stay below the variance, i.e. a gap
    of under about 22 years; beyond it the index is NaN.
    """
    distance = np.asarray(raw_dist, dtype=float)
    mean = MODEL_A * np.asarray(time_years, dtype=float) + MODEL_B
    prob = mean / MODEL_C**2
    with np.errstate(divide="ignore", invalid="ignore"):
        size = (mean**2 / MODEL_C**2) / (1 - prob)
    valid = np.isfinite(size) & (size > 0) & np.isfinite(prob) & (prob > 0) & (prob < 1)
    index = np.full(distance.shape, np.nan)
    index[valid] = nbinom.sf(distance[valid], size[valid], prob[valid])
    return index


def parse_date(value: str) -> datetime.date:
    """Read ``d/m/yyyy`` (day first), ``yyyy-mm-dd`` or ``d-Mon-yy`` / ``d-Mon-yyyy``.

    A month name may be abbreviated (Mar, Sept) or in full, in any case. A two-digit year
    follows the usual convention: 69-99 are 1969-1999 and 00-68 are 2000-2068.
    """
    text = value.strip()
    if match := _DAY_FIRST.fullmatch(text):
        day, month, year = map(int, match.groups())
    elif match := _ISO.fullmatch(text):
        year, month, day = map(int, match.groups())
    elif (match := _DAY_MONTH_NAME.fullmatch(text)) and (
        months := [i for i, name in enumerate(_MONTH_NAMES, start=1) if name.startswith(match[2].lower())]
    ):
        day, month, year = int(match[1]), months[0], int(match[3])
        if len(match[3]) == 2:
            year += 1900 if year >= 69 else 2000
    else:
        raise ValueError(f"{value!r} is not {DATE_FORMATS}")
    return datetime.date(year, month, day)


def decimal_months(date: datetime.date) -> float:
    """Months since year 0, the day of the month counting as a fraction of its month."""
    days_in_month = calendar.monthrange(date.year, date.month)[1]
    return date.year * 12 + (date.month - 1) + (date.day - 1) / days_in_month


def years_between(first: datetime.date, second: datetime.date) -> float:
    """Time between two dates in years, e.g. 1/7/1995 to 1/12/1997 is 29 months = 2.416667."""
    return abs(decimal_months(second) - decimal_months(first)) / 12


def collection_dates(
    metadata_rows: Sequence[Mapping[str, Any]], column: str = DEFAULT_DATE_COLUMN
) -> dict[str, datetime.date]:
    """Map sample_id to its collection date; samples without one are left out."""
    if not metadata_rows:
        raise ValueError("--iap-model needs --metadata with a collection date column")
    available = metadata.metadata_fields(metadata_rows)
    if column not in available:
        raise ValueError(
            f"'{column}' is not a metadata column (available: {', '.join(available)}); "
            "name the collection date column with --date-column"
        )
    dates: dict[str, datetime.date] = {}
    invalid: list[str] = []
    for row in metadata_rows:
        value = str(row.get(column) or "").strip()
        if not value or value == metadata.MISSING_VALUE:
            continue
        try:
            dates[str(row["sample_id"])] = parse_date(value)
        except ValueError:
            invalid.append(f"{row['sample_id']}: {value!r}")
    if invalid:
        raise ValueError(
            f"Invalid '{column}' value(s), expected {DATE_FORMATS}:\n- " + "\n- ".join(invalid)
        )
    return dates


def _fasta_has_records(path: Path) -> bool:
    with path.open() as handle:
        return any(line.startswith(">") for line in handle)


def write_pairs(
    fasta: Path,
    subtype: str,
    dates: Mapping[str, datetime.date],
    out: Path,
    base_ids: Mapping[str, str] | None = None,
    ambiguity: str = DEFAULT_AMBIGUITY,
) -> PairSummary:
    """Write ``IAP_FIELDS`` for every dated pair in one subtype's clustering FASTA.

    Dates are keyed by metadata sample_id, so renamed mixed-infection variants
    (``SAMPLE__v2``) are looked up through ``base_ids``. Pairs collected on the same day
    are kept, with an empty ``snp_rate`` (undefined for them) but a relatedness index.
    ``ambiguity`` is how IUPAC codes are counted (see the module notes); the model needs a whole
    number, so a fractional ``raw_dist`` is rounded to the nearest for the index only.
    """
    base_ids = base_ids or {}
    summary = PairSummary()
    undated: set[str] = set()
    pairs: list[tuple[PairCount, float]] = []
    if _fasta_has_records(fasta):
        for count in count_pairs(fasta, ambiguity):
            date_1 = dates.get(base_ids.get(count.id_1, count.id_1))
            date_2 = dates.get(base_ids.get(count.id_2, count.id_2))
            if date_1 is None or date_2 is None:
                summary.skipped_pairs += 1
                undated.update(i for i, d in ((count.id_1, date_1), (count.id_2, date_2)) if d is None)
                continue
            pairs.append((count, years_between(date_1, date_2)))
    summary.undated_sequences = len(undated)
    summary.pairs = len(pairs)

    indices = relatedness_index([np.floor(p.distance + 0.5) for p, _ in pairs], [years for _, years in pairs])
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=IAP_FIELDS)
        writer.writeheader()
        for (count, years), index in zip(pairs, indices):
            writer.writerow(
                {
                    "seqid_1": count.id_1,
                    "seqid_2": count.id_2,
                    "raw_dist": _count_text(count.distance),
                    "time_diff_years": f"{years:.6f}",
                    "snp_rate": f"{count.distance / years:.6f}" if years else "",
                    "subtype": subtype,
                    "relatedness_index": "" if np.isnan(index) else f"{index:.10g}",
                    "note": NOTE_TOO_FAR_APART if np.isnan(index) else ("" if years else NOTE_SAME_DAY),
                    "compared_sites": count.compared,
                    "ambiguous_sites": count.ambiguous,
                    "flags": count.flags,
                }
            )
            summary.low_coverage_pairs += FLAG_LOW_COVERAGE in count.flags
            summary.noisy_pairs += FLAG_NOISY in count.flags
    summary.unscored_pairs = int(np.isnan(indices).sum())
    return summary


def merge_pairs(tables: Mapping[str, Path], out: Path) -> int:
    """Concatenate the per-subtype tables (in subtype order) into one; return its row count."""
    rows = 0
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=IAP_FIELDS)
        writer.writeheader()
        for _subtype, table in sorted(tables.items()):
            with table.open(newline="", encoding="utf-8") as source:
                for row in csv.DictReader(source):
                    writer.writerow(row)
                    rows += 1
    return rows


def describe(summary: PairSummary, out: Path) -> str:
    text = f"IAP model: {summary.pairs} dated pair(s) written to {out}"
    if summary.skipped_pairs:
        text += (
            f"; {summary.skipped_pairs} pair(s) skipped because {summary.undated_sequences} "
            "sequence(s) have no collection date"
        )
    if summary.unscored_pairs:
        text += f"; {summary.unscored_pairs} pair(s) are too far apart in time for a relatedness index"
    if summary.low_coverage_pairs or summary.noisy_pairs:
        text += (f"; flagged: {summary.low_coverage_pairs} low_coverage, {summary.noisy_pairs} noisy "
                 "(see the flags column; their counts may mislead)")
    return text
