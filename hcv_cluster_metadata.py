#!/usr/bin/env python3
"""Metadata loading and non-destructive joins for HCV network plots."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


MISSING_VALUE = "(missing)"
AGE_RANGE_FIELD = "age_range"
AGE_RANGE_ORDER = ("0–30", "31–60", "61+")
AGE_RANGE_ALIASES = {
    "0-30": AGE_RANGE_ORDER[0],
    "0–30": AGE_RANGE_ORDER[0],
    "31-60": AGE_RANGE_ORDER[1],
    "31–60": AGE_RANGE_ORDER[1],
    "61+": AGE_RANGE_ORDER[2],
}
RESERVED_METADATA_COLUMNS = frozenset(
    {
        "genotype",
        "cluster_id",
        "source_cluster_id",
        "cluster_size",
        "source",
        "target",
        "distance",
    }
)


class MetadataValidationError(ValueError):
    """Raised when metadata cannot be joined to clustering results safely."""


@dataclass(frozen=True)
class MetadataJoinReport:
    """Sample identifiers that did not match during a metadata-to-node join."""

    missing_metadata_sample_ids: tuple[str, ...]
    unknown_metadata_sample_ids: tuple[str, ...]

    @property
    def has_warnings(self) -> bool:
        return bool(self.missing_metadata_sample_ids or self.unknown_metadata_sample_ids)


def age_to_range(age: Any) -> str:
    """Validate an age and return its ordered plotting bin."""
    if age is None or (isinstance(age, str) and age.strip() in ("", MISSING_VALUE)):
        return MISSING_VALUE
    if isinstance(age, bool):
        raise MetadataValidationError(f"Invalid age {age!r}; age must be a non-negative whole number")
    try:
        numeric_age = float(age)
    except (TypeError, ValueError) as exc:
        raise MetadataValidationError(
            f"Invalid age {age!r}; age must be a non-negative whole number"
        ) from exc
    if not math.isfinite(numeric_age) or numeric_age < 0 or not numeric_age.is_integer():
        raise MetadataValidationError(
            f"Invalid age {age!r}; age must be a non-negative whole number"
        )
    integer_age = int(numeric_age)
    if integer_age <= 30:
        return AGE_RANGE_ORDER[0]
    if integer_age <= 60:
        return AGE_RANGE_ORDER[1]
    return AGE_RANGE_ORDER[2]


def _clean_fieldnames(fieldnames: Sequence[str | None] | None) -> list[str]:
    if not fieldnames:
        raise MetadataValidationError("Metadata CSV has no header row")
    fields = [field.strip() if field is not None else "" for field in fieldnames]
    if any(not field for field in fields):
        raise MetadataValidationError("Metadata CSV contains an empty column name")
    duplicates = sorted({field for field in fields if fields.count(field) > 1})
    if duplicates:
        raise MetadataValidationError(f"Duplicate metadata column name(s): {', '.join(duplicates)}")
    if "sample_id" not in fields:
        raise MetadataValidationError("Metadata CSV requires a 'sample_id' column")
    collisions = sorted(RESERVED_METADATA_COLUMNS.intersection(fields))
    if collisions:
        raise MetadataValidationError(
            "Metadata column name(s) reserved for clustering output: " + ", ".join(collisions)
        )
    return fields


def normalize_metadata_rows(
    rows: Iterable[Mapping[str, Any]],
    *,
    fieldnames: Sequence[str] | None = None,
) -> list[dict[str, Any]]:
    """Validate and normalize metadata rows, including the derived age range.

    Empty values are represented explicitly as ``(missing)``. Arbitrary columns are
    retained, while clustering-output names are rejected to prevent accidental
    overwrites during a join.
    """
    source_rows = [dict(row) for row in rows]
    if fieldnames is None:
        if not source_rows:
            raise MetadataValidationError("Metadata has no rows and no column definition")
        fieldnames = list(dict.fromkeys(field for row in source_rows for field in row))
    fields = _clean_fieldnames(fieldnames)
    if "age" in fields and AGE_RANGE_FIELD not in fields:
        fields.append(AGE_RANGE_FIELD)

    normalized: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    duplicate_ids: set[str] = set()
    errors: list[str] = []

    for row_number, source in enumerate(source_rows, start=2):
        sample_id = str(source.get("sample_id") or "").strip()
        if not sample_id:
            errors.append(f"row {row_number}: sample_id is blank")
            continue
        if sample_id in seen_ids:
            duplicate_ids.add(sample_id)
        seen_ids.add(sample_id)

        row: dict[str, Any] = {"sample_id": sample_id}
        for field in fields:
            if field == "sample_id":
                continue
            value = source.get(field)
            if value is None or (
                isinstance(value, str) and value.strip() in ("", MISSING_VALUE)
            ):
                row[field] = MISSING_VALUE
            else:
                row[field] = value.strip() if isinstance(value, str) else value

        supplied_range = row.get(AGE_RANGE_FIELD, MISSING_VALUE)
        if supplied_range != MISSING_VALUE:
            canonical_range = AGE_RANGE_ALIASES.get(str(supplied_range))
            if canonical_range is None:
                errors.append(
                    f"row {row_number} ({sample_id}): invalid age_range {supplied_range!r}; "
                    "expected 0-30, 31-60, or 61+"
                )
            else:
                supplied_range = canonical_range
                row[AGE_RANGE_FIELD] = canonical_range

        if "age" in fields:
            try:
                inferred_range = age_to_range(row.get("age"))
            except MetadataValidationError as exc:
                errors.append(f"row {row_number} ({sample_id}): {exc}")
            else:
                if (
                    supplied_range != MISSING_VALUE
                    and inferred_range != MISSING_VALUE
                    and supplied_range != inferred_range
                ):
                    errors.append(
                        f"row {row_number} ({sample_id}): age and age_range are inconsistent "
                        f"({inferred_range!r} != {supplied_range!r})"
                    )
                elif supplied_range == MISSING_VALUE:
                    row[AGE_RANGE_FIELD] = inferred_range
        normalized.append(row)

    if duplicate_ids:
        errors.append("duplicate sample_id value(s): " + ", ".join(sorted(duplicate_ids)))
    if errors:
        raise MetadataValidationError("Invalid metadata:\n- " + "\n- ".join(errors))
    return normalized


def load_metadata_csv(path: str | Path) -> list[dict[str, Any]]:
    """Load and validate a UTF-8 CSV metadata file."""
    metadata_path = Path(path)
    try:
        with metadata_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = _clean_fieldnames(reader.fieldnames)
            # DictReader retains the original header spelling; headers are stripped
            # here so benign surrounding whitespace does not change field names.
            rows = []
            for row_number, source in enumerate(reader, start=2):
                if None in source:
                    raise MetadataValidationError(
                        f"Metadata row {row_number} has more values than the header"
                    )
                rows.append(
                    {
                        cleaned: source.get(original)
                        for original, cleaned in zip(reader.fieldnames or (), fields, strict=True)
                    }
                )
    except (OSError, UnicodeError, csv.Error) as exc:
        raise MetadataValidationError(f"Could not read metadata CSV {metadata_path}: {exc}") from exc
    return normalize_metadata_rows(rows, fieldnames=fields)


def metadata_fields(metadata_rows: Sequence[Mapping[str, Any]]) -> tuple[str, ...]:
    """Return metadata fields in stable input order, excluding the join key."""
    fields: list[str] = []
    for row in metadata_rows:
        for field in row:
            if field != "sample_id" and field not in fields:
                fields.append(field)
    return tuple(fields)


def write_metadata_csv(metadata_rows: Sequence[Mapping[str, Any]], path: str | Path) -> Path:
    """Write normalized metadata as UTF-8 CSV, preserving first-seen field order."""
    output_path = Path(path)
    fields = ("sample_id", *metadata_fields(metadata_rows))
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in metadata_rows:
            writer.writerow({field: row.get(field, MISSING_VALUE) for field in fields})
    return output_path


def join_metadata(
    node_rows: Sequence[Mapping[str, Any]],
    metadata_rows: Sequence[Mapping[str, Any]],
    *,
    missing_value: str = MISSING_VALUE,
) -> tuple[list[dict[str, Any]], MetadataJoinReport]:
    """Left-join metadata into node rows without changing graph membership.

    Unknown metadata IDs and graph nodes lacking metadata are returned in a report.
    Existing node values are never silently overwritten with a conflicting value.
    """
    fields = metadata_fields(metadata_rows)
    index: dict[str, Mapping[str, Any]] = {}
    for row in metadata_rows:
        sample_id = str(row.get("sample_id") or "").strip()
        if not sample_id:
            raise MetadataValidationError("Metadata row has a blank sample_id")
        if sample_id in index:
            raise MetadataValidationError(f"Duplicate metadata sample_id: {sample_id}")
        index[sample_id] = row

    node_ids = {str(row["sample_id"]) for row in node_rows}
    joined: list[dict[str, Any]] = []
    missing_ids: list[str] = []
    for source in node_rows:
        row = dict(source)
        sample_id = str(row["sample_id"])
        metadata = index.get(sample_id)
        if metadata is None:
            missing_ids.append(sample_id)
            for field in fields:
                row.setdefault(field, missing_value)
        else:
            for field in fields:
                value = metadata.get(field, missing_value)
                if field in row and row[field] not in (None, "", missing_value) and row[field] != value:
                    raise MetadataValidationError(
                        f"Metadata for {sample_id!r} conflicts with existing node field {field!r}"
                    )
                row[field] = value
        joined.append(row)

    report = MetadataJoinReport(
        missing_metadata_sample_ids=tuple(sorted(missing_ids)),
        unknown_metadata_sample_ids=tuple(sorted(set(index).difference(node_ids))),
    )
    return joined, report
