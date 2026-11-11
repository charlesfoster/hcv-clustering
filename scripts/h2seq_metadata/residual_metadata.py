#!/usr/bin/env python3
"""Reformat a residual-sample CSV into the h2seq_metadata.py output columns.

Expected input columns (extra columns are ignored):
    sample_id, DOC, location, two_by_two, DOB, gender, duplicates
    (age is used only when it cannot be recalculated from DOB and DOC)

The supplied sample_id and duplicates are kept as given; age and its bins are
recalculated, and fields the residual sheet cannot provide are "Unknown"
(except age, which is left blank because hcv-cluster requires a number there).

Example:
    python residual_metadata.py residual.csv -o residual.metadata.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from h2seq_metadata import (
    DUPLICATE_FIELD,
    DUPLICATES_LIST_FIELD,
    SEX_CODES,
    UNKNOWN_LABEL,
    RowProblem,
    age_at,
    age_decade,
    age_group,
    long_date,
    parse_date,
    write_records,
)


def convert_residual_row(row: dict[str, str], notes: list[str]) -> dict[str, str]:
    def value(column: str) -> str:
        return (row.get(column) or "").strip()

    def date_field(column: str):
        """Return (parsed date or None, output text); unreadable dates are kept verbatim."""
        raw = value(column)
        try:
            parsed = parse_date(raw, column)
        except RowProblem as exc:
            notes.append(f"{exc}; kept as given")
            return None, raw
        return parsed, long_date(parsed) if parsed else UNKNOWN_LABEL

    doc, doc_text = date_field("DOC")
    dob, dob_text = date_field("DOB")

    age: int | None = None
    if doc and dob and doc >= dob:
        age = age_at(dob, doc)
    elif value("age").isdigit():
        age = int(value("age"))
        notes.append("age not recalculable from DOB/DOC; using supplied age")
    elif doc and dob:
        notes.append("DOC is before DOB; age set to Unknown")

    gender = value("gender")
    sex = SEX_CODES.get(gender.lower())
    if sex is None and gender:
        notes.append(f"unrecognised gender value {gender!r}; set to {UNKNOWN_LABEL}")

    duplicates = value("duplicates")
    return {
        "sample_id": value("sample_id"),
        "sample_collection_date": doc_text,
        "location": value("location") or UNKNOWN_LABEL,
        "test_centre": UNKNOWN_LABEL,
        "moved_6m": UNKNOWN_LABEL,
        "two_by_two": value("two_by_two") or UNKNOWN_LABEL,
        "DOB": dob_text,
        "sex": sex[1] if sex else UNKNOWN_LABEL,
        # hcv-cluster validates age as a whole number, so an unknown age stays blank.
        "age": "" if age is None else str(age),
        "age_group": UNKNOWN_LABEL if age is None else age_group(age),
        "age_decade": UNKNOWN_LABEL if age is None else age_decade(age),
        "indigenous_status": UNKNOWN_LABEL,
        "injecting_status": UNKNOWN_LABEL,
        DUPLICATE_FIELD: "true" if duplicates else "false",
        DUPLICATES_LIST_FIELD: duplicates,
    }


def read_residual(input_path: Path, encoding: str) -> tuple[list[dict[str, str]], int, list[str]]:
    """Convert residual rows; returns (records, skipped row count, warnings)."""
    warnings: list[str] = []
    records: list[dict[str, str]] = []
    skipped = 0
    with input_path.open("r", encoding=encoding, newline="") as handle:
        reader = csv.DictReader(handle)
        reader.fieldnames = [name.strip() for name in reader.fieldnames or ()]
        if "sample_id" not in reader.fieldnames:
            raise SystemExit(f"error: {input_path} has no 'sample_id' column")
        for name in ("DOC", "location", "two_by_two", "DOB", "gender", "duplicates"):
            if name not in reader.fieldnames:
                warnings.append(f"column not found, output will be {UNKNOWN_LABEL}: {name}")

        seen: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            if not any((cell or "").strip() for cell in row.values() if isinstance(cell, str)):
                continue
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                warnings.append(f"row {row_number}: skipped: blank sample_id")
                skipped += 1
                continue
            if sample_id in seen:
                warnings.append(f"row {row_number} ({sample_id}): sample_id repeated in residual file")
            seen.add(sample_id)
            notes: list[str] = []
            records.append(convert_residual_row(row, notes))
            warnings.extend(f"row {row_number} ({sample_id}): {note}" for note in notes)
    return records, skipped, warnings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Reformat a residual-sample CSV into h2seq_metadata.py columns.")
    parser.add_argument("input", type=Path, help="Residual-sample CSV")
    parser.add_argument(
        "-o", "--output", type=Path,
        help="Output CSV (default: <input stem>.metadata.csv next to the input; '-' for stdout)",
    )
    parser.add_argument("--encoding", default="utf-8-sig", help="Input encoding (default: utf-8-sig)")
    parser.add_argument("-f", "--force", action="store_true", help="Overwrite the output file if it exists")
    args = parser.parse_args(argv)

    if not args.input.is_file():
        print(f"error: input not found: {args.input}", file=sys.stderr)
        return 2
    output_path = args.output or args.input.with_name(f"{args.input.stem}.metadata.csv")
    to_stdout = str(output_path) == "-"
    if not to_stdout and output_path.exists() and not args.force:
        print(f"error: {output_path} exists; pass --force to overwrite", file=sys.stderr)
        return 2

    records, skipped, warnings = read_residual(args.input, args.encoding)
    if to_stdout:
        write_records(records, sys.stdout)
    else:
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            write_records(records, handle)
    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    destination = "stdout" if to_stdout else output_path
    print(f"Wrote {len(records)} row(s) to {destination}; skipped {skipped}.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
