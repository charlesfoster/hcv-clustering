#!/usr/bin/env python3
"""Convert an H2Seq REDCap export into a plotting metadata CSV for hcv-cluster.

Standard library only. Columns are located by header name (first occurrence),
so the export's repeated headers ("Date in:", "Tattooing", ...) are harmless.
Warnings go to stderr and cite the row number and study identifier, never
names or dates of birth.

Values that cannot be derived are written as "Unknown".

Example:
    python h2seq_metadata.py export.csv -o sample_metadata.csv
    python residual_metadata.py residual.csv -o residual.metadata.csv
    python h2seq_metadata.py export.csv --residual-metadata residual.metadata.csv -o sample_metadata.csv
    python h2seq_metadata.py export.csv --location-map locations.csv -o sample_metadata.csv   # locations.csv: CES,Cessnock
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import date, datetime
from pathlib import Path

COL_STUDY_ID = "Study identifier:"
COL_FIRST = "First name initials:"
COL_LAST = "Last name initials:"
COL_DOB = "Date of birth:"
COL_DOB_FALLBACK = "Date of birth"  # interview form's copy, used if the case form's is blank
COL_SEX = "Sex assigned at birth (as per medical record):"
COL_INTERVIEW = "Date of interview (DD/MM/YYYY)"
COL_CENTRE = "Current centre:"
COL_TEST_CENTRE = "HCV test centre:"
COL_MOVED = "Has the participant been incarcerated in any additional centres in the last 6 months?"
INDIGENOUS_PREFIX = "Do you identify as Aboriginal and/or Torres Strait Islander origin? (choice="
COL_ABORIGINAL = INDIGENOUS_PREFIX + "Aboriginal)"
COL_TSI = INDIGENOUS_PREFIX + "Torres Strait Islander)"
COL_NEITHER = INDIGENOUS_PREFIX + "Neither)"
COL_INDIG_UNKNOWN = INDIGENOUS_PREFIX + "Not known/Not stated)"
COL_IDU_EVER = "Have you ever injected drugs?"
COL_IDU_2Y = "Have you injected in the last 2 years?"
COL_IDU_6M = "Have you injected in the last 6 months?"
COL_IDU_1M = "Have you injected in the last month?"
COL_REASON = "Reason for the HCV test"
# The short names used in the cluster report's "How diagnosed" table. A new REDCap option is kept as written (with a warning).
REASON_LABELS = {
    "risk factor triggered screening": "Risk factors",
    "voluntary screening at another time": "Voluntary screening",
    "high intensity testing campaign": "HITC",
    "reception screening": "Reception screening",
    "investigation of symptoms": "Symptomatic",
    "other (please specify)": "Other",
}

DUPLICATE_FIELD = "duplicated"
DUPLICATES_LIST_FIELD = "duplicates"

REQUIRED_COLUMNS = (COL_FIRST, COL_LAST, COL_DOB, COL_SEX, COL_INTERVIEW)

OUTPUT_FIELDS = (
    "sample_id",
    "sample_collection_date",
    "location",
    "test_centre",
    "moved_6m",
    "two_by_two",
    "DOB",
    "sex",
    "age",
    "age_group",
    "age_decade",
    "indigenous_status",
    "injecting_status",
    "how_diagnosed",
    DUPLICATE_FIELD,
    DUPLICATES_LIST_FIELD,
)

DATE_FORMATS = (
    "%d/%m/%Y", "%Y-%m-%d", "%d-%m-%Y", "%d.%m.%Y", "%Y/%m/%d", "%d/%m/%y",
    "%d-%B-%Y", "%d-%b-%Y", "%d %B %Y", "%d %b %Y",
)
TRAILING_TIME = re.compile(r"[ T]\d{1,2}:\d{2}(:\d{2})?$")
UNKNOWN_LABEL = "Unknown"
MONTHS = (
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December",
)
YES = {"yes", "y", "1", "true", "checked"}
NO = {"no", "n", "0", "false", "unchecked"}
UNKNOWN = {"unknown", "not known", "don't know", "dont know", "not stated", "prefer not to say", "declined"}
SEX_CODES = {
    "male": ("M", "Male"),
    "m": ("M", "Male"),
    "female": ("F", "Female"),
    "f": ("F", "Female"),
    "unknown": ("U", "Unknown"),
    "u": ("U", "Unknown"),
}


class RowProblem(ValueError):
    """A value in one row could not be interpreted."""


def parse_date(value: str, column: str) -> date | None:
    value = value.strip()
    if not value:
        return None
    # REDCap datetimes ("2026-01-08 10:30") carry a time we do not need.
    value = TRAILING_TIME.sub("", value)
    for fmt in DATE_FORMATS:
        try:
            return datetime.strptime(value, fmt).date()
        except ValueError:
            continue
    raise RowProblem(f"unrecognised date format in '{column}'")


def parse_yes_no(value: str, column: str) -> bool | None:
    text = value.strip().lower()
    if not text or text in UNKNOWN:
        return None
    if text in YES:
        return True
    if text in NO:
        return False
    raise RowProblem(f"unrecognised yes/no value {value.strip()!r} in '{column}'")


def compact_date(value: date) -> str:
    return value.strftime("%d%m%Y")


def long_date(value: date) -> str:
    # Explicit month names keep the output independent of the system locale.
    return f"{value.day:02d}-{MONTHS[value.month - 1]}-{value.year}"


def initials(value: str) -> str:
    return "".join(char for char in value.upper() if char.isalpha())


def age_at(dob: date, on: date) -> int:
    return on.year - dob.year - ((on.month, on.day) < (dob.month, dob.day))


def age_group(age: int) -> str:
    if age <= 30:
        return "0-30"
    if age <= 60:
        return "31-60"
    return "61+"


def age_decade(age: int) -> str:
    start = (age // 10) * 10
    return f"{start}-{start + 9}"


def indigenous_status(get) -> str:
    def checked(column: str) -> bool:
        return bool(parse_yes_no(get(column), column))

    aboriginal, tsi = checked(COL_ABORIGINAL), checked(COL_TSI)
    if aboriginal and tsi:
        return "Aboriginal and Torres Strait Islander"
    if aboriginal:
        return "Aboriginal"
    if tsi:
        return "Torres Strait Islander"
    if checked(COL_NEITHER):
        return "Neither"
    return UNKNOWN_LABEL


def injecting_status(get) -> str:
    ever = parse_yes_no(get(COL_IDU_EVER), COL_IDU_EVER)
    if ever is False:
        return "Never"
    # Most recent period answered "yes" wins; later questions are skipped when earlier ones are "no".
    for column, label in ((COL_IDU_1M, "Last month"), (COL_IDU_6M, "Last 6 months"), (COL_IDU_2Y, "Last 2 years")):
        if parse_yes_no(get(column), column):
            return label
    if parse_yes_no(get(COL_IDU_2Y), COL_IDU_2Y) is False:
        return "More than 2 years ago"
    if ever:
        return "Ever (recency unknown)"
    return UNKNOWN_LABEL


def how_diagnosed(get, notes: list[str]) -> str:
    raw = get(COL_REASON).strip()
    if not raw:
        return UNKNOWN_LABEL
    short = REASON_LABELS.get(raw.lower())
    if short is None:
        notes.append(f"unrecognised reason for the HCV test {raw!r}; kept as written")
        return raw
    return short


def soft(derive, notes: list[str]) -> str:
    """Derive an optional field; a bad value makes it Unknown instead of dropping the row."""
    try:
        return derive() or UNKNOWN_LABEL
    except RowProblem as exc:
        notes.append(f"{exc}; set to {UNKNOWN_LABEL}")
        return UNKNOWN_LABEL


def yes_no_label(value: bool | None) -> str:
    return UNKNOWN_LABEL if value is None else ("Yes" if value else "No")


def convert_row(get, notes: list[str]) -> dict[str, str]:
    first, last = initials(get(COL_FIRST)), initials(get(COL_LAST))
    dob = parse_date(get(COL_DOB), COL_DOB) or parse_date(get(COL_DOB_FALLBACK), COL_DOB_FALLBACK)
    interview = parse_date(get(COL_INTERVIEW), COL_INTERVIEW)
    sex_raw = get(COL_SEX).strip()
    sex = SEX_CODES.get(sex_raw.lower())

    missing = [
        name
        for name, present in (
            ("first name initials", first),
            ("last name initials", last),
            ("date of birth", dob),
            ("sex", sex_raw),
            ("interview date", interview),
        )
        if not present
    ]
    if missing:
        raise RowProblem("cannot build sample_id; missing " + ", ".join(missing))
    if sex is None:
        raise RowProblem(f"unrecognised sex value {sex_raw!r}; expected Male/Female/Unknown or M/F/U")
    if interview < dob:
        raise RowProblem("interview date is before date of birth")

    age = age_at(dob, interview)
    return {
        "sample_id": f"{first}{last}{compact_date(dob)}{sex[0]}{compact_date(interview)}",
        "sample_collection_date": long_date(interview),
        "location": get(COL_CENTRE).strip() or UNKNOWN_LABEL,
        "test_centre": get(COL_TEST_CENTRE).strip() or UNKNOWN_LABEL,
        "moved_6m": soft(lambda: yes_no_label(parse_yes_no(get(COL_MOVED), COL_MOVED)), notes),
        "two_by_two": first + last,
        "DOB": long_date(dob),
        "sex": sex[1],
        "age": str(age),
        "age_group": age_group(age),
        "age_decade": age_decade(age),
        "indigenous_status": soft(lambda: indigenous_status(get), notes),
        "injecting_status": soft(lambda: injecting_status(get), notes),
        "how_diagnosed": how_diagnosed(get, notes),
    }


def read_export(input_path: Path, encoding: str) -> tuple[list[dict[str, str]], int, list[str]]:
    """Convert REDCap export rows; returns (records, skipped row count, warnings)."""
    warnings: list[str] = []
    with input_path.open("r", encoding=encoding, newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = [name.strip() for name in next(reader)]
        except StopIteration:
            raise SystemExit(f"error: {input_path} is empty")

        index: dict[str, int] = {}
        for position, name in enumerate(header):
            index.setdefault(name, position)
        absent = [name for name in REQUIRED_COLUMNS if name not in index]
        if absent:
            raise SystemExit("error: required column(s) not found: " + "; ".join(absent))
        optional_absent = [
            name
            for name in (COL_CENTRE, COL_TEST_CENTRE, COL_MOVED, COL_ABORIGINAL, COL_TSI, COL_NEITHER,
                         COL_INDIG_UNKNOWN, COL_IDU_EVER, COL_IDU_2Y, COL_IDU_6M, COL_IDU_1M, COL_REASON)
            if name not in index
        ]
        for name in optional_absent:
            warnings.append(f"column not found, output will be {UNKNOWN_LABEL}: {name}")

        converted: list[tuple[dict[str, str], int]] = []
        # sample_id -> [(row number, warning label, identifier)]; identifier is the study ID or "row N"
        rows_by_id: dict[str, list[tuple[int, str, str]]] = {}
        skipped = 0
        for row_number, row in enumerate(reader, start=2):
            if not any(cell.strip() for cell in row):
                continue

            def get(column: str, row=row) -> str:
                position = index.get(column)
                return row[position] if position is not None and position < len(row) else ""

            label = f"row {row_number}"
            study_id = get(COL_STUDY_ID).strip()
            if study_id:
                label += f" ({study_id})"
            notes: list[str] = []
            try:
                record = convert_row(get, notes)
            except RowProblem as exc:
                warnings.append(f"{label}: skipped: {exc}")
                skipped += 1
                continue
            warnings.extend(f"{label}: {note}" for note in notes)
            identifier = study_id or f"row {row_number}"
            rows_by_id.setdefault(record["sample_id"], []).append((row_number, label, identifier))
            converted.append((record, row_number))

    for sample_id, entries in rows_by_id.items():
        if len(entries) > 1:
            warnings.append(
                f"duplicate sample_id {sample_id} on {len(entries)} rows: "
                + ", ".join(label for _, label, _ in entries)
            )

    for record, row_number in converted:
        others = [other for number, _, other in rows_by_id[record["sample_id"]] if number != row_number]
        record[DUPLICATE_FIELD] = "true" if len(rows_by_id[record["sample_id"]]) > 1 else "false"
        record[DUPLICATES_LIST_FIELD] = ";".join(others)
    return [record for record, _ in converted], skipped, warnings


def read_formatted(input_path: Path, encoding: str) -> list[dict[str, str]]:
    """Read a CSV already in OUTPUT_FIELDS form (e.g. from residual_metadata.py) unchanged."""
    with input_path.open("r", encoding=encoding, newline="") as handle:
        reader = csv.DictReader(handle)
        fields = [name.strip() for name in reader.fieldnames or ()]
        missing = [name for name in OUTPUT_FIELDS if name not in fields and name != "how_diagnosed"]  # older files lack it
        if missing:
            raise SystemExit(
                f"error: {input_path} is missing column(s): {', '.join(missing)}. "
                "Reformat it first with residual_metadata.py."
            )
        reader.fieldnames = fields
        return [
            {name: (row.get(name) or "").strip() or (UNKNOWN_LABEL if name == "how_diagnosed" else "") for name in OUTPUT_FIELDS}
            for row in reader
            if any((cell or "").strip() for cell in row.values() if isinstance(cell, str))
        ]


def read_location_map(path: Path, encoding: str) -> dict[str, str]:
    """Read a two-column CSV of centre code -> location name (a header row "code,location" is optional)."""
    mapping: dict[str, str] = {}
    with path.open("r", encoding=encoding, newline="") as handle:
        for number, row in enumerate(csv.reader(handle), start=1):
            cells = [cell.strip() for cell in row]
            if not any(cells) or (number == 1 and cells[0].lower() == "code"):
                continue
            if len(cells) < 2 or not cells[0] or not cells[1]:
                raise SystemExit(f"error: {path} line {number}: expected 'code,location'")
            key = cells[0].lower()
            if mapping.get(key, cells[1]) != cells[1]:
                raise SystemExit(f"error: {path} line {number}: {cells[0]!r} is given two different names")
            mapping[key] = cells[1]
    return mapping


def apply_location_map(records: list[dict[str, str]], mapping: dict[str, str], warnings: list[str]) -> None:
    """Rename the location and test_centre values in place; codes the map lacks are kept and listed once."""
    unmapped: set[str] = set()
    for record in records:
        for field in ("location", "test_centre"):
            value = record.get(field, "")
            if value.lower() in mapping:
                record[field] = mapping[value.lower()]
            elif value and value != UNKNOWN_LABEL and value not in mapping.values():
                unmapped.add(value)
    if unmapped:
        warnings.append("not in --location-map, left as written: " + ", ".join(sorted(unmapped)))


def write_records(records: list[dict[str, str]], output) -> None:
    writer = csv.DictWriter(output, fieldnames=OUTPUT_FIELDS)
    writer.writeheader()
    writer.writerows(records)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert an H2Seq REDCap export into a sample metadata CSV for hcv-cluster plots.",
    )
    parser.add_argument("input", type=Path, help="REDCap export CSV (labels or raw codes)")
    parser.add_argument(
        "-o", "--output", type=Path,
        help="Output CSV (default: <input stem>.metadata.csv next to the input; '-' for stdout)",
    )
    parser.add_argument(
        "--residual-metadata", type=Path, metavar="CSV",
        help="Output of residual_metadata.py (already in this script's columns); its rows are written first",
    )
    parser.add_argument(
        "--location-map", type=Path, metavar="CSV",
        help="Two columns, centre code then location name (e.g. CES,Cessnock), applied to location and test_centre; "
             "codes not in the file are kept and listed in a warning",
    )
    parser.add_argument("--encoding", default="utf-8-sig", help="Input encoding (default: utf-8-sig; try cp1252 if it fails)")
    parser.add_argument("--strict", action="store_true", help="Exit with status 1 if any row was skipped or warned about")
    parser.add_argument("-f", "--force", action="store_true", help="Overwrite the output file if it exists")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    for path in (args.input, args.residual_metadata, args.location_map):
        if path is not None and not path.is_file():
            print(f"error: input not found: {path}", file=sys.stderr)
            return 2

    output_path = args.output or args.input.with_name(f"{args.input.stem}.metadata.csv")
    to_stdout = str(output_path) == "-"
    if not to_stdout and output_path.exists() and not args.force:
        print(f"error: {output_path} exists; pass --force to overwrite", file=sys.stderr)
        return 2

    records, skipped, warnings = read_export(args.input, args.encoding)
    if args.residual_metadata is not None:
        residual = read_formatted(args.residual_metadata, args.encoding)
        overlap = sorted({row["sample_id"] for row in residual} & {row["sample_id"] for row in records})
        if overlap:
            warnings.append("sample_id(s) in both residual and export: " + ", ".join(overlap))
        records = residual + records

    if args.location_map is not None:
        apply_location_map(records, read_location_map(args.location_map, args.encoding), warnings)

    if to_stdout:
        write_records(records, sys.stdout)
    else:
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            write_records(records, handle)

    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    destination = "stdout" if to_stdout else output_path
    print(f"Wrote {len(records)} row(s) to {destination}; skipped {skipped}.", file=sys.stderr)
    return 1 if args.strict and warnings else 0


if __name__ == "__main__":
    sys.exit(main())
