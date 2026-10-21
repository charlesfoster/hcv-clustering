from pathlib import Path

import pytest

from hcv_cluster_metadata import (
    AGE_RANGE_ORDER,
    MISSING_VALUE,
    MetadataValidationError,
    age_to_range,
    join_metadata,
    load_metadata_csv,
    normalize_metadata_rows,
    write_metadata_csv,
)


@pytest.mark.parametrize(
    ("age", "expected"),
    [(0, "0–30"), (30, "0–30"), (31, "31–60"), (60, "31–60"), (61, "61+"), ("", MISSING_VALUE)],
)
def test_age_boundaries(age: object, expected: str) -> None:
    assert age_to_range(age) == expected


@pytest.mark.parametrize("age", ["unknown", -1, 3.5, float("inf")])
def test_invalid_age_is_rejected(age: object) -> None:
    with pytest.raises(MetadataValidationError, match="Invalid age"):
        age_to_range(age)


def test_normalization_preserves_arbitrary_columns_and_derives_age_range() -> None:
    rows = normalize_metadata_rows(
        [
            {"sample_id": " A ", "age": "30", "location": "P1", "subtype": ""},
            {"sample_id": "B", "age": "61", "location": "P2", "indigenous": "yes"},
        ]
    )
    assert rows[0] == {
        "sample_id": "A", "age": "30", "location": "P1", "subtype": MISSING_VALUE,
        "indigenous": MISSING_VALUE, "age_range": "0–30",
    }
    assert rows[1]["age_range"] == "61+"
    assert tuple(value for value in AGE_RANGE_ORDER) == ("0–30", "31–60", "61+")


def test_supplied_age_range_is_validated_and_checked_for_consistency() -> None:
    with pytest.raises(MetadataValidationError, match="inconsistent"):
        normalize_metadata_rows([{"sample_id": "A", "age": "20", "age_range": "61+"}])

    rows = normalize_metadata_rows(
        [
            {"sample_id": "A", "age_range": "0-30"},
            {"sample_id": "B", "age_range": "31-60"},
        ]
    )
    assert [row["age_range"] for row in rows] == ["0–30", "31–60"]

    with pytest.raises(MetadataValidationError, match="invalid age_range"):
        normalize_metadata_rows([{"sample_id": "A", "age_range": "young"}])


@pytest.mark.parametrize(
    "field",
    [
        "genotype",
        "cluster_id",
        "source_cluster_id",
        "cluster_size",
        "source",
        "target",
        "distance",
    ],
)
def test_reserved_columns_are_rejected(field: str) -> None:
    with pytest.raises(MetadataValidationError, match="reserved"):
        normalize_metadata_rows([{"sample_id": "A", field: "x"}])


def test_csv_requires_key_and_unique_ids(tmp_path: Path) -> None:
    missing_key = tmp_path / "missing.csv"
    missing_key.write_text("location\nP1\n", encoding="utf-8")
    with pytest.raises(MetadataValidationError, match="sample_id"):
        load_metadata_csv(missing_key)

    duplicate = tmp_path / "duplicate.csv"
    duplicate.write_text("sample_id,location\nA,P1\nA,P2\n", encoding="utf-8")
    with pytest.raises(MetadataValidationError, match="duplicate sample_id"):
        load_metadata_csv(duplicate)


def test_csv_rejects_rows_with_more_values_than_the_header(tmp_path: Path) -> None:
    metadata = tmp_path / "extra_value.csv"
    metadata.write_text("sample_id,location\nA,P1,unexpected\n", encoding="utf-8")

    with pytest.raises(MetadataValidationError, match="more values than the header"):
        load_metadata_csv(metadata)


def test_join_is_left_join_and_reports_both_unmatched_directions() -> None:
    nodes = [
        {"sample_id": "A", "cluster_id": "C1", "cluster_size": "1"},
        {"sample_id": "B", "cluster_id": "C2", "cluster_size": "1"},
    ]
    metadata = normalize_metadata_rows(
        [{"sample_id": "A", "location": "P1"}, {"sample_id": "UNKNOWN", "location": "P9"}]
    )
    joined, report = join_metadata(nodes, metadata)
    assert [row["sample_id"] for row in joined] == ["A", "B"]
    assert joined[0]["location"] == "P1"
    assert joined[1]["location"] == MISSING_VALUE
    assert report.missing_metadata_sample_ids == ("B",)
    assert report.unknown_metadata_sample_ids == ("UNKNOWN",)
    assert report.has_warnings


def test_join_never_silently_overwrites_existing_node_data() -> None:
    nodes = [{"sample_id": "A", "cluster_id": "C1", "cluster_size": "1", "location": "P1"}]
    metadata = normalize_metadata_rows([{"sample_id": "A", "location": "P2"}])
    with pytest.raises(MetadataValidationError, match="conflicts"):
        join_metadata(nodes, metadata)


def test_write_metadata_csv_preserves_field_order_and_utf8(tmp_path: Path) -> None:
    rows = normalize_metadata_rows(
        [{"sample_id": "A", "location": "Métro", "status": "yes"}],
    )
    path = write_metadata_csv(rows, tmp_path / "normalized.csv")
    assert path.read_text(encoding="utf-8").splitlines() == [
        "sample_id,location,status", "A,Métro,yes",
    ]


def test_normalized_metadata_with_missing_age_round_trips(tmp_path: Path) -> None:
    rows = normalize_metadata_rows([{"sample_id": "A", "age": "", "location": ""}])
    path = write_metadata_csv(rows, tmp_path / "metadata.csv")
    assert load_metadata_csv(path) == rows
