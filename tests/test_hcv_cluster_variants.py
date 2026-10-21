import csv
from pathlib import Path

import pytest

import hcv_cluster_metadata
import hcv_cluster_prep
import hcv_cluster_variants as variants


def _records(*pairs: tuple[str, str]) -> list[hcv_cluster_prep.FastaRecord]:
    return [hcv_cluster_prep.FastaRecord(header, sequence) for header, sequence in pairs]


def test_identical_repeats_are_dropped_and_distinct_sequences_renamed() -> None:
    plan = variants.resolve_duplicate_ids(
        _records(("a", "ACGT"), ("b", "ACGT"), ("a", "acgt"), ("c", "AAAA"), ("c", "CCCC"))
    )
    assert [record.header for record in plan.records] == ["a", "b", "c__v1", "c__v2"]
    assert plan.mixed_sample_ids() == {"c"}
    assert plan.base_sample_ids()["c__v2"] == "c"
    statuses = {variant.input_record: variant.status for variant in plan.variants}
    assert statuses[3] == "dropped_identical_repeat"


def test_unique_input_is_left_untouched() -> None:
    plan = variants.resolve_duplicate_ids(_records(("a", "ACGT"), ("b", "ACGT")))
    assert not plan.has_duplicates
    assert [record.header for record in plan.records] == ["a", "b"]


def test_rename_refuses_to_collide_with_an_existing_id() -> None:
    with pytest.raises(ValueError, match="already exists"):
        variants.resolve_duplicate_ids(_records(("c", "AAAA"), ("c", "CCCC"), ("c__v1", "GGGG")))


def _qc(path: Path, coverage: dict[str, float]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("sample_id", "coverage_fraction"))
        writer.writeheader()
        for sample_id, value in coverage.items():
            writer.writerow({"sample_id": sample_id, "coverage_fraction": f"{value:.6f}"})
    return path


def _same_genotype_plan() -> variants.VariantPlan:
    # v1: truncated but covers the region fully; v2: more complete with Ns inside the region.
    plan = variants.resolve_duplicate_ids(_records(("s", "ACGT" * 10), ("s", "ACGT" * 20 + "N" * 4)))
    variants.record_genotypes(
        plan,
        [
            {"sequence_id": "s__v1", "assigned_genotype": "1a", "assignment_status": "pass"},
            {"sequence_id": "s__v2", "assigned_genotype": "1a", "assignment_status": "pass"},
        ],
    )
    return plan


@pytest.mark.parametrize(
    ("mode", "coverage", "expected_kept"),
    [
        ("region-coverage", {"s__v1": 1.0, "s__v2": 0.8}, "s__v1"),
        ("region-coverage", {"s__v1": 0.9, "s__v2": 0.9}, "s__v2"),  # tie -> most complete
        ("completeness", {"s__v1": 1.0, "s__v2": 0.8}, "s__v2"),
    ],
)
def test_same_genotype_selection(tmp_path: Path, mode: str, coverage: dict[str, float], expected_kept: str) -> None:
    plan = _same_genotype_plan()
    dropped = variants.select_same_genotype(plan, "1a", _qc(tmp_path / "qc.csv", coverage), mode)
    kept = [v.sequence_id for v in plan.variants if v.status == "kept"]
    assert kept == [expected_kept]
    assert len(dropped) == 1


def test_metadata_joins_on_base_sample_id() -> None:
    nodes = [
        {"sample_id": "s__v1", "base_sample_id": "s", "cluster_id": "C1", "cluster_size": 1},
        {"sample_id": "s__v2", "base_sample_id": "s", "cluster_id": "C2", "cluster_size": 1},
    ]
    joined, report = hcv_cluster_metadata.join_metadata(nodes, [{"sample_id": "s", "location": "P1"}])
    assert [row["location"] for row in joined] == ["P1", "P1"]
    assert report.missing_metadata_sample_ids == ()
    assert report.unknown_metadata_sample_ids == ()
