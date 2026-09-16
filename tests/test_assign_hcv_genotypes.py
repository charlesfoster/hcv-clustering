from pathlib import Path
from types import SimpleNamespace

import assign_hcv_genotypes_from_fasta as genotyping
import pytest


def test_normalize_records_strips_gaps_before_length_is_measured() -> None:
    records = {
        "seq1": {"description": "seq1", "sequence": "ACGT--ACGT"},
    }

    normalized = genotyping.normalize_records(records)

    assert normalized["seq1"]["sequence"] == "ACGTACGT"


def test_normalize_records_uppercases_and_converts_u_to_t() -> None:
    records = {
        "seq1": {"description": "seq1", "sequence": "acgu"},
    }

    normalized = genotyping.normalize_records(records)

    assert normalized["seq1"]["sequence"] == "ACGT"


def test_normalize_records_all_gap_sequence_becomes_empty_not_a_crash() -> None:
    records = {
        "seq1": {"description": "seq1", "sequence": "----"},
        "seq2": {"description": "seq2", "sequence": "ACGT"},
    }

    normalized = genotyping.normalize_records(records)

    assert normalized["seq1"]["sequence"] == ""
    assert normalized["seq2"]["sequence"] == "ACGT"


def test_write_records_fasta_round_trips_normalized_sequence(tmp_path: Path) -> None:
    records = genotyping.normalize_records(
        {"seq1": {"description": "seq1 some description", "sequence": "ACGT--NNNNacgt"}}
    )

    fasta_path = tmp_path / "query.fasta"
    genotyping.write_records_fasta(fasta_path, records)

    reparsed = genotyping.parse_fasta(fasta_path)
    assert reparsed["seq1"]["sequence"] == "ACGTNNNNACGT"


def _write_split_paf(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "sample\t9453\t4934\t9307\t+\t3a_D17763.1\t9456\t4935\t9308\t3955\t4253\t60\tAS:i:1153\ttp:A:P",
                "sample\t9453\t501\t1409\t+\t3a_D17763.1\t9456\t501\t1409\t851\t907\t60\tAS:i:346\ttp:A:P",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_aggregate_split_hits_combines_collinear_nonoverlapping_segments(tmp_path: Path) -> None:
    paf = tmp_path / "split.paf"
    _write_split_paf(paf)

    raw_hits = genotyping.parse_paf(paf)["sample"]
    aggregated = genotyping.aggregate_split_hits(raw_hits)

    assert len(aggregated) == 1
    hit = aggregated[0]
    assert hit["alignment_segment_count"] == 2
    assert hit["query_aligned_bases"] == 5281
    assert hit["alignment_score"] == 1499
    assert hit["query_coverage"] == pytest.approx(5281 / 9453)
    assert hit["identity"] == pytest.approx((3955 + 851) / (4253 + 907))


def test_split_alignment_passes_query_coverage_after_aggregation(tmp_path: Path) -> None:
    paf = tmp_path / "split.paf"
    _write_split_paf(paf)
    hits = genotyping.parse_paf(paf)
    records = {"sample": {"description": "sample", "sequence": "N" * 9453}}
    args = SimpleNamespace(
        min_query_coverage=0.50,
        min_identity=0.75,
        close_hit_fraction=0.98,
    )

    row = genotyping.build_assignment_rows(records, hits, args)[0]

    assert row["assignment_status"] == "pass"
    assert row["assigned_genotype"] == "3a"
    assert row["qc_fail_reason"] == ""
    assert row["best_ref"] == "3a_D17763.1"
    assert row["coverage_ref"] == "3a_D17763.1"
    assert row["query_coverage"] == "0.558659"
    assert row["alignment_segment_count"] == 2
    assert row["non_n_bases"] == 0
    assert row["non_n_fraction"] == "0.000000"


def test_coverage_qc_uses_broadest_passing_reference_within_winning_genotype(tmp_path: Path) -> None:
    paf = tmp_path / "multiple_3a_refs.paf"
    paf.write_text(
        "\n".join(
            [
                # D17763 wins genotype assignment by score but covers only 55.9%.
                "sample\t9453\t501\t1409\t+\t3a_D17763.1\t9456\t501\t1409\t851\t907\t60\tAS:i:346\ttp:A:P",
                "sample\t9453\t4934\t9307\t+\t3a_D17763.1\t9456\t4935\t9308\t3955\t4253\t60\tAS:i:1153\ttp:A:P",
                # D28917 is a lower-scoring reference of the winning genotype but
                # provides the representative coverage evidence for QC.
                "sample\t9453\t507\t3610\t+\t3a_D28917.1\t9454\t507\t3611\t2833\t3105\t12\tAS:i:339\ttp:A:P",
                "sample\t9453\t4917\t9307\t+\t3a_D28917.1\t9454\t4918\t9308\t3928\t4272\t0\tAS:i:676\ttp:A:S",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    records = {
        "sample": {
            "description": "sample",
            "sequence": "A" * 7428 + "N" * 2025,
        }
    }
    args = SimpleNamespace(
        min_query_coverage=0.70,
        min_identity=0.75,
        close_hit_fraction=0.98,
    )

    row = genotyping.build_assignment_rows(records, genotyping.parse_paf(paf), args)[0]

    assert row["assignment_status"] == "pass"
    assert row["assigned_genotype"] == "3a"
    assert row["best_ref"] == "3a_D17763.1"
    assert row["best_alignment_score"] == 1499
    assert row["coverage_ref"] == "3a_D28917.1"
    assert row["query_coverage"] == "0.792658"
    assert row["identity"] == "0.916497"
    assert row["alignment_segment_count"] == 2
    assert row["non_n_bases"] == 7428
    assert row["non_n_fraction"] == "0.785782"


def test_split_hit_aggregation_rejects_discordant_target_gaps(tmp_path: Path) -> None:
    paf = tmp_path / "discordant.paf"
    paf.write_text(
        "\n".join(
            [
                "sample\t9000\t0\t3000\t+\t3a_D17763.1\t9500\t0\t3000\t2800\t3000\t60\tAS:i:900\ttp:A:P",
                # Query gap is 1000 nt but target gap is 4000 nt: not one collinear missing block.
                "sample\t9000\t4000\t7000\t+\t3a_D17763.1\t9500\t7000\t9500\t2300\t2500\t60\tAS:i:800\ttp:A:P",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    hit = genotyping.aggregate_split_hits(genotyping.parse_paf(paf)["sample"])[0]

    assert hit["alignment_segment_count"] == 1
    assert hit["query_aligned_bases"] == 3000


def test_split_hit_aggregation_rejects_query_overlap(tmp_path: Path) -> None:
    paf = tmp_path / "overlap.paf"
    paf.write_text(
        "\n".join(
            [
                "sample\t9000\t0\t4000\t+\t3a_D17763.1\t9500\t0\t4000\t3700\t4000\t60\tAS:i:900\ttp:A:P",
                "sample\t9000\t3500\t7000\t+\t3a_D17763.1\t9500\t3500\t7000\t3200\t3500\t60\tAS:i:800\ttp:A:P",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    hit = genotyping.aggregate_split_hits(genotyping.parse_paf(paf)["sample"])[0]

    assert hit["alignment_segment_count"] == 1
    assert hit["query_aligned_bases"] == 4000
