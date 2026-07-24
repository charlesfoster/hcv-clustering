from pathlib import Path

import assign_hcv_genotypes_from_fasta as genotyping


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
