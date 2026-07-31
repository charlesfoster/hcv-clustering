import logging
from pathlib import Path

import hcv_cluster_prep


def test_reuse_cached_alignment_partitions_unchanged_and_new_records(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fasta"
    cache = tmp_path / "prep.aligned.fasta"
    reference.write_text(">REF.1\nACGTACGT\n")
    cache.write_text(">REF.1\nACGTACGT\n>old\nACGT--GT\n")

    reference_aligned, reused, to_align = hcv_cluster_prep.reuse_cached_alignment(
        records=[
            hcv_cluster_prep.FastaRecord("old", "ACGTGT"),
            hcv_cluster_prep.FastaRecord("new", "ACGTACGT"),
        ],
        cache_path=cache,
        reference_fasta=reference,
        reference_accession="REF.1",
        logger=logging.getLogger("test"),
    )

    assert reference_aligned == "ACGTACGT"
    assert reused == {"old": "ACGT--GT"}
    assert to_align == [hcv_cluster_prep.FastaRecord("new", "ACGTACGT")]


def test_reuse_cached_alignment_realigns_a_changed_sample(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fasta"
    cache = tmp_path / "prep.aligned.fasta"
    reference.write_text(">REF.1\nACGTACGT\n")
    cache.write_text(">REF.1\nACGTACGT\n>sample\nACGT--GT\n")

    _, reused, to_align = hcv_cluster_prep.reuse_cached_alignment(
        records=[hcv_cluster_prep.FastaRecord("sample", "ACGTATGT")],
        cache_path=cache,
        reference_fasta=reference,
        reference_accession="REF.1",
        logger=logging.getLogger("test"),
    )

    assert reused == {}
    assert to_align == [hcv_cluster_prep.FastaRecord("sample", "ACGTATGT")]
