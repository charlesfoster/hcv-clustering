import csv
from pathlib import Path

import pytest

from hcv_cluster import prep


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        pytest.param("ACGTNNNNACGT", "ACGT----ACGT", id="replaces-n"),
        pytest.param("ACGTnnnnACGT", "ACGT----ACGT", id="case-insensitive"),
        pytest.param("ACGTRYSWKMACGT", "ACGTRYSWKMACGT", id="keeps-ambiguity-codes"),
        pytest.param("AC--NNGT", "AC----GT", id="keeps-existing-gaps"),
    ],
)
def test_mask_n_as_gap(sequence: str, expected: str) -> None:
    assert prep.mask_n_as_gap(sequence) == expected


def test_write_outputs_masks_n_in_clustering_fasta_but_not_aligned_fasta(tmp_path: Path) -> None:
    reference_aligned = "AAAAAAAAAA"
    query_alignment = {"sample1": "ACGTNNGTAC"}
    extracted = {"sample1": "ACGTNNGTAC"}
    selection = prep.RegionSelection(
        expression="core-e2-nohvr1",
        strategy="fixed",
        segments=(prep.RegionSegment("region", 1, 10, "test", ""),),
    )
    spec = prep.ReferenceSpec(genotype="1a", accession="TEST.1")
    out_prefix = tmp_path / "prep"

    import logging

    prep.write_outputs(
        records=[prep.FastaRecord("sample1", "ACGTNNGTAC")],
        query_alignment=query_alignment,
        reference_aligned=reference_aligned,
        extracted=extracted,
        selection=selection,
        spec=spec,
        out_prefix=out_prefix,
        min_coverage=0.5,
        logger=logging.getLogger("test"),
    )

    clustering_fasta = (tmp_path / "prep.clustering.fasta").read_text()
    assert "ACGT--GTAC" in clustering_fasta
    assert "N" not in clustering_fasta

    aligned_fasta = (tmp_path / "prep.aligned.fasta").read_text()
    assert "ACGTNNGTAC" in aligned_fasta

    with (tmp_path / "prep.qc.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["n_bases"] == "2"
    assert rows[0]["passed_qc"] == "true"
