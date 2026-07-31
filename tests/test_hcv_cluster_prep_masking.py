import csv
from pathlib import Path

import hcv_cluster_prep


def test_mask_n_as_gap_replaces_n_only() -> None:
    assert hcv_cluster_prep.mask_n_as_gap("ACGTNNNNACGT") == "ACGT----ACGT"


def test_mask_n_as_gap_is_case_insensitive() -> None:
    assert hcv_cluster_prep.mask_n_as_gap("ACGTnnnnACGT") == "ACGT----ACGT"


def test_mask_n_as_gap_leaves_real_ambiguity_codes_untouched() -> None:
    assert hcv_cluster_prep.mask_n_as_gap("ACGTRYSWKMACGT") == "ACGTRYSWKMACGT"


def test_mask_n_as_gap_leaves_existing_gaps_untouched() -> None:
    assert hcv_cluster_prep.mask_n_as_gap("AC--NNGT") == "AC----GT"


def test_write_outputs_masks_n_in_clustering_fasta_but_not_aligned_fasta(tmp_path: Path) -> None:
    reference_aligned = "AAAAAAAAAA"
    query_alignment = {"sample1": "ACGTNNGTAC"}
    extracted = {"sample1": "ACGTNNGTAC"}
    selection = hcv_cluster_prep.RegionSelection(
        expression="core-e2-nohvr1",
        strategy="fixed",
        segments=(hcv_cluster_prep.RegionSegment("region", 1, 10, "test", ""),),
    )
    spec = hcv_cluster_prep.ReferenceSpec(genotype="1a", accession="TEST.1")
    out_prefix = tmp_path / "prep"

    import logging

    hcv_cluster_prep.write_outputs(
        records=[hcv_cluster_prep.FastaRecord("sample1", "ACGTNNGTAC")],
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
