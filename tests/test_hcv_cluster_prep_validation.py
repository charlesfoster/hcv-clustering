import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from hcv_cluster import prep
from hcv_cluster.prep import FastaRecord

LOG = logging.getLogger("test")


def _random_dna(seed: int, length: int = 200) -> str:
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _validate(monkeypatch, tmp_path: Path, records, mode="kmer"):
    references = {"1a": _random_dna(1), "1b": _random_dna(2)}
    catalog = {g: SimpleNamespace(genotype=g) for g in references}

    def fake_download(spec, *_args, **_kwargs):
        path = tmp_path / f"{spec.genotype}.fasta"
        path.write_text(f">{spec.genotype}\n{references[spec.genotype]}\n")
        return SimpleNamespace(fasta_path=path)

    monkeypatch.setattr(prep, "download_reference", fake_download)
    return prep.validate_input_genotypes(
        records, "1a", mode, 15, tmp_path, catalog, False, None, LOG
    ), references


def test_kmer_conflicts_are_returned_per_sequence_not_raised(monkeypatch, tmp_path: Path) -> None:
    refs_1b = _random_dna(2)
    records = [
        FastaRecord("good", _random_dna(1)),
        FastaRecord("wrong_subtype", refs_1b),
        FastaRecord("too_short", "ACGT" * 5),
    ]
    rejected, _ = _validate(monkeypatch, tmp_path, records)
    assert set(rejected) == {"wrong_subtype", "too_short"}
    assert "closer to 1b" in rejected["wrong_subtype"]


def test_header_genotype_tag_mismatch_still_aborts(monkeypatch, tmp_path: Path) -> None:
    records = [FastaRecord("s1 genotype=1b", "ACGT" * 20)]
    with pytest.raises(prep.HcvPrepError, match="does not match"):
        _validate(monkeypatch, tmp_path, records)
