import json
import logging
from pathlib import Path

import pytest

import hcv_cluster_prep


def write_cache(
    tmp_path: Path,
    reference: str = "ACGTACGT",
    cached_rows: dict[str, str] | None = None,
    inputs: dict[str, str] | None = None,
    reference_accession: str = "REF.1",
) -> tuple[Path, Path]:
    """Build a cached alignment plus the sidecar prep-align would have written."""
    cached_rows = cached_rows if cached_rows is not None else {"old": "ACGT--GT"}
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(f">{reference_accession}\n{reference}\n")

    cache_path = tmp_path / "prep.aligned.fasta"
    lines = [f">{reference_accession}", reference]
    for header, row in cached_rows.items():
        lines.extend([f">{header}", row])
    cache_path.write_text("\n".join(lines) + "\n")

    # By default the recorded input is the degapped row; callers override this to model
    # a sample whose input differs from what --keeplength kept.
    inputs = inputs if inputs is not None else {h: r.replace("-", "") for h, r in cached_rows.items()}
    hcv_cluster_prep.write_alignment_cache_metadata(
        alignment_path=cache_path,
        records=[hcv_cluster_prep.FastaRecord(h, seq) for h, seq in inputs.items()],
        reference_accession=reference_accession,
        reference_ungapped=reference,
        alignment_length=len(reference),
    )
    return reference_path, cache_path


def reuse(records, reference_path, cache_path, accession="REF.1"):
    return hcv_cluster_prep.reuse_cached_alignment(
        records=records,
        cache_path=cache_path,
        reference_fasta=reference_path,
        reference_accession=accession,
        logger=logging.getLogger("test"),
    )


def test_reuse_cached_alignment_partitions_unchanged_and_new_records(tmp_path: Path) -> None:
    reference_path, cache_path = write_cache(tmp_path)

    reference_aligned, reused, to_align = reuse(
        [
            hcv_cluster_prep.FastaRecord("old", "ACGTGT"),
            hcv_cluster_prep.FastaRecord("new", "ACGTACGT"),
        ],
        reference_path,
        cache_path,
    )

    assert reference_aligned == "ACGTACGT"
    assert reused == {"old": "ACGT--GT"}
    assert to_align == [hcv_cluster_prep.FastaRecord("new", "ACGTACGT")]


def test_reuse_cached_alignment_realigns_a_changed_sample(tmp_path: Path) -> None:
    reference_path, cache_path = write_cache(tmp_path, cached_rows={"sample": "ACGT--GT"})

    _, reused, to_align = reuse(
        [hcv_cluster_prep.FastaRecord("sample", "ACGTATGT")], reference_path, cache_path
    )

    assert reused == {}
    assert to_align == [hcv_cluster_prep.FastaRecord("sample", "ACGTATGT")]


def test_reuse_survives_insertions_dropped_by_keeplength(tmp_path: Path) -> None:
    """The regression this cache design exists for.

    MAFFT runs with --keeplength, so residues inserted relative to the reference are
    absent from the aligned row. Degapping the row therefore cannot reproduce the
    input, and matching that way made every real sample a permanent cache miss.
    """
    reference_path, cache_path = write_cache(
        tmp_path,
        cached_rows={"sample": "ACGT--GT"},
        inputs={"sample": "ACGTTTTTGT"},  # insertion MAFFT discarded
    )

    _, reused, to_align = reuse(
        [hcv_cluster_prep.FastaRecord("sample", "ACGTTTTTGT")], reference_path, cache_path
    )

    assert reused == {"sample": "ACGT--GT"}
    assert to_align == []


def test_reuse_without_sidecar_realigns_everything(tmp_path: Path) -> None:
    reference_path, cache_path = write_cache(tmp_path)
    hcv_cluster_prep.alignment_cache_metadata_path(cache_path).unlink()

    _, reused, to_align = reuse(
        [hcv_cluster_prep.FastaRecord("old", "ACGTGT")], reference_path, cache_path
    )

    assert reused == {}
    assert [record.header for record in to_align] == ["old"]


def test_reuse_rejects_a_sidecar_built_against_another_reference(tmp_path: Path) -> None:
    reference_path, cache_path = write_cache(tmp_path)
    meta_path = hcv_cluster_prep.alignment_cache_metadata_path(cache_path)
    payload = json.loads(meta_path.read_text())
    payload["reference_fingerprint"] = hcv_cluster_prep.sequence_fingerprint("TTTTTTTT")
    meta_path.write_text(json.dumps(payload))

    with pytest.raises(hcv_cluster_prep.HcvPrepError, match="different reference sequence"):
        reuse([hcv_cluster_prep.FastaRecord("old", "ACGTGT")], reference_path, cache_path)


def _unsupported_sidecar_version(tmp_path: Path) -> tuple[Path, Path]:
    reference_path, cache_path = write_cache(tmp_path)
    meta_path = hcv_cluster_prep.alignment_cache_metadata_path(cache_path)
    payload = json.loads(meta_path.read_text())
    payload["version"] = hcv_cluster_prep.ALIGNMENT_CACHE_VERSION + 1
    meta_path.write_text(json.dumps(payload))
    return reference_path, cache_path


def _missing_cache_file(tmp_path: Path) -> tuple[Path, Path]:
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">REF.1\nACGTACGT\n")
    return reference_path, tmp_path / "nope.fasta"


def _cache_missing_reference(tmp_path: Path) -> tuple[Path, Path]:
    reference_path, cache_path = write_cache(tmp_path)
    cache_path.write_text(">other\nACGTACGT\n>old\nACGT--GT\n")
    return reference_path, cache_path


@pytest.mark.parametrize(
    ("make_bad_cache", "match"),
    [
        pytest.param(_unsupported_sidecar_version, "unsupported version", id="unsupported-version"),
        pytest.param(_missing_cache_file, "does not exist", id="missing-cache-file"),
        pytest.param(_cache_missing_reference, "does not contain reference", id="missing-reference"),
    ],
)
def test_reuse_rejects_a_bad_cache(tmp_path: Path, make_bad_cache, match: str) -> None:
    reference_path, cache_path = make_bad_cache(tmp_path)

    with pytest.raises(hcv_cluster_prep.HcvPrepError, match=match):
        reuse([hcv_cluster_prep.FastaRecord("old", "ACGTGT")], reference_path, cache_path)


def test_count_reusable_sequences_reports_hits_without_aligning(tmp_path: Path) -> None:
    _, cache_path = write_cache(
        tmp_path, cached_rows={"a": "ACGT--GT", "b": "ACGTAC--"}
    )
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">a\nACGTGT\n>b\nCHANGED\n>c\nACGTAC\n")

    assert hcv_cluster_prep.count_reusable_sequences(input_fasta, cache_path) == (1, 3)


@pytest.mark.parametrize("sidecar_only_missing", [False, True], ids=["no-cache", "no-sidecar"])
def test_count_reusable_sequences_without_usable_cache_is_zero(
    tmp_path: Path, sidecar_only_missing: bool
) -> None:
    if sidecar_only_missing:
        _, cache_path = write_cache(tmp_path)
        hcv_cluster_prep.alignment_cache_metadata_path(cache_path).unlink()
    else:
        cache_path = tmp_path / "absent.fasta"
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">old\nACGTGT\n")

    assert hcv_cluster_prep.count_reusable_sequences(input_fasta, cache_path) == (0, 1)
