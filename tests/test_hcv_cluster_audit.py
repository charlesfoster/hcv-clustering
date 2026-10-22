import csv
from pathlib import Path

import hcv_cluster_audit as audit_mod
import hcv_cluster_prep
import hcv_cluster_variants
from hcv_cluster_prep import FastaRecord


def _geno(seq_id: str, genotype: str = "", status: str = "pass", reason: str = "") -> dict[str, str]:
    return {
        "sequence_id": seq_id,
        "assigned_genotype": genotype,
        "assignment_status": status,
        "qc_fail_reason": reason,
        "query_coverage": "0.9",
        "identity": "0.8",
    }


def _qc(seq_id: str, passed: bool, coverage: float = 0.9) -> dict[str, str]:
    return {
        "sample_id": seq_id,
        "passed_qc": "true" if passed else "false",
        "coverage_fraction": str(coverage),
        "selected_region": "e1-e2",
        "region_length": "100",
        "non_gap_bases": "90",
        "n_bases": "0",
        "gap_bases": "10",
    }


def _build(records, genotype_rows, **kwargs):
    plan = hcv_cluster_variants.resolve_duplicate_ids(records)
    return plan, audit_mod.build_audit(
        records,
        plan,
        genotype_rows,
        min_query_coverage=0.5,
        min_identity=0.75,
        min_coverage=0.7,
        **kwargs,
    )


def test_each_stage_failure_is_attributed_and_later_stages_not_run() -> None:
    records = [FastaRecord(name, "ACGT") for name in ("ok", "lowid", "norefs", "lowcov", "skip", "kmer")]
    genotype_rows = [
        _geno("kmer", "1a"),
        _geno("ok", "1a"),
        _geno("lowid", status="qc_fail", reason="low_identity"),
        _geno("norefs", "2l"),
        _geno("lowcov", "1a"),
        _geno("skip", "3a"),
    ]
    _plan, audits = _build(
        records,
        genotype_rows,
        excluded_subtypes={
            "2l": "No reference genome configured for subtype 2l",
            "3a": "Subtype 3a excluded by --genotype (1a)",
        },
        qc_by_genotype={"1a": {"ok": _qc("ok", True), "lowcov": _qc("lowcov", False, 0.4)}},
        validation_failures={"1a": {"kmer": "15-mer match is closer to 1b than to assigned 1a"}},
        cluster_rows=[{"sample_id": "ok", "cluster_id": "1a_C0001", "cluster_size": "3"}],
    )
    by_id = {a.sample_id: a.as_row() for a in audits}

    assert by_id["ok"]["overall_status"] == "pass"
    assert by_id["ok"]["cluster_id"] == "1a_C0001"
    assert by_id["lowid"]["failed_stage"] == "genotyping"
    assert "identity 0.8 < 0.75" in by_id["lowid"]["overall_reason"]
    assert by_id["lowid"]["subtype_status"] == "not_run"
    assert by_id["norefs"]["failed_stage"] == "subtype"
    assert "No reference genome" in by_id["norefs"]["overall_reason"]
    assert by_id["lowcov"]["failed_stage"] == "region_qc"
    assert "coverage 0.400 < 0.7" in by_id["lowcov"]["overall_reason"]
    assert by_id["lowcov"]["clustering_status"] == "not_run"
    assert by_id["kmer"]["failed_stage"] == "genotype_validation"
    assert by_id["kmer"]["region_qc_status"] == "not_run"
    assert by_id["skip"]["failed_stage"] == "subtype"
    assert "--genotype" in by_id["skip"]["overall_reason"]


def test_no_alignment_and_missing_genotype_row_have_readable_reasons() -> None:
    records = [FastaRecord("a", "ACGT"), FastaRecord("b", "ACGT")]
    _plan, audits = _build(records, [_geno("a", status="qc_fail", reason="no_alignment")])
    reasons = {a.sample_id: a.as_row()["overall_reason"] for a in audits}
    assert "no alignment" in reasons["a"]
    assert "no genotype assignment" in reasons["b"]


def test_repeats_and_variants_get_one_audit_row_per_input_record(tmp_path: Path) -> None:
    records = [
        FastaRecord("s1", "AAAA"),
        FastaRecord("s1", "AAAA"),  # exact repeat
        FastaRecord("s2", "CCCC"),
        FastaRecord("s2", "CCCG"),  # distinct variant, same genotype
    ]
    plan = hcv_cluster_variants.resolve_duplicate_ids(records)
    genotype_rows = [_geno("s1", "1a"), _geno("s2__v1", "1a"), _geno("s2__v2", "1a")]
    hcv_cluster_variants.record_genotypes(plan, genotype_rows)
    qc_rows = {sid: _qc(sid, True) for sid in ("s1", "s2__v1", "s2__v2")}
    qc_rows["s2__v2"]["coverage_fraction"] = "0.95"
    qc_path = tmp_path / "prep.qc.csv"
    with qc_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qc_rows["s1"]))
        writer.writeheader()
        writer.writerows(qc_rows.values())
    hcv_cluster_variants.select_same_genotype(plan, "1a", qc_path, "region-coverage")

    audits = audit_mod.build_audit(
        records,
        plan,
        genotype_rows,
        min_query_coverage=0.5,
        min_identity=0.75,
        min_coverage=0.7,
        qc_by_genotype={"1a": audit_mod.read_qc_rows(qc_path)},
        cluster_rows=[
            {"sample_id": "s1", "cluster_id": "1a_C0001", "cluster_size": "2"},
            {"sample_id": "s2__v2", "cluster_id": "1a_C0001", "cluster_size": "2"},
        ],
    )
    rows = [a.as_row() for a in audits]
    assert [r["input_record"] for r in rows] == [1, 2, 3, 4]
    assert rows[1]["failed_stage"] == "input"
    assert rows[2]["failed_stage"] == "variant_selection"
    assert rows[2]["base_sample_id"] == "s2" and rows[2]["mixed_infection"] == "true"
    assert rows[3]["overall_status"] == "pass"
    assert audit_mod.summarize(audits) == (
        "2/4 input records survived; excluded at input (1), variant_selection (1)"
    )


def test_outputs_metadata_tables_and_labelled_fasta(tmp_path: Path) -> None:
    records = [FastaRecord("good extra words", "ACGT"), FastaRecord("bad", "TTTT")]
    _plan, audits = _build(
        records,
        [_geno("good", "1a"), _geno("bad", status="qc_fail", reason="low_identity")],
        qc_by_genotype={"1a": {"good": _qc("good", True)}},
        cluster_rows=[{"sample_id": "good", "cluster_id": "1a_C0001", "cluster_size": "1"}],
    )
    metadata = [{"sample_id": "good", "site": "A"}, {"sample_id": "bad", "site": "B"}]
    all_path, surviving_path = audit_mod.write_metadata_tables(
        audits, metadata, tmp_path / "all.csv", tmp_path / "surviving.csv"
    )
    read = lambda p: list(csv.DictReader(p.open()))  # noqa: E731
    assert [r["sample_id"] for r in read(all_path)] == ["good", "bad"]
    assert [r["sample_id"] for r in read(surviving_path)] == ["good"]
    assert read(surviving_path)[0]["site"] == "A"
    assert read(all_path)[1]["failed_stage"] == "genotyping"

    fasta = audit_mod.write_surviving_fasta(audits, records, tmp_path / "survivors.fasta")
    parsed = hcv_cluster_prep.read_fasta_records(fasta)
    assert [r.header for r in parsed] == ["good|1a|1a_C0001 extra words"]


def test_metadata_tables_are_written_without_user_metadata(tmp_path: Path) -> None:
    _plan, audits = _build([FastaRecord("only", "ACGT")], [_geno("only", status="qc_fail", reason="low_identity")])
    all_path, surviving_path = audit_mod.write_metadata_tables(
        audits, [], tmp_path / "all.csv", tmp_path / "surviving.csv"
    )
    assert len(list(csv.DictReader(all_path.open()))) == 1
    assert list(csv.DictReader(surviving_path.open())) == []
