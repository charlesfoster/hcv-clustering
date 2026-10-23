"""The stage commands hand files to one another; run the chain the Nextflow workflow runs.

Alignment and genotyping need external tools, so their outputs are written by hand; every
other stage runs for real (SNP distances need no binary).
"""
import csv
from pathlib import Path

from hcv_cluster import stages

BASE = "ACGTACGTAC" * 12  # 120 nt


def _mutate(seq: str, *positions: int) -> str:
    chars = list(seq)
    for p in positions:
        chars[p] = "A" if chars[p] != "A" else "C"
    return "".join(chars)


def _fasta(path: Path, records: dict[str, str]) -> Path:
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records.items()))
    return path


def _rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open()))


def _run(*argv: str) -> None:
    assert stages.main(list(argv)) == 0


def test_stage_chain_produces_labelled_clusters_and_a_complete_audit(tmp_path: Path) -> None:
    # a/b/c are near-identical (one cluster); d is distant; "a" appears twice (mixed infection);
    # "dup" is an exact repeat; "junk" fails genotyping.
    seqs = {"b": _mutate(BASE, 1), "c": _mutate(BASE, 2), "d": _mutate(BASE, *range(0, 60, 2))}
    records = [("a", BASE), ("a", _mutate(BASE, 50)), ("b", seqs["b"]), ("c", seqs["c"]),
               ("d", seqs["d"]), ("dup", BASE), ("dup", BASE), ("junk", "ACGT" * 10)]
    inp = tmp_path / "input.fasta"
    inp.write_text("".join(f">{h}\n{s}\n" for h, s in records))

    _run("resolve-ids", "--input", str(inp), "--out-fasta", str(tmp_path / "resolved.fasta"),
         "--out-variants", str(tmp_path / "variants.csv"))
    resolved = [line[1:].strip() for line in (tmp_path / "resolved.fasta").read_text().splitlines() if line.startswith(">")]
    assert resolved == ["a__v1", "a__v2", "b", "c", "d", "dup", "junk"]

    # Stand-ins for the genotype and align stages (they wrap minimap2 and MAFFT).
    with (tmp_path / "genotypes.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence_id", "assigned_genotype", "assignment_status", "qc_fail_reason"])
        writer.writeheader()
        for sid in ("a__v1", "a__v2", "b", "c", "d", "dup"):
            writer.writerow({"sequence_id": sid, "assigned_genotype": "1a", "assignment_status": "pass", "qc_fail_reason": ""})
        writer.writerow({"sequence_id": "junk", "assigned_genotype": "", "assignment_status": "qc_fail", "qc_fail_reason": "no_alignment"})
    aligned = {"a__v1": BASE, "a__v2": _mutate(BASE, 50), **seqs, "dup": BASE}
    clustering = _fasta(tmp_path / "prep.clustering.fasta", aligned)
    with (tmp_path / "prep.qc.csv").open("w", newline="") as handle:
        fields = ["sample_id", "passed_qc", "coverage_fraction", "selected_region", "region_length", "non_gap_bases", "n_bases", "gap_bases"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for sid in aligned:
            writer.writerow({"sample_id": sid, "passed_qc": "true", "coverage_fraction": "0.95" if sid == "a__v1" else "0.9",
                             "selected_region": "e1-e2", "region_length": "120", "non_gap_bases": "114", "n_bases": "0", "gap_bases": "6"})

    _run("select-variants", "--variants", str(tmp_path / "variants.csv"), "--genotypes", str(tmp_path / "genotypes.csv"),
         "--genotype", "1a", "--qc", str(tmp_path / "prep.qc.csv"), "--clustering-fasta", str(clustering),
         "--out-fasta", str(tmp_path / "selected.fasta"), "--out-decisions", str(tmp_path / "decisions.csv"))
    kept = [line[1:] for line in (tmp_path / "selected.fasta").read_text().splitlines() if line.startswith(">")]
    assert "a__v1" in kept and "a__v2" not in kept  # better-covered variant wins

    _run("cluster", "--method", "snp", "--input", str(tmp_path / "selected.fasta"), "--region", "e1-e2",
         "--out-distances", str(tmp_path / "snp.csv"), "--out-links", str(tmp_path / "links.csv"),
         "--out-clusters", str(tmp_path / "clusters_1a.csv"), "--threshold", "0.05")
    _run("merge", "--clusters", f"1a={tmp_path / 'clusters_1a.csv'}", "--links", f"1a={tmp_path / 'links.csv'}",
         "--variants", str(tmp_path / "variants.csv"),
         "--out-clusters", str(tmp_path / "clusters.csv"), "--out-links", str(tmp_path / "merged_links.csv"))

    registry = tmp_path / "registry"
    _run("label", "--clusters", str(tmp_path / "clusters.csv"), "--registry", str(registry), "--update-registry",
         "--registry-out", str(tmp_path / "registry_out"), "--metric", "snp", "--region", "e1-e2", "--threshold", "0.05",
         "--label-date", "2026-10-23", "--out-clusters", str(tmp_path / "clusters.labelled.csv"),
         "--out-labels", str(tmp_path / "cluster_labels.csv"), "--out-changes", str(tmp_path / "label_changes.csv"))
    assert not registry.exists()  # the staged input registry is never modified
    labelled = {r["sample_id"]: r["cluster_label"] for r in _rows(tmp_path / "clusters.labelled.csv")}
    assert labelled["a__v1"] == labelled["b"] == labelled["c"] == "HCV-2610-001"
    assert labelled["d"] == ""  # singleton
    assert (tmp_path / "registry_out" / "membership.csv").exists()

    _run("audit", "--input", str(inp), "--resolved-fasta", str(tmp_path / "resolved.fasta"),
         "--variants", str(tmp_path / "variants.csv"), "--genotypes", str(tmp_path / "genotypes.csv"),
         "--qc", f"1a={tmp_path / 'prep.qc.csv'}", "--decisions", f"1a={tmp_path / 'decisions.csv'}",
         "--clusters", str(tmp_path / "clusters.labelled.csv"), "--outdir", str(tmp_path / "out"))
    audit = {(r["sample_id"], r["input_record"]): r for r in _rows(tmp_path / "out" / "sample_audit.csv")}
    outcome = {sid: (r["overall_status"], r["failed_stage"]) for (sid, _), r in audit.items()}
    assert len(audit) == len(records)  # one row per input record
    assert outcome["a__v2"] == ("fail", "variant_selection")
    assert outcome["junk"] == ("fail", "genotyping")
    assert outcome["b"] == ("pass", "")
    assert sum(1 for r in audit.values() if r["failed_stage"] == "input") == 1  # the repeat
    survivors = (tmp_path / "out" / "surviving_samples.fasta").read_text()
    assert ">b|1a|HCV-2610-001" in survivors and ">d|1a|singleton" in survivors
