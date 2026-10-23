"""`hcv-cluster run` and the Nextflow workflow must produce the same tables.

Both call the same Python stages, so this guards the wiring (file hand-offs, naming,
channel logic), which is where the two could drift. Needs the real tools and Nextflow:
run it with `pixi run -e nf test`. It is skipped when Nextflow is not on PATH.
"""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "assets" / "test" / "example.fasta"

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None or shutil.which("tn93") is None,
    reason="needs nextflow and the clustering tools (pixi run -e nf test)",
)

# TN93 runs multithreaded, so row order within a table is not stable; compare as sets.
COMPARED = [
    "sample_audit.csv", "clusters.csv", "clusters.snp.csv", "links.csv", "links.snp.csv",
    "metadata.all_samples.csv", "metadata.surviving_samples.csv", "sequence_variants.csv",
    "cluster_labels.csv", "genotypes.csv",
    "by_genotype/1a/clusters.csv", "by_genotype/1a/snp_clusters.csv", "by_genotype/1a/tn93.csv",
    "by_genotype/1a/prep.qc.csv", "by_genotype/1a/variant_decisions.csv",
]


def _rows(path: Path) -> list[tuple]:
    with path.open(newline="") as handle:
        return sorted(tuple(row) for row in csv.reader(handle))


def _write_metadata(path: Path) -> Path:
    """Metadata in which clusterA_1 and clusterA_2 are one patient (sized-by-patient path)."""
    ids = dict.fromkeys(line[1:].strip() for line in FIXTURE.read_text().splitlines() if line.startswith(">"))
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "patient_index", "sample_index", "site"])
        for sample in ids:
            patient = "P_A12" if sample in ("clusterA_1", "clusterA_2") else f"P_{sample}"
            writer.writerow([sample, patient, f"S_{sample}", "X"])
    return path


def test_cli_and_workflow_agree(tmp_path: Path) -> None:
    metadata = _write_metadata(tmp_path / "metadata.csv")
    shared = ["--distance", "both", "--label-date", "2026-10-23"]
    cli_out, nf_out = tmp_path / "cli", tmp_path / "nf"
    subprocess.run(
        [sys.executable, "-m", "hcv_cluster", "run", "-i", str(FIXTURE), "-o", str(cli_out),
         "--registry", str(tmp_path / "reg_cli"), "--update-registry",
         "--metadata", str(metadata), "--patient-column", "patient_index", "--sample-column", "sample_index", *shared],
        check=True, capture_output=True, cwd=tmp_path, env={**os.environ, "PYTHONPATH": str(ROOT)},
    )
    subprocess.run(
        ["nextflow", "run", str(ROOT / "main.nf"), "--input", str(FIXTURE), "--outdir", str(nf_out),
         "--registry", str(tmp_path / "reg_nf"), "--update_registry", "true",
         "--metadata", str(metadata), "--patient_column", "patient_index", "--sample_column", "sample_index",
         "--distance", "both", "--label_date", "2026-10-23"],
        check=True, capture_output=True, cwd=tmp_path, env={**os.environ, "NXF_WORK": str(tmp_path / "work")},
    )
    for name in COMPARED:
        assert _rows(cli_out / name) == _rows(nf_out / name), name
    assert (cli_out / "surviving_samples.fasta").read_text() == (nf_out / "surviving_samples.fasta").read_text()
    assert _rows(tmp_path / "reg_cli" / "membership.csv") == _rows(tmp_path / "reg_nf" / "membership.csv")
