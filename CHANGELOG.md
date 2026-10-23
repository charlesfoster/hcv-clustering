# Changelog

All notable changes to hcv-cluster. The format follows [Keep a Changelog](https://keepachangelog.com/);
see [docs/RELEASING.md](docs/RELEASING.md) for what each kind of version bump means. Anything that can
change which sequences cluster together for the same input is called out as **results-affecting**.

## [0.2.0] - unreleased

First versioned release, and the first that can run as a Nextflow pipeline.

### Added
- **Nextflow workflow** (`main.nf`, `workflows/`, `modules/`) with a reusable `HCV_CLUSTER` subworkflow, a
  Docker image built from `pixi.lock`, `-resume`, alignment reuse (`--reuse_alignments`) and network plots.
  Its tables are identical to `hcv-cluster run` (guarded by `tests/test_nextflow_parity.py`).
- **Sample audit** (`sample_audit.csv`): one row per input record with pass/fail and the reason at every
  stage (input, genotyping, subtype, genotype validation, region QC, variant selection, clustering), plus
  `metadata.all_samples.csv`, `metadata.surviving_samples.csv` and `surviving_samples.fasta`
  (headers suffixed `|genotype|cluster`).
- **Persistent cluster labels** (`--registry`, `--update-registry`): labels such as `HCV-2610-001` survive
  between runs, merges keep the oldest label, splits are recorded, singletons are never labelled. Several
  clustering methods can be tracked in one registry.
- **Cluster sizes by patient** (`--patient-column`); `cluster_sequences`, `cluster_samples` (`--sample-column`) and `cluster_patients` columns.
- **Mixed-infection handling**: repeated FASTA IDs are resolved into variants instead of failing.
- **`run_info.json`** with every run: pipeline version and commit, tool versions, reference data, parameters.
- `hcv-stage` commands (`resolve-ids`, `plan-genotypes`, `select-variants`, `cluster`, `merge`, `label`,
  `audit`, `plot`, `provenance`) so any workflow engine can run the pipeline step by step.
- Region-specific default thresholds with a documented evidence basis (`docs/threshold_rationale.md`),
  HVR1-masked regions, and an absolute SNP-count linking option.
- Reference records for all 84 catalogued subtypes ship with the package; no network access is needed.

### Changed
- **Results-affecting:** default threshold is now chosen by region and distance metric (0.03 for the default
  `e1-e2`), not by genotype; default region is `e1-e2`; TN93 ambiguity handling defaults to `average`;
  minimap2 preset `asm20`.
- A sequence whose k-mer profile conflicts with its assigned subtype is excluded (and audited) instead of
  aborting the whole run.
- Same-sample variants of one genotype compete only if they passed region QC.
- Tools updated: TN93 1.1, minimap2 2.31, biopython 1.88 (outputs verified identical to the previous
  versions on real and synthetic data).
- Code reorganised into the `hcv_cluster` package (`python -m hcv_cluster`); the old flat scripts and their
  import names are gone.

### Fixed
- A subtype without a reference genome no longer stops the run; it is skipped and audited.
- Reference cache files are no longer rewritten on every run, so a read-only reference directory works.
