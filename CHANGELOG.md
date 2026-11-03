# Changelog

All notable changes to hcv-cluster. The format follows [Keep a Changelog](https://keepachangelog.com/);
see [docs/RELEASING.md](docs/RELEASING.md) for what each kind of version bump means. Anything that can
change which sequences cluster together for the same input is called out as **results-affecting**.

## [0.2.0] - unreleased

First versioned release, and the first that can run as a Nextflow pipeline.

### Added
- **IAP model** (`--iap-model`, **Run IAP model** in the GUI, `--iap_model` in Nextflow, `hcv-stage iap-pairs`):
  `iap_pairs.csv` with the raw SNP count, years between collection dates, SNPs per year and relatedness index
  (a negative-binomial model, ported from the R notebook) for every dated pair within a subtype, plus a `subtype` column and a `note` column explaining empty values (`too_far_apart`, `same_day`). One table per subtype (`by_genotype/<subtype>/`) and a merged one. Needs a metadata collection date column (`--date-column`, default `sample_collection_date`; dates as `d/m/yyyy`, `yyyy-mm-dd` or `d-Mon-yy`), and
  aborts for any region but the E1-E2 (HVR1 included) the model is calibrated for.
- **IAP ambiguity codes and flags:** a site where either sequence carries an IUPAC ambiguity code is counted as the chance the two
  differ (`--iap-ambiguity average`, the default; `skip` ignores such sites), the closest stand-in for the majority-base calls the
  model was fitted on. `iap_pairs.csv` gains `compared_sites`, `ambiguous_sites` and `flags` (`low_coverage`: the pair compares
  under 95% of the most complete sequence's sites; `noisy`: a sequence has 10 or more ambiguity codes), carried into the IAP link tables.
- **IAP clusters** (`--iap-alpha`, default 0.05): a second clustering from the relatedness index alone
  (`clusters.iap.csv`, `links.iap.csv`), the main links annotated with the index (`links.iap_information.csv`),
  and `clusters.comparison.csv` (identical / contained / split / absent, in both directions). The IAP clusters
  have their own label track (`iap`) and prefix (`--iap-label-prefix`, default `HCV_IAP`) in the same registry and their own network plots (`iap_network.*`); the GUI
  viewer switches between main and IAP clusters.
- **One-page progress report** (`hcv-report RUN_DIR`): a standalone command, run when a report is wanted (or at the
  end of a run with `--generate-report`; off unless asked for). Writes `report.html`, an A4 summary of a finished run: the sequence funnel, cluster
  statistics, age, sex, Indigenous status, injecting and diagnosis breakdowns from the metadata, collection months,
  a by-location table, a map of the five largest clusters and registry changes since the last report.
  `--since PREV_RUN_DIR` flags new and grown clusters and new sequences relative to the run behind the last report.
  The time charts show the last 36 months by default (`--chart-months`, `--chart-from`); the largest-clusters
  table shows each cluster's first and last sample and whether it is active, and `--active-only` limits the map,
  timelines and table to active clusters. Counts only: no sample identifiers.
- **Report on demand after a run** (`--generate-report`, a tick box in the GUI): writes `report.html` and, where Chrome
  or Chromium is installed, `report.pdf`, with download buttons in the GUI; the HTML has its own Save as PDF button,
  so any browser can make the PDF. A report that fails never fails the run.
  The H2Seq logo and title are the defaults. A run that saves its labels also saves what the report showed in the registry
  (`last_report.json`), so the next report shows growth since it with no `--since` and no old folder to remember.
- **Lab settings** (`~/.hcv-cluster/config.toml`): `registry`, `generate_report` and `threads` apply to the
  command line and the GUI. `hcv-init` writes the file by asking questions. The GUI shows the registry location, a tick box to save the run's
  labels and one to run without a registry, and blocks a run until the registry location is right.
- **Registry safety:** a registry named in the settings must already exist (or be declared new with `--new-registry`);
  a folder that is not a registry, or whose files are damaged, is refused; saving is refused when the input holds fewer
  than half of the samples the registry labels (`--allow-registry-mismatch`); every save is preceded by a backup
  in `backups/<run id>/`. `--update-registry` can now be turned off with `--no-update-registry`.
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
  `audit`, `plot`, `iap-check`, `iap-pairs`, `iap-merge`, `iap-support`, `provenance`) so any workflow engine can run the pipeline step by step.
- Region-specific default thresholds with a documented evidence basis (`docs/threshold_rationale.md`),
  HVR1-masked regions, and an absolute SNP-count linking option.
- Reference records for all 84 catalogued subtypes ship with the package; no network access is needed.

### Changed
- The pixi tasks are now just `hcv-cluster`, `hcv-report`, `hcv-init`, `hcv-cluster-gui` and `test`. `hcv-stage`, `hcv-workflow`
  and `hcv-nextflow` are gone: run `pixi run python -m hcv_cluster.stages` or `.workflow`, and
  `pixi run -e nf nextflow run main.nf`. The `hcv-stage` and `hcv-workflow` commands of a pip install are unchanged.
- **Results-affecting:** default threshold is now chosen by region and distance metric (0.045 for the default
  `e1-e2`), not by genotype; default region is `e1-e2`; TN93 ambiguity handling defaults to `average`;
  minimap2 preset `asm20`.
- A sequence whose k-mer profile conflicts with its assigned subtype is excluded (and audited) instead of
  aborting the whole run.
- Same-sample variants of one genotype compete only if they passed region QC.
- Tools updated: TN93 1.1, minimap2 2.31, biopython 1.88 (outputs verified identical to the previous
  versions on real and synthetic data).
- Added `numpy` and `scipy` (the IAP relatedness index). The container lock also moves `libgcc`/`libgfortran`
  from 15 to 16, which scipy needs; no scientific tool version changed.
- Code reorganised into the `hcv_cluster` package (`python -m hcv_cluster`); the old flat scripts and their
  import names are gone.

### Fixed
- The layout of a small connected component in a network plot depended on Python's per-process hash seed, so
  the same data could be laid out differently in two runs (nodes swapped or mirrored). It is now deterministic.
- With `--plot_network`, the Nextflow workflow replaced each `by_genotype/<subtype>` folder with one holding
  only the plots; it now publishes just the plot files.
- A subtype without a reference genome no longer stops the run; it is skipped and audited.
- Reference cache files are no longer rewritten on every run, so a read-only reference directory works.
