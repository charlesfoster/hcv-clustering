# HCV Clustering Workflow

End-to-end HCV genomic clustering for transmission surveillance. Given a multi-FASTA of consensus sequences, the workflow assigns genotypes, extracts a homologous aligned region per genotype, computes pairwise distances, and defines clusters of likely recent transmission using a threshold-based connected-component approach.

---

## Requirements

[Pixi](https://pixi.sh) is required to manage dependencies. Install it with:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Or via Homebrew on macOS:

```bash
brew install pixi
```

Restart your terminal after installation, then verify with `pixi --version`.

### Supported platforms

Linux (x86-64 and arm64) and macOS (Apple silicon and Intel). The alignment,
genotyping and distance tools (MAFFT, minimap2, TN93) have no native Windows
builds in conda-forge or bioconda, so Windows is not supported directly — run
the workflow under [WSL2](https://learn.microsoft.com/windows/wsl/install),
where it resolves as ordinary Linux. The Streamlit GUI is reachable from a
Windows browser at `localhost` when started inside WSL2.

## Setup

Get the code:

```bash
git clone https://github.com/charlesfoster/hcv-clustering.git
cd hcv-clustering
```

Then install all dependencies into an isolated environment:

```bash
pixi install
```

Every command below is run from inside the `hcv-clustering` directory.

---

## Quick start

### Step 1 — Reference sequences (nothing to do)

The GenBank reference records used for alignment coordinates ship with the repository (`hcv_cluster/data/genbank`), so the first run needs no network access. `pixi run hcv-cluster refs` is only needed if you add a subtype with `--reference-map` and want to download it ahead of time, or force a re-download with `--force-download`.

### Step 2 — Run the clustering pipeline

```bash
pixi run hcv-cluster run -i samples.fasta
```

Genotypes are detected automatically. The default analysis uses full `e1-e2` at a provisional, conservative threshold of 0.03. The same 0.03 operational default is used for `e1-e2-nohvr1`; established alternatives retain their region-specific values (`core-e2` 0.045, `core-e2-nohvr1` 0.03, and `ns5b` 0.015). See [`docs/threshold_rationale.md`](docs/threshold_rationale.md) for the evidence and limitations. Results are written to `results/` by default.

For recurring runs against a growing input FASTA, write to a new output
directory while reusing reference-anchored alignments from an earlier results
directory:

```bash
pixi run hcv-cluster run -i samples.fasta -o results_new \
  --reuse-alignments results_previous
```

Unchanged samples are read from each genotype's existing
`results_previous/by_genotype/<genotype>/prep.aligned.fasta`; only new or
changed sequences are sent to MAFFT. Region extraction and clustering outputs
are rebuilt in `results_new` so the selected settings and new cluster
connections are applied consistently.

A sample counts as unchanged when its input sequence matches the fingerprint
recorded in the sidecar `prep.aligned.fasta.meta.json` written alongside each
alignment. The run reports what the cache achieved, so a cache that is not
helping is visible rather than silent:

```
--> Reusing 812/840 cached alignments from results_previous/...; 28 to align
```

An alignment produced before the sidecar existed cannot be matched against, and
every sequence is realigned (with a warning). Re-run once to write one.

### Reporting differences as whole SNPs

Every SNP-distance run reports the absolute number of differing bases alongside
the p-distance. `snp.csv` and the SNP `links.csv` both carry `snp_count` and
`comparable_sites`, and each subtype prints a plain-English summary:

```
--> SNP   (threshold 0.03): 4 cluster(s), 11 singleton(s), largest: 6
          linked pairs differ by 0-62 SNPs (median 9) over 2098-2157 comparable sites
```

To cluster on that number directly — "link anything within 30 SNPs" — use:

```bash
pixi run hcv-cluster run -i samples.fasta --distance snp --snp-count-threshold 30
```

This is offered because people reason in whole SNPs, but it is **not** the
default and is not evidence-based. `comparable_sites` differs between pairs
(N and gap positions are skipped), so the same count means different divergence
for different pairs, and no HCV clustering threshold in the literature is
defined this way. See [`docs/threshold_rationale.md`](docs/threshold_rationale.md).

Common options:

```bash
pixi run hcv-cluster run \
  -i samples.fasta \
  -o my_results \
  -t 0.03 \
  -d both \
  -r e1-e2 \
  -T 4
```

| Short | Long | Default | Description |
|-------|------|---------|-------------|
| `-i` | `--input` | required | Input multi-FASTA |
| `-o` | `--outdir` | `results` | Output directory |
| `-t` | `--threshold` | region-dependent | Maximum distance for cluster linking |
| `-d` | `--distance` | `tn93` | `tn93`, `snp`, or `both` |
| `-r` | `--region` | `e1-e2` | Genomic region expression (see below) |
| `-T` | `--threads` | `1` | MAFFT alignment threads |
| `-v` | `--verbose` | off | Debug logging |

For all advanced tuning options (minimap2/MAFFT/TN93 settings, QC thresholds, etc.):

```bash
pixi run hcv-cluster run --help-advanced
```

To also save cluster network plots (per genotype, and per metric if `--distance both`):

```bash
pixi run hcv-cluster run -i samples.fasta --plot-network png --plot-hide-singletons
```

`--plot-network` accepts `png` (static image), `html` (interactive, opens in a browser and works offline), `svg` (editable vector image), `both` (PNG + HTML, retained for backward compatibility), or `all` (PNG + HTML + SVG). Per-genotype files are written alongside the other outputs as `by_genotype/<genotype>/network.*` (and `snp_network.*` for the SNP metric under `--distance both`). `--plot-hide-singletons` restricts plots to sequences that fall in a multi-member cluster.

Plot scope is independent of the file format:

```bash
# Keep the default per-genotype plots and also create one combined plot.
pixi run hcv-cluster run -i samples.fasta \
  --plot-network all \
  --plot-scope both \
  --plot-combined-layout by-genotype
```

`--plot-scope` accepts `per-genotype` (the default), `combined`, or `both`. A combined plot is a **presentation-only union** of the genotype networks: distances and clusters are still calculated separately within each genotype, and the workflow never creates cross-genotype links. `--plot-combined-layout packed` packs disconnected clusters together; `by-genotype` places those packed clusters in genotype-specific regions.

### Repeated FASTA IDs and mixed infections

The same sample ID may appear on several FASTA records, for example when a mixed infection yields more than one consensus sequence:

- **Identical sequences** under one ID are exact repeats; only the first is kept.
- **Distinct sequences** under one ID are renamed `ID__v1`, `ID__v2`, ... in input order, so every step keeps a unique key.
- Variants assigned **different genotypes** are all kept, each in its own genotype network.
- Variants sharing a **genotype** are reduced to one per sample, so a person never appears twice in one network. `--duplicate-selection region-coverage` (the default) keeps the variant with the best coverage of the clustering region, falling back to the most complete sequence (most unambiguous A/C/G/T bases) when coverage ties; `--duplicate-selection completeness` keeps the most complete sequence.

Cluster tables carry `base_sample_id` (the original ID) and `mixed_infection` (`true` when a sample had more than one distinct sequence). Metadata is joined on `base_sample_id`, so a single metadata row per person applies to every variant, and both fields can be used for plot encodings or hover. `sequence_variants.csv` records each repeated record's new name, genotype, completeness, region coverage, and whether it was kept or why it was dropped.

### Sample metadata and plot appearance

Supply an optional UTF-8 CSV keyed by the exact FASTA sample identifier:

```csv
sample_id,age,location,indigenous_status,injecting_status,subtype
sample_001,27,Prison A,yes,current,1a
sample_002,44,Prison B,no,former,1a
sample_003,68,Prison A,yes,never,3a
```

```bash
pixi run hcv-cluster run -i samples.fasta \
  --metadata sample_metadata.csv \
  --plot-network all \
  --plot-scope both \
  --plot-color-by location \
  --plot-symbol-by age_range \
  --plot-outline-by indigenous_status \
  --plot-center-by injecting_status \
  --plot-node-spacing expanded \
  --plot-hover-field subtype
```

`sample_id` is required, must be unique, and must match FASTA identifiers. All other columns are retained as arbitrary plotting metadata, so future fields do not require a schema change. Blank values and plotted samples without a metadata row appear as `(missing)`; rows that do not appear in the plotted results are reported. Duplicate/blank identifiers, malformed ages, invalid headers, and reserved clustering-output column names are rejected with an actionable error. Metadata changes presentation only and cannot change distances, links, or cluster membership. A normalized copy is saved as `metadata.csv` in the result directory for later GUI rendering.

When an `age` column is supplied, the workflow validates it as a non-negative whole number and derives the ordered `age_range` categories `0–30`, `31–60`, and `61+`. You may instead supply `age_range` directly; both ASCII forms (`0-30`, `31-60`) and en-dash forms are accepted and normalized. If both columns are present, they must agree.

Colour, shape, outline, and a small centre mark are independent channels and can be used simultaneously through `--plot-color-by`, `--plot-symbol-by`, `--plot-outline-by`, and `--plot-center-by`. The recommended epidemiological mapping is location by colour, `age_range` by shape, indigenous status by outline, and injecting status by centre mark. These channels retain a common node footprint and are easier to distinguish than small changes in node area. Node size remains available through the advanced `--plot-size-by` channel for numeric or intrinsically ordered data; for another categorical field, repeat `--plot-size-order VALUE` from smallest to largest. Shape, outline, and centre marks work best for small numbers of categories. Repeat `--plot-hover-field FIELD` to select several hover fields. Sample ID and cluster details are always included in interactive hover text, while static PNG/SVG files deliberately do not draw sample-name labels.

For a clearer publication view, repeat the same fixed network in up to four colour panels:

```bash
pixi run hcv-cluster run -i samples.fasta \
  --metadata sample_metadata.csv \
  --plot-network svg \
  --plot-scope combined \
  --plot-small-multiple-field location \
  --plot-small-multiple-field age_range \
  --plot-small-multiple-field injecting_status \
  --plot-small-multiple-field indigenous_status
```

Each small-multiple panel has identical node positions, so topology can be compared directly while colour carries only one meaning per panel. Small multiples cannot be combined with composite node-encoding flags. `--plot-node-spacing` accepts `compact`, `normal` (default), or `expanded`; it applies deterministic collision removal before components are packed.

### Cluster sizes by patient

Several sequences from one person (a re-test, a serial sample) are not evidence of transmission. If your metadata has a column identifying the person, name it with `--patient-column`:

```bash
pixi run hcv-cluster run -i samples.fasta --metadata metadata.csv --patient-column patient_id
```

`cluster_size` is then the number of distinct **patients**, so a cluster whose sequences all come from one person behaves like a singleton: it is not drawn as a cluster, gets no persistent label, and the sample audit says why. `clusters.csv` always reports `cluster_sequences`, `cluster_samples` and `cluster_patients` as well, with or without this option. `--sample-column` names a column identifying the physical sample (one patient can give several), so that `cluster_samples` counts distinct samples. A sample with no patient value counts as its own patient. Clustering itself is unchanged: sequences are still linked by distance, and only the counting differs.

### GUI

A Streamlit GUI wraps the same pipeline for users who prefer not to use the command line:

```bash
pixi run hcv-cluster-gui
```

This opens a browser tab where you can upload an input FASTA, set the common options (threshold, distance metric, region, threads), cache reference sequences, and browse/download the resulting cluster tables. It calls the same code as `hcv-cluster` — nothing is duplicated or reimplemented, so both stay in sync automatically.

The FASTA uploader has an optional metadata CSV uploader beside it. The GUI validates the metadata before the run and uses the same `--metadata` pipeline path as the CLI.

After a run, choose **All genotypes** or an individual genotype (and a metric, if `--distance both`) and click **View / update network**. The all-genotype view can either pack all disconnected clusters together or group them spatially by genotype. This toggle changes only the display: clustering remains genotype-stratified and there are no cross-genotype distance comparisons or links.

The displayed metadata can then be changed repeatedly without rerunning clustering. **Composite nodes** provides independent dropdowns for colour, shape, outline, centre mark, and optional advanced size. **Apply epidemiology view** maps location to fill colour, age range to shape, indigenous status to outline, and injecting status to a centre mark when those fields are available; **Reset to cluster view** restores the familiar cluster colours. **Small multiples** shows up to four copies of the same network at once, using identical positions and one colour-coded metadata field per panel. Additional fields can always be selected for hover.

The **Node spacing** selector offers compact, normal, and expanded collision-aware layouts. Positions are cached for the selected graph, layout, singleton setting, and spacing, so changing metadata does not make nodes jump around. Sample ID and cluster details always remain in hover. High-cardinality encodings produce readability warnings rather than silently changing the data.

Tick **Hide singletons** to restrict the plot (and its stats) to sequences in a multi-member cluster. **Show cluster ID labels** adds one text label beside each multi-member cluster (rather than one per node); singleton IDs are omitted to avoid clutter. These labels are retained as editable text in SVG exports and are repeated in each small-multiple panel. **Save image (PNG)** downloads a raster copy, while **Save editable vector (SVG)** downloads an Illustrator-editable vector containing paths and text rather than a flattened bitmap. Plotly's SVG grouping is preserved, although it is not a hand-authored Illustrator layer hierarchy. Sample names remain hover-only and are therefore not printed as labels in either static format. Use the **Exit** button in the sidebar to shut the server down cleanly from the browser instead of returning to the terminal.

---

## Running with Nextflow

The same pipeline runs as a Nextflow workflow, for servers, HPC clusters, or embedding in a larger pipeline. Nextflow only wires files together: every scientific decision (thresholds, regions, supported subtypes, sample bookkeeping) is made by the same Python code that `hcv-cluster run` uses, so the two give identical tables. `tests/test_nextflow_parity.py` checks this (`pixi run -e nf test`).

```bash
# Local, tools from the pixi environment
pixi run -e nf nextflow run main.nf --input samples.fasta --outdir results

# In a container built from pixi.lock (identical tool versions)
docker build -t hcv-cluster .
nextflow run main.nf -profile docker --input samples.fasta --outdir results

# Quick check that an installation works (uses a small synthetic data set)
nextflow run main.nf -profile test,docker --outdir test_results
```

Parameters mirror the CLI flags with underscores (`--distance both`, `--region core-e2`, `--min_coverage 0.8`, `--genotype 1a`, `--registry label_registry --update_registry true`); see `nextflow.config` for the full list. Reference records ship in `hcv_cluster/data/genbank`, so no network access is needed. Use `-resume` to skip steps whose inputs have not changed, and `--reuse_alignments <previous outdir>` to re-align only new or changed sequences (the same behaviour as `--reuse-alignments` on the CLI). Profiles: `standard` (tools on `PATH`), `docker`, `singularity`, `slurm`.

Network plots use the same options as the CLI with underscores: `--plot_network html --plot_scope both --plot_color_by location`, where the colour field is a column of your `--metadata` CSV. Interactive HTML works everywhere, including the container; PNG and SVG additionally need a local Chrome, which the container does not include, so they are skipped there with a warning.

Releases and versioning are described in [docs/RELEASING.md](docs/RELEASING.md) and [CHANGELOG.md](CHANGELOG.md).

To use the clustering inside another Nextflow pipeline, include the workflow rather than running `main.nf` (a complete example is in `examples/embed/`; also load `conf/hcv_cluster.config` for the parameter defaults):

```groovy
include { HCV_CLUSTER } from './hcv-cluster/workflows/hcv_cluster.nf'
HCV_CLUSTER(fasta_channel, metadata_channel)   // metadata: a CSV, or assets/NO_FILE
```

The stages are also available as `pixi run hcv-stage <stage>` (or `python -m hcv_cluster.stages <stage>`) (`resolve-ids`, `plan-genotypes`, `select-variants`, `cluster`, `merge`, `label`, `audit`, `plot`), alongside `python -m hcv_cluster.genotyping` and `python -m hcv_cluster.prep prep-align`, for running steps by hand or from another workflow engine.

**Persistent labels in a workflow.** A process cannot modify its input directory, so the `LABEL` step reads the registry and, only with `--update_registry true`, publishes the updated files back to the `--registry` directory. Run labelling-enabled pipelines one at a time against a given registry. The platform that launches runs is responsible for serialising them and for keeping the registry directory backed up.

**Environments.** `pixi.toml` pins the scientific tools once, in the `core` feature, and three environments share that pin: `default` (development, tests, GUI), `nf` (adds Nextflow) and `pipeline` (only what the container installs). Bump tool versions deliberately, then re-check results; do not let re-locking move them. TN93 1.1 adds exact speed-ups for large data sets that give identical links to the default mode (checked on 3,000 sequences); enable them with `--extra-tn93-args "-S"` (or `--extra_tn93_args '-S'` in Nextflow).

---

## Region expressions

The `-r/--region` argument accepts a flexible expression describing which part of the genome to extract and align. The default `e1-e2` retains full E1 and E2, including HVR1. Use `e1-e2-nohvr1` to retain all of E1 and E2 except the first 81 nt of E2. The E1–E2 threshold of 0.03 is an operational, conservative starting point rather than a cutoff validated for this exact full-length region; sensitivity analysis remains recommended.

**Individual regions** (in genome order):

| Name | Description |
|------|-------------|
| `core` | Core/capsid protein |
| `e1` | Envelope glycoprotein 1 |
| `e2` | Envelope glycoprotein 2 |
| `p7` | Ion channel |
| `ns2` | NS2 protease |
| `ns3` | NS3 helicase/protease |
| `ns4a` | NS4A cofactor |
| `ns4b` | NS4B membrane protein |
| `ns5a` | NS5A replication complex |
| `ns5b` | NS5B RNA-dependent RNA polymerase |

**Named HVR1 variants** (each is its own top-level region name):

| Name | Description | Default threshold |
|------|-------------|-------------------|
| `e1-e2` | E1 through E2, **HVR1 included** — the default region | 0.03 (provisional) |
| `e1-e2-nohvr1` | E1 through E2, **HVR1 masked** | 0.03 (provisional) |
| `core-e2-nohvr1` | Core through E2, **HVR1 masked** | 0.03 |
| `core-e2` | Core through E2, **HVR1 included**, contiguous | 0.045 |

**Named presets:**

| Preset | Equivalent | Description |
|--------|-----------|-------------|
| `structural` | `core-e2-nohvr1` | Core through end of E2, HVR1 masked |
| `envelope` | `e1-e2` | Both envelope glycoproteins (default region) |
| `nonstructural` | `ns2-ns5b` | NS2 through NS5B |
| `cds` | — | Whole coding sequence (polyprotein) |

**Range syntax** — `first-last` selects all regions from `first` through `last` (inclusive):

```
e1-e2          E1 + E2
ns3-ns5b       NS3 through NS5B
core-ns5b      Entire polyprotein by region boundaries
```

**Union syntax** — `a+b` includes both expressions independently:

```
core-e2-nohvr1+ns3        Default structural region plus NS3
e1-e2+ns5a                Envelope glycoproteins plus NS5A
core-e2-nohvr1+ns3+ns5b   Three separate segments
```

Ranges and unions can be combined. Overlapping or adjacent segments are merged automatically.

---

## Outputs

```
results/
  genotypes.csv                         genotype assignment for every input sequence
  genotype_fastas/
    hcv_1a.fasta                        raw sequences for each detected passing genotype
    hcv_3a.fasta
  by_genotype/
    1a/
      prep.clustering.fasta             aligned, region-trimmed sequences (input to distance step)
      prep.aligned.fasta                full reference-anchored alignment including reference
      prep.aligned.fasta.meta.json      input fingerprints enabling --reuse-alignments
      prep.qc.csv                       per-sequence QC metrics (coverage, N-fraction, pass/fail)
      tn93.csv                          pairwise TN93 distances (if --distance tn93 or both)
      snp.csv                           pairwise SNP distances + snp_count (if --distance snp or both)
      links.csv                         pairs within threshold (primary metric)
      clusters.csv                      cluster membership for this genotype
      snp_links.csv                     SNP-based pairs (if --distance both)
      snp_clusters.csv                  SNP-based clusters (if --distance both)
    3a/
      ...
  sample_audit.csv                      one row per input FASTA record: overall pass/fail, the stage that excluded it, and why
  metadata.all_samples.csv              pipeline columns (genotype, status, cluster) + your metadata, every analysed sequence
  metadata.surviving_samples.csv        the same, surviving sequences only
  surviving_samples.fasta               original sequences of survivors; headers suffixed |genotype|cluster_id
  clusters.csv                          merged cluster table, all genotypes, genotype-qualified IDs (+ cluster_label with --registry)
  cluster_labels.csv                    labels present in this run (with --registry)
  label_changes.csv                     labels created, merged or split by this run (with --registry)
  sequence_variants.csv                 repeated FASTA IDs: renames and keep/drop decisions (only if present)
  input.variants_resolved.fasta         input with repeated IDs resolved (only if present)
  clusters.snp.csv                      merged SNP clusters (if --distance both)
  links.csv                             merged edge list, all genotypes (source, target, distance)
  links.snp.csv                         merged SNP edge list, with snp_count (if --distance both)
```

Cluster IDs in `clusters.csv` are genotype-qualified (e.g. `1a_C0001`, `3a_C0001`) so per-genotype component numbers remain unique after merging. `links.csv` is a plain concatenation of every genotype's edges — sample IDs already match `clusters.csv` directly, since distances are never computed across genotypes (see below), so there's nothing to merge conflict on.

### Persistent cluster labels

Run-local IDs such as `1a_C0001` are renumbered whenever samples are added. For labels that survive between runs (and can anchor reports), point every run at a registry directory:

```bash
pixi run hcv-cluster run -i samples.fasta -o results --registry label_registry --update-registry
```

Network plots and the GUI name clusters by their persistent label (singletons and `UNASSIGNED` clusters keep their run-local ID). The GUI run tab has matching registry fields.

Without `--update-registry` the registry is only read: known clusters show their saved label, and clusters that would be newly named show `UNASSIGNED`. Use that for exploratory runs; only commit runs you want reports to be based on.

Labels look like `HCV-2608-001`: the prefix (`--label-prefix`), the **year and month the label was created** (YYMM — *not* when the outbreak began or samples were collected), and a counter within that month. Rules:

- A cluster keeps its label as it grows, even though its run-local ID changes.
- Singletons are never labelled; a sequence gets a label only once it is within threshold of another sequence.
- If clusters join, the **oldest** label survives and the others are retired with `merged_into` recorded.
- If a cluster falls apart (samples removed, QC changed), the largest piece keeps the label and the other pieces get new labels with `split_from` recorded.
- Labels are only valid for the clustering that created them. The registry records the metric, region and threshold, and refuses a run that uses different ones — use a separate registry directory for a different definition.
- Samples are tracked by FASTA ID (original ID for mixed-infection variants) and genotype. Samples missing from a later run keep their membership.

Registry files (plain CSV, safe to diff or keep in git): `labels.csv`, `membership.csv`, `events.csv` (append-only history of every create/join/merge/split), `runs.csv`, `tracks.csv`. Registries contain sample IDs, so give them the same access controls as the input data. The registry tracks the primary distance metric (TN93 unless `--distance snp`). Other clustering methods can be tracked as separate *tracks* in the same registry, sharing one label numbering, through `hcv_cluster_labels.assign_labels`.

### Sample audit

`sample_audit.csv` accounts for every record in the input FASTA, so a missing sample can always be traced. Each row has `overall_status` (`pass`/`fail`), `failed_stage` and `overall_reason`, followed by a `<stage>_status` / `<stage>_reason` pair for each stage, in the order a sequence meets them:

| Stage | Fails when |
|---|---|
| `input` | exact repeat of an earlier record with the same FASTA ID (the first copy is analysed) |
| `genotyping` | no panel alignment, or callable query coverage / identity below `--min-query-coverage` / `--min-identity` |
| `subtype` | subtype excluded by `--genotype`, or no reference genome configured for it |
| `genotype_validation` | the k-mer check finds the sequence closer to a different subtype than the one assigned (only that sample is excluded; a mismatching `genotype=` header tag is still a hard error) |
| `region_qc` | clustering-region coverage below `--min-coverage` |
| `variant_selection` | a different same-genotype variant of the same sample was kept (only variants that passed `region_qc` compete; see repeated FASTA IDs) |
| `clustering` | survivor missing from the cluster table (should not happen) |

A failed stage sets every later stage to `not_run`. Singletons count as survivors: they passed every stage and simply have no sequence within the threshold. The audit is also written when a run stops early (e.g. nothing passes genotyping). The input metadata is split into `metadata.all_samples.csv` and `metadata.surviving_samples.csv`; `metadata.csv` remains the normalized copy of the file you supplied, which the GUI reloads.

### MicrobeTrace

The top-level `clusters.csv` (nodes) and `links.csv` (edges) can be dragged into MicrobeTrace together for a single network spanning all genotypes — color nodes by the `genotype` column to keep each genotype's disconnected components visually distinct. Per-genotype `by_genotype/<gt>/clusters.csv` + `by_genotype/<gt>/links.csv` work the same way if you only want one genotype's network. Swap in `clusters.snp.csv`/`links.snp.csv` (or the per-genotype `snp_clusters.csv`/`snp_links.csv`) for the SNP-based network when `--distance both`.

---

## Why sequences are stratified by genotype

Distances and clusters are computed **separately for each genotype**. Mixing genotypes before distance calculation produces biologically meaningless results. See [`docs/subtype_stratification.md`](docs/subtype_stratification.md) for a full explanation.

---

## Advanced: stage-by-stage commands

The `hcv-workflow` entry point exposes each pipeline stage independently for debugging or custom pipelines.

**Assign genotypes:**

```bash
pixi run hcv-workflow genotype \
  --input samples.fasta \
  --output-csv results/genotypes.csv \
  --outdir results
```

Competitive genotyping uses minimap2's `asm20` preset so divergent envelope-region
sequence remains detectable. Query coverage is calculated over callable (non-`N`)
bases after combining compatible, non-overlapping, collinear segments to the same
reference; leading, trailing, and internal `N` padding therefore does not lower the
coverage fraction. Segments on different strands, overlapping query spans,
rearranged target spans, or strongly discordant query/target gaps are not combined.

Genotype assignment uses the highest-scoring reference. Coverage QC then uses the
broadest compatible passing alignment to a reference of that winning genotype,
recorded in `coverage_ref`. Identity QC uses minimap2's gap-compressed divergence
(`de`) so a long absent genome block counts as one gap event rather than thousands
of nucleotide mismatches; `block_identity` retains the conventional PAF
matches/alignment-block value for audit. `genotypes.csv` also records the old
whole-query coordinate-span coverage as `query_span_coverage`, the number of
combined segments in `alignment_segment_count`, and the input sequence's
`non_n_bases` and `non_n_fraction`.

**Prepare a genotype-specific clustering FASTA:**

```bash
pixi run hcv-workflow prep \
  --input samples.fasta \
  --genotype 1a \
  --out-prefix results/hcv_1a
```

**Pre-cache reference records:**

```bash
pixi run hcv-workflow prep-refs --genotype all
```

**Run TN93 on an aligned FASTA:**

```bash
pixi run hcv-workflow tn93 \
  --input results/hcv_1a.clustering.fasta \
  --output results/tn93.csv \
  --tn93-threshold 1.0
```

**Create linked pairs from a distance table:**

```bash
pixi run hcv-workflow link \
  --distances results/tn93.csv \
  --output results/links.csv \
  --threshold 0.03
```

**Build clusters from linked pairs:**

```bash
pixi run hcv-workflow cluster \
  --links results/links.csv \
  --nodes-fasta results/hcv_1a.clustering.fasta \
  --output results/clusters.csv
```

---

## Code layout

```
hcv_cluster/            the Python package (also pip-installable: pyproject.toml)
  cli.py                hcv-cluster run / refs         (python -m hcv_cluster)
  stages.py             file-in/file-out stages shared by the CLI and Nextflow
  genotyping.py         minimap2 subtype assignment
  prep.py               MAFFT alignment and region extraction
  methods.py            region thresholds, SNP distance
  plotting.py           network plots (CLI flags and the plot stage)
  variants.py audit.py labels.py metadata.py viz.py gui.py workflow.py
  data/                 genotype panel, ICTV table, bundled GenBank references
streamlit_app.py        GUI entry point
main.nf workflows/ modules/ conf/ nextflow.config Dockerfile    the Nextflow workflow
```

## Tests

```bash
pixi run test            # unit and stage-chain tests
pixi run -e nf test      # also runs the CLI-vs-Nextflow parity test
```
