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

### Step 1 — Cache reference sequences (first run only)

Download and cache the genotype-specific reference GenBank records used for alignment coordinate extraction:

```bash
pixi run hcv-cluster refs
```

References are downloaded from NCBI and stored in `refs/`. This step only needs to be repeated if you add a new genotype via `--reference-map` or want to force a re-download with `--force-download`.

### Step 2 — Run the clustering pipeline

```bash
pixi run hcv-cluster run -i samples.fasta
```

Genotypes are detected automatically. The clustering threshold default is region-dependent (0.03 for `core-e2-nohvr1`, the default region; 0.045 for `core-e2` with HVR1 included; 0.015 for `ns5b`; 0.03 as a starting point for any other region — see [`docs/threshold_rationale.md`](docs/threshold_rationale.md) for the evidence behind these values). Results are written to `results/` by default.

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
  -r core-e2-nohvr1 \
  -T 4
```

| Short | Long | Default | Description |
|-------|------|---------|-------------|
| `-i` | `--input` | required | Input multi-FASTA |
| `-o` | `--outdir` | `results` | Output directory |
| `-t` | `--threshold` | region-dependent | Maximum distance for cluster linking |
| `-d` | `--distance` | `tn93` | `tn93`, `snp`, or `both` |
| `-r` | `--region` | `core-e2-nohvr1` | Genomic region expression (see below) |
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
  --plot-symbol-by injecting_status \
  --plot-size-by age_range \
  --plot-outline-by indigenous_status \
  --plot-hover-field subtype
```

`sample_id` is required, must be unique, and must match FASTA identifiers. All other columns are retained as arbitrary plotting metadata, so future fields do not require a schema change. Blank values and plotted samples without a metadata row appear as `(missing)`; rows that do not appear in the plotted results are reported. Duplicate/blank identifiers, malformed ages, invalid headers, and reserved clustering-output column names are rejected with an actionable error. Metadata changes presentation only and cannot change distances, links, or cluster membership. A normalized copy is saved as `metadata.csv` in the result directory for later GUI rendering.

When an `age` column is supplied, the workflow validates it as a non-negative whole number and derives the ordered `age_range` categories `0–30`, `31–60`, and `61+`. You may instead supply `age_range` directly; both ASCII forms (`0-30`, `31-60`) and en-dash forms are accepted and normalized. If both columns are present, they must agree.

Colour, shape, size, and outline are independent channels and can be used simultaneously. `--plot-color-by`, `--plot-symbol-by`, `--plot-size-by`, and `--plot-outline-by` each accept a metadata field; `genotype` is also available. Shape and outline work best for a small number of categories. Size accepts numeric or intrinsically ordered data such as `age_range`; for another categorical field, repeat `--plot-size-order VALUE` from smallest to largest. Repeat `--plot-hover-field FIELD` to select several hover fields. Sample ID and cluster details are always included in interactive hover text, while static PNG/SVG files deliberately do not draw sample-name labels.

### GUI

A Streamlit GUI wraps the same pipeline for users who prefer not to use the command line:

```bash
pixi run hcv-cluster-gui
```

This opens a browser tab where you can upload an input FASTA, set the common options (threshold, distance metric, region, threads), cache reference sequences, and browse/download the resulting cluster tables. It calls the same code as `hcv-cluster` — nothing is duplicated or reimplemented, so both stay in sync automatically.

The FASTA uploader has an optional metadata CSV uploader beside it. The GUI validates the metadata before the run and uses the same `--metadata` pipeline path as the CLI.

After a run, choose **All genotypes** or an individual genotype (and a metric, if `--distance both`) and click **View / update network**. The all-genotype view can either pack all disconnected clusters together or group them spatially by genotype. This toggle changes only the display: clustering remains genotype-stratified and there are no cross-genotype distance comparisons or links.

The displayed metadata can then be changed repeatedly without rerunning clustering. Independent dropdowns map fields to node colour, shape, size, and outline, so several metadata types can be visible at once; additional fields can be selected for hover. **Apply epidemiology view** chooses conservative suggestions from the available fields, while **Reset to cluster view** restores the familiar cluster colours. Sample ID and cluster details always remain in hover. Layout positions are cached for the selected graph, layout, and singleton setting, so changing metadata does not make nodes jump around. High-cardinality encodings produce readability warnings rather than silently changing the data.

Tick **Hide singletons** to restrict the plot (and its stats) to sequences in a multi-member cluster. **Save image (PNG)** downloads a raster copy, while **Save editable vector (SVG)** downloads an Illustrator-editable vector containing paths and text rather than a flattened bitmap. Plotly's SVG grouping is preserved, although it is not a hand-authored Illustrator layer hierarchy. Sample names remain hover-only and are therefore not printed as labels in either static format. Use the **Exit** button in the sidebar to shut the server down cleanly from the browser instead of returning to the terminal.

---

## Region expressions

The `-r/--region` argument accepts a flexible expression describing which part of the genome to extract and align. The default `core-e2-nohvr1` covers Core through E2 with HVR1 (the hypervariable N-terminus of E2) masked out — the standard choice for HCV transmission surveillance, and the region the default threshold (0.03) is actually calibrated for. See [`docs/threshold_rationale.md`](docs/threshold_rationale.md) for why HVR1 is masked by default and how to get the full (HVR1-included) region instead.

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

**Named Core-E2 variants** (not selectable as individual genes; each is its own top-level region name):

| Name | Description | Default threshold |
|------|-------------|-------------------|
| `core-e2-nohvr1` | Core through E2, **HVR1 masked** — the default region | 0.03 |
| `core-e2` | Core through E2, **HVR1 included**, contiguous | 0.045 |

**Named presets:**

| Preset | Equivalent | Description |
|--------|-----------|-------------|
| `structural` | `core-e2-nohvr1` | Core through end of E2, HVR1 masked (default) |
| `envelope` | `e1-e2` | Both envelope glycoproteins |
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
  clusters.csv                          merged cluster table, all genotypes, genotype-qualified IDs
  clusters.snp.csv                      merged SNP clusters (if --distance both)
  links.csv                             merged edge list, all genotypes (source, target, distance)
  links.snp.csv                         merged SNP edge list, with snp_count (if --distance both)
```

Cluster IDs in `clusters.csv` are genotype-qualified (e.g. `1a_C0001`, `3a_C0001`) so per-genotype component numbers remain unique after merging. `links.csv` is a plain concatenation of every genotype's edges — sample IDs already match `clusters.csv` directly, since distances are never computed across genotypes (see below), so there's nothing to merge conflict on.

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

## Tests

```bash
pixi run test
```
