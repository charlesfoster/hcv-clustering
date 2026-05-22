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

## Setup

Install all dependencies into an isolated environment:

```bash
pixi install
```

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

Genotypes are detected automatically. The clustering threshold defaults are genotype-aware (1a/1b: 0.012; 2a/2b/3a and unrecognised genotypes: 0.015). Results are written to `results/` by default.

Common options:

```bash
pixi run hcv-cluster run \
  -i samples.fasta \
  -o my_results \
  -t 0.015 \
  -d both \
  -r core-e2 \
  -T 4
```

| Short | Long | Default | Description |
|-------|------|---------|-------------|
| `-i` | `--input` | required | Input multi-FASTA |
| `-o` | `--outdir` | `results` | Output directory |
| `-t` | `--threshold` | genotype-aware | Maximum distance for cluster linking |
| `-d` | `--distance` | `tn93` | `tn93`, `snp`, or `both` |
| `-r` | `--region` | `core-e2` | Genomic region expression (see below) |
| `-T` | `--threads` | `1` | MAFFT alignment threads |
| `-v` | `--verbose` | off | Debug logging |

For all advanced tuning options (minimap2/MAFFT/TN93 settings, QC thresholds, etc.):

```bash
pixi run hcv-cluster run --help-advanced
```

---

## Region expressions

The `-r/--region` argument accepts a flexible expression describing which part of the genome to extract and align. The default `core-e2` covers the complete structural region and is the standard choice for HCV transmission surveillance.

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

**Named presets:**

| Preset | Equivalent | Description |
|--------|-----------|-------------|
| `structural` | `core-e2` | Core through end of E2 |
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
core-e2+ns3        Structural region plus NS3
e1-e2+ns5a         Envelope glycoproteins plus NS5A
core-e2+ns3+ns5b   Three separate segments
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
      prep.qc.csv                       per-sequence QC metrics (coverage, N-fraction, pass/fail)
      tn93.csv                          pairwise TN93 distances (if --distance tn93 or both)
      snp.csv                           pairwise SNP distances (if --distance snp or both)
      links.csv                         pairs within threshold (primary metric)
      clusters.csv                      cluster membership for this genotype
      snp_links.csv                     SNP-based pairs (if --distance both)
      snp_clusters.csv                  SNP-based clusters (if --distance both)
    3a/
      ...
  clusters.csv                          merged cluster table, all genotypes, genotype-qualified IDs
  clusters.snp.csv                      merged SNP clusters (if --distance both)
```

Cluster IDs in `clusters.csv` are genotype-qualified (e.g. `1a_C0001`, `3a_C0001`) so per-genotype component numbers remain unique after merging.

### MicrobeTrace

`links.csv` (edges: `source`, `target`, `distance`) and `clusters.csv` (nodes: `sample_id`, `cluster_id`, `cluster_size`, `genotype`) can be imported directly into MicrobeTrace for network visualisation.

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
  --threshold 0.015
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
