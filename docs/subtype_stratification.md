# Why HCV Sequences Must Be Stratified by Subtype Before Distance Calculation and Clustering

## The one-sentence version

> Comparing an HCV genotype 1a sequence to a genotype 3a sequence is like comparing two different species — the genetic differences reflect ancient evolutionary divergence, not shared transmission — so mixing them produces meaningless distances and biologically impossible clusters.

---

## Background

The clustering workflow assigns each input sequence to a genotype/subtype (e.g. 1a, 3a) using competitive minimap2 alignment, then processes each subtype **independently**: alignment → region extraction → pairwise distances → threshold linking → connected-component clusters. This document explains why mixing subtypes in any of those steps would produce incorrect results.

---

## Reason 1: The alignment step is subtype-specific — mixing breaks coordinate extraction

The preprocessing step uses `MAFFT --addfragments` to anchor query sequences against a **subtype-specific reference genome**. This anchor is what allows the tool to extract a biologically homologous region (e.g. core–E2) at consistent coordinates across all sequences.

If genotype 1a and 3a sequences were forced through the same reference, the 3a sequences — only ~75–80% identical to 1a at the nucleotide level — would align poorly or misalign. The extracted region boundaries would shift, producing sequences that are not actually homologous to each other. Any distances computed from non-homologous sequences are meaningless.

---

## Reason 2: Inter-subtype distances swamp intra-subtype transmission signal

Within a subtype, contemporary sequences differ by roughly **1–5%** in the core–E2 region. Between subtypes, sequences differ by **20–35%**. Running pairwise distances on a mixed dataset produces a strongly bimodal distribution:

```
1a_sample_A vs 1a_sample_B  →  0.008   ← within-subtype: transmission signal
1a_sample_A vs 3a_sample_X  →  0.27    ← between-subtype: evolutionary divergence
```

The clustering threshold (e.g. 0.03 for the default `core-e2-nohvr1` region; see [`docs/threshold_rationale.md`](threshold_rationale.md)) is calibrated specifically for the within-subtype distribution. In a mixed matrix it is either irrelevant (no cross-subtype pairs will link, so the threshold is redundant) or — if there are alignment artefacts producing spuriously low cross-subtype distances — it becomes a source of false positive clusters.

---

## Reason 3: The TN93 evolutionary model is invalid across subtypes

TN93 (Tamura-Nei 1993) corrects for multiple substitutions using a model that assumes the sequences being compared share the **same evolutionary parameters**: similar base composition, the same transition/transversion ratio, and a common ancestor within the timeframe of interest.

Genotype 1a and 3a diverged hundreds of years ago and have measurably different base compositions and substitution patterns. Applying a single TN93 model across them produces a distance that is neither biologically meaningful nor comparable to within-subtype distances computed separately. It is a statistical artefact, not an evolutionary estimate.

---

## Reason 4: Cluster transitivity means one wrong link contaminates everything

Cluster assignment uses connected components: if A links to B, and B links to C, then A, B, and C are assigned to the same cluster — regardless of whether A and C are genetically close to each other.

If a mixed-subtype distance matrix produced **even a single spurious low-distance pair** between a 1a sequence and a 3a sequence (through alignment artefact, model breakdown, or sequencing error), that one link would **merge the entire 1a cluster with the entire 3a cluster**. The result would be a large, biologically impossible cluster reported as a transmission network.

---

## Reason 5: A cluster is a biological claim — cross-subtype clusters are impossible

A transmission cluster is a statement: *"these patients are likely connected in a recent transmission network."* HCV subtype is fixed for the lifetime of an infection; a patient infected with genotype 1a cannot transmit genotype 3a and vice versa. A cluster that spans subtypes is not a transmission network — it is a computational error. Any public health action taken on the basis of a spurious cross-subtype cluster would be misguided.

---

## Summary table

| Step | Why stratification is required |
|------|-------------------------------|
| MAFFT alignment | Reference coordinates are subtype-specific; mixed alignment produces misaligned, non-homologous sequences |
| Region extraction | core–E2 boundaries occupy different absolute positions in each subtype's genome |
| TN93 distance | Model parameters are invalid across subtypes; estimates are not comparable between subtypes |
| SNP distance | Raw SNP counts between subtypes reflect ancient fixed differences, not recent variation |
| Threshold linking | Threshold is calibrated to within-subtype variation; cross-subtype distances violate this assumption |
| Connected components | Transitivity means one false inter-subtype link collapses entire transmission networks into one spurious cluster |
| Biological interpretation | HCV subtype cannot change during infection; cross-subtype clusters are biologically impossible |

---

## In practice

The workflow handles this automatically. The `genotype` stage assigns each sequence to a subtype and splits the input FASTA accordingly. Each downstream step (alignment, distances, clusters) runs on a single-subtype FASTA. The final merged `clusters.csv` uses genotype-qualified cluster IDs (e.g. `1a_C0001`, `3a_C0001`) so that per-subtype component numbering remains unambiguous after merging.

Users should not override this behaviour by manually combining sequences from different subtypes into a single alignment or distance run.
