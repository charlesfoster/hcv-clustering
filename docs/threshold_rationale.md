# Why the Default Clustering Thresholds Are What They Are

## The one-sentence version

> Default region is `core-e2-nohvr1` (Core-through-E2 with HVR1 masked) at TN93 ≤0.03, because that's the exact coupled specification used by the closest matching precedent — Bartlett et al. 2017, an Australian TN93 pairwise/connected-components network on the same region, the same distance metric, and the same clustering algorithm this pipeline uses; the full HVR1-inclusive `core-e2` region is kept as an explicit alternative at 0.045, matched to Lamoury et al. 2015's HVR1-inclusive Core-E2 value, since threshold and HVR1-handling are not independently swappable.

---

## Defaults

This pipeline offers two named Core-E2 variants, each paired with the threshold validated for it:

- **`core-e2-nohvr1` (default): TN93 ≤0.03**, HVR1 masked.
- **`core-e2` (alternative): TN93 ≤0.045**, HVR1 included, unchanged/contiguous.

---

## What the literature says

| Region | Source | Population / metric | Threshold | Notes |
|---|---|---|---|---|
| Core-early-E2, HVR1 **masked** (H77 nt 347–1750 minus HVR1; matches this pipeline's default `core-e2-nohvr1`) | **Bartlett et al. 2017**, *J Viral Hepat* | Australian recent-HCV cohort (38% G1a, 48% G3a), **TN93** distance, **connected-components/single-linkage network** | **0.03** (sensitivity: 0.02, 0.025) | The strongest available precedent: same distance metric and same clustering algorithm as this pipeline, applied to (near-)the same region. Directly verified — see citation details below. |
| Core-E2, HVR1 **included** (matches this pipeline's `core-e2` alternative) | Lamoury et al. 2015, *PLoS ONE* | GT1a, p-distance, 90% bootstrap + ClusterPicker | **0.045** | Tree-based/p-distance, not TN93 — an extrapolation, not a direct match. Used only because no HVR1-inclusive TN93 precedent exists. |
| Core-E2, HVR1 **excluded** | Lamoury et al. 2015 | GT1a, p-distance, 90% bootstrap + ClusterPicker | 0.030 | Consistent with Bartlett et al.'s 0.03 despite the different metric/algorithm — corroborating, not primary, evidence for `core-e2-nohvr1`. |
| NS5B | Lamoury et al. 2015 | GT1a, p-distance | 0.015 | Nothing found challenges this value. |
| Partial E1 (H77 nt 943–1288) | Rose et al. 2018, *Infect Genet Evol* | High-risk PWID cohort with longitudinal follow-up, **HIV-TRACE** (TN93-based) and PhyloPart, reported as an "optimal threshold" | **0.039** | Not currently wired into `REGION_THRESHOLDS` (this pipeline has no distinct "partial E1" region), but documented here as a fallback if you're working from E1-only amplicon data lacking Core. Do not apply 0.03/0.045 to an E1-only region, and don't merge E1-only clusters with Core-E2 clusters. We confirmed the region, tools, and threshold value directly; we have not independently verified the exact optimization procedure behind "optimal," so treat that detail as reported rather than re-derived. |
| Whole genome | Rodrigo et al. 2017 (InC3 study), *J Viral Hepat* | Early-infection full genomes, IDU cohorts (incl. Australia) | Mean patristic distance **0.01** under ML tree, 95% bootstrap | Tree-based, not TN93 — **not transferable** to this pipeline's method. Included for its important caveat below, not as a threshold source. |
| Core-E2 minus HVR1 (1104 bp) | Bartlett et al. 2019, *J Int AIDS Soc* (ATAHC/ANZ cohort) | Australia/NZ, ClusterPicker + RAxML (GTR-model tree distance), 90% bootstrap | 0.05, chosen after sensitivity-testing 1.5%–5% | Corroborating evidence that thresholds well above 0.015 are appropriate for this region. Bartlett et al. 2017 is treated as primary evidence instead, since its method (TN93 + connected components) matches this pipeline directly, vs. this study's tree+bootstrap approach. |
| E1-HVR1 | Bretaña et al. 2015, *Emerg Infect Dis* (HITS-p prisoners) | Australian PWID/prisoners, PhyloPart percentile-based patristic distance | 1a: 0.034, 3a: 0.022 subs/site | Different region/method; the only genotype-specific Australian number found. 3a's cutoff is *lower* than 1a's here. |
| HIV *pol* gene (different pathogen, cited for comparison) | Weaver et al. 2024, *bioRxiv* preprint (AUTO-TUNE) | HIV-TRACE standard default | 0.015 | The well-known HIV-TRACE convention. Included to make an explicit point: this is an HIV-*pol*-specific convention, not HCV Core-E2 evidence — don't reuse it for HCV without independent support. |

---

## HVR1 masking: what it is and why it's the default

**Definition used**: HVR1 is the N-terminal 27 amino acids (81 nt) of E2, immediately following the E1/E2 cleavage site — H77 polyprotein residues 384–410, with E2 itself starting at residue 384 (Prentoe & Bukh 2018, *Front Immunol*, "Hypervariable Region 1 in Envelope Protein 2 of Hepatitis C Virus"). `hcv_cluster_prep.HVR1_LENGTH_NT = 81` implements this directly: it masks the first 81 nt of the annotated E2 region for every genotype/reference, rather than using genotype-specific coordinates — the position (immediately after E1/E2) is structurally conserved across genotypes even though HVR1's *sequence* is not.

**Why masking is the default, not just an option**: HVR1 is the fastest-evolving part of the HCV genome (rapid immune-escape-driven divergence). Including it measurably shifts the "equivalent" clustering threshold upward — Lamoury et al. 2015 found the Core-E2 cutoff moved from 0.030 (HVR1 excluded) to 0.045 (HVR1 included) using the same method on the same samples. **The threshold and the region definition are a coupled specification, not independently swappable.** You cannot take Bartlett et al. 2017's 0.03 and apply it to an HVR1-inclusive alignment, or take Lamoury's 0.045 and apply it to an HVR1-masked one, and expect either to remain valid.

**What masking does NOT mean**: HVR1 is not discarded from your sequence data — the full aligned FASTA (`prep.aligned.fasta`) still contains it. Only the *clustering* sequence (`prep.clustering.fasta`, and hence the TN93/SNP distance calculation) excludes it by default. If your priority is very-recent-transmission or outbreak-level resolution rather than broader population connectivity, an HVR1-inclusive analysis (or ideally HVR1 deep-sequencing with a haplotype-based method — see Campo et al., *outside this pipeline's scope*) may be more informative; that's exactly what `--region core-e2` gives you, alongside its own matched threshold.

**Genotype 3a caution**: Rodrigo et al.'s InC3 whole-genome study found that Core-partial-E2 (without HVR1) plus NS5B reproduced full-genome clustering well for genotype 1a, but *not* for 3a — several 3a clusters visible only via NS3/NS5A were missed. Neither Bartlett paper independently validated 3a specifically (their cohorts included substantial 3a representation and applied one uniform rule, but did not derive/optimize the threshold separately by genotype). **Treat borderline 3a links with more caution than 1a**, and cross-check against whole-genome phylogenetics where available, particularly for consequential surveillance decisions.

---

## p-distance vs. TN93 — does the metric mismatch matter?

Bartlett et al. 2017 uses TN93 directly, so this caveat applies only to the `core-e2` (HVR1-included) alternative, anchored on Lamoury et al.'s p-distance value:

- TN93 corrects for multiple/unobserved substitutions, transition/transversion bias, and base composition. It is always ≥ p-distance for the same pair, and the gap widens with divergence — but at the low divergence levels relevant to a clustering threshold (single-digit %), the correction is modest.
- Net effect: if anything, this pushes towards a slightly higher TN93-equivalent than 0.045, not lower. We have not computed an exact p-distance→TN93 conversion (that requires transition/transversion ratios we don't have for these datasets), so we kept Lamoury's own reported value rather than invent a falsely precise adjustment.

---

## Region-to-threshold mapping

```python
REGION_THRESHOLDS = {
    "core-e2-nohvr1": 0.03,   # default region
    "core-e2": 0.045,         # HVR1 included, explicit alternative
    "ns5b": 0.015,
}
```

- **`core-e2-nohvr1` (default) → 0.03.** Bartlett et al. 2017: TN93, connected components, Australian 1a/3a-dominated cohort — the closest methodological match available.
- **`core-e2` (alternative, HVR1 included) → 0.045.** Lamoury et al. 2015's HVR1-inclusive Core-E2 p-distance value; see the metric caveat above.
- **`ns5b` → 0.015.** Lamoury et al. 2015.
- **Any other region → falls back to 0.03** (the default region's value) as the best-evidenced anchor, but this is an extrapolation. The CLI prints a warning when this fallback is used without an explicit `--threshold`:

  ```
  WARNING: no HCV-specific clustering threshold evidence for region '<region>'; using the
  core-E2 default (0.03) as a starting point. See docs/threshold_rationale.md and consider
  passing --threshold explicitly.
  ```

  If you have E1-only amplicon data, use `--threshold 0.039` with an E1-restricted `--region` (Rose et al. 2018) rather than this fallback, and keep the resulting clusters separate from any Core-E2 network — don't merge them.

## Why there's no genotype-specific threshold split

A single threshold is used across genotypes, not a per-genotype split:

1. Bartlett et al. 2017 — the primary evidence — applies one uniform threshold across a cohort dominated by both 1a and 3a; it does not derive or validate a separate per-genotype cutoff.
2. Lamoury et al.'s region-comparison numbers are GT1a-only, so there's no multi-genotype Core-E2 data to derive a split from.
3. The one genotype-specific number available for an Australian PWID population (Bretaña et al., HITS-p, different region/method) has genotype 3a's cutoff *lower* than 1a's — evidence doesn't support 3a needing a more permissive (higher) threshold.
4. Instead of a numeric split, genotype 3a gets a documented *qualitative* caution (see above), which is what the evidence actually supports — treat results more cautiously, not "apply a different number."

`--threshold` still overrides everything, per-run, if you have genotype- or population-specific reasons to deviate.

---

## References

1. Bartlett SR, Wertheim JO, Bull RA, Matthews GV, Lamoury FMJ, Scheffler K, Hellard M, Maher L, Dore GJ, Lloyd AR, Applegate TL, Grebely J. "A molecular transmission network of recent hepatitis C infection in people with and without HIV: Implications for targeted treatment strategies." *Journal of Viral Hepatitis*. 2017;24(5):404–411. doi:[10.1111/jvh.12652](https://doi.org/10.1111/jvh.12652)
2. Lamoury FMJ, Jacka B, Bartlett S, Bull RA, Wong A, Amin J, Schinkel J, Poon AF, Matthews GV, Grebely J, Dore GJ, Applegate TL. "The Influence of Hepatitis C Virus Genetic Region on Phylogenetic Clustering Analysis." *PLoS ONE*. 2015;10(7):e0131437. doi:[10.1371/journal.pone.0131437](https://doi.org/10.1371/journal.pone.0131437)
3. Rose R, Lamers SL, Massaccesi G, Osburn W, Ray SC, Thomas DL, Cox AL, Laeyendecker O. "Complex patterns of Hepatitis-C virus longitudinal clustering in a high-risk population." *Infection, Genetics and Evolution*. 2018;58:77–82. doi:[10.1016/j.meegid.2017.12.015](https://doi.org/10.1016/j.meegid.2017.12.015)
4. Rodrigo C, Eltahla AA, Bull RA, Luciani F, Grebely J, Dore GJ, Applegate T, Page K, Bruneau J, Morris MD, Cox AL, Osburn W, Kim AY, Shoukry NH, Lauer GM, Maher L, Schinkel J, Prins M, Hellard M, Lloyd AR; InC3 Study Group. "Phylogenetic analysis of full-length, early infection, hepatitis C virus genomes among people with intravenous drug use: the InC3 Study." *Journal of Viral Hepatitis*. 2017;24(1):43–52. doi:[10.1111/jvh.12616](https://doi.org/10.1111/jvh.12616)
5. Bartlett SR, Applegate TL, Jacka BP, Martinello M, Lamoury FMJ, Danta M, Bradshaw D, Shaw D, Lloyd AR, Hellard M, Dore GJ, Matthews GV, Grebely J. "A latent class approach to identify multi-risk profiles associated with phylogenetic clustering of recent hepatitis C virus infection in Australia and New Zealand from 2004 to 2015." *Journal of the International AIDS Society*. 2019;22(2):e25222. doi:[10.1002/jia2.25222](https://doi.org/10.1002/jia2.25222)
6. Bretaña NA, Boelen L, Bull R, Teutsch S, White PA, Lloyd AR, Luciani F. "Transmission of Hepatitis C Virus among Prisoners, Australia, 2005–2012." *Emerging Infectious Diseases*. 2015;21(5):765–774. doi:[10.3201/eid2105.141832](https://doi.org/10.3201/eid2105.141832)
7. Weaver S, Dávila-Conn V, Ji D, Verdonk H, Ávila-Ríos S, Leigh Brown AJ, Wertheim JO, Kosakovsky Pond SL. "AUTO-TUNE: Selecting the Distance Threshold for Inferring HIV Transmission Clusters." *bioRxiv* [preprint]. 2024. doi:[10.1101/2024.03.11.584522](https://doi.org/10.1101/2024.03.11.584522)
8. Prentoe J, Bukh J. "Hypervariable Region 1 in Envelope Protein 2 of Hepatitis C Virus: A Linchpin in Neutralizing Antibody Evasion and Viral Entry." *Frontiers in Immunology*. 2018;9:2146. doi:[10.3389/fimmu.2018.02146](https://doi.org/10.3389/fimmu.2018.02146) — HVR1 boundary definition (H77 polyprotein residues 384–410, 27 aa / 81 nt), consistently reported across multiple independent HCV E2/HVR1 structural reviews.

See also [`docs/subtype_stratification.md`](subtype_stratification.md) for why thresholds are applied per-genotype-*stratum* (separate distance matrices) even though the threshold *value* itself is not genotype-split.
