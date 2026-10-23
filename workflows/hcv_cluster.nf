include { RESOLVE_IDS; GENOTYPE; PREP_ALIGN; SELECT_VARIANTS; CLUSTER; MERGE; LABEL; AUDIT; PLOT; PROVENANCE } from '../modules/hcv_cluster.nf'

/*
 * HCV genomic clustering.
 *
 *   take:  fasta    consensus sequences (one file)
 *          metadata optional metadata CSV (a file, or the assets/NO_FILE placeholder)
 *   emit:  clusters merged, labelled (if a registry is given) primary cluster table
 *          audit    one row per input record: what happened and why
 */
workflow HCV_CLUSTER {
    take:
    fasta
    metadata

    main:
    RESOLVE_IDS(fasta)
    GENOTYPE(RESOLVE_IDS.out.fasta)

    // One channel item per subtype that can be analysed (plan-genotypes set the rest aside).
    // With --reuse_alignments <previous outdir>, a subtype's earlier alignment (and its
    // fingerprint file) is staged so that only new or changed sequences go to MAFFT.
    genotype_fastas = GENOTYPE.out.fastas
        .flatten()
        .map { f ->
            def gt = f.baseName.replaceFirst(/^hcv_/, '')
            def cached = params.reuse_alignments
                ? ["prep.aligned.fasta", "prep.aligned.fasta.meta.json"]
                    .collect { file("${params.reuse_alignments}/by_genotype/${gt}/${it}") }
                : []
            tuple(gt, f, cached.every { it.exists() } ? cached : [])
        }

    PREP_ALIGN(genotype_fastas)
    SELECT_VARIANTS(
        PREP_ALIGN.out.clustering.join(PREP_ALIGN.out.qc),
        RESOLVE_IDS.out.variants,
        GENOTYPE.out.genotypes,
    )

    // distance: tn93 | snp | both (tn93 is then primary and snp is reported alongside).
    def methods = params.distance == 'both' ? ['tn93', 'snp'] : [params.distance]
    CLUSTER(SELECT_VARIANTS.out.clustering.combine(Channel.fromList(methods)), RESOLVE_IDS.out.variants, metadata)

    // Group by method -> one merge each. ifEmpty keeps the audit running when nothing clustered.
    tables = CLUSTER.out.clusters.mix(CLUSTER.out.links)
        .groupTuple()
        .map { method, files -> tuple(method, files.flatten()) }
    merge_in = Channel.fromList(methods).map { m -> tuple(m, m) }
        .join(tables, remainder: true).map { m, _m, files -> tuple(m, files ?: []) }
    MERGE(merge_in, RESOLVE_IDS.out.variants, metadata)

    primary = MERGE.out.clusters.filter { method, f -> method == methods[0] }
    if (params.registry) {
        registry = Channel.value(file(params.registry).exists() ? file(params.registry) : file("${projectDir}/assets/NO_REGISTRY"))
        LABEL(primary, registry)
        final_clusters = LABEL.out.clusters.map { m, f -> f }
    } else {
        final_clusters = primary.map { m, f -> f }
    }

    AUDIT(
        fasta,
        RESOLVE_IDS.out.fasta,
        RESOLVE_IDS.out.variants,
        GENOTYPE.out.genotypes,
        PREP_ALIGN.out.qc.map { gt, f -> f }.collect().ifEmpty([]),
        PREP_ALIGN.out.validation.collect().ifEmpty([]),
        SELECT_VARIANTS.out.decisions.collect().ifEmpty([]),
        final_clusters,
        metadata,
    )

    // Optional network plots (--plot_network html|png|svg|both|all).
    if (params.plot_network && params.plot_network != 'none') {
        primary_links = MERGE.out.links.filter { it.name == "${methods[0]}.links.csv" }
        secondary_tables = params.distance == 'both'
            ? MERGE.out.clusters.filter { m, f -> m == 'snp' }.map { m, f -> f }
                .mix(MERGE.out.links.filter { it.name == 'snp.links.csv' }).collect()
            : Channel.value([])
        PLOT(
            final_clusters,
            primary_links,
            secondary_tables,
            CLUSTER.out.clusters.mix(CLUSTER.out.links).map { m, f -> f }.collect().ifEmpty([]),
            metadata,
        )
    }

    PROVENANCE()

    emit:
    clusters = final_clusters
    audit    = AUDIT.out[0]
}
