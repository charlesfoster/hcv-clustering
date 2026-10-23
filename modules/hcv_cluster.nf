// Processes for the HCV clustering pipeline.
//
// These contain no scientific logic. Thresholds, region handling, genotype support and
// sample bookkeeping are all decided by the Python package in hcv_cluster/ (stages,
// prep, genotyping), the same code that
// `hcv-cluster run` uses, so the two ways of running cannot drift apart.

// The package is run from params.code_dir (the repository, or /opt/hcv-cluster in the container).
def py(module) { "PYTHONPATH=${params.code_dir} python -m hcv_cluster.${module}" }
def stage(name) { "${py('stages')} ${name}" }
def optionally(flag, value) { value != null && value != "" ? "${flag} ${value}" : "" }

// The primary metric owns clusters.csv / links.csv; with --distance both, snp is secondary.
// Written inline in closures below, because saveAs closures cannot see process inputs.
// primary = params.distance == 'snp' ? 'snp' : 'tn93'

process RESOLVE_IDS {
    input:
    path fasta

    output:
    path "resolved.fasta", emit: fasta
    path "variants.csv",   emit: variants

    script:
    """
    ${stage('resolve-ids')} --input ${fasta} --out-fasta resolved.fasta --out-variants variants.csv
    """
}

process GENOTYPE {
    publishDir "${params.outdir}", mode: "copy", saveAs: { fn -> fn == "genotypes.csv" ? fn : null }

    input:
    path fasta

    output:
    path "genotypes.csv",                emit: genotypes
    path "genotype_fastas/hcv_*.fasta",  emit: fastas, optional: true

    script:
    """
    ${py('genotyping')} \\
        --input ${fasta} \\
        --output-csv genotypes.csv \\
        --outdir genotype_fastas \\
        ${optionally('--panel-fasta', params.panel_fasta)} \\
        --minimap2 ${params.minimap2} \\
        --preset ${params.minimap2_preset} \\
        --max-secondary ${params.max_secondary} \\
        --min-query-coverage ${params.min_query_coverage} \\
        --min-identity ${params.min_identity} \\
        --close-hit-fraction ${params.close_hit_fraction} \\
        --write-genotype-fastas \\
        ${params.extra_minimap2_args ? "--extra-minimap2-args '${params.extra_minimap2_args}'" : ""}

    # Subtypes without a reference (or not requested) are set aside, with a warning.
    ${stage('plan-genotypes')} --genotypes genotypes.csv --fasta-dir genotype_fastas \\
        ${(params.genotype ? [params.genotype].flatten() : []).collect { "--genotype ${it}" }.join(' ')} \\
        ${optionally('--reference-map', params.reference_map)}
    """
}

process PREP_ALIGN {
    tag { genotype }
    publishDir path: { "${params.outdir}/by_genotype/${genotype}" }, mode: "copy",
        saveAs: { fn -> fn.endsWith(".prep.clustering.fasta") ? null : fn.replaceFirst(/^[^.]+\./, "") }

    input:
    tuple val(genotype), path(genotype_fasta), path(cache, stageAs: "cache/*")

    output:
    tuple val(genotype), path("${genotype}.prep.clustering.fasta"), emit: clustering
    tuple val(genotype), path("${genotype}.prep.qc.csv"),           emit: qc
    path "${genotype}.prep.genotype_validation.csv",                emit: validation, optional: true
    path "${genotype}.prep.aligned.fasta*",                         emit: aligned

    script:
    """
    ${py('prep')} prep-align \\
        --input ${genotype_fasta} \\
        --genotype ${genotype} \\
        --out-prefix ${genotype}.prep \\
        --min-coverage ${params.min_coverage} \\
        --region-strategy ${params.region_strategy} \\
        --mafft ${params.mafft} \\
        --threads ${task.cpus} \\
        --genotype-validation ${params.genotype_validation} \\
        --genotype-validation-kmer-size ${params.genotype_validation_kmer_size} \\
        ${cache ? '--cached-alignment cache/prep.aligned.fasta' : ''} \\
        ${optionally('--region', params.region)} \\
        ${optionally('--refs-dir', params.refs_dir)} \\
        ${optionally('--reference-map', params.reference_map)} \\
        ${optionally('--ncbi-email', params.ncbi_email)}
    """
}

process SELECT_VARIANTS {
    tag { genotype }
    publishDir path: { "${params.outdir}/by_genotype/${genotype}" }, mode: "copy",
        saveAs: { fn -> fn.endsWith(".selected.clustering.fasta") ? "prep.clustering.fasta" : "variant_decisions.csv" }

    input:
    tuple val(genotype), path(clustering), path(qc)
    path variants
    path genotypes

    output:
    tuple val(genotype), path("${genotype}.selected.clustering.fasta"), emit: clustering
    path "${genotype}.decisions.csv",                                   emit: decisions

    script:
    """
    ${stage('select-variants')} \\
        --variants ${variants} --genotypes ${genotypes} --genotype ${genotype} \\
        --qc ${qc} --clustering-fasta ${clustering} --mode ${params.duplicate_selection} \\
        --out-fasta ${genotype}.selected.clustering.fasta --out-decisions ${genotype}.decisions.csv
    """
}

process CLUSTER {
    tag { "${genotype}:${method}" }
    // Files are <genotype>.<method>.<kind>.csv; the primary method keeps the plain names.
    publishDir path: { "${params.outdir}/by_genotype/${genotype}" }, mode: "copy", saveAs: { fn ->
        def parts = fn.tokenize(".")
        def m = parts[-3]
        def kind = parts[-2]
        if (kind == "distances") return "${m}.csv"
        m == (params.distance == "snp" ? "snp" : "tn93") ? "${kind}.csv" : "${m}_${kind}.csv"
    }

    input:
    tuple val(genotype), path(clustering_fasta), val(method)
    path variants

    output:
    tuple val(method), path("${genotype}.${method}.clusters.csv"), emit: clusters
    tuple val(method), path("${genotype}.${method}.links.csv"),    emit: links
    path "${genotype}.${method}.distances.csv",                    emit: distances

    script:
    """
    ${stage('cluster')} \\
        --method ${method} --input ${clustering_fasta} \\
        --out-distances ${genotype}.${method}.distances.csv \\
        --out-links ${genotype}.${method}.links.csv \\
        --out-clusters ${genotype}.${method}.clusters.csv \\
        --variants ${variants} \\
        --region ${params.region ?: 'e1-e2'} \\
        ${optionally('--threshold', params.threshold)} \\
        ${optionally('--snp-count-threshold', params.snp_count_threshold)} \\
        --tn93 ${params.tn93} \\
        --ambiguities ${params.ambiguities} \\
        ${optionally('--ambiguity-fraction', params.ambiguity_fraction)} \\
        ${optionally('--min-overlap', params.min_overlap)} \\
        ${params.extra_tn93_args ? "--extra-tn93-args '${params.extra_tn93_args}'" : ""}
    """
}

process MERGE {
    tag { method }
    // Files are <method>.clusters.csv / <method>.links.csv. The primary method owns
    // clusters.csv/links.csv; a secondary one gets a method suffix (clusters.snp.csv).
    publishDir "${params.outdir}", mode: "copy", saveAs: { fn ->
        def parts = fn.tokenize(".")
        def (m, kind) = [parts[0], parts[1]]
        if (m != (params.distance == "snp" ? "snp" : "tn93")) return "${kind}.${m}.csv"
        kind == "clusters" && params.registry ? null : "${kind}.csv"  // LABEL publishes labelled clusters
    }

    input:
    tuple val(method), path("tables/*")
    path variants

    output:
    tuple val(method), path("${method}.clusters.csv"), emit: clusters
    path "${method}.links.csv",                         emit: links

    script:
    """
    ${stage('merge')} --tables-dir tables --method ${method} --variants ${variants} \\
        --out-clusters ${method}.clusters.csv --out-links ${method}.links.csv
    """
}

process LABEL {
    publishDir "${params.outdir}", mode: "copy", saveAs: { fn ->
        fn.startsWith("registry_out/") ? null : (fn == "clusters.labelled.csv" ? "clusters.csv" : fn)
    }
    // Write the updated registry back to where the caller keeps it.
    publishDir "${params.registry}", mode: "copy", overwrite: true,
        saveAs: { fn -> fn.startsWith("registry_out/") ? fn.replaceFirst("registry_out/", "") : null },
        enabled: params.update_registry && params.registry

    input:
    tuple val(method), path(clusters)
    path registry

    output:
    tuple val(method), path("clusters.labelled.csv"), emit: clusters
    path "cluster_labels.csv",                         emit: labels
    path "label_changes.csv",                          emit: changes
    path "registry_out/*",                             emit: registry, optional: true

    script:
    """
    ${stage('label')} \\
        --clusters ${clusters} --registry ${registry} \\
        ${params.update_registry ? '--update-registry --registry-out registry_out' : ''} \\
        --metric ${method} --region ${params.region ?: 'e1-e2'} \\
        --region-strategy ${params.region_strategy} \\
        ${optionally('--threshold', params.threshold)} \\
        ${optionally('--snp-count-threshold', params.snp_count_threshold)} \\
        --ambiguities ${params.ambiguities} \\
        ${optionally('--ambiguity-fraction', params.ambiguity_fraction)} \\
        ${optionally('--min-overlap', params.min_overlap)} \\
        ${params.extra_tn93_args ? "--extra-tn93-args '${params.extra_tn93_args}'" : ""} \\
        --label-prefix ${params.label_prefix} \\
        ${optionally('--label-date', params.label_date)} \\
        --out-clusters clusters.labelled.csv
    """
}

process AUDIT {
    publishDir "${params.outdir}", mode: "copy"

    input:
    path input_fasta
    path resolved_fasta
    path variants
    path genotypes
    path "qc/*"
    path "validation/*"
    path "decisions/*"
    path clusters
    path metadata

    output:
    path "sample_audit.csv"
    path "metadata.all_samples.csv"
    path "metadata.surviving_samples.csv"
    path "surviving_samples.fasta"
    path "sequence_variants.csv", optional: true  // final keep/drop decisions, only if IDs repeat

    script:
    def metadataArg = metadata.name == "NO_FILE" ? "" : "--metadata ${metadata}"
    // Files are staged as <genotype>.<suffix>; the stage wants GENOTYPE=PATH pairs.
    """
    pairs() { for f in "\$1"/*; do [ -e "\$f" ] && echo "\$(basename "\$f" | sed 's/\\..*//')=\$f"; done; }
    ${stage('audit')} \\
        --input ${input_fasta} --resolved-fasta ${resolved_fasta} --variants ${variants} \\
        --genotypes ${genotypes} --clusters ${clusters} ${metadataArg} \\
        --qc \$(pairs qc) --validation \$(pairs validation) --decisions \$(pairs decisions) \\
        --min-query-coverage ${params.min_query_coverage} --min-identity ${params.min_identity} \\
        --min-coverage ${params.min_coverage} \\
        ${(params.genotype ? [params.genotype].flatten() : []).collect { "--genotype ${it}" }.join(' ')} \\
        ${optionally('--reference-map', params.reference_map)} \\
        --outdir .
    """
}
