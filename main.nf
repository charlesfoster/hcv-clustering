#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HCV_CLUSTER } from './workflows/hcv_cluster.nf'

/*
 * HCV genomic clustering for transmission surveillance.
 *
 *   nextflow run main.nf --input sequences.fasta --outdir results
 *
 * All science lives in the Python stages; this file only wires files together. To use the
 * clustering inside a larger pipeline, include the HCV_CLUSTER workflow instead of running
 * this entry point:
 *
 *   include { HCV_CLUSTER } from './path/to/hcv-cluster/workflows/hcv_cluster.nf'
 */
workflow {
    if (!params.input) {
        error "Provide an input multi-FASTA with --input <file>"
    }
    if (params.distance !in ['tn93', 'snp', 'both']) {
        error "--distance must be tn93, snp or both (got '${params.distance}')"
    }
    if (params.patient_column && !params.metadata) {
        error "--patient_column requires --metadata <csv> containing that column"
    }
    if (params.update_registry && !params.registry) {
        error "--update_registry requires --registry <directory>"
    }

    fasta    = Channel.fromPath(params.input, checkIfExists: true).first()
    metadata = Channel.fromPath(params.metadata ?: "${projectDir}/assets/NO_FILE", checkIfExists: true).first()

    HCV_CLUSTER(fasta, metadata)
}
