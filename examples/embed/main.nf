nextflow.enable.dsl = 2

// Example of embedding hcv-cluster in a larger Nextflow pipeline.
//
// hcv-cluster is added to the wrapper repository at ./hcv-cluster (a git submodule or a
// vendored checkout of a release tag), then its HCV_CLUSTER workflow is called like any
// other. Run:  nextflow run main.nf --source sequences.fasta --outdir results
include { HCV_CLUSTER } from './hcv-cluster/workflows/hcv_cluster.nf'

process PREPARE_INPUTS {
    input:
    path source
    output:
    path "sequences.fasta", emit: fasta
    path "metadata.csv", emit: metadata
    script:
    """
    cp ${source} sequences.fasta
    printf 'sample_id,patient_index\\nclusterA_1,P1\\nclusterA_2,P1\\n' > metadata.csv
    """
}

process SUMMARISE {
    publishDir "${params.outdir}/wrapper", mode: "copy"
    input:
    path clusters
    path audit
    output:
    path "summary.txt"
    script:
    """
    echo "clusters: \$(tail -n +2 ${clusters} | wc -l) rows; audit: \$(tail -n +2 ${audit} | wc -l) records" > summary.txt
    """
}

workflow {
    PREPARE_INPUTS(Channel.fromPath(params.source))
    HCV_CLUSTER(PREPARE_INPUTS.out.fasta, PREPARE_INPUTS.out.metadata)
    SUMMARISE(HCV_CLUSTER.out.clusters, HCV_CLUSTER.out.audit)
}
