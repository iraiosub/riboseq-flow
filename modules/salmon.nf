#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process QUANTIFY_SALMON {
    tag "${sample_id}"
    label 'process_medium'

    conda "bioconda::salmon=1.9.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/salmon:1.9.0--h7e5ed60_1' :
    //     'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1' }"

    publishDir "${params.outdir}/salmon_quant", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(transcriptome_bam), path(transcriptome_bai)
        path(transcript_fasta)
        path(gtf)


    output:
        tuple val(sample_id), path("${sample_id}"), emit: results
        tuple val(sample_id), path("*_meta_info.json"), emit: log

    script:

        """
        salmon quant \
            --geneMap $gtf \
            --threads ${task.cpus} \
            -l SF \
            -t $transcript_fasta \
            -a $transcriptome_bam \
            -o ${sample_id}

        if [ -f ${sample_id}/aux_info/meta_info.json ]; then
            cp ${sample_id}/aux_info/meta_info.json "${sample_id}_meta_info.json"
        fi
        """
}


process TXIMPORT_SALMON {
    label "process_medium"

    conda '/camp/lab/luscomben/home/users/iosubi/projects/riboseq_nf/riboseq/env.yml'

    publishDir "${params.outdir}/salmon_quant", mode: 'copy', overwrite: true

    input:
    path(salmon_tables) // path ("salmon_quant/*")
    path(gtf)

    output:
    path("*gene_counts.tsv"), emit: counts_gene
    path("*gene_tpm.tsv.gz"), emit: tpm_gene
    path("*gene_counts_length_scaled.tsv.gz"), emit: counts_gene_length_scaled
    path("*gene_counts_scaled.tsv"), emit: counts_gene_scaled

    path("*transcript_tpm.tsv"), emit: tpm_transcript
    path("*transcript_counts.tsv"), emit: counts_transcript



    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """

    get_salmon_count_tables.R --salmon_results INPUT_TXIMPORT --gtf $gtf

    """
}
    
