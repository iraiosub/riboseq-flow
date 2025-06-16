#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process RIBOLOCO {
    tag "${sample_id}"
    label 'process_high_memory'

    container 'iraiosub/riboloco_env:latest'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path transcript_fa
    path transcript_info

    output:
    tuple val(sample_id), path("${sample_id}.csv.gz"),  emit: results
    tuple val(sample_id), path("*.summary.csv.gz"),     emit: summary
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    def flank = task.ext.riboloco_orf_flank ?: 300
    """
    riboloco_lite.py \\
        --bam ${bam} \\
        --fasta ${transcript_fa} \\
        --info ${transcript_info} \\
        --output ${prefix} \\
        --flank ${flank} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}



process ANALYSE_RIBOLOCO {
    tag "${sample_id}"
    label 'process_medium'

    container 'iraiosub/riboloco_env:latest'

    input:
    tuple val(sample_id), path(riboloco)
    path transcript_info

    output:
    tuple val(sample_id), path("*.footprint_types_per_orf.csv.gz"),                  emit: footprint_types_per_orf
    tuple val(sample_id), path("*.annotated_fractions.csv.gz"),                      emit: annotated_fractions
    tuple val(sample_id), path(".kld_div.csv.gz"), path(".kld_div_decrease.csv.gz"), emit: kld
    tuple val(sample_id), path("*.pdf"),                                             emit: plots

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    def lengths = task.ext.riboloco_lengths ?: "27:30"
    def min_footprints = task.ext.riboloco_min_footprints ?: 10
    def min_unique_footprint_positions = task.ext.riboloco_min_unique_footprint_positions ?: 3

    """
    analyse_riboloco.R \\
        --input ${riboloco} \\
        --transcript_info ${transcript_info} \\
        --lengths ${lengths} \\
        --min_footprints ${min_footprints} \\
        --min_unique_footprint_positions ${min_unique_footprint_positions} \\
        --output ${prefix} \\
        ${args}
    """
}


// RIBOLOCO_UNMIXING

// RIBOLOCO_UNMIXING_ANALYSIS

