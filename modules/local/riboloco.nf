#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process RIBOLOCO {
    tag "${sample_id}"
    label 'process_high_memory'

    // conda "bioconda::pysam=0.23.3 bioconda::pandas=2.3.0"
    container 'iraiosub/riboloco_env:latest'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path transcript_fa
    path transcript_info

    output:
    tuple val(sample_id), path("*.csv.gz"),         emit: results
    tuple val(sample_id), path("*.summary.csv.gz"), emit: summary
    path "versions.yml",                            emit: versions

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

    stub:
    def prefix = task.ext.prefix ?: "${sample_id}"
    """
    touch ${prefix}.csv.gz
    touch ${prefix}.summary.csv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: 0.23.3
        pandas: 2.3.0
    END_VERSIONS
    """
}