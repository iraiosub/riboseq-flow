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
        python: \$(python -c "import sys; print(sys.version.split()[0])")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}


process ANALYSE_RIBOLOCO {
    tag "${sample_id}"
    label 'process_medium'

    container 'iraiosub/analyse_riboloco:latest'

    input:
    tuple val(sample_id), path(riboloco)
    path transcript_info

    output:
    tuple val(sample_id), path("*.footprint_types_per_orf.csv.gz"),                  emit: footprint_types_per_orf
    tuple val(sample_id), path("*.annotated_fractions.csv.gz"),                      emit: annotated_fractions
    tuple val(sample_id), path(".kld_div.csv.gz"), path(".kld_div_decrease.csv.gz"), emit: kld, optional: true
    tuple val(sample_id), path("*.pdf"),                                             emit: plots, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    def lengths = task.ext.riboloco_lengths ?: "27:30"
    def min_footprints = task.ext.riboloco_min_footprints ?: 10
    def min_unique_footprint_positions = task.ext.riboloco_min_unique_footprint_positions ?: 3
    def gene_names = task.ext.riboloco_gene_names ?: "SCN"

    """
    analyse_riboloco_output.R \\
        --input ${riboloco} \\
        --transcript_info ${transcript_info} \\
        --lengths ${lengths} \\
        --min_footprints ${min_footprints} \\
        --min_unique_footprint_positions ${min_unique_footprint_positions} \\
        --gene_names ${gene_names} \\
        --output ${prefix} \\
        ${args}
    """
}


process RIBOLOCO_UNMIXING {
    tag "${sample_id}"
    label 'process_high'

    container 'iraiosub/analyse_riboloco:latest'

    input:
        tuple val(sample_id), path(expected_dist), path(all_footprints)

        output:
        tuple val(sample_id), path("*.unmixing_results.csv"), emit: results, optional: true
        path "versions.yml",                                  emit: versions

        when:
        task.ext.when == null || task.ext.when

        script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${sample_id}"
        def bootstrap_number = task.ext.bootstrap_number ?: 10
        def min_count = task.ext.riboloco_min_footprints ?: 10

        """
        linear_unmixing.py \\
            --expected_dist ${expected_dist} \\
            --all_footprints ${all_footprints} \\
            --output ${prefix}.unmixing_results.csv \\
            --bootstrap_number ${bootstrap_number} \\
            --min_count ${min_count} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            pandas: \$(python -c "import pandas; print(pandas.__version__)")
            numpy: \$(python -c "import numpy; print(numpy.__version__)")
            scipy: \$(python -c "import scipy; print(scipy.__version__)")
            matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        END_VERSIONS
        """
}



