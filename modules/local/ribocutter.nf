#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process RIBOCUTTER {

    tag "${sample_id}"
    label 'process_single'

    // conda 'bioconda::ribocutter=0.1.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribocutter:0.1.1--pyh5e36f6f_0' :
        'quay.io/biocontainers/ribocutter:0.1.1--pyh5e36f6f_0' }"

    publishDir "${params.outdir}/ribocutter", pattern: "*.csv"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.csv"), emit: guides
    path("versions.yml"), emit: versions
    // tuple val(sample_id), path("*stats.csv"), emit: stats


    script:

    args = "-g " + params.guide_number
    args += " -r " + params.max_reads
    args += " " + params.extra_ribocutter_args


    if (!params.min_length) {

        min_read_length_arg = ""
        suffix = ""

    } else {

        min_read_length_arg = " --min_read_length " + params.min_length
        suffix = ".min" + params.min_length
    }


        """

        ribocutter -i $reads -o ${sample_id}${suffix} $min_read_length_arg $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //')
            ribocutter: \$(python -c 'from pkg_resources import get_distribution; print(get_distribution("ribocutter").version)')
        END_VERSIONS

        """
}


process GET_PROPORTION_TARGETED {

    label 'process_single'

    container 'iraiosub/nf-riboseq-qc@sha256:719e18799ff01b3071cb2187fcae78efaafe53faa604ca4a0bf224663f0cefba'

    publishDir "${params.outdir}/ribocutter", pattern: "*.{pdf,tsv}"

    input:
    path(guides)

    output:
    path("*.pdf"), emit: percentage_targeted
    path("ribocutter_mqc.tsv"), emit: ribocutter_mqc
    path("versions.yml"), emit: versions


    script:

        """

        INPUT=`echo $guides | sed 's/ /,/g'`

        echo \$INPUT

        analyse_ribocutter_guides.R -g \$INPUT
        software_versions.R "${task.process}" tidyverse data.table optparse dplyr

        """
}
