#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process RIBOSEQ_QC {
    
    tag "${sample_id}"
    label 'process_medium'

    // conda '/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/riboseq_nf/riboseq/env.yml'
    container 'iraiosub/mapping-length:latest'

    publishDir "${params.outdir}/riboseq_qc", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(bam), path(bai), path(before_dedup_length_analysis), path(after_premap_length_analysis), path(after_dedup_length_analysis)
        path(transcript_info)
    
    output:
        tuple val(sample_id), path("*.qc_results.tsv.gz"), emit: qc
        path("*.qc_results.pdf"), emit: plots
        tuple val(sample_id), path("*_fq_length_mqc.tsv"), emit: fq_length_distr
        tuple val(sample_id), path("*_useful_length_mqc.tsv"), emit: useful_length_distr

    script:
        
        expected_length_range = params.expected_length

        if (params.with_umi) 

            """
            riboseq_qc.R -b $bam -t $transcript_info -o ${sample_id} --after_premap $after_premap_length_analysis --before_dedup $before_dedup_length_analysis --after_dedup $after_dedup_length_analysis --expected_length $expected_length_range
            """
        
        else 
            """
            riboseq_qc.R -b $bam -t $transcript_info -o ${sample_id} --after_premap $after_premap_length_analysis --before_dedup $before_dedup_length_analysis --expected_length $expected_length_range
            """

}


process SUMMARISE_RIBOSEQ_QC {
    
    tag "${workflow.runName}"
    label 'process_low'

    // conda '/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/riboseq_nf/riboseq/env.yml'
    container 'iraiosub/mapping-length:latest'

    publishDir "${params.outdir}/riboseq_qc", mode: 'copy', overwrite: true

    input:
        path(qc_tables)
        path(fq_length_tables)
        path(useful_length_tables)
    
    output:
        path("qc_summary.pdf"), emit: plot
        path("pcoding_percentage_mqc.tsv"), emit: pcoding_percentage_mqc
        path("expected_length_mqc.tsv"), emit: expected_length_mqc
        path("duplication_mqc.tsv"), emit: duplication_mqc
        path("starting_length_mqc.tsv"), emit: length_mqc
        path("useful_length_mqc.tsv"), emit: useful_length_mqc


    script:

        """
        INPUT_QC=`echo $qc_tables | sed 's/ /,/g'`
        INPUT_FQ_LEN=`echo $fq_length_tables | sed 's/ /,/g'`
        INPUT_USEFUL_LEN=`echo $useful_length_tables | sed 's/ /,/g'`

        riboseq_qc_summary.R -i \$INPUT_QC -l \$INPUT_FQ_LEN -u \$INPUT_USEFUL_LEN
        """

}


process TRACK_READS {
    
    tag "${workflow.runName}"
    label 'process_low'

    // conda '/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/riboseq_nf/riboseq/env.yml'
    // container 'iraiosub/mapping-length:latest'

    publishDir "${params.outdir}/riboseq_qc/sankey", mode: 'copy', overwrite: true

    input:
        path(logs)
    
    output:
        path("*_sankey.html"), emit: sankey

    script:

        """
        plot_sankey.R -d .
        """


}