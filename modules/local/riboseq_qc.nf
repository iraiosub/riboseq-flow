#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process RIBOSEQ_QC {
    
    tag "${sample_id}"
    label 'process_medium'

    conda '/camp/lab/luscomben/home/users/iosubi/projects/riboseq_nf/riboseq/env.yml'

    publishDir "${params.outdir}/riboseq_qc", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(bam), path(bai), path(before_dedup_length_analysis), path(after_premap_length_analysis), path(after_dedup_length_analysis)
        path(transcript_info)
    
    output:
        tuple val(sample_id), path("*.qc_results.tsv.gz"), emit: qc
        path("*.qc_results.pdf"), emit: plots

    script:
            
        """
        riboseq_qc.R -b $bam -t $transcript_info -o ${sample_id} --after_premap $after_premap_length_analysis --before_dedup $before_dedup_length_analysis --after_dedup $after_dedup_length_analysis
        """

}



process SUMMARISE_RIBOSEQ_QC {
    
    tag "${workflow.runName}"
    label 'process_low'

    conda '/camp/lab/luscomben/home/users/iosubi/projects/riboseq_nf/riboseq/env.yml'

    publishDir "${params.outdir}/riboseq_qc", mode: 'copy', overwrite: true

    input:
        path(qc_tables)
    
    output:
        path("qc_summary.tsv.gz"), emit: qc_summary
        path("qc_summary.pdf"), emit: plots

    script:

        """
        INPUT_QC=`echo $qc_tables | sed 's/ /,/g'`
        
        riboseq_qc_summary.R -i \$INPUT_QC
        """

}