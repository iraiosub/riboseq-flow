#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process RIBOSEQ_QC {
    
    tag "${sample_id}"
    label 'process_medium'

    // conda '/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/riboseq_nf/riboseq/env.yml'
    container 'iraiosub/riboseq-qc:latest'

    publishDir "${params.outdir}/riboseq_qc", pattern: "*.qc_results.*", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(bam), path(bai), path(before_dedup_length_analysis), path(after_premap_length_analysis), path(after_dedup_length_analysis)
        path(transcript_info)
    
    output:
        tuple val(sample_id), path("*.qc_results.tsv.gz"), emit: qc
        path("*.qc_results.pdf"), emit: plots
        tuple val(sample_id), path("*_fq_length_mqc.tsv"), emit: fq_length_distr
        tuple val(sample_id), path("*_useful_length_mqc.tsv"), emit: useful_length_distr
        tuple val(sample_id), path("*_region_counts_mqc.tsv"), emit: region_counts
        tuple val(sample_id), path("*_start_dist_mqc.tsv"), emit: start_dist
        tuple val(sample_id), path("*_frame_mqc.tsv"), emit: frame_counts
        

    script:
        
        expected_length_range = params.expected_length

        if (params.with_umi && !params.skip_premap) 

            """
            riboseq_qc.R -b $bam -t $transcript_info -o ${sample_id} --after_premap $after_premap_length_analysis --before_dedup $before_dedup_length_analysis --after_dedup $after_dedup_length_analysis --expected_length $expected_length_range
            """
        
        else if (!params.with_umi && !params.skip_premap) 

            """
            riboseq_qc.R -b $bam -t $transcript_info -o ${sample_id} --after_premap $after_premap_length_analysis --before_dedup $before_dedup_length_analysis --expected_length $expected_length_range
            """

        else if (params.with_umi && params.skip_premap) 
            
            """
            riboseq_qc.R -b $bam -t $transcript_info -o ${sample_id} --before_dedup $before_dedup_length_analysis --after_dedup $after_dedup_length_analysis --expected_length $expected_length_range
            """

        else 
            """
            riboseq_qc.R -b $bam -t $transcript_info -o ${sample_id} --before_dedup $before_dedup_length_analysis --expected_length $expected_length_range
            """

}       


process SUMMARISE_RIBOSEQ_QC {
    
    tag "${workflow.runName}"
    label 'process_low'

    container 'iraiosub/riboseq-qc:latest'

    publishDir "${params.outdir}/riboseq_qc", pattern: "*.pdf", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/riboseq_qc/multiqc_tables", pattern: "*_mqc.tsv", mode: 'copy', overwrite: true

    input:
        path(qc_tables)
        path(fq_length_tables)
        path(useful_length_tables)
        path(region_counts_tables)
        path(start_dist_tables)
        path(mapping_counts_tables)
        path(frame_counts_tables)
    
    output:
        path("qc_summary.pdf"), emit: plot
        path("*_mqc.tsv"), emit: mqc


    script:

        """
        INPUT_QC=`echo $qc_tables | sed 's/ /,/g'`
        INPUT_FQ_LEN=`echo $fq_length_tables | sed 's/ /,/g'`
        INPUT_USEFUL_LEN=`echo $useful_length_tables | sed 's/ /,/g'`
        INPUT_REGION=`echo $region_counts_tables | sed 's/ /,/g'`
        INPUT_START_DIST=`echo $start_dist_tables | sed 's/ /,/g'`
        INPUT_MAPPING_COUNTS=`echo $mapping_counts_tables | sed 's/ /,/g'`
        INPUT_FRAME=`echo $frame_counts_tables | sed 's/ /,/g'`
        
        riboseq_qc_summary.R -i \$INPUT_QC -l \$INPUT_FQ_LEN -u \$INPUT_USEFUL_LEN -r \$INPUT_REGION --start_dist_list \$INPUT_START_DIST -m \$INPUT_MAPPING_COUNTS -f \$INPUT_FRAME
        rm $fq_length_tables $useful_length_tables $region_counts_tables $start_dist_tables $mapping_counts_tables $frame_counts_tables
        
        """

}


process TRACK_READS {
    
    tag "${sample_id}"
    label 'process_low'

    container 'iraiosub/riboseq-qc:latest'

    publishDir "${params.outdir}/riboseq_qc/read_fate", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(preprocess_log), path(premap_log), path(map_log), path(dedup_log), path(qc_log)
    
    output:
        tuple val(sample_id), path("*_sankey.html"), path("*_sankey_files"), emit: sankey
        tuple val(sample_id), path("*_mapping_counts_mqc.tsv"), emit: mapping_counts

    script:

        """
        plot_sankey.R -d .
        """


}

process PCA {
 
    label 'process_single'

    container 'iraiosub/nf-riboseq:latest'

    publishDir "${params.outdir}/riboseq_qc/pca", mode: 'copy', overwrite: true
    
    input:
    path(featurecounts_table)
    path(cds_table)
    path(cds_window_table)
    path(transcript_info)

    output:
    path("pca.pdf"), emit: pca, optional: true
    path("pca_mqc.png"), emit: pca_mqc
    // path("longest_cds_coverage_psite.tsv.gz"), emit: longest_cds_counts, optional: true
    // path("*nt_coverage_psite.tsv.gz"),  emit: longest_cds_window_counts, optional: true
    path("*rlog.tsv.gz"), emit: rlog, optional: true

    script:

        """
        plot_pca.R --featurecounts $featurecounts_table --cds $cds_table --cds_window $cds_window_table --transcript_info $transcript_info
            
        """

}

