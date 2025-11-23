#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

/*
PREPARE INPUTS
*/

// Check mandatory parameters
if (!params.input) {
    exit 1, 'Input samplesheet not specified!'
} else {
    ch_input = Channel
            .fromPath( params.input )
            .splitCsv(header:true)
            .map { row -> [ row.sample, file(row.fastq, checkIfExists: true) ] }
            // .view()
}


// Check if genome exists in the config file
if (params.org  && !params.genomes.containsKey(params.org)) {
    exit 1, "The provided genome '${params.org}' is not available. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

/*
PREPARE GENOME CHANNELS
*/

if(params.org) {

    params.fasta = params.genomes[ params.org ].fasta
    params.gtf = params.genomes[ params.org ].gtf
    params.star_index = params.genomes[ params.org ].star_index
    params.contaminants_fasta = params.genomes[ params.org ].contaminants_fasta
    params.transcript_info = params.genomes[ params.org ].transcript_info
    params.transcript_fasta = params.genomes[ params.org ].transcript_fasta

}  else {

    if(!params.fasta ) { exit 1, '--fasta is not specified.' }
    if(!params.gtf ) { exit 1, '--gtf is not specified.' }
    if(!params.contaminants_fasta && !params.skip_premap ) {exit 1, '--contaminants_fasta is not specified.' }

}


ch_genome_fasta = file(params.fasta, checkIfExists: true)
ch_contaminants_fasta = file(params.contaminants_fasta, checkIfExists: true)
ch_genome_gtf = file(params.gtf, checkIfExists: true)

/*
MultiQC config
*/
// ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)

/*
Channel for OPTIONAL INPUT - contains a tuple with sample ids as keys, but no path/file
*/

ch_optional = Channel
            .fromPath( params.input )
            .splitCsv(header:true)
            .map { row -> [ row.sample, [] ] }


/*
SUBWORKFLOWS
*/

include { PREPARE_RIBOSEQ_REFERENCE } from './subworkflows/prepare_reference.nf'
include { PREPROCESS_READS } from './subworkflows/preprocess_reads.nf'
include { DEDUPLICATE } from './subworkflows/dedup.nf'
include { MAPPING_LENGTH_ANALYSES } from './subworkflows/mapping_length_analyses.nf'
include { GET_GENE_LEVEL_COUNTS } from './subworkflows/gene_level_counts.nf'
include { RUN_RIBOCUTTER } from './subworkflows/ribocutter_analysis.nf'

/*
MODULES
*/

include { FASTQC } from './modules/nf-core/fastqc/main'
include { PREMAP } from './modules/local/premap.nf'
include { MAP } from './modules/local/map.nf'
include { RIBOSEQ_QC } from './modules/local/riboseq_qc.nf'
include { SUMMARISE_RIBOSEQ_QC } from './modules/local/riboseq_qc.nf'
include { TRACK_READS } from './modules/local/riboseq_qc.nf'
include { IDENTIFY_PSITES } from './modules/local/ribowaltz.nf'
include { RUST_QC } from './modules/local/ribowaltz.nf'
include { GET_COVERAGE_TRACKS } from './modules/local/get_tracks.nf'
include { GET_PSITE_TRACKS } from './modules/local/ribowaltz.nf'
include { PCA } from './modules/local/riboseq_qc.nf'
include { MULTIQC } from './modules/local/multiqc.nf'


/*
MAIN WORKFLOW
*/

workflow RIBOSEQ {

    // Prepare annotation
    PREPARE_RIBOSEQ_REFERENCE(
        ch_genome_fasta,
        ch_genome_gtf,
        ch_contaminants_fasta
    )

    // Test to asses channel type
    // def test = PREPARE_RIBOSEQ_REFERENCE.out.genome_gtf
    //                 .map { it[1] }

    // println "genome_gtf_main class: ${test.getClass()}"

    // Extract UMIs and/or trim adapters and filter on min length, then run FASTQC
    PREPROCESS_READS(ch_input)
    FASTQC(PREPROCESS_READS.out.fastq)

    // Run ribocutter on trimmed but not length filtered (for ts_trimming, on trimmed but not rGrGrG-cut or length filtered)
    if (!params.skip_ribocutter) {
        RUN_RIBOCUTTER(
            PREPROCESS_READS.out.trimmed_fastq
            )
    }

    // Align reads
    if (!params.skip_premap) {

        // Premap to the contaminant RNA genome
        PREMAP(
            PREPROCESS_READS.out.fastq,
            PREPARE_RIBOSEQ_REFERENCE.out.contaminants_bowtie2_index.collect()
        )

        MAP(
            PREMAP.out.unmapped,
            PREPARE_RIBOSEQ_REFERENCE.out.genome_star_index.collect()
        )

    } else {

        MAP(
            PREPROCESS_READS.out.fastq,
            PREPARE_RIBOSEQ_REFERENCE.out.genome_star_index.collect()
        )
    }

    if (params.with_umi) {

        // Remove duplicate reads from BAM file based on UMIs
        DEDUPLICATE(
            MAP.out.genome_bam,
            MAP.out.transcriptome_bam
        )
    }


    // riboseq QC
    if (!params.skip_qc) {
        // Mapping length analysis, visualisation and read fate tracking

        if (!params.skip_premap && params.with_umi) {

            MAPPING_LENGTH_ANALYSES(
                MAP.out.genome_bam,
                PREPROCESS_READS.out.cut_fastq,
                DEDUPLICATE.out.dedup_genome_bam,
                PREMAP.out.unmapped
            )

            RIBOSEQ_QC(
                DEDUPLICATE.out.dedup_transcriptome_bam
                    .join(MAPPING_LENGTH_ANALYSES.out.before_dedup_length_analysis)
                    .join(MAPPING_LENGTH_ANALYSES.out.after_premap_length_analysis)
                    .join(MAPPING_LENGTH_ANALYSES.out.after_dedup_length_analysis),
                PREPARE_RIBOSEQ_REFERENCE.out.transcript_info.collect()
            )

            // Track read fate through pipeline
            TRACK_READS(
                PREPROCESS_READS.out.logs
                    .join(PREMAP.out.log)
                    .join(MAP.out.log)
                    .join(DEDUPLICATE.out.dedup_genome_log)
                    .join(RIBOSEQ_QC.out.qc)
            )

        } else if (params.skip_premap && params.with_umi) {

            MAPPING_LENGTH_ANALYSES(
                MAP.out.genome_bam,
                PREPROCESS_READS.out.cut_fastq,
                DEDUPLICATE.out.dedup_genome_bam,
                Channel.empty()
            )

            RIBOSEQ_QC(
                DEDUPLICATE.out.dedup_transcriptome_bam
                    .join(MAPPING_LENGTH_ANALYSES.out.before_dedup_length_analysis)
                    .join(MAPPING_LENGTH_ANALYSES.out.after_premap_length_analysis)
                    .join(MAPPING_LENGTH_ANALYSES.out.after_dedup_length_analysis),
                PREPARE_RIBOSEQ_REFERENCE.out.transcript_info.collect()
            )

            TRACK_READS(
                PREPROCESS_READS.out.logs
                    .join(MAP.out.log)
                    .join(DEDUPLICATE.out.dedup_genome_log)
                    .join(RIBOSEQ_QC.out.qc)
                    .join(ch_optional)
            )


        } else if (!params.skip_premap && !params.with_umi) {

            MAPPING_LENGTH_ANALYSES(
                MAP.out.genome_bam,
                PREPROCESS_READS.out.cut_fastq,
                Channel.empty(),
                PREMAP.out.unmapped
            )

            RIBOSEQ_QC(
                MAP.out.transcriptome_bam
                    .join(MAPPING_LENGTH_ANALYSES.out.before_dedup_length_analysis)
                    .join(MAPPING_LENGTH_ANALYSES.out.after_premap_length_analysis)
                    .join(MAPPING_LENGTH_ANALYSES.out.after_dedup_length_analysis),
                PREPARE_RIBOSEQ_REFERENCE.out.transcript_info.collect()
            )

            TRACK_READS(
                PREPROCESS_READS.out.logs
                    .join(PREMAP.out.log)
                    .join(MAP.out.log)
                    .join(RIBOSEQ_QC.out.qc)
                    .join(ch_optional)
            )

        } else if (params.skip_premap && !params.with_umi) {

            MAPPING_LENGTH_ANALYSES(
                MAP.out.genome_bam,
                PREPROCESS_READS.out.cut_fastq,
                Channel.empty(),
                Channel.empty()
            )

            RIBOSEQ_QC(
                MAP.out.transcriptome_bam
                    .join(MAPPING_LENGTH_ANALYSES.out.before_dedup_length_analysis)
                    .join(MAPPING_LENGTH_ANALYSES.out.after_premap_length_analysis)
                    .join(MAPPING_LENGTH_ANALYSES.out.after_dedup_length_analysis),
                PREPARE_RIBOSEQ_REFERENCE.out.transcript_info.collect()
            )

            TRACK_READS(
                PREPROCESS_READS.out.logs
                    .join(MAP.out.log)
                    .join(RIBOSEQ_QC.out.qc)
                    .join(ch_optional)
                    .join(ch_optional)
            )

        }

        ch_merge_qc = RIBOSEQ_QC.out.qc
            .map { [ it[1] ] }
            .collect()

        ch_merge_read_length = RIBOSEQ_QC.out.fq_length_distr
            .map { [ it[1] ] }
            .collect()

        ch_merge_useful_length = RIBOSEQ_QC.out.useful_length_distr
            .map { [ it[1] ] }
            .collect()

        ch_merge_region_counts = RIBOSEQ_QC.out.region_counts
            .map { [ it[1] ] }
            .collect()

        ch_merge_start_dist = RIBOSEQ_QC.out.start_dist
            .map { [ it[1] ] }
            .collect()

        ch_merge_frame = RIBOSEQ_QC.out.frame_counts
            .map { [ it[1] ] }
            .collect()

        ch_merge_mapping_counts = TRACK_READS.out.mapping_counts
            .map { [ it[1] ] }
            .collect() 

        ch_merge_lemgth_filter = TRACK_READS.out.length_filter
            .map { [ it[1] ] }
            .collect() 

        SUMMARISE_RIBOSEQ_QC(
            ch_merge_qc,
            ch_merge_read_length,
            ch_merge_useful_length,
            ch_merge_region_counts,
            ch_merge_start_dist,
            ch_merge_mapping_counts,
            ch_merge_frame,
            ch_merge_lemgth_filter
        )

    }

    // Get gene-level counts from BAM alignments using featureCounts and coverage tracks using deepTools
    if (params.with_umi) {
        GET_GENE_LEVEL_COUNTS(
            DEDUPLICATE.out.dedup_genome_bam,
            PREPARE_RIBOSEQ_REFERENCE.out.genome_gtf.map{ it[1] }
        )

        GET_COVERAGE_TRACKS(DEDUPLICATE.out.dedup_genome_bam)

    } else {
        GET_GENE_LEVEL_COUNTS(
            MAP.out.genome_bam,
            PREPARE_RIBOSEQ_REFERENCE.out.genome_gtf.map{ it[1] }
        )

        GET_COVERAGE_TRACKS(MAP.out.genome_bam)

    }

    // Get P sites
    if (!params.skip_psite) {

         if (params.with_umi) {

            IDENTIFY_PSITES(
                DEDUPLICATE.out.dedup_transcriptome_bam
                    .map { [it[1]] }
                    .collect(),
                PREPARE_RIBOSEQ_REFERENCE.out.genome_gtf
                    .map { it[1] },
                PREPARE_RIBOSEQ_REFERENCE.out.genome_fasta
                    .map { it[1] },
                PREPARE_RIBOSEQ_REFERENCE.out.transcript_info
            )


        } else {

            IDENTIFY_PSITES(
                MAP.out.transcriptome_bam
                    .map { [it[1]] }
                    .collect(),
                PREPARE_RIBOSEQ_REFERENCE.out.genome_gtf
                    .map { it[1] },
                PREPARE_RIBOSEQ_REFERENCE.out.genome_fasta
                    .map { it[1] },
                PREPARE_RIBOSEQ_REFERENCE.out.transcript_info
            )


        }

        // Convert codon coverage from transcriptomic to genomic coordinates
        if (params.genomic_codon_coverage) {
        GET_PSITE_TRACKS(
            IDENTIFY_PSITES.out.psites.flatten(),
            PREPARE_RIBOSEQ_REFERENCE.out.genome_gtf.map{ it[1] },
            PREPARE_RIBOSEQ_REFERENCE.out.genome_fai
        )


        }

    }

    if (!params.skip_psite & !params.skip_qc) {

        RUST_QC(
            PREPARE_RIBOSEQ_REFERENCE.out.transcript_info_fa,
            PREPARE_RIBOSEQ_REFERENCE.out.transcript_info,
            IDENTIFY_PSITES.out.psites.flatten()
        )
    }


    // PCA on gene-level RPF counts and transcript-level P-sites
    if (!params.skip_psite) {

        PCA(
            GET_GENE_LEVEL_COUNTS.out.merged_counts_table,
            IDENTIFY_PSITES.out.cds_coverage,
            IDENTIFY_PSITES.out.cds_window_coverage,
            PREPARE_RIBOSEQ_REFERENCE.out.transcript_info
        )

    } else {

        PCA(
            GET_GENE_LEVEL_COUNTS.out.merged_counts_table,
            Channel.empty(),
            Channel.empty(),
            PREPARE_RIBOSEQ_REFERENCE.out.transcript_info
            )
    }

    // Run MULTIQC
    // if ribocutter skipped, use empty channel in multiqc

    if (params.skip_ribocutter) {

        ch_ribocutter = Channel.empty()

    } else {

        ch_ribocutter = RUN_RIBOCUTTER.out.ribocutter_mqc
    }

    if (!params.skip_premap && !params.skip_qc) {

        ch_logs = FASTQC.out.html.join(FASTQC.out.zip)
            .map { [it[1], it[2]] }
            .collect()
            .mix(
                PREMAP.out.log.map { it[1] }.collect(),
                MAP.out.log.map { it[1] }.collect(),
                PCA.out.pca_mqc,
                SUMMARISE_RIBOSEQ_QC.out.mqc,
                ch_ribocutter
            )
            .collect()

    } else if (!params.skip_premap && params.skip_qc) {

        ch_logs = FASTQC.out.html.join(FASTQC.out.zip)
            .map { [it[1], it[2]] }
            .collect()
            .mix(
                PREMAP.out.log.map { it[1] }.collect(),
                MAP.out.log.map { it[1] }.collect(),
                PCA.out.pca_mqc,
                ch_ribocutter
            )
            .collect()

    } else if (params.skip_premap && params.skip_qc) {

        ch_logs = FASTQC.out.html.join(FASTQC.out.zip)
            .map { [it[1], it[2]] }
            .collect()
            .mix(
                MAP.out.log.map { it[1] }.collect(),
                PCA.out.pca_mqc,
                ch_ribocutter
            )
            .collect()

    } else {

        ch_logs = FASTQC.out.html.join(FASTQC.out.zip)
            .map { [it[1], it[2]] }
            .collect()
            .mix(
                MAP.out.log.map { it[1] }.collect(),
                PCA.out.pca_mqc,
                SUMMARISE_RIBOSEQ_QC.out.mqc,
                ch_ribocutter
            )
            .collect()
    }

    MULTIQC(ch_logs, ch_multiqc_config)

}


workflow {
    RIBOSEQ ()
}

/*
COMPLETION EVENTS
*/

workflow.onComplete {

    if (workflow.success) {
        log.info "-\033[0;34m[riboseq]\033[0;32m Pipeline completed successfully\033[0m-\n"
    } else {
        log.info "-\033[0;34m[riboseq]\033[1;91m Pipeline completed with errors\033[0m\n"
    }

}