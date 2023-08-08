#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// include { BOWTIE2_BUILD } from '../modules/local/bowtie2.nf'
include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main'   
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main'
            
workflow GENERATE_REFERENCE_INDEX {

    take:
    smallrna_fasta
    genome_fasta
    genome_gtf


    main:
    
    // Generate small RNA index
    if(!params.skip_premap) {

        BOWTIE2_BUILD(smallrna_fasta)
        smallrna_bowtie2_index = BOWTIE2_BUILD.out.index
        
    } else {

        smallrna_bowtie2_index = Channel.empty()
    }
   
    // Generate genome index

    if (!params.star_index) {

        STAR_GENOMEGENERATE(genome_fasta.map{ it[1] }, genome_gtf.map{ it[1] })
        genome_star_index = STAR_GENOMEGENERATE.out.index

    } else {

        genome_star_index = params.star_index
    }


    emit:

    smallrna_bowtie2_index
    genome_star_index


}