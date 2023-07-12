#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { GENERATE_BOWTIE_INDEX } from '../modules/local/reference.nf'
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main'
// include { GENERATE_GENOME_STAR_INDEX } from '../modules/local/reference.nf'

            
workflow GENERATE_REFERENCE_INDEX {

    take:
    smallrna_fasta
    genome_fasta
    genome_gtf


    main:
    
    // Generate small RNA index

    if(!params.skip_premap) {

        GENERATE_BOWTIE_INDEX(smallrna_fasta)
        smallrna_bowtie2_index = GENERATE_BOWTIE_INDEX.out.smallrna_index
        
    } else {

        smallrna_bowtie2_index = Channel.empty()
    }
   
    // Generate genome index

    if (!params.star_index) {

        STAR_GENOMEGENERATE(genome_fasta, genome_gtf)
        genome_star_index = STAR_GENOMEGENERATE.out.index

    } else {

        genome_star_index = params.star_index
    }

    

    emit:

    smallrna_bowtie2_index
    genome_star_index


}