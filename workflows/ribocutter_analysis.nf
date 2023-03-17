#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { RIBOCUTTER as RIBOCUTTER } from '../modules/local/ribocutter.nf' addParams(min_length: params.min_read_length)
include { RIBOCUTTER as RIBOCUTTER_MIN23 } from '../modules/local/ribocutter.nf' addParams(min_length: params.min_23_read_length)

// Design ...

workflow RUN_RIBOCUTTER {

    take:
    reads

    main:

    // Deduplication of sequences aligned to genome
    RIBOCUTTER(
        reads
    )

    // Deduplication of sequences aligned to transcriptome
    RIBOCUTTER_MIN23(
        reads
    )

    emit:

    ribocutter_all_lengths = RIBOCUTTER.out.guides
    ribocutter_min23 = RIBOCUTTER_MIN23.out.guides


}
