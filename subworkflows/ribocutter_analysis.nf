#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { RIBOCUTTER as RIBOCUTTER_DEFAULT } from '../modules/local/ribocutter.nf' addParams(min_length: params.min_read_length)
include { RIBOCUTTER as RIBOCUTTER_MIN23 } from '../modules/local/ribocutter.nf' addParams(min_length: params.min_23_read_length)
include { GET_PROPORTION_TARGETED } from '../modules/local/ribocutter.nf'

workflow RUN_RIBOCUTTER {

    take:
    reads

    main:

    // All reads (without length filtering or rGrGrG trimmed in case of --ts_trimming)
    RIBOCUTTER_DEFAULT(
        reads
    )

    // All reads (without length filtering or rGrGrG trimmed in case of --ts_trimming) used as input, but min length 23 for ribocutter
    RIBOCUTTER_MIN23(
        reads
    )

    GET_PROPORTION_TARGETED(RIBOCUTTER_DEFAULT.out.guides.map { [ it[1] ] }.collect().mix(RIBOCUTTER_MIN23.out.guides.map { [ it[1] ] }.collect()).collect())

    emit:

    ribocutter_guides = RIBOCUTTER_DEFAULT.out.guides
    ribocutter_guides_min23 = RIBOCUTTER_MIN23.out.guides
    ribocutter_mqc = GET_PROPORTION_TARGETED.out.ribocutter_mqc

}
