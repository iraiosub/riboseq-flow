#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREPARE_PSITE_TRACKS; GET_PSITE_TRACKS } from '../modules/local/ribowaltz.nf'
include { BEDTOOLS_SORT as PSITE_BEDTOOLS_SORT } from '../modules/nf-core/bedtools/sort/main'
include { UCSC_BEDGRAPHTOBIGWIG as PSITE_BEDGRAPHTOBIGWIG } from '../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow GENOMIC_PSITE_TRACKS {
    take:
    psites
    gtf
    fai
    transcript_info

    main:
    if (!(params.track_format in ['bigwig', 'bedgraph'])) {
        error("P-site track_format must be 'bigwig' or 'bedgraph'")
    }

    PREPARE_PSITE_TRACKS(gtf, fai, transcript_info)
    ch_psites = psites.map { table -> tuple([id: table.name.replaceFirst(/\.psite\.tsv\.gz$/, '')], table) }
    GET_PSITE_TRACKS(ch_psites, PREPARE_PSITE_TRACKS.out.annotation, PREPARE_PSITE_TRACKS.out.sizes)

    ch_strands = GET_PSITE_TRACKS.out.forward.map { meta, bedgraph -> tuple(meta + [strand: 'forward'], bedgraph) }
        .mix(GET_PSITE_TRACKS.out.reverse.map { meta, bedgraph -> tuple(meta + [strand: 'reverse'], bedgraph) })
    PSITE_BEDTOOLS_SORT(ch_strands, [])

    ch_tracks = PSITE_BEDTOOLS_SORT.out.sorted
    ch_versions = PREPARE_PSITE_TRACKS.out.versions.mix(GET_PSITE_TRACKS.out.versions, PSITE_BEDTOOLS_SORT.out.versions)
    if (params.track_format == 'bigwig') {
        PSITE_BEDGRAPHTOBIGWIG(PSITE_BEDTOOLS_SORT.out.sorted, PREPARE_PSITE_TRACKS.out.sizes)
        ch_tracks = PSITE_BEDGRAPHTOBIGWIG.out.bigwig
        ch_versions = ch_versions.mix(PSITE_BEDGRAPHTOBIGWIG.out.versions)
    }

    emit:
    tracks = ch_tracks
    bedgraph = PSITE_BEDTOOLS_SORT.out.sorted
    mapping_stats = GET_PSITE_TRACKS.out.mapping_stats
    versions = ch_versions
}
