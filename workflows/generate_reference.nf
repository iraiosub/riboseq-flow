#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GENERATE_GENOME_STAR_INDEX {
    tag "genome_index"

    cpus 8
    memory '32G'
    time '24h'

    input:
        path(genome_fa)
        path(annotation_gtf)

    output:
        path("genome_index"), emit: genome_index

    script:
    """
    zcat $annotation_gtf > ${annotation_gtf.getSimpleName()}.gtf
    mkdir genome_index
    STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir genome_index \
        --genomeFastaFiles $genome_fa \
        --sjdbGTFfile ${annotation_gff.getSimpleName()}.gtf \
        --sjdbGTFfeatureExon exon \
        --sjdbOverhang 100 \
        --limitOutSJcollapsed 2000000 \
        --limitGenomeGenerateRAM=200000000000
    """
}


    
//     --outFileNamePrefix {params.outdir} \
        
    

workflow GENERATE_REFERENCE_INDEX {

    main:
    
    // Generate small RNA index
    

    // Generate genome index
    GENERATE_GENOME_STAR_INDEX(
        run_annotation.out.genome_fa,
        run_annotation.out.genome_gff
    )
    

    emit:

    genome_star_index = GENERATE_GENOME_STAR_INDEX.out.genome_index

}