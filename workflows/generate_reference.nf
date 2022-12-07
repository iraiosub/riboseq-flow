#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GENERATE_SMALL_RNA_BOWTIE_INDEX {
    tag "smallrna_index"

    cpus 4
    memory '16G'
    time '4h'

    input:
        path(small_rna)

    output:
        path("smallrna_index"), emit: smallrna_index

    script:
    """
    bowtie2-build --threads ${task.cpus} $small_rna smallrna_index

    """
}



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


            
workflow GENERATE_REFERENCE_INDEX {

    main:
    
    // Generate small RNA index
    GENERATE_SMALL_RNA_BOWTIE_INDEX(
        smallrna_genome
    )

    // // Generate genome index
    // GENERATE_GENOME_STAR_INDEX(
    //     fasta, gtf
    // )
    

    emit:

    smallrna_bowtie2_index = GENERATE_SMALL_RNA_BOWTIE_INDEX.out.smallrna_index
    // genome_star_index = GENERATE_GENOME_STAR_INDEX.out.genome_index

}