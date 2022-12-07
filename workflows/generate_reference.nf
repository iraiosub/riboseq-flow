#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GENERATE_SMALL_RNA_BOWTIE_INDEX {
    tag "$smallrna_fasta"
    conda '/camp/home/iosubi/miniconda3/envs/riboseq_nf_env'

    cpus 4
    memory '16G'
    time '4h'

    input:
        path(smallrna_fasta)

    output:
        path("${smallrna_fasta.simpleName}.*.bt2"), emit: smallrna_index

    script:
    """
    bowtie2-build --threads ${task.cpus} $smallrna_fasta ${smallrna_fasta.simpleName}

    """
}



// process GENERATE_GENOME_STAR_INDEX {
//     tag "genome_index"
//     conda '/camp/home/iosubi/miniconda3/envs/riboseq_nf_env'

//     cpus 8
//     memory '32G'
//     time '24h'

//     input:
//         path(genome_fa)
//         path(annotation_gtf)

//     output:
//         path("genome_index"), emit: genome_index

//     script:
//     """
//     zcat $annotation_gtf > ${annotation_gtf.getSimpleName()}.gtf
//     mkdir genome_index
//     STAR --runThreadN ${task.cpus} \
//         --runMode genomeGenerate \
//         --genomeDir genome_index \
//         --genomeFastaFiles $genome_fa \
//         --sjdbGTFfile ${annotation_gff.getSimpleName()}.gtf \
//         --sjdbGTFfeatureExon exon \
//         --sjdbOverhang 100 \
//         --limitOutSJcollapsed 2000000 \
//         --limitGenomeGenerateRAM=200000000000
//     """
// }

//ch_smallrna_fasta = Channel.fromPath(params.smallrna_genome, checkIfExists: true)
            
workflow GENERATE_REFERENCE_INDEX {

    take:
    smallrna_fasta

    main:
    
    // Generate small RNA index
    GENERATE_SMALL_RNA_BOWTIE_INDEX(
        smallrna_fasta
    )

    // // Generate genome index
    // GENERATE_GENOME_STAR_INDEX(
    //     fasta, gtf
    // )
    

    emit:

    smallrna_bowtie2_index = GENERATE_SMALL_RNA_BOWTIE_INDEX.out.smallrna_index
    // genome_star_index = GENERATE_GENOME_STAR_INDEX.out.genome_index

}
