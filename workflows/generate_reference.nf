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



process GENERATE_GENOME_STAR_INDEX {
    tag "$star_index"
    conda '/camp/home/iosubi/miniconda3/envs/riboseq_nf_env'

    cpus 4
    memory '64G'
    time '24h'

    publishDir "${params.outdir}/star_index", pattern: "$star_index", mode: 'copy', overwrite: true

    input:
    path(genome_fasta)
    path(genome_gtf)

    output:
    path("$star_index"), emit: star_index

    script:
    
    //zcat $annotation_gtf > ${annotation_gtf.getSimpleName()}.gtf
    //mkdir genome_index
    
    // cat $genome_gtf > ${genome_gtf.getSimpleName()}.gtf
    //--limitGenomeGenerateRAM=2000000000000

    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''

    """
    samtools faidx $genome_fasta
    NUM_BASES=`awk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${genome_fasta}.fai`

    
    mkdir star_index
    STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir star_index \
        --genomeFastaFiles $genome_fasta \
        --genomeSAindexNbases \$NUM_BASES \
        --sjdbGTFfile $genome_gtf \
        --sjdbGTFfeatureExon exon \
        --sjdbOverhang 100 \
        --limitOutSJcollapsed 2000000 \
        --limitGenomeGenerateRAM=200000000000
        

    """
}


//ch_smallrna_fasta = Channel.fromPath(params.smallrna_genome, checkIfExists: true)
            
workflow GENERATE_REFERENCE_INDEX {

    // label "high_memory"
    // cpus 8
    // memory '128G'
    // time '24h'

    take:
    smallrna_fasta
    genome_fasta
    genome_gtf

    main:
    
    // Generate small RNA index
    GENERATE_SMALL_RNA_BOWTIE_INDEX(
        smallrna_fasta
    )

    // Generate genome index
    GENERATE_GENOME_STAR_INDEX(
        genome_fasta, genome_gtf
    )
    

    emit:

    smallrna_bowtie2_index = GENERATE_SMALL_RNA_BOWTIE_INDEX.out.smallrna_index
    genome_star_index = GENERATE_GENOME_STAR_INDEX.out.star_index

}
