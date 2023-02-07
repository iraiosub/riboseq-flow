#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GUNZIP_FASTA {
    tag "$fasta"

    conda "conda-forge::sed=4.7"

    input:
    path fasta

    output:
    path "*.fa", emit: fasta

    script:
    def gzipped = fasta.toString().endsWith('.gz')
    def outfile = gzipped ? file(fasta.baseName).baseName : fasta.baseName
    // def command = gzipped ? 'zcat' : 'cat'

    if (gzipped)
    """
    zcat $fasta > ${outfile}.fa
    """
    else
    """
    touch $fasta
    """

}


process GENERATE_SMALL_RNA_BOWTIE_INDEX {
    tag "$smallrna_fasta"
    label 'process_medium'
    // cpus 4
    // memory '16G'
    // time '4h'

    conda 'bioconda::bowtie2=2.5.0'

    input:
        path(smallrna_fasta)

    output:
        path("${smallrna_fasta.simpleName}.*.bt2"), emit: smallrna_index

    script:
    """
    bowtie2-build --threads $task.cpus $smallrna_fasta ${smallrna_fasta.simpleName}
    
    """
}



process GENERATE_GENOME_STAR_INDEX {
    tag "$genome_fasta"
    label 'process_high'
    // cpus 8
    // memory '64G'
    // time '8h'

    conda 'bioconda::star=2.7.10a bioconda::samtools=1.16.1 conda-forge::gawk=5.1.0'
    // container 'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:a7908dfb0485a80ca94e4d17b0ac991532e4e989-0'

    if (params.save_reference) {

        publishDir "${params.outdir}/star_index", mode: 'copy', overwrite: true
    }

    input:
    path(genome_fasta)
    path(genome_gtf)

    output:
    path("star_index"), emit: star_index

    script:

    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''

    // --limitGenomeGenerateRAM 8369034848

    """
    samtools faidx $genome_fasta
    NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${genome_fasta}.fai`
    
    mkdir star_index
    STAR \
        --runThreadN $task.cpus \
        --runMode genomeGenerate \
        --genomeDir star_index/ \
        --genomeFastaFiles $genome_fasta \
        --sjdbGTFfile $genome_gtf \
        --sjdbGTFfeatureExon exon \
        --genomeSAindexNbases \$NUM_BASES \
        --sjdbOverhang 100 \
        $memory
              
    """
}
            
workflow GENERATE_REFERENCE_INDEX {

    take:
    smallrna_fasta
    genome_fasta
    genome_gtf


    main:
    
    // Generate small RNA index

    if(!params.skip_premap) {
        
        GENERATE_SMALL_RNA_BOWTIE_INDEX(smallrna_fasta)
        smallrna_bowtie2_index = GENERATE_SMALL_RNA_BOWTIE_INDEX.out.smallrna_index
        
    } else {

        smallrna_bowtie2_index = Channel.empty()
    }
   
    // Generate genome index

    if (!params.star_index) {

        GENERATE_GENOME_STAR_INDEX(genome_fasta, genome_gtf)
        genome_star_index = GENERATE_GENOME_STAR_INDEX.out.star_index

    } else {

        genome_star_index = params.star_index
    }

    

    emit:

    smallrna_bowtie2_index
    genome_star_index


}