#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_TRANSCRIPT_INFO {
    
    tag "$gtf"
    label 'process_single'

    // conda '/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/riboseq_nf/riboseq/env2.yml'
    container 'iraiosub/nf-riboseq:latest'

    publishDir "${params.outdir}/riboseq_qc", mode: 'copy', overwrite: true

    input:
        path(gtf)

    output:
        path("*.longest_cds.transcript_info.tsv"), emit: transcript_info

    script:
    """
    ${workflow.projectDir}/bin/get_transcript_info.R -g $gtf -o 'Homo sapiens'
    
    """
}


process GENERATE_BOWTIE_INDEX {
    
    tag "$smallrna_fasta"
    label 'process_medium'


    // conda "bioconda::bowtie2=2.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0' :
        'quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0' }"

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

    // conda "bioconda::star=2.7.10a bioconda::samtools=1.16.1 conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

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