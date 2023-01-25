# riboseq - A Nextflow pipeline to perform Ribo-seq data analysis

## Table of contents

1. [Introduction](#introduction)
2. [Pipeline summary](#pipeline-summary)
3. [Quick start (test the pipeline)](#quick-start-testing)
4. [Quick start (run the pipeline)](#quick-start-running)
5. [Pipeline parameters](#pipeline-parameters)
6. [Pipeline outputs](#pipeline-outputs)

## Introduction

riboseq is a Nextflow DSL2 pipeline for the analysis of Ribo-seq data.

## Pipeline summary

1. UMI extraction ([`UMI-tools`](https://umi-tools.readthedocs.io/en/latest/)) (Optional)
2. Adapter and quality trimming ([`Cutadapt`](https://cutadapt.readthedocs.io)) (Optional)
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Premapping to remove small RNA mapping reads ([`bowtie2`](), [`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/)) (Optional)
5. Mapping to the genome and transcriptome ([`STAR`](https://github.com/alexdobin/STAR))
6. UMI-based deduplication ([`UMI-tools`](https://umi-tools.readthedocs.io/en/latest/),[`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/)) (Optional)
7. Extensive quality control ([`mapping_length_analysis`](https://pypi.org/project/mapping-length-analysis/),[`R`](https://www.r-project.org/))
8. Gene-level RPF quantification

## Quick start (test the pipeline with a minimal dataset)

1. Ensure `Nextflow` and `Docker`/`Singularity`/`Conda` are installed on your system.

**Note:** The pipeline has so far been tested on the Crick HPC (CAMP) with Nextflow version `21.10.3`.

2. Pull the main version of the pipeline from the GitHub repository:

```
nextflow pull ulelab/riboseq -r dev
```

3. Run the pipeline on the provided test dataset:

```
nextflow run ulelab/riboseq -r dev -profile test,conda,crick
```

4. Review the results

**Note:** The test can conveniently be run using the script [`run.sh`](https://github.com/ulelab/riboseq/blob/dev/run.sh) provided.


## Quick start (run the pipeline)

1. Ensure `Nextflow` and `Docker` or `Singularity` are installed on your system
2. Pull the main version of the pipeline from the GitHub repository:

```
nextflow pull ulelab/riboseq -r dev
```

3. You will need to create a samplesheet `samplesheet.csv` with information about the samples you would like to analyse before running the pipeline. It has to be a comma-separated file with 2 columns, and a header row as shown in the example below. 


```
sample,fastq
sample1,/path/to/file1.fastq.gz
sample2,/path/to/file2.fastq.gz
sample3,/path/to/file3.fastq.gz
```

4. Run the pipeline. The typical command for running the pipeline is as follows (the minimum parameters have been specified):

```
nextflow run ulelab/riboseq -r dev \
-profile conda,crick \
--input samplesheet.csv \
--org GRCh38
```

## Pipeline parameters

### Core Nextflow arguments

- `-profile`: specifies a configuration profile. Profiles can give configuration presets for different compute environments. Options are `test`, `docker`, `docker`,`singularity` and `crick` depending on the system being used and resources available. Others can be found at [nf-core](https://github.com/nf-core/configs).
- `-resume`: specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files'contents as well.

### General parameters

- `--input` specifies the input sample sheet
- `--outdir` specifies the output results directory
    - default: `./results`
- `--tracedir` specifies the pipeline run trace directory
    - default: `./results/pipeline_info`

### Genome parameters

- `--org` specifies the organism (options are currently: `GRCh38`, `GRCm38`).
If `--org` is specified, all annotations will be loaded from the paths in the genomes config file.

If `--org` is not specified, the user must provide paths to all required annotation files, using the following parameters:

- `--fasta` path to FASTA genome file
- `--gtf` path to GTF annotation file
- `--smallrna_fasta` path to the small RNA FASTA genome file (Required if pre-mapping is enabled)
- `--star_index` path to directory of pre-built STAR index (Optional)
- `--save_reference` if generated by the pipeline, save the STAR index in the results directory

### UMI options

- `--with_umi` enables UMI-based read deduplication
- `--skip_umi_extract` skips UMI extraction from the read in case UMIs have been moved to the headers in advance
- `--umi_extract_method` specify method to extract the UMI barcode, either `string` or `regex`
- `--umi_pattern` specifies the UMI barcode pattern
- `--umi_separator` specifies the UMI barcode separator

### Read trimming and filtering options

- `--skip_trimming` skip the adapter and quality trimming step
- `--adapter_threeprime` sequence of 3' adapter (equivalent to -a in `cutadapt`) (default: `AGATCGGAAGAGC`)
- `--adapter_fiveprime` sequence of 5' adapter (equivalent to -g in `cutadapt`)
- `--times_trimmed` number of times a read will be adaptor trimmed (default: `1`)
- `--min_readlength` minimum read length after trimming (default `20`)
- `--min_quality` cutoff value for trimming low-quality ends from reads (default `20`)

If you prepared your library using a TS (template-switch protocol) you may use the convenience option
- `--ts_trimming` which automatically trims the poly-A from the 3'end of reads, and the first 3 nucleotides corresponding to rGrGrG, without the need for the user to specify the adapter sequences. 


### Read alignment options

- `--skip_premap` skips premapping to the small RNA genome

### Optional pipeline modules

- `--skip_qc` skips mapping_length_analysis and generation of QC plots and MultiQC report

## Pipeline outputs

The pipeline outputs results in a number of subfolders:

```
.
├── fastqc
├── premapped
├── mapped
├── deduplicated
├── mapping_length_analysis
├── riboseq_qc
└── pipeline_info
```

### Files

- `fastqc` contains FastQC reports
- `premapped` contains files resulting from alignment to the small RNA genome (smallrna_genome):
    - `*.bam` contains read alignments to the small RNA genome in BAM format
    - `*.bam.seqs.gz` contains the mapping locations and sequences of reads mapped to the small RNA genome
    - `*.unmapped.fastq.gz` contains the sequencing reads that did not map to the small RNA genome, in FASTQ format
    - `*.premap.log` is the bowtie2 log file
- `mapped` contains files resulting from alignment to the genome and transcriptome:
    - `*.Aligned.sortedByCoord.out.bam` contains read alignments to the genome in BAM format
    - `*.Aligned.toTranscriptome.out.bam` contains read alignments to the transcriptome in BAM format
    - `*.Log.final.out` is the STAR log output file
- `deduplicated` contains files resulting after deduplication based on genomic or transcriptomic location and UMIs:
    - `*.genome.dedup.sorted.bam` contains the UMI deduplicated alignments to the genome in BAM format
    - `*.genome.dedup.bed.gz` contains the the UMI deduplicated alignments to the genome in BED format
    - `*.transcriptome.dedup.sorted.bam` contains the UMI deduplicated alignments to the transcriptome in BAM format
    - `*.transcriptome.dedup.bed.gz` contains the the UMI deduplicated alignments to the transcriptomie in BED format
- `mapping_length_analysis` contains csv files with number of raw and mapped reads by length:
    - `*.after_premap.csv` 
    - `*.before_dedup.csv` 
    - `*.after_dedup.csv`
- `riboseq_qc` contains quality-control plots informing on mapping lengths, periodicity, dist around start and stop codons, rRNA proportion, duplication, fraction of useful reads etc
- `pipeline_info` contains the execution reports, traces and timelines generated by Nextflow:
    - `execution_report.html`
    - `execution_timeline.html`
    - `execution_trace.txt`
