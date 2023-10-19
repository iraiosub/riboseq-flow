# riboseq - A Nextflow pipeline to perform Ribo-seq data analysis

## Table of contents

1. [Introduction](#introduction)
2. [Pipeline summary](#pipeline-summary)
3. [Quick start (test the pipeline)](#quick-start-testing)
4. [Quick start (run the pipeline on your data)](#quick-start-running)
5. [Pipeline parameters](#pipeline-parameters)
6. [Pipeline outputs](#pipeline-outputs)
7. [Pre-download container images](#pre-download-container-images)
8. [Authors and contact](#authors-and-contact)
9. [Issues and contributions](#issues-and-contributions)

## Introduction

riboseq is a Nextflow DSL2 pipeline for the analysis of Ribo-seq data.

## Pipeline summary

1. UMI extraction ([`UMI-tools`](https://umi-tools.readthedocs.io/en/latest/)) (Optional)
2. Adapter and quality trimming, read length filtering ([`Cutadapt`](https://cutadapt.readthedocs.io)) (Optional)
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Premapping to remove small RNA mapping reads ([`bowtie2`](), [`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/)) (Optional)
5. Mapping to the genome and transcriptome ([`STAR`](https://github.com/alexdobin/STAR))
6. UMI-based deduplication ([`UMI-tools`](https://umi-tools.readthedocs.io/en/latest/),[`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/)) (Optional)
7. Extensive quality control ([`mapping_length_analysis`](https://pypi.org/project/mapping-length-analysis/),[`R`](https://www.r-project.org/)) (Optional)
8. Gene-level RPF quantification ([`FeatureCounts`](https://subread.sourceforge.net/))
9. P-site identification, CDS occupancy and P-site diagnostics ([`riboWaltz`](https://github.com/LabTranslationalArchitectomics/riboWaltz/)) (Optional)
10. Principal Component Analysis (PCA) using gene-level read counts and P-site counts over CDS regions
11. Generation of coverage tracks ([`deepTools`](https://deeptools.readthedocs.io/en/develop/))
12. Design of sgRNA templates to deplete unwanted abundant contaminants ([`Ribocutter`](https://www.biorxiv.org/content/10.1101/2021.07.14.451473v1.full)) (Optional)
13. MultiQC report of reads QC and mapping statistics ([`MultiQC`](https://multiqc.info/))

![Pipeline summary](img/dataflow.png "Pipeline summary")


## Quick start (test the pipeline with a minimal dataset)

1. Ensure `Nextflow` and `Docker`/`Singularity` are installed on your system.

**Note:** The pipeline has been tested on with `Nextflow` versions `21.10.3`, `22.10.3` and `23.04.2`. `Singularity` version `3.6.4` is required for the time being to run the pipeline.

2. Pull the desired version of the pipeline from the GitHub repository:

```
nextflow pull iraiosub/riboseq -r main
```

3. Run the pipeline on the provided test dataset:

```
nextflow run iraiosub/riboseq -r main -profile test,singularity,crick
```

4. Review the results

**Note:** An example script for running the test is provided [`run.sh`](https://github.com/iraiosub/riboseq/blob/main/run.sh).


## Quick start (run the pipeline on your data)

1. Ensure `Nextflow` and `Docker`/`Singularity` are installed on your system
2. Pull the main version of the pipeline from the GitHub repository:

```
nextflow pull iraiosub/riboseq -r main
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
nextflow run iraiosub/riboseq -r main \
-profile singularity,crick \
-resume \
--input samplesheet.csv \
--org GRCh38 \
--strandedness forward
```

## Pipeline parameters

### Core Nextflow arguments

- `-profile`: specifies a configuration profile. Profiles can give configuration presets for different compute environments. Options are `test`, `docker`,`singularity` and `crick` depending on the system being used and resources available. Others can be found at [nf-core](https://github.com/nf-core/configs).
- `-resume`: specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files'contents as well.

### General parameters

- `--input` specifies the input sample sheet
- `--outdir` specifies the output results directory
    - default: `./results`
- `--tracedir` specifies the pipeline run trace directory
    - default: `./results/pipeline_info`

### Genome parameters

- `--org` specifies the organism (options are currently: `GRCh38`, `GRCm39`).
If `--org` is specified, all annotation files will be loaded from the paths in the [genomes.config](https://github.com/iraiosub/riboseq/blob/main/conf/genomes.config) file.

If `--org` is not specified, the user must provide paths to all required annotation files, using the following parameters:

- `--fasta` path to FASTA genome file
- `--fai` path to FASTA genome index file
- `--gtf` path to GTF annotation file
- `--smallrna_fasta` path to the small RNA FASTA genome file (Required if pre-mapping is enabled)
- `--star_index` path to directory of pre-built STAR index (Optional). If not provided, the pipeline will generate the STAR index.
- `--transcript_info` path to TSV file with a single representative transcript for each gene, with information on CDS start, length, end and transcript length (Optional). If not provided, the pipeline will generate it from the GTF file, selecting the longest CDS transcript as representative for each gene. An example illustrating the format can be found [here](https://github.com/iraiosub/riboseq/tree/dev-nf/assets/transcript_info). 


### Tool specific parameters

The pipeline allows the user to set preferred parameter values, or the option to skip pipeline steps, as detailed below. 
Most parameters have default values, which will be used by the pipeline unless the user overrides them by adding the appropriate options to the run script. 
Where a default value is missing, the user must provide an appropriate value.

#### UMI options

- `--with_umi` enables UMI-based read deduplication. Use if you used UMIs in your protocol. By default, not enabled.
- `--skip_umi_extract` if using UMIs, skips UMI extraction from the read in case UMIs have been moved to the headers in advance
- `--umi_extract_method` specify method to extract the UMI barcode (options: `string` (default) or `regex`)
- `--umi_pattern` specifies the UMI barcode pattern,  e.g. 'NNNNN' indicates that the first 5 nucleotides of the read are from the UMI.
- `--umi_separator` specifies the UMI barcode separator (default: `_`; `rbc:` if Ultraplex was used)

#### Read trimming and filtering options

- `--skip_trimming` skip the adapter and quality trimming and length filtering step
- `--adapter_threeprime` sequence of 3' adapter (equivalent to -a in `cutadapt`) (default: `AGATCGGAAGAGC`)
- `--adapter_fiveprime` sequence of 5' adapter (equivalent to -g in `cutadapt`)
- `--times_trimmed` number of times a read will be adaptor trimmed (default: `1`)
- `--min_readlength` minimum read length after trimming (default `20`)
- `--min_quality` cutoff value for trimming low-quality ends from reads (default `20`)

If you prepared your library using a TS (template-switching-based protocol) you may use the option
- `--ts_trimming` which automatically trims the first 3 nucleotides corresponding to three untemplated bases at the 5′ end of the read (e.g. corresponding to rGrGrG), and the poly-A from the 3'end of reads, without the need for the user to specify adapter sequences. Should the user desire to trim a different sequence from the 3' end of the read, they can additionally set the `--ts_adapter_threeprime` (default: `A{11}`) to a different sequence.

#### Ribocutter options

- `--skip_ribocutter` skips Ribocutter
- `--guide_number` number of guides to design (default: `50`)
- `--max_reads` maximum number of reads analysed (default: `1000000`)
- `--min_read_length` minimum read length threshold for reads to be analysed
- `--ribocutter_args` string specifying additional ribocutter arguments

#### Read alignment options

- `--skip_premap` skips premapping to the small RNA genome
- `--star_args` string specifying additional STAR arguments

#### Gene-level quantification options

- `--strandedness` specifies the library strandedness (options: `forward`, `reverse` or `unstranded`) (Required)

#### Ribo-seq quality control options

- `--skip_qc` skips mapping_length_analysis and generation of riboseq QC plots
- `--expected_length` expected read lengths range. Used to report the proportion of reads of expected lengths in the aligned reads, and for the generation of riboseq QC plots (default `26:32`).

#### P-site identification and quantification options

P-sites are identified with [`riboWaltz`](https://github.com/LabTranslationalArchitectomics/riboWaltz/). It is strongly recommended to check the riboWaltz method to ensure the approach is suitable for your data. 
By default, the reads must be in length bins that satisfy periodicity to be used for P-site offset calculations. 
Additionally, the user can specify the following options:

- `--skip_psite` skips P-site identification and riboWaltz diagnotics plots
- `--length_range` specifies the range of RPF lengths used for P-site identification (default `26:32`)
- `--periodicity_threshold` specifies the periodicity threshold for read lengths to be considered for P-site identification (default `50`)
- `--psite_method` specifies method used for P-site offsets identification (options: `ribowaltz` (default) or `global_max_5end`).
     - For `ribowaltz` P-site offsets are defined using [`riboWaltz`](https://github.com/LabTranslationalArchitectomics/riboWaltz/).
     - For `global_max_5end` P-site offsets are defined by the distances between the first nucleotide of the translation initiation site and the nucleotide corresponding to the global maximum of the read length-specific 5'end profiles (the first step of the riboWaltz method). Compared to riboWaltz-defined offsets, only the 5' extremities of the reads are considered for calculation and no further offset correction is performed after read-length global maximum identification.
- `exclude_start` specifies the number of nucleotides 3' from the start codon to be excluded from CDS P-site quantification (default `42`, i.e. exclude the first 14 codons)
- `exclude_stop` specifies the number of nucleotides 5' from the stop codon to be excluded from CDS P-site quantification (default `27`, i.e. exclude the last 9 codons)

P-sites and information are identified and reported using transcriptomic coordinates. For the analyses, a representative transcript is selected for each gene. The selection is based on the following hierarchy: CDS length > total length > number of exons > 5'UTR length > 3'UTR length.
The selection is performed automatically by the pipeline using the information in the provided GTF file, and stored in `*.longest_cds.transcript_info.tsv`


#### Coverage tracks options

- `--bin_size` spcifies bin size for calculating coverage (i.e. the number of nt per bin). Bins are short consecutive counting windows of a defined size. (default `1`)
- `--track_format` specifies output file type. Either “bigwig” or “bedgraph” (default `bigwig`)

## Pipeline outputs

The pipeline outputs results in a number of subfolders:

```
.
├── annotation
├── preprocessed
├── fastqc
├── premapped
├── mapped
├── deduplicated
├── riboseq_qc
├── featurecounts
├── psites
├── coverage_tracks
├── ribocutter
├── multiqc
└── pipeline_info
```

### Files

- `annotation` contains information on the representative transcript per gene used for riboseq QC and P-site analyses, as well as bowtie2 and STAR indexes used by the pipeline
    - `*.longest_cds.transcript_info.tsv` a TSV file with a single representative transcript for each gene, with information on CDS start, length and end
    - `bowtie2` and `star` folders with indexes used for aligning reads
- `preprocessed` contains reads that have been pre-processed according to user settings for UMI extraction and trimming and filtering options
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
- `riboseq_qc` contains quality-control plots informing on mapping lengths, frame, distribution around start and stop codons, rRNA proportion, duplication, fraction of useful reads, PCA plots
    - `mapping_length_analysis` contains csv files with number of raw and mapped reads by length:
        - `*.after_premap.csv` 
        - `*.before_dedup.csv` 
        - `*.after_dedup.csv`
    - `multiqc_tables` contains tsv files with sample summary metrics for multiQC
    - `read_fate` contains sample-specific html files tracking read fate through the pipeline steps, a visualisation that helps understanding useful reads yield and troubleshooting. 
    - `pca` contains PCA plots and rlog-normalised count tables.
- `featurecounts` contains gene-level quantification of the UMI deduplicated alignments to the genome
- `psites` contains P-sites information, codon coverage and CDS coverage tables, and ribowaltz diagnostic plots:
    - `psite_offset.tsv.gz` contains P-site offsets for each read-length for all samples
    - `*psite.tsv.gz` contains sample-specific P-site information for each read
    - `ribowaltz_qc` folder containing P-site diagnostic plots generated by RiboWaltz
    - `*.coverage_psite.tsv.gz` contains P-site counts over the CDS of representative transcripts, or over a CDS window excluding a spcified region downstream start codons and upstream stop codons
    - `codon_coverage_psite.tsv.gz` contains codon-level P-site counts for each transcript
    - `codon_coverage_rpf.tsv.gz` contains codon-level RPF counts for each transcript
    - `offset_plot` contains plots that detail P-site assignment using the ribowaltz method
- `coverage_tracks` contains coverage track files in bigwig format
- `ribocutter` contains results of ribocutter run on the pre-processed reads of minimum length defined by user, and minimum length of 23 nt, respectively
- `multiqc` contains a MultiQC report summarising the FastQC reports, premapping and mapping logs
- `pipeline_info` contains the execution reports, traces and timelines generated by Nextflow:
    - `execution_report.html`
    - `execution_timeline.html`
    - `execution_trace.txt`

## Pre-download container images

When a Nextflow pipeline requires multiple Docker images, it can sometimes fail to pull them, leading to pipeline execution failures. In such scenarios, pre-downloading the container images to a designated location on your system can prevent this. 

Follow these steps to pre-download and cache the necessary images:
1. Identify the desired location on your system where you want to store the container images.
2. Set the 'NXF_SINGULARITY_CACHEDIR' environment variable to point to this chosen location.
For example, you can add the following line to your shell profile or run script:

```
export NXF_SINGULARITY_CACHEDIR=/path/to/image/cache
```

3. Run the code below to pre-download and cache the required Docker images:

```
#!/bin/sh

# A script to pre-download singularity images required by iraiosub/riboseq pipeline
ml Singularity/3.6.4

# change dir to your NFX_SINGULARITY_CACHEDIR path
cd /path/to/image/cache

singularity pull iraiosub-nf-riboseq-latest.img docker://iraiosub/nf-riboseq
singularity pull iraiosub-mapping-length-latest.img docker://iraiosub/nf-riboseq-qc
singularity pull iraiosub-nf-riboseq-dedup-latest.img docker://iraiosub/nf-riboseq-dedup
```

### Authors and contact

This DSL2 Nextflow pipeline is written and maintained by Ira Iosub in Prof. Jernej Ule's lab at The Francis Crick Institute. It is based on a snakemake pipeline in collaboration with its original author, Oscar Wilkins. 
Contact email: `ira.iosub@crick.ac.uk`

### Issues and contributions

riboseq is under active development by Ira Iosub. For queries related to the pipeline, raise an issue on github.
If you are interested in building more functionality or want to get involved please reach out.
