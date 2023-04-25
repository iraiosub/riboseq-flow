# riboseq - A Nextflow pipeline to perform Ribo-seq data analysis

## Table of contents

1. [Introduction](#introduction)
2. [Pipeline summary](#pipeline-summary)
3. [Quick start (test the pipeline)](#quick-start-testing)
4. [Quick start (run the pipeline)](#quick-start-running)
5. [Pipeline parameters](#pipeline-parameters)
6. [Pipeline outputs](#pipeline-outputs)
7. [Authors and contact](#authors-contact)

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
10. PCA
11. Design of sgRNA templates to deplete unwanted abundant contaminants ([`Ribocutter`](https://www.biorxiv.org/content/10.1101/2021.07.14.451473v1.full)) (Optional)
12. MultiQC report of reads QC and mapping statistics ([`MultiQC`](https://multiqc.info/))

![Pipeline summary](img/dataflow.png "Pipeline summary")


## Quick start (test the pipeline with a minimal dataset)

1. Ensure `Nextflow` and `Docker`/`Singularity` are installed on your system.

**Note:** The pipeline has so far been tested on the Crick HPC (NEMO) with `Nextflow` version `21.10.3`. `Singularity` version `3.6.4` is required for the time being to run the pipeline.

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


## Quick start (run the pipeline)

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
--input samplesheet.csv \
--org GRCh38 \
-resume
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
If `--org` is specified, all annotations will be loaded from the paths in the [genomes.config](https://github.com/iraiosub/riboseq/blob/main/conf/genomes.config) file.

If `--org` is not specified, the user must provide paths to all required annotation files, using the following parameters:

- `--fasta` path to FASTA genome file
- `--gtf` path to GTF annotation file
- `--smallrna_fasta` path to the small RNA FASTA genome file (Required if pre-mapping is enabled)
- `--star_index` path to directory of pre-built STAR index (Optional). If not provided, the pipeline will generate the STAR index.
- `--save_reference` if generated by the pipeline, save the STAR index in the results directory


### Tool specific parameters

The pipeline allows the user to set preferred parameter values, or the option to skip pipeline steps, as detailed below. 
Most parameters have default values, which will be used by the pipeline unless the user overrides them by adding the appropriate options to the run script. 
Where a default value is missing, the user must provide an appropriate value.

#### UMI options

- `--with_umi` enables UMI-based read deduplication (default: `true`). Set to `false` if you didn't use UMIs in your protocol.
- `--skip_umi_extract` skips UMI extraction from the read in case UMIs have been moved to the headers in advance
- `--umi_extract_method` specify method to extract the UMI barcode (options: `string` (default) or `regex`)
- `--umi_pattern` specifies the UMI barcode pattern,  e.g. 'NNNNN' indicates that the first 5 nucleotides of the read are from the UMI.
- `--umi_separator` specifies the UMI barcode separator (default: `_`; `rbc:` if Ultraplex was used)

#### Read trimming and filtering options

- `--skip_preprocessing` skip the adapter and quality trimming and length filtering step
- `--adapter_threeprime` sequence of 3' adapter (equivalent to -a in `cutadapt`) (default: `AGATCGGAAGAGC`)
- `--adapter_fiveprime` sequence of 5' adapter (equivalent to -g in `cutadapt`)
- `--times_trimmed` number of times a read will be adaptor trimmed (default: `1`)
- `--min_readlength` minimum read length after trimming (default `20`)
- `--min_quality` cutoff value for trimming low-quality ends from reads (default `20`)

If you prepared your library using a TS (template-switching-based protocol) you may use the option
- `--ts_trimming` which automatically trims the first 3 nucleotides corresponding to three untemplated bases at the 5′ end of the read (e.g. corresponding to rGrGrG), and the poly-A from the 3'end of reads, without the need for the user to specify adapter sequences. Should the user desire to trim a different sequence from the 3' end of the read, they can additionally set the `--ts_adapter_threeprime` (default: `A{11}`) to a different sequence.

#### Ribocutter options

- `--guide_number` number of guides to design (default: `50`)
- `--max_reads` maximum number of reads analysed (default: `1000000`)
- `--min_read_length` minimum read length threshold for reads to be analysed
- `--ribocutter_args` string specifying additional ribocutter arguments

#### Read alignment options

- `--skip_premap` skips premapping to the small RNA genome
- `--star_args` string specifying additional STAR arguments

#### Gene-level quantification options

- `--strandedness` specifies the library strandedness (options: `forward`, `reverse` or `unstranded`) (Required)

#### P-site identification and quantification options

P-sites are identified with [`riboWaltz`](https://github.com/LabTranslationalArchitectomics/riboWaltz/). It is strongly recommended to check the riboWaltz method to ensure the approach is suitable for your data. 
By default, the reads must be in length bins that satisfy periodicity to be used for P-site offset calculations.
Additionally, the user can specify the following options:

- `--length_range` specifies the range of RPF lengths used for P-site identification (default `26:31`)
- `--psite_method` specifies method used for P-site offsets identification (options: `ribowaltz` (default) or `global_max_5end`).
     - For `ribowaltz` P-site offsets are defined using [`riboWaltz`](https://github.com/LabTranslationalArchitectomics/riboWaltz/).
     - For `global_max_5end` P-site offsets are defined by the distances between the first nucleotide of the translation initiation site and the nucleotide corresponding to the global maximum of the read length-specific 5'end profiles (the first step of the riboWaltz method). Compared to riboWaltz-defined offsets, only the 5' extremities of the reads are considered for calculation and no further offset correction is performed after read-length global maximum identification.
- `exclude_start` specifies the number of nucleotides 3' from the start codon to be excluded from CDS P-site quantification (default `42`, i.e. exclude the first 14 codons)
- `exclude_stop` specifies the number of nucleotides 5' from the stop codon to be excluded from CDS P-site quantification (default `27`, i.e. exclude the last 9 codons)

#### Optional pipeline modules

- `--skip_qc` skips mapping_length_analysis and generation of riboseq QC plots
- `--skip_ribocutter` skips Ribocutter
- `--skip_psite` skips P-site identification and riboWaltz diagnotics plots


## Pipeline outputs

The pipeline outputs results in a number of subfolders:

```
.
├── preprocessed
├── premapped
├── mapped
├── deduplicated
├── mapping_length_analysis
├── riboseq_qc
├── featurecounts
├── psites
├── ribocutter
├── multiqc
└── pipeline_info
```

### Files

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
- `mapping_length_analysis` contains csv files with number of raw and mapped reads by length:
    - `*.after_premap.csv` 
    - `*.before_dedup.csv` 
    - `*.after_dedup.csv`
- `riboseq_qc` contains quality-control plots informing on mapping lengths, frame, distribution around start and stop codons, rRNA proportion, duplication, fraction of useful reads, PCA plots etc
- `featurecounts` contains gene-level quantification of the UMI deduplicated alignments to the genome
- `psites` contains P-sites information, codon coverage and CDS coverage tables, and ribowaltz diagnostic plots
- `ribocutter` contains results of ribocutter run on the pre-processed reads of minimum length defined by user, and minimum length of 23 nt, respectively
- `multiqc` contains a MultiQC report summarising the FastQC reports, premapping and mapping logs
- `pipeline_info` contains the execution reports, traces and timelines generated by Nextflow:
    - `execution_report.html`
    - `execution_timeline.html`
    - `execution_trace.txt`


### Authors and contact

This DSL2 Nextflow pipeline is written and maintained by Ira Iosub in Prof. Jernej Ule's lab at The Francis Crick Institute. It is based on a snakemake pipeline in collaboration with its original author, Oscar Wilkins. 
To raise any issues or comments with the pipeline, please email `ira.iosub@crick.ac.uk` or raise an issue on github.
