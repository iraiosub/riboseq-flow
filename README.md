[![DOI](https://zenodo.org/badge/552979223.svg)](https://zenodo.org/doi/10.5281/zenodo.10372020)

# riboseq-flow - A Nextflow DSL2 pipeline to perform ribo-seq data analysis

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
10. [Contributing guidelines](#contributing-guidelines)

## Introduction

riboseq-flow is a Nextflow DSL2 pipeline for the analysis and quality control of ribo-seq data.

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

1. Ensure `Nextflow`(version `21.10.3` or later) and `Docker` or`Singularity` (version `3.6.4` or later) are installed on your system.
Nextflow installation instructions can be found [here](https://nf-co.re/docs/usage/installation).
We recommend using Nextflow with `Java 17.0.9` or later.

**Note:** The pipeline has been tested on with `Nextflow` versions `21.10.3`, `22.10.3`, `23.04.2` and `23.10.0`. 

2. Pull the desired version of the pipeline from the GitHub repository:

```
nextflow pull iraiosub/riboseq-flow -r v1.0.1
```

3. Run the pipeline on the provided test dataset:

Using Singularity:

```
nextflow run iraiosub/riboseq-flow -r v1.0.1 -profile test,singularity
```

or using Docker:

```
nextflow run iraiosub/riboseq-flow -r v1.0.1 -profile test,docker
```

4. Check succesful execution.


## Quick start (run the pipeline on your data)

1. Ensure `Nextflow` and `Docker`/`Singularity` are installed on your system.
2. Pull the desired version of the pipeline from the GitHub repository:

```
nextflow pull iraiosub/riboseq-flow -r v1.0.1
```

3. Create a samplesheet `samplesheet.csv` with information about the samples you would like to analyse before running the pipeline. It has to be a comma-separated file with 2 columns, and a header row as shown in the example below. 

**Note:** Only single-end read data can be used as input; if you used paired-end sequencing make sure the correct read is used fro the analysis.

```
sample,fastq
sample1,/path/to/file1.fastq.gz
sample2,/path/to/file2.fastq.gz
sample3,/path/to/file3.fastq.gz
```

4. Run the pipeline. The typical command for running the pipeline is as follows (the minimum parameters have been specified):

```
nextflow run iraiosub/riboseq-flow -r v1.0.1 \
-profile singularity,crick \
-resume \
--input samplesheet.csv \
--fasta /path/to/fasta \
--gtf /path/to/gtf \
--smallrna_fasta /path/to/smallrna_fasta \
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

The pipeline is compatible with any well-annotated organism for which a FASTA genone file and GTF annotation (and ideally rRNA and other contaminat sequences) are available. Recommended sources for these files are [GENCODE](https://www.gencodegenes.org/) and [Ensembl](https://www.ensembl.org/index.html).
**Important:** The GTF file must include UTR annotations, and the format should follow the standards set by Ensembl or GENCODE.


Simplified Option for Human and Mouse:

Use the `--org` flag to automatically download and set up reference files for human or mouse genomes, eliminating the need to manually provide them.
Available options for `--org` are `GRCh38` (human) and `GRCm39` (mouse).


When `--org` is specified, all annotation files are sourced from the paths in the [genomes.config](https://github.com/iraiosub/riboseq/blob/main/conf/genomes.config) file.

Manual Annotation File Specification:


If `--org` is not specified, the user needs to provide full paths to all required annotation files:

- `--fasta` path to the FASTA genome file
- `--gtf` path to the GTF annotation file
- `--smallrna_fasta` path to the FASTA file of abundant RNA contaminants (like rRNA). (Required if pre-mapping is enabled)
- `--star_index` path to directory of pre-built STAR index (Optional). If not provided, the pipeline will generate the STAR index.
- `--transcript_info` path to TSV file with a single representative transcript for each gene, with information on CDS start, length, end and transcript length (Optional).  If not provided, the pipeline will create it using the GTF file, selecting the longest CDS transcript per gene. The transcript IDs in this file must match those in the GTF file. An example file format can be found [here](https://github.com/iraiosub/riboseq/blob/main/assets/transcript_info/gencode.v44.primary_assembly.annotation.longest_cds.transcript_info.tsv). 
- `--transcript_fasta` path to transcripts FASTA (full sequence, including CDS and UTRs) and matches the `transcript_info` file. (Required if supplying the `transcript_info` file.). If not provided, the pipeline will generate this file for the selected representative transcripts. 


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
- `--adapter_threeprime` sequence of 3' adapter (equivalent to `-a` in `cutadapt`) (Required if trimming enabled)
- `--adapter_fiveprime` sequence of 5' adapter (equivalent to `-g` in `cutadapt`)
- `--times_trimmed` number of times a read will be adaptor trimmed (default: `1`)
- `--cut_end` number of nucleotides to be trimmed from 5' or 3' end of read (equivalent to `-u` in `cutadapt`). Supply positive value for trimming from the 5' end, and negative value for trimming from the 3'end. (default `0`, no nt are trimmed). 
Important: This step is perfomed after adapter trimming, and after UMIs have been moved to the read header.
- `--min_quality` cutoff value for trimming low-quality ends from reads (default `10`)
- `--min_readlength` minimum read length after trimming (default `20`)

If you prepared your library using a TS (template-switching-based protocol) and you haven't trimmed the non-templated nucleotides and adaptors before running this pipeline, you may use the option
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
- `--expected_length` expected read lengths range. Used to report the proportion of reads of expected lengths in the aligned reads, for the generation of riboseq QC plots, and for specifying the range of RPF lengths used for P-site identification  (default `26:32`). 

**Important:**  The `--expected_length` parameter does not filter footprints based on this length range for any other analyses, including alignment, gene-level quantification or track data generation.

#### P-site identification and quantification options

P-sites are identified with [`riboWaltz`](https://github.com/LabTranslationalArchitectomics/riboWaltz/). It is strongly recommended to check the riboWaltz method to ensure the approach is suitable for your data. 
By default, the reads must be in length bins that satisfy periodicity to be used for P-site offset calculations. 
Additionally, the user can specify the following options:

- `--skip_psite` skips P-site identification and riboWaltz diagnotics plots
- `--expected_length` expected read lengths range. Used to report the proportion of reads of expected lengths in the aligned reads, for the generation of riboseq QC plots, and for specifying the range of RPF lengths used for P-site identification  (default `26:32`). 
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
    - `*.transcriptome.dedup.sorted.bam` contains the UMI deduplicated alignments to the transcriptome in BAM format
- `riboseq_qc` contains quality-control plots informing on mapping lengths, frame, distribution around start and stop codons, rRNA proportion, duplication, fraction of useful reads, PCA plots
    - `mapping_length_analysis` contains csv files with number of raw and mapped reads by length:
        - `*.after_premap.csv` 
        - `*.before_dedup.csv` 
        - `*.after_dedup.csv`
    - `multiqc_tables` contains tsv files with sample summary metrics for multiQC
    - `read_fate` contains sample-specific html files tracking read fate through the pipeline steps, a visualisation that helps understanding useful reads yield and troubleshooting. 
    - `pca` contains PCA plots and rlog-normalised count tables. Only produced if 4 samples or more are analysed.
    - `rust_analysis` contains [`RUST`](https://www.nature.com/articles/ncomms12915) metafootprint analysis, with plots showing the Kullback–Leibler divergence (K–L) profiles stratified by footprint length, using the inferred P-sites.
- `featurecounts` contains gene-level quantification of the UMI deduplicated alignments to the genome
- `psites` contains P-sites information, codon coverage and CDS coverage tables, and ribowaltz diagnostic plots:
    - `psite_offset.tsv.gz` contains P-site offsets for each read-length for all samples
    - `*psite.tsv.gz` contains sample-specific P-site information for each read
    - `ribowaltz_qc` folder containing P-site diagnostic plots generated by RiboWaltz
    - `*.coverage_psite.tsv.gz` contains P-site counts over the CDS of representative transcripts, or over a CDS window excluding a spcified region downstream start codons and upstream stop codons
    - `codon_coverage_psite.tsv.gz` contains codon-level P-site counts for each transcript
    - `codon_coverage_rpf.tsv.gz` contains codon-level RPF counts for each transcript
    - `offset_plot` contains plots that detail P-site assignment using the ribowaltz method
- `coverage_tracks` contains track files in BED and bigWig or bedGraph format
    - `*.bigwig` contains the the UMI deduplicated genome coverage tracks in bigWig format
    - `*.bedgraph` contains the the UMI deduplicated genome coverage tracks in bedGraph format
    - `*.genome.dedup.bed.gz` contains the the UMI deduplicated alignments to the genome in BED format
    - `*.transcriptome.dedup.bed.gz` contains the the UMI deduplicated alignments to the transcriptomie in BED format
    - `psite` folder containing BED files with the genomic coordinates and counts of inferred P-sites
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
3. Make sure you have Singularity installed. Please use the same version you intend to use for running the pipeline. 
4. Run the code below to pre-download and cache the required Docker images:

```
#!/bin/sh

# A script to pre-download singularity images required by iraiosub/riboseq-flow pipeline

# change dir to your NXF_SINGULARITY_CACHEDIR path
cd /path/to/image/cache

singularity pull iraiosub-nf-riboseq-latest.img docker://iraiosub/nf-riboseq
singularity pull iraiosub-mapping-length-latest.img docker://iraiosub/nf-riboseq-qc
singularity pull iraiosub-nf-riboseq-dedup-latest.img docker://iraiosub/nf-riboseq-dedup
```

### Authors and contact

riboseq-flow is written and maintained by Ira Iosub in Prof. Jernej Ule's lab at The Francis Crick Institute. It is based on a Snakemake pipeline in collaboration with its original author, Oscar Wilkins. 
Contact email: `ira.iosub@crick.ac.uk`

### Issues and contributions

riboseq-flow is under active development by Ira Iosub. For queries related to the pipeline, raise an issue on GitHub.
If you are interested in building more functionality or want to get involved please reach out.

### Contributing guidelines

We welcome contributions to riboseq-flow.

If you wish to make an addition or change to the pipeline, please follow these steps:

1. Open an issue to detail the proposed fix or feature and select the appropriate label.
2. Create a new branch based on the `dev` branch, with a short, descriptive name e.g. `feat-colours` for making changes to a color palette
3. Modify the code exclusively on this new branch and mention the relavant issue in the commit messages.
4. When your modifications are complete, submit a pull request to the `dev` branch describing the changes. 
5. Request a review from iraiosub on your pull request.
6. The pull-request will trigger a workflow execution on GitHub Actions for continuous integration (CI) of the pipeline.
This is designed to automatically test riboseq-flow whenever a pull request is made to the main or dev branches of the repository. It ensures that the pipeline runs correctly in an Ubuntu environment, helping to catch any issues or errors early in the development process. 


