#!/usr/bin/sh
#SBATCH --job-name=nf-fetch-ngs
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu

# A script with steps to obtain test data used by the riboseq-flow test profle

WORKDIR=/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/ingolia_p1/data

# Dwnload raw fastq from SRA
# SRR23242345,RNase I ribosome profiling by OTTR of HEK293T cell lysate: cell culture dish 1 (SRR23242345), specified in ids.csv

module purge
ml Nextflow/23.04.2
ml Singularity/3.6.4
ml Graphviz/2.47.2-GCCcore-10.3.0

export NXF_SINGULARITY_CACHEDIR=/nemo/lab/ulej/home/shared/singularity
export NXF_HOME=/nemo/lab/ulej/home/users/luscomben/users/iosubi/.nextflow

nextflow pull nf-core/fetchngs -r 1.10.1

nextflow run nf-core/fetchngs -r 1.10.1 \
   -profile singularity,crick \
   --input ids.csv \
   --outdir $WORKDIR/fastq

cd $WORKDIR/fastq/fastq

# Subsample 30k reads
seqtk sample -s100 SRX19188681_SRR23242345.fastq.gz 30000 > subsampled_SRX19188681_SRR23242345.fastq
gzip subsampled_SRX19188681_SRR23242345.fastq