#!/usr/bin/sh
#SBATCH --job-name=nf-ribo-seq-test
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu

module purge
# ml Java/11.0.2
# ml Nextflow/22.10.3
ml Nextflow/21.10.3
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

export NXF_SINGULARITY_CACHEDIR=/nemo/lab/ulej/home/users/luscomben/users/iosubi/nfcore/riboseq
# export  NXF_CONDA_CACHEDIR=/nemo/lab/ulej/home/users/luscomben/users/iosubi/nfcore/riboseq
export NXF_HOME=/nemo/lab/ulej/home/users/luscomben/users/iosubi/.nextflow # for Nextflow versions > 22.x

nextflow pull iraiosub/riboseq -r dev

nextflow run iraiosub/riboseq -r dev \
-profile crick,singularity,test \
--input /camp/lab/ulej/home/shared/riboseq_workshop/samplesheet.csv \
-resume