#!/usr/bin/sh
#SBATCH --job-name=nf-ribo-seq-test
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu

module purge
ml Nextflow/21.10.3
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

export NXF_SINGULARITY_CACHEDIR=/camp/lab/luscomben/home/users/iosubi/nfcore/riboseq
export  NXF_CONDA_CACHEDIR=/camp/lab/luscomben/home/users/iosubi/nfcore/riboseq
export NXF_HOME=/camp/lab/luscomben/home/users/iosubi/.nextflow # for Nextflow versions > 22.x

nextflow pull ulelab/riboseq -r dev

nextflow run ulelab/riboseq -r dev \
-profile conda,crick \
-resume \
--skip_umi_extract \
--with_umi \
--skip_trimming \
--umi_separator 'rbc:' \
--org GRCh38 \
--input ./data/samplesheet.csv \
--outdir results_full_test
