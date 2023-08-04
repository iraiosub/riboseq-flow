#!/usr/bin/sh
#SBATCH --job-name=nf-ribo-seq-test
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu

module purge
# ml Java/11.0.2
# ml Nextflow/21.10.3
#ml Nextflow/22.10.3
ml Nextflow/23.04.2
ml Singularity/3.6.4
# ml Graphviz/2.38.0-foss-2016b
ml Graphviz/2.47.2-GCCcore-10.3.0

export NXF_SINGULARITY_CACHEDIR=/nemo/lab/ulej/home/shared/singularity
export NXF_HOME=/nemo/lab/ulej/home/users/luscomben/users/iosubi/.nextflow

nextflow pull iraiosub/riboseq -r dev-nf

nextflow run iraiosub/riboseq -r dev-nf \
-profile crick \
-resume \
--skip_umi_extract \
--with_umi \
--skip_trimming \
--umi_separator 'rbc:' \
--strandedness forward \
--org GRCh38 \
--input ./data/samplesheet.csv \
--outdir results_full_test_dev_nf
