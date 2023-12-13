#!/usr/bin/sh
#SBATCH --job-name=riboseq-flow-test
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu

# conda activate env-nf
module purge
# ml Java/11.0.2
# ml Nextflow/22.10.3
ml Nextflow/21.10.3
# ml Nextflow/23.10.0
ml Singularity/3.6.4
# ml Graphviz/2.38.0-foss-2016b

export NXF_SINGULARITY_CACHEDIR=/nemo/lab/ulej/home/shared/singularity
export NXF_HOME=/nemo/lab/ulej/home/users/luscomben/users/iosubi/.nextflow

nextflow pull iraiosub/riboseq-flow -r v1.0.1

nextflow run iraiosub/riboseq-flow -r v1.0.1 \
-profile test,singularity