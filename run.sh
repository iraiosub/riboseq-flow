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

# export NXF_SINGULARITY_CACHEDIR=/camp/lab/luscomben/home/users/iosubi/nfcore/riboseq
export  NXF_CONDA_CACHEDIR=/camp/lab/ulej/home/users/luscomben/users/iosubi/nfcore/riboseq
# export NXF_HOME=/camp/lab/luscomben/home/users/iosubi/.nextflow # for Nextflow versions > 22.x

nextflow pull ulelab/riboseq -r feat-ribocutter

nextflow run ulelab/riboseq -r feat-ribocutter \
-profile conda,crick,test \
-resume \
--min_read_length 23 \
--extra_ribocutter_args '--max_read_length 60'
