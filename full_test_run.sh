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

nextflow pull ulelab/riboseq -r dev

nextflow run ulelab/riboseq -r dev \
-profile conda,crick \
-resume \
--org GRCh38 \
--input $GITHUBDIR/data/samplesheet.csv \
--outdir $GITHUBDIR/results_full_test