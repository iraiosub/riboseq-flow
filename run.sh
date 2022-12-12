#!/usr/bin/sh

ml purge
ml Nextflow/21.10.3
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

cd /camp/project/proj-luscombe-ule/working/ira-jure/RiboSeq_NextFlow

nextflow run main.nf