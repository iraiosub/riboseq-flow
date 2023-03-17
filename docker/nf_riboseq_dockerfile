FROM nfcore/base:1.12.1
MAINTAINER Ira Iosub <ira.iosub@gmail.com>
LABEL authors="Ira Iosub" \
      description="Docker image containing software requirements for the riboseq pipeline"

# Install the conda environment
COPY env2.yml /
RUN conda env create --quiet -f /env2.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-riboseq/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-riboseq > nf-riboseq.yml
