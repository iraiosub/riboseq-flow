FROM nfcore/base:1.12.1
MAINTAINER Ira Iosub <ira.iosub@gmail.com>
LABEL authors="Ira Iosub" \
      description="Docker image containing software requirements for mapping length analysis of riboseq pipeline"

# Install the conda environment
COPY env.yml /
RUN conda env create --quiet -f /env.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/mapping-length/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name mapping-length > mapping-length.yml
