FROM nfcore/base:1.12.1
MAINTAINER Ira Iosub <ira.iosub@gmail.com>
LABEL authors="Ira Iosub" \
      description="Docker image containing software requirements for riboloco"

# Install the conda environment
COPY analyse_riboloco.yml /
RUN conda env create --quiet -f /analyse_riboloco.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/analyse_riboloco/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name analyse_riboloco > analyse_riboloco.yml