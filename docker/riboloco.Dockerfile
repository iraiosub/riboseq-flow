FROM nfcore/base:1.12.1
MAINTAINER Ira Iosub <ira.iosub@gmail.com>
LABEL authors="Ira Iosub" \
      description="Docker image containing software requirements for riboloco"

# Install the conda environment
COPY riboloco.yml /
RUN conda env create --quiet -f /riboloco.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/riboloco_env/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name riboloco_env > riboloco_env.yml
