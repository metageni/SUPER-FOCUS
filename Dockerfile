# conda setup
FROM continuumio/miniconda3
RUN conda create -n env python=3.6
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

# creating conda environment
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# installs super-focus using bioconda
RUN conda install super-focus

# Entrypoint
CMD /bin/bash
