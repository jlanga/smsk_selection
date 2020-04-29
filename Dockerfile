FROM ubuntu:18.04

SHELL ["/bin/bash", "--login", "-c"]

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# apt packages
RUN apt update && apt install --yes \
    autoconf \
    autoconf-archive \
    automake \
    autotools-dev \
    build-essential \
    curl \
    git \
    cmake \
    libarmadillo-dev \
    libatlas-cpp-0.6-dev \
    libboost-all-dev \
    liblapack-dev \
    libnlopt0 \
    libnlopt-dev \
    libopenmpi-dev \
    libopenblas-dev \
    libtool \
    openmpi-bin \
    wget \
&& rm -rf /var/lib/apt/lists/*

# RUN bash .travis/before_install.sh

RUN wget --continue https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -u -b -p "$HOME"/miniconda3
RUN conda clean --all --yes


RUN \
    conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install --yes --channel bioconda snakemake=5.15

# ADD . /root/smsk_orthofinder/

# WORKDIR /root/smsk_orthofinder/

# Non-conda software
# RUN mkdir -p bin/

# phyx
# RUN \
#     pushd src/phyx-1.01/src && \
#     ./configure && \
#     make -j 4 && \
#     cp px* ../../../bin/ && \
#     popd

# # guidance
# RUN \
#     pushd src/guidance.v2.02/ && \
#     make -j 4 && \
#     popd

# # fastcodeml
# RUN \
#     pushd src/fastcodeml && \
#     `# MATH_LIB_NAMES="openblas;lapack" cmake . -DUSE_LAPACK:BOOL=ON` \
#     cmake . && \
#     make -j 4 && \
#     cp fast ../../bin/ && \
#     popd

# RUN snakemake --use-conda --jobs 4 --show-failed-logs busco
# RUN snakemake --use-conda --jobs 4 --show-failed-logs cdhit
# RUN snakemake --use-conda --jobs 4 --show-failed-logs orthofinder
# RUN snakemake --use-conda --jobs 4 --show-failed-logs homologs
# RUN snakemake --use-conda --jobs 4 --show-failed-logs tree
# RUN snakemake --use-conda --jobs 4 --show-failed-logs selection_guidance


