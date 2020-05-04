FROM ubuntu:18.04

SHELL ["/bin/bash", "--login", "-c"]

# apt packages
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update \
&& DEBIAN_FRONTEND="noninteractive" apt install --yes \
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
    tzdata \
    wget \
    xvfb \
&& rm -rf /var/lib/apt/lists/*

ENV miniconda=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/miniconda3/bin:$PATH"

RUN \
    wget --quiet --continue $miniconda && \
    bash Miniconda3-latest-Linux-x86_64.sh -bfp /opt/miniconda3 && \
    conda clean --all --yes && \
    rm Miniconda3-latest-Linux-x86_64.sh

RUN \
    conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install --yes --channel bioconda snakemake-minimal=5.15 pandas && \
    conda clean --all --yes

RUN \
    mkdir /.conda \
    && chmod ugo+rwx /.conda
