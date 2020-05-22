FROM ubuntu:18.04

# apt packages

ENV DEBIAN_FRONTEND="noninteractive" 

RUN apt-get update \
&& apt-get install --yes --no-install-recommends \
    autoconf=2.69-11 \
    autoconf-archive=20170928-2 \
    automake=1:1.15.1-3ubuntu2 \
    autotools-dev=20180224.1 \
    build-essential=12.4ubuntu1 \
    ca-certificates=20180409 \
    cmake=3.10.2-1ubuntu2.18.04.1 \
    curl=7.58.0-2ubuntu3.8 \
    git=1:2.17.1-1ubuntu0.7 \
    libarmadillo-dev=1:8.400.0+dfsg-2 \
    libatlas-cpp-0.6-dev=0.6.3-4ubuntu1 \
    libboost-all-dev=1.65.1.0ubuntu1 \
    liblapack-dev=3.7.1-4ubuntu1 \
    libnlopt0=2.4.2+dfsg-4 \
    libnlopt-dev=2.4.2+dfsg-4 \
    libopenmpi-dev=2.1.1-8 \
    libopenblas-dev=0.2.20+ds-4 \
    libtool=2.4.6-2 \
    openmpi-bin=2.1.1-8 \
    tzdata=2019c-0ubuntu0.18.04 \
    wget=1.19.4-1ubuntu2.2 \
    xvfb=2:1.19.6-1ubuntu4.4 \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

SHELL ["/bin/bash", "-c"]

ENV miniconda=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/miniconda3/bin:$PATH"

RUN wget --no-check-certificate --quiet --continue $miniconda \
&& bash Miniconda3-latest-Linux-x86_64.sh -bfp /opt/miniconda3 \
&& conda clean --all --yes \
&& rm Miniconda3-latest-Linux-x86_64.sh

RUN conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& conda install \
    --yes --channel bioconda \
    snakemake-minimal=5.15 \
    pandas=1.0.3 \
&& conda clean --all --yes

RUN mkdir /.conda \
&& chmod ugo+rwx /.conda
