FROM ubuntu:16.04

LABEL version="1.5"

# Update distro and install basic tools
ENV APT_PACKAGES \
    build-essential \
    libopenmpi-dev \
    openmpi-bin \
    wget

RUN apt update && \
    apt install --yes ${APT_PACKAGES} && \
    rm -rf /var/lib/apt/lists/*


# Install exabayes
RUN \
    wget https://cme.h-its.org/exelixis/resource/download/software/exabayes-1.5.tar.gz && \
    tar xvf exabayes-1.5.tar.gz && \
    cd exabayes-1.5/ && \
    ./configure --enable-mpi && \
    make -j 4 && \
    make install
