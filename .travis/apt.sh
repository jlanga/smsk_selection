#!/usr/bin/env bash
apt update
apt install --yes \
    build-essential \
    git \
    curl \
    autotools-dev \
    autoconf \
    autoconf-archive \
    automake \
    cmake \
    libboost-all-dev \
    libtool \
    liblapack-dev \
    libatlas-cpp-0.6-dev \
    libnlopt-dev \
    libnlopt0 \
    libarmadillo-dev \
    libopenmpi-dev \
    libopenblas-dev \
    openmpi-bin \
    wget
