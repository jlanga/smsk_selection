#!/usr/bin/env bash
set -euxo pipefail

# conda
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install --quiet --yes --channel conda-forge mamba
mamba install --yes --channel bioconda snakemake
snakemake --use-conda --conda-create-envs-only -j 4
conda clean --all --yes

# Non-conda software
mkdir -p bin/

# phyx
pushd src/phyx-1.01/src
./configure
make -j 4
cp px* ../../../bin/
popd

# guidance
pushd src/guidance.v2.02/
make -j 4
popd

# fastcodeml
pushd src/fastcodeml
# MATH_LIB_NAMES="openblas;lapack" cmake . -DUSE_LAPACK:BOOL=ON
cmake .
make -j 4
cp fast ../../bin/
popd