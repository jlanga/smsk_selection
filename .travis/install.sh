#!/usr/bin/env bash
set -euxo pipefail

# conda
export PATH="$HOME"/miniconda3/bin:"$PATH"
echo "$PATH"
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install --quiet --yes --channel bioconda snakemake=5.15
snakemake --use-conda --conda-create-envs-only -j 4
conda clean --all --yes

# Non-conda software
mkdir -p bin/

# phyx
pushd src/phyx-1.01/src
./configure
make -j
cp px* ../../../bin/
popd

# guidance
pushd src/guidance.v2.02/
make -j
popd

# fastcodeml
pushd src/fastcodeml
# MATH_LIB_NAMES="openblas;lapack" cmake . -DUSE_LAPACK:BOOL=ON
cmake .
make -j
cp fast ../../bin/
popd
