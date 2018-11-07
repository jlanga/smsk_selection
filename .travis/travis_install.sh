#!/usr/bin/env bash
set -euo pipefail

export PATH="$HOME"/miniconda3_"$TRAVIS_OS_NAME"/bin:"$PATH"
echo "$PATH"

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda install --yes --channel bioconda snakemake
conda clean --all --yes
