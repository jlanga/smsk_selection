#!/usr/bin/env bash
set -euo pipefail

export PATH="$HOME"/miniconda3_"$TRAVIS_OS_NAME"/bin:"$PATH"
echo "$PATH"

echo "local_repodata_ttl: 1800" >> ~/.condarc

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda install --yes --channel bioconda snakemake
conda clean --all --yes

snakemake --use-conda
