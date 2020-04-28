#!/usr/bin/env bash
set -euxo pipefail

snakemake --cleanup-conda
conda clean --all --yes
rm -rf results/
rm -rf "$CONDA_PATH"/{locks,pkgs,var,conda-meta/history}
