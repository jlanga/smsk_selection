#!/usr/bin/env bash
set -euxo pipefail

echo "$PATH"
python --version
export PATH="$HOME"/miniconda3/bin:"$PATH"
snakemake --use-conda --jobs 4 --show-failed-logs
