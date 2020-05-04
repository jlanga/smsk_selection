#!/usr/bin/env bash
set -euo pipefail

docker build -t smsk_orthofinder .

docker run \
    --interactive \
    --tty \
    --volume "$(pwd):$(pwd)" \
    --user "$(id -u)":"$(id -g)" \
    --workdir "$(pwd)" \
    --name "$(date +%Y%m%d-%H%M%S)" \
    smsk_orthofinder \
    bash -c 'export PATH=/opt/miniconda3/bin:$PATH; snakemake --use-conda --keep-going --keep-incomplete --show-failed-logs -j'
