#!/usr/bin/env bash
set -euo pipefail

docker run \
    --interactive \
    --tty \
    --volume "$(pwd):$(pwd)" \
    `#--user "$(id -u)":"$(id -g)"` \
    --workdir "$(pwd)" \
    smsk_orthofinder snakemake --use-conda -j
