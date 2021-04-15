#!/usr/bin/env bash
set -euxo pipefail

docker build -t smsk_orthofinder .

command="export PATH=/opt/miniconda3/bin:$PATH; "
command+="snakemake --use-conda --keep-going --keep-incomplete --show-failed-logs "
command+="$*"
docker run \
    --interactive \
    --tty \
    --volume "$(pwd):$(pwd)" \
    --user "$(id -u)":"$(id -g)" \
    --workdir "$(pwd)" \
    --name "$(date +%Y%m%d-%H%M%S)" \
    smsk_orthofinder \
    bash -c "$command"
