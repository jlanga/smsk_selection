#!/usr/bin/env bash
set -euxo pipefail


# Install miniconda

url="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
if [[ -d "$HOME"/miniconda3/bin ]]; then
    echo "miniconda already installed."
else
    echo "Installing miniconda."
    mkdir -p "$HOME"/download
    wget \
        --continue \
        --output-document "$HOME"/download/miniconda3.sh \
        $url
    chmod +x "$HOME"/download/miniconda3.sh
    "$HOME"/download/miniconda3.sh \
        -u \
        -b \
        -p "$HOME"/miniconda3
    "$HOME"/miniconda3/bin/conda clean --all --yes
    echo "local_repodata_ttl: 1800" >> ~/.condarc
fi
