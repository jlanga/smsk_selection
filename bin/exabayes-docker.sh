#!/usr/bin/env bash
set -euo pipefail

docker run \
    --rm \
    --volume "$PWD":"$PWD" \
    --user "$UID":"$GID" \
    exabayes:1.5  \
    exabayes "$@"