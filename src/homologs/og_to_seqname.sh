#!/usr/bin/env bash
set -euo pipefail


process_orthogroup(){

    fullname=$1
    filename=$(basename -- "$fullname")
    filename="${filename%.*}"

    for seqname in $(grep ^">" "$fullname" | tr -d ">"); do
        echo -e "$filename\t$seqname"
    done

}

export -f process_orthogroup    

find "$1" -name "*.fa" -type f \
| sort --version-sort \
| xargs -P 1 -I "{}" bash -c "process_orthogroup {}"
