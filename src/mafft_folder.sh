#!/usr/bin/env bash
set -euo pipefail


in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5


mkdir -p "$out_dir"

for filein in "$in_dir"/*."$in_ext" ; do

    cluster_id=$(basename -- "$filein" ."$in_ext")

    fileout="$out_dir/${cluster_id}.$out_ext"

    if [[ -f "$fileout" ]] ; then

        >&2 echo "$fileout exits. Skipping."
    
    else

        mafft --amino --auto --thread "$threads" "$filein" > "$fileout"

    fi

done
