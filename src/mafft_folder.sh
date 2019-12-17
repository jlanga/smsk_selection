#!/usr/bin/env bash
set -euo pipefail


in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5


mkdir -p "$out_dir"

#for filein in "$in_dir"/*."$in_ext" ; do

#    cluster_id=$(basename -- "$filein" ."$in_ext")

#    fileout="$out_dir/${cluster_id}.$out_ext"

#    if [[ -f "$fileout" ]] ; then

#        >&2 echo "$fileout exits. Skipping."

#    else

#        mafft --amino --auto --maxiterate 1000 --thread "$threads" "$filein" > "$fileout"

#    fi

#done

parallel \
    --jobs "$threads" \
    mafft --amino --genafpair --maxiterate 1000 {} ">" "$out_dir"/{/.}."$out_ext" \
::: "$in_dir"/*."$in_ext"
