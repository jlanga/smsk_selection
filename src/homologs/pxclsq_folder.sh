#!/usr/bin/env bash
set -euo pipefail

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
min_occupancy=$5


mkdir -p "$out_dir"

for filein in "$in_dir"/*."$in_ext" ; do

    cluster_id=$(basename -- "$filein" ."$in_ext")

    fileout="$out_dir/${cluster_id}.$out_ext"

    if [[ -f "$fileout" ]] ; then

        >&2 echo "$fileout exits. Skipping."
    
    else

        pxclsq \
            --aminoacid \
            --prop "$min_occupancy" \
            --seqf "$filein" \
            --outf "$fileout"-pht \
            --verbose \
        1>&2
    
        python2 src/homologs/trim_alignment_by_lengths.py \
            "$fileout"-pht \
            10 \
        > "$fileout"
    
        rm "$fileout"-pht
        
    fi

done


rm phyx.logfile