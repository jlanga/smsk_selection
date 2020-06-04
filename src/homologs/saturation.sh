#!/usr/bin/env bash

set -euo pipefail


ete3_m0(){

    # Run ete3 evol to compute only M0
    # From it, get the average omega for all tree

    tree=$1; shift
    msa=$1; shift
    output_dir=$1; shift
    output_txt=$1; shift

    ete3 evol \
        --resume \
        -t "$tree" \
        --alg "$msa" \
        -o "$output_dir" \
        --cpu 1 \
        --resume \
        --codeml_binary "$PWD/bin/codeml" \
        --model M0 \
    > "$output_txt"
    
}


ete3_bfree(){

    # Run ete3 evol to compute the aerage 

    tree=$1; shift
    msa=$1; shift
    output_dir=$1; shift
    output_txt=$1; shift
    species=$1; shift

    transcripts=$(python src/homologs/get_definers.py "$tree" "$species")

    ete3 evol \
        --resume \
        -t "$tree" \
        --alg "$msa" \
        -o "$output_dir" \
        --mark "$transcripts" \
        --cpu 1 \
        --resume \
        --codeml_binary "$PWD/bin/codeml" \
        --models b_free \
    > "$output_txt" 
}


run_models(){

    species_tree=$1; shift
    msa=$1; shift
    output_dir=$1; shift
    species=$1; shift

    mkdir --parents "$output_dir"

    orthogroup=$(basename "$msa" .fa)

    transcript_tree="$output_dir/$orthogroup.nwk" 

    python src/homologs/shape_tree_as_msa.py \
        --tree "$species_tree" \
        --msa "$msa" \
        --output "$transcript_tree"

    m0out="$output_dir/$orthogroup.m0.txt"
    bfreeout="$output_dir/$orthogroup.bfree.txt"

    ete3_m0 "$transcript_tree" "$msa" "$output_dir" "$m0out"
    ete3_bfree "$transcript_tree" "$msa" "$output_dir" "$bfreeout" "$species"

}


get_m0(){
    m0_out=$1; shift

    grep ^"tree length for dS" "$m0_out" \
    | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'
}


get_bfree(){
    bfree_out=$1; shift
    target_species=$1; shift

    tempfile=$(mktemp)

    tree=$(grep ^"dS tree:" -A 1 "$bfree_out" | tail -n+2)
    echo "$tree" > "$tempfile"
    
    branch_length=$(python src/homologs/extract_branch_length.py \
        --tree "$tempfile" \
        --species "$target_species")
    
    # rm -f "$tempfile"
    echo "$branch_length"

}



JOBS="1"


print_help(){
    echo "run_fastcodeml_folder.sh"
    echo "-h|--help: print this message"
    echo "-t|--tree: input species tree (newick). Mandatory"
    echo "-m|--input-msa-folder: folder with all codon alignments in fasta format. Mandatory"
    echo "-t|--target-species: species used for the foreground branch. Mandatory"
    echo "-o|--output: File with the results. Mandatory"
    echo "-j|--jobs: Number of jobs to run in parallel. Default: 1"
}

while (( "$#" )); do

    case "$1" in
        -t|--input-tree)
            TREE=$2
            shift 2
            ;;
        -m|--input-msa-folder)
            MSA_FOLDER=$2
            shift 2
            ;;
        -s|--target-species)
            TARGET_SPECIES=$2
            shift 2
            ;;
        -o|--output-file)
            OUTPUT_FILE=$2
            shift 2
            ;;
        -j|--jobs)
            JOBS=$2
            shift 2
            ;;
        -d|--output-folder)
            OUTPUT_FOLDER=$2
            shift 2
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "ERROR: Unupported flag $1" >&2
            print_help >&2
            exit 1
            ;;
    esac
done

mkdir --parents "$OUTPUT_FOLDER"

export -f ete3_m0
export -f ete3_bfree
export -f run_models

find "$MSA_FOLDER" -type f -name "*.fa" \
| sort --version-sort \
| parallel \
    --jobs "$JOBS" \
    run_models \
        "$TREE" \
        "$MSA_FOLDER/{/.}.fa" \
        "$OUTPUT_FOLDER/{/.}" \
        "$TARGET_SPECIES"



echo -e "orthogroup\tds_m0\tds_bfree" > "$OUTPUT_FILE"

for orthogroup in "$OUTPUT_FOLDER"/* ; do

    orthogroup=$(basename "$orthogroup")

    ds_m0=$(get_m0 "$OUTPUT_FOLDER/$orthogroup"/M0*/out)
    ds_bfree=$(get_bfree "$OUTPUT_FOLDER/$orthogroup"/b_free*/out "$TARGET_SPECIES")

    echo -e "$orthogroup\t$ds_m0\t$ds_bfree" >> "$OUTPUT_FILE"

done
