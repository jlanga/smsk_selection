#!/usr/bin/env bash

set -euo pipefail

MIN_BACKGROUND="2"
MIN_FOREGROUND="2"
OMEGA_ZEROS="0.5,1.0,1.5"
JOBS="1"
FASTCODEML_BINARY="fast"

print_help(){
    echo "run_fastcodeml_folder.sh"
    echo "-h|--help: print this message"
    echo "-t|--tree: input species tree (newick). Mandatory"
    echo "-m|--input-msa-folder: folder with all codon alignments in fasta format. Mandatory"
    echo "-w|--omega-zeros: starting omega values. Default= '0.5,1.0,1.5'"
    echo "-t|--target-species: species used for the foreground branch. Mandatory"
    echo "-b|--min-background: Minimum number of species to use as the background branch. Default: 2"
    echo "-f|--min-foreground: Minimum number of species to use as the foreground branch. Default: 2"
    echo "-o|--output-prefix: Output prefix of the files generated. Mandatory"
    echo "-p|--output-pvalues: File where to write the tsv file with all the p-values"
    echo "-c|--fastcodeml-binary: Path to the fastcodeml binrary. Default: 'fast'"
    echo "-j|--jobs: Number of jobs to run in parallel. Default: 1"
}

while (( "$#" )); do

    case "$1" in
        -t|--tree)
            TREE=$2
            shift 2
            ;;
        -m|--input-msa-folder)
            INPUT_MSA_FOLDER=$2
            shift 2
            ;;
        -w|--omega-zeros)
            OMEGA_ZEROS=$2
            shift 2
            ;;
        -s|--target-species)
            TARGET_SPECIES=$2
            shift 2
            ;;
        -b|--min-background)
            MIN_BACKGROUND=$2
            shift 2
            ;;
        -f|--min-foreground)
            MIN_FOREGROUND=$2
            shift 2
            ;;
        -o|--output-folder)
            OUTPUT_FOLDER=$2
            shift 2
            ;;
        -c|--fastcodeml-binary)
            FASTCODEML_BINARY=$2
            shift 2
            ;;
        -j|--jobs)
            JOBS=$2
            shift 2
            ;;
        -p|--output-pvalues)
            OUTPUT_PVALUES=$2
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

find "$INPUT_MSA_FOLDER" -type f -name "*.fa" \
| sort -V \
| parallel \
    --jobs "$JOBS" \
    --keep-order \
    python src/homologs/run_fastcodeml.py \
        --tree "$TREE" \
        --msa {} \
        --omega-zeros "$OMEGA_ZEROS" \
        --target-species "$TARGET_SPECIES" \
        --min-foreground "$MIN_FOREGROUND" \
        --min-background "$MIN_BACKGROUND" \
        --output-prefix "$OUTPUT_FOLDER/{/.}" \
        --fastcodeml-binary "$FASTCODEML_BINARY"

printf "prefix\tomega_zero\tlnl0\tlnl1\tl\tp_value\n" > "$OUTPUT_PVALUES"
find "$OUTPUT_FOLDER" -name "*.tsv" -type f \
| sort -V \
| xargs --max-lines=1 tail -n+2 \
| sort  -V -k 1,1 -k2,2 >> "$OUTPUT_PVALUES"
