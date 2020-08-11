#!/usr/bin/env python3

import argparse
import multiprocessing as mp
import tempfile
import os
from subprocess import run


from helpers import \
    fasta_to_dict, \
    fix_dir_path, \
    process_folders_parallel


bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases] + ['---'] 
aminoacids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-'
codon_table = dict(zip(codons, aminoacids))
del bases, codons, aminoacids

def translate(string):
    """Translate a nucleotide string to amino acid"""
    n = len(string)
    string = string.upper()
    if n % 3 != 0:
        sys.exit(f"ERROR: Incorrect number of characters: {n}")
    translated = ""
    for i in range(0, n, 3):
        codon = string[i:i+3]
        if codon not in codon_table:
            translated += "-"
        else:
            translated += codon_table[codon]
    return translated

def translate_fasta(fasta_in, fasta_out):
    """Translate an entire fasta aligment"""
    with open(fasta_out, "w") as f_out:
        for name, sequence in fasta_to_dict(fasta_in).items():
            f_out.write(f">{name}\n{translate(sequence)}\n")


def run_trimal(fasta_in, fasta_out, fasta_gapless):
    """Run trimal.
    
    At least 60% of the col is not a gap, simmilariy score < 0.001
    """
    command = [
        "trimal",
        "-in", fasta_in,
        "-out", fasta_out,
        "-gt", "0.6",
        "-st", "0.001",
        "-backtrans", fasta_gapless
    ]
    run(command)


def remove_gaps(fasta_in, fasta_out):
    """Remove all gaps in the alignment"""
    with open(fasta_out, "w") as f_out:
        for name, sequence in fasta_to_dict(fasta_in).items():
            f_out.write(f">{name}\n{sequence.replace('-', '')}\n")


def run_pipeline(raw_fn, trimmed_fn):
    """
    Align CDS with trimal (translate | trim)
    """
    translated = tempfile.NamedTemporaryFile()
    gapless = tempfile.NamedTemporaryFile()

    translate_fasta(fasta_in=raw_fn, fasta_out=translated.name)
    remove_gaps(fasta_in=raw_fn, fasta_out=gapless.name)
    run_trimal(
        fasta_in=translated.name,
        fasta_out=trimmed_fn,
        fasta_gapless=gapless.name
    )
    
    


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run trimal over an entire folder of codon alignments"
    )
    parser.add_argument(
        '-i', '--input-folder',
        help="Input folder. All files must be *.fa",
        required=True
    )
    parser.add_argument(
        '-o', '--output-folder',
        help="Output folder. All files will be *.fa",
        required=True
    )
    parser.add_argument(
        '-t', '--threads',
        help="Number of threads to use",
        default=1,
        type=int,
        required=False
    )
    return  vars(parser.parse_args())



if __name__ == '__main__':

    ARGS = parse_arguments()

    process_folders_parallel(
        fix_dir_path(ARGS["input_folder"]), "fa",
        fix_dir_path(ARGS["output_folder"]), "fa",
        run_pipeline,
        ncpus=ARGS["threads"]
    )
