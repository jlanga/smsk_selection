#!/usr/bin/env python3

"""
Script to perform MSA with t-coffee and associated tools to obtain high quality
alignments:

- align with t-coffe: t_coffee_msa mafftgins_msa
"""

import os
import sys

from Bio import AlignIO

def run_tcoffee_mcoffee(filename_in, filename_out):
    """Run t_coffee with the methods in `methods"""
    methods = ["muscle_msa", "mafftgins_msa", "t_coffee_msa", "kalign_msa"]
    methods_str = " ".join(methods)
    command = f"t_coffee " + \
        f"{filename_in} " + \
        f"-method {methods_str} " + \
        f"-output=aln " + \
        f"-outfile {filename_out} " + \
        f"2>&1"
    sys.stderr.write(command)
    os.system(command)


def run_tcoffee_eval(filename_in, filename_out):
    """Run t_coffee in eval mode"""
    command = f"t_coffee " + \
        f"-infile {filename_in} " + \
        "-evaluate " + \
        "-output=score_ascii " + \
        f"-outfile {filename_out} " + \
        "2>&1"
    sys.stderr.write(command)
    os.system(command)


def parse_scores(filename):
    """
    Parse the scores from the ascii_results from t_coffee -evaluate
    """
    scores = []
    with open(filename, "r") as f_in:
        for line in f_in.readlines():
            if line.startswith("cons") and ":" not in line:
                line = line[4:].strip()
                line.replace(" ", "")
                scores += [int(x) for x in line]
    return scores


def keep_columns(msa, columns_to_keep):
    """Trim the msa, preserving the columns to keep"""
    msa_trimmed = msa[:, :0]
    for column_number in columns_to_keep:
        msa_trimmed += msa[:, column_number:column_number+1]
    return msa_trimmed


def get_gappy_columns(msa, threshold=0.5):
    """Get the column numbers that contain gaps in at least `threshold`
    proportion"""
    columns_to_remove = []
    number_of_sequences = len(msa)
    for column_number in range(msa):
        column = msa[:, column_number]
        number_of_gaps = column.count("-")
        if number_of_gaps / number_of_sequences >= threshold:
            columns_to_remove.append(column_number)
    return columns_to_remove


def get_high_score_columns(scores, min_score):
    """From the array of scores, return the indexes of those above min_score"""
    high_score_columns = list()
    for i, score in enumerate(scores):
        if score >= min_score:
            high_score_columns.append(i)
    return high_score_columns


def get_highly_occupied_columns(msa, min_ratio=0.5):
    """Return the columns that are occupied at least min_ratio"""
    highly_occupied_columns = []
    number_of_sequences = len(msa)
    number_of_columns = msa.get_alignment_length()
    for column_number in range(number_of_columns):
        column_str = msa[:, column_number]
        gap_ratio = column_str.count("-") / number_of_sequences
        if gap_ratio <= min_ratio:
            highly_occupied_columns.append(column_number)
    return highly_occupied_columns


def run_max_align(filename_in, filename_out):
    """Execute maxalign"""
    maxalign = "perl src/maxalign.pl"
    command = f"{maxalign} -p {filename_in} > {filename_out}"
    sys.stderr.write(command)
    os.system(command)


def run_pipeline(filename_in):
    """Run the refinement step per file"""
    orthogroup_id = filename_in.split("/")[-1].split(".")[0]
    output_folder = "/".join(filename_in.split("/")[0:-1])

    run_tcoffee_mcoffee(
        filename_in=filename_in,
        filename_out=f"{output_folder}/{orthogroup_id}.aln"
    )

    run_tcoffee_eval(
        filename_in=f"{output_folder}/{orthogroup_id}.aln",
        filename_out=f"{output_folder}/{orthogroup_id}.cons"
    )

    scores = parse_scores(filename=f"{output_folder}/{orthogroup_id}.cons")
    hq_scores = get_high_score_columns(scores, 9)
    msa = AlignIO.read(
        handle=f"{output_folder}/{orthogroup_id}.aln",
        format="clustal"
    )
    msa_hq = keep_columns(msa, hq_scores)
    highly_occupied_columns = get_highly_occupied_columns(msa_hq, 0.5)
    msa_vhq = keep_columns(msa_hq, highly_occupied_columns)
    AlignIO.write(
        alignments=msa_vhq,
        handle=f"{output_folder}/{orthogroup_id}.tcoffee.fa",
        format="fasta"
    )
    run_max_align(
        filename_in=f"{output_folder}/{orthogroup_id}.tcoffee.fa",
        filename_out=f"{output_folder}/{orthogroup_id}.maxalign.fa"
    )


if __name__ == '__main__':
    run_pipeline(sys.argv[1])
