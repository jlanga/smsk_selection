#!/usr/bin/env python2.7

import os
import sys
import multiprocessing as mp

from Bio import AlignIO
from Bio import SeqIO


def parse_scores(filename):
    """Parse the scores from the ascii_results from t_coffee -evaluate"""
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
    """
    Get the column numbers that contain gaps in at least `threshold` proportion
    """
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


def filter_tcoffee_files(aln_in, cons_in, fasta_out):

    scores = parse_scores(filename=cons_in)
    hq_scores = get_high_score_columns(scores, 9)
    msa = AlignIO.read(
        handle=aln_in,
        format="clustal"
    )
    msa_hq = keep_columns(msa, hq_scores)

    highly_occupied_columns = get_highly_occupied_columns(msa_hq, 0.5)
    msa_vhq = keep_columns(msa_hq, highly_occupied_columns)

    AlignIO.write(
        alignments=msa_vhq,
        handle=fasta_out,
        format="fasta"
    )


if __name__ == '__main__':

    if len(sys.argv) != 4:
        sys.stderr.write(
            "ERROR: Incorrect number of parameters: \n" + \
            "python filter_tcoffee.py tcoffee_aln tcoffee_cons out_aln\n"
        )
        sys.exit(1)
    
    ALN_IN = sys.argv[1]
    CONS_IN = sys.argv[2]
    ALN_OUT = sys.argv[3]

    filter_tcoffee_files(ALN_IN, CONS_IN, ALN_OUT)