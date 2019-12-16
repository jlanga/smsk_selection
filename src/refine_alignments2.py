#!/usr/bin/env python3

"""
Script to perform MSA with t-coffee and associated tools to obtain high quality
alignments:

- align with t-coffe: muscle_msa, t_coffee_msa, mafftgins_msa, kalign_msa
- evaluate alignment
- keep colums with score = 9 (max)
- remove columns where there is a gap in at leas 50% of the times
- run maxalign
- return return the final msa in aa and nt formats
"""

import os
import sys
import multiprocessing as mp

from Bio import AlignIO
from Bio import SeqIO

def run_tcoffee_mcoffee(filename_in, filename_out):
    """Run t_coffee with the methods in `methods"""
    methods = ["muscle_msa", "mafftgins_msa", "t_coffee_msa", "kalign_msa"]
    methods_str = " ".join(methods)
    command = (
        "t_coffee "  + filename_in + " "
        "-method " + methods_str + " "
        "-output=aln "
        "-outfile " + filename_out + " "
        "-quiet "
        "2>&1"
    )
    sys.stderr.write(command + "\n")
    os.system(command)


def run_tcoffee_eval(filename_in, filename_out):
    """Run t_coffee in eval mode"""
    command = (
        "t_coffee "
        "-infile " + filename_in + " "
        "-evaluate "
        "-output=score_ascii "
        "-outfile " + filename_out + " "
        "-quiet "
        "2>&1"
    )
    sys.stderr.write(command + "\n")
    os.system(command)


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


def run_max_align(filename_in, filename_out):
    """Execute maxalign"""
    maxalign = "perl src/maxalign.pl"
    command = maxalign + " -p " + filename_in + " > " + filename_out
    sys.stderr.write(command + "\n")
    os.system(command)


def subset_cds(msa_in, fasta_out, cds_dict):
    """Slice cds_dict with the sequences in msa_in"""
    msa = AlignIO.read(msa_in, format="fasta")
    seq_names = [msa[i].name for i in range(len(msa))]
    cds_sequences = {name: cds_dict[name] for name in seq_names}
    with open(fasta_out, "w") as f_in:
        for name, sequence in cds_sequences.items():
            f_in.write(">{name}\n{sequence}\n".format(
                name=name,
                sequence=sequence
            ))


def run_max_align_cds(filename_in_pep, filename_in_cds, filename_out):
    """Execute maxalign, returning the cds"""
    maxalign = "perl src/maxalign.pl"
    command = maxalign + " -p " + filename_in_pep + " " + filename_in_cds + " > " + filename_out
    sys.stderr.write(command + "\n")
    os.system(command)


def run_pipeline(filename_in, cds_dict):
    """Run the refinement step per file"""
    orthogroup_id = filename_in.split("/")[-1].split(".")[0]
    output_folder = "/".join(filename_in.split("/")[0:-1])

    run_tcoffee_mcoffee(
        filename_in=filename_in,
        filename_out=output_folder + "/" + orthogroup_id + ".aln"
    )

    run_tcoffee_eval(
        filename_in=output_folder + "/" + orthogroup_id + ".aln",
        filename_out=output_folder + "/" + orthogroup_id + ".cons"
    )

    scores = parse_scores(filename=output_folder + "/" + orthogroup_id + ".cons")
    hq_scores = get_high_score_columns(scores, 9)
    msa = AlignIO.read(
        handle=output_folder + "/" + orthogroup_id + ".aln",
        format="clustal"
    )
    msa_hq = keep_columns(msa, hq_scores)

    highly_occupied_columns = get_highly_occupied_columns(msa_hq, 0.5)
    msa_vhq = keep_columns(msa_hq, highly_occupied_columns)

    AlignIO.write(
        alignments=msa_vhq,
        handle=output_folder + "/" + orthogroup_id + ".tcoffee.fa",
        format="fasta"
    )

    run_max_align(
        filename_in=output_folder + "/" + orthogroup_id + ".tcoffee.fa",
        filename_out=output_folder + "/" + orthogroup_id + ".maxalign.fa",
    )

    subset_cds(
        msa_in=output_folder + "/" + orthogroup_id + ".maxalign.fa",
        fasta_out=output_folder + "/" + orthogroup_id + ".cds",
        cds_dict=cds_dict
    )

    run_max_align_cds(
        filename_in_pep=output_folder + "/" + orthogroup_id + ".tcoffee.fa",
        filename_in_cds=output_folder + "/" + orthogroup_id + ".cds",
        filename_out=output_folder + "/" + orthogroup_id + ".maxalign.cds",
    )

    #os.remove(orthogroup_id + ".dnd")


def run_pipeline_starmap(params):
    """Wrapper of run_pipeline to run in parallel"""
    filename_in, cds_dict = params[0:2]
    run_pipeline(filename_in, cds_dict)


if __name__ == '__main__':

    if len(sys.argv) != 5:
        sys.stderr.write(
            "ERROR: Incorrect number of parameters: \n" + \
            "python refine_alignments.py indir in_extension all.cds cores\n"
        )
        exit(1)

    IN_DIR = sys.argv[1]
    IN_EXT = sys.argv[2]
    CDS = sys.argv[3]
    CORES = int(sys.argv[4])

    if not IN_DIR.endswith("/"):
        IN_DIR += "/"

    CDS_DICT = {
        identifier.split()[0]: sequence
        for identifier, sequence in SeqIO.FastaIO.SimpleFastaParser(open(CDS, "r"))
    }

    IN_FILES = tuple(
        IN_DIR + in_file
        for in_file in os.listdir(IN_DIR)
        if in_file.endswith(IN_EXT)
    )

    PARAMS = zip(
        IN_FILES,
        (CDS_DICT for i in IN_FILES)
    )
    #print(PARAMS)
    POOL = mp.Pool(CORES)
    POOL.map(run_pipeline_starmap, PARAMS)
