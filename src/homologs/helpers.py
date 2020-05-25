#!/usr/bin/env python

"""Set of helper functions for the orthofinder/positive selection pipeline"""

import os
import multiprocessing as mp

from Bio.SeqIO.FastaIO import SimpleFastaParser

def fix_dir_path(path):
    """Add a / at the end if necessary"""
    return path if path.endswith("/") else path + "/"


def fasta_to_dict(filename):
    """Compose dict {seqname: sequence} from a fasta file"""
    with open(filename, "r") as handle:
        return {
            identifier.split()[0]: sequence
            for identifier, sequence in SimpleFastaParser(handle)
        }


def process_folders(in_dir, in_ext, out_dir, out_ext, process_function):
    """Apply function for every file in dir_in that has ext_in as extension,
    and store it in dir_out with extension ext_out. The same filenames are used.

    The function has to have the shape function(filein, fileout)
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for in_file in os.listdir(in_dir):
        if in_file.endswith(in_ext):
            prefix = in_file.split("/")[-1].split(in_ext)[0]
            out_file = out_dir + prefix + out_ext
            process_function(in_dir + in_file, out_file)


def process_folders_parallel(in_dir, in_ext, out_dir, out_ext, process_function, ncpus=1):
    """Apply a function in parallel over entire folders"""

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    prefixes = [
        in_file.split("/")[-1].split(in_ext)[0]
        for in_file in os.listdir(in_dir)
        if in_file.endswith(in_ext)
    ]
    file_ins = (in_dir + prefix + in_ext for prefix in prefixes)
    file_outs = (out_dir + prefix + out_ext for prefix in prefixes)

    pool = mp.Pool(ncpus)
    to_process = tuple(zip(file_ins, file_outs))
    pool.starmap(process_function, to_process)
    pool.close()
    pool.join()
