#!/usr/bin/env python

"""parse_ete3_evol_folder.py: parse the output file of ete3 evol to get the 
models and p-values
"""

import sys
import os

from helpers import \
    fix_dir_path, \
    process_folders


def parse_ete3_evol_file(file_in):
    """Parse the output from a single output file from ete3 evol"""
    orthogroup = ".".join(file_in.split("/")[-1].split(".")[0:-3])
    omega = ".".join(file_in.split("/")[-1].split(".")[-3:-1])
    with open(file_in, "r") as f_in:
        lines = [line.strip() for line in f_in.readlines()]
    for numberline, line in enumerate(lines):
        if line[0:3] == "LRT":
            lrt_line = numberline
        if line[0:17] == "SUMMARY BY MODEL":
            summary_line = numberline
    block = lines[lrt_line + 4: summary_line -1]
    assert(block)
    for line in block:
        model1, model2, pvalue = [x.strip() for x in line.split("|")]
        model1 = model1.split(".")[0]
        model2 = model2.split(".")[0]
        pvalue = pvalue.strip("*")
        sys.stdout.write(
            f"{orthogroup}\t{omega}\t{model1}\t{model2}\t{pvalue}\n"
        )
    

def parse_ete3_evol_folder(in_dir, in_ext):
    """Parse an entire folder"""
    sys.stdout.write("orthogroup\tstarting_omega\tmodel1\tmodel2\tp-value\n")
    for file in sorted(os.listdir(in_dir)):
        if file.endswith(in_ext):
            parse_ete3_evol_file(in_dir + "/" + file)


if __name__ == '__main__':

    if len(sys.argv) != 3:
        sys.stderr.write(
            "ERROR! Incorrect number of files\n"
            "python parse_ete3_evol_folder.py folder_in extension_in\n"
        )
        sys.exit(1)

    IN_DIR = fix_dir_path(sys.argv[1])
    IN_EXT = sys.argv[2]

    parse_ete3_evol_folder(IN_DIR, IN_EXT)
