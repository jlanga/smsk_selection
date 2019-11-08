#!/usr/bin/env python3
"""
extract_orthologs.py orthologs.tsv indir outdir in_extension out_extension
Process all fasta files in indir that have extension .in_extension and group
them clustered as in orthologs.tsv, writting them in outdir with extension
out_extension
seqtk necessary
input fasta files have to be named species.in_extension
"""

import os
import sys

import pandas as pd





def prepare_files_and_folders(orthogroup_ids, outdir, out_extension):
    """prepare_files_and_folders: helper function
    Create out folder and create/flush output files
    """
    os.system("mkdir -p {outdir}".format(outdir=outdir))
    for orthogroup_id in orthogroup_ids:
        os.system("cat /dev/null > {outdir}/{orthogroup_id}.{extension}".format(
            outdir=outdir,
            orthogroup_id=orthogroup_id,
            extension=out_extension
        ))


def write_orthogroups(orthologs, indir, inextension, outdir, outextension):
    """write_orthogroups: separate the fasta files as in the orthologs file
    """

    species = orthologs.columns.tolist()
    orthogroup_ids = orthologs.index.tolist()

    prepare_files_and_folders(orthogroup_ids, outdir, outextension)

    for orthogroup_id in orthologs.index:
        for single_species in species:
            transcripts = orthologs.loc[orthogroup_id, single_species]
            if pd.isna(transcripts):
                continue
            transcripts = transcripts.replace(",", "")
            os.system(
                "echo {transcripts} | seqtk subseq {infile} - >> {outfile}".format(
                    infile=indir + "/" + single_species + "." + inextension,
                    transcripts=transcripts,
                    outfile=outdir + "/" + orthogroup_id + "." + outextension
                )
            )


if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit("./extract_orthologs.py Orthologs.csv indir inext outdir outext")

    ORTHOLOGS_FN = sys.argv[1]
    INDIR = sys.argv[2]
    INEXT = sys.argv[3]
    OUTDIR = sys.argv[4]
    OUTEXT = sys.argv[5]

    ORTHOLOGS = pd.read_csv(ORTHOLOGS_FN, sep="\t", index_col=0)
    ORTHOLOGS.rename(columns={"Unnamed: 0": "orthogroup_ids"}, inplace=True)
    ORTHOGROUP_IDS = ORTHOLOGS.index.tolist()

    write_orthogroups(ORTHOLOGS, INDIR, INEXT, OUTDIR, OUTEXT)
