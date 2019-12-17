#!/usr/bin/env python2
"""
extract_orthologs.py orthologs.tsv in_dir outdir in_extension out_extension
Process all fasta files in in_dir that have extension .in_extension and group
them clustered as in orthologs.tsv, writting them in outdir with extension
out_extension
input fasta files have to be named species.in_extension
"""
import sys

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

def get_species_sequences(orthogroups, in_dir, in_ext):
    """Build {species: {identifier: sequence}} dict
    """
    species_list = orthogroups.columns
    dict_of_sequences = {}
    for species in species_list:
        filename = in_dir + "/" + species + "." + in_ext
        with open(filename, "r") as handle:
            dict_of_sequences[species] = {
                seq_id.split()[0]: sequence
                for seq_id, sequence in SimpleFastaParser(handle)
            }
    return dict_of_sequences


def assign_sequences_to_orthogroups(dict_of_sequences, orthogroups):
    """Reformat the dict_of_sequences as {orthogroup_id : {seq_id, seq}}
    """
    dict_of_orthogroups = {orthogroup_id: {} for orthogroup_id in orthogroups.index}

    for orthogroup_id in orthogroups.index:
        for species_id in orthogroups.columns:
            elements = orthogroups\
                .loc[orthogroup_id, species_id]
            if pd.isna(elements):
                continue
            elements = elements\
                .replace(",", "")\
                .split(" ")
            for element in elements:
                dict_of_orthogroups[orthogroup_id][element] = \
                    dict_of_sequences[species_id][element]

    return dict_of_orthogroups


def write_orthogroups(dict_of_orthogroups, out_dir, out_ext):
    """Write dict_of_orthogroups to files"""
    for orthogroup_id, dict_of_sequences in dict_of_orthogroups.items():
        filename = "%s/%s.%s" % (out_dir, orthogroup_id, out_ext)
        with open(filename, "w") as handle:
            for identifier, sequence in dict_of_sequences.items():
                handle.write(">%s\n%s\n" % (identifier, sequence))


if __name__ == '__main__':

    if len(sys.argv) != 6:
        sys.exit(
            "python2 extract_orthologs.py Orthogroups.csv in_dir in_ext out_dir "
            "out_ext"
        )

    ORTHOGROUPS_FN = sys.argv[1]
    IN_DIR = sys.argv[2]
    IN_EXT = sys.argv[3]
    OUT_DIR = sys.argv[4]
    OUT_EXT = sys.argv[5]

    ORTHOGROUPS = pd\
        .read_csv(ORTHOGROUPS_FN, sep="\t", index_col=0)

    SPECIES_TO_SEQUENCES = get_species_sequences(
        orthogroups=ORTHOGROUPS,
        in_dir=IN_DIR,
        in_ext=IN_EXT
    )

    ORTHOGROUPS_TO_SEQUENCES = assign_sequences_to_orthogroups(
        dict_of_sequences=SPECIES_TO_SEQUENCES,
        orthogroups=ORTHOGROUPS
    )

    write_orthogroups(
        ORTHOGROUPS_TO_SEQUENCES,
        out_dir=OUT_DIR,
        out_ext=OUT_EXT
    )
