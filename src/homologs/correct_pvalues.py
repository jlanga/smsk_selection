#!/usr/bin/env python

"""correct_pvalues.py: get the results from fastcodeml and apply BH correction 
gene by gene for the results from ete3/codeml and fastcodeml
"""

import argparse

import pandas as pd
from statsmodels.stats.multitest import multipletests


def process_pvalues(tsv_ete3, tsv_fastcodeml, tsv_output, threshold_pvalue=0.05):
    """Compute BH correction for p-values"""

    # Read tables
    fastcodeml = pd.read_csv(
            filepath_or_buffer=tsv_fastcodeml, sep="\t"
        )\
        [["orthogroup", "omega_zero", "pvalue"]]\
        .assign(program="fastcodeml")

    ete3 = pd.read_csv(
            filepath_or_buffer=tsv_ete3,
            sep="\t"
        )\
        [["orthogroup", "omega_zero", "pvalue"]]

    # Filter ete3 results that already were negative
    fastcodeml_orthogroups = fastcodeml.orthogroup.values

    ete3 = ete3\
        .loc[ete3.orthogroup.isin(fastcodeml_orthogroups)]\
        .assign(program="ete3")\
        .reset_index(drop=True)

    pvalues = pd\
        .concat([ete3, fastcodeml])\
        .sort_values(by=["orthogroup", "program", "omega_zero"])\
        .reset_index(drop=True)

    del fastcodeml_orthogroups, ete3, fastcodeml

    # Compute Benjamini-Hochberch
    bh = pvalues\
        .groupby("orthogroup")\
        .apply(
            lambda x: multipletests(
                x.pvalue,
                alpha=threshold_pvalue, method="fdr_bh"
            )[1]  # Just take the p-values
        )\
        .explode()\
        .rename("bh")

    pvalues["bh"] = bh.values
    del bh

    selected = pvalues\
        .groupby("orthogroup")\
        .apply(
            lambda x: all(x.bh < threshold_pvalue)
        )\
        .rename("selected")\
        .to_frame()\
        .reset_index()

    pvalues = pd.merge(pvalues, selected)
    del selected

    pvalues.to_csv(tsv_output, "\t", index=False)



def parse_arguments():
    """parser for correct_pvalues.py"""
    parser = argparse.ArgumentParser(
        description="correct_pvalues.py: join the results from ete3 and "
        "fastcodeml, perform Benjamini-Hochberch correction per orthogroup, and"
        "mark which p-values are under the target threshold"
    )
    parser.add_argument(
        '-e', '--ete3',
        help='Path to ete3 p-values (TSV)',
        required=True
    )
    parser.add_argument(
        '-f', '--fastcodeml',
        help='Path to fastcodeml p-values (TSV)',
        required=True
    )
    parser.add_argument(
        '-o', '--output',
        help='Path to output (TSV)',
        required=True
    )
    parser.add_argument(
        '-p', '--pvalue',
        help='Threshold p-value for BH-correction and selection',
        default=0.05,
        type=float
    )
    return vars(parser.parse_args())


if __name__ == '__main__':

    ARGS = parse_arguments()
    process_pvalues(
        tsv_ete3=ARGS["ete3"],
        tsv_fastcodeml=ARGS["fastcodeml"],
        tsv_output=ARGS["output"],
        threshold_pvalue=ARGS["pvalue"]
    )
