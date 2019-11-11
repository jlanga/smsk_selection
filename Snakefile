import pandas as pd
import yaml

from snakemake.utils import min_version
min_version("5.7.4")

shell.prefix("set -euo pipefail;")

params = yaml.safe_load(open("params.yml", "r"))
features = yaml.safe_load(open("features.yml", "r"))
samples = pd.read_csv("samples.tsv", sep="\t").set_index("species")

singularity: "docker://continuumio/miniconda3:4.4.10"

SPECIES = samples.index.tolist()
N_SPECIES = len(SPECIES)

snakefiles = "src/snakefiles/"

include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "generic.py"
include: snakefiles + "raw.py"
include: snakefiles + "download.py"
include: snakefiles + "db.py"
include: snakefiles + "busco.py"
include: snakefiles + "transdecoder.py"
include: snakefiles + "cdhit.py"
include: snakefiles + "tidy.py"
include: snakefiles + "orthofinder.py"
#include: snakefiles + "homologs.py"

rule all:
    input:
        # DOWNLOAD + "uniref90.fa.gz",
        # DOWNLOAD + "swissprot.fa.gz",
        # DOWNLOAD + "pfama.hmm.gz",
        # expand(
        #     TRANSDECODER + "{species}.pep",
        #     species=species
        # ),
        # expand(
        #     TAG + "{species}.{extension}",
        #     species=species,
        #     extension="pep cds".split()
        # ),
        #OF_GROUPS + "Orthogroups.csv",
        rules.orthofinder.input,
        OF_SEQUENCES_PEP,
        OF_SEQUENCES_CDS
