import pandas as pd
import yaml

from snakemake.utils import min_version
min_version("5.0")

shell.prefix("set -euo pipefail;")

params = yaml.load(open("params.yml", "r"))
features = yaml.load(open("features.yml", "r"))
samples = pd.read_table("samples.tsv").set_index("species")

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
include: snakefiles + "tag.py"
include: snakefiles + "filterlen.py"
include: snakefiles + "orthofinder.py"
# include: snakefiles + "orthogroups"

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
        ORTHOFINDER + "clean.ok",
        rules.busco.input
