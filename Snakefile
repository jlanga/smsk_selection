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

include: snakefiles + "folders.smk"
include: snakefiles + "clean.smk"
include: snakefiles + "generic.smk"
include: snakefiles + "raw.smk"
include: snakefiles + "download.smk"
include: snakefiles + "db.smk"
include: snakefiles + "busco.smk"
include: snakefiles + "transdecoder.smk"
include: snakefiles + "cdhit.smk"
include: snakefiles + "tidy.smk"
include: snakefiles + "orthofinder.smk"
include: snakefiles + "homologs.smk"

rule all:
    input:
        rules.orthofinder.input,
        rules.homologs_round1_treeshrink.output
