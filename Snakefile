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
include: snakefiles + "transdecoder.py"
include: snakefiles + "tag.py"
include: snakefiles + "filterlen.py"
include: snakefiles + "orthofinder.py"
# include: snakefiles + "orthogroups"

rule all:
    input:
        # download + "uniref90.fa.gz",
        # download + "swissprot.fa.gz",
        # download + "pfama.hmm.gz",
        # expand(
        #     transdecoder + "{species}.pep",
        #     species = species
        # ),
        # expand(
        #     tag + "{species}.{extension}",
        #     species = species,
        #     extension = "pep cds".split()
        # ),
        orthofinder + "prepare.txt"
        # orthofinder + "Orthogroups.csv",
        #expand(
        #    orthofinder + "{species}.pep.fai",
        #    species = species
        #),
        #expand(
        #    orthofinder + "{species}.cds.fai",
        #    species = species
        #),
        # orthogroups + "extract_orthogroups_pep.txt",
        # orthogroups + "extract_orthogroups_cds.txt",
        # orthogroups + "remove_stops_pep.txt",
        # orthogroups + "trim_pep.txt"
        # dynamic(
        #    orthogroups + "fastcodeml/{orthogroup_id}.fastcodeml"
        # )
