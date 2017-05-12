shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

species = [x for x in config["species"]]
n_species = len(species)

snakefiles = "bin/snakefiles/"

include: snakefiles + "folders"
include: snakefiles + "clean"
include: snakefiles + "raw"
include: snakefiles + "download"
include: snakefiles + "db"
include: snakefiles + "transdecoder"
include: snakefiles + "tag"
include: snakefiles + "orthofinder"
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
        #orthofinder + "prepare.txt"
        orthofinder + "Orthogroups.csv",
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
