def get_models(wildcards):
    return params["selection"]["foreground_branches"][wildcards.group]["models"]

def get_species(wildcards):
    return params["selection"]["foreground_branches"][wildcards.group]["species"]


rule selection_trees_group:
    input:
        tree = TREE + "exabayes/ExaBayes.rooted.nwk",
        msa_folder = HOMOLOGS_REFINE2 + "maxalign",
    output:
        tree_folder = directory(SELECTION + "{group}/trees")
    log: SELECTION + "{group}/trees.log"
    benchmark: SELECTION + "{group}/trees.bmk"
    conda: "selection.yml"
    params:
        species = get_species,
        minimum_foreground_species = params["selection"]["minimum_foreground_species"],
        minimum_background_species = params["selection"]["minimum_background_species"]
    shell:
        """
        python src/homologs/ete3_evol_prepare_folder.py \
            {input.tree} \
            {input.msa_folder} fa \
            {output.tree_folder} nwk \
            {params.species} \
            {params.minimum_foreground_species} \
            {params.minimum_background_species} \
        2> {log} 1>&2
        """

rule selection_trees:
    input:
        expand(
            SELECTION + "trees_{group}",
            group = params["selection"]["foreground_branches"]
        )


rule selection_ete3_group:
    input:
        msa_folder = HOMOLOGS_REFINE2 + "maxalign",
        tree_folder = SELECTION + "{group}/trees"
    output:
        ete3_folder = directory(SELECTION + "{group}/ete3"),
        results = SELECTION + "{group}/ete3.tsv"
    log: SELECTION + "{group}/ete3.log"
    benchmark: SELECTION + "{group}/ete3.bmk"
    conda: "selection.yml"
    threads: MAX_THREADS
    params:
        models = get_models,
        species = get_species
    shell:
        """
        bash src/homologs/ete3_evol_folder.sh \
            {input.tree_folder} \
            {input.msa_folder} \
            {output.ete3_folder} \
            {threads} \
            "{params.models}" \
            {params.species} \
        2> {log} 1>&2

        python src/homologs/parse_ete3_evol_folder.py \
            {output.ete3_folder}/values txt \
        > {output.results} 2>> {log}
        """


rule selection_ete3:
    input:
        expand(
            SELECTION + "{group}/ete3.tsv",
            group = params["selection"]["foreground_branches"]
        )



# rule selection_fasta_group:
#     input:
#         fasta = all.cds
#         tree_folder = SELECTION(groups
#         tsv = SELECTION + "{group}/ete3.tsv"
#     output:
#         fasta = directory(SELECTION + "{group}/fasta")
#     log: SELECTION + "{group}/fasta.log"
#     benchmark: SELECTION + "{group}/fasta.bmk"
#     conda: "selection.yml"
