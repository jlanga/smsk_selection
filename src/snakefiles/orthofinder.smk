rule orthofinder_parse_ids:
    """
    Remove the comments from fasta
    """
    input: CDHIT + "{species}.pep"
    output: ORTHOFINDER + "{species}.fasta"
    shell:
        "cut -f 1 -d \" \" < {input} > {output}"


rule orthofinder_link_all:
    input:
        fastas = expand(
            ORTHOFINDER + "{species}.fasta",
            species=SPECIES
        )


rule orthofinder_groups:
    """
    Make diamond dbs, run diamond, cluster orthogroups
    """
    input:
        fastas = expand(
            ORTHOFINDER + "{species}.fasta",
            species=SPECIES
        )
    output: folder = directory(ORTHOFINDER + "groups/")
    threads: 32
    log: ORTHOFINDER + "groups.log"
    benchmark: ORTHOFINDER + "groups.bmk"
    conda: "orthofinder.yml"
    shell:
        """
        orthofinder \
            --algthreads {threads} \
            --threads {threads} \
            --inflation 1.5 \
            --fasta {ORTHOFINDER} \
            --only-groups \
            --name groups \
        2> {log} 1>&2

        ln --symbolic --relative {OF_GROUPS} {output.folder}
        """


rule orthofinder_trees:
    input: ORTHOFINDER + "groups/"
    output: 
        msa = directory(ORTHOFINDER + "msa"),
        trees = directory(ORTHOFINDER + "gene_trees")
    params:
        trees_results = OF_TREES,
        msa = OF_TREES + "MultipleSequenceAlignments",
        trees = OF_TREES + "Gene_Trees"
    threads: 32
    log: ORTHOFINDER + "trees.log"
    benchmark: ORTHOFINDER + "trees.bmk"
    conda: "orthofinder.yml"
    shell:
        """
        orthofinder \
            --only-trees \
            --algthreads 4 \
            --from-groups {OF_GROUPS} \
            --method msa \
            --threads {threads} \
            --name trees \
        2> {log} 1>&2

        ln --symbolic --relative --force {params.trees} {output.trees}
        ln --symbolic --relative --force {params.msa} {output.msa}
        """


checkpoint orthofinder_orthologues:
    input:
        ORTHOFINDER + "msa",
        ORTHOFINDER + "gene_trees"
    output: directory(ORTHOFINDER + "resolved_gene_trees")
    params: OF_ORTHOLOGUES + "Resolved_Gene_Trees/"
    log: ORTHOFINDER + "orthologues.bmk"
    benchmark: ORTHOFINDER + "orthologues.bmk"
    threads: 32
    conda: "orthofinder.yml"
    shell:
        """
        orthofinder \
            --algthreads 4 \
            --threads {threads} \
            --from-trees {OF_TREES} \
            --name orthologues \
        2> {log} 1>&2

        ln --symbolic --force --relative {params} {output}
        """

rule orthofinder:
    input: ORTHOFINDER + "resolved_gene_trees"