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
    output: touch(ORTHOFINDER + "groups.ok")
    threads: MAX_THREADS
    log: ORTHOFINDER + "groups.log"
    benchmark: ORTHOFINDER + "groups.bmk"
    conda: "orthofinder.yml"
    shell:
        """
        # Sometimes complains because it exists
        rm -rf {ORTHOFINDER}/OrthoFinder/Results_groups
        
        orthofinder \
            --algthreads {threads} \
            --threads {threads} \
            --inflation 1.5 \
            --fasta {ORTHOFINDER} \
            --only-groups \
            --name groups \
        2> {log} 1>&2
      """


rule orthofinder_trees:
    input: ORTHOFINDER + "groups.ok"
    output: touch(ORTHOFINDER + "trees.ok")
    params: OF_GROUPS
    threads: MAX_THREADS
    log: ORTHOFINDER + "trees.log"
    benchmark: ORTHOFINDER + "trees.bmk"
    conda: "orthofinder.yml"
    shell:
        """
        # Sometimes complains because it exists
        rm -rf {ORTHOFINDER}/OrthoFinder/Results_trees

        orthofinder \
            --only-trees \
            --algthreads 4 \
            --from-groups {params} \
            --method msa \
            --threads {threads} \
            --name trees \
        2> {log} 1>&2
        """


rule orthofinder_orthologues:
    input: ORTHOFINDER + "trees.ok"
    output: touch(ORTHOFINDER + "orthologues.ok")
    params: OF_TREES
    log: ORTHOFINDER + "orthologues.log"
    benchmark: ORTHOFINDER + "orthologues.bmk"
    threads: MAX_THREADS
    conda: "orthofinder.yml"
    shell:
        """
        # Sometimes complains because it exists
        rm -rf {ORTHOFINDER}/OrthoFinder/Results_orthologues

        orthofinder \
            --algthreads {threads} \
            --threads {threads} \
            --from-trees {params} \
            --name orthologues \
        2> {log} 1>&2
        """

rule orthofinder:
    input: ORTHOFINDER + "orthologues.ok"
