CHUNKS_ORTHO = params["orthofinder"]["number_of_chunks"]

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


# rule orthofinder_prepare:
#     """
#     Split fasta files, rename species and sequences and prepare blast databases.
#     """
#     input:
#         fastas = expand(
#             ORTHOFINDER + "{species}.fasta",
#             species=SPECIES
#         )
#     output:
#         txt = touch(ORTHOFINDER + "prepare.ok"),
#     log:
#         ORTHOFINDER + "prepare.log"
#     benchmark:
#         ORTHOFINDER + "prepare.json"
#     conda:
#         "orthofinder.yml"
#     shell:
#         """
#         orthofinder \
#             --fasta {ORTHOFINDER} \
#             --search diamond \
#             --only-prepare \
#             --name prepare \
#         2> {log} 1>&2
#         """
        # mv {params.temp_dir2}/* {params.fasta_dir}/
        # rm --recursive --force {params.temp_dir1}
        # """


# rule orthofinder_blastp_split:
#     """
#     Split the headers from transdecoder_longest_orfs into multiple files
#     """
#     input:
#         fai = ancient(ORTHOFINDER + "Species{species_number}.fa.fai")
#     output:
#         expand(
#             ORTHOFINDER + "{{species_number}}/chunks/ids_{chunk_id}.tsv",
#             chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)]
#         )
#     params:
#         folder = ORTHOFINDER,
#         number_of_chunks = CHUNKS_ORTHO,
#         species = "{species_number}"
#     log:
#         ORTHOFINDER + "{species_number}/split.log"
#     benchmark:
#         ORTHOFINDER + "{species_number}/split.json"
#     conda:
#         "orthofinder.yml"
#     shell:
#         "split "
#             "--number l/{params.number_of_chunks} "
#             "--numeric-suffixes "
#             "--suffix-length 5 "
#             "--additional-suffix .tsv "
#             "{input.fai} "
#             "{params.folder}/{params.species}/chunks/ids_ "
#         "2> {log}"



# rule orthofinder_blastp:
#     """
#     Run blastp of each chunk
#     """
#     input:
#         fasta = ORTHOFINDER + "Species{species_number}.fa",
#         fai   = ancient(ORTHOFINDER + "Species{species_number}.fa.fai"),
#         chunk = ORTHOFINDER + "{species_number}/chunks/ids_{chunk_id}.tsv",
#         db    = ORTHOFINDER + "diamondDBSpecies{database_number}.dmnd",
#     output:
#         tsv = ORTHOFINDER + "{species_number}/{database_number}/blastp_{chunk_id}.tsv"
#     log:
#         ORTHOFINDER + "{species_number}/{database_number}/blastp_{chunk_id}.log"
#     benchmark:
#         ORTHOFINDER + "{species_number}/{database_number}/blastp_{chunk_id}.json"
#     conda:
#         "orthofinder.yml"
#     shell:
#         "cut -f 1 {input.chunk} "
#         "| seqtk subseq {input.fasta} - "
#         "| diamond blastp "
#             "--db {input.db} "
#             "--outfmt 6 "
#             "--evalue 0.001 "
#             "--out {output.tsv} "
#             "--threads {threads} "
#         "2> {log} 1>&2"



# rule orthofinder_blastp_merge:
#     """
#     Merge results from the different blastps
#     """
#     input:
#         expand(
#             ORTHOFINDER + "{{species_number}}/{{database_number}}/blastp_{chunk_id}.tsv",
#             chunk_id = ['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)]
#         )
#     output:
#         tsv = ORTHOFINDER + "Blast{species_number}_{database_number}.txt"
#     log:
#         ORTHOFINDER + "{species_number}/{database_number}/blastp_merge.log"
#     benchmark:
#         ORTHOFINDER + "{species_number}/{database_number}/blastp_merge.json"
#     conda:
#         "orthofinder.yml"
#     shell:
#         "cat {input} > {output} 2> {log}"

# rule orthofinder_search:
#     input: ORTHOFINDER + "prepare.ok"
#     output: touch(ORTHOFINDER + "search.ok")
#     shell:
#         """
#         orthofinder \
#             --algthreads {threads} \
#             --threads {threads} \
#             --blast 
#             --fasta 
#             --only
#         """


rule orthofinder_groups:
    """
    Join blastp results, normalize bitscores and run mcl.
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