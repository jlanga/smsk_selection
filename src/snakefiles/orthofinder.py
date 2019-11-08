CHUNKS_ORTHO = params["orthofinder"]["number_of_chunks"]

rule orthofinder_search_parse_ids:
    """
    Parse the input transcriptomes and split the fasta header by " "
    """
    input: CDHIT + "{species}.pep"
    output: OF_SEARCH + "{species}.fasta"
    shell:
        """
        cut -f 1 -d \" \" < {input} > {output}
        """


rule orthofinder_search_link_all:
    input:
        fastas = expand(
            OF_SEARCH + "{species}.fasta",
            species=SPECIES
        )


rule orthofinder_search_prepare:
    """
    Split fasta files, rename species and sequences and prepare blast databases.
    """
    input:
        fastas = expand(
            OF_SEARCH + "{species}.fasta",
            species=SPECIES
        )
    output:
        txt = touch(ORTHOFINDER + "prepare.txt"),
        fastas = expand(
            OF_SEARCH + "Species{species_number}.fa",
            species_number=[x for x in range(0, N_SPECIES)]
        ),
        db = expand(
            OF_SEARCH + "diamondDBSpecies{database_number}.dmnd",
            database_number=[x for x in range(0, N_SPECIES)]
        )
    params:
        fasta_dir = OF_SEARCH,
        temp_dir1 = OF_SEARCH + "Results_*",
        temp_dir2 = OF_SEARCH + "Results_*/WorkingDirectory/"
    log:
        OF_SEARCH + "prepare.log"
    benchmark:
        OF_SEARCH + "prepare.json"
    conda:
        "orthofinder.yml"
    shell:
        """
        orthofinder \
            --fasta {params.fasta_dir} \
            --search diamond \
            --only-prepare \
        2> {log} 1>&2

        mv {params.temp_dir2}/* {params.fasta_dir}/
        rm --recursive --force {params.temp_dir1}
        """


rule orthofinder_search_split:
    """
    Split the headers from transdecoder_longest_orfs into multiple files
    """
    input:
        fai = ancient(OF_SEARCH + "Species{species_number}.fa.fai")
    output:
        expand(
            OF_SEARCH + "{{species_number}}/chunks/ids_{chunk_id}.tsv",
            chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)]
        )
    params:
        folder = OF_SEARCH,
        number_of_chunks = CHUNKS_ORTHO,
        species = "{species_number}"
    log:
        OF_SEARCH + "{species_number}/split.log"
    benchmark:
        OF_SEARCH + "{species_number}/split.json"
    conda:
        "orthofinder.yml"
    shell:
        "split "
            "--number l/{params.number_of_chunks} "
            "--numeric-suffixes "
            "--suffix-length 5 "
            "--additional-suffix .tsv "
            "{input.fai} "
            "{params.folder}/{params.species}/chunks/ids_ "
        "2> {log}"



rule orthofinder_search_blastp:
    """
    Run blastp of each chunk
    """
    input:
        fasta = OF_SEARCH + "Species{species_number}.fa",
        fai   = ancient(OF_SEARCH + "Species{species_number}.fa.fai"),
        chunk = OF_SEARCH + "{species_number}/chunks/ids_{chunk_id}.tsv",
        db    = OF_SEARCH + "diamondDBSpecies{database_number}.dmnd",
    output:
        tsv = OF_SEARCH + "{species_number}/{database_number}/blastp_{chunk_id}.tsv"
    log:
        OF_SEARCH + "{species_number}/{database_number}/blastp_{chunk_id}.log"
    benchmark:
        OF_SEARCH + "{species_number}/{database_number}/blastp_{chunk_id}.json"
    conda:
        "orthofinder.yml"
    shell:
        "cut -f 1 {input.chunk} "
        "| xargs samtools faidx {input.fasta}"
        "| diamond blastp "
            "--db {input.db} "
            "--outfmt 6 "
            "--evalue 0.001 "
            "--out {output.tsv} "
            "--threads {threads} "
        "2> {log} 1>&2"



rule orthofinder_search_blastp_merge:
    """
    Merge results from the different blastps
    """
    input:
        expand(
            OF_SEARCH + "{{species_number}}/{{database_number}}/blastp_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)]
        )
    output:
        tsv = OF_SEARCH + "Blast{species_number}_{database_number}.txt"
    log:
        OF_SEARCH + "{species_number}/{database_number}/blastp_merge.log"
    benchmark:
        OF_SEARCH + "{species_number}/{database_number}/blastp_merge.json"
    conda:
        "orthofinder.yml"
    shell:
        "cat {input} > {output} 2> {log}"


rule orthofinder_groups:
    """
    Join blastp results, normalize bitscores and run mcl.
    """
    input:
        tsv = expand(
            OF_SEARCH + "Blast{database_number}_{species_number}.txt",
            species_number = [x for x in range(0,N_SPECIES)],
            database_number = [x for x in range(0,N_SPECIES)]
        )
    output:
        clusters_inflation = OF_GROUPS + "clusters_OrthoFinder_v2.2.7_I1.5.txt",
        clusters_pairs = OF_GROUPS + "clusters_OrthoFinder_v2.2.7_I1.5.txt_id_pairs.txt",
        graph = OF_GROUPS + "OrthoFinder_v2.2.7_graph.txt",
        ortho_csv = OF_GROUPS + "Orthogroups.csv",
        gene_count = OF_GROUPS + "Orthogroups.GeneCount.csv",
        species_overlaps = OF_GROUPS + "Orthogroups_SpeciesOverlaps.csv",
        ortho_txt = OF_GROUPS + "Orthogroups.txt",
        unassigned = protected(OF_GROUPS + "Orthogroups_UnassignedGenes.csv"),
        single_copy = protected(OF_GROUPS + "SingleCopyOrthogroups.txt"),
        statistics_overall = protected(OF_GROUPS + "Statistics_Overall.csv"),
        statistics_species = protected(OF_GROUPS + "Statistics_PerSpecies.csv")
    params:
        fasta_dir = OF_SEARCH,
        inflation = params["orthofinder"]["mcl_inflation"],
        cluster_json = OF_SEARCH + "cluster.json",
        cluster_log = OF_SEARCH + "cluster.log",
        clusters_inflation = OF_SEARCH + "clusters_OrthoFinder_v2.2.7_I1.5.txt",
        clusters_pairs = OF_SEARCH + "clusters_OrthoFinder_v2.2.7_I1.5.txt_id_pairs.txt",
        graph = OF_SEARCH + "OrthoFinder_v2.2.7_graph.txt",
        ortho_csv = OF_SEARCH + "Orthogroups.csv",
        gene_count = OF_SEARCH + "Orthogroups.GeneCount.csv",
        species_overlaps = OF_SEARCH + "Orthogroups_SpeciesOverlaps.csv",
        ortho_txt = OF_SEARCH + "Orthogroups.txt",
        unassigned = OF_SEARCH + "Orthogroups_UnassignedGenes.csv",
        single_copy = OF_SEARCH + "SingleCopyOrthogroups.txt",
        statistics_overall = OF_SEARCH + "Statistics_Overall.csv",
        statistics_species = OF_SEARCH + "Statistics_PerSpecies.csv"
    threads:
        8 # There is no reason to go beyond this value
    log:
        OF_GROUPS + "groups.log"
    benchmark:
        OF_GROUPS + "groups.json"
    conda:
        "orthofinder.yml"
    shell:
        """
        orthofinder \
            --algthreads {threads} \
            --inflation 1.5 \
            --blast {params.fasta_dir} \
            --only-groups \
        2> {log} 1>&2


        mv {params.clusters_inflation} {output.clusters_inflation}
        mv {params.clusters_pairs} {output.clusters_pairs}
        mv {params.graph} {output.graph}
        mv {params.ortho_csv} {output.ortho_csv}
        mv {params.gene_count} {output.gene_count}
        mv {params.species_overlaps} {output.species_overlaps}
        mv {params.ortho_txt} {output.ortho_txt}
        mv {params.unassigned} {output.unassigned}
        mv {params.single_copy} {output.single_copy}
        mv {params.statistics_overall} {output.statistics_overall}
        mv {params.statistics_species} {output.statistics_species}
        """
