CHUNKS_ORTHO = params["orthofinder"]["number_of_chunks"]

rule orthofinder_link:
    """
    Link proteomes into another folder since orthofinder populates such folder.
    """
    input: filterlen + "{species}.pep"
    output: temp(orthofinder + "{species}.fasta")
    shell: "ln --symbolic $(readlink --canonicalize {input}) {output}"


rule orthofinder_prepare:
    """
    Split fasta files, rename species and sequences and prepare blast databases.
    """
    input:
        fastas = expand(
            orthofinder + "{species}.fasta",
            species=SPECIES
        )
    output:
        txt = touch(orthofinder + "prepare.txt"),
        fastas = expand(
            orthofinder + "Species{species_number}.fa",
            species_number = [x for x in range(0, N_SPECIES)]
        ),
        db = touch(
            expand(
                orthofinder + "BlastDBSpecies{database_number}",
                database_number = [x for x in range(0, N_SPECIES)]
            )
        )
    params:
        fasta_dir = orthofinder,
        temp_dir1 = orthofinder + "Results_*",
        temp_dir2 = orthofinder + "Results_*/WorkingDirectory/"
    log:
        orthofinder + "prepare.log"
    benchmark:
        orthofinder + "prepare.json"
    conda:
        "orthofinder.yml"
    shell:
        """
        orthofinder --fasta {params.fasta_dir} --only-prepare 2> {log} 1>&2

        mv {params.temp_dir2}/* {params.fasta_dir}/
        rm --recursive --force {params.temp_dir1}
        """


rule orthofinder_split:
    """
    Split the headers from transdecoder_longest_orfs into multiple files
    """
    input:
        fai = orthofinder + "Species{species_number}.fa.fai"
    output:
        expand(
            orthofinder + "{{species_number}}/chunks/ids_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)]
        )
    params:
        number_of_chunks = CHUNKS_ORTHO,
        species = "{species_number}"
    log:
        orthofinder + "{species_number}/split.log"
    benchmark:
        orthofinder + "{species_number}/split.json"
    conda:
        "orthofinder.yml"
    shell:
        "split "
            "--number l/{params.number_of_chunks} "
            "--numeric-suffixes "
            "--suffix-length 5 "
            "--additional-suffix .tsv "
            "{input.fai} "
            "{orthofinder}/{params.species}/chunks/ids_ "
        "2> {log}"



rule orthofinder_blastp:
    """
    Run blastp of each chunk
    """
    input:
        fasta = orthofinder + "Species{species_number}.fa",
        fai   = orthofinder + "Species{species_number}.fa.fai",
        chunk = orthofinder + "{species_number}/chunks/ids_{chunk_id}.tsv",
        db    = orthofinder + "BlastDBSpecies{database_number}",
    output:
        tsv = orthofinder + "{species_number}/{database_number}/blastp_{chunk_id}.tsv"
    log:
        orthofinder + "{species_number}/{database_number}/blastp_{chunk_id}.log"
    benchmark:
        orthofinder + "{species_number}/{database_number}/blastp_{chunk_id}.json"
    conda:
        "orthofinder.yml"
    shell:
        "cut -f 1 {input.chunk} "
        "| xargs samtools faidx {input.fasta} "
        "| blastp "
            "-db {input.db} "
            "-max_target_seqs 1 "
            "-outfmt 6 "
            "-evalue 1e-5 "
            "-out {output.tsv} "
        "2> {log} 1>&2"



rule orthofinder_blastp_merge:
    """
    Merge results from the different blastps
    """
    input:
        expand(
            orthofinder + "{{species_number}}/{{database_number}}/blastp_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)]
        )
    output:
        tsv = orthofinder + "Blast{species_number}_{database_number}.txt"
    log:
        orthofinder + "{species_number}/{database_number}/blastp_merge.log"
    benchmark:
        orthofinder + "{species_number}/{database_number}/blastp_merge.json"
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
                orthofinder + "Blast{database_number}_{species_number}.txt",
                species_number = [x for x in range(0,N_SPECIES)],
                database_number = [x for x in range(0,N_SPECIES)]
            ),
    output:
        cluster_json = protected(orthofinder + "cluster.json"),
        cluster_log = protected(orthofinder + "cluster.log"),
        clusters_inflation = protected(orthofinder + "clusters_OrthoFinder_v2.2.7_I1.5.txt"),
        clusters_pairs = protected(
            orthofinder + "clusters_OrthoFinder_v2.2.7_I1.5.txt_id_pairs.txt"
        ),
        graph = protected(orthofinder + "OrthoFinder_v2.2.7_graph.txt"),
        ortho_csv = protected(orthofinder + "Orthogroups.csv"),
        gene_count = protected(orthofinder + "Orthogroups.GeneCount.csv"),
        species_overlaps = protected(orthofinder + "Orthogroups_SpeciesOverlaps.csv"),
        ortho_txt = protected(orthofinder + "Orthogroups.txt"),
        unassigned = protected(orthofinder + "Orthogroups_UnassignedGenes.csv"),
        sigle_copy = protected(orthofinder + "SingleCopyOrthogroups.txt"),
        statistics_overall = protected(orthofinder + "Statistics_Overall.csv"),
        statistics_species = protected(orthofinder + "Statistics_PerSpecies.csv")
    params:
        fasta_dir = orthofinder,
        inflation = params["orthofinder"]["mcl_inflation"]
    threads:
        8 # There is no reason to go beyond this value
    log:
        orthofinder + "cluster.log"
    benchmark:
        orthofinder + "cluster.json"
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
        """


rule orthofhinder_msa:
    input:
        orthogroups =
    conda: "orthofinder.yml"
    shell:


rule orthofinder:
    input: orthofinder + "Orthogroups.csv"
