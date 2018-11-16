CHUNKS_ORTHO = params["orthofinder"]["number_of_chunks"]

rule orthofinder_link:
    """
    Link proteomes into another folder since ORTHOFINDER populates such folder.
    """
    input: FILTERLEN + "{species}.pep"
    output: ORTHOFINDER + "{species}.fasta"
    shell: "ln --symbolic $(readlink --canonicalize {input}) {output}"


rule orthofinder_link_all:
    input:
        fastas = expand(
            ORTHOFINDER + "{species}.fasta",
            species=SPECIES
        )


rule orthofinder_prepare:
    """
    Split fasta files, rename species and sequences and prepare blast databases.
    """
    input:
        fastas = expand(
            ORTHOFINDER + "{species}.fasta",
            species=SPECIES
        )
    output:
        txt = touch(ORTHOFINDER + "prepare.txt"),
        fastas = expand(
            ORTHOFINDER + "Species{species_number}.fa",
            species_number=[x for x in range(0, N_SPECIES)]
        ),
        db = expand(
            ORTHOFINDER + "diamondDBSpecies{database_number}.dmnd",
            database_number=[x for x in range(0, N_SPECIES)]
        )
    params:
        fasta_dir = ORTHOFINDER,
        temp_dir1 = ORTHOFINDER + "Results_*",
        temp_dir2 = ORTHOFINDER + "Results_*/WorkingDirectory/"
    log:
        ORTHOFINDER + "prepare.log"
    benchmark:
        ORTHOFINDER + "prepare.json"
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


rule orthofinder_split:
    """
    Split the headers from transdecoder_longest_orfs into multiple files
    """
    input:
        fai = ancient(ORTHOFINDER + "Species{species_number}.fa.fai")
    output:
        expand(
            ORTHOFINDER + "{{species_number}}/chunks/ids_{chunk_id}.tsv",
            chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)]
        )
    params:
        folder = ORTHOFINDER,
        number_of_chunks = CHUNKS_ORTHO,
        species = "{species_number}"
    log:
        ORTHOFINDER + "{species_number}/split.log"
    benchmark:
        ORTHOFINDER + "{species_number}/split.json"
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



rule orthofinder_blastp:
    """
    Run blastp of each chunk
    """
    input:
        fasta = ORTHOFINDER + "Species{species_number}.fa",
        fai   = ancient(ORTHOFINDER + "Species{species_number}.fa.fai"),
        chunk = ORTHOFINDER + "{species_number}/chunks/ids_{chunk_id}.tsv",
        db    = ORTHOFINDER + "diamondDBSpecies{database_number}.dmnd",
    output:
        tsv = ORTHOFINDER + "{species_number}/{database_number}/blastp_{chunk_id}.tsv"
    log:
        ORTHOFINDER + "{species_number}/{database_number}/blastp_{chunk_id}.log"
    benchmark:
        ORTHOFINDER + "{species_number}/{database_number}/blastp_{chunk_id}.json"
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



rule orthofinder_blastp_merge:
    """
    Merge results from the different blastps
    """
    input:
        expand(
            ORTHOFINDER + "{{species_number}}/{{database_number}}/blastp_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)]
        )
    output:
        tsv = ORTHOFINDER + "Blast{species_number}_{database_number}.txt"
    log:
        ORTHOFINDER + "{species_number}/{database_number}/blastp_merge.log"
    benchmark:
        ORTHOFINDER + "{species_number}/{database_number}/blastp_merge.json"
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
            ORTHOFINDER + "Blast{database_number}_{species_number}.txt",
            species_number = [x for x in range(0,N_SPECIES)],
            database_number = [x for x in range(0,N_SPECIES)]
        )
    output:
        touch(ORTHOFINDER + "groups.ok")
        # cluster_json = ORTHOFINDER + "cluster.json",
        # cluster_log = ORTHOFINDER + "cluster.log",
        # clusters_inflation = ORTHOFINDER + "clusters_OrthoFinder_v2.2.7_I1.5.txt",
        # clusters_pairs = ORTHOFINDER + "clusters_OrthoFinder_v2.2.7_I1.5.txt_id_pairs.txt",
        # graph = ORTHOFINDER + "OrthoFinder_v2.2.7_graph.txt",
        # ortho_csv = ORTHOFINDER + "Orthogroups.csv",
        # gene_count = ORTHOFINDER + "Orthogroups.GeneCount.csv",
        # species_overlaps = ORTHOFINDER + "Orthogroups_SpeciesOverlaps.csv",
        # ortho_txt = ORTHOFINDER + "Orthogroups.txt",
        # unassigned = protected(ORTHOFINDER + "Orthogroups_UnassignedGenes.csv"),
        # sigle_copy = protected(ORTHOFINDER + "SingleCopyOrthogroups.txt"),
        # statistics_overall = protected(ORTHOFINDER + "Statistics_Overall.csv"),
        # statistics_species = protected(ORTHOFINDER + "Statistics_PerSpecies.csv")
    params:
        fasta_dir = ORTHOFINDER,
        inflation = params["orthofinder"]["mcl_inflation"]
    threads:
        8 # There is no reason to go beyond this value
    log:
        ORTHOFINDER + "groups.log"
    benchmark:
        ORTHOFINDER + "groups.json"
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


rule orthofinder_trees:
    input: rules.orthofinder_groups.output
    output:
        touch(ORTHOFINDER + "trees.ok")
    params:
        orthofinder_dir = ORTHOFINDER,
        tree_program = params["orthofinder"]["tree_program"],
        msa_program = params["orthofinder"]["msa_program"]
    conda:
        "orthofinder.yml"
    log:
        ORTHOFINDER + "trees.log"
    benchmark:
        ORTHOFINDER + "trees.bmk"
    threads: 32
    shell:
        """
        orthofinder \
            --from-groups {params.orthofinder_dir} \
            --only-trees \
            --method msa \
            --msa_program {params.msa_program} \
            --tree_program {params.tree_program} \
            --algthreads {threads} \
            --threads {threads} \
        2> {log} 1>&2
        """

rule orthofinder_orthologues:
    input: ORTHOFINDER + "trees.ok"
    output: touch(ORTHOFINDER + "orthologues.ok")
    conda: "orthofinder.yml"
    params:
        orthofinder_dir = ORTHOFINDER,
    log: ORTHOFINDER + "orthologues.log"
    threads: 64
    shell:
        """
        orthofinder \
            --from-trees {params.orthofinder_dir}/Orthologues_*/ \
            --algthreads {threads} \
            --threads {threads} \
        2> {log} 1>&2
        """


rule orthofinder_clean:
    input: ORTHOFINDER + "orthologues.ok"
    output: touch(ORTHOFINDER + "clean.ok")
    params:
        orthofinder_dir = ORTHOFINDER,
        n_species = N_SPECIES - 1
    log: ORTHOFINDER + "clean.log"
    shell:
        """
        pushd {params.orthofinder_dir} 2> {log} 1>&2

        mkdir search/
        for i in {{0..{params.n_species}}}; do
            mv ${{i}} search/
        done

        mv \
            diamondDB* \
            Species*.fa* \
            SequenceIDs.txt \
            SpeciesIDs.txt \
            search/

        mkdir groups/
        mv \
            clusters_OrthoFinder_* \
            OrthoFinder_*_graph.txt \
            Orthogroups* \
            SingleCopyOrthogroups.txt \
            Statistics_Overall.csv \
            Statistics_PerSpecies.csv \
            groups/

        mv \
            Orthologues_*/Alignments \
            Orthologues_*/Gene_Trees \
            Orthologues_*/Sequences \
            .

        mkdir species_tree/
        mv Orthologues_*/SpeciesTree* species_tree

        mv Orthologues_*/New_Analysis_*/* .

        rm -rf Orthologues_*
        """
