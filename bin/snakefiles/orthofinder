rule orthofinder_prepare_link_species:
    """
    Link proteomes into another folder since orthofinder populates such folder.
    """
    input:
        pep= tag + "{species}.pep"
    output:
        fasta= temp(
            orthofinder + "{species}.fasta"
        )
    threads:
        1
    log:
        orthofinder + "prepare_link_{species}.log"
    benchmark:
        orthofinder + "prepare_link_{species}.json"
    shell:
        "ln --symbolic $(readlink -f {input.pep}) {output.fasta} 2> {log}"



rule orthofinder_prepare:
    """
    Split fasta files, rename species and sequences and prepare blast databases.
    """
    input:
        fastas = expand(
            orthofinder + "{species}.fasta",
            species = config["species"]
        )
    output:
        txt = touch(
            orthofinder + "prepare.txt"
        ),
        fastas = expand(
            orthofinder + "Species{species_number}.fa",
            species_number = [x for x in range(0,n_species)]
        ),
        db = touch(
            expand(
                orthofinder + "BlastDBSpecies{database_number}",
                database_number = [x for x in range(0,n_species)]
            )
        )
    params:
        fasta_dir = orthofinder
    threads:
        1
    log:
        orthofinder + "prepare.log"
    benchmark:
        orthofinder + "prepare.json"
    shell:
        "docker run "
            "--rm "
            "--volume `pwd`:`pwd` "
            "--workdir `pwd` "
            "--user `id -u $USER`:`id -g $USER` "
            "jlanga/orthofinder-docker "
            "orthofinder.py "
                "--fasta {params.fasta_dir} "
                "--only-prepare "
        "2> {log} 1>&2 ; "
        "mv "
            "{params.fasta_dir}/*/WorkingDirectory/* "
            "{params.fasta_dir}/ "
        "2>> {log} 1>&2"



rule orthofinder_index_fasta_species:
    input:
        fasta = orthofinder + "Species{species_number}.fa"
    output:
        fai = orthofinder + "Species{species_number}.fa.fai"
    log:
        orthofinder + "{species_number}/index_fasta.log"
    shell:
        "samtools faidx {input.fasta} 2> {log}"



rule orthofinder_split_species:
    """
    Split the headers from transdecoder_longest_orfs into multiple files
    """
    input:
        fai = orthofinder + "Species{species_number}.fa.fai"
    output:
        expand(
            orthofinder + "{{species_number}}/chunks/ids_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, config["params"]["orthofinder"]["number_of_chunks"])]
        )
    params:
        number_of_chunks = config["params"]["orthofinder"]["number_of_chunks"],
        species = "{species_number}"
    log:
        orthofinder + "{species_number}/split.log"
    benchmark:
        orthofinder + "{species_number}/split.json"
    shell:
        "split "
            "--number l/{params.number_of_chunks} "
            "--numeric-suffixes "
            "--suffix-length 5 "
            "--additional-suffix .tsv "
            "{input.fai} "
            "{orthofinder}/{params.species}/chunks/ids_ "
        "2> {log}"



rule orthofinder_blastp_species_chunk:
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



rule orthofinder_blastp_merge_species_database_merge:
    """
    Merge results from the different blastps
    """
    input:
        expand(
            orthofinder + "{{species_number}}/{{database_number}}/blastp_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, config["params"]["orthofinder"]["number_of_chunks"])]
        )
    output:
        tsv = orthofinder + "Blast{species_number}_{database_number}.txt"
    log:
        orthofinder + "{species_number}/{database_number}/blastp_merge.log"
    benchmark:
        orthofinder + "{species_number}/{database_number}/blastp_merge.json"
    shell:
        "cat {input} > {output} 2> {log}"



rule orthofinder_groups:
    """
    Join blastp results, normalize bitscores and run mcl.
    """
    input:
        tsv = expand(
                orthofinder + "Blast{database_number}_{species_number}.txt",
                species_number = [x for x in range(0,n_species)],
                database_number = [x for x in range(0,n_species)]
            ),
    output:
        orthologous_groups_csv = protected(
            orthofinder + "Orthogroups.csv"
        ),
        orthologous_groups_txt = protected(
            orthofinder + "Orthogroups.txt"
        )
    params:
        fasta_dir = orthofinder,
        inflation = config["params"]["orthofinder"]["mcl_inflation"]
    threads:
        8 # There is no reason to go beyond this value
    log:
        orthofinder + "cluster.log"
    benchmark:
        orthofinder + "cluster.json"
    shell:
        "docker run "
            "--rm "
            "--volume `pwd`:`pwd` "
            "--workdir `pwd` "
            "--user `id -u $USER`:`id -g $USER` "
            "jlanga/orthofinder-docker "
            "orthofinder.py "
                "--algthreads {threads} "
                "--inflation 1.5 "
                "--blast {params.fasta_dir} "
        "2> {log} 1>&2"
