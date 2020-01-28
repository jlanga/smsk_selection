rule tree_fourfold_degenerate_sites:
    """
    Get codons that are fourfold-degenerate
    """
    input: HOMOLOGS_REFINE2 + "maxalign_cds"
    output: directory(TREE + "fourfold_degenerate_sites")
    log: TREE + "fourfold_degenerate_sites.log"
    benchmark: TREE + "fourfold_degenerate_sites.bmk"
    conda: "tree.yml"
    threads: MAX_THREADS
    shell:
        """
        mkdir --parents {output}

        (find {input} -name "*.fa" -type f \
        | sort -V \
        | parallel --keep-order \
            --jobs {threads} \
            python2.7 src/homologs/extract_4d_sites_cds.py \
                {input}/{{/.}}.fa \
                {output}/{{/.}}.fa \
        ) 2> {log}
        """

rule tree_supermatrix_prepare:
    """
    Strip the transcript id from the sequence names - leave only the species id.

    concatenate_matrices_phyx.py requires the files to be *.aln-cln
    """
    input: TREE + "fourfold_degenerate_sites"
    output: directory(TREE + "supermatrix_prepare")
    log: TREE + "supermatrix_prepare.log"
    benchmark: TREE + "supermatrix_prepare.bmk"
    threads: MAX_THREADS
    conda: "tree.yml"
    shell:
        """
        mkdir --parent {output}

        (find {input} -name "*.fa" \
        | sort -V \
        | parallel \
            --jobs {threads} \
            cut -f 1 -d @ \
                "<" {input}/{{/.}}.fa \
                ">" {output}/{{/.}}.aln-cln \
        ) 2> {log} 1>&2
        """



rule tree_supermatrix:
    """
    Concatenate the alignments in a fasta file and create the partition file for
    raxml.
    """
    input: TREE + "supermatrix_prepare"
    output:
        fasta = TREE + "supermatrix.fa",
        nex = TREE + "supermatrix.nex",
        model = TREE + "supermatrix.model",
        phy = TREE + "supermatrix.phy",
        stats = TREE + "supermatrix_stats.txt"
    log: TREE + "supermatrix.log"
    benchmark: TREE + "supermatrix.bmk"
    threads: MAX_THREADS
    params:
        min_length = params["tree"]["supermatrix"]["min_length"],
        min_taxa = params["tree"]["supermatrix"]["min_taxa"],
        output_prefix = TREE + "supermatrix",
        stats_tmp = TREE + "supermatrix_taxon_occupancy_stats"
    conda: "tree.yml"
    shell:
        """
        PATH="bin:$PATH"

        python2.7 src/pdc2/scripts/concatenate_matrices_phyx.py \
            {input} \
            {params.min_length} \
            {params.min_taxa} \
            {params.output_prefix} \
        2> {log} 1>&2

        mv {params.stats_tmp} {output.stats}
        """


rule tree_phyx_trim_cols:
    """
    Remove columns that are not
    """
    input: TREE + "supermatrix.fa"
    output: TREE + "supermatrix_hq.fa"
    log: TREE + "supermatrix_hq.log"
    benchmark: TREE + "supermatrix_hq.bmk"
    params:
        params["tree"]["min_occupation"]  # proportion of data that is required to be present
    conda: "tree.yml"
    shell:
        """
        ./bin/pxclsq \
            --seqf {input} \
            --outf {output} \
            --prop {params} \
        2> {log} 1>&2

        rm phyx.logfile
        """


rule tree_modeltest:
    input: TREE + "supermatrix_hq.fa"
    output: TREE + "supermatrix_hq.modeltest-ng.out"
    params: TREE + "supermatrix_hq.modeltest-ng"
    log: TREE + "modeltestng.log"
    benchmark: TREE + "modeltest.bmk"
    threads: MAX_THREADS
    conda: "tree.yml"
    shell:
        """
        modeltest-ng \
            --datatype nt \
            --input {input} \
            --output {params} \
            --processes {threads} \
            --rngseed 1 \
        2> /dev/null 1>&2
        """


rule tree_raxmlng:
    input: 
        fasta = TREE + "supermatrix_hq.fa",
        model = TREE + "supermatrix_hq.modeltest-ng.out"
    output: TREE + "supermatrix_hq.raxml.bestTree"
    log: TREE + "raxmlmg.log"
    benchmark: TREE + "raxmlmg.bmk"
    threads: params["tree"]["raxml"]["threads"]
    params:
        bootstraps = params["tree"]["raxml"]["bootstrap_replicates"],
        prefix = TREE + "supermatrix_hq"
    conda: "tree.yml"
    shell:
        """
        model=$(grep raxml-ng {input.model} \
            | tail -1 \
            | grep -Eo -- '--model [A-Za-z0-9]+' \
            | cut -f 2 -d " " \
        )

        raxml-ng \
            --all \
            --msa {input.fasta} \
            --model $model \
            --bs-trees {params.bootstraps} \
            --prefix {params.prefix} \
            --threads {threads} \
        2> {log} 1>&2
        """

rule tree_exabayes_config:
    output: TREE + "exabayes_config.txt"
    log: TREE + "exabayes_config.log"
    benchmark: TREE + "exabayes_config.bmk"
    conda: "tree.yml"
    params:
        mcmc_runs = params["tree"]["exabayes"]["mcmc_runs"],
        coupled_chains = params["tree"]["exabayes"]["coupled_chains"],
        generations = params["tree"]["exabayes"]["generations"],
        sampling_frequency = params["tree"]["exabayes"]["sampling_frequency"]
    shell:
        """
        cat > {output} <<ENDOFTEXT
begin RUN;
    numRuns {params.mcmc_runs}
    numGen {params.generations}
    samplingFreq {params.sampling_frequency}
    numCoupledChains {params.coupled_chains}
end;
ENDOFTEXT
        """


    # """
    # exabayes \
    #     -f {input.phy} \
    #     -s 12345 \
    #     -n exabayes \
    #     .m DNA \
    #     -z \
    #     -w {TREE} \
    #     -C 
    # """

rule tree:
    input: rules.tree_raxmlng.output