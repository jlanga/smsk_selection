rule tree_4ds:
    """
    Get codons that are fourfold-degenerate
    """
    input: HOMOLOGS_REFINE2 + "maxalign"
    output: directory(TREE + "4ds")
    log: TREE + "4ds.log"
    benchmark: TREE + "4ds.bmk"
    conda: "tree.yml"
    shell:
        """
        mkdir --parents {output}

        python2.7 src/homologs/extract_4d_sites_cds.py \
            {input} fa \
            {output} fa \
        2> {log}
        """

rule tree_prepare:
    """
    Strip the transcript id from the sequence names - leave only the species id.

    concatenate_matrices_phyx.py requires the files to be *.aln-cln
    """
    input: TREE + "4ds"
    output: directory(TREE + "prepare")
    log: TREE + "prepare.log"
    benchmark: TREE + "prepare.bmk"
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
    input: TREE + "prepare"
    output:
        fasta = TREE + "supermatrix/supermatrix.fa",
        nex = TREE + "supermatrix/supermatrix.nex",
        model = TREE + "supermatrix/supermatrix.model",
        phy = TREE + "supermatrix/supermatrix.phy",
        stats = TREE + "supermatrix/supermatrix_stats.txt"
    log: TREE + "supermatrix.log"
    benchmark: TREE + "supermatrix.bmk"
    threads: MAX_THREADS
    params:
        min_length = params["tree"]["supermatrix"]["min_length"],
        min_taxa = params["tree"]["supermatrix"]["min_taxa"],
        output_prefix = TREE + "supermatrix/supermatrix",
        stats_tmp = TREE + "supermatrix/supermatrix_taxon_occupancy_stats"
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


rule tree_trim:
    """
    Remove columns that are not very occupied
    """
    input:
        fasta = TREE + "supermatrix/supermatrix.fa"
    output:
        fasta = TREE + "trim/supermatrix.fa",  # For raxml
        phy = TREE + "trim/supermatrix.phy"  # For exabayes
    log: TREE + "trim.log"
    benchmark: TREE + "trim.bmk"
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

        python2.7 src/homologs/fasta_to_phy.py \
        < {output.fasta} \
        > {output.phy} \
        2>> {log}

        rm phyx.logfile
        """


rule tree_modeltestng:
    input: TREE + "trim/supermatrix.fa"
    output: TREE + "modeltestng/modeltestng.out"
    params: TREE + "modeltestng/modeltestng"
    log: TREE + "modeltestng.log"
    benchmark: TREE + "modeltestng.bmk"
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
        fasta = TREE + "trim/supermatrix.fa",
        model = TREE + "modeltestng/modeltestng.out"
    output:
        best_tree = TREE + "raxmlng/supermatrix.raxml.bestTree",
        rooted_tree = TREE + "raxmlng/supermatrix.raxml.rooted.nwk"
    log: TREE + "raxmlng.log"
    benchmark: TREE + "raxmlng.bmk"
    threads: params["tree"]["raxml"]["threads"]
    params:
        bootstraps = params["tree"]["raxml"]["bootstrap_replicates"],
        prefix = TREE + "raxmlng/supermatrix",
        outgroup = params["tree"]["outgroup"]
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

        ./bin/pxrr \
            --outgroups {params.outgroup} \
            --treef {output.best_tree} \
            --outf {output.rooted_tree} \
        2>> {log} 1>&2
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


rule tree_exabayes:
    input:
        config_file = TREE + "exabayes_config.txt",
        phy = TREE + "trim/supermatrix.phy",
        tree = TREE + "raxmlng/supermatrix.raxml.bestTree"
    output:
        sdsf = TREE + "exabayes/ExaBayes.sdsf.txt",
        prostprocparam = TREE + "exabayes/ExaBayes_parameterStatistics.txt",
        consensus = TREE + "exabayes/ExaBayes.consensus.nwk",
        rooted = TREE + "exabayes/ExaBayes.rooted.nwk"
    log: TREE + "exabayes.log"
    benchmark: TREE + "exabayes.bmk"
    conda: "tree.yml"
    threads: params["tree"]["exabayes"]["threads"]
    params:
        out_dir = TREE + "exabayes",
        exabayes = params["tree"]["executables"]["exabayes"],
        sdsf = params["tree"]["executables"]["sdsf"],
        consense = params["tree"]["executables"]["consense"],
        postprocparam = params["tree"]["executables"]["postprocparam"],
        outgroup = params["tree"]["outgroup"]
    shell:
        """
        {params.exabayes} \
            -f {input.phy} \
            -s 12345 \
            -n txt \
            -t {input.tree} \
            -m DNA \
            -c {input.config_file} \
            -w {params.out_dir} \
        2> {log} 1>&2

        topologies=$(find {params.out_dir} -name "ExaBayes_topologies.run*")
        {params.sdsf} -f $topologies > {output.sdsf} 2>&1

        parameters=$(find {params.out_dir} -name "ExaBayes_parameters.run*")
        {params.postprocparam} -n txt -f $parameters 2>>{log} 1>&2
        mv ExaBayes_parameterStatistics.txt {output.prostprocparam}

        {params.consense} -n consensus -f $topologies 2>>{log} 1>&2
        mv ExaBayes_ConsensusExtendedMajorityRuleNewick.consensus {output.consensus}
        rm ExaBayes_ConsensusExtendedMajorityRuleNexus.consensus

        ./bin/pxrr \
            --outgroups {params.outgroup} \
            --treef {output.consensus} \
            --outf {output.rooted} \
        2>> {log} 1>&2
        """

rule tree:
    input:
        rules.tree_raxmlng.output,
        rules.tree_exabayes.output