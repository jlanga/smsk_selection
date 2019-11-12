rule homologs_round1_mafft:
    input:
        OF_SEQUENCES + "{orthogroup_id}.pep"
    output:
        HOMOLOGS_R1 + "{orthogroup_id}.mafft.pep"
    log:
        HOMOLOGS_R1 + "{orthogroup_id}.mafft.pep.log"
    benchmark:
        HOMOLOGS_R1 + "{orthogroup_id}.mafft.pep.bmk"
    conda:
        "homologs.yml"
    shell:
        """
        if [ -s {input} ]; then
            mafft --genafpair --maxiterate 1000 {input} > {output} 2> {log}
        else
            touch {output} 2> {log}
        fi
        """




rule homologs_round1_trimal_pep:
    input: HOMOLOGS_R1 + "{orthogroup_id}.mafft.pep"
    output: HOMOLOGS_R1 + "{orthogroup_id}.trimal.pep"
    log: HOMOLOGS_R1 + "{orthogroup_id}.trimal.log"
    benchmark: HOMOLOGS_R1 + "{orthogroup_id}.trimal.bmk"
    conda: "homologs.yml"
    shell:
        """
        if [ -s {input} ]; then
            trimal -in {input} -out {output} -automated1 2> {log} 1>&2
        else
            touch {output} 2> {log} 1>&2
        fi
        """


rule homologs_round1_trimal_cds:
    input:
        cds_raw = OF_SEQUENCES + "{orthogroup_id}.cds",
        pep_aligned = HOMOLOGS_R1 + "{orthogroup_id}.mafft.pep"
    output: HOMOLOGS_R1 + "{orthogroup_id}.trimal.cds"
    threads: 1
    log: HOMOLOGS_R1 + "{orthogroup_id}.trimal.cds.log"
    benchmark: HOMOLOGS_R1 + "{orthogroup_id}.trimal.cds.bmk"
    conda: "homologs.yml"
    shell:
        """
        if [ -s {input.cds_raw} ] && [ -s {input.pep_aligned} ]; then
            trimal \
                -in {input.pep_aligned} \
                -backtrans {input.cds_raw} \
                -out {output} \
            2> {log} 1>&2
        else
            touch {output} 2> {log} 1>&2
        fi
        """


rule homologs_round1_modeltestng:
    input: HOMOLOGS_R1 + "{orthogroup_id}.trimal.pep"
    output:
        ckp = HOMOLOGS_R1 + "{orthogroup_id}.modeltest-ng.ckp",
        out = HOMOLOGS_R1 + "{orthogroup_id}.modeltest-ng.out",
        tree = HOMOLOGS_R1 + "{orthogroup_id}.modeltest-ng.tree"
    params:
        prefix = HOMOLOGS_R1 + "{orthogroup_id}.modeltest-ng"
    threads: 1
    benchmark: HOMOLOGS_R1 + "{orthogroup_id}.modeltest-ng.bmk"
    conda: "homologs.yml"
    shell:
        """
        if [ -s {input} ]; then
            modeltest-ng \
                --datatype aa \
                --input {input} \
                --template raxml \
                --out {params.prefix} \
            2> /dev/null 1>&2
        else
            touch {output.ckp} {output.out} {output.tree}
        fi
        """


rule homologs_round1_raxmlng:
    input:
        alignment = HOMOLOGS_R1 + "{orthogroup_id}.trimal.pep",
        model_file = HOMOLOGS_R1 + "{orthogroup_id}.modeltest-ng.out"
    output:
        best_model = touch(HOMOLOGS_R1 + "{orthogroup_id}.raxml.bestModel"),
        best_tree = touch(HOMOLOGS_R1 + "{orthogroup_id}.raxml.bestTree"),
        ml_trees = touch(HOMOLOGS_R1 + "{orthogroup_id}.raxml.mlTrees"),
        rba = touch(HOMOLOGS_R1 + "{orthogroup_id}.raxml.rba"),
        start_tree = touch(HOMOLOGS_R1 + "{orthogroup_id}.raxml.startTree")
    params:
        prefix = HOMOLOGS_R1 + "{orthogroup_id}"
    threads: 1
    log: HOMOLOGS_R1 + "{orthogroup_id}.raxml.log"
    benchmark: HOMOLOGS_R1 + "{orthogroup_id}.raxml.bmk"
    conda: "homologs.yml"
    shell:
        """
        if [ -s {input.alignment} ] && \
            [ -s {input.model_file} ] && \
            [ $(grep -c ^">" {input.alignment}) -ge 4 ]
        then
            MODEL=$(\
                grep raxml-ng {input.model_file} \
                | tail -1 \
                | grep -Po "model [A-Z0-9\+\-]+" \
                | sed 's/^model //g' \
            )
            raxml-ng \
                --msa {input.alignment} \
                --model $MODEL \
                --prefix {params.prefix} \
            2> /dev/null 1>&2
        fi
        """


rule homologs_round1_trim_tips:
    input: HOMOLOGS_R1 + "{orthogroup_id}.raxml.bestTree"
    output: touch(HOMOLOGS_R1 + "{orthogroup_id}.trimmed_tips.nwk")
    params:
        relative_cutoff = params["homologs"]["trim_tips"]["relative_cutoff"],
        absolute_cutoff = params["homologs"]["trim_tips"]["absolute_cutoff"]
    threads: 1
    log: HOMOLOGS_R1 + "{orthogroup_id}.trimmed_tips.log"
    benchmark: HOMOLOGS_R1 + "{orthogroup_id}.trimmed_tips.bmk"
    conda: "homologs.yml"
    shell:
        "(python2 src/pdc/trim_tips.py "
            "{input} "
            "{output} "
            "{params.relative_cutoff} "
            "{params.absolute_cutoff} || true)"
        "2> {log} 1>&2"


rule homologs_round1_mask_tips_by_taxon_id:
    input:
        trimmed_tree = HOMOLOGS_R1 + "{orthogroup_id}.trimmed_tips.nwk",
        alignment = HOMOLOGS_R1 + "{orthogroup_id}.trimal.pep"
    output: touch(HOMOLOGS_R1 + "{orthogroup_id}.masked_tips.nwk")
    params:
        mask_tips = params["homologs"]["mask_tips"]["mask_paraphyletic"]
    threads: 1
    log: HOMOLOGS_R1 + "{orthogroup_id}.masked_tips.log",
    benchmark: HOMOLOGS_R1 + "{orthogroup_id}.masked_tips.bmk"
    conda: "homologs.yml"
    shell:
        "(python2 src/pdc/mask_tips_by_taxonID_transcripts.py "
            "{input.trimmed_tree} "
            "{output} "
            "{input.alignment} "
            "{params.mask_tips} || true)"
        "2> {log} 1>&2"


rule homologs_round1_cut_internal_long_branches:
    input:
        masked_tree = HOMOLOGS_R1 + "{orthogroup_id}.masked_tips.nwk",
        original_tree = HOMOLOGS_R1 + "{orthogroup_id}.raxml.bestTree",
    output: touch(HOMOLOGS_R1 + "{orthogroup_id}.final.nwk")
    params:
        internal_branch_cutoff = params["homologs"]["cut_internal_long_branches"]["internal_branch_cutoff"],
        minimum_taxa = params["homologs"]["cut_internal_long_branches"]["minimum_taxa"]
    threads: 1
    log: HOMOLOGS_R1 + "{orthogroup_id}.cut_internal_long_branches.log"
    benchmark: HOMOLOGS_R1 + "{orthogroup_id}.cut_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        "(python2 src/pdc/cut_long_internal_branches.py "
            "{input.masked_tree} "
            "{output} "
            "{params.internal_branch_cutoff} "
            "{params.minimum_taxa} "
            "{input.original_tree} || true)"
        "2> {log} 1>&2"


rule homologs_round1_to_fasta:
    input:
        tree = HOMOLOGS_R1 + "{orthogroup_id}.final.nwk",
        cds = HOMOLOGS_R1 + "{orthogroup_id}.trimal.cds",
        pep = HOMOLOGS_R1 + "{orthogroup_id}.trimal.pep"
    output:
        cds = HOMOLOGS_R1 + "{orthogroup_id}.final.cds",
        pep = HOMOLOGS_R1 + "{orthogroup_id}.final.pep"
    threads: 1
    log:  HOMOLOGS_R1 + "{orthogroup_id}.to_fasta.log"
    benchmark: HOMOLOGS_R1 + "{orthogroup_id}.to_fasta.bmk"
    conda: "homologs.yml"
    shell:
        "python2 src/newick_to_fasta.py "
            "{input.tree} "
            "{input.cds} "
            "{input.pep} "
            "{output.cds} "
            "{output.pep} "
        "2> {log} 1>&2"


def aggregate_homologs_round1_files(wildcards):

    checkpoint_cds = checkpoints.orthofinder_sequences.get(**wildcards).output[0]
    files_cds = expand(
        HOMOLOGS_R1 + "{i}.final.cds",
        i=glob_wildcards(os.path.join(checkpoint_cds, "{i}.cds")).i
    )

    checkpoint_pep = checkpoints.orthofinder_sequences.get(**wildcards).output[0]
    files_pep = expand(
        HOMOLOGS_R1 + "{i}.final.pep",
        i=glob_wildcards(os.path.join(checkpoint_pep, "{i}.pep")).i
    )
    return files_cds + files_pep


rule homologs_round1:
    input:
         aggregate_homologs_round1_files


################################################################################
# ROUND 2
################################################################################

rule homologs_round2_mafft:
    input: HOMOLOGS_R1 + "{orthogroup_id}.final.pep"
    output: HOMOLOGS_R2 + "{orthogroup_id}.mafft.pep"
    log: HOMOLOGS_R2 + "{orthogroup_id}.mafft.pep.log"
    benchmark: HOMOLOGS_R2 + "{orthogroup_id}.mafft.pep.bmk"
    conda: "homologs.yml"
    shell:
        """
        if [ -s {input} ]
        then
            mafft --genafpair --maxiterate 1000 {input} > {output} 2> {log}
        else
            touch {output} 2> {log} 1>&2
        fi
        """


rule homologs_round2_trimal_pep:
    input: HOMOLOGS_R2 + "{orthogroup_id}.mafft.pep"
    output: HOMOLOGS_R2 + "{orthogroup_id}.trimal.pep"
    log: HOMOLOGS_R2 + "{orthogroup_id}.trimal.log"
    benchmark: HOMOLOGS_R2 + "{orthogroup_id}.trimal.bmk"
    threads: 1
    conda: "homologs.yml"
    shell:
        """
        if [ -s {input} ]
        then
            trimal -in {input} -out {output} -automated1 2> {log} 1>&2
        else
            touch {output}
        fi
        """


rule homologs_round2_trimal_cds:
    input:
        cds_raw = OF_SEQUENCES + "{orthogroup_id}.cds",
        pep_aligned = HOMOLOGS_R2 + "{orthogroup_id}.mafft.pep"
    output: HOMOLOGS_R2 + "{orthogroup_id}.trimal.cds"
    threads: 1
    log: HOMOLOGS_R2 + "{orthogroup_id}.trimal.cds.log"
    benchmark: HOMOLOGS_R2 + "{orthogroup_id}.trimal.cds.bmk"
    conda: "homologs.yml"
    shell:
        """
        if [ -s {input.cds_raw} ] && \
            [ -s {input.pep_aligned} ]
        then
            trimal \
                -in {input.pep_aligned} \
                -backtrans {input.cds_raw} \
                -out {output} \
            2> {log} 1>&2
        else
            touch {output}
        fi
        """


rule homologs_round2_modeltestng:
    input: HOMOLOGS_R2 + "{orthogroup_id}.trimal.pep"
    output:
        ckp = HOMOLOGS_R2 + "{orthogroup_id}.modeltest-ng.ckp",
        out = HOMOLOGS_R2 + "{orthogroup_id}.modeltest-ng.out",
        tree = HOMOLOGS_R2 + "{orthogroup_id}.modeltest-ng.tree"
    params:
        prefix = HOMOLOGS_R2 + "{orthogroup_id}.modeltest-ng"
    threads: 1
    benchmark: HOMOLOGS_R2 + "{orthogroup_id}.modeltest-ng.bmk"
    conda: "homologs.yml"
    shell:
        """
        if [ -s {input} ]
        then
            modeltest-ng \
                --datatype aa \
                --input {input} \
                --template raxml \
                --out {params.prefix} \
            2> /dev/null 1>&2
        else
            touch {output}
        fi
        """


rule homologs_round2_raxmlng:
    input:
        alignment = HOMOLOGS_R2 + "{orthogroup_id}.trimal.pep",
        model_file = HOMOLOGS_R2 + "{orthogroup_id}.modeltest-ng.out"
    output:
        best_model = HOMOLOGS_R2 + "{orthogroup_id}.raxml.bestModel",
        best_tree = HOMOLOGS_R2 + "{orthogroup_id}.raxml.bestTree",
        ml_trees = HOMOLOGS_R2 + "{orthogroup_id}.raxml.mlTrees",
        rba = HOMOLOGS_R2 + "{orthogroup_id}.raxml.rba",
        start_tree = HOMOLOGS_R2 + "{orthogroup_id}.raxml.startTree"
    params:
        prefix = HOMOLOGS_R2 + "{orthogroup_id}"
    threads: 1
    log: HOMOLOGS_R2 + "{orthogroup_id}.raxml.log"
    benchmark: HOMOLOGS_R2 + "{orthogroup_id}.raxml.bmk"
    conda: "homologs.yml"
    shell:
        """
        if [ -s {input.alignment} ] && \
            [ -s {input.model_file} ] && \
            [ $(grep -c ^">" {input.alignment}) -ge 4 ]
        then
            MODEL=$(\
                grep raxml-ng {input.model_file} \
                | tail -1 \
                | grep -Po "model [A-Z0-9\+\-]+" \
                | sed 's/^model //g' \
            )
            raxml-ng \
                --msa {input.alignment} \
                --model $MODEL \
                --prefix {params.prefix} \
            2> /dev/null 1>&2
        else
            touch {output}
        fi
        """


rule homologs_round2_trim_tips:
    input: HOMOLOGS_R2 + "{orthogroup_id}.raxml.bestTree"
    output: touch(HOMOLOGS_R2 + "{orthogroup_id}.trimmed_tips.nwk")
    params:
        relative_cutoff = params["homologs"]["trim_tips"]["relative_cutoff"],
        absolute_cutoff = params["homologs"]["trim_tips"]["absolute_cutoff"]
    threads: 1
    log: HOMOLOGS_R2 + "{orthogroup_id}.trimmed_tips.log"
    benchmark: HOMOLOGS_R2 + "{orthogroup_id}.trimmed_tips.bmk"
    conda: "homologs.yml"
    shell:
        "(python2 src/pdc/trim_tips.py "
            "{input} "
            "{output} "
            "{params.relative_cutoff} "
            "{params.absolute_cutoff} "
        "|| true) "
        "2> {log} 1>&2"


rule homologs_round2_mask_tips_by_taxon_id:
    input:
        trimmed_tree = HOMOLOGS_R2 + "{orthogroup_id}.trimmed_tips.nwk",
        alignment = HOMOLOGS_R2 + "{orthogroup_id}.trimal.pep"
    output: touch(HOMOLOGS_R2 + "{orthogroup_id}.masked_tips.nwk")
    params:
        mask_tips = params["homologs"]["mask_tips"]["mask_paraphyletic"]
    threads: 1
    log: HOMOLOGS_R2 + "{orthogroup_id}.masked_tips.log",
    benchmark: HOMOLOGS_R2 + "{orthogroup_id}.masked_tips.bmk"
    conda: "homologs.yml"
    shell:
        "(python2 src/pdc/mask_tips_by_taxonID_transcripts.py "
            "{input.trimmed_tree} "
            "{output} "
            "{input.alignment} "
            "{params.mask_tips} "
        "|| true) "
        "2> {log} 1>&2"


rule homologs_round2_cut_internal_long_branches:
    input:
        masked_tree = HOMOLOGS_R2 + "{orthogroup_id}.masked_tips.nwk",
        original_tree = HOMOLOGS_R2 + "{orthogroup_id}.raxml.bestTree",
    output: touch(HOMOLOGS_R2 + "{orthogroup_id}.final.nwk")
    params:
        internal_branch_cutoff = params["homologs"]["cut_internal_long_branches"]["internal_branch_cutoff"],
        minimum_taxa = params["homologs"]["cut_internal_long_branches"]["minimum_taxa"]
    threads: 1
    log: HOMOLOGS_R2 + "{orthogroup_id}.cut_internal_long_branches.log"
    benchmark: HOMOLOGS_R2 + "{orthogroup_id}.cut_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        "(python2 src/pdc/cut_long_internal_branches.py "
            "{input.masked_tree} "
            "{output} "
            "{params.internal_branch_cutoff} "
            "{params.minimum_taxa} "
            "{input.original_tree} "
        "|| true) "
        "2> {log} 1>&2"


rule homologs_round2_to_fasta:
    input:
        tree = HOMOLOGS_R2 + "{orthogroup_id}.final.nwk",
        cds = HOMOLOGS_R2 + "{orthogroup_id}.trimal.cds",
        pep = HOMOLOGS_R2 + "{orthogroup_id}.trimal.pep"
    output:
        cds = HOMOLOGS_R2 + "{orthogroup_id}.final.cds",
        pep = HOMOLOGS_R2 + "{orthogroup_id}.final.pep"
    threads: 1
    log:  HOMOLOGS_R2 + "{orthogroup_id}.to_fasta.log"
    benchmark: HOMOLOGS_R2 + "{orthogroup_id}.to_fasta.bmk"
    conda: "homologs.yml"
    shell:
        "python2 src/newick_to_fasta.py "
            "{input.tree} "
            "{input.cds} "
            "{input.pep} "
            "{output.cds} "
            "{output.pep} "
        "2> {log} 1>&2"


def aggregate_homologs_round2_files(wildcards):

    checkpoint_cds = checkpoints.orthofinder_sequences.get(**wildcards).output[0]
    files_cds = expand(
        HOMOLOGS_R2 + "{i}.final.cds",
        i=glob_wildcards(os.path.join(checkpoint_cds, "{i}.cds")).i
    )

    checkpoint_pep = checkpoints.orthofinder_sequences.get(**wildcards).output[0]
    files_pep = expand(
        HOMOLOGS_R2 + "{i}.final.pep",
        i=glob_wildcards(os.path.join(checkpoint_pep, "{i}.pep")).i
    )
    return files_cds + files_pep


rule homologs_round2:
    input:
         aggregate_homologs_round2_files




# rule homologs_filter_1to1_orthologs:
# rule homologs_prune_paralogs_mi:
# rule homologs_prune_paralogs_mo:
# rule homologs_prune_paralogs_mt:
