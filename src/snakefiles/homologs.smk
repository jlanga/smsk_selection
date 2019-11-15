rule homologs_merge_pep:
    input: aggregate_orthofinder_pep
    output: HOMOLOGS + "all.pep"
    shell: "cat {input} > {output}"


rule homologs_merge_cds:
    input: aggregate_orthofinder_cds
    output: HOMOLOGS + "all.cds"
    shell: "cat {input} > {output}"


rule homologs_round1_prepare:
    input:
        OF_SEQUENCES + "{orthogroup_id}.pep"
    output:
        HOMOLOGS_R1 + "{orthogroup_id}.fa"
    shell:
        "ln -s $(readlink -f {input}) $(readlink -f {output})"


def aggregate_homologs_round1_prepare(wildcards):
    checkpoint_pep = checkpoints.orthofinder_sequences.get(**wildcards).output[0]
    files_pep = expand(
        HOMOLOGS_R1 + "{i}.fa",
        i=glob_wildcards(os.path.join(checkpoint_pep, "{i}.pep")).i
    )
    return files_pep


rule homologs_round1_fasta_to_tree:
    input: aggregate_homologs_round1_prepare
    output:
        touch(HOMOLOGS + "round1_fasta_to_tree.ok")
    params:
        in_dir = HOMOLOGS_R1
    log: HOMOLOGS + "round1_fasta_to_tree.log"
    benchmark: HOMOLOGS + "round1_fasta_to_tree.bmk"
    threads: 32
    conda: "homologs.yml"
    shell:
        """
        PATH="bin:$PATH"
        python src/pdc3/scripts/fasta_to_tree_pxclsq.py \
            {params.in_dir} \
            {threads} \
            aa \
            y \
        2> {log} 1>&2
        """

rule homologs_round1_treeshrink:
    input: HOMOLOGS + "round1_fasta_to_tree.ok"
    output: touch(HOMOLOGS + "round1_treeshrink.ok")
    log: HOMOLOGS + "round1_treeshrink.log"
    benchmark: HOMOLOGS + "round1_treeshrink.bmk"
    threads: 32
    params:
        in_dir = HOMOLOGS_R1,
        quantile = 0.05
    conda: "homologs.yml"
    shell:
        """
        PATH="bin:$PATH"
        python src/pdc3/scripts/tree_shrink_wrapper.py \
            {params.in_dir} \
            .tre \
            {params.quantile} \
            {params.in_dir} \
        2> {log} 1>&2
        """


rule homologs_round1_mask_tips_by_taxon_id:
    input: HOMOLOGS + "round1_treeshrink.ok"
    output: touch(HOMOLOGS + "round1_mask_tips_by_taxon_id.ok")
    params:
        in_dir = HOMOLOGS_R1,
        mask_tips = params["homologs"]["mask_tips"]["mask_paraphyletic"]
    threads: 1
    log: HOMOLOGS + "round1_mask_tips_by_taxon_id.log",
    benchmark: HOMOLOGS + "round1_mask_tips_by_taxon_id.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/mask_tips_by_taxonID_transcripts.py "
            "{params.in_dir} "
            "{params.in_dir} "
            "{params.mask_tips} "
        "2> {log} 1>&2"


rule homologs_round1_cut_internal_long_branches:
    input: HOMOLOGS + "round1_mask_tips_by_taxon_id.ok"
    output: touch(HOMOLOGS + "round1_cut_internal_long_branches.ok")
    params:
        in_dir = HOMOLOGS_R1,
        internal_branch_cutoff = params["homologs"]["cut_internal_long_branches"]["internal_branch_cutoff"],
        minimum_taxa = params["homologs"]["cut_internal_long_branches"]["minimum_taxa"]
    threads: 1
    log: HOMOLOGS + "round1_cut_internal_long_branches.log"
    benchmark: HOMOLOGS + "homologs_round1_cut_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/cut_long_internal_branches.py "
            "{params.in_dir} "
            ".mm "
            "{params.internal_branch_cutoff} "
            "{params.minimum_taxa} "
            "{params.in_dir} "
        "2> {log} 1>&2"


rule homologs_round1_write_fasta_files_from_trees:
    input:
        fasta = HOMOLOGS + "all.pep",
        ok = HOMOLOGS + "round1_cut_internal_long_branches.ok"
    output:
        touch(HOMOLOGS + "round1_write_fasta_files_from_trees.ok")
    params:
        indir = HOMOLOGS_R1,
        outdir = HOMOLOGS_R2
    log: HOMOLOGS + "round1_write_fasta_files_from_trees.log"
    benchmark: HOMOLOGS + "round1_write_fasta_files_from_trees.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/write_fasta_files_from_trees.py "
            "{input.fasta} "
            "{params.indir} "
            ".subtree "
            "{params.indir} "
        "2> {log} 1>&2"


rule homologs_round1:
    input:
        rules.homologs_round1_write_fasta_files_from_trees.output


################################################################################
# ROUND 2
################################################################################

rule homologs_round2_prepare:
    input: HOMOLOGS + "round1_write_fasta_files_from_trees.ok"
    output: touch(HOMOLOGS + "round2_prepare.ok")
    shell:
        "mkdir -p {HOMOLOGS_R2}; "
        "ln $(readlink -f {HOMOLOGS_R1}/*rr.fa) {HOMOLOGS_R2}"

rule homologs_round2_fasta_to_tree:
    input: HOMOLOGS + "round2_prepare.ok"
    output:
        touch(HOMOLOGS + "round2_fasta_to_tree.ok")
    params:
        in_dir = HOMOLOGS_R2
    log: HOMOLOGS + "round2_fasta_to_tree.log"
    benchmark: HOMOLOGS + "round2_fasta_to_tree.bmk"
    threads: 32
    conda: "homologs.yml"
    shell:
        """
        PATH="bin:$PATH"
        python src/pdc3/scripts/fasta_to_tree_pxclsq.py \
            {params.in_dir} \
            {threads} \
            aa \
            y \
        2> {log} 1>&2
        """

rule homologs_round2_treeshrink:
    input: HOMOLOGS + "round2_fasta_to_tree.ok"
    output: touch(HOMOLOGS + "round2_treeshrink.ok")
    log: HOMOLOGS + "round2_treeshrink.log"
    benchmark: HOMOLOGS + "round2_treeshrink.bmk"
    threads: 32
    params:
        in_dir = HOMOLOGS_R2,
        quantile = 0.05
    conda: "homologs.yml"
    shell:
        """
        PATH="bin:$PATH"
        python src/pdc3/scripts/tree_shrink_wrapper.py \
            {params.in_dir} \
            .tre \
            {params.quantile} \
            {params.in_dir} \
        2> {log} 1>&2
        """


rule homologs_round2_mask_tips_by_taxon_id:
    input: HOMOLOGS + "round2_treeshrink.ok"
    output: touch(HOMOLOGS + "round2_mask_tips_by_taxon_id.ok")
    params:
        in_dir = HOMOLOGS_R2,
        mask_tips = params["homologs"]["mask_tips"]["mask_paraphyletic"]
    threads: 1
    log: HOMOLOGS + "round2_mask_tips_by_taxon_id.log",
    benchmark: HOMOLOGS + "round2_mask_tips_by_taxon_id.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/mask_tips_by_taxonID_transcripts.py "
            "{params.in_dir} "
            "{params.in_dir} "
            "{params.mask_tips} "
        "2> {log} 1>&2"


rule homologs_round2_cut_internal_long_branches:
    input: HOMOLOGS + "round2_mask_tips_by_taxon_id.ok"
    output: touch(HOMOLOGS + "round2_cut_internal_long_branches.ok")
    params:
        in_dir = HOMOLOGS_R2,
        internal_branch_cutoff = params["homologs"]["cut_internal_long_branches"]["internal_branch_cutoff"],
        minimum_taxa = params["homologs"]["cut_internal_long_branches"]["minimum_taxa"]
    threads: 1
    log: HOMOLOGS + "round2_cut_internal_long_branches.log"
    benchmark: HOMOLOGS + "homologs_round2_cut_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/cut_long_internal_branches.py "
            "{params.in_dir} "
            ".mm "
            "{params.internal_branch_cutoff} "
            "{params.minimum_taxa} "
            "{params.in_dir} "
        "2> {log} 1>&2"


checkpoint homologs_round2_write_fasta_files_from_trees:
    input:
        fasta = HOMOLOGS + "all.pep",
        ok = HOMOLOGS + "round2_cut_internal_long_branches.ok"
    output:
        touch(HOMOLOGS + "round2_write_fasta_files_from_trees.ok")
    params:
        indir = HOMOLOGS_R2,
        outdir = HOMOLOGS_R2
    log: HOMOLOGS + "round2_write_fasta_files_from_trees.log"
    benchmark: HOMOLOGS + "round2_write_fasta_files_from_trees.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/write_fasta_files_from_trees.py "
            "{input.fasta} "
            "{params.indir} "
            ".subtree "
            "{params.indir} "
        "2> {log} 1>&2"


rule homologs_round2:
    input:
        rules.homologs_round2_write_fasta_files_from_trees.output

def write_inout(samples, output):
    samples\
        .reset_index()\
        [["inout", "species"]]\
        .to_csv(
            path_or_buf=output,
            sep="\t",
            header=False,
            index=False
        )


rule homologs_create_taxa_inout:
    output: HOMOLOGS + "in_out.tsv"
    log: HOMOLOGS + "in_out.log"
    benchmark: HOMOLOGS + "in_out.bmk"
    run:
        samples\
            .reset_index()\
            [["inout", "species"]]\
            .to_csv(
                path_or_buf=output[0],
                sep='\t',
                index=False,
                header=False
            )

rule homologs_prepare_rt:
    input: rules.homologs_round2.input
    output: touch(HOMOLOGS + "prepare_rt.ok")
    params:
        in_dir = HOMOLOGS_R2,
        out_dir = HOMOLOGS_RT
    log: HOMOLOGS + "prepare_rt.log"
    benchmark: HOMOLOGS + "prepare_rt.bmk"
    shell:
        """
        mkdir -p {params.out_dir}
        ln {params.in_dir}/*.tre {params.out_dir}
        """


rule homologs_prune_paralogs_rt:
    input:
        tsv = HOMOLOGS + "in_out.tsv",
        ok = HOMOLOGS + "prepare_rt.ok"
    output:
        ok = touch(HOMOLOGS + "prune_paralogs_rt.ok")
    params:
        in_dir = HOMOLOGS_R2,
        tree_ending = ".tre",
        out_dir = HOMOLOGS_RT,
        minimum_ingroup_taxa = params["homologs"]["root_to_tip"]["minimum_ingroup_taxa"]
    log: HOMOLOGS + "prune_paralogs_rt.log"
    benchmark: HOMOLOGS + "prune_paralogs_rt.bmk"
    conda: "homologs.yml"
    shell:
        """
        python src/pdc3/scripts/prune_paralogs_RT.py \
            {params.in_dir} \
            {params.tree_ending} \
            {params.out_dir} \
            {params.minimum_ingroup_taxa} \
            {input.tsv} \
        2> {log} 1>&2
        """

rule homologs_rt:
    input: HOMOLOGS + "prune_paralogs_rt.ok"



