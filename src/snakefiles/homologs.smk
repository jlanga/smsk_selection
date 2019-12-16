INTERNAL_BRANCH_CUTOFF = params["homologs"]["cut_internal_long_branches"]["internal_branch_cutoff"],
MINIMUM_TAXA = params["homologs"]["cut_internal_long_branches"]["minimum_taxa"]
MASK_PARAPHYLETIC_TIPS = params["homologs"]["mask_tips"]["mask_paraphyletic"]
MINIMUM_INGROUP_TAXA = params["homologs"]["root_to_tip"]["minimum_ingroup_taxa"]

rule homologs_join_cds:
    input: expand(CDHIT + "{species}.cds", species=SPECIES)
    output: HOMOLOGS + "all.cds"
    conda: "homologs.yml"
    shell:
        """
        cat /dev/null > {output}

        for file in {input}; do
            species=$(basename $file | sed "s/.cds//")
            seqtk seq $file \
            | cut -f 1 -d " " \
            | paste - - \
            | sed "s/^>/>$species@/" \
            | sed 's/ENA|[A-Z0-9]*|//' \
            | tr "\t" "\n" \
            >> {output}
        done
        """

rule homologs_join_pep:
    input: expand(CDHIT + "{species}.pep", species=SPECIES)
    output: HOMOLOGS + "all.pep"
    conda: "homologs.yml"
    shell:
        """
        cat /dev/null > {output}

        for file in {input}; do
            species=$(basename $file | sed "s/.pep//")
            seqtk seq $file \
            | cut -f 1 -d " " \
            | paste - - \
            | sed "s/^>/>$species@/" \
            | sed 's/ENA|[A-Z0-9]*|//' \
            | tr "\t" "\n" \
            >> {output}
        done
        """


rule homologs_correct_leafs:
    input: ORTHOFINDER + "orthologues.ok"
    output: directory(HOMOLOGS + "leafs_corrected")
    log: HOMOLOGS + "correct_leafs.log"
    benchmark: HOMOLOGS + "correct_leafs.bmk"
    conda: "homologs.yml"
    threads: MAX_THREADS
    params:
        gene_trees_folder = OF_RESOLVED_GENE_TREES
    shell:
        """
        python2 src/correct_tree_leaf_names.py \
            {params} \
            _tree.txt \
            {output} \
            .nwk \
        2> {log} 1>&2
        """


rule homologs_round1_tree_to_fasta:
    input:
        pep = HOMOLOGS + "all.pep",
        folder = HOMOLOGS + "leafs_corrected"
    output:
        folder = directory(HOMOLOGS_R1 + "fasta")
    log: HOMOLOGS_R1 + "tree_to_fasta.log"
    benchmark: HOMOLOGS_R1 + "tree_to_fasta.bmk"
    conda: "homologs.yml"
    shell:
        """
        python2 src/tree_to_fasta.py \
            {input.pep} \
            {input.folder} \
            nwk \
            {output.folder} \
            fa \
        2> {log} 1>&2
        """


rule homologs_round1_mafft:
    input: HOMOLOGS_R1 + "fasta"
    output: directory(HOMOLOGS_R1 + "mafft")
    log: HOMOLOGS_R1 + "mafft.log"
    benchmark: HOMOLOGS_R1 + "mafft.bmk"
    conda: "homologs.yml"
    shell:
        """
        bash src/mafft_folder.sh \
            {input} \
            fa \
            {output} \
            fa \
            {threads} \
        2> {log} 1>&2
        """


rule homologs_round1_pxclsq:
    input: HOMOLOGS_R1 + "mafft"
    output: directory(HOMOLOGS_R1 + "pxclsq")
    log: HOMOLOGS_R1 + "pxclsq.log"
    benchmark: HOMOLOGS_R1 + "pxclsq.bmk"
    conda: "homologs.yml"
    params:
        min_occupancy = params["homologs"]["pxclsq"]["min_occupancy"]
    shell:
        """
        PATH="bin:$PATH"

        bash src/pxclsq_folder.sh \
            {input} \
            fa \
            {output} \
            fa \
            {params.min_occupancy} \
        2> {log} 1>&2
        """


rule homologs_round1_raxmlng:
    input: HOMOLOGS_R1 + "pxclsq"
    output: directory(HOMOLOGS_R1 + "raxmlng")
    log: HOMOLOGS_R1 + "raxmlng.log"
    benchmark: HOMOLOGS_R1 + "raxmlng.bmk"
    conda: "homologs.yml"
    threads: MAX_THREADS
    shell:
        """
        bash src/raxmlng_folder.sh \
            {input} \
            fa \
            {output} \
            nwk \
            {threads} \
        2> {log} 1>&2
        """


rule homologs_round1_trim_tips:
    """Trim tips via the old method"""
    input: HOMOLOGS_R1 + "raxmlng"
    output: directory(HOMOLOGS_R1 + "trimmed_tips")
    log: HOMOLOGS_R1 + "trimmed_tips.log"
    benchmark: HOMOLOGS_R1 + "trimmed_tips.bmk"
    params:
        relative_cutoff=params["homologs"]["trim_tips"]["relative_cutoff"],
        absolute_cutoff=params["homologs"]["trim_tips"]["absolute_cutoff"]
    conda: "homologs.yml"
    shell:
        """
        mkdir -p {output}
        cp -R {input}/*.bestTree {output} 2> {log} 1>&2
        rename.ul .raxml.bestTree .tree {output}/*.bestTree


        python2.7 src/pdc2/scripts/trim_tips.py \
            {output} \
            .tree \
            {params.relative_cutoff} \
            {params.absolute_cutoff} \
        2>> {log} 1>&2
        """


rule homologs_round1_mask_tips_by_taxon_id:
    """
    Result:
    """    
    input:
        trimmed_tips = HOMOLOGS_R1 + "trimmed_tips",
        alignments = HOMOLOGS_R1 + "pxclsq"
    output: directory(HOMOLOGS_R1 + "masked_tips")
    params:
        mask_paraphyletic = MASK_PARAPHYLETIC_TIPS
    threads: 1
    log: HOMOLOGS_R1 + "masked_tips.log",
    benchmark: HOMOLOGS_R1 + "masked_tips.bmk"
    conda: "homologs.yml"
    shell:
        """
        mkdir -p {output}
        cp --recursive {input.alignments}/*.fa {output}/
        cp --recursive {input.trimmed_tips}/*.tt {output}/
        rename.ul .fa .aln-cln {output}/*.fa

        python2.7 src/pdc2/scripts/mask_tips_by_taxonID_transcripts.py \
            {output} \
            {output} \
            {params.mask_paraphyletic} \
        2> {log} 1>&2

        rm -rf {output}/*.aln-cln {output}/*.tree.tt
        """


rule homologs_round1_cut_internal_long_branches:
    """
    Result: OG\d+_\d
    """
    input: HOMOLOGS_R1 + "masked_tips"
    output: directory(HOMOLOGS_R1 + "cutted_internal_long_branches")
    params:
        internal_branch_cutoff = INTERNAL_BRANCH_CUTOFF,
        minimum_taxa = MINIMUM_TAXA
    log: HOMOLOGS_R1 + "cutted_internal_long_branches.log"
    benchmark: HOMOLOGS_R1 + "cutted_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        """
        mkdir -p {output} 
        
        python2.7 src/pdc2/scripts/cut_long_internal_branches.py \
            {input} \
            .mm \
            {params.internal_branch_cutoff} \
            {params.minimum_taxa} \
            {output} \
        2> {log} 1>&2
        """


rule homologs_round1:
    input:
        rules.homologs_round1_cut_internal_long_branches.output





################################################################################
# ROUND 2
################################################################################

rule homologs_round2_tree_to_fasta:
    input:
        pep = HOMOLOGS + "all.pep",
        folder = HOMOLOGS_R1 + "cutted_internal_long_branches"
    output:
        folder = directory(HOMOLOGS_R2 + "fasta")
    log: HOMOLOGS_R2 + "tree_to_fasta.log"
    benchmark: HOMOLOGS_R2 + "tree_to_fasta.bmk"
    conda: "homologs.yml"
    shell:
        """
        python2 src/tree_to_fasta.py \
            {input.pep} \
            {input.folder} \
            subtree \
            {output.folder} \
            fa \
        2> {log} 1>&2
        """


rule homologs_round2_mafft:
    input: HOMOLOGS_R2 + "fasta"
    output: directory(HOMOLOGS_R2 + "mafft")
    log: HOMOLOGS_R2 + "mafft.log"
    benchmark: HOMOLOGS_R2 + "mafft.bmk"
    conda: "homologs.yml"
    shell:
        """
        bash src/mafft_folder.sh \
            {input} \
            fa \
            {output} \
            fa \
            {threads} \
        2> {log} 1>&2
        """


rule homologs_round2_pxclsq:
    input: HOMOLOGS_R2 + "mafft"
    output: directory(HOMOLOGS_R2 + "pxclsq")
    log: HOMOLOGS_R2 + "pxclsq.log"
    benchmark: HOMOLOGS_R2 + "pxclsq.bmk"
    conda: "homologs.yml"
    params:
        min_occupancy = params["homologs"]["pxclsq"]["min_occupancy"]
    shell:
        """
        PATH="bin:$PATH"

        bash src/pxclsq_folder.sh \
            {input} \
            fa \
            {output} \
            fa \
            {params.min_occupancy} \
        2> {log} 1>&2
        """


rule homologs_round2_raxmlng:
    input: HOMOLOGS_R2 + "pxclsq"
    output: directory(HOMOLOGS_R2 + "raxmlng")
    log: HOMOLOGS_R2 + "raxmlng.log"
    benchmark: HOMOLOGS_R2 + "raxmlng.bmk"
    conda: "homologs.yml"
    threads: MAX_THREADS
    shell:
        """
        bash src/raxmlng_folder.sh \
            {input} \
            fa \
            {output} \
            nwk \
            {threads} \
        2> {log} 1>&2
        """


rule homologs_round2_trim_tips:
    """Trim tips via the old method"""
    input: HOMOLOGS_R2 + "raxmlng"
    output: directory(HOMOLOGS_R2 + "trimmed_tips")
    log: HOMOLOGS_R2 + "trimmed_tips.log"
    benchmark: HOMOLOGS_R2 + "trimmed_tips.bmk"
    params:
        relative_cutoff=params["homologs"]["trim_tips"]["relative_cutoff"],
        absolute_cutoff=params["homologs"]["trim_tips"]["absolute_cutoff"]
    conda: "homologs.yml"
    shell:
        """
        mkdir -p {output}
        cp -R {input}/*.bestTree {output} 2> {log} 1>&2
        rename.ul .raxml.bestTree .tree {output}/*.bestTree

        python2.7 src/pdc2/scripts/trim_tips.py \
            {output} \
            .tree \
            {params.relative_cutoff} \
            {params.absolute_cutoff} \
        2>> {log} 1>&2
        """


rule homologs_round2_mask_tips_by_taxon_id:
    """
    Result:
    """    
    input:
        trimmed_tips = HOMOLOGS_R2 + "trimmed_tips",
        alignments = HOMOLOGS_R2 + "pxclsq"
    output: directory(HOMOLOGS_R2 + "masked_tips")
    params:
        mask_paraphyletic = MASK_PARAPHYLETIC_TIPS
    log: HOMOLOGS_R2 + "masked_tips.log",
    benchmark: HOMOLOGS_R2 + "masked_tips.bmk"
    conda: "homologs.yml"
    shell:
        """
        mkdir -p {output}
        cp --recursive {input.alignments}/*.fa {output}/
        cp --recursive {input.trimmed_tips}/*.tt {output}/
        rename.ul .fa .aln-cln {output}/*.fa

        python2.7 src/pdc2/scripts/mask_tips_by_taxonID_transcripts.py \
            {output} \
            {output} \
            {params.mask_paraphyletic} \
        2> {log} 1>&2

        rm -rf {output}/*.aln-cln {output}/*.tree.tt
        """


rule homologs_round2_cut_internal_long_branches:
    """
    Result: OG\d+_\d
    """
    input: HOMOLOGS_R2 + "masked_tips"
    output: directory(HOMOLOGS_R2 + "cutted_internal_long_branches")
    params:
        internal_branch_cutoff = INTERNAL_BRANCH_CUTOFF,
        minimum_taxa = MINIMUM_TAXA
    log: HOMOLOGS_R2 + "cutted_internal_long_branches.log"
    benchmark: HOMOLOGS_R2 + "cutted_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        """
        mkdir -p {output} 
        
        python2.7 src/pdc2/scripts/cut_long_internal_branches.py \
            {input} \
            .mm \
            {params.internal_branch_cutoff} \
            {params.minimum_taxa} \
            {output} \
        2> {log} 1>&2
        """


rule homologs_round2:
    input:
        rules.homologs_round2_cut_internal_long_branches.output


##############################################################################
## Prune paralogs
##############################################################################

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

# rule homologs_rt_prepare:
#     input: rules.homologs_round2.input
#     output: touch(HOMOLOGS + "rt_prepare.ok")
#     params:
#         in_dir = HOMOLOGS_R2,
#         out_dir = HOMOLOGS_RT
#     log: HOMOLOGS + "rt_prepare.log"
#     benchmark: HOMOLOGS + "rt_prepare.bmk"
#     shell:
#         """
#         mkdir -p {params.out_dir}
#         ln {params.in_dir}/*.subtree {params.out_dir}
#         """


rule homologs_rt_prune_paralogs:
    input:
        tsv = HOMOLOGS + "in_out.tsv",
        tree_dir = HOMOLOGS_R2 + "cutted_internal_long_branches"
    output: directory(HOMOLOGS_RT)
    params:
        minimum_ingroup_taxa = MINIMUM_INGROUP_TAXA
    log: HOMOLOGS + "root_to_tip.log"
    benchmark: HOMOLOGS + "root_to_tip.bmk"
    conda: "homologs.yml"
    shell:
        """
        mkdir -p {output}

        python2.7 src/pdc2/scripts/prune_paralogs_RT.py \
            {input.tree_dir} \
            subtree \
            {output} \
            {params.minimum_ingroup_taxa} \
            {input.tsv} \
        2> {log} 1>&2
        """

rule homologs_rt:
    input: HOMOLOGS_RT

###############################################################################
## Refinement 1
###############################################################################



rule homologs_refine1:
    input: 
        in_dir = HOMOLOGS_RT,
        pep = HOMOLOGS + "all.pep",
        cds = HOMOLOGS + "all.cds"
    output: directory(HOMOLOGS_REFINE1)
    log: HOMOLOGS + "refine1.log"
    benchmark: HOMOLOGS + "refine1.bmk"
    threads: MAX_THREADS
    conda: "homologs.yml"
    shell:
        """
        bash src/tree_to_fasta.sh \
            {input.pep} \
            {input.in_dir} \
            tre \
            {output} \
            fa \
        2> {log} 1>&2

        python2.7 src/refine_alignments2.py \
            {HOMOLOGS_REFINE1} \
            fa \
            {input.cds} \
            {threads} \
        2> {log} 1>&2

        rm -rf *.dnd
        """


rule homologs_refine2:
    input: 
        in_dir = HOMOLOGS_REFINE1,
        cds = HOMOLOGS + "all.cds"
    output: directory(HOMOLOGS_REFINE2)
    log: HOMOLOGS + "refine2.log"
    benchmark: HOMOLOGS + "refine2.bmk"
    threads: MAX_THREADS
    conda: "homologs.yml"
    shell:
        """
        cp {input.in_dir}/*.maxalign.fa {output}
        rename.ul .maxalign.fa .fa {output}/*.maxalign.fa

        python2.7 src/refine_alignments2.py \
            {HOMOLOGS_REFINE2} \
            fa \
            {input.cds} \
            {threads} \
        2> {log} 1>&2

        rm -rf *.dnd
        """


rule homologs:
    input: HOMOLOGS_REFINE2
