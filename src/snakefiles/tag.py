rule tag_pep_species:
    """
    Add to the beginning of the fasta header the species from which the fasta
    comes, i.e. for Danio rerio: >drer|seq_id
    PEP version
    """
    input: TRANSDECODER + "{species}.pep"
    output: TAG + "{species}.pep"
    params: "{species}"
    conda: "tag.yml"
    shell: "python3 src/fasta_tagger.py {params} < {input} > {output} "


rule tag_cds_species:
    """
    Add to the beginning of the fasta header the species from which the fasta
    comes, i.e. for Danio rerio: >drer|seq_id
    CDS version
    """
    input: TRANSDECODER + "{species}.cds"
    output: TAG + "{species}.cds"
    params: "{species}"
    conda: "tag.yml"
    shell: "python3 src/fasta_tagger.py {params} < {input} > {output}"


rule tag:
    """Perform all TAG tasks"""
    input:
        expand(
            TAG + "{species}.{extension}",
            species=SPECIES,
            extension="pep cds".split()
        )
