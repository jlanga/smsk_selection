rule tag_pep_species:
    """
    Add to the beginning of the fasta header the species from which the fasta
    comes, i.e. for Danio rerio: >drer|seq_id
    PEP version
    """
    input: transdecoder + "{species}.pep"
    output: tag + "{species}.pep"
    params: "{species}"
    shell: "python3 src/fasta_tagger.py {params} < {input} > {output} "


rule tag_cds_species:
    """
    Add to the beginning of the fasta header the species from which the fasta
    comes, i.e. for Danio rerio: >drer|seq_id
    CDS version
    """
    input: transdecoder + "{species}.cds"
    output: tag + "{species}.cds"
    params: "{species}"
    shell: "python3 src/fasta_tagger.py {params} < {input} > {output}"


rule tag:
    """Perform all tag tasks"""
    input:
        expand(
            tag + "{species}.{extension}",
            species=SPECIES,
            extension="pep cds".split()
        )
