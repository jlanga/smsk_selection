rule tag_pep_species:
    """
    Add to the beginning of the fasta header the species from which the fasta
    comes, i.e. for Danio rerio: >drer|seq_id
    PEP version
    """
    input:
        pep = transdecoder + "{species}.pep"
    output:
        pep = tag + "{species}.pep"
    threads:
        1
    params:
        tag = "{species}"
    log:
        tag + "pep_{species}.log"
    benchmark:
        tag + "pep_{species}.json"
    shell:
        "python3 bin/fasta_tagger.py "
            "{params.tag} "
        "< {input.pep} "
        "> {output.pep} "
        "2> {log}"



rule tag_cds_species:
    """
    Add to the beginning of the fasta header the species from which the fasta
    comes, i.e. for Danio rerio: >drer|seq_id
    CDS version
    """
    input:
        cds = transdecoder + "{species}.cds"
    output:
        cds = tag + "{species}.cds"
    threads:
        1
    params:
        tag = "{species}"
    log:
        tag + "cds_{species}.log"
    benchmark:
        tag + "cds_{species}.json"
    shell:
        "python3 bin/fasta_tagger.py "
            "{params.tag} "
        "< {input.cds} "
        "> {output.cds} "
        "2> {log}"
