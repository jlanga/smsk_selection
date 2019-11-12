rule cdhit_pep:
    input: TRANSDECODER + "{species}.pep"
    output: CDHIT + "{species}.pep"
    threads: 1
    benchmark: CDHIT + "{species}.cdhit.bmk"
    log: CDHIT + "{species}.cdhit.log"
    params:
        threshold = params["cd-hit"]["threshold"]
    conda: "cdhit.yml"
    shell:
        "cd-hit "
            "-i {input} "
            "-o {output} "
            "-c {params.threshold} "
            "-n 5 "  # Word length
            "-T {threads} "
        "2> {log} 1>&2"


rule cdhit_filter_cds:
    input:
        cds = TRANSDECODER + "{species}.cds",
        fai = CDHIT + "{species}.pep.fai"
    output:
        cds = CDHIT + "{species}.cds"
    threads: 1
    benchmark: CDHIT + "{species}.cds.bmk"
    log: CDHIT + "{species}.cds.log"
    conda: "cdhit.yml"
    shell:
        "seqtk subseq "
            "{input.cds} "
            "<(cut -f 1 {input.fai}) "
        "> {output.cds} "
        "2> {log}"


rule cdhit:
    input:
        expand(
            CDHIT + "{species}.pep",
            species=SPECIES
        ),
        expand(
            CDHIT + "{species}.cds",
            species=SPECIES
        )