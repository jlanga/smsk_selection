rule generic_index:
    input: "{prefix}"
    output: "{prefix}.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"