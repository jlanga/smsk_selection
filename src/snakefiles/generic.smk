rule index_pep:
    input: "{prefix}.pep"
    output: "{prefix}.pep.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"

rule index_fasta:
    input: "{prefix}.fasta"
    output: "{prefix}.fasta.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"

rule index_cds:
    input: "{prefix}.cds"
    output: "{prefix}.cds.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"

rule index_fa:
    input: "{prefix}.fa"
    output: "{prefix}.fa.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"

rule index_faa:
    input: "{prefix}.faa"
    output: "{prefix}.faa.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"

rule index_fna:
    input: "{prefix}.fna"
    output: "{prefix}.fna.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"
