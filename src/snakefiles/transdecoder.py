CHUNKS = params["transdecoder"]["chunks"]

rule transdecoder_longorfs:
    """
    Predict ORFs by length
    """
    input:
        fasta = RAW + "{species}.fasta",
        tsv = RAW + "{species}.g2t.tsv"
    output:
        "{species}.fasta.transdecoder_dir/longest_orfs.pep"
    log:
        TRANSDECODER + "{species}.longorfs.log"
    benchmark:
        TRANSDECODER + "{species}.longorfs.json"
    conda:
        "transdecoder.yml"
    shell:
        """
        TransDecoder.LongOrfs \
            -t {input.fasta} \
            --gene_trans_map {input.tsv} \
        2> {log} 1>&2
        """


rule transdecoder_split_longest_orfs:
    """
    Split the headers from transdecoder_longest_orfs into multiple files
    """
    input:
        fai = "{species}.fasta.transdecoder_dir/longest_orfs.pep.fai"
    output:
        expand(
            TRANSDECODER + "{species}/chunks/longest_orfs_{chunk_id}.tsv",
            species="{species}",
            chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS)]
        )
    params:
        number_of_chunks = CHUNKS,
        prefix = TRANSDECODER + "{species}/chunks/longest_orfs_"
    log:
        TRANSDECODER + "{species}/split_longest_orfs.log"
    benchmark:
        TRANSDECODER + "{species}/split_longest_orfs.json"
    conda:
        "transdecoder.yml"
    shell:
        """
        split \
            --number l/{params.number_of_chunks} \
            --numeric-suffixes \
            --suffix-length 5 \
            --additional-suffix .tsv \
            {input.fai} \
            {params.prefix} \
        2> {log}
        """

rule transdecoder_hmmscan:
    """
    hmmscan over one chunk
    """
    input:
        pep = "{species}.fasta.transdecoder_dir/longest_orfs.pep",
        fai = "{species}.fasta.transdecoder_dir/longest_orfs.pep.fai",
        chunk = TRANSDECODER + "{species}/chunks/longest_orfs_{chunk_id}.tsv",
        hmm = DB + "Pfam-A.hmm"
    output:
        tsv = TRANSDECODER + "{species}/hmmscan/longest_orfs_{chunk_id}.tsv"
    log:
        TRANSDECODER + "{species}/hmmscan/longest_orfs_{chunk_id}.log"
    benchmark:
        TRANSDECODER + "{species}/hmmscan/longest_orfs_{chunk_id}.json"
    conda:
        "transdecoder.yml"
    shell:
        """
        cut -f 1 {input.chunk} \
        | xargs samtools faidx {input.pep} \
        | hmmscan \
            --domtblout {output.tsv} \
            --cpu {threads} \
            {input.hmm} \
            - \
        2> {log} 1>&2
        """


rule transdecoder_hmmscan_merge:
    """
    Merge hmmscan results into one file
    """
    input:
        expand(
            TRANSDECODER + "{species}/hmmscan/longest_orfs_{chunk_id}.tsv",
            species="{species}",
            chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS)]
        )
    output:
        tsv = TRANSDECODER + "{species}/hmmscan.tsv"
    log:
        TRANSDECODER + "{species}/hmmscan_merge.log"
    benchmark:
        TRANSDECODER + "{species}/hmmscan_merge.json"
    conda:
        "transdecoder.yml"
    shell:
        """cat {input} > {output} 2> {log}"""



rule transdecoder_diamond_chunk:
    """
    Run blastp of each chunk
    """
    input:
        pep = "{species}.fasta.transdecoder_dir/longest_orfs.pep",
        fai = "{species}.fasta.transdecoder_dir/longest_orfs.pep.fai",
        chunk = TRANSDECODER + "{species}/chunks/longest_orfs_{chunk_id}.tsv",
        db = DB + "uniprot_sprot.dmnd"
    output:
        tsv = TRANSDECODER + "{species}/blastp/longest_orfs_{chunk_id}.tsv"
    log:
        TRANSDECODER + "{species}/blastp/longest_orfs_{chunk_id}.log"
    benchmark:
        TRANSDECODER + "{species}/blastp/longest_orfs_{chunk_id}.json"
    conda:
        "transdecoder.yml"
    shell:
        """
        cut -f 1 {input.chunk} \
        | xargs samtools faidx {input.pep} \
        | diamond blastp \
            --db {input.db} \
            --max-target-seqs 1 \
            --outfmt 6 \
            --evalue 1e-5 \
            --out {output.tsv} \
            --threads {threads} \
        2> {log} 1>&2
        """



rule transdecoder_blastp_merge:
    """
    Merge results from the different blastps
    """
    input:
        expand(
            TRANSDECODER + "{species}/blastp/longest_orfs_{chunk_id}.tsv",
            species="{species}",
            chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS)]
        )
    output:
        tsv = TRANSDECODER + "{species}/blastp.tsv"
    log:
        TRANSDECODER + "{species}/blastp_merge.log"
    benchmark:
        TRANSDECODER + "{species}/blastp_merge.json"
    conda:
        "transdecoder.yml"
    shell:
        "cat {input} > {output} 2> {log}"


rule transdecoder_predict:
    """
    Join results from blast and hmmr to predict coding sequences
    """
    input:
        fasta = RAW + "{species}.fasta",
        pfam_tsv = TRANSDECODER + "{species}/hmmscan.tsv",
        blastp_tsv = TRANSDECODER + "{species}/blastp.tsv"
    output:
        bed = TRANSDECODER + "{species}.bed",
        cds = TRANSDECODER + "{species}.cds",
        gff3 = TRANSDECODER + "{species}.gff3",
        pep = TRANSDECODER + "{species}.pep",
    params:
        folder = "{species}.fasta.transdecoder_dir",
        bed = "{species}.fasta.transdecoder.bed",
        cds = "{species}.fasta.transdecoder.cds",
        gff3 = "{species}.fasta.transdecoder.gff3",
        pep = "{species}.fasta.transdecoder.pep",
        checkpoints = "{species}.fasta.transdecoder_dir.__checkpoints"
    threads:
        24
    log:
        TRANSDECODER + "{species}/predict.log"
    benchmark:
        TRANSDECODER + "{species}/predict.json"
    conda:
        "transdecoder.yml"
    shell:
        """
        TransDecoder.Predict \
            -t {input.fasta} \
            --retain_pfam_hits {input.pfam_tsv} \
            --retain_blastp_hits {input.blastp_tsv} \
            --cpu {threads} \
            --no_refine_starts \
        2> {log} 1>&2

        mv {params.bed} {output.bed} 2>> {log} 1>&2
        mv {params.cds} {output.cds} 2>> {log} 1>&2
        mv {params.gff3} {output.gff3} 2>> {log} 1>&2
        mv {params.pep} {output.pep} 2>> {log} 1>&2
        rm -rf {params.folder} {params.checkpoints} 2>> {log} 1>&2
        """


rule transdecoder:
    input:
        expand(
            TRANSDECODER + "{species}.{ending}",
            species=SPECIES,
            ending="pep cds bed gff3".split()
        )
