CHUNKS = params["transdecoder"]["chunks"]

rule transdecoder_transcripts:
    input: HOMOLOGS + "all.cds"
    output: TRANSDECODER + "transcripts.fasta"
    log: TRANSDECODER + "transcripts.log"
    benchmark: TRANSDECODER + "transcripts.bmk"
    conda: "transdecoder.yml"
    shell: "sed \'s/\.p[0-9]*$//g\' < {input} > {output} 2> {log}"



rule transdecoder_g2t:
    input: TRANSDECODER + "transcripts.fasta"
    output: TRANSDECODER + "g2t.tsv"
    log: TRANSDECODER + "g2t.log"
    conda: "transdecoder.yml"
    shell: 
        """
        grep ^'>' {input} \
        | tr -d '>' \
        | awk \'{{print $1 "\t" $1}}\' \
        > {output}
        """



rule transdecoder_longorfs:
    """
    Predict proteins according to the presence of long ORFs
    """
    input:
        fasta = TRANSDECODER + "transcripts.fasta",
        tsv = TRANSDECODER + "g2t.tsv"
    output:
        "transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    log:
        TRANSDECODER + "longorfs.log"
    benchmark:
        TRANSDECODER + "longorfs.bmk"
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
        fai = "transcripts.fasta.transdecoder_dir/longest_orfs.pep.fai"
    output:
        expand(
            TRANSDECODER + "chunks/longest_orfs_{chunk_id}.tsv",
            chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS)]
        )
    params:
        number_of_chunks = CHUNKS
    log:
        TRANSDECODER + "split_longest_orfs.log"
    benchmark:
        TRANSDECODER + "split_longest_orfs.bmk"
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
            {TRANSDECODER}/chunks/longest_orfs_ \
        2> {log}
        """


rule transdecoder_hmmscan_chunk:
    """
    hmmscan over one chunk
    """
    input:
        pep = "transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        fai = "transcripts.fasta.transdecoder_dir/longest_orfs.pep.fai",
        chunk = TRANSDECODER + "chunks/longest_orfs_{chunk_id}.tsv",
        hmm = DB + "Pfam-A.hmm"
    output:
        tsv = TRANSDECODER + "hmmscan/longest_orfs_{chunk_id}.tsv"
    log:
        TRANSDECODER + "hmmscan/longest_orfs_{chunk_id}.log"
    benchmark:
        TRANSDECODER + "hmmscan/longest_orfs_{chunk_id}.bmk"
    conda:
        "transdecoder.yml"
    shell:
        """
        cut -f 1 {input.chunk} \
        | xargs samtools faidx {input.pep} \
        | hmmscan \
            --cpu {threads} \
            --domtblout {output.tsv} \
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
            TRANSDECODER + "hmmscan/longest_orfs_{chunk_id}.tsv",
            chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS)]
        )
    output:
        tsv = TRANSDECODER + "hmmscan.tsv"
    log:
        TRANSDECODER + "hmmscan_merge.log"
    benchmark:
        TRANSDECODER + "hmmscan_merge.bmk"
    conda:
        "transdecoder.yml"
    shell:
        """cat {input} > {output} 2> {log}"""


rule transdecoder_blastp:
    input:
        pep = "transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        db = DB + "uniprot_sprot.dmnd"
    output: TRANSDECODER + "blastp.tsv"
    threads: MAX_THREADS
    log: TRANSDECODER + "blastp.log"
    benchmark: TRANSDECODER + "blastp.bmk"
    conda: "transdecoder.yml"
    shell:
        "diamond blastp "
            "--query {input.pep} "
            "--db {input.db} "
            "--outfmt 6 "
            "--evalue 1e-5 "
            "--out {output} "
            "--threads {threads} "
        "2> {log} 1>&2"


rule transdecoder_predict:
    """
    Join results from blast and hmmr to predict coding sequences
    """
    input:
        fasta = TRANSDECODER + "transcripts.fasta",
        pfam_tsv = TRANSDECODER + "hmmscan.tsv",
        blastp_tsv = TRANSDECODER + "blastp.tsv"
    output:
        bed = TRANSDECODER + "transdecoder.bed",
        cds = TRANSDECODER + "transdecoder.cds",
        gff3 = TRANSDECODER + "transdecoder.gff3",
        pep = TRANSDECODER + "transdecoder.pep",
    params:
        dir = "transcripts.fasta.transdecoder_dir",
        bed = "transcripts.fasta.transdecoder.bed",
        cds = "transcripts.fasta.transdecoder.cds",
        gff3 = "transcripts.fasta.transdecoder.gff3",
        pep = "transcripts.fasta.transdecoder.pep",
        checkpoints = "transcripts.fasta.transdecoder_dir.__checkpoints"
    threads:
        24
    log:
        TRANSDECODER + "predict.log"
    benchmark:
        TRANSDECODER + "predict.bmk"
    conda:
        "transdecoder.yml"
    shell:
        """
        TransDecoder.Predict \
            -t {input.fasta} \
            --retain_pfam_hits {input.pfam_tsv} \
            --retain_blastp_hits {input.blastp_tsv} \
            --cpu {threads} \
        2> {log} 1>&2

        mv {params.bed} {output.bed} 2>> {log} 1>&2
        mv {params.cds} {output.cds} 2>> {log} 1>&2
        mv {params.gff3} {output.gff3} 2>> {log} 1>&2
        mv {params.pep} {output.pep} 2>> {log} 1>&2
        rm -rf {params.dir} {params.checkpoints} 2>> {log} 1>&2
        rm -rf pipeliner.*
        """


rule transdecoder:
    input: rules.transdecoder_predict.output
