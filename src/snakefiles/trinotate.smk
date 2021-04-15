CHUNKS_ANNOTATE = params["trinotate"]["chunks"]






rule trinotate_blastx_assembly:
    """Run diamond blastx over assembly"""
    input:
        fa = TRANSDECODER + "transcripts.fasta",
        db = DB + "uniprot_sprot.dmnd"
    output: TRINOTATE + "blastx.tsv"
    threads: MAX_THREADS
    log: TRINOTATE + "blastx.log"
    benchmark: TRINOTATE + "blastx.bmk"
    conda: "trinotate.yml"
    shell:
        "diamond blastp "
            "--query {input.fa} "
            "--db {input.db} "
            "--outfmt 6 "
            "--evalue 1e-5 "
            "--out {output} "
            "--threads {threads} "
        "2> {log} 1>&2"


rule trinotate_split_proteome:
    """
    Split the headers from transcriptome.pep into multiple files
    """
    input:
        fai = TRANSDECODER + "transdecoder.pep.fai"
    output:
        expand(
            TRINOTATE + "chunks/proteome_{chunk_id}.tsv",
            chunk_id=[
                '{0:05d}'.format(x)
                for x in range(0, CHUNKS_ANNOTATE)
            ]
        )
    params:
        number_of_chunks = CHUNKS_ANNOTATE
    log: TRINOTATE + "split_proteome.log"
    benchmark: TRINOTATE + "split_proteome.bmk"
    conda: "trinotate.yml"
    shell:
        """
        split \
            --number l/{params.number_of_chunks} \
            --numeric-suffixes \
            --suffix-length 5 \
            --additional-suffix .tsv \
            {input.fai} \
            {TRINOTATE}/chunks/proteome_ \
        2> {log}
        """


rule trinotate_blastp_proteome:
    """Run diamond blastp over transdecoder proteome"""
    input:
        pep = TRANSDECODER + "transdecoder.pep",
        db = DB + "uniprot_sprot.dmnd"
    output: TRINOTATE + "blastp.tsv"
    threads: MAX_THREADS
    log: TRINOTATE + "blastp.log"
    benchmark: TRINOTATE + "blastp.bmk"
    conda: "trinotate.yml"
    shell:
        "diamond blastp "
            "--query {input.pep} "
            "--db {input.db} "
            "--outfmt 6 "
            "--evalue 1e-5 "
            "--out {output} "
            "--threads {threads} "
        "2> {log} 1>&2"



rule trinotate_hmmscan_chunk:
    """
    hmmscan over one chunk
    """
    input:
        pep = TRANSDECODER + "transdecoder.pep",
        fai = TRANSDECODER + "transdecoder.pep.fai",
        chunk = TRINOTATE + "chunks/proteome_{chunk_id}.tsv",
        hmm = DB + "Pfam-A.hmm"
    output:
        tsv = TRINOTATE + "hmmscan/proteome_{chunk_id}.tsv"
    log: TRINOTATE + "hmmscan/proteome_{chunk_id}.log"
    benchmark: TRINOTATE + "hmmscan/proteome_{chunk_id}.bmk"
    conda: "trinotate.yml"
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


rule trinotate_hmmscan_merge:
    """
    Merge hmmscan results into one file
    """
    input:
        expand(
            TRINOTATE + "hmmscan/proteome_{chunk_id}.tsv",
            chunk_id=[
                '{0:05d}'.format(x)
                for x in range(0, CHUNKS_ANNOTATE)
            ]
        )
    output:
        tsv = TRINOTATE + "hmmscan.tsv"
    log: TRINOTATE + "hmmscan_merge.log"
    benchmark: TRINOTATE + "hmmscan_merge.bmk"
    conda: "trinotate.yml"
    shell: "cat {input} > {output} 2> {log}"


rule trinotate_create:
    """Build sqlite database"""
    output:
        sqlite = TRINOTATE + "trinotate.sqlite",
        is_created = touch(TRINOTATE + "create.txt")
    log: TRINOTATE + "create.log"
    benchmark: TRINOTATE + "create.bmk"
    conda: "trinotate.yml"
    shell:
        """
        EMBL_dat_to_Trinotate_sqlite_resourceDB.pl \
            --sqlite {output.sqlite} \
            --create \
        2> {log} 1>&2
        """


rule trinotate_init:
    """Initialize DB with genes, transcripts and proteins"""
    input:
        sqlite = ancient(TRINOTATE + "trinotate.sqlite"),
        is_created = TRINOTATE + "create.txt",
        assembly = TRANSDECODER + "transcripts.fasta",
        proteome = TRANSDECODER + "transdecoder.pep",
        g2t = TRANSDECODER + "g2t.tsv"
    output: touch(TRINOTATE + "init.txt")
    log: TRINOTATE + "init.log"
    benchmark: TRINOTATE + "init.bmk"
    conda: "trinotate.yml"
    shell:
        """
        Trinotate {input.sqlite} init \
            --gene_trans_map {input.g2t} \
            --transcript_fasta {input.assembly} \
            --transdecoder_pep {input.proteome} \
        2> {log} 1>&2
        """


rule trinotate_fill:
    input:
        sqlite = ancient(TRINOTATE + "trinotate.sqlite"),
        is_init = TRINOTATE + "init.txt",
        eggnog = DB + "NOG.annotations.tsv.bulk_load",
        go = DB + "go-basic.obo.tab",
        uniprot_index = DB + "trinotate.UniprotIndex",
        taxonomy_index = DB + "trinotate.TaxonomyIndex",
        pfam = DB + "Pfam-A.hmm.gz.pfam_sqlite_bulk_load"
    output: touch(TRINOTATE + "fill.txt")
    log: TRINOTATE + "fill.log"
    benchmark: TRINOTATE + "fill.bmk"
    conda: "trinotate.yml"
    shell:
        """
        EMBL_dat_to_Trinotate_sqlite_resourceDB.pl \
            --sqlite {input.sqlite} \
            --eggnog {input.eggnog} \
            --go_obo_tab {input.go} \
            --uniprot_index {input.uniprot_index} \
            --taxonomy_index {input.taxonomy_index} \
            --pfam {input.pfam} \
        2> {log} 1>&2
        """


rule trinotate_load:
    input:
        sqlite = ancient(TRINOTATE + "trinotate.sqlite"),
        is_filled = TRINOTATE + "fill.txt",
        blastx = TRINOTATE + "blastx.tsv",
        blastp = TRINOTATE + "blastp.tsv",
        pfam = TRINOTATE + "hmmscan.tsv"
    output: touch(TRINOTATE + "load.txt")
    log: TRINOTATE + "load.log"
    benchmark: TRINOTATE + "load.bmk"
    conda: "trinotate.yml"
    shell:
        """
        Trinotate \
            {input.sqlite} \
            LOAD_swissprot_blastp \
            {input.blastp} \
        2> {log} 1>&2

        Trinotate \
            {input.sqlite} \
            LOAD_pfam \
            {input.pfam} \
        2>> {log} 1>&2

        Trinotate \
            {input.sqlite} \
            LOAD_swissprot_blastx \
            {input.blastx} \
        2>> {log} 1>&2
        """


rule trinotate_report:
    input:
        sqlite = ancient(TRINOTATE + "trinotate.sqlite"),
        is_loaded = TRINOTATE + "load.txt"
    output: TRINOTATE + "trinotate.tsv"
    log: TRINOTATE + "report.log"
    benchmark: TRINOTATE + "report.bmk"
    conda: "trinotate.yml"
    params:
        evalue = params["trinotate"]["evalue"],
        pfam_cutoff = params["trinotate"]["pfam_cutoff"]
    shell:
        """
        Trinotate {input.sqlite} report \
            -E {params.evalue} \
            --pfam_cutoff {params.pfam_cutoff} \
        > {output} \
        2> {log}
        """


rule trinotate:
    input: rules.trinotate_report.output
