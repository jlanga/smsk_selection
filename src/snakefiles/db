rule db_makeblastdb_uniref90:
    input:
        fa_gz = download + "uniref90.fa.gz"
    output:
        db = touch(
            db + "uniref90"
            )
    threads:
        1
    log:
        db + "makeblastdb_uniref90.log"
    benchmark:
        db + "makeblastdb_uniref90.json"
    shell:
        "gzip "
            "--decompress "
            "--stdout "
            "{input.fa_gz} | "
        "makeblastdb "
            "-dbtype prot "
            "-title {output.db} "
            "-out {output.db} "
            "-parse_seqids "
        "2> {log} 1>&2"



rule db_makeblastdb_swissprot:
    input:
        fa_gz = download + "swissprot.fa.gz"
    output:
        db = touch(
            db + "swissprot"
            )
    threads:
        1
    log:
        db + "makeblastdb_swissprot.log"
    benchmark:
        db + "makeblastdb_swissprot.json"
    shell:
        "gzip "
            "--decompress "
            "--stdout "
            "{input.fa_gz} | "
        "makeblastdb "
            "-dbtype prot "
            "-title {output.db} "
            "-out {output.db} "
            "-parse_seqids "
        "2> {log} 1>&2"




rule db_hmmpress_pfama: ##data/db/Pfam-A.hmm - Pfam database
    input:
        hmm_gz = download + "pfama.hmm.gz"
    output:
        hmm = db + "pfama.hmm" # Be careful! The database and the original hmm file share the same name!
    threads:
        1
    log:
        db + "hmmpress_pfama.log"
    benchmark:
        db + "hmmpress_pfama.json"
    shell:
        "pigz --decompress --keep --stdout {input.hmm_gz} > {output.hmm} 2> {log} ; "
        "hmmpress {output.hmm} 2>> {log} 1>&2"
