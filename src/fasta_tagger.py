#!/usr/bin/env python3

if __name__ == "__main__":

    import sys

    if len(sys.argv) != 2:
        exit(
            """ERROR!: Incorrect number of arguments
            Usage: python fasta_tagger.py tag < file_in > file_out
            """
        )

# Usage: python tagger tag file_in.fasta file_out.fasta
# python tagger eenc_ eenc_contigs2.fasta eenc_contigs3.fasta

    tag   = sys.argv[1]

    for line in sys.stdin:
        if line.startswith(">"):
            sys.stdout.write(">" + tag + "|" +  line[1:])
        else:
            sys.stdout.write(line)
    sys.stdin.close()
    sys.stdout.close()
