#!/usr/bin/env python3

'''
Removes stop codons (*) from the end of a sequence in fasta format.
python3 remove_stop_codons.py < in.fasta > out.fasta
'''

import sys  # To send error messages, take from STDIN and drop to STDOUT

__author__ = "Jorge Langa"
__copyright__ = "Copyright 2016, Jorge Langa"
__credits__ = ["Jorge Langa"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Jorge Langa"
__email__ = "jorge.langa.arranz@gmail.com"
__status__ = "Prototype"


if __name__ == "__main__":

    # Check CLI arguments
    if len(sys.argv) != 1:
        usage = "{} < in.fasta > out.fasta".format(sys.argv[0])
        exit("Incorrect number of inputs:\n" + usage)
    # Main execution - Convert
    for line in sys.stdin:
        sys.stdout.write(line.replace("*\n", "\n"))
