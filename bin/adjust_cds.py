from Bio import SeqIO
import sys

# Check sys.argv

if len(sys.argv) != 4:
    exit("ERROR!: Incorrect number of arguments\nUsage: python adjust_cds.py fasta_in fasta_out code\n- is allowed for pipes!")

fasta_in  = sys.argv[1]
fasta_out = sys.argv[2]
code      = sys.argv[3]

if len(code) > 4 or len(code) < 1:
    exit("ERROR! The code should have a length between 1 to 4 characters. OrthoMCL asks for it.")

# If recibing from pipes:

if fasta_in == "-":
    fasta_in = sys.stdin

if fasta_out == "-":
    fasta_out = sys.stdout


# Main function

def create_compliant_cds(filename_in, filename_out, code):
    '''Create a fasta file with the headers modified a la orthomclAdjustFasta:
           - No description
           - "|" translated to "_" in the id fields.
           - A 1 to 4 letter code for an organism
           - A "|" between the four letter code and the id of the CDS
       Input:
           - Filename for the raw fasta
           - Filename for with the modified fasta
           - Code
       Ouptut:
           - A file will be (over)written 
    '''
    records = SeqIO.parse(filename_in,"fasta")
    new_records = []
    for record in records:
        record.description = ""
        header = record.id
        header = header.replace("|", "_")
        record.id = code + "|" + header
        new_records.append(record)
    SeqIO.write(new_records, filename_out, "fasta")


# Execute:

create_compliant_cds(fasta_in, fasta_out, code)
