from Bio import SeqIO
import sys

if len(sys.argv) != 5:
    exit("ERROR!: Enter four files:\nfilter_fasta_as_orthoMCL.py fasta_to_filter goodProteins_from_orthomcl goodProteins.fasta poorProteins.fasta")

fasta_in = sys.argv[1]
fasta_template= sys.argv[2]
fasta_good = sys.argv[3]
fasta_poor = sys.argv[4]


def filter_fasta(fasta_to_filter, fasta_template, fasta_filtered_good, fasta_filtered_poor):
    '''
    Filter a fasta file according to the presence of the IDs on a template fasta file, and write a third fasta file
    
    Input:
        - Fasta to filter
        - Fasta that contains the set of IDs that are used to filter (goodProteins.fasta from orthoMCL)
        - Good filtered CDS
        - Poor filtered CDS
    '''
    
    
    # Get the IDs from the fasta with the pattern
    template_records = SeqIO.parse(fasta_template, "fasta") 
    template_ids = {record.id for record in template_records}
    del template_records
    
    # Subset the fasta_to_filter according if the ID is in the pattern fasta
    filtered_good_records = []
    filtered_poor_records = []
    for record in SeqIO.parse(fasta_to_filter, "fasta"):
        if record.id in template_ids:
            filtered_good_records.append(record)
        else:
            filtered_poor_records.append(record)
    
    # Write the filtered files
    SeqIO.write(filtered_good_records, fasta_filtered_good, "fasta")
    SeqIO.write(filtered_poor_records, fasta_filtered_poor, "fasta")
    


# Main
filter_fasta(
    fasta_in,
    fasta_template,
    fasta_good,
    fasta_poor
    )
