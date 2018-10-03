


rule fastcodeml


# Extract single copy orthologs. Remember that we have in the array FISHES the codes of each species:
find data/groups_protein -delete
python scripts/extract_single_copy_orthologs.py             \
    --out_dir   data/groups_protein                         \
    --infasta   data/proteomes_filtered/goodProteins.fasta  \
    --groups    data/groups/groups.txt                      \
    --species   ${FISHES[*]}



# Remove stop codons
mkdir -p data/groups_protein_no_stop

parallel                                        \
    python scripts/remove_stop_codons.py        \
        data/groups_protein/{/.}.fasta          \
        data/groups_protein_no_stop/{/.}.fasta  \
::: $(find data/groups_protein -name "group_?????.fasta" | sort )

mv data/groups_protein_no_stop/* data/groups_protein/
rmdir data/groups_protein_no_stop



## Use adjust_cds.py to create a fasta with the header fields as orthomclAdjustFasta does
mkdir -p data/filtered_cds
for sample in ${FISHES[*]} ; do
    python scripts/adjust_cds.py            \
        data/proteomes/${sample}.cds        \
        data/filtered_cds/${sample}.cds    \
        $sample
done


# Put all the cds together
cat data/filtered_cds/{char,eenc,erin,ssco,tthy}.cds \
> data/filtered_cds/all.cds
rm data/filtered_cds/{char,eenc,erin,ssco,tthy}.cds

# Filter all.cds according to goodProteins and poorProteins:
python scripts/filter_fasta_as_orthomcl.py      \
    data/filtered_cds/all.cds                   \
    data/proteomes_filtered/goodProteins.fasta  \
    data/filtered_cds/good_cds.fasta            \
    data/filtered_cds/poor_cds.fasta

# Extract CDS from the groups file:
find data/groups_cds -delete
python scripts/extract_single_copy_orthologs.py     \
    --out_dir   data/groups_cds                     \
    --infasta   data/filtered_cds/good_cds.fasta    \
    --groups    data/groups/groups.txt              \
    --species   ${FISHES[*]}
