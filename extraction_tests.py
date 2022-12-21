from Bio import Entrez

Entrez.email = 'adam.j.p.hunt@warwick.ac.uk'

handle = Entrez.efetch(db='pubmed', id='IEZ16_RS25390', rettype='fasta')

# handle = Entrez.efetch(db='nuccore', id='34577062', rettype='fasta')
print(handle.read())

handle = Entrez.efetch(db='protein', id='41018013', rettype='fasta')

print(handle.read())

# rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt ./
# grep -E 'Streptomyces.*coelicolor' assembly_summary_genbank.txt

# grep -E 'Streptomyces.*coelicolor' assembly_summary_genbank.txt | cut -f 20 > ftp_links.txt
# head ftp_links.txt

# awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_protein_files.sh

# head download_protein_files.sh

# source download_protein_files.sh

# gunzip GCA_000203835.1_ASM20383v1_protein.faa.gz

# ls

# u2186477@Adams-MacBook-Air Adam_test % grep -E 'Streptomyces.*coelicolor' assembly_summary_genbank.txt | > ftp_links.txt
# u2186477@Adams-MacBook-Air Adam_test % head ftp_links.txt

# select for the supercluster - save 10 genes (so they can be added later) and constantly check for relevant genes
# if [gene name] = [list of gene names], supercluster=True / if x in (,):

# get all streptomyces genome links - save to file then open one at a time, find cluster, run hmmsearch

# some lack protein fasta - ignore sadly
# consider how to include poorly assembled genomes - take +/- 10 regions around each relevant gene

# head -n 1 ftp_links.txt


fasta_sequences = SeqIO.parse(open(sequence_filename), 'fasta')
for index, fasta in enumerate(fasta_sequences):
    # write fasta to a file
    fasta_name = f">{fasta.name}\n"


    # Create input file
    filename = f"{fasta_file_directory.name}{os.sep}{fasta.name}.fasta"
    print(f"\nProcessing number {index} - {fasta_name.strip()}")
    with open(filename, 'w') as file:
        file.write(fasta_name)
        file.write(data_to_save)


#hmmsearch RNA_pol_Rpb1_3.hmm GCA_008124905.1_ASM812490v1_protein.faa > rnap.out



#ncbi-blast-2.13.0+/bin/update_blastdb.pl --passive --decompress DNA-directed_RNA_polymerase_subunit_beta
#ncbi-blast-2.13.0+/bin/blastp -h
#ncbi - blast - 2.13
#.0 + / bin / update_blastdb.pl - -showall *


        with open(output_filename, 'r') as file:
            content = file.read()
            if content.find("[No hits detected that satisfy reporting thresholds]") != -1:
                # the file does not contain a match
                empty_file_list.append(output_filename)
            else:
                # Add to list of files with a hit
                hits_file_list.append(output_filename)

    # Remove empty files
    print("\nRemoving output files with no hits")
    for file in empty_file_list:
        print(f'\tRemoving {file}')
        os.remove(file)

    # List files with hits
    print('\nOutput Files with hits detected:')
    for file in hits_file_list:
        print(f"\t{file}")



