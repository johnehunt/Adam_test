

# add automatic download
# add automatic hmmsearch


import subprocess
from domain_search_test_2 import run_multiple_sequences
import os
import datetime
from pathlib import Path
import shutil
from Bio import SeqIO

# run this in command line the first time
# rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt ./

# grep -E 'Streptomyces.*' assembly_summary_genbank.txt | cut -f 20 > ftp_links.txt
# grep -E 'Streptomyces.*' assembly_summary_genbank.txt | cut -f 20 > ftp_links.txt
# awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_protein_files.sh
# rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/765/GCA_000009765.2_ASM976v2/GCA_000009765.2_ASM976v2_protein.faa.gz ./ # what source does
# source download_protein_files.sh
# head -n 1 download_protein_files.sh
# gunzip *.faa.gz


print("Defining globals")
working_genome = ""
supercluster = []
gene = ""

START_INDEX = 1
END_INDEX = 1500

def delete_directory(dir_name):
    dir_path = Path(dir_name)
    for file in dir_path.iterdir():
        if os.path.isfile(file):
            os.remove(file)
        else:
            delete_directory(file.absolute())
    os.rmdir(dir_name)

def clean_up():
    """Look to see if files and directories need to be removed
    to clean up environment before starting to run the program"""
    file_to_delete = "Summary.txt"
    if os.path.exists(file_to_delete):
        print(f"Deleting {file_to_delete}")
        os.remove(file_to_delete)
    directories = [f for f in os.listdir('.') if os.path.isdir(f)]
    directories = list(filter(lambda f: f.startswith("hits") or f.startswith("output"), directories))
    print(directories)
    for directory in directories:
        print(f'Deleting {directory}')
        delete_directory(directory)

def genome_fetch(species="\'Streptomyces.*\'",
                 source="assembly_summary_genbank.txt",
                 output="ftp_links.txt"):
    cmd = f"grep -E {species} {source} | cut -f 20 > {output}"
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")

def fetch_source(source="ftp_links.txt",
                 output="download_protein_files.sh",
                 bracket_1="{FS=OFS=\"/\";filesuffix=\"protein.faa.gz\"}", # change the suffix here
                 bracket_2="{ftpdir=$0;asm=$10;file=asm\"_\"filesuffix;print \"rsync -t -v \"ftpdir,file\" ./\"}"):
    cmd = f"awk 'BEGIN{bracket_1}{bracket_2}' {source} | sed 's/https/rsync/g' > {output}"
    print(cmd)
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")

def run_hmmsearch(operation="hmmsearch",
            domains="RNA_pol_Rpb1_3.hmm", # RNA_pol_Rpb1_3.hmm --- for normal supercluster search / LAL.hmm / MftR.hmm / Ribosomal_protein_S12.hmm / SecY.hmm
            target="GCA_000009765.2_ASM976v2_protein.faa",
            output="rnap.out"):
    cmd = f"{operation} {domains} {target} > {output}"
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")

def rnap_search():
    print("Running rnap_search")
    with open("rnap.out", 'r') as file:
        # search = file.readlines()
        gene_line = 0
        word_region = ""
        subject = ""
        gene = ""
        region = []
        line_count = 0
        for line in file:
            line_count = line_count + 1
            # print(line)
            target = "    E-value  score  bias    E-value  score  bias    exp  N  Sequence   Description"
            if line.find(target) != -1:
                gene_line = line_count + 2
                pos = 0
                for word in line:
                    subject = subject + word
                    region.append(pos) # region = [region + pos]
                    if word == " ":
                        subject = ""
                        region = []
                    if subject == "Sequence":
                        additional = pos + 2
                        region.append(additional)
                        word_region = region
                    pos = pos + 1
            if line_count == gene_line:
                x = 0
                gene = ""
                for word in line:
                    for item in word_region:
                        if x == item:
                            gene = gene + word
                    x = x + 1
        print(gene)
    return gene



def supercluster_region(gene):
    print("Running supercluster_region")
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    count = 0
    supercluster = []
    found = False
    for index, fasta in enumerate(fasta_sequences):
        if gene == fasta.id:
            print("Match")
            rnap = count
            low_boundary = rnap - 20
            upper_boundary = rnap + 20
            found = True
            for i in range(low_boundary, upper_boundary):
                supercluster.append(i)
        else:
            count = count + 1
    if found == False:
        print("Failed")
        return found
    return supercluster

def supercluster_extraction(supercluster):
    print("Running supercluster_extraction")
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    supercluster_count = 0
    supercluster_output = open("supercluster_output.fasta", "w")
    for index, fasta in enumerate(fasta_sequences):
        for gene in supercluster:
            seq = ""
            if gene == supercluster_count:
                # print(fasta.id, fasta.description)
                # print(fasta.seq)
                for i in range(1, len(fasta.seq)):
                    seq = seq + fasta.seq[i]
                supercluster_output.write(">" + fasta.id + fasta.description + "\n" + seq + "\n")
        supercluster_count = supercluster_count + 1
    supercluster_output.close()
    print("done")
    return(supercluster_output)


def clear_out_empty_files():
    print("Checking for empty files!")
    dir_path = Path(".")
    for file in dir_path.iterdir():
        if os.path.isfile(file):
            if os.path.getsize(file) == 0:
                print(f"Removig empty file {file}")
                os.remove(file)



def main():
    global working_genome, supercluster, gene
    print("=" * 25)
    print("Cleaning up files")
    clean_up()
    print("-" * 25)
    print('Starting - Supercluster Search')
    print("=" * 25)

    genome_fetch(species="\'Pantoea.*\'", output="ftp_links.txt") #change back to salmonella or Streptomyces or Burkholderia or Nocardiopsis or Planobispora or Mycobacterium or Rhodococcus or Gordonia or Sphaerisporangium or Actinomadura or Kocuria or Corynebacterium or Nocardioides
    fetch_source(source="ftp_links.txt", output="download_protein_files.sh")

    # delete directory when done
    # save file name
    # save positive result fastas
    # turn commands into functions to be called

    with open("download_protein_files.sh", 'r') as file:
        count = START_INDEX # 1636 dont change range change count - set to one to start from beginning
        working_genome = ""
        supercluster = []
        gene = ""
        today = datetime.datetime.now()
        date_suffix = f"-{today.year}-{today.month}-{today.day}"
        hit_regions_directory = Path(f"hits{date_suffix}")
        hit_regions_directory.mkdir(exist_ok=True)
        for genome in range(1, END_INDEX):
            search_genome = ""
            print("working")
            cmd = f"cat download_protein_files.sh | head -n {count} | tail -1"
            # cmd = f"head -n {count} download_protein_files.sh"
            print(cmd)
            search_genome = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
            print("here")
            print(f"search_genome cmd: {search_genome}") # change this from a list later
            print("There")
            subprocess.run(search_genome, shell=True)
            cmd = "gunzip *.faa.gz"
            subprocess.run(cmd, shell=True)

            genome_name = ""
            for letter in search_genome:
                genome_name = genome_name + letter
                if letter == "/":
                    genome_name = ""
                if genome_name.endswith(".faa"):
                    working_genome = genome_name
            print(f"working genome = {working_genome}")

            working_genome_zip = f"{working_genome}.gz"
            print(f"Checking for {working_genome_zip}")
            if os.path.exists(working_genome_zip):
                print(f"Removing {working_genome_zip}")
                os.remove(working_genome_zip)
            else:
                print(f"Could not find {working_genome_zip}")
            print("-" * 25)
            print("still working")

            supercluter_copy = f"hits{date_suffix}/{working_genome}.fasta"
            run_hmmsearch(target=working_genome)
            gene = rnap_search()
            if os.path.exists(working_genome):
                supercluster = supercluster_region(gene)
                if supercluster != False:
                    summary_filename = "Summary.txt"
                    with open(summary_filename, 'a') as file:
                        file.write(f"{working_genome} ")
                        # file.write("\n")
                    supercluster_output = supercluster_extraction(supercluster)
                    os.remove(working_genome)
                    shutil.copy("supercluster_output.fasta", supercluter_copy)
                    run_multiple_sequences("supercluster_output.fasta")
                    fasta_files = Path(f"fasta_files{date_suffix}")
                    for fasta_file in fasta_files.iterdir():
                        print(f'Attempting to Delete - {fasta_file}')
                        if os.path.exists(fasta_file):
                            print(f'Deleting {fasta_file}')
                            os.remove(fasta_file)
                        else:
                            print(f"FILE seems to have disappeared {fasta_file}")
                    os.rmdir(fasta_files)
                    with open(summary_filename, 'r') as file:
                        for line in file:
                            if line.endswith("NO HITS DETECTED\n"):
                                print(line)
                                index = line.index(" NO HITS DETECTED")
                                filename = line[:index]
                                file_to_delete = f"{hit_regions_directory.absolute()}{os.sep}{filename}.fasta"
                                if os.path.exists(file_to_delete):
                                    print(f"Deleting {file_to_delete}")
                                    os.remove(file_to_delete)

            count = count + 1

if __name__ == "__main__":
    main()
    clear_out_empty_files()
