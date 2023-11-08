import subprocess
#from domain_search_test_2 import run_multiple_sequences
import os
#import datetime
from pathlib import Path
import shutil
#from Bio import SeqIO


START_INDEX = 1
END_INDEX = 100

def genome_fetch(species="\'Streptomyces.*\'",
                 source="assembly_summary_genbank.txt",
                 output="ftp_links.txt"):
    cmd = f"grep -E {species} {source} | cut -f 20 > {output}"
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")

def fetch_source(source="ftp_links.txt",
                 output="download_protein_files.sh",
                 bracket_1="{FS=OFS=\"/\";filesuffix=\"genomic.gbff.gz\"}", # change the suffix here
                 bracket_2="{ftpdir=$0;asm=$10;file=asm\"_\"filesuffix;print \"rsync -t -v \"ftpdir,file\" ./\"}"):
    cmd = f"awk 'BEGIN{bracket_1}{bracket_2}' {source} | sed 's/https/rsync/g' > {output}"
    print(cmd)
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")





def main():
    global working_genome, supercluster, gene
    print("=" * 25)
    print("Cleaning up files")
    print("-" * 25)
    print('Starting - Supercluster Search')
    print("=" * 25)

    genome_fetch(species="\'Streptomyces.*\'", output="ftp_links.txt") #change back to salmonella or Streptomyces or Burkholderia or Nocardiopsis or Planobispora or Mycobacterium or Rhodococcus or Gordonia or Sphaerisporangium or Actinomadura or Kocuria or Corynebacterium or Nocardioides
    fetch_source(source="ftp_links.txt", output="download_protein_files.sh")

with open("download_protein_files.sh", 'r') as file:
    genome_directory = Path(f"HTP_antismash/genomes")
    genome_directory.mkdir(exist_ok=True)
    count = START_INDEX
    for genome in range(1, END_INDEX):
        print(f'count {count} START')
        cmd = f"cat download_protein_files.sh | head -n {count} | tail -1"
        # cmd = f"head -n {count} download_protein_files.sh"
        print(cmd)
        search_genome = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
        # print("here")
        print(f"search_genome cmd: {search_genome}")  # change this from a list later
        # print("There")
        subprocess.run(search_genome, shell=True)
        cmd = "gunzip *.gbff.gz"
        subprocess.run(cmd, shell=True)

        genome_name = ""
        for letter in search_genome:
            genome_name = genome_name + letter
            if letter == "/":
                genome_name = ""
            if genome_name.endswith(".gbff"):
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

        if os.path.exists(working_genome):
            shutil.move(working_genome, genome_directory)

        count = count + 1
