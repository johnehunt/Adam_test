import subprocess
import os
from pathlib import Path
#import shutil
#from Bio import SeqIO
import antismash

# run this in command line the first time
# rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt ./

print("Defining globals")
working_genome = ""
supercluster = []
gene = ""

START_INDEX = 1
END_INDEX = 3

def delete_directory(dir_name):
    dir_path = Path(dir_name)
    for file in dir_path.iterdir():
        if os.path.isfile(file):
            os.remove(file)
        else:
            delete_directory(file.absolute())
    os.rmdir(dir_name)

def clean_up(file_to_delete):
    if os.path.exists(file_to_delete):
        print(f"Deleting {file_to_delete}")
        os.remove(file_to_delete)
    directories = [f for f in os.listdir('.') if os.path.isdir(f)]
    directories = list(filter(lambda f: f.startswith("hits") or f.startswith("output") or f.startswith("fasta_file_store"), directories))
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
                 bracket_1="{FS=OFS=\"/\";filesuffix=\"genomic.gbff.gz\"}", # change the suffix here
                 bracket_2="{ftpdir=$0;asm=$10;file=asm\"_\"filesuffix;print \"rsync -t -v \"ftpdir,file\" ./\"}"):
    cmd = f"awk 'BEGIN{bracket_1}{bracket_2}' {source} | sed 's/https/rsync/g' > {output}"
    print(cmd)
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")


def main():
    global working_genome, supercluster, gene # change?
    print("=" * 25)
    print("Cleaning up files")
    #clean_up() # specify
    print("-" * 25)
    print('Starting - Supercluster Search')
    print("=" * 25)

    genome_fetch(species="\'Streptomyces.*\'", output="ftp_links.txt") #change back to salmonella or Streptomyces or Burkholderia or Nocardiopsis or Planobispora or Mycobacterium or Rhodococcus or Gordonia or Sphaerisporangium or Actinomadura or Kocuria or Corynebacterium or Nocardioides
    fetch_source(source="ftp_links.txt", output="download_protein_files.sh")

    #source = /opt/anaconda3/envs/antismash
    #cmd = "conda init antismash"
    #subprocess.run(cmd, shell=True)

    # delete directory when done
    # save file name
    # save positive result fastas
    # turn commands into functions to be called

    with open("download_protein_files.sh", 'r') as file:
        hit_regions_directory = Path(f"hits_antismash_output")
        hit_regions_directory.mkdir(exist_ok=True)
        count = START_INDEX
        for genome in range(1, END_INDEX):
            print(f'count {count} START')
            working_genome = ""
            cmd = f"cat download_protein_files.sh | head -n {count} | tail -1"
            # cmd = f"head -n {count} download_protein_files.sh"
            print(cmd)
            search_genome = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
            print(f"search_genome cmd: {search_genome}")  # change this from a list later
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
                print(working_genome)
                #cmd = f"antismash {working_genome}"
                #subprocess.run(cmd, shell = True)
                antismash(working_genome)


if __name__ == "__main__":
    main()
