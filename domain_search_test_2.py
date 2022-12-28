

#to do
#create one folder per genome
#incorporate cblaster
#create a summary file for ease of analysis

import os
import subprocess
from pathlib import Path

from Bio import SeqIO
import datetime
from Bio import Entrez

def run_multiple_sequences(sequence_filename):
    empty_file_list = []
    hits_file_list = []

    today = datetime.datetime.now()
    date_suffix = f"-{today.year}-{today.month}-{today.day}"
    # Create a directory for all the fasta files
    fasta_file_directory = Path(f"fasta_files{date_suffix}")
    fasta_file_directory.mkdir(exist_ok=True)
    # Create a directory for output
    output_directory = Path(f"output{date_suffix}")
    output_directory.mkdir(exist_ok=True)
    # Load sequences from file
    fasta_sequences = SeqIO.parse(open(sequence_filename), 'fasta')
    for index, fasta in enumerate(fasta_sequences):
        # write fasta to a file
        fasta_name = f">{fasta.name}\n"
        data_to_save = str(fasta.seq)

        # Create input file
        filename = f"{fasta_file_directory.name}{os.sep}{fasta.name}.fasta"
        print(f"\nProcessing number {index} - {fasta_name.strip()}")
        with open(filename, 'w') as file:
            file.write(fasta_name)
            file.write(data_to_save)

        # Now process data in file
        output_filename = f"{output_directory}{os.sep}{fasta.name}.output.txt"
        run_hmmscan(output=output_filename,
                    protein_file=filename)

        # Now look to see if output file has any data in it
        # It assumes 'No hits detected' in the file indicates the file is not useful
        with open(output_filename, 'r') as file:
            content = file.read()
            if content.find("[No hits detected that satisfy reporting thresholds]") != -1:
                # the file does not contain a match
                empty_file_list.append(output_filename)
            else:
                # Add to list of files with a hit
                hits_file_list.append(output_filename)

    # create a summary output
    count = 0
    summary_filename = "Summary.txt"
    for hits in hits_file_list:
        count = count + 1
    if count != 0:
        # write genome name and successful gene hits to a summary file
        with open(summary_filename, 'a') as file:
            file.write(f"Genome hits: {count}")
            file.write("\n")
    else:
        with open(summary_filename, 'a') as file:
            file.write(f"NO HITS DETECTED")
            file.write("\n")

    # Remove empty files
    print("\nRemoving output files with no hits")
    for file in empty_file_list:
        print(f'\tRemoving {file}')
        os.remove(file)

    # List files with hits
    print('\nOutput Files with hits detected:')
    for file in hits_file_list:
        print(f"\t{file}")


def run_hmmscan(operation="hmmscan",
                output="output.txt",
                options="--noali",
                domains="test",
                protein_file="test.fasta"):
    cmd = f"{operation} -o {output} {options} {domains} {protein_file}"
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")


if __name__ == "__main__":
    # # alternative use the defaults and just change one value
    # run_hmmscan(protein_file="test.fasta")
    #
    # run_hmmscan(output="output2.txt",
    #             protein_file="test.fasta")

    run_multiple_sequences("caelestamide.fasta")


def run_cblaster(operation="cblaster search",
                query="protein_expression.fasta",
                 # taxinomy is an optional field input in the form of {-eg "txid1902[orgn]"}
                database="refseq_protein",
                hits="2",
                output="summary.csv",
                genes="100",
                range="20000",
                save="session.json",
                max_distance="50000"):
    cmd = f"{operation} -qf {query} -db {database} -mh {hits} -o {output} -ig -mic {genes} -g {range} -s {save} -md {max_distance}"
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")
    print("Extracting clusters")
    extract = f"cblaster extract_clusters session.json -o example_directory"
    subprocess.run(extract, shell=True)

# cblaster currently gives locus tag as protein_ID which prevents hmmscan - download datasets and serch with locus tag?

# cblaster search -qf protein_expression.fasta -db refseq_protein -eq "txid1902[orgn]" -mh 7 -o summary.cvs -ig -mic 100 -g 20000 -b binary.csv -s session.json -md 50000 -p plot.html
# cblaster extract_clusters session.json -o example_directory
# cblaster plot_clusters session.json -o plot.html

# either use the website links to access the fasta files rather than cblaster
# or code and get genbanks that have the locus tag like in cblaster


On Wed, Dec 28, 2022 at 3:52 PM Adam Hunt <adam.j.hunt2@gmail.com> wrote:

    # add automatic download
    # add automatic hmmsearch


    import subprocess
    from domain_search_test_2 import run_multiple_sequences
    import os
    import datetime
    from pathlib import Path
    import shutil
    import domain_search_test_2

    # run this in command line the first time
    # rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt ./

    # grep - E 'Streptomyces.*' assembly_summary_genbank.txt | cut -f 20 > ftp_links.txt
    # awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_protein_files.sh
    # rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/765/GCA_000009765.2_ASM976v2/GCA_000009765.2_ASM976v2_protein.faa.gz ./ # what source does
    # source download_protein_files.sh
    # head -n 1 download_protein_files.sh
    # gunzip *.faa.gz


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


    clean_up()

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
                domains="RNA_pol_Rpb1_3.hmm",
                target="GCA_000009765.2_ASM976v2_protein.faa",
                output="rnap.out"):
        cmd = f"{operation} {domains} {target} > {output}"
        print(f"Running -> '{cmd}'")
        subprocess.run(cmd, shell=True)
        print("Command completed")


    print('Starting')

    genome_fetch(species="\'Streptomyces.*\'", output="ftp_links.txt")
    fetch_source(source="ftp_links.txt", output="download_protein_files.sh")




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

    from Bio import SeqIO

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







    # delete directory when done
    # save file name
    # save positive result fastas
    # turn commands into functions to be called


    with open("download_protein_files.sh", 'r') as file:
        count = 1635 # 1636 dont change range change count - set to one to start from beginning
        working_genome = ""
        supercluster = []
        gene = ""
        today = datetime.datetime.now()
        date_suffix = f"-{today.year}-{today.month}-{today.day}"
        hit_regions_directory = Path(f"hits{date_suffix}")
        hit_regions_directory.mkdir(exist_ok=True)
        for genome in range(1, 1000):
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

