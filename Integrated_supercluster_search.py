

# add automatic download
# add automatic hmmsearch
# navigate based on accession number

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
END_INDEX = 1000

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



def supercluster_region(gene, working_genome):
    print("Running supercluster_region")
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    count = 0
    supercluster = []
    found = False
    gene_hit_output = f'RPOB_hits.fasta' #change based on target gene - Ribosomal_S10_hits.fasta
    for index, fasta in enumerate(fasta_sequences):
        if gene == fasta.id:
            #print(f' who thought this was a good idea? {len(fasta.seq)}') # delete this
            if len(fasta.seq) < 4000: #150 for R_S10
                print("Match")
                rnap = count
                low_boundary = rnap - 20 # changed from 20 (35 for Rhodococcus - 20 for Strep?)
                upper_boundary = rnap + 20 # changed from 20
                found = True
                for i in range(low_boundary, upper_boundary):
                    supercluster.append(i)
            # reintroduce this next part to create a copy of the gene cluster - eventually could be a tick box whether this is activated?
            with open(gene_hit_output, 'a') as file:
                file.write(f'>{fasta.id} \n {fasta.seq} \n')
        else:
            count = count + 1
    if found == False:
        print("Failed")
        return found
    return supercluster

def gene_location(gene, working_genome):
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    count = 0
    rnap = 0
    gene_hit = ''
    for index, fasta in enumerate(fasta_sequences):
        if gene == fasta.id:
            rnap = count
            gene_hit = gene
        count = count + 1
    gene_on_contig = f'{gene_hit} {rnap}/{count}'
    return gene_on_contig

def supercluster_extraction(supercluster, working_genome):
    print("Running supercluster_extraction")
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    supercluster_count = 0
    supercluster_output = open("supercluster_output.fasta", "w")
    for index, fasta in enumerate(fasta_sequences):
        for gene in supercluster:
            seq = ""
            if gene == supercluster_count:
                #print(f'this is working {fasta.id}, {fasta.description}')
                #print(fasta.seq)
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

    genome_fetch(species="\'Streptomyces.*\'", output="ftp_links.txt") #change back to salmonella or Streptomyces or Burkholderia or Nocardiopsis or Planobispora or Mycobacterium or Rhodococcus or Gordonia or Sphaerisporangium or Actinomadura or Kocuria or Corynebacterium or Nocardioides
    fetch_source(source="ftp_links.txt", output="download_protein_files.sh")

    # delete directory when done
    # save file name
    # save positive result fastas
    # turn commands into functions to be called

    with open("download_protein_files.sh", 'r') as file:
        today = datetime.datetime.now()
        date_suffix = f"-{today.year}-{today.month}-{today.day}"
        hit_regions_directory = Path(f"hits{date_suffix}")
        hit_regions_directory.mkdir(exist_ok=True)
        count = START_INDEX
        for genome in range(1, END_INDEX):
            print(f'count {count} START')
            working_genome = ""
            supercluster = []
            gene = ""
            search_genome = ""
            #print("working")
            cmd = f"cat download_protein_files.sh | head -n {count} | tail -1"
            # cmd = f"head -n {count} download_protein_files.sh"
            print(cmd)
            search_genome = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
            #print("here")
            print(f"search_genome cmd: {search_genome}")  # change this from a list later
            #print("There")
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

            fasta_file_store = Path(f"fasta_file_store")
            fasta_file_store.mkdir(exist_ok=True)
            fasta_file_contigs = Path(f"{fasta_file_store}{os.sep}fasta_file_contig_{working_genome}")
            fasta_file_contigs.mkdir(exist_ok=True)
            iteration = 1
            filename = f"{fasta_file_contigs}{os.sep}fasta_rewrite{iteration}.fasta"

            working_genome_zip = f"{working_genome}.gz"
            print(f"Checking for {working_genome_zip}")
            if os.path.exists(working_genome_zip):
                print(f"Removing {working_genome_zip}")
                os.remove(working_genome_zip)
            else:
                print(f"Could not find {working_genome_zip}")
                print("-" * 25)
            print("still working")
            #count = count + 1

            if os.path.exists(working_genome):
                print(working_genome)
                with open(working_genome, "r") as genbank_file:
                    #content = file.read()
                    genome_type = ''
                    print(f'genome check')
                    print(genome_type)
                    line_number = 0
                    record = False
                    gene_found = False
                    set_next = False
                    translation = ""
                    single_line_test = False
                    single_line = False
                    no_skip_needed = True
                    skip = False
                    for line in genbank_file:
                        ignore_first_iter = False
                        if skip == True:
                            ignore_first_iter = True
                        record = False
                        skip = False
                        if set_next == False:
                            translation = ""
                            translation_len = 0
                        if set_next == True:
                            record = True
                            for AA in line:
                                if AA == '_':
                                    record = False
                                    skip = True
                                    no_skip_needed = False
                                if skip == False:
                                    if AA == '"' and ignore_first_iter == False:
                                        set_next = False
                                    if AA == '"' and ignore_first_iter == True:
                                        ignore_first_iter = False
                                    if record == True and skip == False:
                                        if AA != '"' and AA != ' ':
                                            translation = translation + AA
                                            translation_len = translation_len + 1
                        if gene_found == True:
                            single_line_test = False
                            single_line = False
                            for AA in line:
                                if AA == '_':
                                    record = False
                                    skip = True
                                    gene_found = False
                                    set_next = True
                                if skip == False:
                                    if AA == '"' and single_line_test is True:
                                        single_line = True
                                        set_next = False
                                    if AA == '"' and single_line is False:
                                        record = True
                                        gene_found = False
                                        set_next = True
                                        single_line_test = True
                                    if set_next == True:
                                        if AA != '"' and AA != ' ':
                                            translation = translation + AA
                                            translation_len = translation_len + 1
                        if set_next == False and translation != "":
                            fasta = ''
                            for i in range(0, translation_len):
                                fasta = fasta + translation[i]
                            scribe = False
                            fasta_remove = ''
                            for position in fasta:
                                if scribe == True or no_skip_needed == True:
                                    fasta_remove = fasta_remove + position
                                if position == '=':
                                    scribe = True
                            fasta_out = fasta_remove.replace("\n", "").strip()
                            with open(filename, 'a') as file:
                                file.write(f'{fasta_out}')
                                file.write("\n")
                        gene_name = ""
                        gene_sequence = ""
                        if line_number == 1:
                            print(line)
                            for letter in line:
                                if letter == ',':
                                    record = True
                                if record == True and letter != ',':
                                    genome_type = genome_type + letter
                        counter = 0
                        contig_check = ''
                        for pos in line:
                            contig_check = contig_check + pos
                            if contig_check == "ORIGIN":
                                iteration = iteration + 1
                                filename = f"{fasta_file_contigs}{os.sep}fasta_rewrite{iteration}.fasta"
                            if pos == '"':
                                counter = counter + 1
                                if counter == 2:
                                    set_name = gene_name
                            if pos != '"':
                                gene_name = gene_name + pos
                            if gene_name == '                     /protein_id=':
                                gene_name = ''
                                gene_found = True
                        if gene_found == True:
                            with open(filename, 'a') as file:
                                file.write(f'>{set_name}')
                                file.write("\n")
                            #print(f'>{set_name}')
                        line_number = line_number + 1
                    # print(f'take 2')
                    # print(genome_type)

            supercluter_copy = f"hits{date_suffix}/{working_genome}.fasta"
            for contig in range(1, iteration):
                target_contig = f"{fasta_file_contigs}{os.sep}fasta_rewrite{contig}.fasta"
                run_hmmsearch(target=target_contig)
                gene = rnap_search()
                if os.path.exists(target_contig):
                    supercluster = supercluster_region(gene, target_contig)
                    if supercluster != False:
                        gene_position = gene_location(gene, target_contig)
                        summary_filename = "Summary.txt"
                        with open(summary_filename, 'a') as file:
                            file.write(f"{target_contig} ")
                            file.write(f"{gene_position} ")
                            # file.write("\n")
                        supercluster_output = supercluster_extraction(supercluster, target_contig)
                        os.remove(target_contig)
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
                        dir_name = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test"
                        dir_name_specified = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/" # add [fasta_file_store] ?
                        for fasta_file in fasta_file_store.iterdir():
                            fasta_dir = os.path.join(dir_name_specified, fasta_file)
                            fasta_dir_path = Path(fasta_dir)
                            print(f'fasta dir = {fasta_dir}')
                            for contigs in fasta_dir_path.iterdir():
                                if os.path.exists(contigs):
                                    os.remove(contigs)
                            os.rmdir(fasta_dir_path)

                        find_files = os.listdir(dir_name)
                        for file in find_files:
                            if file.endswith(".gbff"):
                                os.remove(os.path.join(dir_name, file))
                        # put downloaded files into new folder?
                        with open(summary_filename, 'r') as file:
                            for line in file:
                                if line.endswith("NO HITS DETECTED\n"):
                                    print(line)
                                    index = line.index(" NO HITS DETECTED")
                                    filename = line[:index]
                                    file_to_delete = f"{hit_regions_directory.absolute()}{os.sep}{filename}.gbff"
                                    if os.path.exists(file_to_delete):
                                        print(f"Deleting {file_to_delete}")
                                        os.remove(file_to_delete)
            print(f'count {count} END')
            count = count + 1
            if os.path.exists(working_genome):
                os.remove(working_genome)
            today = datetime.datetime.now()
            date_suffix = f"-{today.year}-{today.month}-{today.day}"
            output_directory = Path(f"output{date_suffix}")
            if os.path.exists(output_directory):
                for output_to_remove in output_directory.iterdir():
                    if os.path.exists(output_to_remove):
                        os.remove(output_to_remove)


if __name__ == "__main__":
    main()
    clear_out_empty_files()


# add code to ignore files without protein data
# delete unncessary files - stop output writing, delete fasta files, make sure gbff files are deleted