import subprocess
import os
import datetime
from pathlib import Path

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

print(f'starting')
genome_fetch(species="\'Rhodococcus.*\'", output="ftp_links.txt")
fetch_source(source="ftp_links.txt", output="download_protein_files.sh")


with open("download_protein_files.sh", 'r') as file:
    count = 4
    working_genome = ""
    supercluster = []
    gene = ""
    today = datetime.datetime.now()
    date_suffix = f"-{today.year}-{today.month}-{today.day}"
    hit_regions_directory = Path(f"hits{date_suffix}")
    hit_regions_directory.mkdir(exist_ok=True)
    for genome in range(1, 2):
        #search_genome = ""
        print("working")
        cmd = f"cat download_protein_files.sh | head -n {count} | tail -1"
        # cmd = f"head -n {count} download_protein_files.sh"
        print(cmd)
        search_genome = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
        print("here")
        print(f"search_genome cmd: {search_genome}") # change this from a list later
        print("There")
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
        count = count + 1

        print(working_genome)
        with open(working_genome, "r") as genbank_file:
            content = file.read()
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
                #print(f'do we skip? {no_skip_needed}')
                #print(f'{line} {set_next} {skip}') # delete this
                record = False
                skip = False
                #print(skip)
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
                            #print(f'SKIPPING {line}')  # delete this
                        if skip == False:
                            #print(f"skip F") # delete this
                            #print(f'progress? {ignore_first_iter}')
                            if AA == '"' and ignore_first_iter == False:
                                #print(f'This is the error {line}')
                                set_next = False
                            if AA == '"' and ignore_first_iter == True:
                                #print(f'This is not the error {line}')
                                ignore_first_iter = False
                            if record == True and skip == False:
                                if AA != '"' and AA != ' ':
                                    translation = translation + AA
                                    translation_len = translation_len + 1
                    #print(f'scribe {translation}')
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
                            #print(f'single line error? {single_line} {single_line_test}')
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
                    #print(f'prepping') # delete this
                    fasta = ''
                    for i in range(0, translation_len):
                        fasta = fasta + translation[i]
                    scribe = False
                    fasta_remove = ''
                    for position in fasta:
                        if scribe == True or no_skip_needed == True:
                            fasta_remove = fasta_remove + position
                            #print(f'updating fasta {fasta_remove}')  # delete this
                        if position == '=':
                            scribe = True
                    fasta_out = fasta_remove.replace("\n", "").strip()
                    #print(f'output {fasta_out}') # delete this
                    with open(filename, 'a') as file:
                        file.write(f'{fasta_out}')
                        file.write("\n")
                    #print(f'{fasta_out}')
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
                    print(f'>{set_name}')
                line_number = line_number + 1
            #print(f'take 2')
            #print(genome_type)


print(f'COMPLETE')

# use iterations to inform how many scans to run
#if gene_name = ____ skip