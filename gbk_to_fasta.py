import subprocess
from domain_search_test_2 import run_multiple_sequences
import os
import datetime
from pathlib import Path
import shutil
from Bio import SeqIO

files_to_convert = Path(f'Partial_cluster_test_gbk') #gbk_to_convert
print(f'directory {files_to_convert}')
iteration = 1
fasta_files = Path(f"{files_to_convert}")
fasta_files.mkdir(exist_ok=True)

for gbks in files_to_convert.iterdir():
    filename = f'{gbks}.fasta'
    print(f'converting {gbks}')
    if gbks != "gbk_to_convert/.DS_Store":
        with open(gbks, "r") as genbank_file:
            print(f'filename {gbks}')
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
                #print(f'line number {line}')
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
                    print('five')
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
                    print(f'is this working???')
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
                        print(f'writing {gbks}')
                        file.write(f'{fasta_out}')
                        file.write("\n")
                gene_name = ""
                gene_sequence = ""
                if line_number == 1:
                    #print(line)
                    genome_type = ''
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
                        filename = f"{gbks}{os.sep}fasta_rewrite{iteration}.fasta"
                    if pos == '"':
                        counter = counter + 1
                        if counter == 2:
                            set_name = gene_name
                    if pos != '"':
                        gene_name = gene_name + pos
                    if gene_name == '                     /product=': #change this as needed
                        gene_name = ''
                        gene_found = True
                if gene_found == True:
                    with open(filename, 'a') as file:
                        file.write(f'>{set_name}')
                        file.write("\n")
                    # print(f'>{set_name}')
                line_number = line_number + 1
            # print(f'take 2')
            # print(genome_type)