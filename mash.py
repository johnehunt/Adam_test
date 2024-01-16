import subprocess
from domain_search_test_2 import run_multiple_sequences
import os
import datetime
from pathlib import Path
import shutil
from Bio import SeqIO

# for query cluster
input_dir = Path("/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/GBK")
output_dir = Path("/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/fasta")
# for reference clusters
input_dir = Path("/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/MIBIG")
output_dir = Path("/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/MIBIG_fasta")
# silence as appropriate


#for gbk in input_dir.iterdir():
for gbk in input_dir.glob("*.gbk"):
    identifier = ''
    gene_copy = False
    gbk_str = str(gbk)
    gbk_rewrite = ''
    for letter in gbk_str:
        gbk_rewrite = gbk_rewrite + letter
        if letter == '/':
            gbk_rewrite = ''
    gbk_rewrite = gbk_rewrite[:-4]
    print(gbk_rewrite)
    gbk_rewrite = f'{gbk_rewrite}.fasta'
    output_write = os.path.join(f'{output_dir}/{gbk_rewrite}')
    with open(gbk, "r") as gbk_read:
        line_count = 1
        DNA_copy = ''
        for line in gbk_read:
            line_copy = ''
            if line_count == 2:
                identifier = line
            if gene_copy == True:
                line_copy = ''.join(filter(str.isalpha, line))
            if gene_copy == True:
                DNA_copy = f'{DNA_copy}{line_copy}\n'
            if line == "ORIGIN\n":
                gene_copy = True
            line_count = line_count + 1
        edit = ''
        for point in identifier:
            edit = edit + point
            if edit == "DEFINITION  ":
                edit = ''
        with open(output_write, "w") as file:
            file.write(f'>{edit}')
            file.write(f'{DNA_copy}')
