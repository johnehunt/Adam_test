

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
