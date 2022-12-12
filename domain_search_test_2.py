import os
import subprocess
from pathlib import Path

from Bio import SeqIO
import datetime


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

    # Remove empty files
    print("Removing output files with no hits")
    for file in empty_file_list:
        print(f'\tRemoving {file}')
        os.remove(file)

    # List files with hits
    print('Output Files with hits detected:')
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
