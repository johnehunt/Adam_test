import os
import subprocess
from pathlib import Path

from Bio import SeqIO


def run_multiple_sequences(sequence_filename):
    # Create a directory for all the fasta files
    fasta_file_directory = Path("fasta_files")
    fasta_file_directory.mkdir(exist_ok=True)
    # Create a directory for output
    output_directory = Path("output")
    output_directory.mkdir(exist_ok=True)
    # Load sequences from file
    fasta_sequences = SeqIO.parse(open(sequence_filename), 'fasta')
    for index, fasta in enumerate(fasta_sequences):
        # write fasta to a file
        print(f"Processing sequence {index}")
        fasta_name = f">{fasta.name}\n"
        data_to_save = str(fasta.seq)

        # Create input file
        filename = f"{fasta_file_directory.name}{os.sep}{fasta.name}.fasta"
        print(f"processing\n{fasta_name}{data_to_save}")
        with open(filename, 'w') as file:
            file.write(fasta_name)
            file.write(data_to_save)

        # Now process data in file
        output_filename = f"{output_directory}{os.sep}{fasta.name}.output.txt"
        run_hmmscan(output=output_filename,
                    protein_file=filename)


def run_hmmscan(operation="hmmscan",
                output="output.txt",
                options="--noali",
                domains="test",
                protein_file="test.fasta"):
    cmd = f"{operation} -o {output} {options} {domains} {protein_file}"

    print(f"Running -> '{cmd}'")

    subprocess.call(cmd, shell=True)

    print("Command completed")


if __name__ == "__main__":
    operation = "hmmscan"
    output = "output.txt"
    options = "--noali"
    domains = "test"
    protein_file = "test.fasta"

    # run_hmmscan(operation, output, options, domains, protein_file)
    #
    # # alternative use the defaults and just change one value
    # run_hmmscan(protein_file="test.fasta")
    #
    # run_hmmscan(output="output2.txt",
    #             protein_file="test.fasta")

    run_multiple_sequences("caelestamide.fasta")
