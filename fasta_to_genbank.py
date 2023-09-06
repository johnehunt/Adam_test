from Bio import SeqIO

input_handle = open("test.fasta", "r")
output_handle = open("test.gb", "w")

sequences = list(SeqIO.parse(input_handle, "fasta"))

# assign molecule type
for seq in sequences:
  seq.annotations['molecule_type'] = 'DNA'

count = SeqIO.write(sequences, output_handle, "genbank")

output_handle.close()
input_handle.close()
print("Converted {} records".format(count))