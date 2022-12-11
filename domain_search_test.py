from contextlib import redirect_stdout
from Bio import SeqIO
from pyhmmer import hmmsearch
from pyhmmer import hmmpress
from pyhmmer import hmmer
from pyhmmer import phmmer
import subprocess
import pyhmmer
from pyhmmer.hmmer import hmmscan


fasta_sequences = SeqIO.parse(open("caelestamide.fasta"),'fasta')
with open('output.txt', 'w') as out_file:
    for fasta in fasta_sequences:
        #x = '1'
        #with redirect_stdout('output.txt'):
        #out_file.write(x)
        #name, sequence = fasta.id, str(fasta.seq)
        #new_sequence = hmmscan(sequence, test) #run hmmscan here
        #write_fasta('query.fasta')

hmmbuild thioesterase.hmm thioesterase.sto
hmmbuild p450.hmm p450.sto
hmmbuild condensation.hmm condensation.sto
cat condensation.hmm p450.hmm thioesterase.hmm > test
hmmpress test
hmmscan test test.fasta
hmmscan -o 'output.txt' --noali  test test.fasta




