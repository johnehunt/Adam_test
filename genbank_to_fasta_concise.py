from Bio import SeqIO

for genbank in gbk_to_convert:
    count = SeqIO.convert("cor6_6.gb", "genbank", "cor6_6.fasta", "fasta")
    print("Converted %i records" % count)