


def line_finder(file, line):
    with open(f"{file}", 'r') as input:
        pos = 1
        line_limit = line + 10
        for line in input:
            if pos < line_limit:
                print(f'{pos} {line}')
            pos = pos + 1

line_finder(f'fasta_file_store/fasta_file_contig_GCA_000010105.1_ASM1010v1_genomic.gbff/fasta_rewrite1.fasta', 3270)

# 3270 MGSVIKKRRKRMSKKKHRKLLRRTRVQRRKLGKgene1776752..1777885/locus_tag=RER_16360
# why is that added on the end???
#>BAH32585.1
#MKVQPSVKKICEKCKVIRRNGRVMVICENLRHKQRQGgene2059000..2059374/gene=rpsM
#Occurs >1