import subprocess
from Bio import SeqIO
import os

def rnap_search():
    print("Running rnap_search")
    with open("rnap.out", 'r') as file:
        # search = file.readlines()
        result = 0
        gene_line = 0
        word_region = ""
        subject = ""
        region = []
        line_count = 0
        for line in file:
            line_count = line_count + 1
            # print(line)
            target = "    E-value  score  bias    E-value  score  bias    exp  N  Sequence   Description"
            if line.find(target) != -1:
                gene_line = line_count + 2
                pos = 0
                for word in line:
                    subject = subject + word
                    region.append(pos) # region = [region + pos]
                    if word == " ":
                        subject = ""
                        region = []
                    if subject == "Sequence":
                        additional = pos + 2
                        region.append(additional)
                        word_region = region
                    pos = pos + 1
            if line_count == gene_line:
                output = ""
                result = result + 1
                next_line = gene_line + 1
                x = 0
                gene = ""
                for word in line:
                    output = output + word
                    for item in word_region:
                        if x == item:
                            gene = gene + word
                    x = x + 1
                with open(result_data, 'a') as file:
                    file.write(f'{output}')
                    file.write("\n")
            result_check = ""
            while result_check == " ":
                x = 0
                gene = ""
                for word in line:
                    for item in word_region:
                        if x == item:
                            gene = gene + word
                            if gene.length == 1:
                                result_check = gene
                                if gene != " ":
                                    result = result + 1
                                    next_line = next_line + 1
    return (result)
#check for any data under sequence until there is a space for number of hits?

def genome_fetch(species="\'Streptomyces.*\'",
                 source="assembly_summary_genbank.txt",
                 output="ftp_links.txt"):
    cmd = f"grep -E {species} {source} | cut -f 20 > {output}"
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")

def fetch_source(source="ftp_links.txt",
                 output="download_protein_files.sh",
                 bracket_1="{FS=OFS=\"/\";filesuffix=\"protein.faa.gz\"}", # change the suffix here
                 bracket_2="{ftpdir=$0;asm=$10;file=asm\"_\"filesuffix;print \"rsync -t -v \"ftpdir,file\" ./\"}"):
    cmd = f"awk 'BEGIN{bracket_1}{bracket_2}' {source} | sed 's/https/rsync/g' > {output}"
    print(cmd)
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")

def run_hmmsearch(operation="hmmsearch",
            domains="SPASM.hmm", #RNA_pol_Rpb1_3.hmm or tryptophan_halogenase.hmm or LAL.hmm or MftR.hmm or SPASM.hmm
            target="GCA_000009765.2_ASM976v2_protein.faa",
            output="rnap.out"):
    cmd = f"{operation} {domains} {target} > {output}"
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")

def supercluster_region(gene):
    print("Running gene check")
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    found = False
    for index, fasta in enumerate(fasta_sequences):
        if gene == fasta.id:
            print("Match")
            found = True
            rpoB_in_genome = rpoB_in_genome + 1
    if found == False:
        print("Failed")
        return found
    return supercluster

genome_fetch(species="\'Streptomyces.*\'", output="ftp_links.txt") # Streptomyces or Rhodococcus
fetch_source(source="ftp_links.txt", output="download_protein_files.sh")

result_data = "Gene_count.txt"
with open(result_data, 'a') as file:
    file.write(f'    E-value  score  bias    E-value  score  bias    exp  N  Sequence   Description')
    file.write("\n")
working_genome = ""
supercluster = []
rpoB = []
count = 0

def main():
    with open("download_protein_files.sh", 'r') as file:
        count = 0
        working_genome = ""
        rpoB = []
        for genome in range(1, 1500):
            search_genome = ""
            print("working")
            cmd = f"cat download_protein_files.sh | head -n {count} | tail -1"
            print(cmd)
            search_genome = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
            print("here")
            print(f"search_genome cmd: {search_genome}")
            print("There")
            subprocess.run(search_genome, shell=True)
            cmd = "gunzip *.faa.gz"
            subprocess.run(cmd, shell=True)
            genome_name = ""
            for letter in search_genome:
                genome_name = genome_name + letter
                if letter == "/":
                    genome_name = ""
                if genome_name.endswith(".faa"):
                    working_genome = genome_name
            print(f"working genome = {working_genome}")
            working_genome_zip = f"{working_genome}.gz"
            print(f"Checking for {working_genome_zip}")
            if os.path.exists(working_genome_zip):
                print(f"Removing {working_genome_zip}")
                os.remove(working_genome_zip)
            else:
                print(f"Could not find {working_genome_zip}")
            print("-" * 25)
            print("still working")
            run_hmmsearch(target=working_genome)
            result = rnap_search()
            add = False
            if os.path.exists(working_genome):
                print(f'printing working geome {working_genome}')
                os.remove(working_genome)
                add = True
            if add == True:
                rpoB.append(result)
                print(result)
            count = count + 1
        Zero = 0
        One = 0
        Two = 0
        Three = 0
        Four_or_more = 0
        for hits in rpoB:
            greater = True
            if hits == 0:
                Zero = Zero + 1
                greater = False
            if hits == 1:
                One = One + 1
                greater = False
            if hits == 2:
                Two = Two + 1
                greater = False
            if hits == 3:
                Three = Three + 1
                greater = False
            if greater is True:
                Four_or_more = Four_or_more + 1
        print(
            f'Zero rpoB: {Zero} /n One rpoB: {One} /n Two rpoB: {Two} /n Three rpoB: {Three} /n Four or more rpoB: {Four_or_more}')


if __name__ == "__main__":
    main()




# is rnap.out needed?







# write output into a file (with "    E-value  score  bias    E-value  score  bias    exp  N  Sequence   Description" as a header)


