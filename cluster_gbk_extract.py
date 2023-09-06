import subprocess
import os
import datetime
from pathlib import Path
import shutil
from Bio import SeqIO
from string import digits

START_INDEX = 1
END_INDEX = 10

def delete_directory(dir_name):
    dir_path = Path(dir_name)
    for file in dir_path.iterdir():
        if os.path.isfile(file):
            os.remove(file)
        else:
            delete_directory(file.absolute())
    os.rmdir(dir_name)


def genome_fetch(species="\'Rhodococcus.*\'",
                 source="assembly_summary_genbank.txt",
                 output="ftp_links.txt"):
    cmd = f"grep -E {species} {source} | cut -f 20 > {output}"
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")


def fetch_source(source="ftp_links.txt",
                 output="download_protein_files.sh",
                 bracket_1="{FS=OFS=\"/\";filesuffix=\"genomic.gbff.gz\"}", # change the suffix here
                 bracket_2="{ftpdir=$0;asm=$10;file=asm\"_\"filesuffix;print \"rsync -t -v \"ftpdir,file\" ./\"}"):
    cmd = f"awk 'BEGIN{bracket_1}{bracket_2}' {source} | sed 's/https/rsync/g' > {output}"
    print(cmd)
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")


def run_hmmsearch(operation="hmmsearch",
            domains="RNA_pol_Rpb1_3.hmm", # RNA_pol_Rpb1_3.hmm --- for normal supercluster search / LAL.hmm / MftR.hmm / Ribosomal_protein_S12.hmm / SecY.hmm
            target="GCA_000009765.2_ASM976v2_protein.faa",
            output="rnap.out"):
    cmd = f"{operation} {domains} {target} > {output}"
    print(f"Running -> '{cmd}'")
    subprocess.run(cmd, shell=True)
    print("Command completed")


def rnap_search():
    print("Running rnap_search")
    with open("rnap.out", 'r') as file:
        # search = file.readlines()
        gene_line = 0
        word_region = ""
        subject = ""
        gene = ""
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
                x = 0
                gene = ""
                for word in line:
                    for item in word_region:
                        if x == item:
                            gene = gene + word
                    x = x + 1
        print(gene)
    return gene


def supercluster_region(gene, working_genome):
    print("Running supercluster_region")
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    count = 0
    supercluster = []
    found = False
    gene_hit_output = f'RPOB_hits.fasta' #change based on target gene - Ribosomal_S10_hits.fasta
    for index, fasta in enumerate(fasta_sequences):
        if gene == fasta.id:
            #print(f' who thought this was a good idea? {len(fasta.seq)}') # delete this
            if len(fasta.seq) < 4000: #150 for R_S10
                print("Match")
                rnap = count
                low_boundary = rnap - 35 # changed from 20 (35 for Rhodococcus - 20 for Strep?)
                upper_boundary = rnap + 35 # changed from 20
                found = True
                for i in range(low_boundary, upper_boundary):
                    supercluster.append(i)
            # reintroduce this next part to create a copy of the gene cluster - eventually could be a tick box whether this is activated?
            with open(gene_hit_output, 'a') as file:
                file.write(f'>{fasta.id} \n {fasta.seq} \n')
        else:
            count = count + 1
    if found == False:
        print("Failed")
        return found
    return supercluster


def supercluster_region_genes(gene, working_genome):
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    count = 0
    supercluster_genes = ''
    found = False
    gene_hit_output = f'RPOB_hits.fasta' #change based on target gene - Ribosomal_S10_hits.fasta
    for index, fasta in enumerate(fasta_sequences):
        if gene == fasta.id:
            if len(fasta.seq) < 4000: #150 for R_S10
                print("Match")
                rnap = count
                low_boundary = rnap - 35 # changed from 20 (35 for Rhodococcus - 20 for Strep?)
                upper_boundary = rnap + 35 # changed from 20
                found = True
                for i in range(low_boundary, upper_boundary):
                    supercluster_genes = supercluster_genes + fasta.id
        else:
            count = count + 1
    if found == False:
        print("Failed")
        return found
    return supercluster_genes


def gene_location(gene, working_genome):
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    count = 0
    rnap = 0
    gene_hit = ''
    for index, fasta in enumerate(fasta_sequences):
        if gene == fasta.id:
            rnap = count
            gene_hit = gene
        count = count + 1
    gene_on_contig = f'{gene_hit} {rnap}/{count}'
    return gene_on_contig


def supercluster_extraction(supercluster, working_genome):
    print("Running supercluster_extraction")
    fasta_sequences = SeqIO.parse(open(working_genome), 'fasta')
    supercluster_count = 0
    supercluster_output = open("supercluster_output.fasta", "w")
    for index, fasta in enumerate(fasta_sequences):
        for gene in supercluster:
            seq = ""
            if gene == supercluster_count:
                #print(f'this is working {fasta.id}, {fasta.description}')
                #print(fasta.seq)
                for i in range(1, len(fasta.seq)):
                    seq = seq + fasta.seq[i]
                supercluster_output.write(">" + fasta.id + fasta.description + "\n" + seq + "\n")
        supercluster_count = supercluster_count + 1
    supercluster_output.close()
    print("done")
    return(supercluster_output)







def main():
    global working_genome, supercluster, gene
    print("=" * 25)
    print("Cleaning up files")
    print("-" * 25)
    print('Starting - Supercluster Search')
    print("=" * 25)

    genome_fetch(species="\'Rhodococcus.*\'", output="ftp_links.txt") #change back to salmonella or Streptomyces or Burkholderia or Nocardiopsis or Planobispora or Mycobacterium or Rhodococcus or Gordonia or Sphaerisporangium or Actinomadura or Kocuria or Corynebacterium or Nocardioides
    fetch_source(source="ftp_links.txt", output="download_protein_files.sh")

    with open("download_protein_files.sh", 'r') as file:
        today = datetime.datetime.now()
        date_suffix = f"-{today.year}-{today.month}-{today.day}"
        hit_regions_directory = Path(f"hits{date_suffix}")
        hit_regions_directory.mkdir(exist_ok=True)
        count = START_INDEX
        for genome in range(1, END_INDEX):
            print(f'count {count} START')
            working_genome = ""
            supercluster = []
            gene = ""
            cmd = f"cat download_protein_files.sh | head -n {count} | tail -1"
            # cmd = f"head -n {count} download_protein_files.sh"
            print(cmd)
            search_genome = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
            print(f"search_genome cmd: {search_genome}")
            subprocess.run(search_genome, shell=True)
            cmd = "gunzip *.gbff.gz"
            subprocess.run(cmd, shell=True)

            genome_name = ""
            for letter in search_genome:
                genome_name = genome_name + letter
                if letter == "/":
                    genome_name = ""
                if genome_name.endswith(".gbff"):
                    working_genome = genome_name
                    print(f"working genome = {working_genome}")

            fasta_file_store = Path(f"fasta_file_store")
            fasta_file_store.mkdir(exist_ok=True)
            fasta_file_contigs = Path(f"{fasta_file_store}{os.sep}fasta_file_contig_{working_genome}")
            fasta_file_contigs.mkdir(exist_ok=True)
            iteration = 1
            filename = f"{fasta_file_contigs}{os.sep}fasta_rewrite{iteration}.fasta"

            working_genome_zip = f"{working_genome}.gz"
            print(f"Checking for {working_genome_zip}")
            if os.path.exists(working_genome_zip):
                print(f"Removing {working_genome_zip}")
                os.remove(working_genome_zip)
            else:
                print(f"Could not find {working_genome_zip}")
                print("-" * 25)
            print("still working")
            #count = count + 1

            if os.path.exists(working_genome):
                print(working_genome)
                with open(working_genome, "r") as genbank_file:
                    #content = file.read()
                    genome_type = ''
                    print(f'genome check')
                    print(genome_type)
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
                                file.write(f'{fasta_out}')
                                file.write("\n")
                        gene_name = ""
                        gene_sequence = ""
                        if line_number == 1:
                            print(line)
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
                                filename = f"{fasta_file_contigs}{os.sep}fasta_rewrite{iteration}.fasta"
                            if pos == '"':
                                counter = counter + 1
                                if counter == 2:
                                    set_name = gene_name
                            if pos != '"':
                                gene_name = gene_name + pos
                            if gene_name == '                     /protein_id=':
                                gene_name = ''
                                gene_found = True
                        if gene_found == True:
                            with open(filename, 'a') as file:
                                file.write(f'>{set_name}')
                                file.write("\n")
                            #print(f'>{set_name}')
                        line_number = line_number + 1
                    # print(f'take 2')
                    # print(genome_type)

            supercluter_copy = f"hits{date_suffix}/{working_genome}.fasta"
            for contig in range(1, iteration):
                target_contig = f"{fasta_file_contigs}{os.sep}fasta_rewrite{contig}.fasta"
                run_hmmsearch(target=target_contig)
                gene = rnap_search()
                if os.path.exists(target_contig):
                    supercluster = supercluster_region(gene, target_contig)
                    supercluster_genes = supercluster_region_genes(gene, target_contig)
                    if supercluster != False:
                        gene_position = gene_location(gene, target_contig)
                        summary_filename = "Summary.txt"
                        with open(summary_filename, 'a') as file:
                            file.write(f"{target_contig} ")
                            file.write(f"{gene_position} ")
                            # file.write("\n")
                        supercluster_output = supercluster_extraction(supercluster, target_contig)
                        #os.remove(target_contig) move to later!!!
                        shutil.copy("supercluster_output.fasta", supercluter_copy)

                        print(f'made it here!')
                        # current issue is the writing of fasta file stuff only???

                        target_gbk = working_genome_zip[:-3]
                        iterations = 0
                        seq = ''
                        save = False
                        fasta_number = ''
                        new_target_contig = ''
                        for letter in target_contig:
                            seq = seq + letter
                            if letter == '/':
                                seq = ''
                            if save == True:
                                fasta_number = letter
                            if seq == 'fasta_rewrite':
                                save = True
                        copy_target = False
                        with open(target_gbk, 'r') as target:
                            for line in target:
                                seq = ''
                                for letter in line:
                                    seq = seq + letter
                                    print(f'{seq}')
                                    if seq == "LOCUS":
                                        iterations = iterations + 1
                                    if iterations == fasta_number and seq == 'LOCUS':
                                        new_target_contig = f'{working_genome}{fasta_number}.gbk'
                                        copy_target = True
                                    if iterations != fasta_number:
                                        copy_target = False
                                if copy_target == True:
                                    with open(new_target_contig, 'a') as file:
                                        file.write(f'{seq}\n')

                        # write some code to re-extract the sequence as gbff not fasta!!!!
                        name_count = 0
                        genome_title = ''
                        name_add = True
                        for letter in working_genome:
                            if letter == '_':
                                name_count = name_count + 1
                            if name_count == 2:
                                name_add = False
                            if name_add == True:
                                genome_title = genome_title + letter
                        cluster_genbank = f'{genome_title}.genbank'
                        open(cluster_genbank, 'w')
                        record = True
                        supercluster_record = False
                        dna_start_record = False
                        dna_end_record = False
                        dna_start = ''
                        dna_end = ''
                        dna_ordered = ''
                        dna_collect = ''
                        dna_prepare = ''
                        dna = ''
                        dna_record = ''
                        with open(new_target_contig, "r") as target_contig:
                            for info in target_contig:
                                dna_write = False
                                seq = ''
                                check_gene = ''
                                should_be_true = False
                                record = True
                                for letter in info:
                                    seq = seq + letter
                                    if seq == "            ##Genome-Assembly-Data-START##":
                                        record = False
                                    if seq == '                     /protein_id=':
                                        check_gene = check_gene + letter
                                        should_be_true = True
                                    if seq == '                     /protein_id=' and record == False:
                                        dna_start_record = True
                                    if seq == '                     /protein_id=' and record == True:
                                        dna_start_record = False
                                        dna_end_record = True
                                    if seq == '     gene            ' and dna_start_record == True:
                                        dna_write = True
                                    if letter == '.' and dna_write == True and dna_start_record == True:
                                        dna_write = False
                                    if dna_write == True and dna_start_record == True:
                                        dna_start = dna_start + letter
                                    if dna_start == 'complement(' and dna_start_record == True:
                                        dna_start = ''
                                    if seq == '     gene            ' and dna_end_record == True:
                                        dna_prepare = True
                                    if letter == '.' and dna_prepare == True and dna_end_record == True:
                                        dna_write = True
                                    if letter == ' ' or ')' and dna_write == True and dna_end_record == True:
                                        dna_write = False
                                    if dna_write == True and dna_end_record == True:
                                        dna_end = dna_end + letter
                                    if dna_end == 'complement(' and dna_end_record == True:
                                        dna_end = ''
                                for supercluster_gene in supercluster_genes:
                                    if check_gene == supercluster_gene:
                                        supercluster_record = True
                                if supercluster_record != True and should_be_true == True:
                                    record = False
                                if record == True or supercluster_record == True:
                                    with open(cluster_genbank, 'a') as file:
                                        file.write(f'{seq}')
                                if seq == 'ORIGIN':
                                    dna_collect = True
                                    with open(cluster_genbank, 'a') as file:
                                        file.write(f'ORIGIN\n')
                                if dna_collect == True and letter != ' ':
                                    dna = dna + letter
                            with open(cluster_genbank, 'a') as file:
                                file.write(f'\n')
                            remove_digits = str.maketrans('', '', digits)
                            dna_only = dna.translate(remove_digits)
                        nucleotide_count = 0
                        region = False
                        for nucleotide in dna_only:
                            if nucleotide_count == dna_start:
                                region = True
                            if nucleotide_count == dna_end:
                                region = False
                            if region == True:
                                dna_record = dna_record + nucleotide
                            nucleotide_count = nucleotide_count + 1
                        nucleotide_count = 0
                        for nucleotide in dna_record:
                            if nucleotide_count == 0 or nucleotide_count % 60 == 0:
                                marker = nucleotide_count + 1
                                addition = f'{marker} {nucleotide}'
                                dna_ordered = dna_ordered + addition
                            if nucleotide_count % 10 == 0 and nucleotide_count % 60 != 0:
                                addition = f' {nucleotide}'
                                dna_ordered = dna_ordered + addition
                            if nucleotide_count % 10 != 0:
                                dna_ordered = dna_ordered + nucleotide
                            end_line_check = nucleotide_count + 1
                            if end_line_check % 60 == 0:
                                with open(cluster_genbank, 'a') as file:
                                    file.write(f'{dna_ordered}\n')
                                dna_ordered = ''




                        dir_name = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test"
                        dir_name_specified = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/" # add [fasta_file_store] ?
                        for fasta_file in fasta_file_store.iterdir():
                            fasta_dir = os.path.join(dir_name_specified, fasta_file)
                            fasta_dir_path = Path(fasta_dir)
                            print(f'fasta dir = {fasta_dir}')
                            for contigs in fasta_dir_path.iterdir():
                                if os.path.exists(contigs):
                                    os.remove(contigs)
                            os.rmdir(fasta_dir_path)

                        find_files = os.listdir(dir_name)
                        for file in find_files:
                            if file.endswith(".gbff"):
                                os.remove(os.path.join(dir_name, file))
                        # put downloaded files into new folder?
                        with open(summary_filename, 'r') as file:
                            for line in file:
                                if line.endswith("NO HITS DETECTED\n"):
                                    print(line)
                                    index = line.index(" NO HITS DETECTED")
                                    filename = line[:index]
                                    file_to_delete = f"{hit_regions_directory.absolute()}{os.sep}{filename}.gbff"
                                    if os.path.exists(file_to_delete):
                                        print(f"Deleting {file_to_delete}")
                                        os.remove(file_to_delete)
            print(f'count {count} END')
            count = count + 1
            if os.path.exists(working_genome):
                os.remove(working_genome)
            today = datetime.datetime.now()
            date_suffix = f"-{today.year}-{today.month}-{today.day}"
            output_directory = Path(f"output{date_suffix}")
            if os.path.exists(output_directory):
                for output_to_remove in output_directory.iterdir():
                    if os.path.exists(output_to_remove):
                        os.remove(output_to_remove)



if __name__ == "__main__":
    main()

# hmmscan for rna polymerase
# read through file - if match to cluster then turn on record
# when no matches - turn off record
# capture dna with command line code (store dna region in memory)