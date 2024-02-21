import csv
import os
from os import path
from pandas import DataFrame
from Bio import SeqIO
import pandas
import subprocess
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from pathlib import Path
#import category_encoders as ce

# fetch ref cluster


# separate genes in gbk
def gbk_to_fasta(gbk, csv_out): # can add fasta for gbk to fasta
    protein_id = []
    protein = []
    with open(gbk, 'r') as file:
        for records in SeqIO.parse(file, "genbank"):
            for seq_feature in records.features:
                if seq_feature.type == "CDS":
                    protein_id.append({seq_feature.qualifiers["protein_id"][0]})
                    protein.append({seq_feature.qualifiers["translation"][0]})
                    #gene = (f'>{seq_feature.qualifiers["protein_id"][0]}\n{seq_feature.qualifiers["translation"][0]}\n')
                    #print(gene)
                    #with open(fasta, 'a') as file:
                    #    file.write(gene)
    df = DataFrame({'protein_id': protein_id, 'Translation': protein})
    df.to_csv(csv_out, index=False)

# hmmscan gene
#/Users/u2186477/Documents/PhD/Year\ 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/hmmer-3.4/src/hmmpress
def run_hmmscan(operation="hmmscan",
                output="output.txt",
                output_format="--tblout",
                options="--noali",
                domains="Pfam-A.hmm",
                protein_file="protein"):
    #cmd = f"{operation} -o {output} {options} {domains} {protein_file}"
    cmd = f"{operation} {output_format} {output}  {options} {domains} {protein_file}"
    subprocess.run(cmd, shell=True)
    cmd = f"{operation} --domtblout domain_out.txt  {options} {domains} {protein_file}"
    subprocess.run(cmd, shell=True)


# process HMM output
def sort_hmmscan(output, x, domain_output, BGC, domain_csv_out):
    domains = []
    with open(output, 'r') as file:
        for line in file:
            line = line.split()
            if not line[0].startswith('#'):
                gene = line[2]
                print(f'Domain {line[0]} ({line[1]}) occurs {line[17]} time/s with a confidence of {line[4]}')
                for i in range(0, int(line[17])):
                    if float(line[4]) < float(0.0000000001): # set inclusion threshold here
                        domains.append(line[1])
    domain_ref = []
    domain_pos = []
    with open(domain_output, 'r') as file:
        for line in file:
            line = line.split()
            if not line[0].startswith('#'):
                domain_ref.append(line[1])
                domain_pos.append(int(line[19]))
    if len(domain_ref) > 0:
        #domain_pos, domain_ref = (list(t) for t in zip(*sorted(zip(domain_pos, domain_ref))))
        domain_arranged = [x for _, x in sorted(zip(domain_pos, domain_ref))]
        domain_ordered = []
        for i in range(0, len(domain_arranged)):
            match = False
            for hit in domains:
                if domain_arranged[i] == hit:
                    match = True
            if match == True:
                domain_ordered.append(domain_arranged[i])
        y = x + 1
        count = 1
        if os.path.exists(domain_csv_out):
            df_gbk = pandas.read_csv(domain_csv_out)
            repeats = []
            for domain in domain_ordered:
                skip = False
                if df_gbk.isin([domain]).any().any():
                    skip = True
                if domain in repeats:
                    skip = True
                if skip == False:
                    repeats.append(domain)
                    with open('domain_record.txt', 'a') as file:
                        file.write(f'{domain} ')
        else:
            for domain in domain_ordered:
                with open('domain_record.txt', 'a') as file:
                    file.write(f'{domain} ')

        if os.path.exists(domain_csv_out):
            for hit in domain_ordered:
                gene_info = f'{BGC},{gene},{y},{count},{hit}\n'
                with open(domain_csv_out, 'a') as df_gbk:
                    df_gbk.write(gene_info)
                count += 1
        else:
            hit = domain_ordered[0]
            df_gbk = DataFrame({'BGC': BGC, 'Gene': gene, 'BGC_position': y, 'Domain_position': count, 'domain_id': hit}, index=[0])
            df_gbk.to_csv(domain_csv_out, index=False)
            if len(domains) > 1:
                skip = 0
                for hit in domain_ordered:
                    if skip != 0:
                        gene_info = f'{BGC},{gene},{y},{count},{hit}\n'
                        with open(domain_csv_out, 'a') as df_gbk:
                            df_gbk.write(gene_info)
                    skip += 1
                    count += 1


# find potential matches in ref cluster
#def ref_domain_compare(df, ref_cluster, query_cluster_domains):
    #ref_cluster_domains = df[df[cluster_id].isin([ref_cluster])]


# sequence similarity comparison for best match - set threshold
#def ref_seq_compare(ref_csv, ref_protein_id, query_csv, query_protein_id):

    # flip sequence and do both comparisons to find right orientation


# subcluster assessment




# split cluster






#run the program
def process(in_gbk, BGC, domain_csv):
    domain_csv_out = domain_csv
    #in_gbk = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/BGC0000455.gbk"
    csv_out = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/van_genes.csv"
    #out_fasta = 'Van.fasta'
    gbk_to_fasta(in_gbk, csv_out) # out_fasta
    df = pandas.read_csv(csv_out)
    pandas.set_option('display.max_colwidth', None)
    x = 0
    #BGC = 'BGC0000455'
    for gene in range(0, len(df)):
        gene_to_extract = df.iloc[x, 0]
        gene_extracted = df[df["protein_id"].isin([gene_to_extract])]
        #print(gene_extracted)
        gene_seq = gene_extracted['Translation'].to_string(index=False)
        trans_table = gene_seq.maketrans("", "", "'{}")
        gene_seq = gene_seq.translate(trans_table)
        gene_name = gene_to_extract.translate(trans_table)
        if os.path.exists('seq_storage.csv'):
            gene_info = f'{BGC},{gene_name},{gene_seq}\n'
            with open('seq_storage.csv', 'a') as df_store:
                df_store.write(gene_info)
        if not os.path.exists('seq_storage.csv'):
            df_store = DataFrame({'Cluster': BGC, 'protein_id': gene_name, 'Translation': gene_seq}, index=[0])
            df_store.to_csv('seq_storage.csv', index=False)

        #gene_seq.translate(string.maketrans("", "", ), "''{}")
        #print(gene_seq)
        with open('protein_file.fasta', 'w') as file:
            file.write(f'>{gene_name}\n{gene_seq}')
        run_hmmscan(protein_file='protein_file.fasta')
        sort_hmmscan('output.txt', x, 'domain_out.txt', BGC, domain_csv_out)
        x += 1

def oh_ec(domains_csv, df_name):
    domain_db = pandas.read_csv(domains_csv)
    ct = ColumnTransformer([('oh_enc', OneHotEncoder(sparse_output=False), [4]), ], remainder='passthrough')
    domain_db_transformed = ct.fit_transform(domain_db)
    print(domain_db_transformed)
    count_ON = True
    Pfam = []
    for row in domain_db_transformed:
        count = 0
        gene_info = []
        domain_info = []
        pfam_record = True
        for column in row:
            count += 1
            Pfam.append(count)
            if column != 0.0 and column != 1.0:
                pfam_record = False
            if pfam_record == False:
                gene_info.append(column)
            if pfam_record == True:
                domain_info.append(column)
                count_set = count
        if count_ON == True:
            df_transformed = DataFrame({'BGC': gene_info[0], 'Gene': gene_info[1], 'BGC_position': gene_info[2], 'Domain_position': gene_info[3],}, index=[0])
            df_transformed.to_csv(df_name, index=False)
            transformed_csv = pandas.read_csv(df_name)
            with open('domain_record.txt', 'r') as file:
                domain_list = file.read()
                domain_list = domain_list.split(" ")
            for i in range(0, count_set):
                #Pfam_in = f'Pfam{i+1}'
                Pfam_in = domain_list[i]
                transformed_csv[Pfam_in] = transformed_csv['Domain_position']
                #transformed_csv.to_csv('domains_encoded.csv', index=False)
                if (i+1) != count_set:
                    #domain_add = f'{str(domain_info[i])},'
                    domain_add = domain_info[i]
                if (i+1) == count_set:
                    #domain_add = f'{str(domain_info[i])}\n'
                    domain_add = domain_info[i]
                #print(domain_add)
                transformed_csv[Pfam_in] = [domain_add]
                transformed_csv = transformed_csv.replace(Pfam_in, domain_add, regex=True)
                #df_transformed.at[0, Pfam_in] = domain_add
                #df_transformed = df_transformed.assign(Pfam_in = domain_add)
                #df_transformed.replace(to_replace='1', value=domain_add, regex=False)
                #with open('domains_encoded.csv', 'a') as df_transformed:
                #    df_transformed.write(domain_add)
                #df_transformed.assign(Pfam_in=domain_add)
            transformed_csv.to_csv(df_name, index=False)
        if count_ON == False:
            input_info = gene_info + domain_info
            with open(df_name, 'a') as df_transformed:
                row_append = csv.writer(df_transformed)
                row_append.writerow(input_info)
        count_ON = False
    #df_transformed = DataFrame.to_csv(domain_db_transformed)
    #print(df_transformed)

def seq_compare(hit_list, reference): # return list with similarity and with best hits



def cluster_compare(query):
    df_ref = pandas.read_csv('domains_ref_encoded.csv')
    df_test = pandas.read_csv(query)
    gbk_len = len(df_test)
    ref_len = len(df_ref)
    gbk_domains = []
    ref_gene = []
    query_gene = []
    query_name = []
    domain_list = []
    for line in range(0, gbk_len):
        x = 0
        domain = ''
        gene_test = df_test.iloc[line]
        #print(gene_test)
        #print(gene_test.isin([1.0]).any())
        for column in gene_test:
            #print(column)
            if column == 1.0 and x > 3:
                gbk_domains.append(df_test.columns[x])
                query_gene.append(gene_test[1])
            x += 1

    for row in range(0, ref_len):
        y = 0
        gene_ref = df_ref.iloc[row]
        potential_hits = []  # store for seq comparison
        for col in gene_ref:
            domain_test = ''
            if col == 1.0 and y > 3:
                domain_test = df_ref.columns[y]
                z = 0
                for domain_compare in gbk_domains:
                    #print(domain_compare)
                    #print(gbk_domains[z])
                    if domain_compare == domain_test:
                        #print(f'{domain_compare} {domain_test} {gene_ref[1]} {query_gene[z]}')
                        ref_gene.append(gene_ref[1])
                        potential_hits.append(query_gene[z])
                        #query_name.append(seq_compare(potential_hits, gene_ref[1]))
                        #query_name.append(query_gene[z])
                        domain_list.append(domain_compare) # wrong information going into df
                    z += 1
            y += 1
    print(query_name)
    print(ref_gene)
    print(len(query_name))
    print(len(ref_gene))
    df_cluster_compare = DataFrame({'Gene': query_name, 'Ref': ref_gene, 'Domain': domain_list}) # still not working and moving a bit too slow - change approach once it works
    print(df_cluster_compare)

def build_dataset(in_gbk, out_name, ref):
    file_path = Path(in_gbk)
    df_name = out_name
    for gbk in file_path.iterdir():
        if ref == 'True':
            domain_csv = 'gbk_domains.csv'
        else:
            domain_csv = 'gbk_query_domains.csv'
        BGC = str(gbk)
        BGC = BGC.replace('/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/MIBIG_TEST/', '')
        BGC = BGC.replace('/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/Query/', '')
        BGC = BGC.replace('.gbk', '')
        process(gbk, BGC, domain_csv)
        if ref == 'False':
            oh_ec('gbk_query_domains.csv', 'domains_encoded.csv')
            cluster_compare('domains_encoded.csv')
    if ref == 'True':
        oh_ec('gbk_domains.csv', df_name)

# fetch query cluster
def main():
    if os.path.exists('domains_ref_encoded.csv') is not True:
        in_gbk = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/MIBIG_TEST"
        build_dataset(in_gbk, 'domains_ref_encoded.csv', 'True')
    in_gbk = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/Query"
    build_dataset(in_gbk, 'domains_encoded.csv', 'False')

# add AA sequence to csv file for future comparisons



#    df = DataFrame({'protein_id': protein_id, 'Translation': protein, 'BGC':   , 'Position':   , DOMAINS})



if __name__ == "__main__":
    main()