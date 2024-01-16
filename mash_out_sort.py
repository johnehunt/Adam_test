from pandas import DataFrame

#use bigscape for classification of clusters - classify reference or query? Both?

# add check to only count the relevant part of the output (deleted start for now) - maybe include a print in the sh that can be used as reference?
# need tofix the final data so the correct info is copied
# if I do all of the MIBIG at once how do I tell what each one is?
# how do I process so much data?
# do all vs reference and all vs all

mash_out = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial/output.txt"
with open(mash_out, "r") as file:
    cluster_test = []
    reference_clusters = []
    distance_calc = []
    P_value_calc = []
    matching_hashes_calc = []
    for line in file:
        test = ''
        target_copy = False
        reference_copy = False
        mash_dist = False
        P_value = False
        matching_hashes = False
        P_value_off = False
        mash_dist_off = False
        target_cluster = ''
        reference_cluster = ''
        P_value_copy = ''
        mash_dist_copy = ''
        matching_hashes_copy = ''
        for letter in line:
            test = test + letter
            if letter != '	' and matching_hashes == True:
                matching_hashes_copy = matching_hashes_copy + letter
            if letter != '	' and P_value == True:
                P_value_copy = P_value_copy + letter
                P_value_off = True
            if letter == '	' and P_value_off == True:
                P_value = False
                matching_hashes = True
            if letter != "	" and mash_dist == True:
                mash_dist_copy = mash_dist_copy + letter
                mash_dist_off = True
            if letter == "	" and mash_dist_off == True:
                mash_dist = False
                if len(P_value_copy) < 1:
                    P_value = True
            if test == "fasta/":
                target_copy = True
            if test == 'MIBIG_fasta/':
                reference_copy = True
            if letter == "	":
                if target_copy == True:
                    target_cluster = test
                    target_cluster = target_cluster[:-7]
                if reference_copy == True:
                    reference_cluster = test
                    reference_cluster = reference_cluster[:-17]
                    mash_dist = True
                target_copy = False
                reference_copy = False
            if letter == "/":
                test = ''
        cluster_test.append(target_cluster)
        reference_clusters.append(reference_cluster)
        distance_calc.append(mash_dist_copy)
        P_value_calc.append(P_value_copy)
        matching_hashes_copy = matching_hashes_copy[:-1]
        matching_hashes_calc.append(matching_hashes_copy)

    #print(f'cluster: {cluster_test}\n MIBIG reference: {reference_clusters}\n Mash-distance: {distance_calc}\n P-value: {P_value_calc}\n Matching-hashes: {matching_hashes_calc}')
    df = DataFrame({'cluster': cluster_test, 'MIBIG reference': reference_clusters, 'Mash-distance': distance_calc, 'P-value': P_value_calc, 'Matching-hashes': matching_hashes_calc})
    print(df)
    #df.to_excel('test.xlsx', sheet_name='sheet1', index=False)
    df.to_csv('mash_result.csv', sep='\t')