import os
from pathlib import Path
import shutil

# delete all empty directories (with the correct name)
def remove_empty_dirs(x):
    x = int(x)
    file = f'output_directory{x}'
    file_path = Path(f'HTP_antismash/{file}')
    if os.path.exists(file_path):
        #print(file)
        if len(os.listdir(file_path)) != 0:
            return(x)
        if len(os.listdir(file_path)) == 0:
            print(f'Removing {file}')
            os.rmdir(file_path)

def remove_dirs(file):
    file_path = Path(f'HTP_antismash/{file}')
    if os.path.exists(file_path):
        #print(file)
        #os.rmdir(file_path)
        shutil.rmtree(file_path, ignore_errors=True)

# delete all GBK files
def remove_gbff(file):
    if os.path.exists(file):
        print(file)
        os.remove(file)

def manage_output(file):
    clusters = Path('HTP_antismash/clusters')
    file_path = Path(f'HTP_antismash/{file}')
    for out in file_path.iterdir():
        #print(out)
        region_check = ''
        region_true = False
        for letter in str(out):
            region_check = region_check + letter
            if region_check == 'region':
                region_true = True
            if letter == '.':
                region_check = ''
        if region_true == True:
            shutil.move(out, clusters)




# move antismash data to single file


if __name__ == "__main__":

    antismash_folder = Path('HTP_antismash')
    genomes_success = []

    x = 0
    for file in antismash_folder.iterdir():
        to_save = remove_empty_dirs(x)
        genomes_success.append(to_save)
        x = x + 1

    antismash_gbff = Path('HTP_antismash/genomes')
    for file in antismash_gbff.iterdir():
        remove_gbff(file)

    genomes_success_filtered = list(filter(lambda item: item is not None, genomes_success))
    genomes_success_len = len(genomes_success_filtered)
    genomes_success_end = genomes_success_filtered[(genomes_success_len - 1)]

    y = 0
    #for file in antismash_folder.iterdir():
    for file in range(0, (genomes_success_end+1)):
        for pos in range(0, (genomes_success_len)):
            if y == genomes_success_filtered[pos]:
                x = genomes_success_filtered[pos]
                x = x
                folder_itr = f'output_directory{x}'
                print(folder_itr)
                target_file = (f'{antismash_folder}/{folder_itr}')
                if os.path.exists(target_file):
                    manage_output(folder_itr)
                if os.path.exists(target_file):
                    remove_dirs(folder_itr)
        y = y + 1



