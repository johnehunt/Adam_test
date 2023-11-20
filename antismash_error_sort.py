import os


with open("antismash_error.txt", 'r') as file:
    missed_genomes = []
    for line in file:
        target = False
        copy = False
        word = ''
        word_record = ''
        turn_off = False
        for letter in line:
            word = word + letter
            if letter == ' ' and len(word_record) > 0:
                # missed_genomes = list.append(missed_genomes, word_record)
                missed_genomes.append(word_record)
            if letter == ' ':
                word = ''
            if len(word) > 0 and copy == True:
                word_record = word_record + letter
                turn_off = True
            if word == 'ERROR':
                target = True
            if word == 'Record' and target == True:
                copy = True
            if turn_off == True and letter == ' ':
                copy = False
                word_record = ''
    print(missed_genomes)
    print(missed_genomes[0])
    print(len(missed_genomes))


directory = '/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/HTP_antismash/genomes'
files = os.listdir(directory)
for docs in files:
    word = ''
    name = False
    prepare = False
    name_record = ''
    file_edit = ''
    file_path = os.path.join(directory, docs)
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            first_line = file.readline()
            for letter in first_line:
                word = word + letter
                if word == 'LOCUS       ':
                    name = True
                if prepare == True and letter == ' ':
                    name = False
                if name == True and letter == ' ':
                    prepare = True
                if name == True and letter != ' ':
                    name_record = name_record + letter
            name_record = f'{name_record}.1'
            for ref in missed_genomes:
                if name_record == ref:
                    file_edit = file_path
        word_check = ''
        genome_edited = False
        if len(file_edit) > 0:
            for letter in file_edit:
                word_check = word_check + letter
                if letter == '/':
                    word_check = ''
                if word_check == 'genomes':
                    genome_edited = True
                if genome_edited == True and letter == ' ':
                    genome_name = word_check
            print(word_check)
    except UnicodeDecodeError:
        print(f"Error reading {file}. Ignoring and moving to the next file.")
        continue
