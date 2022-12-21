with open("rnap.out", 'r') as file:
    # search = file.readlines()
    gene = False
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
            gene = True
            print(line_count)
            gene_line = line_count + 2
            print(gene_line)
            print("success")
            pos = 0
            for word in line:
                subject = subject + word
                region.append(pos) # region = [region + pos]
                if word == " ":
                    print(subject)
                    print(region)
                    subject = ""
                    region = []
                if subject == "Sequence":
                    additional = pos + 2
                    region.append(additional)
                    print(region)
                    print("check")
                    word_region = region
                    print(word_region)
                pos = pos + 1
        if line_count == gene_line:
            print(line)
            print(region)
            print(word_region)
            x = 0
            gene = ""
            for word in line:
                for item in word_region:
                    if x == item:
                        gene = gene + word
                        print("Yes")
                x = x + 1
        print(gene)
    if gene == False:
        print("Failed")