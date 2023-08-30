


class Genome_hit:

    hit_list = []

    def __repr__(self):
        return f"Genome_hit('{self.name}', {self.hits})"

    def __init__(self, name, hits):
        self.name = name
        self.hits = hits

def sort_list(result):
    lst = list()
    pos = 0
    with open(result) as file:
        for line in file:
                line = line.split()
                if line[3] == "DETECTED":
                    line[3] = 0
                lst.append(Genome_hit(line[0], int(line[3])))
                # print(f"{line[0]} {line[3]}")
                pos = pos + 1
        lst.sort(key=lambda a: (a.hits, a.name))
    filename = f"Sorted_output.fasta"
    with open(filename, 'w') as file:
        for output in range(1, pos):
            # print(lst[output])
            file.write(f'{lst[output]}\n')
    print(f"COMPLETE")
    return lst


lst = sort_list(result="Hits_to_sort.txt")



