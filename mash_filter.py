import pandas as pd
import csv

mash_path = "/Users/u2186477/Documents/PhD/Year 1/Python/cblaster/project_reference/attempt_2/hmmer-3.3.1/Adam_test/mash_trial"

input_path = f'{mash_path}/mash_result_ref_control.csv'
output_path = f'{mash_path}/mash_result_filtered_control.csv'

with open(input_path, 'r') as file:
    cvs_reader = csv.reader(file)
    count = 0
    numerical = 0
    for row in cvs_reader:
        if count == 0:
            with open(output_path, 'w') as out_file:
                writer = csv.writer(out_file)
                writer.writerow(row)
        if 0 < count: # if 0 < count < 1000:
            record = False
            decimal_points = ''
            decimal_count = 0
            first_point = True
            less_than_one = False
            for point in row[3]:
                if first_point == True and point == '0':
                    first_point = False
                    less_than_one = True
                if first_point == True and point == '1':
                    first_point = False
                if record == True and point != '-':
                    decimal_points = decimal_points + point
                if point == 'e':
                    record = True
                    decimal_count_done = decimal_count
                decimal_count = decimal_count + 1
            if record == True:
                int(decimal_points)
                decimal_point = int(decimal_points) + int(decimal_count_done) - 2
                #print(f'{row[3]} {decimal_point} {decimal_points} {decimal_count_done}')
                start_text_start = '"{:.'
                start_text_end = 'f}"'
                start_text = f'{start_text_start}{decimal_point}{start_text_end}'
                #print("{:.20f}".format(float(row[3])))
                #print("{:.11f}".format(float(row[3])))
                #print(f'{start_text}.format(float({row[3]}))')
                #print(start_text.format(float(row[3])))
                numerical = start_text.format(float(row[3]))
                numerical = numerical.replace('"', '')
            if record == False and less_than_one == False:
                numerical = int(row[3])
            if less_than_one == True and record == False:
                numerical = float(row[3])
            # if not an int convert and do to decimal determined by the number after e
            # print("{:.20f}".format(float(row[3])))
            #if float(numerical) != 1: # less than or not 1?
            if float(numerical) < float(0.0000001): # e-7 because this is just beyond 1/1000
            #if float(numerical) < float(0.00000000000001):  # e-15
            #if float(numerical) < float(0.0000000000000000001):  # e-20
                print('========================================')
                print(f'cluster name is {row[0]}')
                print(f'reference cluster is {row[1]}')
                print(f'p-value score {row[3]} or {numerical}')
                print(f'similarity is {row[4]}')
                with open(output_path, 'a') as out_file:
                    writer = csv.writer(out_file)
                    writer.writerow(row)
        count = count + 1



#data = pd.read_csv(output_path)
# df =  df[df.name != "1"]
# print "{:.16f}".format(float("1.70000043572e-05"))
#writer.writerow(row)
#denote reference clusters


# use controls made up of partial clusters - track real and percieved VAN interactions on figure
# create network of reference clusters - might be complicating results