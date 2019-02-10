#Created by Austin Pederson 02/08/2019 10:34
import random
import datetime
import csv
import re

#Question 1
def read_fasta(file_name, text=True):
    if text:
        with open(file_name, "r") as fasta:
            content = fasta.readlines()
            header = content[0].strip('\n')
            seq = ""
            for tmp in content[1:]:
                seq += tmp.strip('\n')
            fasta_dict = {header: seq}

    else:
        with open(file_name, "rb") as fasta:
            content = fasta.read().splitlines()
            seq = bytearray()
            for tmp in content[1:]:
                seq += tmp
            fasta_dict = {content[0]: bytes(seq)}

    return fasta_dict


#Question 2
def calculate_gc_content(seq):

    if type(seq) == str:
        g = seq.count("G")
        c = seq.count("C")
        a = seq.count("A")
        t = seq.count("T")
    else:
        g = seq.count(b'G')
        c = seq.count(b'C')
        a = seq.count(b'A')
        t = seq.count(b'T')

    return (g + c)/(a + t + c + g)


#Question 3
def k_mer(sequence, k):
    mer_dict = {}
    seq_len = len(sequence)

    for pos in range((seq_len-k)+1):
        mer = sequence[pos:pos+k]
        if(not mer_dict.get(mer)):
            mer_dict[mer] = 1
        else:
            mer_dict[mer] += 1

    return mer_dict


#Question 4
def jaccard_dist(mer_dict1, mer_dict2):
    union = {}
    intersection = {}
    int_count = 0
    un_count = 0

    for mer in mer_dict1:

        if mer_dict2.get(mer):
            count = min(mer_dict1[mer], mer_dict2[mer])
            intersection[mer] = count
            int_count += count

        if mer not in union:
            count = max(mer_dict1[mer], mer_dict2[mer])
            union[mer] = count
            un_count += count

    return int_count/un_count


#Question 5
def read_samples(sequence, mer_length, read_count):
    samples_dict = {}
    sample_len = len(sequence)
    for num in range(read_count):
        pos = random.randint(0, sample_len)
        header = str(num + 1) + "|" + str(pos)
        samples_dict[header] = sequence[pos:pos+mer_length]


def main():
    seq_004722_bin = read_fasta("NC_004722.fasta", False)
    seq_003997_bin = read_fasta("NC_003997.fasta", False)
    seq_004722_str = read_fasta("NC_004722.fasta")

    header_004722_bin = next(iter(seq_004722_bin.keys()))
    header_003997_bin = next(iter(seq_003997_bin.keys()))
    header_004722_str = next(iter(seq_004722_str.keys()))

    k1 = k_mer(seq_003997_bin[header_003997_bin], 4)
    k2 = k_mer(seq_004722_bin[header_004722_bin], 4)

    jaccard_dist(k1, k2)

"""     #Question 6

    #part1
    with open("q6_p1.csv", "w") as csvfile:
        writer = csv.writer(csvfile, delimiter = ",")
        writer.writerow(["sequence_name", "text", "replicate", "seconds"])
        
        for i in range(10):
            tstart = datetime.datetime.now()
            read_fasta("NC_004722.fasta", False)
            t_bin = datetime.datetime.now() - tstart

            tstart = datetime.datetime.now()
            #read_fasta("NC_003997.fasta")
            t_str = datetime.datetime.now() - tstart
    
            writer.writerow([header_004722_bin, False, i, t_bin])
            writer.writerow([header_004722_str, True, i, t_str])

    #part2
    with open("q6_p2.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sequence_name", "text", "replicate", "k"])
        for k in [4, 8, 12]:
            for i in range(10):
                tstart = datetime.datetime.now()
                k_mer(seq_004722_bin[header_004722_bin], k)
                t_bin = datetime.datetime.now() - tstart

                tstart = datetime.datetime.now()
                k_mer(seq_004722_str[header_004722_str], k)
                t_str = datetime.datetime.now() - tstart
                
                writer.writerow([header_004722_bin, False, i, t_bin, k])
                writer.writerow([header_004722_str, True, i, t_str, k])
            

   #part3
    with open("q6_p3.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sequence_name", "text", "replicate", "mer_length", "read_count"])
        for read_count in [100000, 500000, 1000000, 10000000, 25000000, 50000000]:
            for i in range(10):
                tstart = datetime.datetime.now()
                read_samples(seq_004722_bin[header_004722_bin], 50, read_count)
                t_bin = datetime.datetime.now() - tstart

                tstart = datetime.datetime.now()
                read_samples(seq_004722_str[header_004722_str], 50, read_count)
                t_str = datetime.datetime.now() - tstart
                
                writer.writerow([header_004722_bin, False, i, t_bin, 50, read_count])
                writer.writerow([header_004722_str, True, i, t_str, 50, read_count]) """



if __name__ == "__main__":
    main()