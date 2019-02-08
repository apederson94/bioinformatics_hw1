#Created by Austin Pederson 02/08/2019 10:34
import random

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

def k_mer(sequence, k):
    mer_dict = {}
    seq_len = len(sequence)

    for pos in range(seq_len):
        mer = sequence[pos:pos+k]
        if(mer_dict.get(mer)):
            continue
        else:
            mer_dict[mer] = sequence.count(mer)

    return mer_dict

#TODO: this
def jaccard_dist(mer_dict1, mer_dict2):
    union = []
    intersection = []
    for mer in mer_dict1:

        if mer_dict2.get(mer):
            intersection.append(mer)

        if mer not in union:
            union.append(mer)

    #print(union)
    #print("\n")
    #print(intersection)
    
def read_samples(sequence, mer_length, read_count):
    samples_dict = {}
    sample_len = len(sequence)
    for num in range(read_count):
        pos = random.randint(0, sample_len)
        header = str(num + 1) + "|" + str(pos)
        samples_dict[header] = sequence[pos:pos+mer_length]
    print(samples_dict)

def main():
    seq_003997_bin = read_fasta("NC_003997.fasta", False)
    seq_003997_str = read_fasta("NC_003997.fasta")
    seq_004722_bin = read_fasta("NC_004722.fasta", False)
    seq_004722_str = read_fasta("NC_004722.fasta")

    header_003997_bin = next(iter(seq_003997_bin.keys()))
    header_003997_str = next(iter(seq_003997_str.keys()))
    header_004722_bin = next(iter(seq_004722_bin.keys()))
    header_004722_str = next(iter(seq_004722_str.keys()))
    
    calculate_gc_content(seq_003997_bin[header_003997_bin])

    k1 = k_mer(seq_003997_bin[header_003997_bin], 4)
    k2 = k_mer(seq_003997_str[header_003997_str], 4)


    #jaccard_dist(k1, k2)

    #read_samples(seq_003997_bin[header_003997_bin], 45, 25)



if __name__ == "__main__":
    main()