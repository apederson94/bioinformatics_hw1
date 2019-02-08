#Created by Austin Pederson 02/08/2019 10:34

#Question 1
def read_fasta(file_name, text=True):
    with open(file_name, "r") as fasta:
        content = fasta.readlines()
        header = content[0]
        seq = ""

        for string in content[1:]:
            seq = seq + string.strip("\n")

        if text:
            fasta_dict = {"header": content[0], "sequence": seq}
        else:
            fasta_dict = {"header": bytes(content[0], 'utf-8'), "sequence": bytes(seq, 'utf-8')}
                
    return fasta_dict

def calculate_gc_content(dict):
    seq = dict["sequence"]

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
    mer_len = 0
    seq_is_str = type(sequence) == str
    if seq_is_str:
        mer = ""
    else:
        mer = bytearray()

    for num in range(k-1):
        sequence = sequence[num:]
        for symbol in sequence:

            if seq_is_str:
                mer += symbol
            else:
                mer.append(symbol)

            mer_len += 1

            if mer_len == k:
                if not seq_is_str:
                    mer = bytes(mer)
                if(mer_dict.get(mer)):
                    mer_dict[mer] += 1
                else:
                    mer_dict[mer] = 1

                mer_len = 0

                if seq_is_str:
                    mer = ""
                else:
                    mer = bytearray()

    return mer_dict
    
def main():
    seq_003997_bin = read_fasta("NC_003997.fasta", False)
    seq_003997_str = read_fasta("NC_003997.fasta")
    seq_004722_bin = read_fasta("NC_004722.fasta", False)
    seq_004722_str = read_fasta("NC_004722.fasta")

    gc_003997 = calculate_gc_content(seq_003997_bin)
    gc_004722 = calculate_gc_content(seq_004722_bin)

    k1 = k_mer(seq_003997_bin["sequence"], 4)
    k2 = k_mer(seq_004722_bin["sequence"], 4)

    


    



if __name__ == "__main__":
    main()