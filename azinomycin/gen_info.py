# Provide basic information about the DNA sequence such as length, GC content, number of 'N' base call.

def read_dna(file_path):
    with open(file_path, 'r') as file:
        seq = file.read().replace('\n', '').upper()
    return seq

def main():
    file_path = input("Enter the path to your sequence file: ")
    seq = read_dna(file_path)

    # Calculate length of DNA
    dna_len = len(seq)

    # Calculate GC content
    gc_count = seq.count('G') + seq.count('C')
    gc_content = (gc_count / len(seq)) * 100

    # Calculate number of 'N' base call
    n_count = seq.count('N')
  
    # Print the results
    print(f"DNA Sequence Length: {dna_len} bp")
    print(f"GC Content: {gc_content:.3f}%")
    print(f"Count of 'N' bases: {n_count}")

if __name__ == "__main__":
    main()
