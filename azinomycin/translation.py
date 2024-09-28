from Bio.Seq import Seq

def read_dna(file_path):
    with open(file_path, 'r') as file:
        seq = file.read().replace('\n', '').upper()
    return seq

def translate(dna_seq):
    dna_obj = Seq(dna_seq)
    protein_seq = dna_obj.translate()
    return protein_seq

def main():
    file_path = input("Enter the path to your DNA sequence file: ")
    seq = read_dna(file_path)

    # Translate to protein sequence
    protein_seq = translate(seq)

    # Print the result
    print(f"Translated Amino Acid Sequence:\n{protein_seq}")

if __name__ == "__main__":
    main()
