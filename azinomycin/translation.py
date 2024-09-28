from Bio import SeqIO
from Bio.Seq import Seq

def read_fa(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(record)
    return sequences

def main():
    file_path = input("Enter the path to your FASTA file: ")

    # Read sequences from the FASTA file
    sequences = read_fasta_file(file_path)

    # Translate each DNA sequence to protein and print results
    for record in sequences:
        seq_id = record.id
        dna_seq = str(record.seq)
        dna_obj = Seq(dna_seq)
        protein_seq = dna_obj.translate()
        
        print(f"Sequence ID: {seq_id}")
        print(f"DNA Sequence: {dna_seq}")
        print(f"Protein Sequence: {protein_seq}")
        print("-" * 50)

if __name__ == "__main__":
    main()