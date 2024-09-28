# Provide basic information about the DNA sequence such as length, GC content, number of 'N' base call.

dfrom Bio import SeqIO
from prettytable import PrettyTable

def read_fa(file_path)
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(record)
    return sequences

def gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = (gc_count / len(sequence)) * 100
    return gc_content

def count_n_bases(sequence):
    return sequence.count('N')

def create_table(sequences):
    table = PrettyTable()
    table.field_names = ["Sequence ID", "Length", "GC Content (%)", "N Base Count"]

    for record in sequences:
        seq_id = record.id
        sequence = str(record.seq)
        seq_len = len(sequence)
        gc_content = gc_content(sequence)
        n_count = count_n_bases(sequence)
        table.add_row([seq_id, seq_len, f"{gc_content:.2f}", n_count])
    
    return table

def main():
    file_path = input("Enter the path to your FASTA file: ")

    sequences = read_fa(file_path)
    summary_table = create_table(sequences)

    # Print the summary table
    print(summary_table)

if __name__ == "__main__":
    main()
