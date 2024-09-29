# Merge multiple sequence .fasta file into one single .fasta file

from Bio import SeqIO

def merge_fa(fasta_files, output):
    with open(output, 'w') as out_file:
        for file in fasta_files:
            with open(file, 'r') as in_file:
                for record in SeqIO.parse(in_file, 'fasta'):
                    SeqIO.write(record, out_file, 'fasta')

def main():
    input = input("Enter the paths of the FASTA files to merge, separated by commas: ").split(',')
    res = input("Enter the output FASTA file name (e.g., merged_seq.fasta): ")
    
    merge_fa(input, res)

if __name__ == "__main__":
    main()