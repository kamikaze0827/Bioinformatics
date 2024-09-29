# BLASTP analysis

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

def read_fa(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(record)
    return sequences

def run_blastp(sequence):
    print(f"Submitting sequence {sequence.id} for BLASTP analysis...")
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence.seq)
    return result_handle

def filter_res(result_handle, e_threshold=0.01):
    blast_records = NCBIXML.parse(result_handle)
    filtered_results = []

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < e_threshold:
                    filtered_results.append({
                        'sequence_id': blast_record.query_id,
                        'hit_id': alignment.hit_id,
                        'hit_def': alignment.hit_def,
                        'e_value': hsp.expect,
                        'alignment_length': hsp.align_length,
                        'identity': hsp.identities,
                    })

    return filtered_results

def main():
    path = input("Enter the path to your FASTA file (with protein sequences): ")
    res = input("Enter the output file path for the filtered BLASTP report (e.g., 'blastp_report.txt'): ")
    sequences = read_fa(file_path)

    # Open output file to write filtered BLAST results
    with open(res, 'w') as out_file:
        out_file.write("Sequence_ID\tHit_ID\tHit_Definition\tE-value\tAlignment_Length\tIdentity\n")

        # Submit each sequence to BLASTP, filter and write results
        for record in sequences:
            result_handle = run_blastp(record)
            filtered_results = filter_res(result_handle)

            for result in filtered_results:
                out_file.write(f"{result['sequence_id']}\t{result['hit_id']}\t{result['hit_def']}\t{result['e_value']:.4e}\t{result['alignment_length']}\t{result['identity']}\n")

            print(f"Filtered results for sequence {record.id} saved.")

    print(f"BLASTP analysis complete. Filtered results saved to {res}")

if __name__ == "__main__":
    main()