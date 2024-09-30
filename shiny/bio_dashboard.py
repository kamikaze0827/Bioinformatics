# Create a Shiny dashboard to summarize bioinformatics analysis

import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from shiny import App, reactive, ui, render

# Function to extract general information of the DNA sequence
def gen_info(sequences):
    metrics = []
    for record in sequences:
        seq_length = len(record.seq)
        gc_content = (record.seq.count('G') + record.seq.count('C')) / seq_length * 100 if seq_length > 0 else 0
        n_count = record.seq.count('N')
        
        metrics.append({
            'Sequence_ID': record.id,
            'Length': seq_length,
            'GC_Content': gc_content,
            'N_Count': n_count
        })
    return pd.DataFrame(metrics)

# Function to run BLASTN analysis 
def blastn(fasta_file, e_threshold=0.01):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    results = []

    for record in sequences:
        print(f"Submitting sequence {record.id} for BLASTN analysis...")
        result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < e_threshold:
                        results.append({
                            'Sequence_ID': blast_record.query_id,
                            'Hit_ID': alignment.hit_id,
                            'Hit_Definition': alignment.hit_def,
                            'E-value': hsp.expect,
                            'Alignment_Length': hsp.align_length,
                            'Identity': hsp.identities,
                        })
    
    return pd.DataFrame(results)

# Define the UI layout
app_ui = ui.page_fluid(
    ui.panel_title("BLASTN Analysis Summary Dashboard"),
    
    # File upload input
    ui.input_file("fasta_file", "Upload FASTA file", multiple=False, accept=[".fasta", ".fa"]),
    
    # Action button to run BLAST analysis
    ui.input_action_button("run_blast", "Run BLAST Analysis"),
    
    # Data table to show results
    ui.output_table("blast_results"),
    
    # Show summary statistics
    ui.output_text("summary_stats")
)

# Define the server-side logic
def server(input, output, session):
    
    # Placeholder for the filtered results
    filtered_blast_results = reactive.Value(None)
    
    # Initialize the output for blast_results
    @output
    @render.table
    def blast_results():
        res = filtered_blast_results.get()
        if res is not None and not res.empty:
            return res
        else:
            return pd.DataFrame()
        
    # Initialize the output for summary_stats
    @output
    @render.text
    def summary_stats():
        res = filtered_blast_results.get()
        if res is not None and not res.empty:
            num_sequences = len(res)
            avg_identity = res["Identity"].mean() 
            return f"Number of Sequences: {num_sequences}, Average Identity: {avg_identity:.2f}"
        else:
            return "No valid results to display."
    
    # Function to run BLAST analysis when the button is clicked
    @reactive.Effect
    @reactive.event(input.run_blast)
    def on_blast_run():
        fasta_file = input.fasta_file()
        
        if fasta_file:
            print("Running BLAST analysis and generating summary...")

            fasta_path = fasta_file[0]['datapath']
            blast_res = blastn(fasta_path)

            if isinstance(blast_res, pd.DataFrame):
                filtered_blast_results.set(blast_res)
            else:
                print("The result is not a valid DataFrame.")
                filtered_blast_results.set(pd.DataFrame())
        else:
            filtered_blast_results.set(pd.DataFrame())
    
# Create the Shiny App
app = App(app_ui, server)