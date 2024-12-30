from Bio import Entrez, SeqIO, Seq
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "ntwalisolomon9@gmail.com"  # Replace with your email

def fetch_top_similar_sequences(query, db="nt", retmax=10):
    try:
        # Perform BLAST search
        result_handle = NCBIWWW.qblast("blastn", db, query)
        
        # Parse BLAST results
        blast_records = NCBIXML.read(result_handle)
        
        sequences = []
        for alignment in blast_records.alignments[:retmax]:
            for hsp in alignment.hsps:
                sequence_data = {
                    "title": alignment.title,
                    "sequence": hsp.sbjct,
                    "accession": alignment.accession,
                    "e_value": hsp.expect,
                    "percent_similarity": (hsp.identities / hsp.align_length) * 100
                }
                sequences.append(sequence_data)
        
        # Sort sequences by percent similarity in descending order
        sequences.sort(key=lambda x: x['percent_similarity'], reverse=True)
        
        return sequences
    except Exception as e:
        return {"error": str(e)}

def calculate_percent_identity(alignment):
    matches = sum(1 for a, b in zip(alignment[0], alignment[1]) if a == b)
    length = min(len(alignment[0]), len(alignment[1]))
    return (matches / length) * 100

def color_code_alignment(alignment):
    html_alignment = ""
    for a, b in zip(alignment[0], alignment[1]):
        if a == b:
            html_alignment += f'<span style="color: green;">{a}</span>'
        else:
            html_alignment += f'<span style="color: red;">{a}</span>'
    return html_alignment

def perform_pairwise_alignment(seq1, seq2, mode='global'):
    try:
        aligner = PairwiseAligner()
        aligner.mode = mode  # Set alignment mode based on user input

        # Convert sequences to Seq objects
        seq1 = Seq(seq1)
        seq2 = Seq(seq2)

        alignments = aligner.align(seq1, seq2)
        
        # Get the best alignment
        best_alignment = alignments[0]
        percent_identity = calculate_percent_identity(best_alignment)
        colored_alignment = color_code_alignment(best_alignment)
        
        result = {
            "score": best_alignment.score,
            "alignment": colored_alignment,
            "percent_identity": percent_identity
        }
    
        return [result]
    except Exception as e:
        print(f"Error: {e}")
        return [{"error": "Unable to perform alignment."}]