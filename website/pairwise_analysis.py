from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def perform_pairwise_alignment(seq1, seq2, mode='global', seq_type='nucleotide'):
    try:
        aligner = PairwiseAligner()
        aligner.mode = mode  # Set alignment mode based on user input

        # Convert sequences to Seq objects based on sequence type
        if seq_type == 'protein':
            seq1 = Seq(seq1)
            seq2 = Seq(seq2)
        else:
            seq1 = Seq(seq1)
            seq2 = Seq(seq2)

        alignments = aligner.align(seq1, seq2)
        
        # Get the best alignment
        best_alignment = alignments[0]
        result = f"Score: {best_alignment.score}\n{best_alignment}"
    
        return [result]
    except Exception as e:
        print(f"Error: {e}")
        return ["Error: Unable to perform alignment."]