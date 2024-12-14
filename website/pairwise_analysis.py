from Bio.Align import PairwiseAligner

def perform_pairwise_alignment(seq1, seq2):
    try:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Choose 'global' or 'local' alignment mode
        alignments = aligner.align(seq1, seq2)
        
        # Get the best alignment
        best_alignment = alignments[0]
        result = f"Score: {best_alignment.score}\n{best_alignment}"
    
        return [result]
    except Exception as e:
        print(f"Error: {e}")
        return ["Error: Unable to perform alignment."]