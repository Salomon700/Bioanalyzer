from flask import Blueprint, render_template, request, jsonify
from flask_login import login_required, current_user
from .pairwise_analysis import perform_pairwise_alignment
from .multiple_alignment import perform_multiple_alignment
from .phylogenetic_analysis import perform_phylogenetic_analysis

routes = Blueprint('routes', __name__)

@routes.route('/')
def home():
    return render_template("home.html")

@routes.route('/pairwise-alignment', methods=['GET', 'POST'])
@login_required
def pairwise_alignment():
    if request.method == 'POST':
        data = request.get_json()
        seq1 = data.get('seq1')
        seq2 = data.get('seq2')
        results = perform_pairwise_alignment(seq1, seq2)
        if "Error" in results[0]:
            return jsonify({"error": results[0]}), 500
        return jsonify(results)
    return render_template('pairwise_alignment.html')

@routes.route('/multiple-alignment', methods=['GET', 'POST'])
@login_required
def multiple_alignment():
    if request.method == 'POST':
        data = request.get_json()
        sequences = data.get('sequences')
        input_file = 'input.fasta'
        output_file = 'output.aln'

        with open(input_file, 'w') as f:
            f.write(sequences)

        alignment_result = perform_multiple_alignment(input_file, output_file)
        if "Error" in alignment_result:
            return jsonify({"error": alignment_result}), 500
        return jsonify({'alignment': alignment_result})
    return render_template('multiple_alignment.html')

@routes.route('/phylogenetic-tree', methods=['GET', 'POST'])
@login_required
def phylogenetic_tree():
    if request.method == 'POST':
        data = request.get_json()
        newick_data = data.get('newick')
        tree_ascii = perform_phylogenetic_analysis(newick_data)
        return jsonify({'phylogenetic_tree': tree_ascii})
    return render_template('phylogenetic_tree.html')

@routes.route('/health', methods=['GET'])
def health_check():
    return jsonify({"status": "UP"}), 200
