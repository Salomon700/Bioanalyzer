from flask import Blueprint, render_template, request, jsonify, flash, redirect
from flask_login import login_required, current_user
from .pairwise_analysis import perform_pairwise_alignment, fetch_top_similar_sequences

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
        mode = data.get('mode', 'global')
        results = perform_pairwise_alignment(seq1, seq2, mode)
        if "error" in results[0]:
            return jsonify(results[0]), 500
        return jsonify(results[0])
    return render_template('pairwise_alignment.html')

@routes.route('/fetch-similar-sequences', methods=['GET', 'POST'])
@login_required
def fetch_similar_sequences():
    if request.method == 'POST':
        data = request.get_json()
        query = data.get('query')
        database = data.get('database', 'nt')
        results = fetch_top_similar_sequences(query, db=database)
        if "error" in results:
            return jsonify(results), 500
        return jsonify(results)
    return render_template('fetch_similar_sequences.html')

@routes.route('/about')
def about():
    return render_template('about.html')

@routes.route('/contact')
def contact():
    return render_template('contact.html')