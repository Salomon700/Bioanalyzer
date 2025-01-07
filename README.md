# Bioinformatics Web App

## Overview

This Bioinformatics Web App provides tools for performing pairwise sequence alignment, multiple sequence alignment, and phylogenetic tree construction. The app also includes user authentication to manage access to these features.

## Features

- **User Authentication**: Secure login and registration system.(Functional)
- **Pairwise Sequence Alignment**: Align two sequences to identify regions of similarity. (Functional)
- **Multiple Sequence Alignment**: Align multiple sequences to identify conserved regions. (Not yet Functional)
- **Phylogenetic Tree Construction**: Generate phylogenetic trees from sequence data. (Not yet Functional)

## Technologies Used

- **Flask**: Web framework for Python.
- **SQLAlchemy**: ORM for database management.
- **Flask-Login**: User session management.
- **Clustal Omega**: Tool for multiple sequence alignment.
- **Bootstrap**: Front-end framework for responsive design.

## Installation

### Prerequisites

- Python 3.x
- Flask
- Flask-Login
- Flask-Migrate
- Flask-SQLAlchemy
- biopython

### Steps

1. **Clone the repository**:

   ```sh
   git clone https://github.com/yourusername/bioinformatics-web-app.git
   cd bioinformatics-web-app
   ```

2. **Set up a virtual environment**:

   ```sh
   python3 -m venv venv
   source venv/bin/activate  # On Windows use `venv\Scripts\activate`
   ```

3. **Install the required packages**:

   ```sh
   pip install -r requirements.txt
   ```

4. **Set up the database**:

   ```sh
   flask db init
   flask db migrate -m "Initial migration."
   flask db upgrade
   ```

5. **Run the application**:
   ```sh
   python main.py
   ```

## Usage

1. **Access the web app**:
   Open your web browser and go to `http://127.0.0.1:5000`.

2. **Register and log in**:
   Create a new account or log in with an existing account.

3. **Use the tools**:
   - Navigate to the Pairwise Alignment, Multiple Sequence Alignment, or Phylogenetic Tree pages to use the respective tools.

## Acknowledgements

- [Flask](https://flask.palletsprojects.com/)
- [Bootstrap](https://getbootstrap.com/)
- [Biopython](https://biopython.org/)
- [Tech with Tim](https://github.com/techwithtim/Flask-Web-App-Tutorial)
