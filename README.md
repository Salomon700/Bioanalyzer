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
- SQLAlchemy
- Flask-Login
- Clustal Omega

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

4. **Install Clustal Omega**:

   - **Ubuntu**:
     ```sh
     sudo apt-get update
     sudo apt-get install clustalo
     ```
   - **macOS**:
     ```sh
     brew install clustal-omega
     ```
   - **Windows**:
     - Download the Clustal Omega executable from [Clustal Omega download page](http://www.clustal.org/omega/).
     - Extract the executable to a directory and add the directory to your system's PATH.

5. **Set up the database**:

   ```sh
   flask db init
   flask db migrate -m "Initial migration."
   flask db upgrade
   ```

6. **Run the application**:
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

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgements

- [Clustal Omega](http://www.clustal.org/omega/)
- [Flask](https://flask.palletsprojects.com/)
- [Bootstrap](https://getbootstrap.com/)
