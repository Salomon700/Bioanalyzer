import subprocess
import os

def perform_multiple_alignment(input_file: str, output_file: str):
    # Run Clustal Omega for MSA using subprocess
    command = [
        "clustalo",  # Ensure 'clustalo' is in your PATH
        "--infile", input_file,
        "--outfile", output_file,
        "--verbose",
        "--auto"
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stderr:
        raise Exception(f"ClustalOmega Error: {result.stderr}")
    with open(output_file, "r") as f:
        return f.read()