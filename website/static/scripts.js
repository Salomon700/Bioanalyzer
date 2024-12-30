// Pairwise alignment
async function performPairwiseAlignment() {
  const seq1 = document.getElementById("sequence1").value;
  const seq2 = document.getElementById("sequence2").value;
  const alignmentMode = document.getElementById("alignmentMode").value;

  try {
    const response = await fetch("/pairwise-alignment", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        seq1: seq1,
        seq2: seq2,
        mode: alignmentMode,
      }),
    });

    if (!response.ok) {
      const errorData = await response.json();
      document.getElementById("results").innerHTML = errorData.error;
      return;
    }

    const result = await response.json();
    document.getElementById("results").innerHTML = `
      <p><strong>Score:</strong> ${result.score}</p>
      <p><strong>Percent Identity:</strong> ${result.percent_identity}%</p>
      <p><strong>Alignment:</strong></p>
      <pre>${result.alignment}</pre>
    `;
  } catch (error) {
    document.getElementById("results").innerHTML = `Error: ${error.message}`;
  }
}

//sequence alignment ncbi
async function fetchSimilarSequences() {
  const query = document.getElementById("sequenceQuery").value;
  const resultsDiv = document.getElementById("similar-sequences-results");
  const spinner = document.getElementById("loading-spinner");
  const downloadButton = document.getElementById("download-button");

  // Show the loading spinner
  spinner.classList.remove("d-none");
  resultsDiv.innerHTML = "";
  downloadButton.classList.add("d-none");

  try {
    const response = await fetch("/fetch-similar-sequences", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({ query }),
    });

    // Hide the loading spinner
    spinner.classList.add("d-none");

    if (!response.ok) {
      const errorData = await response.json();
      resultsDiv.innerHTML = errorData.error;
      return;
    }

    const results = await response.json();
    resultsDiv.innerHTML = "";

    let sequencesText = "";

    results.forEach((result, index) => {
      const resultHtml = `
        <div class="result-item">
          <h5>Result ${index + 1}</h5>
          <p><strong>Title:</strong> ${result.title}</p>
          <p><strong>Accession:</strong> ${result.accession}</p>
          <p><strong>E-value:</strong> ${result.e_value}</p>
          <p><strong>Percent Similarity:</strong> ${result.percent_similarity.toFixed(
            2
          )}%</p>
          <pre>${result.sequence}</pre>
        </div>
        <hr>
      `;
      resultsDiv.innerHTML += resultHtml;

      sequencesText += `>Result ${index + 1}\n`;
      sequencesText += `Title: ${result.title}\n`;
      sequencesText += `Accession: ${result.accession}\n`;
      sequencesText += `E-value: ${result.e_value}\n`;
      sequencesText += `Percent Similarity: ${result.percent_similarity.toFixed(
        2
      )}%\n`;
      sequencesText += `${result.sequence}\n\n`;
    });

    // Show the download button
    downloadButton.classList.remove("d-none");

    // Store the sequences text in the button's dataset
    downloadButton.dataset.sequences = sequencesText;
  } catch (error) {
    // Hide the loading spinner in case of an error
    spinner.classList.add("d-none");
    resultsDiv.innerHTML = `Error: ${error.message}`;
  }
}

function downloadSequences() {
  const downloadButton = document.getElementById("download-button");
  const sequencesText = downloadButton.dataset.sequences;

  const blob = new Blob([sequencesText], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "sequences.txt";
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// Perform Multiple Sequence Alignment

// Generate Phylogenetic Tree
