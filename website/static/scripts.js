// Perform Pairwise Alignment
async function performPairwiseAlignment() {
  const seq1 = document.getElementById("sequence1").value;
  const seq2 = document.getElementById("sequence2").value;

  try {
    const response = await fetch("/pairwise-alignment", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({ seq1, seq2 }),
    });

    if (!response.ok) {
      const errorData = await response.json();
      document.getElementById("results").textContent = errorData.error;
      return;
    }

    const result = await response.json();
    document.getElementById("results").textContent = result.join("\n\n");
  } catch (error) {
    document.getElementById("results").textContent = `Error: ${error.message}`;
  }
}

// Perform Multiple Sequence Alignment
async function performMultipleAlignment() {
  const sequences = document.getElementById("sequences").value;

  // Clear previous results
  document.getElementById("msa-results").textContent = "Loading...";

  try {
    const response = await fetch("/multiple-alignment", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sequences }),
    });

    if (!response.ok) {
      const errorData = await response.json();
      document.getElementById("msa-results").textContent = errorData.error;
      return;
    }

    const data = await response.json();
    const resultsDiv = document.getElementById("msa-results");
    if (data.alignment) {
      resultsDiv.innerText = data.alignment;
    } else {
      resultsDiv.innerText =
        "Error: Unable to perform alignment. Please check your input.";
    }
  } catch (error) {
    document.getElementById(
      "msa-results"
    ).textContent = `Error: ${error.message}`;
  }
}

// Generate Phylogenetic Tree
async function generatePhylogeneticTree() {
  const sequences = document.getElementById("tree-sequences").value;

  const response = await fetch("/phylogenetic-tree", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ newick: sequences }),
  });

  const data = await response.json();
  const resultsDiv = document.getElementById("tree-results");
  if (data.phylogenetic_tree) {
    resultsDiv.innerText = data.phylogenetic_tree;
  } else {
    resultsDiv.innerText =
      "Error: Unable to generate tree. Please check your input.";
  }
}
