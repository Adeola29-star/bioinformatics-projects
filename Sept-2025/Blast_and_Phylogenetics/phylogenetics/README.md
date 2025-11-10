## Phylogenetic Analysis of Hemoglobin Alpha Sequences


---

###  Introduction

Title: Phylogenetic Analysis of Hemoglobin Alpha

Description:
This notebook demonstrates the workflow for constructing a phylogenetic tree from protein sequences using Biopython. The analysis includes sequence retrieval, alignment, distance calculation, tree construction, and visualization.


---

### Fetching Sequences

Goal: Download protein sequences from NCBI Entrez.

Description:

Use the Entrez module to fetch sequences by accession numbers.

Save sequences in FASTA format for local processing.

Example sequences include human, mouse, chicken, frog, and zebrafish hemoglobin alpha proteins.



---

### Renaming Sequences

Goal: Add species information to sequence IDs.

Description:

Assign each sequence a unique ID with species prefix.

Update descriptions to include species names.

Save renamed sequences to a new FASTA file.



---

### Multiple Sequence Alignment

Goal: Align sequences to identify conserved regions.

Description:

Perform multiple sequence alignment using MUSCLE.

Save aligned sequences in FASTA format.

Alignments allow for downstream distance calculations and tree construction.



---

### Parsing Aligned Sequences

Goal: Load and inspect aligned sequences.

Description:

Use Biopython AlignIO to parse the aligned FASTA file.

Ensure all sequence IDs are unique.

Count total sequences and alignment length for reference.



---

### Distance Matrix Calculation

Goal: Measure sequence similarity for tree construction.

Description:

Use the identity model to calculate the fraction of identical positions between sequences.

Generate a distance matrix from the alignment.



---

### Phylogenetic Tree Construction

Goal: Build a Neighbor-Joining phylogenetic tree.

Description:

Use the distance matrix to construct a tree using the Neighbor-Joining method.

Inspect the tree structure in ASCII format in the terminal.



---

### Tree Visualization

Goal: Graphical representation of the phylogenetic tree.

Description:

Read the Newick tree file.

Plot the tree using matplotlib for clear visualization.

Save the tree as a high-resolution PNG for reports or presentations.



---

### Outputs

Files generated in this workflow:

hemoglobin_alpha_renamed.fasta — sequences with species info

Aligned_hemoglobin_alpha_renamed.fasta — multiple sequence alignment

hemoglobin_alpha_renamed_tree.nwk — Newick tree

hemoglobin_tree.jpg — graphical tree


```python

```
