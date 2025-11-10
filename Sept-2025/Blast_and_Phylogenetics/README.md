# Phylogenetics, BLAST & Sequence Analysis

### Overview

This folder contains scripts and data for sequence analysis, including:

- Phylogenetic analysis of protein sequences


- Offline BLASTX analysis of multiple sequences


- Mini-Project pipeline for sequence cleaning, GC calculation, transcription, translation, ORF detection, and alignments



The scripts use Biopython, Pandas, Matplotlib, and MUSCLE for multiple sequence alignment.


---

### Contents

- Phylogenetic Analysis

Retrieves protein sequences from NCBI (Entrez)

Renames sequences with species information

Performs multiple sequence alignment (MSA) using MUSCLE

Computes distance matrix and constructs a phylogenetic tree (Neighbor Joining)

Visualizes tree using Matplotlib and saves as PNG

Saves tree in Newick format for future use


- BLAST Analysis

Performs offline BLASTX on multiple sequences using a local SwissProt database

Parses XML results into a Pandas DataFrame

Extracts hits, alignment length, E-values, and partial alignments

Queries with no hits are recorded as None


- Mini-Project: Sequence Analysis Pipeline

Parses and cleans raw FASTA sequences

Computes sequence metrics:

Sequence length

GC content (%)

RNA sequence (transcription)

Protein sequence (translation)

Longest ORF per reading frame


Stores results in a Pandas DataFrame

Plots GC content distribution using Matplotlib

Performs pairwise alignment of sequences

Performs multiple sequence alignment (MSA) using MUSCLE

---

### Requirements

Python 

Biopython

Pandas

Matplotlib

MUSCLE executable (for multiple sequence alignment)



---

### Notes

Place supporting data files in their respective subfolders (Phylogenetic, BLAST, Mini_Project)

All scripts are modular and can be adapted for different sequences or projects

Ensure MUSCLE executable paths are correctly set in the scripts


```python

```
