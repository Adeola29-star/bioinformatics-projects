# BioPython Practice – September 2025

This folder documents my hands-on learning and experiments with *BioPython*, covering key bioinformatics operations such as sequence manipulation, file parsing, alignments, and multiple sequence analysis (MSA).  
It builds upon my foundational biology and Python knowledge from previous modules.

---

## 📂 Folder Overview

| File | Description |
|------|--------------|
| *01_dna_sequence_basics.ipynb* | Basic Biopython usage for creating and manipulating Seq objects. Demonstrates complement, reverse complement, and reading sequences from FASTA files using SeqIO. |
| *02_create_fasta.ipynb* | Script for generating a FASTA file (globin.fasta) containing both mRNA and protein sequences for testing. |
| *03_msa_muscle.ipynb* | Multiple Sequence Alignment (MSA) using *MUSCLE* through a subprocess. Parses the aligned sequences with Biopython’s AlignIO module. |
| *04_my_biopython_project.ipynb* | Mini-project combining custom Python functions with Biopython utilities. Includes sequence cleaning, GC content calculation, transcription, translation, DataFrame creation (Pandas), and plotting (Matplotlib). |
| *05_pairwise_alignments.ipynb* | Demonstrates *global* and *local* alignments using PairwiseAligner, with custom match/mismatch/gap scores. Also integrates MUSCLE alignment runs for practice. |
| *06_protein_alignment_matrices.ipynb* | Compares alignment results using *BLOSUM62* and *PAM250* substitution matrices for protein sequences under both global and local alignment modes. |
| *data/* | Contains all related FASTA files (e.g., globin.fasta, clean.fasta, Aligned.fasta, etc.) used in the above scripts. |

---

## 🧠 Key Concepts Covered

- Working with *Seq* objects and SeqIO  
- Parsing and writing *FASTA* files  
- *GC content* and sequence cleaning  
- *DNA → RNA → Protein* transcription & translation  
- *Pairwise alignment* and *substitution matrices*  
- *Multiple sequence alignment* (MSA) with MUSCLE  
- Using *Pandas* and *Matplotlib* for analysis and visualization  

---

## ⚙️ Tools & Libraries Used

- *Python* 
- *BioPython*
- *Pandas*
- *Matplotlib*
- *Subprocess* (for MUSCLE)
- *MUSCLE* executable (installed locally)

---

## 📊 Example Workflow

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq

results_for_pandas= []
new_records= []
for seqRecord in SeqIO.parse("practice3.fasta", "fasta"):
    clean_sequence= clean(str(seqRecord.seq))
