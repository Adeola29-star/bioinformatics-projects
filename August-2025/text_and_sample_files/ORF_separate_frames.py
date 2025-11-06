#!/usr/bin/env python
# coding: utf-8

# Exercise: ORFs in All 3 Reading Frames
# 
# Goal:
# Write a Python function that finds all ORFs in all three forward reading frames (frame 0, 1, and 2) of a DNA sequence.
# 
# Details:
# 
# Use your existing find_orfs function as a base.
# 
# Modify it so it accepts a parameter for the frame (0, 1, or 2).
# 
# For each frame, scan the DNA starting at that frame position (i.e., start at index 0 for frame 0, index 1 for frame 1, index 2 for frame 2).
# 
# Collect and print all ORFs found in each frame separately, showing start and end positions relative to the full sequence (1-based indexing).

# In[1]:


def fasta_parsing(practice3txt):
    fasta_dict= {}
    with open("practice3.txt", "r") as file:
        header= None
        sequence_lines= []
        for line in file:
            line= line.strip()
            if line.startswith('>'):
                if header:
                    fasta_dict[header]= "".join(sequence_lines)
                header= line[1:]
                sequence_lines= []

            else:
                sequence_lines.append(line)
        if header:
            fasta_dict[header]= "".join(sequence_lines)
    return fasta_dict

def clean(seq):
    cleaned_sequence= "".join([base for base in seq.upper() if base in 'ATGC'])
    return cleaned_sequence



def find_orf_frames(dna):
    
    start_codon= 'ATG'
    stop_codons= {'TAG', 'TGA', 'TAA'}
    orfs= []

    for frame in range (3):
        i= frame
        while i < len(dna)-2:
            codon= dna[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(dna)-2, 3):
                    stop_codon= dna[j:j+3]
                    if stop_codon in stop_codons:
                        orf= dna[i:j+3]
                        orfs.append((frame+1, i+1, j+3, orf))
                        i= j+3
                        break
                else:
                    i += 3
            else:
                i += 3
    return orfs


sequences= fasta_parsing("practice3.txt")

clean_sequences= {}
for header, seq in sequences.items():
    clean_sequences[header]= clean(seq)

for header, seq in clean_sequences.items():
    print(f"Header: {header}")
    orf_frames= find_orf_frames(seq)
    if orf_frames:
        for frame, s_index, e_index, orf_sequence in orf_frames:
            print(f"Frame_index: {frame}\nStart index: {s_index}\nEnd index: {e_index}\n",
              f"ORF sequence: {orf_sequence}\n")
    else:
        print("No ORFs found.\n")
          


# In[ ]:





# In[ ]:




