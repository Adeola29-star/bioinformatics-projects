#!/usr/bin/env python
# coding: utf-8

# In[2]:


# 13-08-2025
""" Load in a fasta file, parsed it, cleaned it up, and save the cleaned
dictionary in a file
"""
import random
def load(fasta_file):
    fasta_dict= {}
    with open(fasta_file, "r") as file:
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
    cleaned_seq= "".join([base for base in seq.upper() if base in 'ATGCU'])
    return cleaned_seq


def reverse_complement(seq):
    DNA_complement= {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    RNA_complement= {'A':'U', 'U':'A', 'G':'C', 'C':'G'}
    complement_seq=""
    
    if 'U' in seq:    # Detect sequence type
        seq_type= "RNA"
    else:
        seq_type= "DNA"
        
    for base in seq.upper():
        if seq_type == "RNA":
            complement_seq+= RNA_complement[base]
        else:
            complement_seq+= DNA_complement[base]
    reverse_complement_seq= complement_seq[::-1]
    return reverse_complement_seq


def gc_whole_sequences(seq):
    gc= 0
    total= len(seq)
    for base in seq.upper():
        if base in 'GC':
            gc += 1
    gc_content= round(gc/total * 100, 2) if total > 0 else 0
    return gc_content

def gc_sliding_window(seq, window_size):

    gc_content_window= []
    seq= seq.upper()
    for i in range(len(seq)- window_size + 1):
        window= seq[i:i+window_size]
        gc_window= 0
        for base in window:
            if base in 'GC':
                gc_window += 1
        gc_content_window.append(round(gc_window/window_size * 100, 2)) 
    return gc_content_window
            

def orf_finder(seq):
    
    start_codon= 'AUG'
    stop_codons= ['UAA', 'UGA', 'UAG']
    orfs= []
    longest_orf= ""
    longest_frame= None
    for frame in range(3):
        i= frame
        while i < len(seq)-2:
            codon= seq[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(seq)-2, 3):
                    stop_codon = seq[j:j+3]
                    if stop_codon in stop_codons:
                        orf= seq[i:j+3]
                        if len(orf) > len(longest_orf):
                            longest_orf= orf
                            longest_frame= frame + 1 # frames are numbered 1-3
                        orfs.append((frame+1, orf))
                        i= j+3
                else:
                    i += 3
            else:
                i += 3
    return orfs, longest_orf, longest_frame

def transcription(seq):
    rna= ""
    for base in seq:
        if base == 'T':
            rna += 'U'
        else:
            rna += base
    return rna
    

def translation(seq):
    
    codon_map = {
    'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
    'UUA': 'Leucine', 'UUG': 'Leucine',
    'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
    'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine',
    'AUG': 'Methionine',  # Start codon
    'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',

    'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
    'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
    'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
    'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',

    'UAU': 'Tyrosine', 'UAC': 'Tyrosine',
    'CAU': 'Histidine', 'CAC': 'Histidine',
    'CAA': 'Glutamine', 'CAG': 'Glutamine',
    'AAU': 'Asparagine', 'AAC': 'Asparagine',
    'AAA': 'Lysine', 'AAG': 'Lysine',
    'GAU': 'Aspartic acid', 'GAC': 'Aspartic acid',
    'GAA': 'Glutamic acid', 'GAG': 'Glutamic acid',

    'UGU': 'Cysteine', 'UGC': 'Cysteine',
    'UGG': 'Tryptophan',
    'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
    'AGU': 'Serine', 'AGC': 'Serine',
    'AGA': 'Arginine', 'AGG': 'Arginine',
    'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine',
    'UAG': 'Stop', 'UGA':'Stop', 'UAA':'Stop' # stop codons
}


    aminoacids=[]
    codons= [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    for codon in codons:
        aa = codon_map.get(codon, '?')
        if aa == 'Stop':
            break
        aminoacids.append(aa)
    protein = "-".join(aminoacids)
    return protein


def unique_mutation(seq, n):
    seq_list= list(seq)
    bases= ('A','T','G','C')
    mutated_positions= set()
    while len(mutated_positions) < n:
        position= random.randint(0, len(seq)-1)

        if position not in mutated_positions:
            original_base= seq_list[position]

            possible_outcomes=[]
            for base in bases:
                if base != original_base:
                    possible_outcomes.append(base)
            new_base= random.choice(possible_outcomes)
            seq_list[position] = new_base
            mutated_positions.add(position)
    mutated_sequence= "".join(seq_list)
    return mutated_sequence

def single_point_mutation(seq, pos_num, new_base):
    seq_list= list(seq)
    index= pos_num - 1
    if pos_num > len(seq) or new_base not in 'ATGC':
        return "Invalid information"
    seq_list[index] = new_base
    mutated_seq= "".join(seq_list)
    return mutated_seq, index
        
if __name__ == "__main__":
    
    fasta_file= "practice3.txt"
    sequences= load("practice3.txt")
    clean_dict= {}
    for header, seq in sequences.items():
        clean_dict[header]= clean(seq)

    with open("cleaned_fasta_file.fasta", "w") as file:
        for header, seq in clean_dict.items():
            file.write(f">{header}\n{seq}\n")
    
  
    for header, seq in clean_dict.items():
        print(f">{header}\n{seq}\n")

        gc_whole= gc_whole_sequences(seq)
        print(f"GC Content(Whole): {gc_whole}\n")

        gc_s_window= gc_sliding_window(seq, window_size=5)
        print(f"GC Content(Window): {gc_s_window}\n")

        r_complement= reverse_complement(seq)
        print(f"Reverse Complement: {r_complement}\n")

        rna = transcription(seq)
        print(f"RNA: {rna}\n")
        
        _, longest_orf, _= orf_finder(rna)
        print(f"Longest ORF:{longest_orf}\n")

        protein= translation(longest_orf)
        print(f"Protein: {protein}\n")

        mutated = unique_mutation(seq, 3)
        print(f"Mutated Sequence: {mutated}\n")

        single_pm= single_point_mutation(seq, 5, "C")
        print(f"Single Point Mutation: {single_pm}\n")
        
            


# In[ ]:





# In[ ]:




