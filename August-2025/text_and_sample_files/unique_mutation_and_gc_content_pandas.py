#!/usr/bin/env python
# coding: utf-8

# In[11]:


import pandas as pd
import random

def fasta_parsing(practice3txt):
    fasta_dict= {}
    with open("practice3.txt") as file:
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
    cleaned_seq= "".join([base for base in seq.upper() if base in 'ATCG'])
    return cleaned_seq

def calculate_gc(seq):
    gc= 0
    for base in seq.upper():
        if base in 'GC':
            gc += 1
    gc_percent= round(gc/len(seq) * 100, 2) if len(seq) > 0 else 0
    return gc_percent

def unique_mutation(seq, n):
    seq_list= list(seq)
    bases= ('A', 'T', 'C', 'G')
    mutated_position= set()
    while len(mutated_position) < n:
        position= random.randint(0, len(seq)-1)
        if position not in mutated_position:
            original_base= seq_list[position]
            
            possible_outcomes= []
            for base in bases:
                if base != original_base:
                    possible_outcomes.append(base)
            new_base= random.choice(possible_outcomes)
            seq_list[position]= new_base
            mutated_position.add(position)
    mutated_sequence= "".join(seq_list)
    return mutated_sequence

sequences= fasta_parsing("practice3.txt")
n = int(input("Enter the numbers of mutations you want to occur:"))

cleaned_sequences= {}
for header, seq in sequences.items():
    cleaned_sequences[header]= clean(seq)


data= []
for header, seq in cleaned_sequences.items():
    mutated_sequence= unique_mutation(seq, n)
    data.append({
        "Header": header,
        "Sequence": seq,
        "Length": len(seq),
        "GC_Content(%)": calculate_gc(seq),
        "Mutated_Sequence": mutated_sequence,
        "Mutated_GC_Content(%)": calculate_gc(mutated_sequence)
    })
        
            
df= pd.DataFrame(data)
print(df)
            
    
    


# In[ ]:




