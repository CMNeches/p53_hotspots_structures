#!/usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#load in the data
df = pd.read_csv('../data/Germline.csv')
df = df[df["Short_topo"] == "COLON"]

#we are defining the codon array (Python starts counting at 0)
codon_array = np.arange(1,394)

#df is the dataset that we have (ProtDescription is a column, we have already filtered for missense), and we are converting that column (or series) into a numpy array.  This shows the slow way
#all_mutations = df["ProtDescription"].to_numpy()
 
#here is the fast way
all_codons = df["Codon_number"].to_numpy()
print(all_codons)

#this is called a list comprehension (everything within the square brackets), which is a one line for loop.  For example, when c=1, the code is all_codons ==1, which means that we are looking for all of the places within all of the codons where the value is equal to one.  .sum() is the sum of all of the comparisons
mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])

#here is how you loop through the mutation types
types = ["missense", "FS", "nonsense", "silent", "splice", "other", "intronic", "large del"]
for t in types:
    temporary_df = df[df["Effect"] == t]
    codon_array = np.arange(1,394)
    all_codons = temporary_df["Codon_number"].to_numpy()
    mutation_counts = np.array([(all_codons == c).sum() for c in codon_array])
    
    domain_dict = {"TAD1" : [1, 41, 'b'], 
                   "TAD2" : [41, 62, 'r'],
                   "PR" : [65, 93, 'g'],
                   "DBD" : [100, 300, 'c'],
                   "OD": [323, 355, 'm'],
                   "C-terminus" : [364, 393, 'y']}

    #discovering novel hotspot locations
    #mutation counts is the number of mutations at each codon, and then codon array is an array from 1 to 393, argsort sorts an array in ascending order.  That is why we sort on the negative of the mutation count array (we take the top 10 mutations)
    for k,v in domain_dict.items():
        start, end, _ = domain_dict[k]
        new_codon_array = codon_array[start : end]
        new_mutation_counts = mutation_counts[start : end]
        hotspots = new_codon_array[(-new_mutation_counts).argsort()][:10]
        freqs = (new_mutation_counts[(-new_mutation_counts).argsort()]/new_mutation_counts.sum())[:10]
    
        print(t, k, hotspots)    
        print(t, k, freqs)
        print()

